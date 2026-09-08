program dust_backend_smoke
  use snrt_dust_ir
  use snrt_runtime_backend
  use mpi_mod
  use iso_c_binding, only: c_float
  implicit none
  integer,parameter::nc=1031,ng=2,nd=2
  type(dust_ir_table)::table
  type(dust_ir_diagnostics)::ref_diag,trial_diag
  real(dust_dp)::rays(3,nd),weights(nd),rho(nc),heat(nc),cap(nc)
  real(dust_dp)::ref(ng,nd,nc),trial(ng,nd,nc),ref_t(nc),trial_t(nc),ref_p(ng,nc),trial_p(ng,nc)
  real(dust_dp)::ref_e(nc),trial_e(nc),error
  real(dust_dp)::ghosts(ng,nd,1)
  integer::remote(6,nc)
  logical::blocked(6,nc)
  integer::links(6,nc),i,mode,ierr,info
  real(c_float)::photons(4,3,2),saved_photons(4,3,2)
  real(dust_dp)::scatter_weight(4),scatter_density(2),scatter_sigma(3),expected,original_sum
  integer::g,d
  real(dust_dp)::gas_e(2),gas_c(2),dust_e(2),dust_t(2),exchange(2),saved_dust(2)
  real(dust_dp)::coupled_gas(nc),coupled_cv(nc),coupled_k(nc),coupled_q(nc),old_gas(nc),old_dust(nc),total_change
  real(dust_dp)::cell_u(4,nc),test_power(4),test_band(2,4),test_rate(2,nc)
  real(dust_dp)::one_rate(2,1),one_t(1),one_e(1),one_q(1)
  call MPI_INIT(info)
  call snrt_backend_initialize(ierr)
  if(ierr/=0)then
     write(*,*)'BACKEND_INIT_REJECTED',ierr
     call MPI_ABORT(MPI_COMM_WORLD,2,info)
  endif
  gas_e=[8d-20,2d-20];gas_c=1d-21;dust_e=[2d-20,8d-20];saved_dust=dust_e
  dust_t=0;exchange=0
  call snrt_runtime_dust_exchange(gas_e,gas_c,dust_e,[1d3,1d3],[1d3,1d3], &
       [10d0,20d0,50d0,100d0],[1d-23,2d-23,5d-23,1d-22],1d-21,.5d0,1d8,10d0,dust_t,exchange,ierr)
  if(ierr/=0)stop 30
  if(exchange(1)<=0.or.exchange(2)>=0.or.any(gas_e-exchange<0))stop 31
  if(maxval(abs((dust_e-saved_dust)-exchange))>1d-32)stop 32
  if(maxval(abs(gas_e-exchange+dust_e-gas_e-saved_dust))>1d-32)stop 33
  saved_dust=dust_e
  call snrt_runtime_dust_exchange(gas_e,gas_c,dust_e,[1d3,1d3],[1d3,1d3], &
       [10d0,20d0,50d0,100d0],[1d-23,2d-23,5d-23,1d-22],-1d-21,.5d0,1d8,10d0,dust_t,exchange,ierr)
  if(ierr==0.or.any(dust_e/=saved_dust))stop 34
  write(*,*)'GAS_DUST_EXCHANGE_FORTRAN_CONSERVATION_ROLLBACK_PASS'
  scatter_weight=[1d0,2d0,3d0,4d0];scatter_density=[0d0,2d0];scatter_sigma=[0d0,.15d0,50d0]
  do i=1,2
     do g=1,3
        photons(:,g,i)=[real(i*g,c_float),0.0_c_float,0.0_c_float,0.0_c_float]
     enddo
  enddo
  saved_photons=photons
  call snrt_runtime_isotropic_scatter(photons,scatter_weight,scatter_density,scatter_sigma,1d0,ierr)
  if(ierr/=0)stop 20
  do i=1,2
     do g=1,3
        original_sum=sum(real(saved_photons(:,g,i),dust_dp))
        do d=1,4
           expected=real(saved_photons(d,g,i),dust_dp)*exp(-scatter_density(i)*scatter_sigma(g)) + &
                original_sum*scatter_weight(d)/10*(1-exp(-scatter_density(i)*scatter_sigma(g)))
           if(abs(real(photons(d,g,i),dust_dp)-expected)>3d-7)stop 21
        enddo
        if(abs(sum(real(photons(:,g,i),dust_dp))-original_sum)>3d-7)stop 22
     enddo
  enddo
  saved_photons=photons;scatter_density(2)=-1
  call snrt_runtime_isotropic_scatter(photons,scatter_weight,scatter_density,scatter_sigma,1d0,ierr)
  if(ierr==0.or.any(photons/=saved_photons))stop 23
  write(*,*)'SCATTER_FORTRAN_LAYOUT_ANALYTIC_ROLLBACK_PASS'
  rays(:,1)=[1d0,0d0,0d0];rays(:,2)=[-1d0,0d0,0d0];weights=.5d0;links=0
  remote=0;remote(1,1)=1;ghosts=1d-25;blocked=.false.;blocked(2,nc)=.true.
  do i=1,nc
     rho(i)=.5d0+real(i,dust_dp)/nc
     if(i>1)links(1,i)=i-1
     if(i<nc)links(2,i)=i+1
  enddo
  rho(1)=0;cap=1d-24
  do mode=0,1
     if(mode==0)then
        call snrt_dust_ir_initialize(table,[.001d0,.01d0],[.001d0,.01d0],[1d-21,1d-12], &
             [10d0,20d0,50d0,100d0],10d0,ierr)
        ref_e=20d0*cap
     else
        call snrt_dust_ir_initialize(table,[.001d0,.01d0],[.001d0,.01d0],[1d-21,1d-12], &
             [10d0,20d0,50d0,100d0],10d0,ierr,[1d-23,2d-23,8d-23,2d-22])
        ref_e=rho*2d-23
     endif
     if(ierr/=0)stop 3
     heat=1d-30*rho;ref=1d-25;ref_t=20;ref_p=0
     trial=ref;trial_t=ref_t;trial_p=ref_p;trial_e=ref_e
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          ref,ref_t,ref_p,ref_diag,ierr,1d-10,128,ref_e,cap,ghosts,remote,blocked)
     if(ierr/=0)stop 4
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap,ghosts,remote,blocked, &
          material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
          absorb_dispatch=snrt_runtime_ir_absorb)
     if(ierr/=0)stop 5
     error=max(maxval(abs(ref-trial))/maxval(ref),maxval(abs(ref_t-trial_t))/maxval(ref_t), &
          maxval(abs(ref_p-trial_p))/maxval(ref_p),maxval(abs(ref_e-trial_e))/maxval(ref_e))
     if(error>1d-10)stop 6
     if(ref_diag%escaped_erg/=trial_diag%escaped_erg.or.ref_diag%interface_erg/=trial_diag%interface_erg)stop 8
     if(trial_diag%balance_relative>1d-10.or.trial_diag%local_relative>1d-10)stop 9
     write(*,'(A,I0,A,ES14.6)')'DUST_BACKEND_FORTRAN_PARITY material_u=',mode,' relative=',error
     ! Error after entering the material callback, not just outer shape checks.
     heat=1d99
     ref=trial;ref_e=trial_e;ref_t=trial_t;ref_p=trial_p
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap,ghosts,remote,blocked, &
          material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
          absorb_dispatch=snrt_runtime_ir_absorb)
     if(ierr==0.or.any(trial/=ref).or.any(trial_e/=ref_e).or.any(trial_p/=ref_p).or.any(trial_t/=ref_t))stop 7
  enddo
  ! A stiff gas reservoir can cool through dust into IR without driving the
  ! finite-capacity grain beyond its table. The split post-IR kick fails here.
  coupled_cv=1d-24;coupled_gas=1d-21;coupled_k=1d-29*rho;coupled_q=0
  heat=1d-30*rho;trial=1d-25;trial_t=20;trial_p=0;trial_e=rho*2d-23
  old_gas=coupled_gas;old_dust=trial_e;ref=trial
  call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
       trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap, &
       material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
       absorb_dispatch=snrt_runtime_ir_absorb,gas_energy=coupled_gas,gas_capacity=coupled_cv, &
       conductance=coupled_k,gas_transfer=coupled_q)
  if(ierr/=0)stop 40
  total_change=(sum(trial-ref)*.5d0+sum(trial_e-old_dust)+sum(coupled_gas-old_gas))*1d36 &
       +trial_diag%escaped_erg-sum(heat)*1d42
  if(abs(total_change)>1d-10*sum(old_gas)*1d36.or.any(coupled_q<0).or.sum(coupled_q)<=0)stop 41
  if(maxval(trial_t)>100.or.minval(coupled_gas)<0)stop 42
  do i=1,nc
     ! Independent variable-speed BE gas equation, with final (not initial) Tgas.
     expected=1d6*coupled_k(i)*sqrt(coupled_gas(i)/old_gas(i))
     expected=expected/(coupled_cv(i)+expected)*(old_gas(i)-coupled_cv(i)*trial_t(i))
     if(abs(expected-coupled_q(i))>1d-11*old_gas(i))stop 44
  enddo
  old_gas=coupled_gas;old_dust=trial_e;ref=trial;ref_t=trial_t;ref_p=trial_p
  coupled_k(nc)=-1
  call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
       trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap, &
       material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
       absorb_dispatch=snrt_runtime_ir_absorb,gas_energy=coupled_gas,gas_capacity=coupled_cv, &
       conductance=coupled_k,gas_transfer=coupled_q)
  if(ierr==0.or.any(coupled_gas/=old_gas).or.any(trial_e/=old_dust).or.any(trial/=ref).or. &
       any(trial_t/=ref_t).or.any(trial_p/=ref_p))stop 43
  write(*,*)'GAS_DUST_IR_VARIABLE_SPEED_JOINT_CONSERVATION_STIFF_ROLLBACK_PASS'
  ! Heterogeneous per-cell U(T), including cells across hybrid batch edges.
  test_power=[1d-30,2d-30,8d-30,2d-29]
  test_band(1,:)=.6d0*test_power;test_band(2,:)=.4d0*test_power
  do i=1,nc
     cell_u(:,i)=[1d-23,2d-23,8d-23,2d-22]*(.5d0+real(i,dust_dp)/nc)
  enddo
  heat=1d-30*rho;old_dust=rho*cell_u(2,:);cap=1d-24
  old_gas=8d-23;coupled_cv=1d-24;coupled_k=1d-29*rho
  do mode=0,1
     if(mode==0)then
        call snrt_runtime_dust_material(heat,rho,old_dust,cap,log([10d0,20d0,50d0,100d0]), &
             test_power,test_band,cell_u(:,1),.true.,1d6,test_power(1),10d0,1d-10, &
             test_rate,trial_t,trial_e,ierr,cell_material_u=cell_u)
     else
        call snrt_runtime_dust_material(heat,rho,old_dust,cap,log([10d0,20d0,50d0,100d0]), &
             test_power,test_band,cell_u(:,1),.true.,1d6,test_power(1),10d0,1d-10, &
             test_rate,trial_t,trial_e,ierr,old_gas,coupled_cv,coupled_k,coupled_q,cell_u)
     endif
     if(ierr/=0)then
        write(*,*)'CELL_MATERIAL_TEST_FAILED mode/status=',mode,ierr
        stop 50
     endif
     do i=2,nc,257
        if(mode==0)then
           call snrt_runtime_dust_material(heat(i:i),rho(i:i),old_dust(i:i),cap(i:i), &
                log([10d0,20d0,50d0,100d0]),test_power,test_band,cell_u(:,i),.true.,1d6, &
                test_power(1),10d0,1d-10,one_rate,one_t,one_e,ierr)
        else
           call snrt_runtime_dust_material(heat(i:i),rho(i:i),old_dust(i:i),cap(i:i), &
                log([10d0,20d0,50d0,100d0]),test_power,test_band,cell_u(:,i),.true.,1d6, &
                test_power(1),10d0,1d-10,one_rate,one_t,one_e,ierr, &
                old_gas(i:i),coupled_cv(i:i),coupled_k(i:i),one_q)
           if(abs(one_q(1)-coupled_q(i))>1d-12*old_gas(i))stop 51
        endif
        if(ierr/=0.or.abs(one_t(1)-trial_t(i))>1d-11*trial_t(i))stop 52
        if(maxval(abs(one_rate(:,1)-test_rate(:,i)))>1d-11*maxval(test_rate(:,i)))stop 53
        if(abs(one_e(1)-trial_e(i))>1d-11*trial_e(i))stop 54
     enddo
  enddo
  write(*,*)'DUST_CELL_MATERIAL_HETEROGENEOUS_HYBRID_LAYOUT_PASS'
  call check_d03_optics()
  write(*,*)'DUST_BACKEND_PARITY_AND_ROLLBACK_PASS'
  call MPI_FINALIZE(info)
contains
  subroutine check_d03_optics()
    use dust_composition_optics
    use dust_composition_material, only: dust_composition_curve
    use snrt_dust_contract
    type(dust_ir_table)::local_table
    type(dust_ir_diagnostics)::local_diag,other_diag
    character(len=2048)::contract_path
    real(dust_dp)::bins(4),pa(d03_ng),ps(d03_ng),pg(d03_ng),depth(d03_ng)
    real(dust_dp)::ia(d03_nir),isc(d03_nir),ig(d03_nir),pure(d03_ng,4)
    real(dust_dp)::spectra(d03_nir,5),pa_changed(d03_ng)
    real(dust_dp)::u(4),tgrid(4),mass_ref,total_before,total_after,parity
    real(dust_dp)::field(d03_nir,2,1),other(d03_nir,2,1),emitted(d03_nir,1),other_emitted(d03_nir,1)
    real(dust_dp)::td(1),ed(1),other_td(1),other_ed(1),heating(1),dust_density(1)
    integer::case_id,status,neighbors(6,1)
    call get_command_argument(1,contract_path)
    if(len_trim(contract_path)==0)return ! Existing no-argument smoke remains valid.
    call snrt_dust_contract_load(trim(contract_path),status)
    if(status/=0)stop 60
    if(.not.d03_optics_binding(snrt_dust_contract_group_edges_ev(1:10), &
         snrt_dust_contract_absorption_mean_energy_ev(1:9), &
         snrt_dust_contract_ir_energy_ev(1:d03_nir),d03_radius_cm,d03_solid_density))stop 61
    pa_changed=d03_primary_ev;pa_changed(9)=pa_changed(9)*1.001d0
    if(d03_optics_binding(d03_edges,pa_changed,d03_ir_ev,d03_radius_cm,d03_solid_density))stop 62
    call d03_mix_opacity([1d0,0d0,0d0,0d0],[5d-7,1d-5],d03_solid_density,1d0, &
         pa,ps,pg,ia,isc,ig,status)
    if(status==0)stop 63
    call d03_mix_opacity([1d0,0d0,0d0,0d0],d03_radius_cm,[2.2d0,3.3d0],1d0, &
         pa,ps,pg,ia,isc,ig,status)
    if(status==0)stop 64
    call d03_mix_opacity([0d0,0d0,0d0,0d0],d03_radius_cm,d03_solid_density,1d0, &
         pa,ps,pg,ia,isc,ig,status)
    if(status/=0.or.any(pa/=0).or.any(ia/=0))stop 65
    tgrid=[10d0,20d0,50d0,100d0];mass_ref=snrt_dust_contract_mass_per_h_g
    dust_density=1d8;heating=1d-30*dust_density;neighbors=1 ! closed reciprocal cell
    do case_id=1,5
       bins=0
       if(case_id<=4)then
          bins(case_id)=1
       else
          bins=[.1d0,.2d0,.3d0,.4d0]
       endif
       call d03_mix_opacity(bins,d03_radius_cm,d03_solid_density,mass_ref,pa,ps,pg,ia,isc,ig,status)
       if(status/=0.or.any(pa<=0).or.any(ps<=0).or.any(abs(pg)>1).or.any(ia<=0))stop 66
       if(pg(9)<.999d0)stop 67 ! Hard X rays are NOT isotropic.
       if(case_id<=4)pure(:,case_id)=pa
       if(case_id==5)then
          if(maxval(abs(pa-matmul(pure,bins))/pa)>1d-13)stop 68
       endif
       call d03_absorption_depth(bins*1d-26,d03_radius_cm,d03_solid_density,1d20,depth,status)
       if(status/=0.or.maxval(abs(depth-pa*(sum(bins)*1d-6/mass_ref))/depth)>1d-13)stop 69
       call dust_composition_curve(tgrid,[sum(bins(1:2)),sum(bins(3:4))],mass_ref,u,status)
       if(status/=0)stop 70
       ! The SAME per-mixture absorption cross section enters IR attenuation
       ! and Kirchhoff emission via the existing native table initializer.
       call snrt_dust_ir_initialize(local_table,d03_ir_ev,snrt_dust_contract_ir_weight_ev(1:d03_nir), &
            ia,tgrid,10d0,status,u)
       if(status/=0)stop 71
       field=0;emitted=0;td=20;ed=dust_density*u(2)
       other=field;other_emitted=emitted;other_td=td;other_ed=ed
       total_before=sum(ed)+1d6*sum(heating)
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,1d6,1d5,dust_density,heating, &
            field,td,emitted,local_diag,status,1d-10,128,ed,[1d0])
       if(status/=0)then
          write(*,*)'D03_REFERENCE_REJECTED',case_id,status
          stop 72
       endif
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,1d6,1d5,dust_density,heating, &
            other,other_td,other_emitted,other_diag,status,1d-10,128,other_ed,[1d0], &
            material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
            absorb_dispatch=snrt_runtime_ir_absorb)
       if(status/=0)stop 73
       total_after=sum(ed)+sum(field)*.5d0
       if(abs(total_after-total_before)>1d-10*total_before.or.local_diag%escaped_erg/=0)stop 74
       parity=max(maxval(abs(other-field))/maxval(field),maxval(abs(other_ed-ed))/maxval(ed), &
            maxval(abs(other_td-td))/maxval(td),maxval(abs(other_emitted-emitted))/maxval(emitted))
       if(parity>1d-10.or.other_diag%balance_relative>1d-10)stop 75
       spectra(:,case_id)=emitted(:,1)/sum(emitted)
       ! In the background-subtracted formulation, Td=Tbath, no primary and
       ! no excess IR is an exact stationary state, not spurious bath heating.
       field=0;td=10;ed=dust_density*u(1);emitted=0
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,1d6,1d5,dust_density,[0d0], &
            field,td,emitted,local_diag,status,1d-10,128,ed,[1d0], &
            material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
            absorb_dispatch=snrt_runtime_ir_absorb)
       if(status/=0.or.any(abs(ed-dust_density*u(1))>1d-12*ed))stop 76
       if(sum(field)>1d-12*sum(ed))stop 77
       write(*,'(A,I0,A,ES12.4)')'D03_NATIVE_IR_REFERENCE_BACKEND case=',case_id,' error=',parity
    enddo
    if(sum(abs(spectra(:,1)-spectra(:,3)))<.01d0)stop 78
    write(*,*)'D03_GRAIN_OPACITY_MASS_UNITS_BINDING_KIRCHHOFF_CONSERVATION_PASS'
    call check_d03_mixtures()
  end subroutine

  subroutine check_d03_mixtures()
    use dust_composition_optics
    use dust_composition_material, only: dust_composition_curve
    use snrt_dust_contract, only: snrt_dust_contract_ir_weight_ev,snrt_dust_contract_mass_per_h_g
    type(dust_ir_table)::mixed_table,single_table
    type(dust_ir_diagnostics)::diag,scalar_diag
    real(dust_dp)::pa(d03_ng,4),ps(d03_ng,4),pg(d03_ng,4),ia(d03_nir,4),isc(d03_nir,4),ig(d03_nir,4)
    real(dust_dp)::f(4,nc),uc(4,nc),tr(4),r(nc),h(nc),ed(nc),old_ed(nc),td(nc)
    real(dust_dp)::field(d03_nir,2,nc),emitted(d03_nir,nc),before(d03_nir,2,nc)
    real(dust_dp)::scalar_field(d03_nir,2,1),scalar_emitted(d03_nir,1),scalar_t(1),scalar_e(1)
    real(dust_dp)::sigma(d03_nir),gas(nc),cv(nc),k(nc),q(nc),old_g(nc),balance
    real(dust_dp)::ps_cell(d03_ng,nc),flux,alpha
    real(dust_dp)::ir_sigma(d03_nir,nc),ir_weight(2),mean_ir
    real(c_float)::beam(2,d03_ng,nc)
    integer::neighbors(6,nc),status,j,g,phase
    tr=[10d0,20d0,50d0,100d0]
    call d03_opacity_basis(snrt_dust_contract_mass_per_h_g,pa,ps,pg,ia,isc,ig,status)
    if(status/=0)stop 80
    do j=1,nc
       call d03_cell_weights([real(j,dust_dp),real(nc-j,dust_dp),3d2,1d2],f(:,j),status)
       call dust_composition_curve(tr,[sum(f(1:2,j)),sum(f(3:4,j))], &
            snrt_dust_contract_mass_per_h_g,uc(:,j),status)
       if(status/=0)stop 81
       neighbors(:,j)=j
    enddo
    call snrt_dust_ir_initialize(mixed_table,d03_ir_ev,snrt_dust_contract_ir_weight_ev(1:d03_nir), &
         ia(:,1),tr,10d0,status,uc(:,1),optical_sigma=ia)
    if(status/=0)stop 82
    r=1d8;r(1)=0;h=1d-30*r
    do phase=0,1
       ed=r*uc(2,:);old_ed=ed;td=20;field=0;before=field;emitted=0
       cv=1d-14;gas=cv*70;old_g=gas;k=r*1d-31;q=0
       if(phase==0)then
          call snrt_dust_ir_advance(mixed_table,rays,weights,neighbors,1d12,1d6,1d5,r,h,field,td, &
               emitted,diag,status,1d-10,128,ed,spread(1d0,1,nc), &
               material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
               absorb_dispatch=snrt_runtime_ir_absorb,cell_material_u=uc,cell_weights=f)
       else
          call snrt_dust_ir_advance(mixed_table,rays,weights,neighbors,1d12,1d6,1d5,r,h,field,td, &
               emitted,diag,status,1d-10,128,ed,spread(1d0,1,nc), &
               material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
               absorb_dispatch=snrt_runtime_ir_absorb,cell_material_u=uc,cell_weights=f, &
               gas_energy=gas,gas_capacity=cv,conductance=k,gas_transfer=q)
       endif
       if(status/=0)then
          write(*,*)'D03_MIXTURE_REJECT',phase,status
          stop 83
       endif
       balance=abs(sum(ed-old_ed)+sum(field)*.5d0+sum(gas-old_g)-1d6*sum(h))
       if(balance>1d-10*(sum(old_ed)+sum(old_g)+1d6*sum(h)))stop 84
       if(phase==0)then
          do j=2,nc,257
             sigma=matmul(ia,f(:,j))
             call snrt_dust_ir_initialize(single_table,d03_ir_ev,snrt_dust_contract_ir_weight_ev(1:d03_nir), &
                  sigma,tr,10d0,status,uc(:,j))
             if(status/=0)stop 85
             scalar_field=0;scalar_t=20;scalar_e=old_ed(j);scalar_emitted=0
             call snrt_dust_ir_advance(single_table,rays,weights,reshape([1,1,1,1,1,1],[6,1]),1d12,1d6,1d5, &
                  r(j:j),h(j:j),scalar_field,scalar_t,scalar_emitted,scalar_diag,status,1d-10,128,scalar_e,[1d0])
             if(status/=0.or.maxval(abs(scalar_field(:,:,1)-field(:,:,j)))>1d-10*maxval(field(:,:,j)))stop 86
             if(abs(scalar_e(1)-ed(j))>1d-10*ed(j).or.abs(scalar_t(1)-td(j))>1d-10*td(j))stop 87
          enddo
       else
          do j=1,nc
             alpha=1d6*k(j)*sqrt(gas(j)/old_g(j))
             flux=alpha/(cv(j)+alpha)*(old_g(j)-cv(j)*td(j))
             if(abs(flux-q(j))>1d-10*old_g(j))stop 88
          enddo
       endif
    enddo
    ps_cell=matmul(ps-pg,f);beam=0;beam(1,:,:)=1
    call snrt_runtime_isotropic_scatter(beam,weights,r,pa(:,1),1d16,status,ps_cell)
    if(status/=0)stop 89
    do j=1,nc
       do g=1,d03_ng
          if(abs(sum(real(beam(:,g,j),dust_dp))-1)>3d-7)stop 90
          flux=exp(-1d16*r(j)*ps_cell(g,j))
          if(abs(real(beam(1,g,j)-beam(2,g,j),dust_dp)-flux)>3d-7)stop 91
       enddo
    enddo
    ! FP64 IR layout differs from primary. Unequal quadrature weights test
    ! conservation of weighted energy, not an incorrect sum of intensities.
    ir_sigma=matmul(isc-ig,f);ir_weight=[.3d0,.7d0]
    field(:,1,:)=3d-15;field(:,2,:)=1d-15;before=field
    call snrt_runtime_ir_scatter(field,ir_weight,r,ir_sigma,1d16,status)
    if(status/=0)stop 93
    do j=1,nc
       do g=1,d03_nir
          mean_ir=sum(before(g,:,j)*ir_weight)
          flux=exp(-1d16*r(j)*ir_sigma(g,j))
          if(abs(sum(field(g,:,j)*ir_weight)-mean_ir)>2d-14*mean_ir)stop 94
          if(maxval(abs(field(g,:,j)-(mean_ir+(before(g,:,j)-mean_ir)*flux)))>2d-14*mean_ir)stop 95
       enddo
    enddo
    before=field;ir_sigma(1,nc)=-1
    call snrt_runtime_ir_scatter(field,ir_weight,r,ir_sigma,1d16,status)
    if(status==0.or.any(field/=before))stop 96
    ! Wrong mixture cannot alter the radiation/material transaction.
    before=field;old_ed=ed;old_g=gas;f(1,2)=-1
    call snrt_dust_ir_advance(mixed_table,rays,weights,neighbors,1d12,1d6,1d5,r,h,field,td, &
         emitted,diag,status,1d-10,128,ed,spread(1d0,1,nc), &
         material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
         absorb_dispatch=snrt_runtime_ir_absorb,cell_material_u=uc,cell_weights=f)
    if(status==0.or.any(field/=before).or.any(ed/=old_ed))stop 92
    write(*,*)'D03_HETEROGENEOUS_OPTICS_IR_GAS_HYBRID_FLUX_MOMENT_ROLLBACK_PASS'
    write(*,*)'D03_IR_TRANSPORT_SCATTER_WEIGHTED_ENERGY_ROLLBACK_PASS'
  end subroutine
end program

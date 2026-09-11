program dust_backend_smoke
  use dust_sublimation_material
  use dust_iron_material
  use dust_iron_radiation
  use dust_iron_optics
  use snrt_dust_ir
  use snrt_runtime_backend
  use dust_multifluid
  use mpi_mod
  use iso_c_binding, only: c_float
  use snrt_spectral_contract, only: snrt_d03_band_enabled
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
  real(dust_dp),parameter::fe_test_reference_mass=1.398d-26
  real(dust_dp)::moving_cold_bins(6),moving_cold_basis(fe_nir,fe_nt,6)
  call MPI_INIT(info)
  if(snrt_d03_band_enabled())then
     call check_d03_spectrum()
     call MPI_FINALIZE(info)
     stop
  endif
  call check_fe_cold_live_cell()
  call snrt_backend_initialize(ierr)
  if(ierr/=0)then
     write(*,*)'BACKEND_INIT_REJECTED',ierr
     call MPI_ABORT(MPI_COMM_WORLD,2,info)
  endif
  call check_primary_momentum()
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
          absorb_dispatch=snrt_runtime_ir_absorb,thin_reabsorption=.true.)
     if(ierr/=0)stop 5
     error=max(maxval(abs(ref-trial))/maxval(ref),maxval(abs(ref_t-trial_t))/maxval(ref_t), &
          maxval(abs(ref_p-trial_p))/maxval(ref_p),maxval(abs(ref_e-trial_e))/maxval(ref_e))
     if(error>1d-10)stop 6
     if(trial_diag%iterations>=ref_diag%iterations)stop 164
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
  call check_radiative_sublimation()
  call check_iron_radiation()
  call check_six_component_ir()
  call check_pah_mixed()
  write(*,*)'DUST_BACKEND_PARITY_AND_ROLLBACK_PASS'
  call MPI_FINALIZE(info)
contains
  subroutine check_pah_charged_mixed()
    use dust_pah_mixed
    use dust_pah_live_model
    use dust_mass_physics, only: dust_mass_enabled,dust_iron_model,dust_pah_model,dust_pah_nbin, &
         dust_pah_solid_charge,dust_pah_molecule_g
    use snrt_dust_contract
    integer,parameter::n=2,np=2*dust_pah_nbin
    type(dust_ir_table)::tab
    type(dust_ir_diagnostics)::diag
    real(dust_dp)::field(fe_nir,2,n),pop(np,n),old_pop(np,n),bins(6,n),primary(9,n)
    real(dust_dp)::ed(n),old_ed(n),td(n),ph(fe_nir,n),gas(n),old_gas(n),cv(n),q(n),ne(n),old_ne(n)
    real(dust_dp)::ghost(fe_nir,2,0),rays(3,2),weights(2),hi,residual,saved_field(fe_nir,2,n)
    real(dust_dp),parameter::electron_cv=1.5d0*1.380649d-16
    integer::neighbors(6,n),remote(6,n),status,i
    logical::blocked(6,n)
    dust_mass_enabled=.true.;dust_iron_model='none';dust_pah_model='pah_charge_fixed_h_v1'
    pop=0;pop(dust_pah_nbin+1,1)=dust_pah_molecule_g
    if(abs(dust_pah_solid_charge(pop(:,1),1.66d-24)/1.66d-24-1d0)>1d-14)stop 348
    call pah_live_prepare(status)
    if(status/=0.or.size(pah_level)/=np)stop 340
    call snrt_dust_ir_initialize(tab,snrt_dust_contract_ir_energy_ev(1:fe_nir), &
         snrt_dust_contract_ir_weight_ev(1:fe_nir),snrt_dust_contract_ir_absorption_per_h_cm2(1:fe_nir), &
         snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature),10d0,status)
    if(status/=0)stop 341
    bins=0;bins(1:4,1)=1d-28;ed=0;td=20
    call iron_mixture_enthalpy(20d0,[2d-28,2d-28,0d0],ed(1),hi,status)
    if(status/=0)stop 342
    pop=0;pop(1,:)=1d-6;field=0;ph=0;cv=1d-20;gas=300*cv;ne=1d-4
    primary=0;primary(3,1)=1d-20
    rays=0;rays(1,:)=[1d0,-1d0];weights=.5d0
    neighbors=0;neighbors(1:2,1)=2;neighbors(1:2,2)=1;remote=0;blocked=.false.
    old_pop=pop;old_ed=ed;old_gas=gas;old_ne=ne;diag=dust_ir_diagnostics()
    call pah_mixed_advance(tab,rays,weights,neighbors,1d12,1d0,3d8,bins,primary,field,pop,ed,td,ph,diag,status, &
         ghost,remote,blocked,gas,cv,[0d0,0d0],q,gas_electrons=ne,electron_capacity=electron_cv)
    if(status/=0.or.ne(1)<=old_ne(1).or.sum(pop(dust_pah_nbin+1:,1))<=0.or.sum(field)<=0)stop 343
    residual=sum(ed-old_ed)+sum(gas-old_gas)+.5d0*sum(field)-sum(primary)
    do i=1,n
       residual=residual+dot_product(pah_level,pop(:,i)-old_pop(:,i))
       if(abs(sum(pop(dust_pah_nbin+1:,i))-(ne(i)-old_ne(i)))>1d-17)stop 344
    enddo
    if(abs(residual)/sum(primary)>1d-8.or.any(abs(gas-old_gas+q)>1d-30))stop 345
    write(*,'(A,3ES22.14)')'PAH_CHARGED_MIXED balance,cations,gas_heat=', &
         residual/sum(primary),sum(pop(dust_pah_nbin+1:,1)),sum(gas-old_gas)
    old_ne=ne;old_pop=pop;old_gas=gas;old_ed=ed;saved_field=field
    primary=0;primary(9,1)=1d-20
    call pah_mixed_advance(tab,rays,weights,neighbors,1d12,1d0,3d8,bins,primary,field,pop,ed,td,ph,diag,status, &
         ghost,remote,blocked,gas,cv,[0d0,0d0],q,gas_electrons=ne,electron_capacity=electron_cv)
    if(status==0.or.any(ne/=old_ne).or.any(pop/=old_pop).or.any(gas/=old_gas).or.any(ed/=old_ed))stop 346
    if(any(field/=saved_field))stop 347
    write(*,*)'PAH_CHARGED_MIXED_PRIMARY_PARTITION_GAS_IR_ROLLBACK_PASS'
  end subroutine

  subroutine check_d03_spectrum()
    use iso_c_binding, only: c_int,c_double
    use dust_composition_optics, only: d03_band_ev,d03_band_abs,d03_band_transport,d03_edges
    use snrt_spectral_contract, only: snrt_fe_band_enabled,snrt_grain_band_bins
    use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment
    integer,parameter::cells=2,dirs=2,bands=9
    real(c_float)::q(cells,dirs,bands),direction(3,dirs),tau(cells,bands),stau(cells,bands,3)
    real(c_float)::dtau(cells,bands),budget(cells,3),hh(cells,bands,3),dd(cells,bands)
    real(c_float)::returned(cells,bands),raw(cells,bands),ag(cells,bands),at(cells),saved_q(cells,dirs,bands)
    real(c_double)::shift(cells,dirs,bands),reference(bands),he(cells,bands,3),de(cells,bands)
    real(c_double)::em(cells,bands,3),columns(cells,3),xi(cells),deposition(cells,8)
    real(c_double),allocatable::grains(:,:),kabs(:,:,:),ksca(:,:,:)
    real(c_double)::weights(dirs),before_n,before_e,after_n,after_e,mean,ta,ts,expected
    real(c_double)::saved_shift(cells,dirs,bands),nerror,eerror,maximum_n,maximum_e
    integer(c_int)::neighbors(6,cells),rc
    integer::mode,c,g,d,status,nb,last_mode
    logical::iron,reject
    iron=snrt_fe_band_enabled();nb=snrt_grain_band_bins();last_mode=4
    allocate(grains(cells,nb))
    if(iron)then
       kabs=fe_six_band_abs;ksca=fe_six_band_transport;last_mode=9
    else
       kabs=d03_band_abs;ksca=d03_band_transport
    endif
    call snrt_secondary_tables_load_from_environment(status)
    if(status/=0)stop 400
    call snrt_backend_initialize(status)
    if(status/=0)stop 401
    direction(:,1)=[1.0_c_float,0.0_c_float,0.0_c_float];direction(:,2)=-direction(:,1)
    neighbors(:,1)=1;neighbors(:,2)=2;weights=[.3d0,.7d0]
    reference=.5d0*(d03_edges(:bands)+d03_edges(2:));maximum_n=0;maximum_e=0
    do mode=0,last_mode
       tau=0;stau=0;dtau=0;columns=0;xi=.1d0;budget=1000
       grains=0
       if(mode>0)then
          grains(1,1:4)=[1d-5,2d-5,3d-5,4d-5]
          if(iron)grains(1,5:6)=[5d-5,6d-5]
          grains(2,:)=grains(1,:)*2
       endif
       if(iron.and.(mode==2.or.mode==6))grains(:,5:6)=0
       if(mode==5)grains(:,1:4)=0
       if(mode==2)then
          columns=1d18;budget=.05_c_float
       endif
       if(mode==4)grains=1d5
       q=0;shift=0
       do g=1,bands
          do c=1,cells
             if(iron.and.mode/=0.and.mode/=2.and.mode/=6.and.g>1)cycle
             q(c,1,g)=4
             mean=d03_edges(g)
             if(mode==2)mean=.3d0*d03_edges(g)+.7d0*d03_edges(g+1)
             if(mode==9)mean=.1d0
             shift(c,1,g)=(mean-reference(g))*real(q(c,1,g),c_double)
          enddo
       enddo
       if(mode==7.or.mode==8)then
          q=0;shift=0;q(:,1,2)=4
          mean=2.3664319132398464d0
          if(mode==8)mean=d03_edges(3)
          shift(:,1,2)=(mean-reference(2))*4
       endif
       saved_q=q;saved_shift=shift;before_n=sum(real(q,c_double));before_e=0
       do g=1,bands
          before_e=before_e+sum(reference(g)*real(q(:,:,g),c_double)+shift(:,:,g))
       enddo
       if(mode==3)grains(2,nb)=-1
       reject=mode==3.or.mode==7.or.mode==8
       he=-17;de=-17;deposition=-17;hh=-17;dd=-17
       rc=snrt_runtime_species_dust_step(q,direction,neighbors,tau,stau,dtau,budget,hh,dd,returned,raw,ag,at, &
            int(cells,c_int),int(cells,c_int),int(dirs,c_int),int(bands,c_int),0.0_c_float, &
            shift=shift,reference_ev=reference,hhe_energy=he,dust_energy=de,dust_energy_moment=em, &
            species_columns=columns,secondary_xi=xi,band_deposition=deposition, &
            grain_columns=grains,angular_weights=weights)
       if(reject)then
          if((mode==7.or.mode==8).and.rc/=10)stop 412
          if(mode==3.and.rc/=2)stop 413
          if(rc==0.or.any(q/=saved_q).or.any(shift/=saved_shift).or.any(he/=-17).or.any(de/=-17))stop 402
          if(any(hh/=-17).or.any(dd/=-17).or.any(deposition/=-17))stop 403
          cycle
       endif
       if(rc/=0)then
          write(*,*)'D03 spectrum mode/status',mode,rc
          stop 404
       endif
       after_n=sum(real(q,c_double))+sum(real(hh,c_double))+sum(real(dd,c_double))
       after_e=sum(he)+sum(de)
       do g=1,bands
          after_e=after_e+sum(reference(g)*real(q(:,:,g),c_double)+shift(:,:,g))
       enddo
       nerror=abs(after_n/before_n-1);eerror=abs(after_e/before_e-1)
       maximum_n=max(maximum_n,nerror);maximum_e=max(maximum_e,eerror)
       if(nerror>3d-7.or.eerror>2d-12.or.any(budget<0))stop 405
       if(mode==0.and.(any(q/=saved_q).or.any(shift/=saved_shift).or.any(de/=0)))stop 406
       if(mode==1.or.mode==5.or.mode==6)then
          do c=1,cells
             do g=1,bands
                if(saved_q(c,1,g)==0)cycle
                ta=sum(grains(c,:)*kabs(1,g,:));ts=sum(grains(c,:)*ksca(1,g,:))
                do d=1,dirs
                   expected=4*exp(-ta)*weights(d)*(1-exp(-ts))
                   if(d==1)expected=expected+4*exp(-ta-ts)
                   if(abs(real(q(c,d,g),c_double)-expected)>4d-7)stop 407
                enddo
                expected=4*d03_edges(g)*(1-exp(-ta))
                if(abs(de(c,g)-expected)>2d-12*max(expected,1d-30))stop 408
             enddo
          enddo
       endif
       if(mode==2)then
          if(sum(hh)<=0.or.sum(de)<=0.or.sum(deposition(:,4:8))<=0)stop 409
          if(any(sum(real(hh,c_double),dim=2)>.0500001d0))stop 410
       endif
       if(mode==4.and.(any(q/=0).or.any(shift/=0)))stop 411
    enddo
    write(*,'(A,2ES24.16)')'D03_ACTUAL_NODE_CONSERVATION max_N,max_E=',maximum_n,maximum_e
    write(*,*)'D03_SPECTRAL_ABSORPTION_SCATTER_FS2010_FINITE_ATOMS_ROLLBACK_PASS'
    if(iron)write(*,*)'FE_SPECTRAL_SIX_BIN_PURE_MIXED_ZERO_FE_HARD_TAIL_ROLLBACK_PASS'
  end subroutine

  subroutine check_primary_momentum()
    use iso_c_binding, only: c_int,c_float,c_double
    type(dust_fv_layout)::layout
    real(c_float)::q(2,2,1),direction(3,2),tau(2,1),stau(2,1,3),dtau(2,1),budget(2,3)
    real(c_float)::hhe(2,1,3),absorbed_dust(2,1),returned(2,1),raw(2,1),group(2,1),total(2)
    real(c_double)::moment(2,1,3),saved_moment(2,1,3),u(13),v(13),heat(2),fraction(2,1)
    real(c_double)::expected,gas_rho(0:2),velocity(3,0:2),thermal,kinetic,pressure,sound
    integer(c_int)::neighbor(6,2),rc
    integer::status,cell
    q(:,1,1)=4;q(:,2,1)=[0.0_c_float,4.0_c_float]
    direction(:,1)=[1.0_c_float,0.0_c_float,0.0_c_float];direction(:,2)=-direction(:,1)
    neighbor(:,1)=1;neighbor(:,2)=2
    tau=.5_c_float;dtau=tau;stau=0;budget=0;moment=-42
    rc=snrt_runtime_species_dust_step(q,direction,neighbor,tau,stau,dtau,budget,hhe,absorbed_dust,returned,raw, &
         group,total,2_c_int,2_c_int,2_c_int,1_c_int,0.0_c_float,moment)
    if(rc/=0.or.abs(moment(1,1,1)-absorbed_dust(1,1))>1d-12.or.any(moment(2,1,:)/=0))stop 110
    call dust_fv_initialize(2,[integer::],[real(c_double)::],layout,status)
    if(status/=0)stop 111
    u=0;u(1)=1.03d0;u(5)=1.5d0;u(6)=.01d0;u(10)=.02d0;fraction(:,1)=[.25d0,.75d0]
    do cell=1,2
       v=u;heat=-42
       call dust_fv_absorption_kick(layout,u,5d0/3,real(absorbed_dust(cell,:),c_double), &
            transpose(moment(cell,:,:)),[1d-12],fraction,1d6,1d-3,1d-13,v,heat,status)
       expected=moment(cell,1,1)*1d-6/(2.99792458d10*1d-13)
       if(status/=0.or.abs(v(2)-expected)>1d-16)stop 112
       if(abs(v(5)-u(5)+sum(heat)-absorbed_dust(cell,1)*1d-3)>1d-14)stop 113
       call dust_fv_decode(layout,v,5d0/3,gas_rho,velocity,thermal,kinetic,pressure,sound,status)
       if(status/=0.or.maxval(abs(velocity(:,0)))>1d-15.or.abs(thermal-1.5d0)>1d-14)stop 114
    enddo
    ! Fortran/C boundary rejects a bad shape without publishing the output.
    saved_moment=moment
    rc=snrt_runtime_species_dust_step(q,direction,neighbor,tau,stau,dtau,budget,hhe,absorbed_dust,returned,raw, &
         group,total,2_c_int,2_c_int,2_c_int,1_c_int,0.0_c_float,moment(:,:,1:2))
    if(rc==0.or.any(moment/=saved_moment))stop 115
    write(*,*)'SNRT_PRIMARY_BACKEND_TO_DUST_FV_MOMENTUM_WORK_HEAT_PASS'
  end subroutine

  subroutine check_pah_mixed()
    use dust_pah_mixed
    use dust_pah_live_model
    use dust_mass_physics, only: dust_pah_nbin,dust_mass_enabled,dust_iron_model
    use snrt_dust_contract
    integer,parameter::n=2
    type(dust_ir_table)::tab
    type(dust_ir_diagnostics)::diag
    real(dust_dp)::field(fe_nir,2,n),pop(dust_pah_nbin,n),old_pop(dust_pah_nbin,n)
    real(dust_dp)::bins(6,n),primary(9,n),ed(n),old_ed(n),td(n),ph(fe_nir,n),gas(n),cv(n),coupling(n),q(n)
    real(dust_dp)::saved_field(fe_nir,2,n),saved_ph(fe_nir,n),saved_ed(n),saved_td(n),hi,residual
    real(dust_dp)::ghost(fe_nir,2,0),rays(3,2),weights(2)
    integer::neighbors(6,n),remote(6,n),status,i
    logical::blocked(6,n)
    character(len=2048)::path
    call get_environment_variable('SNRT_TEST_PAH_CHARGED',path)
    if(trim(path)=='1')then
       call check_pah_charged_mixed()
       return
    endif
    call get_environment_variable('SNRT_PAH_NEUTRAL_TABLE',path)
    if(len_trim(path)==0.or..not.snrt_dust_contract_loaded)return
    dust_mass_enabled=.true.;dust_iron_model='fe_electric_compare_v1'
    call pah_live_prepare(status)
    if(status/=0)stop 301
    call snrt_dust_ir_initialize(tab,snrt_dust_contract_ir_energy_ev(1:fe_nir), &
         snrt_dust_contract_ir_weight_ev(1:fe_nir),snrt_dust_contract_ir_absorption_per_h_cm2(1:fe_nir), &
         snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature),10d0,status)
    if(status/=0)stop 302
    bins(:,1)=1d-28;bins(:,2)=0
    ed=0;td=20
    call iron_mixture_enthalpy(20d0,[2d-28,2d-28,2d-28],ed(1),hi,status)
    if(status/=0)stop 303
    old_ed=ed;field=0;ph=0;pop=0;pop(1,:)=1d-8;old_pop=pop
    primary=0;primary(2,:)=1d-22;gas=1d-12;cv=1d-14;coupling=0
    rays=0;rays(1,:)=[1d0,-1d0];weights=.5d0
    neighbors(:,1)=2;neighbors(:,2)=1;remote=0;blocked=.false.
    call pah_mixed_advance(tab,rays,weights,neighbors,1d12,1d0,3d8,bins,primary,field,pop,ed,td,ph,diag,status, &
         ghost,remote,blocked,gas,cv,coupling,q)
    if(status/=0.or.any(pop<0).or.any(abs(sum(pop,dim=1)-1d-8)>1d-20))stop 304
    residual=sum(ed-old_ed)+sum(field)*.5d0-sum(primary)
    do i=1,n
       residual=residual+dot_product(pah_level,pop(:,i)-old_pop(:,i))
    enddo
    if(abs(residual)>2d-9*sum(primary).or.any(gas/=1d-12))stop 305
    if(sum(pop(2:,:))<=0.or.sum(ph)<=0.or.ed(2)/=0)stop 306
    write(*,'(A,ES16.8)')'PAH_MIXED_AND_PURE_PRIMARY_IR_ENERGY_RESIDUAL=',abs(residual)/sum(primary)
    saved_field=field;saved_ph=ph;saved_ed=ed;saved_td=td;old_pop=pop
    primary=0;primary(9,:)=1d-22
    call pah_mixed_advance(tab,rays,weights,neighbors,1d12,1d0,3d8,bins,primary,field,pop,ed,td,ph,diag,status, &
         ghost,remote,blocked,gas,cv,coupling,q)
    if(status==0.or.any(field/=saved_field).or.any(ph/=saved_ph).or.any(ed/=saved_ed).or. &
         any(td/=saved_td).or.any(pop/=old_pop))stop 307
    write(*,*)'PAH_MIXED_PRIMARY_PARTITION_PURE_PAH_AND_HARD_PHOTON_ROLLBACK_PASS'
    dust_iron_model='none';dust_mass_enabled=.false.
  end subroutine

  subroutine check_six_component_ir()
    ! A physical electric/eddy BASE, not completed magnetic Fe or live AMR.
    integer,parameter::ncell=4
    type(dust_ir_table)::six_table
    type(dust_ir_diagnostics)::diag
    real(dust_dp)::pa(fe_ng,6),ps(fe_ng,6),psg(fe_ng,6),ia(fe_nir,6),isc(fe_nir,6),isg(fe_nir,6)
    real(dust_dp)::fraction(6,ncell),density(ncell),bins(6),m(3),ed(ncell),td(ncell),heat(ncell)
    real(dust_dp)::material(fe_nt,ncell),hi,opacity(fe_nir,ncell),field(fe_nir,2,ncell),photons(fe_nir,ncell)
    real(dust_dp)::old_field(fe_nir,2,ncell),old_ed(ncell),old_td(ncell),old_photons(fe_nir,ncell)
    real(dust_dp)::bad_edges(size(fe_edges)),balance,delta_scatter,dt
    real(dust_dp)::gas(ncell),gas_before(ncell),gas_capacity(ncell),coupling(ncell),transfer(ncell)
    integer::neighbors(6,ncell),i,j,status
    if(fe_full_optics_admitted)stop 220
    if(.not.fe_base_binding(fe_edges,fe_primary_ev,fe_ir_ev,fe_radius_cm,fe_density(1)))stop 221
    bad_edges=fe_edges;bad_edges(9)=bad_edges(9)+1
    if(fe_base_binding(bad_edges,fe_primary_ev,fe_ir_ev,fe_radius_cm,fe_density(1)))stop 222
    call fe_six_opacity_basis(fe_test_reference_mass,pa,ps,psg,ia,isc,isg,status)
    if(status/=0.or.any(pa<=0).or.any(ia<=0).or.any(ps<0).or.any(abs(psg)>ps))stop 223
    if(any(psg(9,5:6)/ps(9,5:6)<.99d0))stop 224
    fraction(:,1)=[.5d0,.5d0,0d0,0d0,0d0,0d0]
    fraction(:,2)=[0d0,0d0,.5d0,.5d0,0d0,0d0]
    fraction(:,3)=[0d0,0d0,0d0,0d0,.5d0,.5d0]
    fraction(:,4)=[.1d0,.1d0,.15d0,.15d0,.25d0,.25d0]
    density=[1d-26,2d-26,1d-26,1d-26]/fe_test_reference_mass
    do i=1,ncell
       bins=fraction(:,i)*fe_test_reference_mass
       m=[sum(bins(1:2)),sum(bins(3:4)),sum(bins(5:6))]
       do j=1,fe_nt
          call iron_mixture_enthalpy(fe_temperature(j),m,material(j,i),hi,status)
          if(status/=0)stop 225
       enddo
       call iron_mixture_enthalpy(40d0,m,ed(i),hi,status)
       if(status/=0)stop 226
       ed(i)=ed(i)*density(i)
       ! 2.366eV optical photons, below Fe's photoelectric threshold: this
       ! exercise does not assign an absorbed keV photoelectron to grain heat.
       heat(i)=density(i)*dot_product(pa(2,:),fraction(:,i))*1d9*fe_primary_ev(2)*1.602176634d-12
    enddo
    call snrt_dust_ir_initialize(six_table,fe_ir_ev,fe_ir_weight_ev,sum(ia,dim=2)/6,fe_temperature,10d0, &
         status,sum(material,dim=2)/ncell,ia)
    if(status/=0)stop 227
    neighbors=0
    do i=1,ncell
       if(i>1)neighbors(1,i)=i-1
       if(i<ncell)neighbors(2,i)=i+1
    enddo
    field=0;photons=0;td=40;old_ed=ed;dt=1d8
    diag=dust_ir_diagnostics()
    call snrt_dust_ir_advance(six_table,rays,weights,neighbors,1d17,dt,1d8,density,heat,field,td, &
         photons,diag,status,1d-9,128,ed,spread(1d0,1,ncell),material_dispatch=fe_test_dispatch, &
         transport_dispatch=snrt_runtime_ir_transport,absorb_dispatch=snrt_runtime_ir_absorb, &
         cell_material_u=material,cell_weights=fraction,thin_reabsorption=.true.)
    if(status/=0.or.any(ed<0).or.any(field<0).or.any(photons<0))then
       write(*,*)'FE_SIX_COMPONENT_IR_REJECT status=',status
       stop 228
    endif
    balance=(sum(ed-old_ed)+sum(field(:,1,:)*weights(1))+sum(field(:,2,:)*weights(2)) &
         +(diag%escaped_erg-diag%primary_erg)/1d51)/max(sum(old_ed),sum(heat)*dt)
    if(abs(balance)>1d-9)stop 229
    old_field=field;opacity=matmul(isc-isg,fraction)
    call snrt_runtime_ir_scatter(field,weights,density,opacity,dt,status)
    if(status/=0)stop 230
    delta_scatter=maxval(abs(sum(field,dim=2)-sum(old_field,dim=2)))/maxval(old_field)
    if(delta_scatter>1d-13)stop 231
    write(*,'(A,6ES16.8)')'FE_SIX_COMPONENT_IR T(4),energy,scatter=',td,balance,delta_scatter
    old_field=field;old_ed=ed;old_td=td;old_photons=photons
    ! The old four-component material backend must reject six components.
    ! Generalizing the Fortran bank must not silently send six into its ABI.
    call snrt_dust_ir_advance(six_table,rays,weights,neighbors,1d17,dt,1d8,density,heat,field,td, &
         photons,diag,status,1d-9,128,ed,spread(1d0,1,ncell),material_dispatch=snrt_runtime_dust_material, &
         cell_material_u=material,cell_weights=fraction)
    if(status==0.or.any(field/=old_field).or.any(ed/=old_ed).or.any(td/=old_td))stop 232
    if(any(photons/=old_photons))stop 233
    gas=1d-14;gas_capacity=1d-16;coupling=1d-28;transfer=0;gas_before=gas
    call snrt_dust_ir_advance(six_table,rays,weights,neighbors,1d17,dt,1d8,density,heat,field,td, &
         photons,diag,status,1d-9,128,ed,spread(1d0,1,ncell),material_dispatch=fe_test_dispatch, &
         transport_dispatch=snrt_runtime_ir_transport,absorb_dispatch=snrt_runtime_ir_absorb, &
         cell_material_u=material,cell_weights=fraction,thin_reabsorption=.true., &
         gas_energy=gas,gas_capacity=gas_capacity,conductance=coupling,gas_transfer=transfer)
    if(status/=0.or.any(gas<0).or.any(transfer<=0))stop 234
    ! Scale by exchanged/radiated energy, not by the dominant gas reservoir.
    balance=(sum(ed-old_ed)+sum(gas-gas_before) &
         +sum((field(:,1,:)-old_field(:,1,:))*weights(1)) &
         +sum((field(:,2,:)-old_field(:,2,:))*weights(2)) &
         +(diag%escaped_erg-diag%primary_erg)/1d51)/max(sum(old_ed),sum(transfer),sum(heat)*dt)
    if(abs(balance)>1d-9.or.diag%balance_relative>1d-9)stop 235
    write(*,'(A,ES16.8)')'FE_SIX_COMPONENT_IR_GAS_COUPLED_BALANCE=',balance
    write(*,*)'FE_REAL_ELECTRIC_BASE_SIX_COMPONENT_IR_TRANSPORT_SCATTER_ROLLBACK_PASS'
  end subroutine

  subroutine fe_test_dispatch(heating,density,old_energy,capacity,log_t,power,band, &
       material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
       gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
    real(dust_dp),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
    real(dust_dp),intent(in)::material_u(:),dt,background,bath,tolerance
    logical,intent(in)::use_u
    real(dust_dp),intent(out)::rate(:,:),temperature(:),next_energy(:)
    integer,intent(out)::ierr
    real(dust_dp),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(dust_dp),optional,intent(out)::gas_transfer(:)
    real(dust_dp),optional,intent(in)::cell_material_u(:,:),cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
    real(dust_dp)::bins(6,size(density)),f(4,size(density))
    ierr=dust_err_config;rate=0;temperature=bath;next_energy=old_energy
    if(present(gas_transfer))gas_transfer=0
    if(.not.use_u.or..not.present(cell_weights).or..not.present(basis_band))return
    if(size(log_t)/=fe_nt.or.size(cell_weights,1)/=6)return
    if(any(log_t/=log(fe_temperature)))return
    bins=cell_weights*spread(density*fe_test_reference_mass,1,6);f=0
    call iron_radiative_batch(fe_temperature,bath,bins,old_energy,heating,dt, &
         fe_test_reference_mass,basis_band,2,next_energy,temperature,f,rate,ierr, &
         gas_energy,gas_capacity,conductance,gas_transfer)
  end subroutine

  subroutine check_iron_radiation()
    ! Synthetic spectrum tests the native closure, NOT admitted Fe opacity.
    real(dust_dp)::nodes(8),basis(2,8,6),m(3),bins(6),u0,u1,td,f(4),r(2),q,h,dt,lo,hi,target
    real(dust_dp)::tt(3),res,eg,cv,kdt
    real(dust_dp)::batch_bins(6,4),batch_old(4),batch_heat(4),batch_e(4),batch_t(4),batch_f(4,4),batch_r(2,4)
    real(dust_dp)::saved_e(4),saved_t(4),saved_f(4,4),saved_r(2,4),batch_q(4),saved_q(4)
    integer::i,j,k,status
    nodes=[5d0,298.15d0,1000d0,1184d0,1665d0,1809d0,2000d0,3000d0]
    tt=[1184d0,1665d0,1809d0]
    do i=1,8
       do k=1,6
          basis(:,i,k)=[.3d0,.7d0]*1d6*log(nodes(i)/5d0)
       enddo
    enddo
    bins=[.1d0,.1d0,.15d0,.15d0,.25d0,.25d0];m=[.2d0,.3d0,.5d0]
    do i=1,2
       td=300d0;dt=1d3
       if(i==2)then
          td=3000d0;dt=1d4
       endif
       call iron_mixture_enthalpy(td,m,u0,hi,status)
       if(status/=0)stop 210
       do j=1,3
          call iron_mixture_enthalpy(tt(j),m,lo,hi,status)
          if(status/=0)stop 211
          target=lo+.37d0*(hi-lo)
          h=(target-u0)/dt+1d6*log(tt(j)/5d0)
          call iron_radiative_cell(nodes,5d0,bins,u0,h,dt,1d0,basis,0d0,1d0,0d0,u1,td,f,r,q,status)
          if(status/=0.or.td/=tt(j).or.abs(f(j+1)-.37d0)>1d-10)then
             write(*,*)'FE_IR_PHASE_REJECT direction/transition/status/T/f=',i,j,status,td,f
             stop 212
          endif
          res=(u1-u0+dt*(sum(r)-h)-q)/max(u0,u1,dt*h)
          if(abs(res)>2d-12.or.q/=0.or.any(r<0))stop 213
       enddo
    enddo
    ! Non-plateau equilibrium heating and stiff gas reservoir at trace dust.
    call iron_mixture_enthalpy(300d0,m,u0,hi,status)
    if(status/=0)stop 236
    batch_bins=spread(bins,2,4);batch_bins(:,4)=0
    batch_old=u0;batch_old(4)=0;batch_heat=0;dt=1d3
    do j=1,3
       call iron_mixture_enthalpy(tt(j),m,lo,hi,status)
       if(status/=0)stop 237
       batch_heat(j)=(lo+.37d0*(hi-lo)-u0)/dt+1d6*log(tt(j)/5d0)
    enddo
    batch_e=-1;batch_t=-1;batch_f=-1;batch_r=-1;batch_q=-1
    call iron_radiative_batch(nodes,5d0,batch_bins,batch_old,batch_heat,dt,1d0,basis,1, &
         batch_e,batch_t,batch_f,batch_r,status,spread(0d0,1,4),spread(1d0,1,4),spread(0d0,1,4),batch_q)
    if(status/=0.or.any(batch_t(:3)/=tt).or.batch_e(4)/=0.or.any(batch_r(:,4)/=0))stop 238
    do j=1,3
       if(abs(batch_f(j+1,j)-.37d0)>1d-10)stop 239
    enddo
    saved_e=batch_e;saved_t=batch_t;saved_f=batch_f;saved_r=batch_r;saved_q=batch_q
    call iron_radiative_batch(nodes,5d0,batch_bins,batch_old,batch_heat,dt,1d0,basis,2, &
         batch_e,batch_t,batch_f,batch_r,status,spread(0d0,1,4),spread(1d0,1,4),spread(0d0,1,4),batch_q)
    if(status/=0.or.any(batch_e/=saved_e).or.any(batch_t/=saved_t).or.any(batch_f/=saved_f))stop 240
    if(any(batch_r/=saved_r).or.any(batch_q/=saved_q))stop 241
    batch_heat(3)=1d99 ! Valid earlier cells must not be committed on this failure.
    call iron_radiative_batch(nodes,5d0,batch_bins,batch_old,batch_heat,dt,1d0,basis,2, &
         batch_e,batch_t,batch_f,batch_r,status,spread(0d0,1,4),spread(1d0,1,4),spread(0d0,1,4),batch_q)
    if(status==0.or.any(batch_e/=saved_e).or.any(batch_t/=saved_t).or.any(batch_f/=saved_f))stop 242
    if(any(batch_r/=saved_r).or.any(batch_q/=saved_q))stop 243
    write(*,*)'FE_OPENMP_BATCH_PHASE_EMPTY_CELL_BITWISE_PARITY_ATOMIC_ROLLBACK_PASS'
    h=1d6*log(300d0/5d0)
    call iron_radiative_cell(nodes,5d0,bins,u0,h,1d3,1d0,basis,0d0,1d0,0d0,u1,td,f,r,q,status)
    if(status/=0.or.abs(td-300d0)>1d-7.or.abs(u1/u0-1)>1d-10)stop 214
    bins=bins*1d-27;m=m*1d-27;u0=u0*1d-27
    eg=1d-12;cv=1d-16;kdt=1d-25;dt=1d15;h=1d-24
    call iron_radiative_cell(nodes,5d0,bins,u0,h,dt,1d0,basis,eg,cv,kdt,u1,td,f,r,q,status)
    if(status/=0.or.eg-q<0.or.u1<0.or.any(r<0))stop 215
    res=(u1-u0+dt*(sum(r)-h)-q)/max(u0,u1,dt*h,dt*sum(r),abs(q))
    if(abs(res)>2d-12)stop 216
    write(*,'(A,3ES16.8)')'FE_IR_STIFF_TRACE T,gas_transfer,closure=',td,q,res
    call iron_radiative_cell(nodes,5d0,bins,u0,1d5,dt,1d0,basis,eg,cv,kdt,u1,td,f,r,q,status)
    if(status==0.or.u1/=u0.or.any(r/=0).or.q/=0.or.any(f/=0))stop 217
    basis(1,4,6)=-1
    call iron_radiative_cell(nodes,5d0,bins,u0,h,dt,1d0,basis,eg,cv,kdt,u1,td,f,r,q,status)
    if(status==0.or.u1/=u0.or.any(r/=0).or.q/=0)stop 218
    write(*,*)'FE_RADIATIVE_PHASE_HEATING_COOLING_GAS_ENERGY_ROLLBACK_PASS'
  end subroutine

! Exact rejected live cell: no RAMSES execution, files or synthetic tolerances.
subroutine check_fe_cold_live_cell
  use, intrinsic :: iso_fortran_env, only: real64
  use dust_iron_radiation
  use dust_iron_optics
  implicit none
  real(real64)::pa(fe_ng,6),ps(fe_ng,6),pg(fe_ng,6),ia(fe_nir,6),isc(fe_nir,6),ig(fe_nir,6)
  real(real64)::basis(fe_nir,fe_nt,6),bins(6),ed,td,phase(4),rate(fe_nir),q,x,occupation,factor,residual
  real(real64),parameter::ref=1.398d-26,old=1.1753624451508176d-24,heat=2.0085971718334332d-40
  real(real64),parameter::dt=4.5617174005087106d14,eg=1.4138734008650590d-16
  real(real64),parameter::cv=1.6622049915806113d-19,k=4.2679263680107152d-38
  integer::i,j,status,n
  real(real64)::kn,cold_old,hi,cold_dt,expected,mass(3),cold_eg,cold_cv,cold_k
  bins=[0d0,1.0748067799354278d-35,4.1498580567538355d-30,4.1498591791261309d-30, &
       8.2997161135076692d-31,8.2997296440224048d-31]
  call fe_six_opacity_basis(ref,pa,ps,pg,ia,isc,ig,status)
  if(status/=0)stop 1
  do j=1,fe_nt
     do i=1,fe_nir
        x=fe_ir_ev(i)/(8.617333262145d-5*fe_temperature(j))
        if(x<1d-3)then
           occupation=1/x-.5d0+x/12-x**3/720
        else
           occupation=exp(-x)/(1-exp(-x))
        endif
        factor=8*acos(-1d0)*(fe_ir_ev(i)*1.602176634d-12)**3/(6.62607015d-27**3*2.99792458d10**2)*1.602176634d-12
        basis(i,j,:)=factor*ia(i,:)*fe_ir_weight_ev(i)*occupation
     enddo
  enddo
  ! Adjacent representable conductances cover compiler/roundoff variations
  ! of the real cell without loosening the final energy tolerance.
  kn=k
  do n=0,32
  call iron_radiative_cell(fe_temperature,10d0,bins,old,heat,dt,ref,basis,eg,cv,kn,ed,td,phase,rate,q,status)
  residual=ed-old+dt*(sum(rate)-heat)-q
  if(n==0)write(*,'(A,I4,4ES25.16)')'COLD_LIVE_REPRO status,T,Ed,Q,residual=',status,td,ed,q,residual
  if(status/=0)stop 2
  if(abs(residual)>2d-12*max(old,ed,dt*sum(rate),abs(q)))stop 3
  kn=nearest(kn,1d0)
  enddo
  write(*,*)'FE_COLD_LIVE_ROUNDOFF_NEIGHBOR_ENERGY_PASS'
  ! Actual Fe/D03 quadrature in a dark field. Both mixed Fe and Fe-free PAH
  ! hosts must cool below the former 5 K bracket without energy clipping.
  cold_dt=1d17
  do n=1,3
     if(n==3)bins(5:6)=0
     mass=[sum(bins(1:2)),sum(bins(3:4)),sum(bins(5:6))]
     call iron_mixture_enthalpy(4d0,mass,cold_old,hi,status)
     if(status/=0)stop 4
     cold_eg=0;cold_cv=1;cold_k=0
     if(n==2)then
        cold_cv=cold_old/4;cold_eg=cold_cv;cold_k=cold_cv/cold_dt
     endif
     call iron_radiative_cell(fe_temperature,5d0,bins,cold_old,0d0,cold_dt,ref,basis, &
          cold_eg,cold_cv,cold_k,ed,td,phase,rate,q,status,absolute_emission=.true.,photon_ev=fe_ir_ev)
     if(status/=0.or.td<=0.or.td>=4.or.ed<=0.or.cold_eg-q<0)stop 5
     residual=ed-cold_old+cold_dt*sum(rate)-q
     if(abs(residual)>2d-12*max(cold_old,ed,cold_dt*sum(rate),abs(q)))stop 6
     do i=1,fe_nir
        x=fe_ir_ev(i)/(8.617333262145d-5*td)
        occupation=0
        if(x<700)occupation=exp(-x)/(1-exp(-x))
        factor=8*acos(-1d0)*(fe_ir_ev(i)*1.602176634d-12)**3/ &
             (6.62607015d-27**3*2.99792458d10**2)*1.602176634d-12
        expected=factor*dot_product(ia(i,:),bins)/ref*fe_ir_weight_ev(i)*occupation
        if(abs(rate(i)-expected)>1d-12*max(sum(rate),tiny(1d0)))stop 7
     enddo
     write(*,'(A,I2,3ES20.10)')'COLD_ABSOLUTE T,E,residual: ',n,td,ed,residual
  enddo
  call iron_radiative_cell(fe_temperature,5d0,bins,0d0,0d0,cold_dt,ref,basis, &
       0d0,1d0,0d0,ed,td,phase,rate,q,status,absolute_emission=.true.,photon_ev=fe_ir_ev)
  if(status/=0.or.ed/=0.or.td/=0.or.any(rate/=0).or.q/=0)stop 8
  write(*,*)'FE_PAH_COLD_ABSOLUTE_PLANCK_GAS_ENERGY_ZERO_STATE_PASS'
  call check_moving_cold_ir(bins,basis,ia)
end subroutine

  subroutine check_moving_cold_ir(bins,basis,sigma)
    real(dust_dp),intent(in)::bins(6),basis(fe_nir,fe_nt,6),sigma(fe_nir,6)
    type(dust_ir_table)::cold_table
    type(dust_ir_diagnostics)::diag
    real(dust_dp)::u(fe_nt,1),mass(3),hi,rays(3,6),w(6),field(fe_nir,6,1),saved(fe_nir,6,1)
    real(dust_dp)::ph(fe_nir,1),ed(1),td(1),old_ed,old_t,density(1),cap(1),primary(1)
    real(dust_dp)::rho(6,1),p(3,6,1),p0(3,6,1),work(6,1),alpha(6,fe_nir,1),total,residual
    integer::links(6,1),j,b,s,mode,status
    logical::blocked(6,1)
    moving_cold_bins=bins;moving_cold_basis=basis
    mass=[sum(bins(1:2)),sum(bins(3:4)),sum(bins(5:6))]
    do j=1,fe_nt
       call iron_mixture_enthalpy(fe_temperature(j),mass,u(j,1),hi,status)
       if(status/=0)stop 260
    enddo
    call snrt_dust_ir_initialize(cold_table,fe_ir_ev,fe_ir_weight_ev, &
         matmul(sigma,bins)/fe_test_reference_mass,fe_temperature,10d0,status,u(:,1))
    if(status/=0)stop 261
    rays=0;w=1d0/6;links=0;blocked=.true.;density=1;cap=1;primary=0;rho(:,1)=bins
    do j=1,3
       rays(j,2*j-1)=1;rays(j,2*j)=-1
    enddo
    do b=1,6
       alpha(b,:,1)=sigma(:,b)*bins(b)/fe_test_reference_mass
    enddo
    ! Exercise both per-cell U(T) and table U(T) prechecks. The second
    ! long step starts below the 10 K bath and below the first 5 K knot.
    do mode=1,2
       call iron_mixture_enthalpy(20d0,mass,ed(1),hi,status)
       if(status/=0)stop 262
       td=20;field=0;ph=0;p=0
       p(1,:,1)=bins*2.99792458d10*1d-4
       do s=1,2
          old_ed=ed(1);old_t=td(1);saved=field;p0=p;work=0
          total=old_ed+sum(field)/6
          if(mode==1)then
             call snrt_dust_ir_advance(cold_table,rays,w,links,1d20,1d17,1d0,density,primary, &
                  field,td,ph,diag,status,1d-9,256,ed,cap,blocked_face=blocked, &
                  material_dispatch=moving_cold_material,cell_material_u=u,thin_reabsorption=.true., &
                  phase_density=rho,phase_momentum=p,phase_absorption=alpha,phase_work=work)
          else
             call snrt_dust_ir_advance(cold_table,rays,w,links,1d20,1d17,1d0,density,primary, &
                  field,td,ph,diag,status,1d-9,256,ed,cap,blocked_face=blocked, &
                  material_dispatch=moving_cold_material,thin_reabsorption=.true., &
                  phase_density=rho,phase_momentum=p,phase_absorption=alpha,phase_work=work)
          endif
          write(*,'(A,3I4,3ES20.10)')'MOVING_COLD_IR mode,step,status,T,E,work=',mode,s,status,td,ed,sum(work)
          if(status/=dust_ok.or.td(1)<=0.or.td(1)>=min(old_t,5d0).or.ed(1)<=0)stop 263
          residual=ed(1)-old_ed+sum(field-saved)/6+sum(work)
          if(abs(residual)>2d-9*total.or.diag%balance_relative>1d-9)stop 264
          if(any(field<0).or.any(ph<0))stop 265
       enddo
       ! Net-bath callers still reject this cold state without publication.
       saved=field;old_ed=ed(1);old_t=td(1)
       call snrt_dust_ir_advance(cold_table,rays,w,links,1d20,1d17,1d0,density,primary, &
            field,td,ph,diag,status,1d-9,256,ed,cap,blocked_face=blocked,material_dispatch=moving_cold_material)
       if(status==dust_ok.or.any(field/=saved).or.ed(1)/=old_ed.or.td(1)/=old_t)stop 266
       ! The moving branch retains the upper-energy guard and rollback.
       ed=2*u(fe_nt,1);old_ed=ed(1);p0=p;work=-123
       call snrt_dust_ir_advance(cold_table,rays,w,links,1d20,1d17,1d0,density,primary, &
            field,td,ph,diag,status,1d-9,256,ed,cap,blocked_face=blocked, &
            material_dispatch=moving_cold_material,phase_density=rho,phase_momentum=p, &
            phase_absorption=alpha,phase_work=work)
       if(status==dust_ok.or.any(field/=saved).or.ed(1)/=old_ed.or.any(p/=p0).or.any(work/=-123))stop 267
    enddo
    write(*,*)'MOVING_COLD_IR_ANALYTIC_PLANCK_LONG_STEP_ENERGY_ROLLBACK_PASS'
  end subroutine

  subroutine moving_cold_material(heating,density,old_energy,capacity,log_t,power,band, &
       material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
       gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
    real(dust_dp),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
    real(dust_dp),intent(in)::material_u(:),dt,background,bath,tolerance
    logical,intent(in)::use_u
    real(dust_dp),intent(out)::rate(:,:),temperature(:),next_energy(:)
    integer,intent(out)::ierr
    real(dust_dp),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(dust_dp),optional,intent(out)::gas_transfer(:)
    real(dust_dp),optional,intent(in)::cell_material_u(:,:),cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
    real(dust_dp)::phase(4),q
    call iron_radiative_cell(fe_temperature,fe_temperature(1),moving_cold_bins,old_energy(1),heating(1),dt, &
         fe_test_reference_mass,moving_cold_basis,0d0,1d0,0d0,next_energy(1),temperature(1),phase,rate(:,1),q,ierr, &
         absolute_emission=.true.,photon_ev=fe_ir_ev)
    if(present(gas_transfer))gas_transfer=q
  end subroutine

  subroutine check_radiative_sublimation()
    use dust_mass_physics
    use dust_composition_material
    real(dust_dp)::nodes(6),basis(2,6,4),curve(6),bins(4),next(4),ed,td,rate(2),q,phase,old
    real(dust_dp)::work_bins(4),work_e,gas,emitted,phases,step,remaining(3),hot_remaining(3)
    integer::k,b,status,r,nsub,local_steps
    dust_mass_enabled=.true.;dust_mass_model='carbon_olivine_2size_v1'
    dust_material_model='dl01_composition_v1';dust_sublimation='gd89_xu25_olivine_v1'
    dust_size_radius_cm=[1d-6,1d-5];dust_size_density=[2.2d0,3.8d0]
    nodes=[5d0,100d0,300d0,1000d0,2000d0,3000d0];bins=[0d0,0d0,.4d0,.6d0]
    call dust_composition_curve(nodes,[0d0,1d0],1d0,curve,status)
    if(status/=0)stop 151
    old=curve(3)
    do b=1,4
       do k=1,6
          basis(:,k,b)=[.4d0,.6d0]*1d-7*nodes(k)**4
       enddo
    enddo
    call dust_radiative_sublimation_cell(nodes,5d0,bins,old,1d8,1d3,1d0,basis, &
         1d12,1d9,0d0,next,ed,td,rate,q,phase,status)
    if(status/=0.or.sum(next)>=sum(bins).or.phase<=0.or.td<=300.or.any(next<0))stop 152
    if(abs(ed+phase-old+sum(rate)*1d3-q-1d11)>1d-10*1d11)stop 153
    write(*,'(A,4ES17.8)')'DUST_RADIATIVE_SUBLIMATION T/lost/phase/balance=',td,sum(bins-next),phase, &
         (ed+phase-old+sum(rate)*1d3-q-1d11)/1d11
    write(*,*)'DUST_RADIATIVE_SUBLIMATION_HEATING_MASS_PHASE_ENERGY_PASS'
    do r=1,3
       nsub=2**(r+2);step=1d3/nsub
       work_bins=bins;work_e=old;gas=1d12;emitted=0;phases=0
       do k=1,nsub
          call dust_radiative_sublimation_cell(nodes,5d0,work_bins,work_e,1d8,step,1d0,basis, &
               gas,1d9,0d0,next,ed,td,rate,q,phase,status)
          if(status/=0)stop 154
          work_bins=next;work_e=ed;gas=gas-q;emitted=emitted+sum(rate)*step;phases=phases+phase
       enddo
       if(abs(work_e+phases+emitted+(gas-1d12)-old-1d11)>1d-9*1d11)stop 155
       remaining(r)=sum(work_bins)
       write(*,'(A,I4,3ES17.8)')'DUST_RADIATIVE_SUBLIMATION_RESOLUTION n/mass/T/balance=',nsub,remaining(r),td, &
            (work_e+phases+emitted+(gas-1d12)-old-1d11)/1d11
    enddo
    if(abs(remaining(3)-remaining(2))>=abs(remaining(2)-remaining(1)))stop 156
    ! Cold, trace grains transferring a much larger gas energy over a long
    ! hydro step: the material residual must scale with the exchanged energy.
    bins=bins*8.3d-30;old=curve(1)*sum(bins)
    call dust_radiative_sublimation_cell(nodes,5d0,bins,old,2.76559726d-39,6.00154024d14,1d0,basis, &
         2.36850979d-16,1.66237703d-19,5.03754669d-38,next,ed,td,rate,q,phase,status)
    if(status/=0.or.q<=1d4*old.or.q>=2.36850979d-16)stop 157
    if(abs(ed+phase-old+(sum(rate)-2.76559726d-39)*6.00154024d14-q)>1d-10*q)stop 158
    write(*,'(A,3ES17.8)')'DUST_SUBLIMATION_GAS_DOMINATED T/Q_over_U/balance=',td,q/old, &
         (ed+phase-old+(sum(rate)-2.76559726d-39)*6.00154024d14-q)/q
    ! A single stiff endpoint loses the initial hot evaporation pulse. Resolve
    ! it locally while still spanning a cosmological-scale caller timestep.
    bins=[0d0,0d0,.4d0,.6d0];old=curve(6)
    do r=1,3
       call dust_radiative_sublimation_evolve(nodes,5d0,bins,old,0d0,1d15,1d0,basis, &
            1d12,1d9,0d0,next,ed,td,rate,q,phase,status,10d0**(-r-2),local_steps)
       if(status/=0)then
          write(*,*)'DUST_SUBLIMATION_TRANSIENT_FAILED',r,status
          stop 159
       endif
       if(sum(next)>=.999d0.or.any(next<0).or.phase<=0)stop 160
       if(abs(ed+phase-old+sum(rate)*1d15-q)>1d-8*old)stop 161
       hot_remaining(r)=sum(next)
       write(*,'(A,I4,I6,3ES17.8)')'DUST_SUBLIMATION_HOT_TRANSIENT r/steps/mass/T/balance=',r,local_steps, &
            hot_remaining(r),td,(ed+phase-old+sum(rate)*1d15-q)/old
    enddo
    if(abs(hot_remaining(3)-hot_remaining(2))>=abs(hot_remaining(2)-hot_remaining(1)))stop 162
    call dust_radiative_sublimation_evolve(nodes,5d0,bins,old,1d99,1d3,1d0,basis, &
         1d12,1d9,0d0,next,ed,td,rate,q,phase,status)
    if(status==0.or.any(next/=bins).or.ed/=old.or.any(rate/=0).or.q/=0.or.phase/=0)stop 163
    write(*,*)'DUST_SUBLIMATION_ADAPTIVE_FAILURE_ROLLBACK_PASS'
    dust_sublimation='none';dust_mass_model='bulk_v1';dust_material_model='fixed_mix';dust_mass_enabled=.false.
  end subroutine
  subroutine check_d03_optics()
    use dust_composition_optics
    use dust_composition_material, only: dust_composition_curve
    use snrt_dust_contract
    type(dust_ir_table)::local_table
    type(dust_ir_diagnostics)::local_diag,other_diag
    character(len=2048)::contract_path
    real(dust_dp)::bins(4),pa(d03_ng),ps(d03_ng),pg(d03_ng),depth(d03_ng)
    real(dust_dp)::ia(d03_nir),isc(d03_nir),ig(d03_nir),pure(d03_ng,4)
    real(dust_dp)::spectra(d03_nir,6),pa_changed(d03_ng)
    real(dust_dp)::u(4),tgrid(4),mass_ref,total_before,total_after,parity,dt_test
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
    do case_id=1,6
       dt_test=1d6
       if(case_id==6)then
          tgrid=[10d0,1000d0,2000d0,3000d0];dt_test=1d-6
       endif
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
       field=0;emitted=0;td=tgrid(2);ed=dust_density*u(2)
       other=field;other_emitted=emitted;other_td=td;other_ed=ed
       total_before=sum(ed)+dt_test*sum(heating)
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,dt_test,1d5,dust_density,heating, &
            field,td,emitted,local_diag,status,1d-10,128,ed,[1d0])
       if(status/=0)then
          write(*,*)'D03_REFERENCE_REJECTED',case_id,status
          stop 72
       endif
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,dt_test,1d5,dust_density,heating, &
            other,other_td,other_emitted,other_diag,status,1d-10,128,other_ed,[1d0], &
            material_dispatch=snrt_runtime_dust_material,transport_dispatch=snrt_runtime_ir_transport, &
            absorb_dispatch=snrt_runtime_ir_absorb)
       if(status/=0)stop 73
       total_after=sum(ed)+sum(field)*.5d0
       if(abs(total_after-total_before)>1d-10*total_before.or.local_diag%escaped_erg/=0)stop 74
       parity=max(maxval(abs(other-field))/maxval(field),maxval(abs(other_ed-ed))/maxval(ed), &
            maxval(abs(other_td-td))/maxval(td),maxval(abs(other_emitted-emitted))/maxval(emitted))
       if(parity>1d-10.or.other_diag%balance_relative>1d-10)stop 75
       if(case_id==6.and.any(td<300d0))stop 97
       spectra(:,case_id)=emitted(:,1)/sum(emitted)
       ! In the background-subtracted formulation, Td=Tbath, no primary and
       ! no excess IR is an exact stationary state, not spurious bath heating.
       field=0;td=10;ed=dust_density*u(1);emitted=0
       call snrt_dust_ir_advance(local_table,rays,weights,neighbors,1d12,dt_test,1d5,dust_density,[0d0], &
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

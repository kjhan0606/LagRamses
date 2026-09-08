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
  write(*,*)'DUST_BACKEND_PARITY_AND_ROLLBACK_PASS'
  call MPI_FINALIZE(info)
end program

! Small state test of the production hook, not a surrogate energy formula.
! The two levels mimic an in-flight coarse state and two fine substeps.
program dust_expansion_smoke
  use amr_parameters, only: dp,amr_block_size,twotondim,cosmo,hydro,pic,ncontrol,ngridmax
  use amr_commons, only: ncoarse,ncpu,active,reception,myid,nstep,nstep_coarse,nstep_coarse_old, &
       n_frw,aexp_frw,hexp_frw,tau_frw,t_frw,aexp,t,dtnew
  use hydro_parameters, only: nvar,idust_energy
  use hydro_commons, only: uold,unew
  use dust_mass_physics, only: dust_mass_enabled
  use dust_mass_runtime, only: dust_expansion_level
  use snrt_dust_live, only: snrt_dust_live_pack,snrt_dust_live_restore,snrt_dust_live_expand
  use snrt_dust_contract, only: snrt_dust_contract_load_from_environment,snrt_dust_contract_number_ir
  use snrt_state, only: snrt_state_restore_cell,snrt_checkpoint_cell_width,snrt_ndirection
  use dust_iron_radiation, only: iron_radiative_cell
  use dust_iron_material, only: iron_mixture_enthalpy
  use snrt_dust_ir, only: dust_ir_table,dust_ir_diagnostics,snrt_dust_ir_initialize,snrt_dust_ir_advance
  use mpi_mod
  implicit none
  real(dp),allocatable::reference(:,:),expected(:,:)
  real(dp),allocatable::ir(:),read_ir(:)
  real(dp)::primary_state(snrt_checkpoint_cell_width)
  integer::level,grid,child,cell,info
  real(dp)::test_bath
  logical::reject_bath=.false.
  call MPI_INIT(info)
  if(nvar<7)stop 1
  ncoarse=1; ncpu=1; amr_block_size=2; idust_energy=7
  allocate(active(2),reception(1,2))
  allocate(uold(1+6*twotondim,nvar),unew(1+6*twotondim,nvar))
  allocate(reference(size(uold,1),nvar),expected(size(uold,1),nvar))
  do cell=1,size(uold,1)
     reference(cell,:)=real(cell,dp)
  enddo
  uold=reference; unew=reference; expected=reference
  do level=1,2
     active(level)%ngrid=1; reception(1,level)%ngrid=1
     allocate(active(level)%igrid(1),reception(1,level)%igrid(1))
     active(level)%igrid=2*level-1
     reception(1,level)%igrid=2*level
  enddo
  dust_mass_enabled=.false.
  call dust_expansion_level(1,2d0)
  if(any(unew/=reference))stop 2
  dust_mass_enabled=.true.
  ! a: .01 -> .011 -> .012. The coarse in-flight energy sees both ratios.
  call dust_expansion_level(1,1.1d0)
  call dust_expansion_level(2,1.1d0)
  call dust_expansion_level(1,12d0/11d0)
  call dust_expansion_level(2,12d0/11d0)
  do grid=1,4
     do child=1,twotondim
        cell=ncoarse+((grid-1)/2)*(twotondim*2)+(child-1)*2+mod(grid-1,2)+1
        expected(cell,idust_energy)=reference(cell,idust_energy)*1.2d0**2
     enddo
  enddo
  if(maxval(abs(unew-expected)/reference)>1d-14)stop 3
  if(any(uold/=reference))stop 4
  ! Changing a by unity must not change any reservoir, including gas energy.
  expected=unew
  call dust_expansion_level(1,1d0)
  call dust_expansion_level(2,1d0)
  if(any(unew/=expected))stop 5
  ! The physical specific grain energy is unchanged after scale_v^2 correction.
  if(abs(unew(2,idust_energy)/reference(2,idust_energy)/1.2d0**2-1d0)>1d-14)stop 6
  ! Also exercise the actual update_time -> every in-flight level call site.
  ! This is an explicit linear clock fixture, not a fitted Friedmann history
  ! or a full AMR/reflux evolution. No main-loop/output branch is entered.
  unew=reference
  call snrt_dust_contract_load_from_environment(info)
  if(info/=0)stop 9
  call check_cmb_material()
  call check_cmb_transport()
  ngridmax=6
  primary_state=0;primary_state(1)=1
  call snrt_state_restore_cell(2,primary_state,info)
  if(info/=0)stop 10
  allocate(ir(snrt_dust_contract_number_ir*snrt_ndirection))
  allocate(read_ir(size(ir)))
  ir=1d-20
  call snrt_dust_live_restore(2,ir,info)
  if(info/=0)stop 11
  myid=1;cosmo=.true.;hydro=.true.;pic=.false.;ncontrol=100000
  nstep=1;nstep_coarse=0;nstep_coarse_old=0
  n_frw=2
  allocate(aexp_frw(0:2),hexp_frw(0:2),tau_frw(0:2),t_frw(0:2))
  tau_frw=[0d0,-1d0,-2d0];aexp_frw=[.012d0,.011d0,.01d0]
  hexp_frw=1d0;t_frw=[1d0,.5d0,0d0]
  t=-2d0;aexp=.01d0;dtnew(2)=1d0
  call update_time(2)
  call update_time(2)
  if(abs(aexp/.012d0-1)>1d-14)stop 7
  if(maxval(abs(unew-expected)/reference)>1d-14.or.any(uold/=reference))stop 8
  call snrt_dust_live_pack(2,read_ir,info)
  if(info/=0.or.maxval(abs(read_ir/ir*1.2d0**3-1d0))>1d-14)stop 12
  ! Invalid conversion is atomic; new state restored at the current epoch
  ! must not inherit a hidden pre-restart multiplier.
  call snrt_dust_live_expand(-1d0,info)
  if(info==0)stop 13
  call snrt_dust_live_pack(2,read_ir,info)
  if(info/=0.or.maxval(abs(read_ir/ir*1.2d0**3-1d0))>1d-14)stop 14
  call snrt_dust_live_restore(2,ir,info)
  if(info/=0)stop 15
  call snrt_dust_live_expand(1d0,info)
  call snrt_dust_live_pack(2,read_ir,info)
  if(info/=0.or.any(read_ir/=ir))stop 16
  write(*,'(A)')'DUST_IR_DILUTION PASS: global once/epoch, a^-3 composition, invalid rollback, restored epoch'
  write(*,'(A)')'DUST_EXPANSION_CLOCK_HOOK PASS: update_time -> both in-flight levels'
  write(*,'(A)')'DUST_EXPANSION PASS: active/reception, two levels, untouched gas/mass/uold/inactive cells'
  call MPI_FINALIZE(info)
contains
  subroutine check_cmb_transport()
    use snrt_dust_contract
    type(dust_ir_table)::table
    type(dust_ir_diagnostics)::diag,saved_diag
    real(dp),allocatable::opacity(:,:),curve(:,:),field(:,:,:),photons(:,:),saved_field(:,:,:),saved_photons(:,:)
    real(dp)::direction(3,1),weight(1),density(1),heat(1),td(1),ed(1),cv(1),gas(1),gas_cv(1),conduct(1),q(1)
    real(dp)::mixture(4,1),mass(3),old,upper,dx,dt,old_gas,residual,saved_td(1),saved_ed(1),saved_gas(1),saved_q(1)
    integer::ng,nt,j,status,neighbors(6,1)
    ng=snrt_dust_contract_number_ir;nt=snrt_dust_contract_number_temperature
    allocate(opacity(ng,4),curve(nt,1),field(ng,1,1),photons(ng,1))
    opacity=spread(snrt_dust_contract_ir_absorption_per_h_cm2(1:ng),2,4)
    mixture(:,1)=[.2d0,.2d0,.3d0,.3d0]
    mass=[.4d0,.6d0,0d0]*1d-26
    density=1d-26/snrt_dust_contract_mass_per_h_g
    do j=1,nt
       call iron_mixture_enthalpy(snrt_dust_contract_temperature_k(j),mass,curve(j,1),upper,status)
       if(status/=0)stop 30
    enddo
    curve=curve/density(1)
    call snrt_dust_ir_initialize(table,snrt_dust_contract_ir_energy_ev(1:ng), &
         snrt_dust_contract_ir_weight_ev(1:ng),opacity(:,1),snrt_dust_contract_temperature_k(1:nt), &
         10d0,status,curve(:,1),opacity)
    if(status/=0)stop 31
    test_bath=272.7d0
    call iron_mixture_enthalpy(20d0,mass,old,upper,status)
    ed=old;td=20;cv=1d-20;gas_cv=1d-20;gas=2d-20;old_gas=gas(1);conduct=1d-30;q=0
    field=0;photons=0;heat=0;direction=0;direction(1,1)=1;weight=1;neighbors=0
    dx=1d20;dt=1d8;diag=dust_ir_diagnostics()
    call snrt_dust_ir_advance(table,direction,weight,neighbors,dx,dt,2.99792458d8,density,heat, &
         field,td,photons,diag,status,1d-9,256,ed,cv,gas_energy=gas,gas_capacity=gas_cv, &
         conductance=conduct,gas_transfer=q,cell_material_u=curve,cell_weights=mixture,bath_dispatch=test_material)
    if(status/=0)then
       write(*,*)'CMB transport rejected:',status
       stop 32
    endif
    residual=(sum(field)+ed(1)-old+gas(1)-old_gas)*dx**3+diag%escaped_erg-diag%background_erg
    if(abs(residual)>1d-9*diag%background_erg.or.diag%background_erg<=0)stop 33
    if(any(field<0).or.any(photons<0).or.gas(1)<=old_gas)stop 34
    saved_field=field;saved_photons=photons;saved_td=td;saved_ed=ed;saved_gas=gas;saved_q=q;saved_diag=diag
    reject_bath=.true.
    call snrt_dust_ir_advance(table,direction,weight,neighbors,dx,dt,2.99792458d8,density,heat, &
         field,td,photons,diag,status,1d-9,256,ed,cv,gas_energy=gas,gas_capacity=gas_cv, &
         conductance=conduct,gas_transfer=q,cell_material_u=curve,cell_weights=mixture,bath_dispatch=test_material)
    reject_bath=.false.
    if(status==0)stop 35
    if(any(field/=saved_field).or.any(photons/=saved_photons).or.any(ed/=saved_ed).or.any(td/=saved_td))stop 36
    if(any(gas/=saved_gas).or.any(q/=saved_q).or.diag%background_erg/=saved_diag%background_erg)stop 37
    write(*,'(A,ES12.4)')'CMB_TRANSPORT PASS: explicit bath closure and rejected-call rollback; relative=', &
         abs(residual)/diag%background_erg
  end subroutine

  subroutine test_material(heating,density,old_energy,log_t,basis_band,cell_weights,dt, &
       gas_energy,gas_capacity,conductance,rate,temperature,next_energy,gas_transfer,background_transfer,ierr)
    use snrt_dust_contract
    real(dp),intent(in)::heating(:),density(:),old_energy(:),log_t(:),basis_band(:,:,:),cell_weights(:,:),dt
    real(dp),intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(dp),intent(out)::rate(:,:),temperature(:),next_energy(:),gas_transfer(:),background_transfer(:)
    integer,intent(out)::ierr
    real(dp)::bins(6),phase(4),bands(size(rate,1),size(log_t),6)
    integer::i
    ! Deliberately fill rejected trial outputs; none may leak into caller state.
    rate=0;temperature=0;next_energy=0;gas_transfer=0;background_transfer=0;ierr=2
    if(reject_bath)return
    bands=0;bands(:,:,1:4)=basis_band
    do i=1,size(density)
       bins=0;bins(1:4)=cell_weights(:,i)*density(i)*snrt_dust_contract_mass_per_h_g
       call iron_radiative_cell(exp(log_t),test_bath,bins,old_energy(i),heating(i),dt, &
            snrt_dust_contract_mass_per_h_g,bands,gas_energy(i),gas_capacity(i),conductance(i), &
            next_energy(i),temperature(i),phase,rate(:,i),gas_transfer(i),ierr, &
            photon_ev=snrt_dust_contract_ir_energy_ev(1:size(rate,1)),cmb_exchange=background_transfer(i))
       if(ierr/=0)return
    enddo
  end subroutine

  subroutine check_cmb_material()
    use snrt_dust_contract
    real(dp),allocatable::basis(:,:,:),rate(:)
    real(dp)::bins(6),mass(3),old,upper,ed,td,phase(4),q,bg,heat,cv,conductance,dt,residual,x,factor
    real(dp)::initial_t(4),bath_t(4),gas_t(4),steps(4)
    integer::nt,ng,i,j,test,status
    nt=snrt_dust_contract_number_temperature;ng=snrt_dust_contract_number_ir
    allocate(basis(ng,nt,6),rate(ng))
    do j=1,nt
       do i=1,ng
          x=snrt_dust_contract_ir_energy_ev(i)/(8.617333262145d-5*snrt_dust_contract_temperature_k(j))
          factor=8*acos(-1d0)*(snrt_dust_contract_ir_energy_ev(i)*1.602176634d-12)**3 / &
               (6.62607015d-27**3*2.99792458d10**2)*1.602176634d-12
          basis(i,j,:)=factor*snrt_dust_contract_ir_absorption_per_h_cm2(i)* &
               snrt_dust_contract_ir_weight_ev(i)*exp(-x)/(1-exp(-x))
       enddo
    enddo
    bins=[2d-27,2d-27,3d-27,3d-27,0d0,0d0]
    mass=[sum(bins(1:2)),sum(bins(3:4)),0d0]
    initial_t=[20d0,2.727d0,30d0,40d0];bath_t=[272.7d0,2.727d0,30d0,20d0]
    gas_t=[0d0,0d0,2d0,0d0];steps=[1d8,1d12,1d8,1d8]
    do test=1,4
       call iron_mixture_enthalpy(initial_t(test),mass,old,upper,status)
       if(status/=0)stop 20
       heat=0;cv=1d-20;conductance=0;dt=steps(test)
       if(test==3)conductance=1d-30
       call iron_radiative_cell(snrt_dust_contract_temperature_k(1:nt),bath_t(test),bins,old,heat,dt, &
            snrt_dust_contract_mass_per_h_g,basis,cv*gas_t(test),cv,conductance,ed,td,phase,rate,q,status, &
            photon_ev=snrt_dust_contract_ir_energy_ev(1:ng),cmb_exchange=bg)
       write(*,'(A,2I4,5ES22.12)')'CMB_MATERIAL case,status,T,Ed,Qgas,Qcmb,IR=', &
            test,status,td,ed,q,bg,dt*sum(rate)
       if(status/=0)stop 21
       residual=ed-old+dt*sum(rate)-q-bg
       if(abs(residual)>3d-12*max(old,ed,abs(q),bg,dt*sum(rate),tiny(1d0)))stop 22
       if(any(rate<0).or.bg<0)stop 23
       if(test==1.and.(ed<=old.or.bg<=0.or.td>bath_t(test)*(1+1d-10)))stop 24
       if(test==2.and.abs(td/bath_t(test)-1)>1d-10)stop 25
       if(test==3.and.(q>=0.or.bg<=0.or.td>=bath_t(test)))stop 26
       if(test==4.and.(sum(rate)<=0.or.bg/=0))stop 27
    enddo
    write(*,'(A)')'CMB_MATERIAL PASS: finite heating, 2.727K equilibrium, colder gas, warm excess; explicit bath ledger'
  end subroutine
end program

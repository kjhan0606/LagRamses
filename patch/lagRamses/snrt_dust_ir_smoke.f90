! Compact numeric-fixture driver for the actual native module, no RAMSES stubs.
program snrt_dust_ir_smoke
  use snrt_dust_ir
  use snrt_dust_coupling
  use snrt_thermochemistry, only: snrt_partition_absorption, snrt_thermochemistry_ok
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(dust_ir_table) :: table
  type(dust_ir_diagnostics) :: diag, before_diag
  real(dust_dp), allocatable :: freq(:), fw(:), sigma(:), nodes(:), direction(:,:), weight(:)
  real(dust_dp), allocatable :: density(:), primary(:), energy(:,:,:), photons(:,:), temperature(:)
  real(dust_dp), allocatable :: before(:,:,:), before_p(:,:), before_t(:), integral(:,:)
  integer, allocatable :: neighbor(:,:)
  real(dust_dp) :: dx,dt,c,cmb,tol,escaped,absorbed,injected,stored,balance
  integer :: ng,nt,nd,nc,steps,i,g,d,it,unit,ierr
  character(len=2048) :: filename
  call get_command_argument(1,filename)
  open(newunit=unit,file=trim(filename),status='old',action='read')
  read(unit,*) ng,nt,nd,nc,steps
  read(unit,*) dx,dt,c,cmb,tol
  allocate(freq(ng),fw(ng),sigma(ng),nodes(nt),direction(3,nd),weight(nd),neighbor(6,nc))
  allocate(density(nc),primary(nc),energy(ng,nd,nc),photons(ng,nc),temperature(nc),integral(ng,nc))
  read(unit,*) freq
  read(unit,*) fw
  read(unit,*) sigma
  read(unit,*) nodes
  read(unit,*) direction
  read(unit,*) weight
  read(unit,*) neighbor
  read(unit,*) density
  read(unit,*) primary
  close(unit)
  call snrt_dust_ir_initialize(table,freq,fw,sigma,nodes,cmb,ierr)
  if (ierr/=dust_ok) stop 1
  energy=0; photons=0; temperature=0; escaped=0; absorbed=0; injected=0
  do it=1,steps
     call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
          energy,temperature,photons,diag,ierr,tol,128)
     if (ierr/=dust_ok) then
        write(*,*) 'NATIVE_DUST_STEP_ERROR',ierr,it
        stop 2
     end if
     escaped=escaped+diag%escaped_erg; absorbed=absorbed+diag%absorbed_erg
     injected=injected+diag%primary_erg
  end do
  before=energy; before_p=photons; before_t=temperature; before_diag=diag
  ! Force a late nonlinear failure, and separately a CFL admission failure.
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,1)
  if (ierr/=dust_err_convergence) stop 3
  call assert_unchanged()
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,100*dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_err_cfl) stop 4
  call assert_unchanged()
  ! A reciprocal-link violation is caught without a partial state commit.
  i=neighbor(2,1)
  if (i<=0) stop 5
  neighbor(1,i)=0
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_err_state) stop 6
  call assert_unchanged()
  neighbor(1,i)=1
  ! Save physical result before zero/weak checks, which have their own state.
  integral=0
  do d=1,nd
     integral=integral+energy(:,d,:)*weight(d)
  end do
  stored=sum(integral)*dx**3
  balance=abs(stored+escaped-injected)/injected
  if (balance>tol) stop 7
  write(*,'(a)') 'NATIVE_DUST_IR_OK'
  write(*,'(es26.16e3)') stored,escaped,absorbed,balance
  write(*,'(es26.16e3)') temperature
  write(*,'(es26.16e3)') integral
  write(*,'(es26.16e3)') photons
  energy=0; photons=0; temperature=0; primary=0
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_ok .or. any(energy/=0) .or. any(photons/=0)) stop 8
  density=1; primary=1d-42
  call snrt_dust_ir_advance(table,direction,weight,neighbor,dx,dt,c,density,primary, &
       energy,temperature,photons,diag,ierr,tol,128)
  if (ierr/=dust_ok .or. maxval(energy)<=0 .or. maxval(photons)<=0) stop 9
  call coupling_checks()
  call transient_checks()
  call halo_checks()
contains
  subroutine halo_checks()
    type(dust_ir_table) :: halo_table
    type(dust_ir_diagnostics) :: whole,part(2)
    real(dust_dp) :: rays(3,2),field(2,2,2),old(2,2,2),split(2,2,2),p(2,2),t(2)
    integer :: links(6,2),remote(6,1),status,i,j
    rays(:,1)=[1d0,0d0,0d0]; rays(:,2)=[-1d0,0d0,0d0]
    links=0; links(1:2,1)=2; links(1:2,2)=1
    remote=0; remote(1:2,1)=1
    old(:,:,1)=.2d0; old(:,:,2)=.7d0; field=old; p=0; t=0
    call snrt_dust_ir_initialize(halo_table,[.001d0,.01d0],[.001d0,.01d0], &
         [1d-21,1d-21],[10d0,20d0],10d0,status)
    if(status/=dust_ok)stop 70
    call snrt_dust_ir_advance(halo_table,rays,[.5d0,.5d0],links,1d0,1d0,.1d0,[0d0,0d0],[0d0,0d0], &
         field,t,p,whole,status,1d-10,128)
    if(status/=dust_ok)stop 71
    split=old; p=0; t=0; links=0
    do i=1,2
       j=3-i
       call snrt_dust_ir_advance(halo_table,rays,[.5d0,.5d0],links(:,i:i),1d0,1d0,.1d0,[0d0],[0d0], &
            split(:,:,i:i),t(i:i),p(:,i:i),part(i),status,1d-10,128, &
            ghost_energy=old(:,:,j:j),ghost_index=remote)
       if(status/=dust_ok)stop 72
    end do
    if(maxval(abs(split-field))>1d-14.or.abs(sum(part%interface_erg))>1d-14)stop 73
    if(any(part%escaped_erg/=0).or.part(1)%interface_erg==0)stop 74
    old=-1
    call snrt_dust_ir_advance(halo_table,rays,[.5d0,.5d0],links(:,1:1),1d0,1d0,.1d0,[0d0],[0d0], &
         split(:,:,1:1),t(1:1),p(:,1:1),whole,status,1d-10,128,ghost_energy=old(:,:,1:1),ghost_index=remote)
    if(status/=dust_err_state.or.any(split/=field))stop 75
    write(*,'(a)')'NATIVE_DUST_IR_HALO_OK split_matches_whole=1 interface_cancels=1 rollback=1'
  end subroutine

  subroutine transient_checks()
    type(dust_ir_table) :: transient_table
    type(dust_ir_diagnostics) :: result
    real(dust_dp) :: rays(3,2), weights(2), field(2,2,1), emitted(2,1), grain_t(1)
    real(dust_dp) :: grain_e(1), capacity(1), previous_e(1), saved_field(2,2,1), saved_emitted(2,1)
    real(dust_dp) :: total,initial,stored_ir,saved_temperature(1)
    real(dust_dp) :: test_energy(2),expected_photons,prefactor,warm,cold
    integer :: links(6,1),status,j
    rays(:,1)=[1d0,0d0,0d0]; rays(:,2)=[-1d0,0d0,0d0]
    weights=.5d0; links=0; field=0; emitted=0
    capacity=1d-24; grain_t=20d0; grain_e=capacity*grain_t; initial=sum(grain_e)*1d36
    call snrt_dust_ir_initialize(transient_table,[.001d0,.01d0],[.001d0,.01d0], &
         [1d-21,1d-21],[10d0,20d0,50d0,100d0],10d0,status)
    if(status/=dust_ok)stop 51
    previous_e=grain_e
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_ok)then
       write(*,*)'TRANSIENT_ERROR',status
       stop 52
    end if
    stored_ir=sum(field)*.5d0*1d36
    total=stored_ir+sum(grain_e)*1d36+result%escaped_erg
    if(abs(total-initial)/initial>1d-10)stop 53
    if(any(grain_e>=previous_e).or.any(grain_t<10d0).or.sum(emitted)<=0)stop 54
    saved_field=field; saved_emitted=emitted; previous_e=grain_e; saved_temperature=grain_t
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,[0d0])
    if(status/=dust_err_state.or.any(field/=saved_field).or.any(grain_e/=previous_e).or. &
         any(emitted/=saved_emitted).or.any(grain_t/=saved_temperature))stop 55
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,dust_energy=grain_e)
    if(status/=dust_err_config.or.any(field/=saved_field).or.any(grain_e/=previous_e))stop 56
    initial=(sum(field)*.5d0+sum(grain_e))*1d36
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[1d-27], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_ok)stop 57
    total=(sum(field)*.5d0+sum(grain_e))*1d36+result%escaped_erg
    if(abs(total-initial-result%primary_erg)/max(initial,result%primary_erg)>1d-10)stop 58
    ! Large capacity approaches fixed temperature, not zero emission. The
    ! physical radiative loss can be much smaller than one material-energy ULP.
    field=0; emitted=0; capacity=1d20; grain_t=20d0; grain_e=capacity*grain_t
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d-6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_ok.or.sum(emitted)<=0.or.sum(field)<=0)stop 59
    test_energy=[.001d0,.01d0]
    do j=1,2
       prefactor=8*acos(-1d0)*(test_energy(j)*1.602176634d-12)**3 / &
            ((6.62607015d-27)**3*(2.99792458d10)**2)*1.602176634d-12*1d-21*test_energy(j)
       warm=1/(exp(test_energy(j)/(8.617333262145d-5*20d0))-1)
       cold=1/(exp(test_energy(j)/(8.617333262145d-5*10d0))-1)
       expected_photons=1d-6*prefactor*(warm-cold)/(test_energy(j)*1.602176634d-12)
       if(abs(emitted(j,1)/expected_photons-1)>1d-12)stop 60
    end do
    field=0; emitted=0; capacity=1d-24; grain_t=10
    grain_e=nearest(capacity*grain_t,-1d0); initial=sum(grain_e)*1d36
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_ok.or.abs(sum(grain_e)*1d36-initial)/initial>1d-13)stop 61
    ! An advected near-bath deficit consumes the declared energy error;
    ! it is neither free background heating nor an unreported clip.
    field=0; emitted=0; grain_t=10d0*(1-2d-12); grain_e=capacity*grain_t
    previous_e=grain_e; initial=sum(grain_e)*1d36
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_ok.or.result%balance_relative<=0.or.result%balance_relative>1d-10)stop 63
    if(abs(result%balance_relative-(sum(grain_e)-sum(previous_e))*1d36/initial)>1d-15)stop 64
    if(any(field/=0).or.any(emitted/=0))stop 65
    grain_t=10d0*(1-2d-9); grain_e=capacity*grain_t; previous_e=grain_e
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_err_state.or.any(grain_e/=previous_e))stop 66
    grain_t=9.99d0; grain_e=capacity*grain_t; previous_e=grain_e
    call snrt_dust_ir_advance(transient_table,rays,weights,links,1d12,1d6,1d5,[1d0],[0d0], &
         field,grain_t,emitted,result,status,1d-10,128,grain_e,capacity)
    if(status/=dust_err_state.or.any(grain_e/=previous_e))stop 62
    write(*,'(a)')'NATIVE_DUST_TRANSIENT_OK cooling=1 heating=1 material_plus_radiation_closure=1 rollback=1'
  end subroutine

  subroutine coupling_checks()
    real(dust_dp) :: sigma_test(3), tau_hhe(3), tau, tau_total
    real(dust_dp) :: available(3), assigned(3), assigned_dust, returned, unassigned
    real(dust_dp) :: reference_available(3), reference_assigned(3), reference_unassigned
    real(dust_dp) :: absorbed_test(2,2), mean_energy(2), heating(2), expected(2)
    real(dust_dp) :: tau_groups(3,3), dust_groups(3), raw_groups(3)
    real(dust_dp) :: initial_available(3), previous_available(3), total_hhe(3)
    real(dust_dp) :: raw, closure, dust_fraction, total_raw
    real(dust_dp) :: before_available(3), bad_heat(2), bad_absorbed(3,1)
    integer :: local_ierr, reference_ierr, group

    sigma_test = [1d-21, 2d-21, 3d-21]
    call snrt_dust_prepare_optical_depth(1d3, .5d0, 2d0, sigma_test(1), tau, local_ierr)
    if (local_ierr /= dust_coupling_ok .or. abs(tau-1d-18)>1d-30) stop 21
    call snrt_dust_total_optical_depth([.2d0,.1d0,.1d0], .2d0, tau_total, local_ierr)
    if (local_ierr /= dust_coupling_ok .or. abs(tau_total-.6d0)>1d-14) stop 22
    call snrt_dust_total_optical_depth([.2d0,.1d0], .2d0, tau_total, local_ierr)
    if (local_ierr /= dust_coupling_err_shape) stop 23

    ! With ample H/He inventory, the allocation follows the component optical
    ! depths and the direct dust fraction is tau_dust/tau_total.
    raw = .5d0
    tau_hhe = [.2d0,.1d0,.1d0]
    available = [10d0,10d0,10d0]
    call snrt_dust_partition_group(raw, tau_hhe, .2d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    closure = raw-sum(assigned)-assigned_dust-returned-unassigned
    if (local_ierr /= dust_coupling_ok .or. abs(closure)>1d-14 .or. &
         abs(assigned_dust-raw*.2d0/.6d0)>1d-14 .or. returned/=0d0 .or. &
         unassigned/=0d0) stop 24

    ! When H/He inventories saturate, only the dust optical-depth fraction of
    ! the H/He excess is transferred.  The rest is returned to photons.
    available = [.05d0,.05d0,.05d0]
    call snrt_dust_partition_group(raw, tau_hhe, .2d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    dust_fraction = 1d0-exp(-.2d0)
    if (local_ierr /= dust_coupling_ok .or. returned<=0d0 .or. &
         abs(assigned_dust-(raw*.2d0/.6d0 + &
         (raw*.4d0/.6d0-.15d0)*dust_fraction))>1d-13 .or. &
         abs(raw-sum(assigned)-assigned_dust-returned-unassigned)>1d-13) stop 25

    ! Several groups consume one shared H/He reservoir.  The adapter must
    ! conserve the component-wise reservoir decrement, not just each scalar
    ! call in isolation.
    tau_groups(:,1) = [.20d0,.10d0,.10d0]
    tau_groups(:,2) = [.10d0,.20d0,.05d0]
    tau_groups(:,3) = [.05d0,.15d0,.25d0]
    dust_groups = [.20d0,.10d0,.30d0]
    raw_groups = [.15d0,.11d0,.09d0]
    initial_available = [.08d0,.07d0,.06d0]
    available = initial_available
    previous_available = available
    total_hhe = 0d0
    total_raw = 0d0
    do group = 1, 3
       call snrt_dust_partition_group(raw_groups(group), tau_groups(:,group), &
            dust_groups(group), available, assigned, assigned_dust, returned, &
            unassigned, local_ierr)
       if (local_ierr /= dust_coupling_ok .or. &
            abs(raw_groups(group)-sum(assigned)-assigned_dust-returned- &
            unassigned)>1d-13 .or. any(available>previous_available+1d-14)) stop 26
       total_hhe = total_hhe + assigned
       total_raw = total_raw + raw_groups(group)
       previous_available = available
    end do
    if (any(abs(initial_available-available-total_hhe)>1d-13) .or. &
         total_raw <= 0d0) stop 27

    ! Zero dust is an exact regression of the existing H/He partition.
    raw = .12d0
    reference_available = [.1d0,.1d0,.1d0]
    available = reference_available
    call snrt_partition_absorption(raw, tau_hhe, reference_available, &
         reference_assigned, reference_ierr, reference_unassigned)
    call snrt_dust_partition_group(raw, tau_hhe, 0d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    if (reference_ierr /= snrt_thermochemistry_ok .or. &
         local_ierr /= dust_coupling_ok .or. &
         any(transfer(available,[0_int64],3) /= &
             transfer(reference_available,[0_int64],3)) .or. &
         any(transfer(assigned,[0_int64],3) /= &
             transfer(reference_assigned,[0_int64],3)) .or. &
         assigned_dust/=0d0 .or. returned/=0d0 .or. unassigned/=0d0) stop 26

    ! No finite H/He inventory is a returned photon, not a fake dust heat
    ! source, and failed calls do not mutate the caller's reservoir.
    available = 0d0
    before_available = available
    call snrt_dust_partition_group(.2d0, tau_hhe, 0d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    if (local_ierr /= dust_coupling_ok .or. returned/=.2d0 .or. &
         assigned_dust/=0d0 .or. unassigned/=0d0 .or. &
         any(transfer(available,[0_int64],3) /= transfer(before_available, &
         [0_int64],3))) stop 27
    available = [.1d0,.1d0,.1d0]
    before_available = available
    call snrt_dust_partition_group(.2d0, tau_hhe, -1d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    if (local_ierr /= dust_coupling_err_input .or. &
         any(transfer(available,[0_int64],3) /= transfer(before_available, &
         [0_int64],3))) stop 28
    call snrt_dust_partition_group(.2d0, [0d0,0d0,0d0], 0d0, available, assigned, &
         assigned_dust, returned, unassigned, local_ierr)
    if (local_ierr /= dust_coupling_err_input) stop 29
    call snrt_dust_prepare_optical_depth(huge(1d0), 1d0, huge(1d0), 1d0, &
         tau, local_ierr)
    if (local_ierr /= dust_coupling_err_input) stop 30

    ! Heating ledger: (photons cm^-3)*(eV)*EV_ERG/dt[s].
    absorbed_test = reshape([1d10,2d10,3d10,4d10], shape(absorbed_test))
    mean_energy = [13.6d0,100d0]
    call snrt_dust_heating_from_absorbed(absorbed_test, mean_energy, 10d0, &
         heating, local_ierr)
    expected = [(1d10*13.6d0+2d10*100d0)*1.602176634d-12/10d0, &
                (3d10*13.6d0+4d10*100d0)*1.602176634d-12/10d0]
    if (local_ierr /= dust_coupling_ok .or. any(abs(heating-expected)>1d-14*expected)) stop 30
    absorbed_test = 0d0
    call snrt_dust_heating_from_absorbed(absorbed_test, mean_energy, 10d0, &
         heating, local_ierr)
    if (local_ierr /= dust_coupling_ok .or. any(heating/=0d0)) stop 31
    absorbed_test = 1d-42
    call snrt_dust_heating_from_absorbed(absorbed_test, mean_energy, 10d0, &
         heating, local_ierr)
    if (local_ierr /= dust_coupling_ok .or. any(.not. (heating>0d0))) stop 32
    bad_absorbed = 0d0
    call snrt_dust_heating_from_absorbed(bad_absorbed, mean_energy, 1d0, &
         bad_heat, local_ierr)
    if (local_ierr /= dust_coupling_err_shape .or. any(bad_heat/=0d0)) stop 33
    write(*,'(a)') 'NATIVE_DUST_COUPLING_OK proportional=1 saturation=1 zero_dust=1 heating=1 invalid=1'
  end subroutine coupling_checks

  subroutine assert_unchanged()
    if (any(transfer(energy,[0_int64],size(energy))/=transfer(before,[0_int64],size(before)))) stop 11
    if (any(transfer(photons,[0_int64],size(photons))/=transfer(before_p,[0_int64],size(before_p)))) stop 12
    if (any(transfer(temperature,[0_int64],size(temperature))/=transfer(before_t,[0_int64],size(before_t)))) stop 13
    ! Compare fields, not derived-type padding bytes.
    if (any(transfer([diag%escaped_erg,diag%absorbed_erg,diag%primary_erg, &
         diag%balance_relative,diag%local_relative],[0_int64],5)/= &
         transfer([before_diag%escaped_erg,before_diag%absorbed_erg,before_diag%primary_erg, &
         before_diag%balance_relative,before_diag%local_relative],[0_int64],5))) stop 14
    if (diag%iterations/=before_diag%iterations) stop 15
  end subroutine
end program

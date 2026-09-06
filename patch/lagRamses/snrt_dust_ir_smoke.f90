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
contains
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

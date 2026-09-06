! Transaction and bounded fixed-point primitives for the native SNRT path.
!
! This module deliberately operates on explicit arrays.  The RAMSES driver is
! responsible only for preparing the arrays and committing the result after a
! collective decision.  No trial is allowed to mutate persistent chemistry or
! radiation state through this module.
module snrt_rt_transaction
  use, intrinsic :: iso_c_binding, only: c_float
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: dp
  implicit none

  private

  integer, parameter, public :: snrt_transaction_ok = 0
  integer, parameter, public :: snrt_transaction_err_shape = 1
  integer, parameter, public :: snrt_transaction_err_state = 2
  integer, parameter, public :: snrt_transaction_err_config = 3
  integer, parameter, public :: snrt_transaction_err_mpi = 4
  integer, parameter, public :: snrt_transaction_contract_version = 1
  integer, parameter, public :: snrt_transaction_max_iterations_limit = 32

  integer, parameter, public :: snrt_failure_none = 0
  integer, parameter, public :: snrt_failure_partition = 1
  integer, parameter, public :: snrt_failure_chemistry = 2
  integer, parameter, public :: snrt_failure_receiver = 3
  integer, parameter, public :: snrt_failure_transport = 4
  integer, parameter, public :: snrt_failure_convergence = 5
  integer, parameter, public :: snrt_failure_unassigned = 6

  type, public :: snrt_rt_iteration_config
     integer :: max_iterations = snrt_transaction_max_iterations_limit
     real(dp) :: fraction_absolute_tolerance = 1.0d-6
     real(dp) :: tau_relative_tolerance = 1.0d-5
     real(dp) :: tau_floor = 1.0d-12
     real(dp) :: relaxation = 0.5d0
     integer :: failure_stage = snrt_failure_none
     integer :: failure_leaf = 0
  end type snrt_rt_iteration_config

  type, public :: snrt_rt_transaction_snapshot
     logical :: active = .false.
     real(c_float), allocatable :: photon_before(:,:,:)
     real(dp), allocatable :: hydrogen_ii_before(:)
     real(dp), allocatable :: helium_ii_before(:)
     real(dp), allocatable :: helium_iii_before(:)
     real(dp), allocatable :: neutral_hydrogen_before(:)
     real(dp), allocatable :: thermal_before(:)
  end type snrt_rt_transaction_snapshot

  public :: snrt_transaction_load_config
  public :: snrt_transaction_begin
  public :: snrt_transaction_restore
  public :: snrt_transaction_commit_level
  public :: snrt_transaction_check_convergence
  public :: snrt_transaction_reduce_decision
  public :: snrt_transaction_reduce_sum
  public :: snrt_transaction_failure_requested
  public :: snrt_transaction_failure_name
  public :: snrt_transaction_error_message

contains

  subroutine snrt_transaction_load_config(production_allowed, config, ierr, message)
    logical, intent(in) :: production_allowed
    type(snrt_rt_iteration_config), intent(out) :: config
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: message
    character(len=1024) :: value
    integer :: length, status, read_status

    config = snrt_rt_iteration_config()
    ierr = snrt_transaction_ok
    message = ''

    call get_environment_variable('SNRT_RT_TX_MAX_ITER', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%max_iterations
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_MAX_ITER is not an integer'
          return
       end if
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_FRACTION_TOL', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%fraction_absolute_tolerance
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_FRACTION_TOL is not a number'
          return
       end if
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_TAU_TOL', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%tau_relative_tolerance
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_TAU_TOL is not a number'
          return
       end if
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_TAU_FLOOR', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%tau_floor
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_TAU_FLOOR is not a number'
          return
       end if
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_RELAXATION', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%relaxation
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_RELAXATION is not a number'
          return
       end if
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_TEST_FAIL_STAGE', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       select case (trim(value(1:length)))
       case ('none')
          config%failure_stage = snrt_failure_none
       case ('partition')
          config%failure_stage = snrt_failure_partition
       case ('chemistry')
          config%failure_stage = snrt_failure_chemistry
       case ('receiver')
          config%failure_stage = snrt_failure_receiver
       case default
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_TEST_FAIL_STAGE is unknown'
          return
       end select
    end if
    value = ''
    call get_environment_variable('SNRT_RT_TX_TEST_FAIL_LEAF', value, length=length, status=status)
    if (status == 0 .and. length > 0) then
       read(value(1:length),*,iostat=read_status) config%failure_leaf
       if (read_status /= 0) then
          ierr = snrt_transaction_err_config
          message = 'SNRT_RT_TX_TEST_FAIL_LEAF is not an integer'
          return
       end if
    end if

    if (config%max_iterations < 1 .or. &
         config%max_iterations > snrt_transaction_max_iterations_limit .or. &
         config%fraction_absolute_tolerance < 0.0d0 .or. &
         config%tau_relative_tolerance < 0.0d0 .or. config%tau_floor <= 0.0d0 .or. &
         config%relaxation <= 0.0d0 .or. config%relaxation > 1.0d0 .or. &
         config%failure_leaf < 0) then
       ierr = snrt_transaction_err_config
       message = 'native RT transaction configuration is outside its declared domain (max 32 trials)'
       return
    end if
    if (config%failure_stage /= snrt_failure_none .and. config%failure_leaf <= 0) then
       ierr = snrt_transaction_err_config
       message = 'failure injection requires a positive leaf index'
       return
    end if
    if (production_allowed .and. config%failure_stage /= snrt_failure_none) then
       ierr = snrt_transaction_err_config
       message = 'failure injection is forbidden under an approved production contract'
       return
    end if
  end subroutine snrt_transaction_load_config

  subroutine snrt_transaction_begin(transaction, persistent_intensity, leaf_slot, &
       hydrogen_ii, helium_ii, helium_iii, neutral_hydrogen, thermal, ierr)
    type(snrt_rt_transaction_snapshot), intent(inout) :: transaction
    real(c_float), intent(in) :: persistent_intensity(:,:,:)
    integer, intent(in) :: leaf_slot(:)
    real(dp), intent(in) :: hydrogen_ii(:), helium_ii(:), helium_iii(:)
    real(dp), intent(in) :: neutral_hydrogen(:), thermal(:)
    integer, intent(out) :: ierr
    integer :: i, islot, nleaf, max_slot

    call snrt_transaction_finalize(transaction)
    ierr = snrt_transaction_ok
    nleaf = size(leaf_slot)
    max_slot = 0
    if (nleaf > 0) max_slot = maxval(leaf_slot)
    if (size(persistent_intensity,1) < 1 .or. size(persistent_intensity,2) < 1 .or. &
         (nleaf > 0 .and. any(leaf_slot < 1)) .or. max_slot > size(persistent_intensity,3) .or. &
         size(hydrogen_ii) < max_slot .or. size(helium_ii) < max_slot .or. &
         size(helium_iii) < max_slot .or. size(neutral_hydrogen) < max_slot .or. &
         size(thermal) < nleaf) then
       ierr = snrt_transaction_err_shape
       return
    end if
    allocate(transaction%photon_before(size(persistent_intensity,1), &
         size(persistent_intensity,2), nleaf), &
         transaction%hydrogen_ii_before(nleaf), transaction%helium_ii_before(nleaf), &
         transaction%helium_iii_before(nleaf), transaction%neutral_hydrogen_before(nleaf), &
         transaction%thermal_before(nleaf))
    do i = 1, nleaf
       islot = leaf_slot(i)
       transaction%photon_before(:,:,i) = persistent_intensity(:,:,islot)
       transaction%hydrogen_ii_before(i) = hydrogen_ii(islot)
       transaction%helium_ii_before(i) = helium_ii(islot)
       transaction%helium_iii_before(i) = helium_iii(islot)
       transaction%neutral_hydrogen_before(i) = neutral_hydrogen(islot)
       transaction%thermal_before(i) = thermal(i)
    end do
    transaction%active = .true.
  end subroutine snrt_transaction_begin

  subroutine snrt_transaction_restore(transaction, persistent_intensity, leaf_slot, &
       hydrogen_ii, helium_ii, helium_iii, neutral_hydrogen, thermal, ierr)
    type(snrt_rt_transaction_snapshot), intent(inout) :: transaction
    real(c_float), intent(inout) :: persistent_intensity(:,:,:)
    integer, intent(in) :: leaf_slot(:)
    real(dp), intent(inout) :: hydrogen_ii(:), helium_ii(:), helium_iii(:)
    real(dp), intent(inout) :: neutral_hydrogen(:), thermal(:)
    integer, intent(out) :: ierr
    integer :: i, islot, nleaf, max_slot

    ierr = snrt_transaction_ok
    nleaf = size(leaf_slot)
    max_slot = 0
    if (nleaf > 0) max_slot = maxval(leaf_slot)
    if (.not. transaction%active .or. size(transaction%photon_before,3) /= nleaf .or. &
         size(hydrogen_ii) < max_slot .or. size(helium_ii) < max_slot .or. &
         size(helium_iii) < max_slot .or. size(neutral_hydrogen) < max_slot .or. &
         size(thermal) < nleaf .or. any(leaf_slot < 1) .or. &
         max_slot > size(persistent_intensity,3)) then
       ierr = snrt_transaction_err_state
       return
    end if
    do i = 1, nleaf
       islot = leaf_slot(i)
       persistent_intensity(:,:,islot) = transaction%photon_before(:,:,i)
       hydrogen_ii(islot) = transaction%hydrogen_ii_before(i)
       helium_ii(islot) = transaction%helium_ii_before(i)
       helium_iii(islot) = transaction%helium_iii_before(i)
       neutral_hydrogen(islot) = transaction%neutral_hydrogen_before(i)
       thermal(i) = transaction%thermal_before(i)
    end do
    call snrt_transaction_finalize(transaction)
  end subroutine snrt_transaction_restore

  subroutine snrt_transaction_commit_level(transaction, persistent_intensity, leaf_slot, &
       persistent_hydrogen_ii, persistent_helium_ii, persistent_helium_iii, &
       persistent_neutral_hydrogen, trial_intensity, coarse_flux_trial, &
       trial_hydrogen_ii, trial_helium_ii, trial_helium_iii, trial_neutral_hydrogen, &
       thermal, trial_thermal, ierr)
    type(snrt_rt_transaction_snapshot), intent(inout) :: transaction
    real(c_float), intent(inout) :: persistent_intensity(:,:,:)
    integer, intent(in) :: leaf_slot(:)
    real(dp), intent(inout) :: persistent_hydrogen_ii(:), persistent_helium_ii(:)
    real(dp), intent(inout) :: persistent_helium_iii(:), persistent_neutral_hydrogen(:)
    real(c_float), intent(in) :: trial_intensity(:,:,:), coarse_flux_trial(:,:,:)
    real(dp), intent(in) :: trial_hydrogen_ii(:), trial_helium_ii(:), trial_helium_iii(:)
    real(dp), intent(in) :: trial_neutral_hydrogen(:), trial_thermal(:)
    real(dp), intent(inout) :: thermal(:)
    integer, intent(out) :: ierr
    integer :: i, islot, nleaf, max_slot

    ierr = snrt_transaction_ok
    nleaf = size(leaf_slot)
    max_slot = 0
    if (nleaf > 0) max_slot = maxval(leaf_slot)
    if (.not. transaction%active .or. size(trial_intensity,1) /= size(persistent_intensity,1) .or. &
         size(trial_intensity,2) /= size(persistent_intensity,2) .or. &
         size(trial_intensity,3) /= nleaf .or. size(coarse_flux_trial,1) /= size(persistent_intensity,1) .or. &
         size(coarse_flux_trial,2) /= size(persistent_intensity,2) .or. &
         size(coarse_flux_trial,3) /= size(persistent_intensity,3) .or. &
         size(trial_hydrogen_ii) < nleaf .or. size(trial_helium_ii) < nleaf .or. &
         size(trial_helium_iii) < nleaf .or. size(trial_neutral_hydrogen) < nleaf .or. &
         size(trial_thermal) < nleaf .or. size(persistent_hydrogen_ii) < size(persistent_intensity,3) .or. &
         size(persistent_helium_ii) < size(persistent_intensity,3) .or. &
         size(persistent_helium_iii) < size(persistent_intensity,3) .or. &
         size(persistent_neutral_hydrogen) < size(persistent_intensity,3) .or. &
         size(thermal) < nleaf .or. any(leaf_slot < 1) .or. &
         max_slot > size(persistent_intensity,3)) then
       ierr = snrt_transaction_err_shape
       return
    end if
    do i = 1, nleaf
       islot = leaf_slot(i)
       persistent_intensity(:,:,islot) = trial_intensity(:,:,i)
       persistent_hydrogen_ii(islot) = trial_hydrogen_ii(i)
       persistent_helium_ii(islot) = trial_helium_ii(i)
       persistent_helium_iii(islot) = trial_helium_iii(i)
       persistent_neutral_hydrogen(islot) = trial_neutral_hydrogen(i)
       thermal(i) = trial_thermal(i)
    end do
    persistent_intensity = persistent_intensity + coarse_flux_trial
    call snrt_transaction_finalize(transaction)
  end subroutine snrt_transaction_commit_level

  subroutine snrt_transaction_check_convergence(current_fraction, target_fraction, &
       current_tau, target_tau, config, residual, converged, ierr)
    real(dp), intent(in) :: current_fraction(:,:), target_fraction(:,:)
    real(c_float), intent(in) :: current_tau(:,:), target_tau(:,:)
    type(snrt_rt_iteration_config), intent(in) :: config
    real(dp), intent(out) :: residual
    logical, intent(out) :: converged
    integer, intent(out) :: ierr
    integer :: i, j
    real(dp) :: fraction_residual, tau_residual, denominator, current_value, target_value

    ierr = snrt_transaction_ok
    residual = 0.0d0
    converged = .false.
    if (size(current_fraction,1) /= size(target_fraction,1) .or. &
         size(current_fraction,2) /= size(target_fraction,2) .or. &
         size(current_tau,1) /= size(target_tau,1) .or. &
         size(current_tau,2) /= size(target_tau,2) .or. &
         size(current_fraction,2) /= 3 .or. config%fraction_absolute_tolerance < 0.0d0 .or. &
         config%tau_relative_tolerance < 0.0d0 .or. config%tau_floor <= 0.0d0) then
       ierr = snrt_transaction_err_shape
       return
    end if
    fraction_residual = 0.0d0
    if (size(current_fraction) > 0) then
       if (any(.not. ieee_is_finite(current_fraction)) .or. &
            any(.not. ieee_is_finite(target_fraction)) .or. &
            any(current_fraction < 0.0d0) .or. any(target_fraction < 0.0d0)) then
          ierr = snrt_transaction_err_state
          return
       end if
       fraction_residual = maxval(abs(target_fraction-current_fraction))
    end if
    if (size(current_tau) > 0) then
       if (any(.not. ieee_is_finite(current_tau)) .or. &
            any(.not. ieee_is_finite(target_tau)) .or. any(current_tau < 0.0_c_float) .or. &
            any(target_tau < 0.0_c_float)) then
          ierr = snrt_transaction_err_state
          return
       end if
    end if
    tau_residual = 0.0d0
    do j = 1, size(current_tau,2)
       do i = 1, size(current_tau,1)
          current_value = max(0.0d0, real(current_tau(i,j),dp))
          target_value = max(0.0d0, real(target_tau(i,j),dp))
          if (max(current_value,target_value) <= config%tau_floor) cycle
          denominator = max(current_value, config%tau_floor)
          tau_residual = max(tau_residual, abs(target_value-current_value)/denominator)
       end do
    end do
    residual = max(fraction_residual, tau_residual)
    converged = fraction_residual <= config%fraction_absolute_tolerance .and. &
         tau_residual <= config%tau_relative_tolerance
  end subroutine snrt_transaction_check_convergence

  subroutine snrt_transaction_reduce_decision(local_failed, local_converged, local_residual, &
       global_failed, global_converged, global_residual, ierr)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    integer, intent(in) :: local_failed, local_converged
    real(dp), intent(in) :: local_residual
    integer, intent(out) :: global_failed, global_converged
    real(dp), intent(out) :: global_residual
    integer, intent(out) :: ierr
    integer :: info

    ierr = snrt_transaction_ok
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_failed, global_failed, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
    if (info /= 0) then
       ierr = snrt_transaction_err_mpi
       return
    end if
    call MPI_ALLREDUCE(local_converged, global_converged, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, info)
    if (info /= 0) then
       ierr = snrt_transaction_err_mpi
       return
    end if
    call MPI_ALLREDUCE(local_residual, global_residual, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         MPI_COMM_WORLD, info)
    if (info /= 0) ierr = snrt_transaction_err_mpi
#else
    global_failed = local_failed
    global_converged = local_converged
    global_residual = local_residual
#endif
  end subroutine snrt_transaction_reduce_decision

  subroutine snrt_transaction_reduce_sum(local_value, global_value, ierr)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    real(dp), intent(in) :: local_value
    real(dp), intent(out) :: global_value
    integer, intent(out) :: ierr
    integer :: info

    ierr = snrt_transaction_ok
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_value, global_value, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, info)
    if (info /= 0) ierr = snrt_transaction_err_mpi
#else
    global_value = local_value
#endif
  end subroutine snrt_transaction_reduce_sum

  logical function snrt_transaction_failure_requested(config, stage, leaf) result(requested)
    type(snrt_rt_iteration_config), intent(in) :: config
    integer, intent(in) :: stage, leaf

    requested = config%failure_stage == stage .and. config%failure_leaf == leaf
  end function snrt_transaction_failure_requested

  function snrt_transaction_failure_name(failure) result(name)
    integer, intent(in) :: failure
    character(len=32) :: name

    select case (failure)
    case (snrt_failure_none)
       name = 'none'
    case (snrt_failure_partition)
       name = 'partition'
    case (snrt_failure_chemistry)
       name = 'chemistry'
    case (snrt_failure_receiver)
       name = 'receiver'
    case (snrt_failure_transport)
       name = 'transport'
    case (snrt_failure_convergence)
       name = 'convergence'
    case (snrt_failure_unassigned)
       name = 'unassigned_absorption'
    case default
       name = 'unknown'
    end select
  end function snrt_transaction_failure_name

  function snrt_transaction_error_message(ierr) result(message)
    integer, intent(in) :: ierr
    character(len=80) :: message

    select case (ierr)
    case (snrt_transaction_ok)
       message = 'ok'
    case (snrt_transaction_err_shape)
       message = 'transaction array shape mismatch'
    case (snrt_transaction_err_state)
       message = 'transaction state is inactive or inconsistent'
    case (snrt_transaction_err_config)
       message = 'invalid native RT transaction configuration'
    case (snrt_transaction_err_mpi)
       message = 'MPI collective decision failed'
    case default
       message = 'unknown native RT transaction error'
    end select
  end function snrt_transaction_error_message

  subroutine snrt_transaction_finalize(transaction)
    type(snrt_rt_transaction_snapshot), intent(inout) :: transaction

    if (allocated(transaction%photon_before)) deallocate(transaction%photon_before)
    if (allocated(transaction%hydrogen_ii_before)) deallocate(transaction%hydrogen_ii_before)
    if (allocated(transaction%helium_ii_before)) deallocate(transaction%helium_ii_before)
    if (allocated(transaction%helium_iii_before)) deallocate(transaction%helium_iii_before)
    if (allocated(transaction%neutral_hydrogen_before)) deallocate(transaction%neutral_hydrogen_before)
    if (allocated(transaction%thermal_before)) deallocate(transaction%thermal_before)
    transaction%active = .false.
  end subroutine snrt_transaction_finalize

end module snrt_rt_transaction

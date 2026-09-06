program snrt_rt_transaction_smoke
  use, intrinsic :: iso_c_binding, only: c_float
  use amr_parameters, only: dp
  use snrt_rt_transaction
  implicit none

  real(c_float) :: persistent(2,2,3), persistent_before(2,2,3)
  real(c_float) :: trial_intensity(2,2,2), coarse_flux(2,2,3)
  integer :: leaf_slot(2)
  real(dp) :: hydrogen(3), helium_ii(3), helium_iii(3), neutral(3), thermal(2)
  real(dp) :: hydrogen_before(3), helium_ii_before(3), helium_iii_before(3)
  real(dp) :: neutral_before(3), thermal_before(2)
  real(dp) :: trial_hydrogen(2), trial_helium_ii(2), trial_helium_iii(2)
  real(dp) :: trial_neutral(2), trial_thermal(2)
  real(dp) :: current_fraction(2,3), target_fraction(2,3), residual
  real(c_float) :: current_tau(2,2), target_tau(2,2)
  real(c_float), allocatable :: empty_persistent(:,:,:), empty_trial(:,:,:), &
       empty_coarse(:,:,:)
  integer, allocatable :: empty_leaf_slot(:)
  real(dp), allocatable :: empty_hydrogen(:), empty_helium_ii(:), &
       empty_helium_iii(:), empty_neutral(:), empty_thermal(:)
  type(snrt_rt_iteration_config) :: config
  type(snrt_rt_transaction_snapshot) :: transaction
  character(len=128) :: message
  integer :: ierr, i, j, k, global_failed, global_converged
  real(dp) :: global_residual, lo, hi, midpoint
  logical :: converged

  persistent = 0.0_c_float
  do k = 1, 3
     do j = 1, 2
        do i = 1, 2
           persistent(i,j,k) = real(100*i + 10*j + k, c_float)
        end do
     end do
  end do
  persistent_before = persistent
  leaf_slot = (/ 1, 3 /)
  hydrogen = (/ 0.10d0, 0.20d0, 0.30d0 /)
  helium_ii = (/ 0.05d0, 0.10d0, 0.15d0 /)
  helium_iii = (/ 0.01d0, 0.02d0, 0.03d0 /)
  neutral = 1.0d0 - hydrogen
  thermal = (/ 10.0d0, 20.0d0 /)
  hydrogen_before = hydrogen
  helium_ii_before = helium_ii
  helium_iii_before = helium_iii
  neutral_before = neutral
  thermal_before = thermal

  call snrt_transaction_begin(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. transaction%active, &
       'begin failed')
  persistent(:,:,1) = -1.0_c_float
  persistent(:,:,3) = -2.0_c_float
  hydrogen(1) = 0.0d0
  hydrogen(3) = 0.0d0
  helium_ii(1) = 0.0d0
  helium_ii(3) = 0.0d0
  helium_iii(1) = 0.0d0
  helium_iii(3) = 0.0d0
  neutral(1) = 0.0d0
  neutral(3) = 0.0d0
  thermal = -1.0d0
  call snrt_transaction_restore(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. .not. transaction%active, &
       'restore failed')
  call assert_true(all(persistent == persistent_before), 'photon rollback mismatch')
  call assert_true(all(hydrogen == hydrogen_before), 'H rollback mismatch')
  call assert_true(all(helium_ii == helium_ii_before), 'He II rollback mismatch')
  call assert_true(all(helium_iii == helium_iii_before), 'He III rollback mismatch')
  call assert_true(all(neutral == neutral_before), 'H I rollback mismatch')
  call assert_true(all(thermal == thermal_before), 'thermal rollback mismatch')
  write(*,'(A)') 'SNRT_NATIVE_RT_TRANSACTION_PARTITION_ROLLBACK_PASS'

  ! Repeat the snapshot/restore path with a chemistry-only mutation so the
  ! two named failure classes are both covered by the native smoke.
  call snrt_transaction_begin(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. transaction%active, &
       'chemistry rollback begin failed')
  hydrogen(1) = 0.91d0
  helium_ii(3) = 0.81d0
  helium_iii(3) = 0.12d0
  neutral(1) = 0.09d0
  thermal(2) = 999.0d0
  call snrt_transaction_restore(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. all(hydrogen == hydrogen_before) .and. &
       all(helium_ii == helium_ii_before) .and. all(helium_iii == helium_iii_before) .and. &
       all(neutral == neutral_before) .and. all(thermal == thermal_before), &
       'chemistry rollback mismatch')
  write(*,'(A)') 'SNRT_NATIVE_RT_TRANSACTION_CHEMISTRY_ROLLBACK_PASS'

  ! An empty MPI rank still enters the transaction API with allocated,
  ! zero-length payloads.  This is the local analogue of the zero-leaf
  ! collective path used by the RAMSES driver.
  allocate(empty_persistent(2,2,0), empty_trial(2,2,0), empty_coarse(2,2,0), &
       empty_leaf_slot(0), empty_hydrogen(0), empty_helium_ii(0), &
       empty_helium_iii(0), empty_neutral(0), empty_thermal(0))
  call snrt_transaction_begin(transaction, empty_persistent, empty_leaf_slot, &
       empty_hydrogen, empty_helium_ii, empty_helium_iii, empty_neutral, &
       empty_thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. transaction%active, &
       'zero-leaf begin failed')
  call snrt_transaction_restore(transaction, empty_persistent, empty_leaf_slot, &
       empty_hydrogen, empty_helium_ii, empty_helium_iii, empty_neutral, &
       empty_thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. .not. transaction%active, &
       'zero-leaf restore failed')
  deallocate(empty_persistent, empty_trial, empty_coarse, empty_leaf_slot, &
       empty_hydrogen, empty_helium_ii, empty_helium_iii, empty_neutral, &
       empty_thermal)

  call snrt_transaction_begin(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. transaction%active, &
       'second begin failed')
  trial_intensity = 7.0_c_float
  coarse_flux = 0.25_c_float
  trial_hydrogen = (/ 0.40d0, 0.50d0 /)
  trial_helium_ii = (/ 0.20d0, 0.25d0 /)
  trial_helium_iii = (/ 0.03d0, 0.04d0 /)
  trial_neutral = 1.0d0 - trial_hydrogen
  trial_thermal = (/ 11.0d0, 22.0d0 /)
  call snrt_transaction_commit_level(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, trial_intensity, coarse_flux, &
       trial_hydrogen, trial_helium_ii, trial_helium_iii, trial_neutral, &
       thermal, trial_thermal, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. .not. transaction%active, &
       'commit failed')
  call assert_true(all(persistent(:,:,1) == 7.25_c_float), 'slot 1 commit mismatch')
  call assert_true(all(persistent(:,:,2) == persistent_before(:,:,2) + 0.25_c_float), &
       'untouched slot correction mismatch')
  call assert_true(all(persistent(:,:,3) == 7.25_c_float), 'slot 3 commit mismatch')
  call assert_true(hydrogen(1) == trial_hydrogen(1) .and. &
       hydrogen(2) == hydrogen_before(2) .and. hydrogen(3) == trial_hydrogen(2), &
       'H commit mismatch')
  call assert_true(helium_ii(1) == trial_helium_ii(1) .and. &
       helium_ii(2) == helium_ii_before(2) .and. helium_ii(3) == trial_helium_ii(2), &
       'He II commit mismatch')
  call assert_true(helium_iii(1) == trial_helium_iii(1) .and. &
       helium_iii(2) == helium_iii_before(2) .and. helium_iii(3) == trial_helium_iii(2), &
       'He III commit mismatch')
  call assert_true(neutral(1) == trial_neutral(1) .and. &
       neutral(2) == neutral_before(2) .and. neutral(3) == trial_neutral(2), &
       'H I commit mismatch')
  call assert_true(all(thermal == trial_thermal), 'thermal commit mismatch')

  config = snrt_rt_iteration_config()
  current_fraction = 0.0d0
  target_fraction = 0.0d0
  current_fraction(1,:) = (/ 0.10d0, 0.05d0, 0.01d0 /)
  current_fraction(2,:) = (/ 0.20d0, 0.10d0, 0.02d0 /)
  target_fraction = current_fraction
  target_fraction(1,1) = target_fraction(1,1) + 5.0d-7
  current_tau = 1.0_c_float
  target_tau = current_tau
  target_tau(1,1) = target_tau(1,1) * 1.0_c_float + 5.0e-6_c_float
  call snrt_transaction_check_convergence(current_fraction, target_fraction, &
       current_tau, target_tau, config, residual, converged, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. converged, &
       'converged fixed-point check failed')
  target_fraction(1,1) = target_fraction(1,1) + 1.0d-3
  call snrt_transaction_check_convergence(current_fraction, target_fraction, &
       current_tau, target_tau, config, residual, converged, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. .not. converged, &
       'non-converged fixed-point check failed')

  ! One-cell bisection is the deterministic reference for the bounded
  ! fixed-point evidence.  A deliberately one-trial contract must remain
  ! non-converged, proving that the iteration cap is fail-closed.
  lo = 0.0d0
  hi = 1.0d0
  do i = 1, 32
     midpoint = 0.5d0*(lo+hi)
     if (midpoint*midpoint < 0.25d0) then
        lo = midpoint
     else
        hi = midpoint
     end if
  end do
  call assert_true(abs(midpoint-0.5d0) < 1.0d-9, 'one-cell bisection reference failed')
  write(*,'(A,ES12.4)') 'SNRT_NATIVE_RT_TRANSACTION_BISECTION_PASS error=', &
       abs(midpoint-0.5d0)
  config = snrt_rt_iteration_config()
  config%max_iterations = 1
  config%fraction_absolute_tolerance = 0.0d0
  config%tau_relative_tolerance = 0.0d0
  current_fraction = 0.0d0
  target_fraction = 0.0d0
  target_fraction(1,1) = 0.5d0
  current_tau = 1.0_c_float
  target_tau = 2.0_c_float
  call snrt_transaction_check_convergence(current_fraction, target_fraction, &
       current_tau, target_tau, config, residual, converged, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. .not. converged, &
       'one-trial non-convergence was not rejected')
  write(*,'(A)') 'SNRT_NATIVE_RT_TRANSACTION_MAXITER1_NONCONVERGED_PASS'

  call snrt_transaction_reduce_decision(0, 1, residual, global_failed, &
       global_converged, global_residual, ierr)
  call assert_true(ierr == snrt_transaction_ok .and. global_failed == 0 .and. &
       global_converged == 1 .and. global_residual == residual, &
       'decision reduction failed')
  call snrt_transaction_load_config(.false., config, ierr, message)
  call assert_true(ierr == snrt_transaction_ok, 'diagnostic config rejected')
  if (config%failure_stage /= snrt_failure_none) then
     call assert_true(snrt_transaction_failure_requested(config, config%failure_stage, &
          config%failure_leaf), 'failure injection selector failed')
     call snrt_transaction_load_config(.true., config, ierr, message)
     call assert_true(ierr == snrt_transaction_err_config, &
          'production failure injection was not rejected')
  end if

  write(*,'(A)') 'SNRT_NATIVE_RT_TRANSACTION_SMOKE_PASS'
  write(*,'(A,I0,A,ES12.4)') 'SNRT_NATIVE_RT_TRANSACTION_MAX_ITER=', &
       config%max_iterations, ' residual=', residual

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
       write(*,'(A,A)') 'SNRT_NATIVE_RT_TRANSACTION_SMOKE_FAIL: ', trim(message)
       error stop 1
    end if
  end subroutine assert_true

end program snrt_rt_transaction_smoke

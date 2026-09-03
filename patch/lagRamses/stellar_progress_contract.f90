! Transactional stellar feedback progress contract.
!
! indtab is the last age whose source was successfully deposited.  A pending
! age is never exported until the caller has completed all field and particle
! updates.  Repeating a committed age is therefore an exact no-op, while an
! interrupted interval can be retried without double deposition.

module stellar_progress_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp
  implicit none

  private

  integer, parameter, public :: progress_ok = 0
  integer, parameter, public :: progress_err_argument = 1
  integer, parameter, public :: progress_err_uninitialized = 2
  integer, parameter, public :: progress_err_busy = 4
  integer, parameter, public :: progress_err_stale = 8

  type, public :: stellar_progress_t
     logical :: initialized = .false.
     logical :: pending = .false.
     real(stellar_dp) :: committed_age_code = 0.0_stellar_dp
     real(stellar_dp) :: pending_age_code = 0.0_stellar_dp
     real(stellar_dp) :: pending_dt_code = 0.0_stellar_dp
  end type stellar_progress_t

  public :: progress_initialize
  public :: progress_begin
  public :: progress_commit
  public :: progress_abort
  public :: progress_export

contains

  subroutine progress_initialize(progress, committed_age_code, ierr)
    type(stellar_progress_t), intent(out) :: progress
    real(stellar_dp), intent(in) :: committed_age_code
    integer, intent(out) :: ierr

    progress%initialized = .false.
    progress%pending = .false.
    progress%committed_age_code = 0.0_stellar_dp
    progress%pending_age_code = 0.0_stellar_dp
    progress%pending_dt_code = 0.0_stellar_dp
    ierr = progress_ok
    if (.not. ieee_is_finite(committed_age_code) .or. &
         committed_age_code < 0.0_stellar_dp) then
       ierr = progress_err_argument
       return
    end if
    progress%initialized = .true.
    progress%committed_age_code = committed_age_code
  end subroutine progress_initialize

  subroutine progress_begin(progress, current_age_code, tolerance, &
       should_deposit, dt_code, ierr)
    type(stellar_progress_t), intent(inout) :: progress
    real(stellar_dp), intent(in) :: current_age_code, tolerance
    logical, intent(out) :: should_deposit
    real(stellar_dp), intent(out) :: dt_code
    integer, intent(out) :: ierr
    real(stellar_dp) :: tol

    should_deposit = .false.
    dt_code = 0.0_stellar_dp
    ierr = progress_ok
    if (.not. progress%initialized) then
       ierr = progress_err_uninitialized
       return
    end if
    if (progress%pending) then
       ierr = progress_err_busy
       return
    end if
    if (.not. ieee_is_finite(current_age_code) .or. &
         .not. ieee_is_finite(tolerance) .or. current_age_code < 0.0_stellar_dp &
         .or. tolerance < 0.0_stellar_dp) then
       ierr = progress_err_argument
       return
    end if

    tol = max(tolerance, 0.0_stellar_dp)
    if (current_age_code < progress%committed_age_code - tol) then
       ierr = progress_err_stale
       return
    end if
    if (current_age_code <= progress%committed_age_code + tol) return

    dt_code = current_age_code - progress%committed_age_code
    if (.not. ieee_is_finite(dt_code) .or. dt_code <= 0.0_stellar_dp) then
       ierr = progress_err_argument
       dt_code = 0.0_stellar_dp
       return
    end if
    progress%pending = .true.
    progress%pending_age_code = current_age_code
    progress%pending_dt_code = dt_code
    should_deposit = .true.
  end subroutine progress_begin

  subroutine progress_commit(progress, ierr)
    type(stellar_progress_t), intent(inout) :: progress
    integer, intent(out) :: ierr

    ierr = progress_ok
    if (.not. progress%initialized .or. .not. progress%pending) then
       ierr = merge(progress_err_uninitialized, progress_err_argument, &
            .not. progress%initialized)
       return
    end if
    progress%committed_age_code = progress%pending_age_code
    progress%pending_age_code = 0.0_stellar_dp
    progress%pending_dt_code = 0.0_stellar_dp
    progress%pending = .false.
  end subroutine progress_commit

  subroutine progress_abort(progress, ierr)
    type(stellar_progress_t), intent(inout) :: progress
    integer, intent(out) :: ierr

    ierr = progress_ok
    if (.not. progress%initialized) then
       ierr = progress_err_uninitialized
       return
    end if
    progress%pending_age_code = 0.0_stellar_dp
    progress%pending_dt_code = 0.0_stellar_dp
    progress%pending = .false.
  end subroutine progress_abort

  subroutine progress_export(progress, committed_age_code, ierr)
    type(stellar_progress_t), intent(in) :: progress
    real(stellar_dp), intent(out) :: committed_age_code
    integer, intent(out) :: ierr

    committed_age_code = 0.0_stellar_dp
    ierr = progress_ok
    if (.not. progress%initialized) then
       ierr = progress_err_uninitialized
       return
    end if
    committed_age_code = progress%committed_age_code
  end subroutine progress_export

end module stellar_progress_contract

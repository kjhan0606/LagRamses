! Restart accounting for the approved zero-remnant SNIa baseline.
!
! For a committed SNIa interval, the persisted particle mass and the generic
! cumulative population ledger are sufficient to recover the already-applied
! SNIa return.  Keep this arithmetic isolated and unit-tested: a future
! nonzero-terminal-remnant model must provide a versioned per-particle debit
! payload instead of reusing this invariant.

module stellar_snia_runtime_accounting
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp
  implicit none
  private

  integer, parameter, public :: snia_accounting_ok = 0
  integer, parameter, public :: snia_accounting_err_argument = 1
  integer, parameter, public :: snia_accounting_err_nonfinite = 2
  integer, parameter, public :: snia_accounting_err_inconsistent = 4

  public :: reconstruct_prior_snia_return

contains

  subroutine reconstruct_prior_snia_return(particle_mass_before_code, &
       generic_interval_return_code, generic_remaining_code, particle_mass_scale, &
       tolerance, prior_snia_returned_code, ierr)
    real(stellar_dp), intent(in) :: particle_mass_before_code
    real(stellar_dp), intent(in) :: generic_interval_return_code
    real(stellar_dp), intent(in) :: generic_remaining_code
    real(stellar_dp), intent(in) :: particle_mass_scale, tolerance
    real(stellar_dp), intent(out) :: prior_snia_returned_code
    integer, intent(out) :: ierr

    real(stellar_dp) :: tol, scale

    prior_snia_returned_code = 0.0_stellar_dp
    ierr = snia_accounting_ok
    if (.not. all(ieee_is_finite((/particle_mass_before_code, &
         generic_interval_return_code, generic_remaining_code, &
         particle_mass_scale, tolerance/))) .or. &
         particle_mass_scale <= 0.0_stellar_dp .or. tolerance < 0.0_stellar_dp) then
       ierr = snia_accounting_err_argument
       return
    end if

    tol = max(tolerance, 1.0e-12_stellar_dp)
    scale = max(particle_mass_scale, abs(particle_mass_before_code), &
         abs(generic_interval_return_code), abs(generic_remaining_code))
    if (particle_mass_before_code < -tol * scale .or. &
         generic_interval_return_code < -tol * scale .or. &
         generic_remaining_code < -tol * scale) then
       ierr = snia_accounting_err_inconsistent
       return
    end if

    prior_snia_returned_code = generic_remaining_code - &
         (particle_mass_before_code - generic_interval_return_code)
    if (.not. ieee_is_finite(prior_snia_returned_code)) then
       prior_snia_returned_code = 0.0_stellar_dp
       ierr = snia_accounting_err_nonfinite
       return
    end if
    if (prior_snia_returned_code < -tol * scale .or. &
         prior_snia_returned_code > generic_remaining_code + tol * scale) then
       prior_snia_returned_code = 0.0_stellar_dp
       ierr = snia_accounting_err_inconsistent
       return
    end if
    prior_snia_returned_code = max(0.0_stellar_dp, prior_snia_returned_code)
  end subroutine reconstruct_prior_snia_return

end module stellar_snia_runtime_accounting

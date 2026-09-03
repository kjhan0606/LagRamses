! Review-only Type-Ia delay-time distribution kernel.
!
! This module supplies the interval mathematics needed by F-P2.  It does not
! select a binary population, a physical delay range, an event normalization,
! or an event-yield source.  The production stellar driver must continue to
! reject SNIa until those inputs and an approval sidecar are present.

module stellar_snia_dtd
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp
  implicit none

  private

  integer, parameter, public :: snia_dtd_ok = 0
  integer, parameter, public :: snia_dtd_err_argument = 1
  integer, parameter, public :: snia_dtd_err_nonfinite = 2
  integer, parameter, public :: snia_dtd_err_normalization = 4

  public :: integrate_snia_dtd_interval
  public :: integrate_snia_dtd_cumulative

contains

  subroutine integrate_snia_dtd_interval(age_old_gyr, age_new_gyr, &
       minimum_delay_gyr, maximum_delay_gyr, power_law_index, &
       events_per_initial_msun, expected_events, ierr)
    real(stellar_dp), intent(in) :: age_old_gyr, age_new_gyr
    real(stellar_dp), intent(in) :: minimum_delay_gyr, maximum_delay_gyr
    real(stellar_dp), intent(in) :: power_law_index
    real(stellar_dp), intent(in) :: events_per_initial_msun
    real(stellar_dp), intent(out) :: expected_events
    integer, intent(out) :: ierr

    real(stellar_dp) :: lower, upper, numerator, denominator
    integer :: integral_ierr

    expected_events = 0.0_stellar_dp
    ierr = snia_dtd_ok
    if (.not. ieee_is_finite(age_old_gyr) .or. &
         .not. ieee_is_finite(age_new_gyr) .or. &
         .not. ieee_is_finite(minimum_delay_gyr) .or. &
         .not. ieee_is_finite(maximum_delay_gyr) .or. &
         .not. ieee_is_finite(power_law_index) .or. &
         .not. ieee_is_finite(events_per_initial_msun) .or. &
         age_old_gyr < 0.0_stellar_dp .or. &
         age_new_gyr < age_old_gyr .or. &
         minimum_delay_gyr <= 0.0_stellar_dp .or. &
         maximum_delay_gyr <= minimum_delay_gyr .or. &
         events_per_initial_msun < 0.0_stellar_dp) then
       ierr = snia_dtd_err_argument
       return
    end if

    call power_law_integral(minimum_delay_gyr, maximum_delay_gyr, &
         power_law_index, denominator, integral_ierr)
    if (integral_ierr /= snia_dtd_ok .or. denominator <= 0.0_stellar_dp) then
       ierr = snia_dtd_err_normalization
       return
    end if

    lower = max(age_old_gyr, minimum_delay_gyr)
    upper = min(age_new_gyr, maximum_delay_gyr)
    if (upper <= lower .or. events_per_initial_msun == 0.0_stellar_dp) return

    call power_law_integral(lower, upper, power_law_index, numerator, &
         integral_ierr)
    if (integral_ierr /= snia_dtd_ok .or. numerator < 0.0_stellar_dp) then
       ierr = snia_dtd_err_normalization
       return
    end if
    expected_events = events_per_initial_msun * numerator / denominator
    if (.not. ieee_is_finite(expected_events) .or. expected_events < 0.0_stellar_dp) then
       expected_events = 0.0_stellar_dp
       ierr = snia_dtd_err_nonfinite
    end if
  end subroutine integrate_snia_dtd_interval

  subroutine integrate_snia_dtd_cumulative(age_gyr, minimum_delay_gyr, &
       maximum_delay_gyr, power_law_index, events_per_initial_msun, &
       expected_events, ierr)
    real(stellar_dp), intent(in) :: age_gyr, minimum_delay_gyr
    real(stellar_dp), intent(in) :: maximum_delay_gyr, power_law_index
    real(stellar_dp), intent(in) :: events_per_initial_msun
    real(stellar_dp), intent(out) :: expected_events
    integer, intent(out) :: ierr

    call integrate_snia_dtd_interval(0.0_stellar_dp, age_gyr, &
         minimum_delay_gyr, maximum_delay_gyr, power_law_index, &
         events_per_initial_msun, expected_events, ierr)
  end subroutine integrate_snia_dtd_cumulative

  subroutine power_law_integral(lower, upper, power_law_index, integral, ierr)
    real(stellar_dp), intent(in) :: lower, upper, power_law_index
    real(stellar_dp), intent(out) :: integral
    integer, intent(out) :: ierr

    real(stellar_dp) :: exponent, log_lower, log_upper, log_ratio, z
    real(stellar_dp) :: exp_minus_one_over_z

    integral = 0.0_stellar_dp
    ierr = snia_dtd_ok
    if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. &
         .not. ieee_is_finite(power_law_index) .or. lower <= 0.0_stellar_dp .or. &
         upper <= lower) then
       ierr = snia_dtd_err_argument
       return
    end if
    exponent = power_law_index + 1.0_stellar_dp
    log_lower = log(lower)
    log_upper = log(upper)
    log_ratio = log_upper - log_lower
    if (exponent == 0.0_stellar_dp) then
       integral = log_ratio
    else
       ! Write the power-law primitive in a form that does not subtract two
       ! nearly equal powers when alpha is close to -1:
       !
       !   (u**e-l**e)/e = l**e log(u/l) [expm1(e log(u/l))/(e log(u/l))].
       !
       ! The local series for expm1(z)/z also avoids loss of significance for
       ! small z.  This is a numerical contract, not a choice of DTD model.
       z = exponent * log_ratio
       if (abs(z) <= 1.0e-4_stellar_dp) then
          exp_minus_one_over_z = 1.0_stellar_dp + z * (0.5_stellar_dp + &
               z * (1.0_stellar_dp / 6.0_stellar_dp + z * &
               (1.0_stellar_dp / 24.0_stellar_dp + z / &
               120.0_stellar_dp)))
       else
          exp_minus_one_over_z = (exp(z) - 1.0_stellar_dp) / z
       end if
       integral = exp(exponent * log_lower) * log_ratio * &
            exp_minus_one_over_z
    end if
    if (.not. ieee_is_finite(integral) .or. integral <= 0.0_stellar_dp) then
       integral = 0.0_stellar_dp
       ierr = snia_dtd_err_normalization
    end if
  end subroutine power_law_integral

end module stellar_snia_dtd

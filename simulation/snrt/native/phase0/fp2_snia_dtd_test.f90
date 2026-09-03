program fp2_snia_dtd_test
  use stellar_enrichment_config, only: stellar_dp
  use stellar_snia_dtd, only: snia_dtd_ok, snia_dtd_err_argument, &
       integrate_snia_dtd_interval, integrate_snia_dtd_cumulative
  implicit none

  integer :: failures, ierr, index
  real(stellar_dp) :: value, cumulative, sum_intervals, expected
  real(stellar_dp), parameter :: tmin = 0.04_stellar_dp
  real(stellar_dp), parameter :: tmax = 13.8_stellar_dp
  real(stellar_dp), parameter :: alpha = -1.0_stellar_dp
  real(stellar_dp), parameter :: normalization = 2.0e-3_stellar_dp
  real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp
  real(stellar_dp), parameter :: near_log_alpha = -0.99999999_stellar_dp
  real(stellar_dp), dimension(6), parameter :: age_edges = &
       (/0.0_stellar_dp, 0.02_stellar_dp, 0.04_stellar_dp, &
         0.1_stellar_dp, 1.0_stellar_dp, 13.8_stellar_dp/)

  failures = 0
  call integrate_snia_dtd_interval(0.0_stellar_dp, 0.02_stellar_dp, tmin, &
       tmax, alpha, normalization, value, ierr)
  call expect(ierr == snia_dtd_ok .and. value == 0.0_stellar_dp, &
       'pre-delay interval is empty', failures)

  call integrate_snia_dtd_interval(tmin, 0.1_stellar_dp, tmin, tmax, alpha, &
       normalization, value, ierr)
  expected = normalization * log(0.1_stellar_dp / tmin) / log(tmax / tmin)
  call expect(ierr == snia_dtd_ok .and. abs(value - expected) < tolerance, &
       'inverse-delay interval has analytic logarithmic integral', failures)

  call integrate_snia_dtd_cumulative(0.0_stellar_dp, tmin, tmax, alpha, &
       normalization, cumulative, ierr)
  call expect(ierr == snia_dtd_ok .and. cumulative == 0.0_stellar_dp, &
       'zero-age cumulative event count is zero', failures)

  call integrate_snia_dtd_cumulative(20.0_stellar_dp, tmin, tmax, alpha, &
       normalization, cumulative, ierr)
  call expect(ierr == snia_dtd_ok .and. abs(cumulative - normalization) < tolerance, &
       'post-support cumulative event count reaches normalization', failures)

  sum_intervals = 0.0_stellar_dp
  do index = 1, size(age_edges) - 1
     call integrate_snia_dtd_interval(age_edges(index), age_edges(index+1), &
          tmin, tmax, alpha, normalization, value, ierr)
     call expect(ierr == snia_dtd_ok, 'adjacent interval remains valid', failures)
     sum_intervals = sum_intervals + value
  end do
  call integrate_snia_dtd_cumulative(13.8_stellar_dp, tmin, tmax, alpha, &
       normalization, cumulative, ierr)
  call expect(abs(sum_intervals - cumulative) < tolerance, &
       'adjacent intervals telescope to cumulative count', failures)

  call integrate_snia_dtd_interval(0.1_stellar_dp, 1.0_stellar_dp, tmin, &
       tmax, alpha, normalization, value, ierr)
  call integrate_snia_dtd_cumulative(1.0_stellar_dp, tmin, tmax, alpha, &
       normalization, cumulative, ierr)
  call integrate_snia_dtd_cumulative(0.1_stellar_dp, tmin, tmax, alpha, &
       normalization, expected, ierr)
  call expect(abs(value + expected - cumulative) < tolerance, &
       'restart split has no hidden endpoint state', failures)

  call integrate_snia_dtd_interval(1.0_stellar_dp, 0.5_stellar_dp, tmin, &
       tmax, alpha, normalization, value, ierr)
  call expect(ierr == snia_dtd_err_argument .and. value == 0.0_stellar_dp, &
       'backward interval is rejected transactionally', failures)

  call integrate_snia_dtd_interval(0.0_stellar_dp, 1.0_stellar_dp, tmin, &
       tmax, alpha, -1.0_stellar_dp, value, ierr)
  call expect(ierr == snia_dtd_err_argument .and. value == 0.0_stellar_dp, &
       'negative event normalization is rejected', failures)

  call integrate_snia_dtd_interval(0.0_stellar_dp, 1.0_stellar_dp, 0.0_stellar_dp, &
       tmax, alpha, normalization, value, ierr)
  call expect(ierr == snia_dtd_err_argument .and. value == 0.0_stellar_dp, &
       'nonpositive minimum delay is rejected', failures)

  call integrate_snia_dtd_interval(0.0_stellar_dp, 1.0_stellar_dp, tmin, &
       tmax, -0.5_stellar_dp, normalization, value, ierr)
  expected = normalization * (sqrt(1.0_stellar_dp) - sqrt(tmin)) / &
       (sqrt(tmax) - sqrt(tmin))
  call expect(ierr == snia_dtd_ok .and. abs(value - expected) < tolerance, &
       'non-log power-law interval has analytic integral', failures)

  call integrate_snia_dtd_interval(tmin, 1.0_stellar_dp, tmin, tmax, near_log_alpha, &
       normalization, value, ierr)
  expected = normalization * near_log_reference(tmin, 1.0_stellar_dp, &
       near_log_alpha + 1.0_stellar_dp) / near_log_reference(tmin, tmax, &
       near_log_alpha + 1.0_stellar_dp)
  call expect(ierr == snia_dtd_ok .and. abs(value - expected) < &
       1.0e-15_stellar_dp, &
       'near-log DTD integral remains stable under cancellation test', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_DTD_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_DTD_TEST_OK'

contains

  subroutine expect(condition, label, failure_count)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failure_count

    if (condition) then
       write(*, '(a)') 'PASS: '//trim(label)
    else
       write(*, '(a)') 'FAIL: '//trim(label)
       failure_count = failure_count + 1
    end if
  end subroutine expect

  real(stellar_dp) function near_log_reference(lower, upper, exponent)
    real(stellar_dp), intent(in) :: lower, upper, exponent
    real(stellar_dp) :: log_lower, log_upper

    log_lower = log(lower)
    log_upper = log(upper)
    near_log_reference = (log_upper - log_lower) + exponent * &
         (log_upper**2 - log_lower**2) / 2.0_stellar_dp + exponent**2 * &
         (log_upper**3 - log_lower**3) / 6.0_stellar_dp + exponent**3 * &
         (log_upper**4 - log_lower**4) / 24.0_stellar_dp
  end function near_log_reference

end program fp2_snia_dtd_test

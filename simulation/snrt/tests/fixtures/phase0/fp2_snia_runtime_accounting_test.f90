program fp2_snia_runtime_accounting_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp
  use stellar_snia_runtime_accounting, only: &
       reconstruct_prior_snia_return, close_effective_snia_return, snia_accounting_ok, &
       snia_accounting_err_inconsistent
  implicit none

  real(stellar_dp) :: prior, tolerance, remaining, invalid
  integer :: ierr, failures

  failures = 0
  tolerance = 1.0e-12_stellar_dp

  ! First interval: initial mass 100, generic return 2, no prior Ia debit.
  call reconstruct_prior_snia_return(100.0_stellar_dp, 2.0_stellar_dp, &
       98.0_stellar_dp, 100.0_stellar_dp, tolerance, prior, ierr)
  call expect(ierr == snia_accounting_ok .and. prior == 0.0_stellar_dp, &
       'first interval reconstructs zero prior SNIa return', failures)

  ! After a committed 1.4-code-mass Ia return, a restart at the next interval
  ! must recover 1.4 even though the generic cumulative ledger advanced again.
  call reconstruct_prior_snia_return(96.6_stellar_dp, 0.5_stellar_dp, &
       97.5_stellar_dp, 100.0_stellar_dp, tolerance, prior, ierr)
  call expect(ierr == snia_accounting_ok, &
       'normal restart accounting accepts persisted state', failures)
  call expect_close(prior, 1.4_stellar_dp, &
       'normal restart recovers the prior SNIa debit', failures)

  ! The same arithmetic is independent of the absolute mass scale.
  call reconstruct_prior_snia_return(0.966_stellar_dp, 0.005_stellar_dp, &
       0.975_stellar_dp, 1.0_stellar_dp, tolerance, prior, ierr)
  call expect_close(prior, 0.014_stellar_dp, &
       'restart accounting scales in code units', failures)

  call reconstruct_prior_snia_return(0.0_stellar_dp, 2.0_stellar_dp, &
       98.0_stellar_dp, 100.0_stellar_dp, tolerance, prior, ierr)
  call expect(ierr == snia_accounting_err_inconsistent, &
       'inconsistent persisted mass is rejected', failures)

  call close_effective_snia_return(100d0,2d0,98d0,0d0,0d0,1.4d0,100d0,tolerance,remaining,ierr)
  call expect(ierr == snia_accounting_ok, 'effective SSP needs no asserted WD inventory', failures)
  call expect_close(remaining,96.6d0,'generic and Ia returns charged once',failures)
  call close_effective_snia_return(96.6d0,.5d0,97.5d0,1.4d0,1.4d0,.7d0,100d0,tolerance,remaining,ierr)
  call expect(ierr == snia_accounting_ok, 'effective SSP resumed interval closes', failures)
  call expect_close(remaining,95.4d0,'resumed return is not double charged',failures)
  call close_effective_snia_return(100d0,2.5d0,97.5d0,0d0,0d0,2.1d0,100d0,tolerance,remaining,ierr)
  call expect_close(remaining,95.4d0,'one interval equals two intervals',failures)
  call close_effective_snia_return(1d-6,2d-8,9.8d-7,0d0,0d0,1.4d-8,1d-6,tolerance,remaining,ierr)
  call expect(ierr == snia_accounting_ok,'small code-unit masses accepted',failures)
  call expect_close(remaining,9.66d-7,'no one-Msun floor in aggregate accounting',failures)
  call close_effective_snia_return(2d0,1d0,3d0,2d0,2d0,1.4d0,100d0,tolerance,remaining,ierr)
  call expect(ierr /= snia_accounting_ok .and. remaining == 0d0,'combined overdraw rejects atomically',failures)
  call close_effective_snia_return(96.6d0,.5d0,97.5d0,1.4d0,0d0,.7d0,100d0,tolerance,remaining,ierr)
  call expect(ierr /= snia_accounting_ok,'duplicate or inconsistent prior DTD debit rejects',failures)
  call close_effective_snia_return(98d0,.5d0,97.5d0,0d0,1.4d0,.7d0,100d0,tolerance,remaining,ierr)
  call expect(ierr /= snia_accounting_ok,'missing prior DTD debit rejects',failures)
  invalid=ieee_value(0d0,ieee_quiet_nan)
  call close_effective_snia_return(100d0,2d0,98d0,0d0,0d0,invalid,100d0,tolerance,remaining,ierr)
  call expect(ierr /= snia_accounting_ok .and. remaining == 0d0,'nonfinite Ia return rejects',failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_RUNTIME_ACCOUNTING_TEST_FAILED failures=', failures
     error stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_RUNTIME_ACCOUNTING_TEST_OK'

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

  subroutine expect_close(actual, expected, label, failure_count)
    real(stellar_dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failure_count
    real(stellar_dp) :: scale

    scale = max(1.0e-300_stellar_dp, abs(actual), abs(expected))
    call expect(abs(actual - expected) <= 1.0e-12_stellar_dp * scale, &
         label, failure_count)
  end subroutine expect_close

end program fp2_snia_runtime_accounting_test

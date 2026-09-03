program fp2_snia_population_contract_test
  use stellar_enrichment_config, only: stellar_dp, population_binary_ssp, &
       stellar_imf_kroupa
  use stellar_snia_population_contract, only: &
       snia_population_realization_t, snia_population_contract_ok, &
       snia_population_contract_err_unapproved, &
       snia_population_contract_err_parameter, &
       snia_population_contract_err_realization, snia_realization_expectation, &
       snia_realization_poisson, snia_metallicity_factor_supplied, &
       validate_snia_population_realization, evaluate_snia_interval_events
  implicit none

  type(snia_population_realization_t) :: realization
  real(stellar_dp) :: expected_events, repeat_events
  real(stellar_dp) :: first_interval, second_interval
  integer :: failures, ierr

  failures = 0
  call validate_snia_population_realization(realization, ierr)
  call expect(ierr == snia_population_contract_err_unapproved, &
       'unapproved population realization is rejected', failures)

  realization%approved = .true.
  realization%population_source_id = 'unit-binary-population'
  realization%population_model_id = population_binary_ssp
  realization%imf_id = stellar_imf_kroupa
  realization%binary_fraction = 0.5_stellar_dp
  realization%imf_conversion_factor = 2.0_stellar_dp
  realization%minimum_delay_gyr = 0.01_stellar_dp
  realization%maximum_delay_gyr = 10.0_stellar_dp
  realization%power_law_index = -1.0_stellar_dp
  realization%events_per_initial_msun = 1.0e-3_stellar_dp
  realization%event_realization_policy = snia_realization_expectation
  realization%metallicity_policy = snia_metallicity_factor_supplied
  realization%source_commit_binding = &
       '0123456789abcdef0123456789abcdef01234567'
  realization%approval_id = 'unit-test-only'
  call validate_snia_population_realization(realization, ierr)
  call expect(ierr == snia_population_contract_ok, &
       'complete approved realization contract is valid', failures)

  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.01_stellar_dp, 0.1_stellar_dp, 0.5_stellar_dp, expected_events, ierr)
  call expect(ierr == snia_population_contract_ok, &
       'interval expectation is evaluated from the realization contract', failures)
  call expect_close(expected_events, 1.0_stellar_dp / 30.0_stellar_dp, &
       'IMF conversion and metallicity factor are applied once', failures)
  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.01_stellar_dp, 0.1_stellar_dp, 0.5_stellar_dp, repeat_events, ierr)
  call expect(ierr == snia_population_contract_ok, &
       'repeated interval evaluation remains valid', failures)
  call expect_close(repeat_events, expected_events, &
       'interval expectation is restart-deterministic', failures)

  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.01_stellar_dp, 0.03_stellar_dp, 0.5_stellar_dp, first_interval, ierr)
  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.03_stellar_dp, 0.1_stellar_dp, 0.5_stellar_dp, second_interval, ierr)
  call expect_close(first_interval + second_interval, expected_events, &
       'adjacent population intervals telescope', failures)
  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.01_stellar_dp, 0.1_stellar_dp, 0.0_stellar_dp, repeat_events, ierr)
  call expect(ierr == snia_population_contract_ok .and. &
       repeat_events == 0.0_stellar_dp, &
       'zero metallicity factor yields zero events explicitly', failures)

  realization%event_realization_policy = snia_realization_poisson
  call evaluate_snia_interval_events(realization, 100.0_stellar_dp, &
       0.01_stellar_dp, 0.1_stellar_dp, 0.5_stellar_dp, repeat_events, ierr)
  call expect(ierr == snia_population_contract_err_realization .and. &
       repeat_events == 0.0_stellar_dp, &
       'Poisson realization is rejected without a seeded RNG contract', failures)
  realization%event_realization_policy = snia_realization_expectation

  realization%source_commit_binding = 'not-a-commit'
  call validate_snia_population_realization(realization, ierr)
  call expect(ierr == snia_population_contract_err_parameter, &
       'source commit binding must be an immutable hex revision', failures)
  realization%source_commit_binding = &
       '0123456789abcdef0123456789abcdef01234567'
  realization%maximum_delay_gyr = 0.005_stellar_dp
  call validate_snia_population_realization(realization, ierr)
  call expect(ierr == snia_population_contract_err_parameter, &
       'non-increasing delay support is rejected', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_POPULATION_CONTRACT_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_POPULATION_CONTRACT_TEST_OK'

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

end program fp2_snia_population_contract_test

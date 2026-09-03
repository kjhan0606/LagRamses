program fp2_snia_physical_contract_test
  use stellar_enrichment_config, only: stellar_dp
  use stellar_snia_physical_contract, only: snia_physical_contract_t, &
       snia_event_budget_t, snia_contract_ok, snia_contract_err_unapproved, &
       snia_contract_err_argument, &
       snia_contract_err_reservoir, snia_contract_err_mass, &
       snia_contract_err_momentum, snia_contract_err_yield, snia_wd_debit_per_event, &
       snia_shortfall_reject_event_interval, &
       snia_momentum_source_frame_vector, snia_momentum_isotropic_zero_vector, &
       snia_momentum_radial_magnitude, validate_snia_physical_contract, &
       resolve_snia_event_momentum, build_snia_event_budget
  implicit none

  type(snia_physical_contract_t) :: contract
  type(snia_event_budget_t) :: budget
  real(stellar_dp) :: momentum(3), bad_direction(3), unit_direction(3)
  integer :: failures, ierr

  failures = 0
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_unapproved, &
       'unapproved physical contract is rejected', failures)

  contract%approved = .true.
  contract%wd_debit_policy = snia_wd_debit_per_event
  contract%shortfall_policy = snia_shortfall_reject_event_interval
  contract%momentum_policy = snia_momentum_source_frame_vector
  contract%returned_mass_per_event = 1.30_stellar_dp
  contract%terminal_remnant_per_event = 0.10_stellar_dp
  contract%wd_debit_per_event = 1.40_stellar_dp
  contract%energy_per_event = 1.0e51_stellar_dp
  contract%momentum_per_event = (/1.0e40_stellar_dp, -2.0e39_stellar_dp, &
       3.0e39_stellar_dp/)
  contract%yield_source_id = 'unit-yield-source'
  contract%yield_source_sha256 = &
       '0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef'
  contract%source_commit_binding = '0123456789abcdef0123456789abcdef01234567'
  contract%conversion_code_sha256 = &
       'abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789'
  contract%approval_id = 'unit-test-only'
  contract%ejected_mass_per_event(1) = 1.0_stellar_dp
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_ok, 'explicit approved contract is valid', failures)

  contract%shortfall_policy = 0
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_argument, &
       'SNIa reservoir shortfall policy must be explicit', failures)
  contract%shortfall_policy = snia_shortfall_reject_event_interval

  call build_snia_event_budget(contract, 2.0_stellar_dp, 3.0_stellar_dp, &
       budget, ierr)
  call expect(ierr == snia_contract_ok, 'event budget fits WD reservoir', failures)
  call expect(abs(budget%wd_reservoir_debit - 2.8_stellar_dp) < 1.0e-12_stellar_dp, &
       'WD reservoir debit scales by event count', failures)
  call expect(abs(budget%returned_mass - 2.6_stellar_dp) < 1.0e-12_stellar_dp, &
       'returned mass scales by event count', failures)
  call expect(abs(budget%terminal_remnant_mass - 0.2_stellar_dp) < 1.0e-12_stellar_dp, &
       'terminal remnant scales by event count', failures)
  call expect(abs(budget%ejected_mass(1) - 2.0_stellar_dp) < 1.0e-12_stellar_dp, &
       'tracked SNIa ejecta scales by event count', failures)
  call expect(all(abs(budget%momentum - 2.0_stellar_dp * &
       contract%momentum_per_event) < 1.0e28_stellar_dp), &
       'signed source-frame momentum is preserved', failures)

  call build_snia_event_budget(contract, 3.0_stellar_dp, 3.0_stellar_dp, &
       budget, ierr)
  call expect(ierr == snia_contract_err_reservoir .and. &
       all(budget%momentum == 0.0_stellar_dp), &
       'WD reservoir overdraw is rejected transactionally', failures)

  contract%returned_mass_per_event = 1.5_stellar_dp
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_mass, &
       'event mass exceeding WD debit is rejected', failures)
  contract%returned_mass_per_event = 1.30_stellar_dp

  contract%terminal_remnant_per_event = 0.0_stellar_dp
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_mass, &
       'WD mass deficit is rejected rather than assigned to living stars', failures)
  contract%terminal_remnant_per_event = 0.10_stellar_dp

  contract%momentum_policy = snia_momentum_isotropic_zero_vector
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_momentum, &
       'isotropic policy rejects nonzero net vector momentum', failures)
  contract%momentum_per_event = 0.0_stellar_dp
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_ok, 'isotropic zero-vector policy is explicit', failures)

  contract%momentum_policy = snia_momentum_radial_magnitude
  contract%radial_momentum_per_event = 5.0e40_stellar_dp
  bad_direction = (/0.0_stellar_dp, 0.0_stellar_dp, 2.0_stellar_dp/)
  unit_direction = (/0.0_stellar_dp, 0.0_stellar_dp, 1.0_stellar_dp/)
  call resolve_snia_event_momentum(contract, bad_direction, momentum, ierr)
  call expect(ierr == snia_contract_err_momentum .and. all(momentum == 0.0_stellar_dp), &
       'radial policy rejects a non-unit deposition direction', failures)
  call build_snia_event_budget(contract, 1.0_stellar_dp, 2.0_stellar_dp, &
       budget, ierr)
  call expect(ierr == snia_contract_err_momentum .and. all(budget%momentum == 0.0_stellar_dp), &
       'radial budget rejects a missing deposition direction', failures)
  call build_snia_event_budget(contract, 2.0_stellar_dp, 3.0_stellar_dp, &
       budget, ierr, unit_direction)
  call expect(ierr == snia_contract_ok, 'radial budget accepts an explicit unit direction', failures)
  call expect(all(abs(budget%momentum - (/0.0_stellar_dp, 0.0_stellar_dp, &
       1.0e41_stellar_dp/)) < 1.0e29_stellar_dp), &
       'radial momentum maps into the declared cell direction', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_PHYSICAL_CONTRACT_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_PHYSICAL_CONTRACT_TEST_OK'

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

end program fp2_snia_physical_contract_test

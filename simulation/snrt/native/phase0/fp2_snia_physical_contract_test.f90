program fp2_snia_physical_contract_test
  use stellar_enrichment_config, only: stellar_dp
  use stellar_snia_physical_contract, only: snia_physical_contract_t, &
       snia_event_budget_t, snia_contract_ok, snia_contract_err_unapproved, &
       snia_contract_err_reservoir, snia_contract_err_mass, &
       snia_contract_err_momentum, snia_wd_debit_per_event, &
       snia_momentum_source_frame_vector, snia_momentum_isotropic_zero_vector, &
       validate_snia_physical_contract, build_snia_event_budget
  implicit none

  type(snia_physical_contract_t) :: contract
  type(snia_event_budget_t) :: budget
  integer :: failures, ierr

  failures = 0
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_unapproved, &
       'unapproved physical contract is rejected', failures)

  contract%approved = .true.
  contract%wd_debit_policy = snia_wd_debit_per_event
  contract%momentum_policy = snia_momentum_source_frame_vector
  contract%returned_mass_per_event = 1.30_stellar_dp
  contract%terminal_remnant_per_event = 0.0_stellar_dp
  contract%wd_debit_per_event = 1.40_stellar_dp
  contract%energy_per_event = 1.0e51_stellar_dp
  contract%momentum_per_event = (/1.0e40_stellar_dp, -2.0e39_stellar_dp, &
       3.0e39_stellar_dp/)
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_ok, 'explicit approved contract is valid', failures)

  call build_snia_event_budget(contract, 2.0_stellar_dp, 3.0_stellar_dp, &
       budget, ierr)
  call expect(ierr == snia_contract_ok, 'event budget fits WD reservoir', failures)
  call expect(abs(budget%wd_reservoir_debit - 2.8_stellar_dp) < 1.0e-12_stellar_dp, &
       'WD reservoir debit scales by event count', failures)
  call expect(abs(budget%returned_mass - 2.6_stellar_dp) < 1.0e-12_stellar_dp, &
       'returned mass scales by event count', failures)
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

  contract%momentum_policy = snia_momentum_isotropic_zero_vector
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_err_momentum, &
       'isotropic policy rejects nonzero net vector momentum', failures)
  contract%momentum_per_event = 0.0_stellar_dp
  call validate_snia_physical_contract(contract, ierr)
  call expect(ierr == snia_contract_ok, 'isotropic zero-vector policy is explicit', failures)

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

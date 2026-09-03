program fp2_snia_event_ledger_test
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       channel_snia
  use stellar_enrichment_contract, only: stellar_source_t
  use stellar_snia_event_ledger, only: snia_event_ledger_ok, &
       snia_event_ledger_err_argument, snia_event_ledger_err_mass, &
       build_snia_event_source
  implicit none

  type(stellar_source_t) :: source
  real(stellar_dp) :: ejecta(n_stellar_elements), net_yield(n_stellar_elements)
  real(stellar_dp) :: momentum(3), expected_events, returned_mass, energy
  real(stellar_dp) :: tracked, scale
  integer :: failures, ierr

  failures = 0
  ejecta = 0.0_stellar_dp
  net_yield = 0.0_stellar_dp
  ejecta(1) = 0.55_stellar_dp
  ejecta(2) = 0.20_stellar_dp
  ejecta(5) = 0.10_stellar_dp
  ejecta(11) = 0.05_stellar_dp
  net_yield(5) = 0.08_stellar_dp
  momentum = (/1.0e40_stellar_dp, -2.0e39_stellar_dp, 3.0e39_stellar_dp/)
  expected_events = 2.5_stellar_dp
  returned_mass = 1.2_stellar_dp
  energy = 1.0e51_stellar_dp

  call build_snia_event_source(expected_events, returned_mass, ejecta, &
       net_yield, energy, momentum, source, ierr)
  call expect(ierr == snia_event_ledger_ok, 'valid event source accepted', failures)
  call expect(abs(source%returned_mass - expected_events * returned_mass) < 1.0e-12_stellar_dp, &
       'returned mass scales by expected event count', failures)
  call expect(abs(sum(source%ejected_mass) - expected_events * sum(ejecta)) < 1.0e-12_stellar_dp, &
       'tracked ejecta scales by expected event count', failures)
  call expect(abs(source%energy - expected_events * energy) < 1.0e36_stellar_dp, &
       'energy scales by expected event count', failures)
  call expect(all(abs(source%momentum - expected_events * momentum) < 1.0e28_stellar_dp), &
       'signed source-frame momentum scales by expected event count', failures)
  call expect(source%channel_returned_mass(channel_snia) == source%returned_mass, &
       'SNIa owns its returned source channel', failures)
  call expect(sum(source%channel_returned_mass) == source%returned_mass, &
       'event source has no duplicate channel return', failures)
  call expect(sum(source%channel_ejected_mass(channel_snia,:)) <= &
       source%channel_returned_mass(channel_snia) + 1.0e-12_stellar_dp, &
       'tracked ejecta remains below returned mass', failures)

  tracked = sum(ejecta)
  scale = max(1.0_stellar_dp, returned_mass, tracked)
  call expect(abs(source%returned_mass - sum(source%ejected_mass) - &
       expected_events * (returned_mass - tracked)) < 1.0e-12_stellar_dp * expected_events * scale, &
       'untracked event residual is preserved', failures)

  call build_snia_event_source(expected_events, returned_mass, ejecta, &
       net_yield, energy, momentum, source, ierr)
  call expect(ierr == snia_event_ledger_ok, 'repeat event interval is deterministic', failures)

  call build_snia_event_source(-1.0_stellar_dp, returned_mass, ejecta, &
       net_yield, energy, momentum, source, ierr)
  call expect(ierr == snia_event_ledger_err_argument .and. source%returned_mass == 0.0_stellar_dp, &
       'negative event count is rejected transactionally', failures)

  ejecta(1) = returned_mass + 1.0_stellar_dp
  call build_snia_event_source(expected_events, returned_mass, ejecta, &
       net_yield, energy, momentum, source, ierr)
  call expect(ierr == snia_event_ledger_err_mass .and. source%returned_mass == 0.0_stellar_dp, &
       'tracked over-return is rejected transactionally', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_EVENT_LEDGER_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_EVENT_LEDGER_TEST_OK'

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

end program fp2_snia_event_ledger_test

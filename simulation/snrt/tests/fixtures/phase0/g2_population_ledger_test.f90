program g2_population_ledger_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp, n_stellar_channels, &
       n_stellar_elements, &
       set_enrichment_defaults, enable_wind, enable_agb, enable_snii, &
       enable_snia, enable_pisn, yield_basis_per_star_cumulative, &
       yield_basis_ssp_cumulative
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, stellar_source_t, untracked_ejecta_mass, &
       generic_metal_ejecta_mass
  use stellar_population_ledger, only: stellar_population_ledger_t, &
       finalize_population_ledger, &
       set_white_dwarf_reservoir, apply_snia_event_budget, &
       compute_unresolved_mass_bucket, &
       population_ledger_ok, population_ledger_err_mass, &
       population_ledger_err_owner, population_ledger_err_nonfinite, &
       population_ledger_err_argument, population_ledger_err_snia
  use stellar_yield_tables, only: stellar_yield_table_t, clear_yield_table
  use stellar_enrichment_driver, only: compute_stellar_source_increment, &
       compute_stellar_cumulative, enrichment_driver_ok, &
       enrichment_driver_err_unsupported, &
       enrichment_driver_err_ledger
  use stellar_ssp_sources, only: integrate_ssp_channel, calculate_imf_normalization, &
       evaluate_imf, ssp_source_err_basis
  use stellar_snia_physical_contract, only: snia_physical_contract_t, &
       snia_event_budget_t, snia_contract_ok, snia_wd_debit_per_event, &
       snia_shortfall_reject_event_interval, &
       snia_momentum_source_frame_vector, build_snia_event_budget
  implicit none

  type(stellar_population_t) :: population
  type(stellar_cumulative_t) :: states(n_stellar_channels)
  type(stellar_population_ledger_t) :: ledger
  type(stellar_cumulative_t) :: driver_states(n_stellar_channels)
  type(stellar_population_ledger_t) :: driver_ledger
  type(stellar_yield_table_t) :: mini_table
  type(stellar_yield_table_t) :: unloaded_table
  logical :: enabled(n_stellar_channels), owners(n_stellar_channels)
  real(stellar_dp) :: channel_mass_min(n_stellar_channels)
  real(stellar_dp) :: channel_mass_max(n_stellar_channels)
  real(stellar_dp) :: chemistry_probe(n_stellar_elements)
  real(stellar_dp) :: normalization, normalization_truncated, mass_integral
  real(stellar_dp) :: unresolved_fraction, unresolved_bucket
  real(stellar_dp) :: overlapping_min(2), overlapping_max(2)
  integer :: ierr, failures, channel, imf, row
  type(stellar_source_t) :: source
  type(snia_physical_contract_t) :: snia_contract
  type(snia_event_budget_t) :: snia_budget

  failures = 0
  call set_enrichment_defaults()
  population%formation_time = 0.0_stellar_dp
  population%initial_mass = 100.0_stellar_dp
  population%current_mass = 100.0_stellar_dp
  population%birth_metallicity = 1.0e-3_stellar_dp
  population%birth_mass_fraction = 0.0_stellar_dp
  population%imf_id = 1
  population%imf_mass_min = 0.08_stellar_dp
  population%imf_mass_max = 120.0_stellar_dp
  population%population_id = 0
  population%binary_fraction = 0.0_stellar_dp
  population%yield_basis_id = yield_basis_per_star_cumulative
  population%pisn_enabled = .false.

  do channel = 1, n_stellar_channels
     states(channel)%ejected_mass = 0.0_stellar_dp
     states(channel)%net_yield = 0.0_stellar_dp
     states(channel)%returned_mass = 0.0_stellar_dp
     states(channel)%remnant_mass = 0.0_stellar_dp
     states(channel)%living_mass = 0.0_stellar_dp
     states(channel)%energy = 0.0_stellar_dp
     states(channel)%momentum = 0.0_stellar_dp
     states(channel)%channel_returned_mass = 0.0_stellar_dp
     states(channel)%channel_energy = 0.0_stellar_dp
     states(channel)%channel_momentum = 0.0_stellar_dp
     states(channel)%channel_ejected_mass = 0.0_stellar_dp
     states(channel)%channel_net_yield = 0.0_stellar_dp
  end do
  enabled = (/ .true., .true., .true., .false., .false. /)
  owners = (/ .false., .true., .true., .false., .true. /)

  call set_state(states(1), 1, 5.0_stellar_dp, 0.0_stellar_dp)
  call set_state(states(2), 2, 15.0_stellar_dp, 2.0_stellar_dp)
  call set_state(states(3), 3, 20.0_stellar_dp, 5.0_stellar_dp)
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(ierr == population_ledger_ok, &
       'channel states close a population mass ledger', failures)
  call expect_close(ledger%returned_mass, 40.0_stellar_dp, &
       'returned mass is summed once by channel', failures)
  call expect_close(ledger%remnant_mass, 7.0_stellar_dp, &
       'terminal remnant mass is retained by owner channels', failures)
  call expect_close(ledger%living_mass, 53.0_stellar_dp, &
       'living mass is derived after channel aggregation', failures)
  call expect_close(ledger%initial_mass - ledger%living_mass - &
       ledger%remnant_mass - ledger%returned_mass, 0.0_stellar_dp, &
       'population mass closure residual is zero', failures)
  call expect_close(ledger%untracked_ejecta_mass, 0.0_stellar_dp, &
       'fully tracked synthetic states have zero untracked residual', failures)
  call set_white_dwarf_reservoir(ledger, 2.0_stellar_dp, 1.0e-12_stellar_dp, ierr)
  call expect(ierr == population_ledger_ok, &
       'explicit WD reservoir is admitted as a remnant subset', failures)
  call set_white_dwarf_reservoir(ledger, 6.0_stellar_dp, 1.0e-12_stellar_dp, ierr)
  call expect(ierr == population_ledger_err_argument .or. &
       ierr == population_ledger_err_snia, &
       'SNII remnant mass cannot fund a WD reservoir', failures)
  call set_white_dwarf_reservoir(ledger, 2.0_stellar_dp, 1.0e-12_stellar_dp, ierr)
  snia_contract%approved = .true.
  snia_contract%wd_debit_policy = snia_wd_debit_per_event
  snia_contract%shortfall_policy = snia_shortfall_reject_event_interval
  snia_contract%momentum_policy = snia_momentum_source_frame_vector
  snia_contract%yield_source_id = 'unit-yield-source'
  snia_contract%yield_source_sha256 = &
       '0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef'
  snia_contract%source_commit_binding = '0123456789abcdef0123456789abcdef01234567'
  snia_contract%conversion_code_sha256 = &
       'abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789'
  snia_contract%approval_id = 'unit-test-only'
  snia_contract%returned_mass_per_event = 1.3_stellar_dp
  snia_contract%terminal_remnant_per_event = 0.1_stellar_dp
  snia_contract%wd_debit_per_event = 1.4_stellar_dp
  snia_contract%energy_per_event = 1.0e51_stellar_dp
  snia_contract%ejected_mass_per_event(1) = 1.0_stellar_dp
  snia_contract%momentum_per_event = 0.0_stellar_dp
  call build_snia_event_budget(snia_contract, 1.0_stellar_dp, &
       ledger%wd_reservoir_mass, snia_budget, ierr)
  call expect(ierr == snia_contract_ok, 'SNIa budget is built from explicit WD reservoir', failures)
  call apply_snia_event_budget(ledger, snia_budget, 1.0e-12_stellar_dp, ierr)
  call expect(ierr == population_ledger_ok, 'SNIa budget updates population ledger transactionally', failures)
  call expect_close(ledger%wd_reservoir_mass, 0.6_stellar_dp, &
       'WD reservoir is debited exactly once', failures)
  call expect_close(ledger%snia_wd_debit_mass, 1.4_stellar_dp, &
       'SNIa WD debit is retained in the ledger', failures)
  call expect_close(ledger%returned_mass, 41.3_stellar_dp, &
       'SNIa returned mass enters aggregate return', failures)
  call expect_close(ledger%remnant_mass, 5.7_stellar_dp, &
       'consumed WD remnant leaves aggregate remnant mass', failures)
  call expect_close(ledger%living_mass, 53.0_stellar_dp, &
       'living mass is recomputed after SNIa debit', failures)
  call expect_close(ledger%channel_remnant_mass(2), 0.6_stellar_dp, &
       'WD debit leaves the owning AGB remnant channel', failures)
  call expect_close(ledger%channel_remnant_mass(4), 0.1_stellar_dp, &
       'terminal SNIa remnant is recorded in its event channel', failures)
  call expect_close(ledger%tracked_ejecta_mass, 41.0_stellar_dp, &
       'SNIa tracked ejecta updates the aggregate ledger', failures)
  call expect_close(ledger%untracked_ejecta_mass, 0.3_stellar_dp, &
       'SNIa untracked returned residual updates the aggregate ledger', failures)
  call expect_close(ledger%channel_tracked_ejecta_mass(4), 1.0_stellar_dp, &
       'SNIa tracked ejecta updates its channel ledger', failures)
  call expect_close(ledger%channel_untracked_ejecta_mass(4), 0.3_stellar_dp, &
       'SNIa untracked residual updates its channel ledger', failures)
  call expect_close(ledger%initial_mass - ledger%living_mass - ledger%remnant_mass - &
       ledger%returned_mass, 0.0_stellar_dp, &
       'SNIa ledger update preserves population mass closure', failures)
  call expect(ledger%unresolved_initial_mass_fraction > 0.1178_stellar_dp .and. &
       ledger%unresolved_initial_mass_fraction < 0.1180_stellar_dp, &
       'unresolved fate bucket is computed from the configured IMF', failures)
  call expect_close(ledger%unresolved_initial_mass, &
       population%initial_mass * ledger%unresolved_initial_mass_fraction, &
       'unresolved fate bucket remains an initial-mass diagnostic', failures)
  overlapping_min = (/0.8_stellar_dp, 0.9_stellar_dp/)
  overlapping_max = (/1.0_stellar_dp, 1.1_stellar_dp/)
  call compute_unresolved_mass_bucket(population, overlapping_min, &
       overlapping_max, 2, unresolved_fraction, unresolved_bucket, ierr)
  call expect(ierr == population_ledger_err_argument, &
       'overlapping unresolved intervals are rejected', failures)
  call compute_unresolved_mass_bucket(population, (/0.8d0, 40.0d0/), &
       (/1.0d0, 120.0d0/), 2, unresolved_fraction, unresolved_bucket, ierr)
  call expect(ierr == population_ledger_ok, &
       'explicit unresolved interval bucket succeeds', failures)

  chemistry_probe = 0.0_stellar_dp
  chemistry_probe(1) = 4.0_stellar_dp
  chemistry_probe(2) = 3.0_stellar_dp
  chemistry_probe(3) = 1.0_stellar_dp
  chemistry_probe(n_stellar_elements) = 0.5_stellar_dp
  call expect_close(untracked_ejecta_mass(10.0_stellar_dp, chemistry_probe), &
       1.5_stellar_dp, 'untracked ejecta are derived without a new table field', &
       failures)
  call expect_close(generic_metal_ejecta_mass(10.0_stellar_dp, &
       chemistry_probe), 3.0_stellar_dp, &
       'generic metal source includes tracked metals plus residual', failures)

  call make_mini_table(mini_table)
  channel_mass_min = (/1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
  channel_mass_max = (/2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0/)
  call compute_stellar_cumulative(mini_table, population, 0.5_stellar_dp, &
       channel_mass_min, channel_mass_max, 16, driver_states, driver_ledger, &
       ierr)
  call expect(ierr == 0, 'SSP driver returns a closed population ledger', failures)
  call expect(driver_states(1)%remnant_mass == 0.0_stellar_dp, &
       'wind channel cannot contribute a terminal remnant', failures)
  call expect(driver_states(2)%remnant_mass > 0.0_stellar_dp .and. &
       driver_states(3)%remnant_mass > 0.0_stellar_dp, &
       'SSP driver retains remnant contributions from terminal channels', failures)
  call expect(driver_ledger%remnant_mass > 0.0_stellar_dp .and. &
       driver_ledger%living_mass >= 0.0_stellar_dp, &
       'SSP driver ledger derives nonnegative living mass', failures)

  call compute_stellar_source_increment(mini_table, population, &
       0.0_stellar_dp, 0.5_stellar_dp, channel_mass_min, channel_mass_max, &
       16, source, ierr, driver_ledger)
  call expect(ierr == 0, &
       'timestep source path executes the cumulative population ledger', failures)
  call expect_close(source%returned_mass, driver_ledger%returned_mass, &
       'age-zero anchored timestep return matches cumulative ledger', failures)

  do row = 1, mini_table%n_rows
     if (mini_table%channel(row) == 1 .and. &
          mini_table%age_gyr(row) > 0.0_stellar_dp) then
        mini_table%remnant_mass(row) = 0.01_stellar_dp
     end if
  end do
  call compute_stellar_source_increment(mini_table, population, &
       0.0_stellar_dp, 0.5_stellar_dp, channel_mass_min, channel_mass_max, &
       16, source, ierr)
  call expect(ierr == enrichment_driver_err_ledger, &
       'timestep source cannot bypass non-owner remnant rejection', failures)
  call expect_close(source%returned_mass, 0.0_stellar_dp, &
       'ledger rejection clears the timestep source transaction', failures)
  do row = 1, mini_table%n_rows
     if (mini_table%channel(row) == 1) then
        mini_table%remnant_mass(row) = 0.0_stellar_dp
     end if
  end do

  do imf = 0, 4
     call calculate_imf_normalization(imf, population%imf_mass_min, &
          population%imf_mass_max, normalization, ierr)
     call expect(ierr == 0 .and. normalization > 0.0_stellar_dp, &
          'configured IMF support has a positive normalization', failures)
     call integrate_imf_mass(imf, population%imf_mass_min, &
          population%imf_mass_max, normalization, mass_integral)
     call expect(abs(mass_integral - 1.0_stellar_dp) < 1.0e-5_stellar_dp, &
          'IMF normalization closes one unit of initial stellar mass', failures)
  end do
  call calculate_imf_normalization(population%imf_id, 1.0_stellar_dp, &
       population%imf_mass_max, normalization_truncated, ierr)
  call calculate_imf_normalization(population%imf_id, population%imf_mass_min, &
       population%imf_mass_max, normalization, ierr)
  call expect(abs(normalization_truncated-normalization) > 1.0e-6_stellar_dp, &
       'changing configured IMF support changes its normalization', failures)
  call expect(abs(evaluate_imf(1.0_stellar_dp,2) - &
       exp(-(log10(1.0_stellar_dp)-log10(0.079_stellar_dp))**2 / &
       (2.0_stellar_dp*0.69_stellar_dp**2))) < 1.0e-14_stellar_dp, &
       'Chabrier IMF is continuous at the one-solar-mass branch', failures)

  population%yield_basis_id = yield_basis_ssp_cumulative
  call integrate_ssp_channel(mini_table, population, 1, 0.5_stellar_dp, &
       1.0_stellar_dp, 2.0_stellar_dp, 16, driver_states(1), ierr)
  call expect(ierr == ssp_source_err_basis, &
       'SSP-normalized source is rejected before a second IMF convolution', failures)
  population%yield_basis_id = yield_basis_per_star_cumulative

  states(1)%remnant_mass = 1.0_stellar_dp
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(iand(ierr, population_ledger_err_owner) /= 0, &
       'non-terminal channel remnant is rejected', failures)
  states(1)%remnant_mass = 0.0_stellar_dp

  states(3)%channel_ejected_mass(3,1) = 19.0_stellar_dp
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(iand(ierr, population_ledger_err_mass) /= 0, &
       'channel ejecta/returned mismatch is rejected', failures)
  states(3)%ejected_mass(1) = 19.0_stellar_dp
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(ierr == population_ledger_ok, &
       'matched tracked ejecta below returned mass are accepted', failures)
  call expect_close(ledger%channel_untracked_ejecta_mass(3), &
       1.0_stellar_dp, 'channel ledger preserves untracked ejecta residual', &
       failures)
  call expect_close(ledger%untracked_ejecta_mass, 1.0_stellar_dp, &
       'population ledger preserves untracked ejecta residual', failures)
  states(3)%ejected_mass(1) = 20.0_stellar_dp
  states(3)%channel_ejected_mass(3,1) = 20.0_stellar_dp

  states(2)%returned_mass = ieee_value(0.0_stellar_dp, ieee_quiet_nan)
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(iand(ierr, population_ledger_err_nonfinite) /= 0, &
       'nonfinite channel ledger input is rejected', failures)
  states(2)%returned_mass = 15.0_stellar_dp

  ! A widened channel-3 candidate window must not be integrated until the
  ! source-node fate consumer exists, even if a caller reaches the driver.
  enable_wind = .false.
  enable_agb = .false.
  enable_snii = .true.
  enable_snia = .false.
  enable_pisn = .false.
  channel_mass_min = (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
  channel_mass_max = (/120.0d0, 8.0d0, 120.0d0, 8.0d0, 260.0d0/)
  call compute_stellar_source_increment(unloaded_table, population, &
       0.0_stellar_dp, 1.0_stellar_dp, channel_mass_min, channel_mass_max, &
       8, source, ierr)
  call expect(ierr == enrichment_driver_err_unsupported, &
       'high-mass SNII window is refused without source-node fate consumer', failures)
  call compute_stellar_cumulative(unloaded_table, population, 1.0_stellar_dp, &
       channel_mass_min, channel_mass_max, 8, driver_states, driver_ledger, ierr)
  call expect(ierr == enrichment_driver_err_unsupported, &
       'cumulative high-mass SNII is refused without source-node fate consumer', failures)

  ! SNIa is intentionally absent from the generic IMF-only driver: the
  ! runtime DTD caller owns its interval convolution and event ledger.
  enable_wind = .false.
  enable_agb = .false.
  enable_snii = .false.
  enable_snia = .true.
  channel_mass_min = (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
  channel_mass_max = (/120.0d0, 8.0d0, 120.0d0, 8.0d0, 260.0d0/)
  call compute_stellar_source_increment(unloaded_table, population, &
       0.0_stellar_dp, 1.0_stellar_dp, channel_mass_min, channel_mass_max, &
       8, source, ierr)
  call expect(ierr == enrichment_driver_ok, &
       'generic driver leaves SNIa interval evaluation to the runtime DTD caller', failures)
  call expect_close(source%returned_mass, 0.0_stellar_dp, &
       'generic SNIa path emits no prompt-table return', failures)

  ! PISN still cannot enter this driver without its explicit fate resolver.
  enable_snia = .false.
  enable_pisn = .true.
  call compute_stellar_source_increment(unloaded_table, population, &
       0.0_stellar_dp, 1.0_stellar_dp, channel_mass_min, channel_mass_max, &
       8, source, ierr)
  call expect(ierr == enrichment_driver_err_unsupported, &
       'PISN activation is refused without the fate resolver', failures)
  call set_enrichment_defaults()
  call clear_yield_table(mini_table)

  if (failures == 0) then
     write(*, '(a)') 'G2_POPULATION_LEDGER_TEST_OK'
  else
     write(*, '(a,i0)') 'G2_POPULATION_LEDGER_TEST_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine integrate_imf_mass(imf_id, mass_min, mass_max, normalization, integral)
    integer, intent(in) :: imf_id
    real(stellar_dp), intent(in) :: mass_min, mass_max, normalization
    real(stellar_dp), intent(out) :: integral
    integer, parameter :: n_bins = 131072
    integer :: bin
    real(stellar_dp) :: dlog_mass, log_mass, mass, dm

    dlog_mass = (log10(mass_max)-log10(mass_min))/real(n_bins,stellar_dp)
    integral = 0.0_stellar_dp
    do bin = 1, n_bins
       log_mass = log10(mass_min) + &
            (real(bin,stellar_dp)-0.5_stellar_dp)*dlog_mass
       mass = 10.0_stellar_dp**log_mass
       dm = mass*log(10.0_stellar_dp)*dlog_mass
       integral = integral + mass*normalization*evaluate_imf(mass,imf_id)*dm
    end do
  end subroutine integrate_imf_mass

  subroutine make_mini_table(table)
    type(stellar_yield_table_t), intent(out) :: table
    integer :: row, channel, im, iz, ia
    real(stellar_dp), parameter :: masses(2) = (/1.0d0, 2.0d0/)
    real(stellar_dp), parameter :: metallicities(2) = (/1.0d-3, 1.0d-2/)
    real(stellar_dp), parameter :: ages(2) = (/0.0d0, 1.0d0/)
    real(stellar_dp) :: returned, age

    table%n_rows = 3 * 2 * 2 * 2
    allocate(table%channel(table%n_rows), table%initial_mass(table%n_rows), &
         table%birth_metallicity(table%n_rows), table%age_gyr(table%n_rows), &
         table%returned_mass(table%n_rows), table%remnant_mass(table%n_rows), &
         table%energy(table%n_rows), table%momentum(table%n_rows,3), &
         table%ejected_mass(table%n_rows,n_stellar_elements), &
         table%net_yield(table%n_rows,n_stellar_elements))
    row = 0
    do channel = 1, 3
       do im = 1, 2
          do iz = 1, 2
             do ia = 1, 2
                row = row + 1
                age = ages(ia)
                returned = age * (0.1_stellar_dp * real(channel, stellar_dp) + &
                     0.05_stellar_dp * masses(im) + &
                     0.2_stellar_dp * metallicities(iz))
                table%channel(row) = channel
                table%initial_mass(row) = masses(im)
                table%birth_metallicity(row) = metallicities(iz)
                table%age_gyr(row) = age
                table%returned_mass(row) = returned
                if (channel == 1) then
                   table%remnant_mass(row) = 0.0_stellar_dp
                else
                   table%remnant_mass(row) = age * &
                        0.05_stellar_dp * real(channel, stellar_dp)
                end if
                table%energy(row) = age * real(channel, stellar_dp) * 1.0e48_stellar_dp
                table%momentum(row,:) = 0.0_stellar_dp
                table%ejected_mass(row,:) = 0.0_stellar_dp
                table%ejected_mass(row,1) = returned
                table%net_yield(row,:) = 0.0_stellar_dp
             end do
          end do
       end do
    end do
    table%loaded = .true.
  end subroutine make_mini_table

  subroutine set_state(state, channel, returned_mass, remnant_mass)
    type(stellar_cumulative_t), intent(inout) :: state
    integer, intent(in) :: channel
    real(stellar_dp), intent(in) :: returned_mass, remnant_mass

    state%returned_mass = returned_mass
    state%remnant_mass = remnant_mass
    state%ejected_mass(1) = returned_mass
    state%channel_returned_mass = 0.0_stellar_dp
    state%channel_returned_mass(channel) = returned_mass
    state%channel_ejected_mass = 0.0_stellar_dp
    state%channel_ejected_mass(channel,1) = returned_mass
  end subroutine set_state

  subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (condition) then
       write(*, '(a)') 'PASS: ' // trim(label)
    else
       failures = failures + 1
       write(*, '(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine expect

  subroutine expect_close(value, expected, label, failures)
    real(stellar_dp), intent(in) :: value, expected
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    call expect(abs(value - expected) <= 1.0e-12_stellar_dp, label, failures)
  end subroutine expect_close

end program g2_population_ledger_test

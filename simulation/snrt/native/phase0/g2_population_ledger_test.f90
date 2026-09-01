program g2_population_ledger_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp, n_stellar_channels, &
       n_stellar_elements, &
       set_enrichment_defaults, enable_wind, enable_agb, enable_snii, &
       enable_snia
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, stellar_source_t, untracked_ejecta_mass, &
       generic_metal_ejecta_mass
  use stellar_population_ledger, only: stellar_population_ledger_t, &
       finalize_population_ledger, &
       population_ledger_ok, population_ledger_err_mass, &
       population_ledger_err_owner, population_ledger_err_nonfinite
  use stellar_yield_tables, only: stellar_yield_table_t, clear_yield_table
  use stellar_enrichment_driver, only: compute_stellar_source_increment, &
       compute_stellar_cumulative, enrichment_driver_err_unsupported
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
  integer :: ierr, failures, channel
  type(stellar_source_t) :: source

  failures = 0
  call set_enrichment_defaults()
  population%formation_time = 0.0_stellar_dp
  population%initial_mass = 100.0_stellar_dp
  population%current_mass = 100.0_stellar_dp
  population%birth_metallicity = 1.0e-3_stellar_dp
  population%birth_mass_fraction = 0.0_stellar_dp
  population%imf_id = 1
  population%population_id = 0
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
  call set_state(states(2), 2, 15.0_stellar_dp, 1.0_stellar_dp)
  call set_state(states(3), 3, 20.0_stellar_dp, 5.0_stellar_dp)
  call finalize_population_ledger(population, states, enabled, owners, &
       1.0e-12_stellar_dp, ledger, ierr)
  call expect(ierr == population_ledger_ok, &
       'channel states close a population mass ledger', failures)
  call expect_close(ledger%returned_mass, 40.0_stellar_dp, &
       'returned mass is summed once by channel', failures)
  call expect_close(ledger%remnant_mass, 6.0_stellar_dp, &
       'terminal remnant mass is retained by owner channels', failures)
  call expect_close(ledger%living_mass, 54.0_stellar_dp, &
       'living mass is derived after channel aggregation', failures)
  call expect_close(ledger%initial_mass - ledger%living_mass - &
       ledger%remnant_mass - ledger%returned_mass, 0.0_stellar_dp, &
       'population mass closure residual is zero', failures)
  call expect_close(ledger%untracked_ejecta_mass, 0.0_stellar_dp, &
       'fully tracked synthetic states have zero untracked residual', failures)

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

  ! An enabled SNIa/PISN switch cannot fall through to ordinary IMF integration.
  enable_wind = .false.
  enable_agb = .false.
  enable_snii = .false.
  enable_snia = .true.
  channel_mass_min = (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
  channel_mass_max = (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
  call compute_stellar_source_increment(unloaded_table, population, &
       0.0_stellar_dp, 1.0_stellar_dp, channel_mass_min, channel_mass_max, &
       8, source, ierr)
  call expect(ierr == enrichment_driver_err_unsupported, &
       'SNIa activation is refused before DTD implementation', failures)
  call set_enrichment_defaults()
  call clear_yield_table(mini_table)

  if (failures == 0) then
     write(*, '(a)') 'G2_POPULATION_LEDGER_TEST_OK'
  else
     write(*, '(a,i0)') 'G2_POPULATION_LEDGER_TEST_FAIL count=', failures
     error stop 1
  end if

contains

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

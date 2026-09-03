program g1_contract_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       channel_snii
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_source_t
  use stellar_yield_tables, only: stellar_yield_table_t, load_yield_table, &
       clear_yield_table, yield_table_ok, set_yield_mass_assignment_mode, &
       yield_mass_assignment_linear, yield_mass_assignment_piecewise_constant, &
       yield_table_err_assignment_mode
  use stellar_yield_interpolation, only: interpolate_yield_row, &
       interpolation_ok, interpolation_err_grid, interpolation_err_assignment_mode
  use stellar_yield_audit, only: audit_yield_table, yield_audit_ok, &
       yield_audit_err_grid, yield_audit_err_duplicate, &
       yield_audit_err_nonfinite, yield_audit_err_mass, &
       yield_audit_err_monotonic
  use stellar_source_increment, only: integrate_ssp_channel_increment, &
       source_increment_ok, source_increment_err_negative
  use stellar_ssp_sources, only: imf_kroupa
  use stellar_native_units, only: code_time_to_age_gyr, &
       code_interval_to_age_gyr, mass_code_to_msun, mass_msun_to_code, &
       energy_erg_to_code, units_ok
  use stellar_ramses_field_map, only: stellar_field_map_t, clear_field_map, &
       validate_field_map
  use stellar_progress_contract, only: stellar_progress_t, &
       progress_initialize, progress_begin, progress_commit, &
       progress_abort, progress_export, progress_ok
  implicit none

  type(stellar_yield_table_t) :: table, incomplete_table, duplicate_table
  type(stellar_yield_table_t) :: nan_table, negative_table
  type(stellar_yield_table_t) :: residual_table, overfull_table
  type(stellar_yield_table_t) :: nonmonotonic_residual_table
  type(stellar_population_t) :: population
  type(stellar_source_t) :: source_a, source_b, source_full
  type(stellar_source_t) :: negative_source
  type(stellar_field_map_t) :: field_map
  type(stellar_progress_t) :: progress, restarted
  real(stellar_dp) :: returned_mass, remnant_mass, energy
  real(stellar_dp) :: momentum(3), ejected_mass(n_stellar_elements)
  real(stellar_dp) :: net_yield(n_stellar_elements)
  real(stellar_dp) :: query_mass, query_z, query_age
  real(stellar_dp) :: expected, dt_gyr, mass_msun, mass_code
  real(stellar_dp) :: native_age, native_dt, native_energy
  real(stellar_dp) :: committed_age
  real(stellar_dp) :: source_tolerance
  integer :: ierr, audit_ierr, failures, i, unit, env_status
  integer :: source_ierr, progress_ierr, units_ierr, map_ierr
  logical :: should_deposit
  character(len=512) :: table_filename, result_filename
  integer, parameter :: n_differential_queries = 6
  real(stellar_dp), parameter :: differential_masses(n_differential_queries) = &
       (/1.0d0, 2.0d0, 1.0d0, 1.5d0, 1.5d0, 1.5d0/)
  real(stellar_dp), parameter :: differential_z(n_differential_queries) = &
       (/1.0d-3, 1.0d-2, 5.5d-3, 1.0d-3, 5.5d-3, 5.5d-3/)
  real(stellar_dp), parameter :: differential_age(n_differential_queries) = &
       (/0.0d0, 1.0d0, 0.5d0, 0.5d0, 0.0d0, 0.5d0/)

  failures = 0
  table_filename = 'g1_synthetic_yields.dat'
  result_filename = 'g1_native_interpolation.txt'
  call get_environment_variable('G1_TABLE_PATH', table_filename, &
       status=env_status)
  if (env_status /= 0 .or. len_trim(table_filename) == 0) then
     table_filename = 'g1_synthetic_yields.dat'
  end if
  call get_environment_variable('G1_NATIVE_RESULT', result_filename, &
       status=env_status)
  if (env_status /= 0 .or. len_trim(result_filename) == 0) then
     result_filename = 'g1_native_interpolation.txt'
  end if

  call write_synthetic_table(trim(table_filename), failures)
  call load_yield_table(trim(table_filename), table, ierr)
  call expect(ierr == yield_table_ok, 'ASCII yield table loads', failures)
  call expect(abs(table%age_gyr(2) - 1.0_stellar_dp) < 1.0e-12_stellar_dp, &
       'age_yr is converted to age_gyr in memory', failures)

  call audit_yield_table(table, 1.0e-10_stellar_dp, audit_ierr, .true.)
  call expect(audit_ierr == yield_audit_ok, &
       'complete finite table passes strict audit', failures)

  residual_table = table
  residual_table%ejected_mass(:,1) = 0.9_stellar_dp * &
       residual_table%returned_mass
  call audit_yield_table(residual_table, 1.0e-10_stellar_dp, audit_ierr, .true.)
  call expect(audit_ierr == yield_audit_ok, &
       'nonnegative untracked ejecta residual passes strict audit', failures)

  overfull_table = table
  overfull_table%ejected_mass(2,1) = &
       1.1_stellar_dp * overfull_table%returned_mass(2)
  call audit_yield_table(overfull_table, 1.0e-10_stellar_dp, audit_ierr)
  call expect(iand(audit_ierr, yield_audit_err_mass) /= 0, &
       'tracked ejecta above returned mass are rejected', failures)

  nonmonotonic_residual_table = table
  nonmonotonic_residual_table%ejected_mass(1,1) = &
       nonmonotonic_residual_table%returned_mass(1) - 0.1_stellar_dp
  nonmonotonic_residual_table%ejected_mass(2,1) = &
       nonmonotonic_residual_table%returned_mass(2) - 0.05_stellar_dp
  call audit_yield_table(nonmonotonic_residual_table, 1.0e-10_stellar_dp, &
       audit_ierr)
  call expect(iand(audit_ierr, yield_audit_err_monotonic) /= 0, &
       'decreasing cumulative untracked ejecta are rejected', failures)

  query_mass = 1.5_stellar_dp
  query_z = 5.5e-3_stellar_dp
  query_age = 0.5_stellar_dp
  call interpolate_yield_row(table, channel_snii, query_mass, query_z, &
       query_age, returned_mass, remnant_mass, energy, momentum, ejected_mass, &
       net_yield, ierr)
  expected = expected_returned(channel_snii, query_mass, query_z, query_age)
  call expect(ierr == interpolation_ok, 'interior trilinear interpolation passes', &
       failures)
  call expect_close(returned_mass, expected, 1.0e-12_stellar_dp, &
       'interpolated returned mass agrees with analytic table', failures)
  call expect_close(ejected_mass(1), returned_mass, 1.0e-12_stellar_dp, &
       'ejected mass closes at the interpolated row', failures)

  call set_yield_mass_assignment_mode(table, yield_mass_assignment_piecewise_constant, ierr)
  call expect(ierr == yield_table_ok, 'piecewise source-cell mode can be selected', failures)
  call interpolate_yield_row(table, channel_snii, query_mass, query_z, query_age, &
       returned_mass, remnant_mass, energy, momentum, ejected_mass, net_yield, ierr)
  expected = expected_returned(channel_snii, 1.0_stellar_dp, query_z, query_age)
  call expect(ierr == interpolation_ok, 'piecewise source-cell query succeeds', failures)
  call expect_close(returned_mass, expected, 1.0e-12_stellar_dp, &
       'piecewise source-cell query selects the left mass node', failures)
  call interpolate_yield_row(table, channel_snii, 2.0_stellar_dp, query_z, query_age, &
       returned_mass, remnant_mass, energy, momentum, ejected_mass, net_yield, ierr)
  expected = expected_returned(channel_snii, 2.0_stellar_dp, query_z, query_age)
  call expect_close(returned_mass, expected, 1.0e-12_stellar_dp, &
       'piecewise source-cell upper edge selects the exact node', failures)
  call set_yield_mass_assignment_mode(table, 99, ierr)
  call expect(ierr == yield_table_err_assignment_mode, &
       'unsupported source-cell mode is rejected', failures)
  table%mass_assignment_mode = 99
  call interpolate_yield_row(table, channel_snii, query_mass, query_z, query_age, &
       returned_mass, remnant_mass, energy, momentum, ejected_mass, net_yield, ierr)
  call expect(ierr == interpolation_err_assignment_mode, &
       'invalid source-cell mode cannot enter interpolation', failures)
  call set_yield_mass_assignment_mode(table, yield_mass_assignment_linear, ierr)
  call expect(ierr == yield_table_ok, 'linear mode can be restored after review test', failures)

  call interpolate_yield_row(table, channel_snii, 0.5_stellar_dp, query_z, &
       query_age, returned_mass, remnant_mass, energy, momentum, ejected_mass, &
       net_yield, ierr)
  call expect(ierr == interpolation_err_grid, &
       'below-domain mass query is rejected without clamping', failures)
  call interpolate_yield_row(table, channel_snii, query_mass, query_z, 1.5d0, &
       returned_mass, remnant_mass, energy, momentum, ejected_mass, net_yield, &
       ierr)
  call expect(ierr == interpolation_err_grid, &
       'above-domain age query is rejected without extrapolation', failures)

  incomplete_table = table
  incomplete_table%n_rows = table%n_rows - 1
  call audit_yield_table(incomplete_table, 1.0e-10_stellar_dp, audit_ierr, .true.)
  call expect(iand(audit_ierr, yield_audit_err_grid) /= 0, &
       'missing Cartesian corner is rejected', failures)

  duplicate_table = table
  duplicate_table%age_gyr(table%n_rows) = &
       duplicate_table%age_gyr(table%n_rows - 1)
  duplicate_table%returned_mass(table%n_rows) = &
       duplicate_table%returned_mass(table%n_rows - 1)
  duplicate_table%ejected_mass(table%n_rows,:) = &
       duplicate_table%ejected_mass(table%n_rows - 1,:)
  call audit_yield_table(duplicate_table, 1.0e-10_stellar_dp, audit_ierr, .true.)
  call expect(iand(audit_ierr, yield_audit_err_duplicate) /= 0, &
       'duplicate coordinate tuple is rejected', failures)

  nan_table = table
  nan_table%energy(1) = ieee_value(0.0_stellar_dp, ieee_quiet_nan)
  call audit_yield_table(nan_table, 1.0e-10_stellar_dp, audit_ierr)
  call expect(iand(audit_ierr, yield_audit_err_nonfinite) /= 0, &
       'nonfinite table value is rejected', failures)

  population%formation_time = 0.0_stellar_dp
  population%initial_mass = 1.0e5_stellar_dp
  population%current_mass = 1.0e5_stellar_dp
  population%birth_metallicity = query_z
  population%birth_mass_fraction = 0.0_stellar_dp
  population%imf_id = imf_kroupa
  population%imf_mass_min = 0.08_stellar_dp
  population%imf_mass_max = 120.0_stellar_dp
  population%population_id = 17
  population%binary_fraction = 0.0_stellar_dp
  population%yield_basis_id = 0
  population%pisn_enabled = .false.
  source_tolerance = 1.0e-10_stellar_dp

  call integrate_ssp_channel_increment(table, population, channel_snii, &
       0.0_stellar_dp, 0.25_stellar_dp, 1.0_stellar_dp, 2.0_stellar_dp, 32, &
       source_a, source_ierr)
  call expect(source_ierr == source_increment_ok, &
       'first cumulative source interval succeeds', failures)
  call integrate_ssp_channel_increment(table, population, channel_snii, &
       0.25_stellar_dp, 1.0_stellar_dp, 1.0_stellar_dp, 2.0_stellar_dp, 32, &
       source_b, source_ierr)
  call expect(source_ierr == source_increment_ok, &
       'second cumulative source interval succeeds', failures)
  call integrate_ssp_channel_increment(table, population, channel_snii, &
       0.0_stellar_dp, 1.0_stellar_dp, 1.0_stellar_dp, 2.0_stellar_dp, 32, &
       source_full, source_ierr)
  call expect(source_ierr == source_increment_ok, &
       'full cumulative source interval succeeds', failures)
  call expect_close(source_a%returned_mass + source_b%returned_mass, &
       source_full%returned_mass, source_tolerance, &
       'variable-width intervals telescope in returned mass', failures)
  call expect_close(source_a%energy + source_b%energy, source_full%energy, &
       source_tolerance, 'variable-width intervals telescope in energy', failures)
  call expect_close(source_a%ejected_mass(1) + source_b%ejected_mass(1), &
       source_full%ejected_mass(1), source_tolerance, &
       'variable-width intervals telescope in ejecta', failures)

  call integrate_ssp_channel_increment(table, population, channel_snii, &
       0.0_stellar_dp, 0.0_stellar_dp, 1.0_stellar_dp, 2.0_stellar_dp, 32, &
       negative_source, source_ierr)
  call expect(source_ierr == source_increment_ok .and. &
       negative_source%returned_mass == 0.0_stellar_dp, &
       'zero-width source interval is an exact zero', failures)

  negative_table = table
  do i = 1, negative_table%n_rows
     if (negative_table%channel(i) /= channel_snii) cycle
     if (abs(negative_table%age_gyr(i) - 1.0_stellar_dp) > 1.0e-12_stellar_dp) cycle
     negative_table%returned_mass(i) = -0.01_stellar_dp
     negative_table%ejected_mass(i,:) = 0.0_stellar_dp
     negative_table%ejected_mass(i,1) = negative_table%returned_mass(i)
  end do
  call integrate_ssp_channel_increment(negative_table, population, channel_snii, &
       0.0_stellar_dp, 1.0_stellar_dp, 1.0_stellar_dp, 2.0_stellar_dp, 32, &
       negative_source, source_ierr)
  call expect(source_ierr == source_increment_err_negative, &
       'negative physical source increment is rejected', failures)

  call code_time_to_age_gyr(2.0_stellar_dp, 1.0d9 * 365.25d0 * 86400.0d0, &
       0.5_stellar_dp, native_age, units_ierr)
  call expect(units_ierr == units_ok, 'code-time conversion succeeds', failures)
  call expect_close(native_age, 8.0_stellar_dp, 1.0e-14_stellar_dp, &
       'code-time conversion retains explicit aexp**2', failures)
  call code_interval_to_age_gyr(0.25_stellar_dp, &
       1.0d9 * 365.25d0 * 86400.0d0, 0.5_stellar_dp, native_dt, units_ierr)
  call expect_close(native_dt, 1.0_stellar_dp, 1.0e-14_stellar_dp, &
       'code interval conversion uses the same convention', failures)
  call energy_erg_to_code(2.0d51, 4.0d49, native_energy, units_ierr)
  call expect(units_ierr == units_ok, 'physical energy conversion succeeds', failures)
  call expect_close(native_energy, 50.0_stellar_dp, 1.0e-14_stellar_dp, &
       'energy conversion is erg/scale_energy with no magic 1e51', failures)
  call mass_code_to_msun(2.0_stellar_dp, 3.0_stellar_dp, mass_msun, units_ierr)
  call mass_msun_to_code(mass_msun, 3.0_stellar_dp, mass_code, units_ierr)
  call expect_close(mass_code, 2.0_stellar_dp, 1.0e-14_stellar_dp, &
       'mass conversion round trip succeeds', failures)

  call clear_field_map(field_map)
  field_map%density_index = 1
  field_map%momentum_index = (/2, 3, 4/)
  field_map%energy_index = 5
  field_map%total_metal_index = 6
  field_map%delayed_cooling_index = 7
  field_map%element_index = (/8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18/)
  call validate_field_map(field_map, 18, 3, map_ierr)
  call expect(map_ierr == 0, 'versioned NVAR=18 field map passes', failures)
  call validate_field_map(field_map, 17, 3, map_ierr)
  call expect(map_ierr /= 0, 'NVAR mismatch is rejected', failures)
  field_map%delayed_cooling_index = 8
  call validate_field_map(field_map, 18, 3, map_ierr)
  call expect(map_ierr /= 0, 'delayed-cooling/element overlap is rejected', failures)

  call progress_initialize(progress, 0.0_stellar_dp, progress_ierr)
  call expect(progress_ierr == progress_ok, 'progress initializes from restart state', &
       failures)
  call progress_begin(progress, 1.0_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call expect(progress_ierr == progress_ok .and. should_deposit .and. &
       dt_gyr == 1.0_stellar_dp, 'new age opens one pending interval', failures)
  call progress_commit(progress, progress_ierr)
  call progress_export(progress, committed_age, progress_ierr)
  call expect(progress_ierr == progress_ok .and. committed_age == 1.0_stellar_dp, &
       'commit exports the deposited age', failures)
  call progress_begin(progress, 1.0_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call expect(progress_ierr == progress_ok .and. .not. should_deposit .and. &
       dt_gyr == 0.0_stellar_dp, 'repeated committed age is an exact no-op', failures)
  call progress_initialize(restarted, committed_age, progress_ierr)
  call progress_begin(restarted, 1.0_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call expect(progress_ierr == progress_ok .and. .not. should_deposit, &
       'restart restores idempotent committed progress', failures)
  call progress_begin(restarted, 2.0_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call progress_abort(restarted, progress_ierr)
  call progress_begin(restarted, 2.0_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call expect(should_deposit .and. dt_gyr == 1.0_stellar_dp, &
       'aborted interval can be retried without advancing progress', failures)
  call progress_commit(restarted, progress_ierr)
  call progress_begin(restarted, 1.5_stellar_dp, 1.0e-12_stellar_dp, &
       should_deposit, dt_gyr, progress_ierr)
  call expect(progress_ierr /= progress_ok, 'stale age is rejected after commit', failures)

  call interpolate_yield_row(table, channel_snii, query_mass, query_z, &
       query_age, returned_mass, remnant_mass, energy, momentum, ejected_mass, &
       net_yield, ierr)
  call expect(ierr == interpolation_ok, 'native differential reference query succeeds', &
       failures)
  open(newunit=unit, file=trim(result_filename), status='replace', &
       action='write', iostat=ierr)
  if (ierr == 0) then
     do i = 1, n_differential_queries
        query_mass = differential_masses(i)
        query_z = differential_z(i)
        query_age = differential_age(i)
        call interpolate_yield_row(table, channel_snii, query_mass, query_z, &
             query_age, returned_mass, remnant_mass, energy, momentum, &
             ejected_mass, net_yield, ierr)
        call expect(ierr == interpolation_ok, 'native differential matrix query succeeds', &
             failures)
        write(unit, '(9(1x,es26.17))') query_mass, query_z, query_age, &
             returned_mass, remnant_mass, energy, momentum(1), ejected_mass(1), &
             net_yield(1)
     end do
     close(unit)
  else
     call expect(.false., 'native interpolation result file opens', failures)
  end if

  call clear_yield_table(table)
  call clear_yield_table(incomplete_table)
  call clear_yield_table(duplicate_table)
  call clear_yield_table(nan_table)
  call clear_yield_table(negative_table)
  call clear_yield_table(residual_table)
  call clear_yield_table(overfull_table)
  call clear_yield_table(nonmonotonic_residual_table)

  if (failures == 0) then
     write(*, '(a)') 'G1_NATIVE_CONTRACT_TEST_OK'
  else
     write(*, '(a,i0)') 'G1_NATIVE_CONTRACT_TEST_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine write_synthetic_table(filename, failures)
    character(len=*), intent(in) :: filename
    integer, intent(inout) :: failures
    integer :: unit, ios, channel, im, iz, ia
    real(stellar_dp) :: masses(2), metallicities(2), ages(2)
    real(stellar_dp) :: returned, remnant, energy, age_yr
    real(stellar_dp) :: momentum(3), ejecta(n_stellar_elements)
    real(stellar_dp) :: net(n_stellar_elements)

    masses = (/1.0_stellar_dp, 2.0_stellar_dp/)
    metallicities = (/1.0e-3_stellar_dp, 1.0e-2_stellar_dp/)
    ages = (/0.0_stellar_dp, 1.0_stellar_dp/)
    open(newunit=unit, file=filename, status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) then
       call expect(.false., 'synthetic table file opens', failures)
       return
    end if
    do channel = 1, 3
       do im = 1, 2
          do iz = 1, 2
             do ia = 1, 2
                returned = expected_returned(channel, masses(im), &
                     metallicities(iz), ages(ia))
                remnant = ages(ia) * 0.05_stellar_dp * &
                     real(channel, stellar_dp)
                energy = expected_energy(channel, ages(ia))
                age_yr = ages(ia) * 1.0e9_stellar_dp
                momentum = 0.0_stellar_dp
                momentum(1) = ages(ia) * &
                     (real(channel, stellar_dp) * 1.0e30_stellar_dp + &
                     ages(ia) * 1.0e29_stellar_dp)
                ejecta = 0.0_stellar_dp
                ejecta(1) = returned
                net = 0.0_stellar_dp
                net(1) = 0.01_stellar_dp * ages(ia)
                write(unit, '(i4,31(1x,es24.16))') channel, masses(im), &
                     metallicities(iz), age_yr, &
                     returned, remnant, energy, momentum, ejecta, net
             end do
          end do
       end do
    end do
    close(unit)
  end subroutine write_synthetic_table

  real(stellar_dp) function expected_returned(channel, mass, metallicity, age)
    integer, intent(in) :: channel
    real(stellar_dp), intent(in) :: mass, metallicity, age

    expected_returned = age * (0.1_stellar_dp * real(channel, stellar_dp) + &
         0.05_stellar_dp * mass + 0.2_stellar_dp * metallicity + &
         0.4_stellar_dp)
  end function expected_returned

  real(stellar_dp) function expected_energy(channel, age)
    integer, intent(in) :: channel
    real(stellar_dp), intent(in) :: age

    expected_energy = age * (real(channel, stellar_dp) * 1.0e48_stellar_dp + &
         1.0e47_stellar_dp)
  end function expected_energy

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

  subroutine expect_close(actual, expected, tolerance, label, failures)
    real(stellar_dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    real(stellar_dp) :: scale

    scale = max(1.0_stellar_dp, abs(actual), abs(expected))
    call expect(abs(actual - expected) <= tolerance * scale, label, failures)
  end subroutine expect_close

end program g1_contract_test

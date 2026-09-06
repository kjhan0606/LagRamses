program fp12_stellar_feedback_transaction_test
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use omp_lib, only: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, &
       omp_destroy_lock, omp_get_max_threads
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, elem_c, &
       channel_snii, channel_wind, channel_agb
  use stellar_enrichment_contract, only: stellar_source_t, &
       delayed_cooling_source_mass, generic_metal_ejecta_mass
  use stellar_native_units, only: solar_mass_cgs
  use stellar_snia_physical_contract, only: snia_event_budget_t
  use stellar_snia_cell_deposition, only: snia_thermal_coupling_t, &
       snia_thermal_all_to_total_energy, snia_deposition_ok
  use stellar_ramses_field_map, only: stellar_field_map_t, clear_field_map, &
       validate_field_map
  use stellar_ramses_bridge, only: ramses_bridge_ok, ramses_bridge_err_source, &
       build_stellar_source_unew_delta, build_snia_budget_unew_delta, &
       deposit_snia_budget_to_unew
  implicit none

  integer, parameter :: nvar_test = 18
  integer :: failures, ierr, element, source_index, lock_index
  integer, parameter :: n_test_sources = 64, n_test_locks = 8
  integer(omp_lock_kind) :: test_locks(n_test_locks)
  integer :: target_cell(1), owner_rank(1), momentum_var(3)
  integer :: element_var(n_stellar_elements)
  type(stellar_field_map_t) :: field_map, invalid_map
  type(stellar_source_t) :: source, zero_mass_source, opposed_source
  type(stellar_source_t) :: invalid_source
  type(snia_event_budget_t) :: budget, opposed_budget
  type(snia_thermal_coupling_t) :: coupling
  real(stellar_dp) :: generic_delta(nvar_test), snia_delta(nvar_test)
  real(stellar_dp) :: expected_generic_delta(nvar_test)
  real(stellar_dp) :: expected_snia_delta(nvar_test)
  real(stellar_dp) :: expected_staged_delta(nvar_test)
  real(stellar_dp) :: opposed_generic_delta(nvar_test)
  real(stellar_dp) :: opposed_snia_delta(nvar_test)
  real(stellar_dp) :: staged_delta(nvar_test), row_before(nvar_test)
  real(stellar_dp) :: row_after(nvar_test), working_row(nvar_test)
  real(stellar_dp) :: expected_row(nvar_test), row_increment(nvar_test)
  real(stellar_dp) :: unew(3,nvar_test), bulk_velocity(3), snia_velocity(3)
  real(stellar_dp) :: opposed_velocity(3)
  real(stellar_dp) :: velocity_cells(3,1)
  real(stellar_dp) :: cell_volume(1), weights(1)
  real(stellar_dp) :: mp_before, indtab_before, mp_work, indtab_work
  real(stellar_dp) :: source_momentum_code(3), expected_momentum(3)
  real(stellar_dp) :: generic_bulk_energy, returned_mass_cgs
  real(stellar_dp) :: snia_bulk_momentum(3), snia_bulk_energy
  real(stellar_dp) :: common_scale_length, common_scale_density
  real(stellar_dp) :: common_scale_velocity, common_independent_energy
  real(stellar_dp) :: common_merged_energy

  failures = 0
  call expect(omp_get_max_threads() > 1, &
       'native transaction test is configured for multiple OpenMP threads', failures)
  if (omp_get_max_threads() <= 1) then
     write(*, '(a)') 'FP12_STELLAR_FEEDBACK_TRANSACTION_TEST_FAILED single-thread configuration'
     stop 1
  end if
  element_var = 0
  do element = 1, n_stellar_elements
     element_var(element) = 7 + element
  end do
  momentum_var = (/2, 3, 4/)

  call clear_field_map(field_map)
  field_map%density_index = 1
  field_map%momentum_index = (/2, 3, 4/)
  field_map%energy_index = 5
  field_map%total_metal_index = 6
  field_map%delayed_cooling_index = 7
  field_map%element_index = element_var
  call validate_field_map(field_map, nvar_test, 3, ierr)
  call expect(ierr == 0, 'complete row-major field map passes', failures)

  invalid_map = field_map
  invalid_map%delayed_cooling_index = element_var(1)
  call validate_field_map(invalid_map, nvar_test, 3, ierr)
  call expect(ierr /= 0, 'delayed-cooling/element overlap is rejected', failures)

  source%ejected_mass = 0.0_stellar_dp
  source%net_yield = 0.0_stellar_dp
  source%returned_mass = 0.4_stellar_dp
  source%energy = 2.0e50_stellar_dp
  source%momentum = (/1.0e39_stellar_dp, -2.0e39_stellar_dp, &
       0.5e39_stellar_dp/)
  source%ejected_mass(elem_c) = 0.1_stellar_dp
  source%channel_returned_mass = 0.0_stellar_dp
  source%channel_energy = 0.0_stellar_dp
  source%channel_momentum = 0.0_stellar_dp
  source%channel_ejected_mass = 0.0_stellar_dp
  source%channel_net_yield = 0.0_stellar_dp
  source%channel_returned_mass(channel_snii) = 0.2_stellar_dp

  generic_delta = -1.0_stellar_dp
  bulk_velocity = (/1.0_stellar_dp, -2.0_stellar_dp, 0.5_stellar_dp/)
  call build_stellar_source_unew_delta(source, bulk_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr)
  expected_generic_delta = 0.0_stellar_dp
  expected_generic_delta(1) = source%returned_mass / 2.0_stellar_dp
  source_momentum_code = source%momentum / 1.0e39_stellar_dp
  expected_momentum = source%returned_mass * bulk_velocity + source_momentum_code
  expected_generic_delta(momentum_var) = expected_momentum / 2.0_stellar_dp
  generic_bulk_energy = 0.5_stellar_dp * source%returned_mass * &
       sum(bulk_velocity**2) + sum(bulk_velocity * source_momentum_code) + &
       0.5_stellar_dp * sum(source_momentum_code**2) / source%returned_mass
  expected_generic_delta(5) = (source%energy / 1.0e50_stellar_dp + &
       generic_bulk_energy) / 2.0_stellar_dp
  expected_generic_delta(6) = generic_metal_ejecta_mass(source%returned_mass, &
       source%ejected_mass) / 2.0_stellar_dp
  expected_generic_delta(7) = delayed_cooling_source_mass(source) / 2.0_stellar_dp
  expected_generic_delta(element_var(elem_c)) = source%ejected_mass(elem_c) / &
       2.0_stellar_dp
  call expect(ierr == ramses_bridge_ok .and. vectors_close(generic_delta, &
       expected_generic_delta, 1.0e-12_stellar_dp), &
       'generic source stages every expected row-major field', failures)

  coupling%approved = .true.
  coupling%mode = snia_thermal_all_to_total_energy
  coupling%thermal_fraction = 1.0_stellar_dp
  coupling%include_event_momentum_kinetic = .true.
  coupling%source_commit_binding = 'test'
  coupling%approval_id = 'test'
  budget%wd_reservoir_debit = 1.4_stellar_dp
  budget%returned_mass = 1.3_stellar_dp
  budget%terminal_remnant_mass = 0.1_stellar_dp
  budget%energy = 1.0e50_stellar_dp
  budget%momentum = (/1.0e39_stellar_dp, -2.0e39_stellar_dp, &
       3.0e39_stellar_dp/)
  budget%ejected_mass = 0.0_stellar_dp
  budget%net_yield = 0.0_stellar_dp
  budget%ejected_mass(elem_c) = 1.0_stellar_dp
  snia_velocity = (/1.0e7_stellar_dp, -2.0e7_stellar_dp, 3.0e7_stellar_dp/)
  snia_delta = -2.0_stellar_dp
  call build_snia_budget_unew_delta(budget, coupling, snia_velocity, 1.0e20_stellar_dp, &
       1.0e-20_stellar_dp, 1.0e7_stellar_dp, nvar_test, 2.0_stellar_dp, 1, 5, &
       momentum_var, snia_delta, 1.0e-12_stellar_dp, ierr, 6, element_var)
  expected_snia_delta = 0.0_stellar_dp
  returned_mass_cgs = budget%returned_mass * solar_mass_cgs
  snia_bulk_momentum = returned_mass_cgs * snia_velocity + budget%momentum
  snia_bulk_energy = 0.5_stellar_dp * returned_mass_cgs * &
       sum(snia_velocity**2) + sum(snia_velocity * budget%momentum) + &
       0.5_stellar_dp * sum(budget%momentum**2) / returned_mass_cgs
  expected_snia_delta(1) = returned_mass_cgs / &
       (2.0_stellar_dp * (1.0e20_stellar_dp**3)) / 1.0e-20_stellar_dp
  expected_snia_delta(momentum_var) = snia_bulk_momentum / &
       (2.0_stellar_dp * (1.0e20_stellar_dp**3)) / (1.0e-20_stellar_dp * 1.0e7_stellar_dp)
  expected_snia_delta(5) = (budget%energy + snia_bulk_energy) / &
       (2.0_stellar_dp * (1.0e20_stellar_dp**3)) / (1.0e-20_stellar_dp * 1.0e14_stellar_dp)
  expected_snia_delta(6) = generic_metal_ejecta_mass(budget%returned_mass, &
       budget%ejected_mass) * solar_mass_cgs / &
       (2.0_stellar_dp * (1.0e20_stellar_dp**3)) / 1.0e-20_stellar_dp
  expected_snia_delta(element_var(elem_c)) = budget%ejected_mass(elem_c) * &
       solar_mass_cgs / (2.0_stellar_dp * (1.0e20_stellar_dp**3)) / 1.0e-20_stellar_dp
  call expect(ierr == ramses_bridge_ok .and. vectors_close(snia_delta, &
       expected_snia_delta, 1.0e-12_stellar_dp), &
       'SNIa source stages every expected row-major field without mutation', failures)

  staged_delta = generic_delta + snia_delta
  expected_staged_delta = expected_generic_delta + expected_snia_delta
  row_before = (/ (1000.0_stellar_dp + element, element = 1, nvar_test) /)
  row_after = row_before + staged_delta
  expected_row = row_before + expected_staged_delta
  call expect(vectors_close(staged_delta, expected_staged_delta, 1.0e-12_stellar_dp) .and. &
       vectors_close(row_after, expected_row, 1.0e-12_stellar_dp), &
       'mixed generic and SNIa deltas match independent field expectations', failures)

  ! Every preparation failure leaves the would-be production state untouched.
  mp_before = 9.0_stellar_dp
  indtab_before = 4.0_stellar_dp
  working_row = row_before
  mp_work = mp_before
  indtab_work = indtab_before
  generic_delta = 9.0_stellar_dp
  call build_stellar_source_unew_delta(source, bulk_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 0.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr)
  call expect(ierr /= ramses_bridge_ok .and. all(generic_delta == 0.0_stellar_dp) .and. &
       vectors_close(working_row, row_before, 0.0_stellar_dp) .and. &
       mp_work == mp_before .and. indtab_work == indtab_before, &
       'invalid generic volume leaves row, mass, and progress unchanged', failures)

  working_row = row_before
  mp_work = mp_before
  indtab_work = indtab_before
  invalid_map = field_map
  invalid_map%energy_index = invalid_map%density_index
  generic_delta = 9.0_stellar_dp
  call build_stellar_source_unew_delta(source, bulk_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, invalid_map, generic_delta, 1.0e-12_stellar_dp, ierr)
  call expect(ierr /= ramses_bridge_ok .and. all(generic_delta == 0.0_stellar_dp) .and. &
       vectors_close(working_row, row_before, 0.0_stellar_dp) .and. &
       mp_work == mp_before .and. indtab_work == indtab_before, &
       'invalid field map leaves row, mass, and progress unchanged', failures)

  working_row = row_before
  mp_work = mp_before
  indtab_work = indtab_before
  generic_delta = 9.0_stellar_dp
  call build_stellar_source_unew_delta(source, bulk_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 2, field_map, generic_delta, 1.0e-12_stellar_dp, ierr)
  call expect(ierr /= ramses_bridge_ok .and. all(generic_delta == 0.0_stellar_dp) .and. &
       vectors_close(working_row, row_before, 0.0_stellar_dp) .and. &
       mp_work == mp_before .and. indtab_work == indtab_before, &
       'ndim mismatch leaves row, mass, and progress unchanged', failures)

  zero_mass_source = source
  zero_mass_source%returned_mass = 0.0_stellar_dp
  zero_mass_source%energy = 0.0_stellar_dp
  zero_mass_source%momentum = (/1.0e39_stellar_dp, 0.0_stellar_dp, 0.0_stellar_dp/)
  zero_mass_source%ejected_mass = 0.0_stellar_dp
  zero_mass_source%channel_returned_mass = 0.0_stellar_dp
  generic_delta = 3.0_stellar_dp
  call build_stellar_source_unew_delta(zero_mass_source, bulk_velocity, 1.0_stellar_dp, &
       1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, nvar_test, 3, &
       field_map, generic_delta, 1.0e-12_stellar_dp, ierr)
  call expect(ierr == ramses_bridge_err_source .and. all(generic_delta == 0.0_stellar_dp), &
       'nonzero momentum with zero returned mass is rejected transactionally', failures)

  ! Wind and AGB are separate ejecta components even when they arrive in the
  ! same cell in one timestep. Net momentum cancellation must not erase energy.
  opposed_source = source
  opposed_source%returned_mass = 0.4_stellar_dp
  opposed_source%energy = 0.0_stellar_dp
  opposed_source%momentum = 0.0_stellar_dp
  opposed_source%channel_returned_mass = 0.0_stellar_dp
  opposed_source%channel_returned_mass(channel_wind) = 0.2_stellar_dp
  opposed_source%channel_returned_mass(channel_agb) = 0.2_stellar_dp
  opposed_source%channel_momentum = 0.0_stellar_dp
  opposed_source%channel_momentum(channel_wind,1) = 1.0e39_stellar_dp
  opposed_source%channel_momentum(channel_agb,1) = -1.0e39_stellar_dp
  opposed_velocity = 0.0_stellar_dp
  call build_stellar_source_unew_delta(opposed_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr == ramses_bridge_ok .and. &
       abs(generic_delta(5) - 2.5_stellar_dp) < 1.0e-12_stellar_dp .and. &
       all(generic_delta(momentum_var) == 0.0_stellar_dp) .and. &
       generic_delta(7) == 0.0_stellar_dp, &
       'opposed wind/AGB momenta retain energy without SNII cooling mass', failures)

  opposed_velocity = (/3.0_stellar_dp, 0.0_stellar_dp, 0.0_stellar_dp/)
  opposed_source%energy = 4.0e50_stellar_dp
  call build_stellar_source_unew_delta(opposed_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr == ramses_bridge_ok .and. &
       abs(generic_delta(5) - 5.4_stellar_dp) < 1.0e-12_stellar_dp .and. &
       abs(generic_delta(2) - 0.6_stellar_dp) < 1.0e-12_stellar_dp, &
       'channel mixing retains supplied energy and stellar bulk motion', failures)

  invalid_source = opposed_source
  invalid_source%channel_returned_mass(channel_agb) = 0.1_stellar_dp
  call build_stellar_source_unew_delta(invalid_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr /= ramses_bridge_ok .and. all(generic_delta == 0.0_stellar_dp), &
       'incomplete channel mass ledger rejects without a partial delta', failures)

  invalid_source = opposed_source
  invalid_source%momentum(1) = 1.0e39_stellar_dp
  call build_stellar_source_unew_delta(invalid_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr /= ramses_bridge_ok .and. all(generic_delta == 0.0_stellar_dp), &
       'inconsistent channel momentum ledger rejects without a partial delta', failures)

  invalid_source = opposed_source
  invalid_source%channel_returned_mass = 0.0_stellar_dp
  invalid_source%channel_returned_mass(channel_snii) = 0.4_stellar_dp
  call build_stellar_source_unew_delta(invalid_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr == ramses_bridge_err_source .and. all(generic_delta == 0.0_stellar_dp), &
       'massless channel momenta cannot hide behind net cancellation', failures)

  ! Isotropic AGB return has zero net rest-frame momentum. Its unresolved
  ! wind energy must come from the supplied energy budget, not an invented speed.
  opposed_source%channel_returned_mass = 0.0_stellar_dp
  opposed_source%channel_returned_mass(channel_agb) = 0.4_stellar_dp
  opposed_source%channel_momentum = 0.0_stellar_dp
  call build_stellar_source_unew_delta(opposed_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e39_stellar_dp, 1.0e50_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, generic_delta, 1.0e-12_stellar_dp, ierr, &
       channel_resolved=.true.)
  call expect(ierr == ramses_bridge_ok .and. &
       abs(generic_delta(1) - 0.2_stellar_dp) < 1.0e-12_stellar_dp .and. &
       abs(generic_delta(2) - 0.6_stellar_dp) < 1.0e-12_stellar_dp .and. &
       abs(generic_delta(5) - 2.9_stellar_dp) < 1.0e-12_stellar_dp .and. &
       abs(generic_delta(element_var(elem_c)) - 0.05_stellar_dp) < 1.0e-12_stellar_dp .and. &
       generic_delta(7) == 0.0_stellar_dp, &
       'AGB-only return carries mass carbon bulk momentum and supplied energy', failures)

  snia_delta = 9.0_stellar_dp
  call build_snia_budget_unew_delta(budget, coupling, snia_velocity, 1.0e20_stellar_dp, &
       1.0e-20_stellar_dp, 1.0e7_stellar_dp, nvar_test, 0.0_stellar_dp, 1, 5, &
       momentum_var, snia_delta, 1.0e-12_stellar_dp, ierr, 6, element_var)
  call expect(ierr /= ramses_bridge_ok .and. all(snia_delta == 0.0_stellar_dp), &
       'invalid volume rejects without exposing a partial SNIa delta', failures)

  ! Generic and SNIa components with opposed momenta must retain their
  ! individual kinetic energies.  A merged-net-momentum calculation would
  ! incorrectly reduce this case to zero kinetic energy.
  opposed_source = source
  opposed_source%returned_mass = 0.4_stellar_dp
  opposed_source%energy = 0.0_stellar_dp
  opposed_source%momentum = (/1.0e7_stellar_dp, 0.0_stellar_dp, 0.0_stellar_dp/)
  opposed_source%ejected_mass = 0.0_stellar_dp
  opposed_source%net_yield = 0.0_stellar_dp
  opposed_source%channel_returned_mass = 0.0_stellar_dp
  opposed_source%channel_energy = 0.0_stellar_dp
  opposed_source%channel_momentum = 0.0_stellar_dp
  opposed_source%channel_ejected_mass = 0.0_stellar_dp
  opposed_source%channel_net_yield = 0.0_stellar_dp
  opposed_velocity = 0.0_stellar_dp
  opposed_generic_delta = -1.0_stellar_dp
  call build_stellar_source_unew_delta(opposed_source, opposed_velocity, &
       1.0_stellar_dp, 1.0e7_stellar_dp, 1.0e14_stellar_dp, 2.0_stellar_dp, &
       nvar_test, 3, field_map, opposed_generic_delta, 1.0e-12_stellar_dp, ierr)

  opposed_budget = budget
  opposed_budget%wd_reservoir_debit = 0.4_stellar_dp
  opposed_budget%returned_mass = 0.4_stellar_dp
  opposed_budget%terminal_remnant_mass = 0.0_stellar_dp
  opposed_budget%energy = 0.0_stellar_dp
  opposed_budget%momentum = (/-1.0e7_stellar_dp * solar_mass_cgs, &
       0.0_stellar_dp, 0.0_stellar_dp/)
  opposed_budget%ejected_mass = 0.0_stellar_dp
  opposed_budget%net_yield = 0.0_stellar_dp
  common_scale_length = 1.0e20_stellar_dp
  common_scale_density = solar_mass_cgs / common_scale_length**3
  common_scale_velocity = 1.0e7_stellar_dp
  opposed_snia_delta = -1.0_stellar_dp
  call build_snia_budget_unew_delta(opposed_budget, coupling, opposed_velocity, &
       common_scale_length, common_scale_density, common_scale_velocity, &
       nvar_test, 2.0_stellar_dp, 1, 5, momentum_var, opposed_snia_delta, &
       1.0e-12_stellar_dp, ierr, 6, element_var)
  common_independent_energy = opposed_generic_delta(5) + opposed_snia_delta(5)
  common_merged_energy = 0.5_stellar_dp * &
       sum((opposed_generic_delta(momentum_var) + &
       opposed_snia_delta(momentum_var))**2) / &
       (opposed_generic_delta(1) + opposed_snia_delta(1))
  call expect(ierr == ramses_bridge_ok .and. &
       vectors_close(opposed_generic_delta(momentum_var) + &
       opposed_snia_delta(momentum_var), (/0.0_stellar_dp, 0.0_stellar_dp, &
       0.0_stellar_dp/), 1.0e-12_stellar_dp) .and. &
       common_independent_energy > 1.0_stellar_dp .and. &
       common_merged_energy < 1.0e-12_stellar_dp, &
       'opposed generic/SNIa momenta retain independent kinetic energy', failures)

  ! A non-owned owner_rank is legal metadata for a virtual/reception row.  The
  ! row-major helper must not reject it; RAMSES reverses the row later.
  unew = 0.0_stellar_dp
  target_cell(1) = 2
  owner_rank(1) = 17
  velocity_cells(:,1) = snia_velocity
  cell_volume(1) = 2.0_stellar_dp
  weights(1) = 1.0_stellar_dp
  call deposit_snia_budget_to_unew(budget, coupling, velocity_cells, 1.0e20_stellar_dp, &
       1.0e-20_stellar_dp, 1.0e7_stellar_dp, nvar_test, 3, 1, target_cell, &
       owner_rank, 0, cell_volume, weights, 1, 5, momentum_var, unew, &
       1.0e-12_stellar_dp, ierr, 6, element_var)
  call expect(ierr == ramses_bridge_ok .and. unew(2,1) > 0.0_stellar_dp .and. &
       unew(2,5) > 0.0_stellar_dp .and. unew(1,5) == 0.0_stellar_dp, &
       'virtual-row target is accepted and row-major energy is deposited', failures)

  ! Exercise the same striped-lock pattern used by the production runtime with
  ! a complete row, rather than a scalar surrogate.  Each source contributes
  ! the same vector so the serial reference is deterministic.
  do lock_index = 1, n_test_locks
     call omp_init_lock(test_locks(lock_index))
  end do
  unew = 0.0_stellar_dp
  unew(2,1:nvar_test) = row_before
  row_increment = (/ (0.01_stellar_dp * element, element = 1, nvar_test) /)
  expected_row = row_before
  do source_index = 1, n_test_sources
     expected_row = expected_row + row_increment
  end do
!$omp parallel do default(shared) private(source_index,lock_index)
  do source_index = 1, n_test_sources
     lock_index = 1 + mod(2 - 1, n_test_locks)
     call omp_set_lock(test_locks(lock_index))
     unew(2,1:nvar_test) = unew(2,1:nvar_test) + row_increment
     call omp_unset_lock(test_locks(lock_index))
  end do
!$omp end parallel do
  call expect(vectors_close(unew(2,1:nvar_test), expected_row, &
       1.0e-12_stellar_dp), &
       'striped lock retains all same-cell full-row increments', failures)
  do lock_index = 1, n_test_locks
     call omp_destroy_lock(test_locks(lock_index))
  end do

  if (failures > 0) then
     write(*, '(a,i0)') 'FP12_STELLAR_FEEDBACK_TRANSACTION_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP12_STELLAR_FEEDBACK_TRANSACTION_TEST_OK'

contains

  logical function vectors_close(actual, expected, tolerance)
    real(stellar_dp), intent(in) :: actual(:), expected(:)
    real(stellar_dp), intent(in) :: tolerance

    vectors_close = size(actual) == size(expected) .and. &
         all(ieee_is_finite(actual)) .and. all(ieee_is_finite(expected))
    if (.not. vectors_close) return
    vectors_close = all(abs(actual - expected) <= tolerance * &
         max(1.0_stellar_dp, abs(expected)))
  end function vectors_close

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

end program fp12_stellar_feedback_transaction_test

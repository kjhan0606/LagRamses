program fp2_snia_cell_deposition_test
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, elem_c
  use stellar_native_units, only: solar_mass_cgs
  use stellar_snia_physical_contract, only: snia_event_budget_t
  use stellar_snia_cell_deposition, only: snia_thermal_coupling_t, &
       snia_cell_increment_t, snia_thermal_all_to_total_energy, &
       snia_deposition_ok, snia_deposition_err_argument, &
       snia_deposition_err_policy, snia_deposition_err_budget, &
       validate_snia_thermal_coupling, build_snia_cell_increment
  implicit none

  type(snia_event_budget_t) :: budget
  type(snia_thermal_coupling_t) :: coupling
  type(snia_cell_increment_t) :: increment
  real(stellar_dp) :: volume_cm3, bulk_velocity_cm_s(3)
  real(stellar_dp) :: expected_mass_density, expected_momentum_density(3)
  real(stellar_dp) :: expected_event_energy_density
  real(stellar_dp) :: expected_bulk_energy_density, returned_mass_cgs
  integer :: failures, ierr

  failures = 0
  call validate_snia_thermal_coupling(coupling, ierr)
  call expect(ierr == snia_deposition_err_policy, &
       'unapproved thermal coupling policy is rejected', failures)

  coupling%approved = .true.
  coupling%mode = snia_thermal_all_to_total_energy
  coupling%thermal_fraction = 0.5_stellar_dp
  call validate_snia_thermal_coupling(coupling, ierr)
  call expect(ierr == snia_deposition_err_policy, &
       'fractional thermal policy without a non-thermal receiver is rejected', failures)

  coupling%thermal_fraction = 1.0_stellar_dp
  coupling%include_event_momentum_kinetic = .true.
  call validate_snia_thermal_coupling(coupling, ierr)
  call expect(ierr == snia_deposition_ok, &
       'approved all-to-total-energy policy is valid', failures)

  budget%wd_reservoir_debit = 1.4_stellar_dp
  budget%returned_mass = 1.3_stellar_dp
  budget%terminal_remnant_mass = 0.1_stellar_dp
  budget%energy = 1.0e51_stellar_dp
  budget%momentum = (/1.0e40_stellar_dp, -2.0e39_stellar_dp, &
       3.0e39_stellar_dp/)
  budget%ejected_mass(elem_c) = 1.0_stellar_dp
  volume_cm3 = 1.0e60_stellar_dp
  bulk_velocity_cm_s = (/1.0e7_stellar_dp, -2.0e7_stellar_dp, &
       3.0e7_stellar_dp/)

  call build_snia_cell_increment(budget, volume_cm3, bulk_velocity_cm_s, &
       coupling, increment, ierr)
  call expect(ierr == snia_deposition_ok, &
       'valid SNIa budget produces a cell increment', failures)
  expected_mass_density = budget%returned_mass * solar_mass_cgs / volume_cm3
  expected_momentum_density = &
       (budget%returned_mass * solar_mass_cgs * bulk_velocity_cm_s + &
       budget%momentum) / volume_cm3
  expected_event_energy_density = budget%energy / volume_cm3
  returned_mass_cgs = budget%returned_mass * solar_mass_cgs
  expected_bulk_energy_density = (0.5_stellar_dp * returned_mass_cgs * &
       sum(bulk_velocity_cm_s**2) + sum(bulk_velocity_cm_s * budget%momentum) + &
       0.5_stellar_dp * sum(budget%momentum**2) / returned_mass_cgs) / volume_cm3
  call expect_close(increment%mass_density, expected_mass_density, &
       'returned mass is converted to cell mass density', failures)
  call expect(all(abs(increment%momentum_density - expected_momentum_density) <= &
       1.0e-12_stellar_dp * max(1.0_stellar_dp, &
       abs(expected_momentum_density))), &
       'bulk and event momentum are combined before cell division', failures)
  call expect_close(increment%event_energy_density, expected_event_energy_density, &
       'event energy is coupled to the total-energy increment', failures)
  call expect_close(increment%bulk_kinetic_energy_density, &
       expected_bulk_energy_density, &
       'returned-mass bulk kinetic energy is included', failures)
  call expect_close(increment%total_energy_density, &
       expected_event_energy_density + expected_bulk_energy_density, &
       'total-energy increment closes event plus bulk energy', failures)
  call expect_close(increment%element_mass_density(elem_c), &
       solar_mass_cgs / volume_cm3, &
       'tracked SNIa ejecta is converted to element mass density', failures)
  call expect_close(increment%total_metal_density, &
       1.3_stellar_dp * solar_mass_cgs / volume_cm3, &
       'total-metal SNIa density follows actual returned ejecta', failures)

  call build_snia_cell_increment(budget, 0.0_stellar_dp, bulk_velocity_cm_s, &
       coupling, increment, ierr)
  call expect(ierr == snia_deposition_err_argument .and. &
       increment%mass_density == 0.0_stellar_dp .and. &
       increment%total_energy_density == 0.0_stellar_dp, &
       'non-positive cell volume is rejected transactionally', failures)

  budget%returned_mass = 1.5_stellar_dp
  call build_snia_cell_increment(budget, volume_cm3, bulk_velocity_cm_s, &
       coupling, increment, ierr)
  call expect(ierr == snia_deposition_err_budget .and. &
       increment%mass_density == 0.0_stellar_dp, &
       'budget mass exceeding the WD debit is rejected', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_CELL_DEPOSITION_TEST_FAILED failures=', failures
     stop 1
  end if
  write(*, '(a)') 'FP2_SNIa_CELL_DEPOSITION_TEST_OK'

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

end program fp2_snia_cell_deposition_test

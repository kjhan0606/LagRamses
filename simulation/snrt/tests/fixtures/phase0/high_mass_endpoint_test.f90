program high_mass_endpoint_test
  use stellar_enrichment_config, only: stellar_dp
  use stellar_yield_audit, only: high_mass_endpoint_t, resolve_high_mass_endpoint, yield_audit_ok
  implicit none
  type(high_mass_endpoint_t) :: raw, result
  real(stellar_dp) :: adjustment
  integer :: ierr

  ! Synthetic complete endpoint; energies are non-bulk source energy, as
  ! required by the native bridge. Directed kinetic energy is separate.
  raw%initial_mass = 60
  raw%wind_mass = 10
  raw%terminal_mass = 20
  raw%remnant_mass = 30
  raw%wind_elements(1) = 9
  raw%terminal_elements(6) = 18
  raw%wind_energy = 1e48_stellar_dp
  raw%terminal_energy = 1e51_stellar_dp
  call resolve_high_mass_endpoint(raw, 'source_consistent', .true., 0d0, result, adjustment, ierr)
  if (ierr /= yield_audit_ok .or. adjustment /= 0d0 .or. result%terminal_mass /= 20d0) stop 1
  call resolve_high_mass_endpoint(raw, 'source_consistent', .false., 0d0, result, adjustment, ierr)
  if (ierr == yield_audit_ok .or. result%initial_mass /= 0d0) stop 2
  raw%remnant_mass = 29.5d0
  call resolve_high_mass_endpoint(raw, 'source_consistent', .true., 0d0, result, adjustment, ierr)
  if (ierr == yield_audit_ok) stop 3
  call resolve_high_mass_endpoint(raw, 'mixed_remnant', .false., .01d0, result, adjustment, ierr)
  if (ierr /= yield_audit_ok .or. adjustment /= .5d0 .or. result%remnant_mass /= 30d0) stop 4
  if (any(result%wind_elements /= raw%wind_elements) .or. &
      any(result%terminal_elements /= raw%terminal_elements)) stop 5
  call resolve_high_mass_endpoint(raw, 'mixed_remnant', .false., .001d0, result, adjustment, ierr)
  if (ierr == yield_audit_ok .or. adjustment /= 0d0 .or. result%wind_mass /= 0d0) stop 6
  call resolve_high_mass_endpoint(raw, 'wind_only_collapse', .false., 0d0, result, adjustment, ierr)
  if (ierr /= yield_audit_ok .or. result%remnant_mass /= 50d0 .or. result%terminal_energy /= 0d0) stop 7
  if (any(result%terminal_elements /= 0d0) .or. result%wind_energy /= raw%wind_energy) stop 8
  raw%terminal_elements(1) = 10
  call resolve_high_mass_endpoint(raw, 'wind_only_collapse', .false., 0d0, result, adjustment, ierr)
  if (ierr == yield_audit_ok) stop 9 ! Invalid raw elements are not hidden by selecting collapse.
  raw%terminal_elements(1) = 0
  raw%wind_mass = 61
  call resolve_high_mass_endpoint(raw, 'wind_only_collapse', .false., 0d0, result, adjustment, ierr)
  if (ierr == yield_audit_ok) stop 10

  ! Actual LC18 integrated wind total for 60 Msun, [Fe/H]=-3, vrot=150.
  ! This checks ONLY the selected mass budget, not missing elemental/time/energy inputs.
  raw = high_mass_endpoint_t()
  raw%initial_mass = 60
  raw%wind_mass = 18.06491167124404d0
  call resolve_high_mass_endpoint(raw, 'wind_only_collapse', .false., 0d0, result, adjustment, ierr)
  if (ierr /= yield_audit_ok .or. abs(result%remnant_mass-41.93508832875596d0) > 1d-12) stop 11
  print *, 'HIGH_MASS_ENDPOINT_TEST_OK'
end program high_mass_endpoint_test

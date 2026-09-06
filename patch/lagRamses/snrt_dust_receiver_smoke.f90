program snrt_dust_receiver_smoke
  use snrt_dust_receiver
  implicit none

  integer, parameter :: nc = 2, ng = 3
  real(snrt_dust_receiver_dp) :: edges(ng + 1), sigma(ng), mean_energy(ng)
  real(snrt_dust_receiver_dp) :: metallicity(nc), dust_to_metal(nc), abundance(nc), bad_abundance(nc)
  real(snrt_dust_receiver_dp) :: n_h(nc), path(nc), tau(nc,ng)
  real(snrt_dust_receiver_dp) :: photons(ng,nc), capacity(nc), old_energy(nc)
  real(snrt_dust_receiver_dp) :: old_temperature(nc), staged_energy(nc)
  real(snrt_dust_receiver_dp) :: staged_temperature(nc), absorbed_energy(nc)
  real(snrt_dust_receiver_dp) :: state_energy(nc), state_temperature(nc)
  real(snrt_dust_receiver_dp) :: state_energy_before(nc), state_temperature_before(nc)
  real(snrt_dust_receiver_dp) :: expected_absorbed(nc)
  integer :: ierr
  character(len=64) :: hash

  edges = [1.0d0, 10.0d0, 100.0d0, 1000.0d0]
  sigma = [1.0d-21, 2.0d-21, 0.0d0]
  mean_energy = [3.0d0, 20.0d0, 200.0d0]
  hash = repeat('a', 64)
  call snrt_dust_validate_opacity_binding(edges, sigma, mean_energy, &
       'candidate_source_sed_matched', 'draine_source_sed', hash, hash, ierr)
  if (ierr /= snrt_dust_receiver_ok) error stop 1
  call snrt_dust_validate_opacity_binding(edges, sigma, mean_energy, &
       'unapproved_status', 'draine_source_sed', hash, hash, ierr)
  if (ierr /= snrt_dust_receiver_err_binding) error stop 2
  call snrt_dust_validate_opacity_binding(edges, sigma, mean_energy, &
       'candidate_source_sed_matched', 'draine_source_sed', 'bad', hash, ierr)
  if (ierr /= snrt_dust_receiver_err_binding) error stop 3

  metallicity = [0.5d0, 2.0d0]
  dust_to_metal = [0.1d0, 0.2d0]
  call snrt_dust_map_cell_abundance(metallicity, dust_to_metal, abundance, ierr)
  if (ierr /= snrt_dust_receiver_ok .or. any(abs(abundance - [0.05d0, 0.4d0]) > 1.0d-15)) error stop 4
  call snrt_dust_map_cell_abundance([-1.0d0, 1.0d0], dust_to_metal, bad_abundance, ierr)
  if (ierr /= snrt_dust_receiver_err_input) error stop 5

  n_h = [1.0d3, 2.0d3]
  path = [1.0d0, 2.0d0]
  call snrt_dust_prepare_cell_optical_depth(n_h, path, abundance, sigma, tau, ierr)
  if (ierr /= snrt_dust_receiver_ok .or. abs(tau(1,1) - 5.0d-20) > 1.0d-32 .or. &
       abs(tau(2,2) - 3.2d-18) > 1.0d-30 .or. any(tau(:,3) /= 0.0d0)) error stop 6

  photons = reshape([1.0d12, 2.0d12, 0.0d0, 0.0d0, 1.0d12, 3.0d12], [ng,nc])
  capacity = [2.0d0, 4.0d0]
  old_energy = [10.0d0, 20.0d0]
  old_temperature = [20.0d0, 30.0d0]
  state_energy = old_energy
  state_temperature = old_temperature
  state_energy_before = state_energy
  state_temperature_before = state_temperature
  expected_absorbed = [(1.0d12*3.0d0 + 2.0d12*20.0d0) * &
       snrt_dust_receiver_ev_to_erg, (1.0d12*20.0d0 + 3.0d12*200.0d0) * &
       snrt_dust_receiver_ev_to_erg]
  call snrt_dust_receiver_stage(photons, mean_energy, 10.0d0, abundance, capacity, &
       old_energy, old_temperature, staged_energy, staged_temperature, absorbed_energy, ierr)
  if (ierr /= snrt_dust_receiver_ok .or. any(abs(absorbed_energy - expected_absorbed) > 1.0d-12) .or. &
       any(abs(staged_energy - old_energy - absorbed_energy) > 1.0d-12) .or. &
       any(abs(staged_temperature - old_temperature - absorbed_energy/capacity) > 1.0d-12)) error stop 7
  if (any(state_energy /= state_energy_before) .or. any(state_temperature /= state_temperature_before)) error stop 8
  call snrt_dust_receiver_commit(state_energy, state_temperature, staged_energy, staged_temperature, ierr)
  if (ierr /= snrt_dust_receiver_ok .or. any(state_energy /= staged_energy) .or. &
       any(state_temperature /= staged_temperature)) error stop 9

  ! Positive dust absorption must never be accepted without a positive cell
  ! abundance.  The persistent state remains unchanged on this failed trial.
  state_energy_before = state_energy
  state_temperature_before = state_temperature
  call snrt_dust_receiver_stage(photons, mean_energy, 10.0d0, [0.0d0, abundance(2)], capacity, &
       state_energy, state_temperature, staged_energy, staged_temperature, absorbed_energy, ierr)
  if (ierr /= snrt_dust_receiver_err_state .or. any(state_energy /= state_energy_before) .or. &
       any(state_temperature /= state_temperature_before)) error stop 10

  ! Zero photons are a valid zero-dust no-op; malformed candidates cannot
  ! enter the persistent state.
  photons = 0.0d0
  call snrt_dust_receiver_stage(photons, mean_energy, 10.0d0, [0.0d0, 0.0d0], capacity, &
       [0.0d0, 0.0d0], old_temperature, staged_energy, staged_temperature, absorbed_energy, ierr)
  if (ierr /= snrt_dust_receiver_ok .or. any(staged_energy /= 0.0d0) .or. &
       any(staged_temperature /= old_temperature)) error stop 11
  call snrt_dust_receiver_commit(state_energy, state_temperature, staged_energy, &
       [-1.0d0, staged_temperature(2)], ierr)
  if (ierr /= snrt_dust_receiver_err_state .or. any(state_energy /= state_energy_before) .or. &
       any(state_temperature /= state_temperature_before)) error stop 12

  write(*,'(a)') 'SNRT_NATIVE_DUST_MAPPING_RECEIVER_OK binding=1 mapping=1 opacity=1 thermal=1 closure=1 rollback=1'
end program snrt_dust_receiver_smoke

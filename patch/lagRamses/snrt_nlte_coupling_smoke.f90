program snrt_nlte_coupling_smoke
  use amr_parameters, only: dp
  use snrt_nlte_coupling, only: snrt_nlte_optical_depth, snrt_nlte_photo_source
  implicit none

  real(dp) :: tau, ionization_increment, heating_rate
  real(dp) :: expected_tau, expected_heat

  call snrt_nlte_optical_depth(0.5d0, 1.0d-3, 1.0d12, 1.0d-18, tau)
  expected_tau = 2.99792458d10 * 1.0d12 * 0.5d0 * 1.0d-3 * 1.0d-18
  if (abs(tau-expected_tau) / expected_tau > 1.0d-13) error stop 1

  call snrt_nlte_photo_source(0.25d0, 1.0d0, 0.5d0, 1.0d-3, 1.0d12, 20.0d0, &
       ionization_increment, heating_rate)
  expected_heat = 0.25d0 * 1.0d-3 * 6.4d0 * 1.602176634d-12 / 1.0d12
  if (abs(ionization_increment-0.25d0) > 1.0d-13) error stop 2
  if (abs(heating_rate-expected_heat) / expected_heat > 1.0d-13) error stop 3

  write(*,'(a,es14.6,a,es14.6)') 'SNRT_NLTE_COUPLING_OK tau=', tau, &
       ' heating_rate=', heating_rate
end program snrt_nlte_coupling_smoke

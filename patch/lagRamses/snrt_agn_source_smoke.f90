program snrt_agn_source_smoke
  use, intrinsic :: iso_c_binding, only: c_float
  use amr_parameters, only: dp
  use snrt_agn_source, only: snrt_agn_photon_budget, snrt_agn_isotropic_packet, &
       snrt_agn_photons_to_density_code, snrt_agn_deposit_isotropic, &
       snrt_agn_deposit_transaction, snrt_c_cgs, snrt_ev_to_erg, &
       snrt_agn_photon_budget_energy, snrt_agn_source_commit
  implicit none

  real(dp) :: luminosity_erg_s, emitted_photons
  real(dp) :: expected_energy_erg, expected_photons
  real(dp) :: angular_weights(4), directional_photons(4)
  real(dp) :: photon_density_code, transaction_deposited
  real(dp) :: emitted_groups(2)
  real(dp) :: pending,reference_photons
  real(c_float) :: state(4,2,1), state_before(4,2,1)
  integer :: ierr

  ! The first positional argument is supplied inflow, not retained BH mass.
  call snrt_agn_photon_budget(2.0d0, 5.0d33, 4.0d0, 0.1d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)

  expected_energy_erg = 0.1d0 * 2.0d0 * 5.0d33 * snrt_c_cgs**2
  expected_photons = 0.25d0 * expected_energy_erg / (20.0d0 * snrt_ev_to_erg)
  if (abs(luminosity_erg_s - expected_energy_erg / 4.0d0) / &
       (expected_energy_erg / 4.0d0) > 1.0d-13) error stop 1
  if (abs(emitted_photons - expected_photons) / expected_photons > 1.0d-13) error stop 2
  call snrt_agn_photon_budget_energy(expected_energy_erg,4d0,0.25d0,20d0, &
       luminosity_erg_s,emitted_photons,ierr)
  if(ierr/=0.or.abs(emitted_photons/expected_photons-1d0)>1d-13)error stop 15
  call snrt_agn_photon_budget_energy(0d0,4d0,0.25d0,20d0,luminosity_erg_s,emitted_photons,ierr)
  if(ierr/=0.or.luminosity_erg_s/=0d0.or.emitted_photons/=0d0)error stop 16
  call snrt_agn_photon_budget_energy(-1d0,4d0,0.25d0,20d0,luminosity_erg_s,emitted_photons,ierr)
  if(ierr==0)error stop 17

  call snrt_agn_photon_budget(0.0d0, 5.0d33, 4.0d0, 0.1d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)
  if (luminosity_erg_s /= 0.0d0 .or. emitted_photons /= 0.0d0) error stop 3

  call snrt_agn_photon_budget(2.0d0, 5.0d33, 4.0d0, 1.0d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)
  if (luminosity_erg_s /= 0.0d0 .or. emitted_photons /= 0.0d0) error stop 14

  angular_weights = (/1.0d0, 2.0d0, -1.0d0, 1.0d0/)
  call snrt_agn_isotropic_packet(12.0d0, angular_weights, directional_photons)
  if (abs(sum(directional_photons) - 12.0d0) > 1.0d-13) error stop 4
  if (maxval(abs(directional_photons - (/3.0d0, 6.0d0, 0.0d0, 3.0d0/))) > &
       1.0d-13) error stop 5

  ! Production deposition requires a non-negative quadrature.  Keep the
  ! negative-weight compatibility check above local to the pure splitter.
  angular_weights = (/1.0d0, 2.0d0, 1.0d0, 1.0d0/)

  call snrt_agn_photons_to_density_code(1.0d63, 125.0d0, 1.0d21, 1.0d-3, &
       photon_density_code)
  if (abs(photon_density_code - 8.0d0) > 1.0d-13) error stop 6

  state = 0.0_c_float
  call snrt_agn_deposit_isotropic(state, 1, 1, 1.0d63, 125.0d0, 1.0d21, &
       1.0d-3, angular_weights, photon_density_code, ierr)
  if (ierr /= 0) error stop 7
  if (abs(sum(real(state,dp)) - 8.0d0) > 1.0d-6) error stop 8

  emitted_groups = (/1.0d63, 2.0d63/)
  call snrt_agn_deposit_transaction(state, 1, emitted_groups, 125.0d0, &
       1.0d21, 1.0d-3, angular_weights, transaction_deposited, ierr)
  if (ierr /= 0) error stop 9
  if (abs(transaction_deposited - 24.0d0) > 1.0d-12) error stop 10
  if (abs(sum(real(state,dp)) - 32.0d0) > 1.0d-6) error stop 11

  state_before = state
  emitted_groups = (/1.0d63, -1.0d0/)
  call snrt_agn_deposit_transaction(state, 1, emitted_groups, 125.0d0, &
       1.0d21, 1.0d-3, angular_weights, transaction_deposited, ierr)
  if (ierr == 0) error stop 12
  if (maxval(abs(real(state,dp) - real(state_before,dp))) > 0.0d0) error stop 13

  ! Accepted-energy receipt survives a failed all-group source transaction.
  pending=1d53; state=0
  call snrt_agn_photon_budget_energy(pending,4d0,0.25d0,20d0,luminosity_erg_s,emitted_groups(1),ierr)
  reference_photons=emitted_groups(1)
  emitted_groups(2)=-1d0
  call snrt_agn_deposit_transaction(state,1,emitted_groups,125d0,1d21,1d-3,angular_weights, &
       transaction_deposited,ierr)
  call snrt_agn_source_commit(pending,ierr==0)
  if(ierr==0.or.pending/=1d53.or.any(state/=0))error stop 18
  ! Retry succeeds. Unrepresented bolometric energy is not renormalized.
  emitted_groups(2)=reference_photons
  call snrt_agn_deposit_transaction(state,1,emitted_groups,125d0,1d21,1d-3,angular_weights, &
       transaction_deposited,ierr)
  call snrt_agn_source_commit(pending,ierr==0)
  if(ierr/=0.or.pending/=0d0)error stop 19
  state_before=state ! The existing RT rollback snapshot is AFTER source commit.
  state=state_before
  call snrt_agn_photon_budget_energy(pending,4d0,0.25d0,20d0,luminosity_erg_s,emitted_photons,ierr)
  if(ierr/=0.or.emitted_photons/=0d0.or.any(state/=state_before))error stop 20
  write(*,'(A)') 'SNRT_ACCEPTED_ENERGY_SOURCE_RETRY_OK'

  write(*,'(a,es14.6,a,es14.6)') 'SNRT_AGN_SOURCE_OK luminosity=', &
       expected_energy_erg / 4.0d0, ' photons=', expected_photons
end program snrt_agn_source_smoke

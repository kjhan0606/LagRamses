program snrt_spectral_contract_smoke
  use amr_parameters, only: dp
  use snrt_agn_source, only: snrt_agn_cgs => snrt_c_cgs, &
       snrt_agn_photon_budget
  use snrt_spectral_contract, only: snrt_ngroups, snrt_nedges, &
       snrt_group_edges_ev, snrt_group_mean_energy_ev, &
       snrt_group_energy_fraction, snrt_group_cross_section_cm2, &
       snrt_group_cross_section_hei_cm2, snrt_group_cross_section_heii_cm2, &
       snrt_group_photoelectron_excess_energy_ev, &
       snrt_group_photoelectron_excess_hei_ev, &
       snrt_group_photoelectron_excess_heii_ev, &
       snrt_group_energy_fraction_sum, snrt_group_unrepresented_energy_fraction, &
       snrt_spectral_contract_load_from_environment, &
       snrt_spectral_contract_validate_values, &
       snrt_spectral_contract_checkpoint_identity_matches, &
       snrt_spectral_contract_loaded, snrt_spectral_contract_runtime_allowed, &
       snrt_spectral_contract_error_name, &
       snrt_spectral_contract_source_id, snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, &
       snrt_spectral_contract_approval_id, snrt_spectral_contract_status, &
       snrt_spectral_contract_group_edges_sha256, &
       snrt_spectral_contract_interval_convention, &
       snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_error_message, &
       snrt_spectral_contract_ok, snrt_spectral_contract_err_edges, &
       snrt_spectral_contract_err_values
  implicit none

  real(dp) :: edges(snrt_nedges), mean_energy(snrt_ngroups)
  real(dp) :: fraction(snrt_ngroups), hi(snrt_ngroups), hei(snrt_ngroups)
  real(dp) :: heii(snrt_ngroups), excess_hi(snrt_ngroups)
  real(dp) :: excess_hei(snrt_ngroups), excess_heii(snrt_ngroups)
  real(dp) :: luminosity, emitted_photons, sum_group_energy, expected_luminosity
  character(len=128) :: mismatched_edges_sha256
  integer :: ierr, failures

  failures = 0
  call snrt_spectral_contract_load_from_environment(ierr)
  write(*,'(a,i0,a,a)') 'loader_status=', ierr, ' (', &
       trim(snrt_spectral_contract_error_name(ierr))//')'
  if (ierr /= snrt_spectral_contract_ok) write(*,'(a,a)') &
       'loader_message=', trim(snrt_spectral_contract_error_message)
  call expect(ierr == snrt_spectral_contract_ok, &
       'reference contract loads from SNRT_GROUP_CONTRACT', failures)
  call expect(snrt_spectral_contract_loaded .and. &
       snrt_spectral_contract_runtime_allowed, &
       'reference contract is runtime-admissible as a control', failures)
  call expect(snrt_ngroups == 9 .and. snrt_nedges == 10, &
       'native state uses nine groups and ten boundaries', failures)
  call expect(maxval(abs(snrt_group_edges_ev - &
       (/0.01d0,1.0d0,5.6d0,11.2d0,13.6d0,24.59d0,54.42d0, &
         500.0d0,2000.0d0,10000.0d0/))) < 1.0d-12, &
       'canonical group edges are exact at native boundary', failures)
  call expect(snrt_group_energy_fraction_sum > 0.0d0 .and. &
       snrt_group_energy_fraction_sum < 1.0d0 .and. &
       abs(snrt_group_energy_fraction_sum + &
       snrt_group_unrepresented_energy_fraction - 1.0d0) < 1.0d-14, &
       'represented and unrepresented source fractions close', failures)
  call expect(snrt_group_cross_section_cm2(1) == 0.0d0 .and. &
       snrt_group_cross_section_cm2(5) > 0.0d0 .and. &
       snrt_group_cross_section_hei_cm2(5) == 0.0d0 .and. &
       snrt_group_cross_section_hei_cm2(6) > 0.0d0 .and. &
       snrt_group_cross_section_heii_cm2(6) == 0.0d0 .and. &
       snrt_group_cross_section_heii_cm2(7) > 0.0d0, &
       'threshold-safe H/He cross-section support is retained', failures)
  call expect(snrt_spectral_contract_checkpoint_identity_matches( &
       snrt_spectral_contract_source_id, snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, &
       snrt_spectral_contract_approval_id, &
       snrt_spectral_contract_group_edges_sha256, &
       snrt_spectral_contract_interval_convention, &
       snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_status), &
       'checkpoint identity matches the loaded spectral contract', failures)
  mismatched_edges_sha256 = snrt_spectral_contract_group_edges_sha256
  mismatched_edges_sha256(1:1) = '0'
  call expect(.not. snrt_spectral_contract_checkpoint_identity_matches( &
       snrt_spectral_contract_source_id, snrt_spectral_contract_source_sha256, &
       snrt_spectral_contract_source_commit_binding, &
       snrt_spectral_contract_approval_id, mismatched_edges_sha256, &
       snrt_spectral_contract_interval_convention, &
       snrt_spectral_contract_fraction_semantics, &
       snrt_spectral_contract_status), &
       'edge digest mismatch is rejected by checkpoint identity binding', failures)

  sum_group_energy = 0.0d0
  do ierr = 1, snrt_ngroups
     call snrt_agn_photon_budget(1.0d0, 1.0d0, 1.0d0, 0.1d0, &
          snrt_group_energy_fraction(ierr), snrt_group_mean_energy_ev(ierr), &
          luminosity, emitted_photons)
     sum_group_energy = sum_group_energy + emitted_photons * &
          snrt_group_mean_energy_ev(ierr) * 1.602176634d-12
  end do
  expected_luminosity = 0.1d0 * snrt_agn_cgs**2 * &
       snrt_group_energy_fraction_sum
  call expect(abs(sum_group_energy-expected_luminosity) / expected_luminosity < 1.0d-13, &
       'source groups conserve the represented Lbol fraction without hidden rescaling', failures)

  edges = snrt_group_edges_ev
  mean_energy = snrt_group_mean_energy_ev
  fraction = snrt_group_energy_fraction
  hi = snrt_group_cross_section_cm2
  hei = snrt_group_cross_section_hei_cm2
  heii = snrt_group_cross_section_heii_cm2
  excess_hi = snrt_group_photoelectron_excess_energy_ev
  excess_hei = snrt_group_photoelectron_excess_hei_ev
  excess_heii = snrt_group_photoelectron_excess_heii_ev

  edges(7) = edges(7) + 1.0d0
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_edges, &
       'canonical edge mismatch is rejected', failures)
  edges = snrt_group_edges_ev

  hi(1) = 1.0d-20
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_values, &
       'sub-threshold H opacity is rejected', failures)
  hi = snrt_group_cross_section_cm2

  mean_energy(5) = 100.0d0
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_values, &
       'out-of-band group representative energy is rejected', failures)

  mean_energy = snrt_group_mean_energy_ev
  excess_hi(5) = 100.0d0
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_values, &
       'excess-energy upper bound is rejected', failures)
  excess_hi = snrt_group_photoelectron_excess_energy_ev

  hi(5) = 1.0d-15
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_values, &
       'cross-section upper bound is rejected', failures)
  hi = snrt_group_cross_section_cm2
  hi(5) = 3.5d-17
  call snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, &
       hei, heii, excess_hi, excess_hei, excess_heii, ierr)
  call expect(ierr == snrt_spectral_contract_err_values, &
       'one-decade cross-section transcription slip is rejected', failures)

  if (failures == 0) then
     write(*,'(a)') 'SNRT_SPECTRAL_CONTRACT_OK'
  else
     write(*,'(a,i0)') 'SNRT_SPECTRAL_CONTRACT_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (condition) then
       write(*,'(a)') 'PASS: ' // trim(label)
    else
       failures = failures + 1
       write(*,'(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine expect

end program snrt_spectral_contract_smoke

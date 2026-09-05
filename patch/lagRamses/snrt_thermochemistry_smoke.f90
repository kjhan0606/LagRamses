program snrt_thermochemistry_smoke
  use amr_parameters, only: dp
  use snrt_thermochemistry, only: &
       snrt_secondary_tables_load_from_environment, snrt_secondary_tables_loaded, &
       snrt_secondary_source_id, snrt_secondary_upstream_commit, &
       snrt_secondary_manifest_sha256, snrt_secondary_fractions, &
       snrt_alpha_hydrogen_case_b, snrt_alpha_helium_ii_case_b, &
       snrt_alpha_helium_ii_radiative_case_b, &
       snrt_alpha_helium_ii_dielectronic_case_b, snrt_alpha_helium_iii_case_b, &
       snrt_partition_absorption, &
       snrt_thermochemistry_advance_cell, snrt_thermochemistry_result, &
       snrt_thermochemistry_ok, snrt_thermochemistry_err_inventory
  implicit none

  integer, parameter :: ngroup = 9
  real(dp), parameter :: expected_fion_200 = 0.09607314283559579d0
  ! The public table row has a 9.41e-6 rounding residual; the native
  ! contract normalizes all five deposition channels before returning them.
  real(dp), parameter :: expected_fheat_200 = 0.8050044271114285d0
  real(dp), parameter :: expected_fexc_200 = 0.09892333418982854d0
  real(dp) :: fheat, fhi, fhei, fheii, fexc, fion
  real(dp) :: fheat_low, fhi_low, fhei_low, fheii_low, fexc_low
  real(dp) :: fheat_high, fhi_high, fhei_high, fheii_high, fexc_high
  real(dp) :: max_delta, max_floor_delta
  real(dp) :: alpha_h, alpha_heii, alpha_heiii, alpha_heii_hot
  real(dp) :: alpha_heii_rad_hot, alpha_heii_dielectronic_hot
  real(dp) :: opacity(3), available(3), partition(3)
  real(dp) :: unassigned
  real(dp) :: absorbed(3,ngroup), excess(3,ngroup)
  real(dp) :: simplex
  type(snrt_thermochemistry_result) :: result
  integer :: ierr, failures

  failures = 0
  call snrt_secondary_tables_load_from_environment(ierr)
  call expect(ierr == snrt_thermochemistry_ok .and. snrt_secondary_tables_loaded, &
       'native FS2010 contract and all fourteen tables load', failures)
  call expect(trim(snrt_secondary_source_id) == &
       'furlanetto_stoever_2010_21cmfast', &
       'native table source identity is pinned', failures)
  call expect(len_trim(snrt_secondary_upstream_commit) == 40 .and. &
       len_trim(snrt_secondary_manifest_sha256) == 64, &
       'native table contract exposes upstream and manifest identities', failures)

  call snrt_secondary_fractions(200.0d0, 0.1d0, fheat, fhi, fhei, fheii, &
       fexc, ierr, fion)
  call expect(ierr == 0, 'native FS2010 interpolation returns successfully', failures)
  call expect(abs(fion-expected_fion_200) < 2.0d-12 .and. &
       abs(fheat-expected_fheat_200) < 2.0d-12 .and. &
       abs(fexc-expected_fexc_200) < 2.0d-12, &
       'native 200 eV, xHII=0.1 interpolation matches pinned reference', failures)
  call expect(abs(fheat+fhi+fhei+fheii+fexc-1.0d0) < 2.0d-15 .and. &
       min(fheat,fhi,fhei,fheii,fexc) >= 0.0d0, &
       'native deposition fractions close and remain non-negative', failures)

  call snrt_secondary_fractions(99.9d0, 0.1d0, fheat_low, fhi_low, fhei_low, &
       fheii_low, fexc_low, ierr)
  call snrt_secondary_fractions(100.1d0, 0.1d0, fheat_high, fhi_high, fhei_high, &
       fheii_high, fexc_high, ierr)
  max_delta = maxval(abs((/fheat_high-fheat_low, fhi_high-fhi_low, &
       fhei_high-fhei_low, fheii_high-fheii_low, fexc_high-fexc_low/)))
  call expect(max_delta < 5.0d-3, 'native 99.9/100.1 eV interpolation is continuous', failures)
  call snrt_secondary_fractions(9.999d0, 0.1d0, fheat_low, fhi_low, fhei_low, &
       fheii_low, fexc_low, ierr)
  call snrt_secondary_fractions(10.001d0, 0.1d0, fheat_high, fhi_high, fhei_high, &
       fheii_high, fexc_high, ierr)
  max_floor_delta = maxval(abs((/fheat_high-fheat_low, fhi_high-fhi_low, &
       fhei_high-fhei_low, fheii_high-fheii_low, fexc_high-fexc_low/)))
  call expect(max_floor_delta < 5.0d-3, 'native 10 eV table-floor transition is bounded', failures)
  call snrt_secondary_fractions(5.0d0, 0.1d0, fheat, fhi, fhei, fheii, fexc, ierr)
  call expect(ierr == 0 .and. fheat == 1.0d0 .and. fhi == 0.0d0 .and. &
       fhei == 0.0d0 .and. fheii == 0.0d0 .and. fexc == 0.0d0, &
       'below-table electron energy is assigned entirely to heat', failures)

  alpha_h = snrt_alpha_hydrogen_case_b(10000.0d0)
  alpha_heii = snrt_alpha_helium_ii_case_b(10000.0d0)
  alpha_heiii = snrt_alpha_helium_iii_case_b(10000.0d0)
  call expect(alpha_h > 0.0d0 .and. alpha_heii > 0.0d0 .and. alpha_heiii > 0.0d0, &
       'native case-B recombination coefficients are positive', failures)
  call expect(abs(alpha_heiii-2.0d0*snrt_alpha_hydrogen_case_b(2500.0d0)) < 1.0d-24, &
       'He III uses exactly 2 alpha_H,B(T/4)', failures)
  call expect(abs(alpha_heii-2.616130035d-13)/2.616130035d-13 < 2.0d-6, &
       'He II case-B coefficient matches the temperature-resolved reference', failures)
  alpha_heii_hot = snrt_alpha_helium_ii_case_b(100000.0d0)
  alpha_heii_rad_hot = snrt_alpha_helium_ii_radiative_case_b(100000.0d0)
  alpha_heii_dielectronic_hot = snrt_alpha_helium_ii_dielectronic_case_b(100000.0d0)
  call expect(abs(alpha_heii_hot-alpha_heii_rad_hot-alpha_heii_dielectronic_hot) / &
       alpha_heii_hot < 1.0d-12 .and. alpha_heii_dielectronic_hot > 0.0d0 .and. &
       alpha_heii_dielectronic_hot/alpha_heii_hot > 1.0d-3, &
       'He II case-B retains the non-negligible dielectronic term at 1e5 K', failures)

  opacity = (/1.0d0, 10.0d0, 1.0d0/)
  available = (/0.10d0, 0.02d0, 0.01d0/)
  call snrt_partition_absorption(0.12d0, opacity, available, partition, ierr)
  call expect(ierr == 0 .and. abs(sum(partition)-0.12d0) < 1.0d-12 .and. &
       minval(available) >= -1.0d-14, &
       'native absorption partition redistributes around species inventory caps', failures)
  opacity = (/1.0d0, 0.0d0, 0.0d0/)
  available = (/0.10d0, 0.02d0, 0.01d0/)
  call snrt_partition_absorption(0.10000001d0, opacity, available, partition, ierr, unassigned)
  call expect(ierr == 0 .and. abs(partition(1)-0.10d0) < 1.0d-12 .and. &
       partition(2) == 0.0d0 .and. partition(3) == 0.0d0 .and. &
       unassigned > 0.0d0, &
       'partition redistribution never assigns a group to an opaque-zero species', failures)
  opacity = (/1.0d0, 0.0d0, 0.0d0/)
  available = (/1.0d-10, 0.0d0, 0.0d0/)
  call snrt_partition_absorption(1.00001d-10, opacity, available, partition, ierr, &
       unassigned, inventory_scale_code=1.0d0)
  call expect(ierr == 0 .and. abs(partition(1)-1.0d-10) < 1.0d-24 .and. &
       unassigned > 0.0d0, &
       'partition tolerance remains tied to the pre-partition cell scale', failures)
  available = (/0.10d0, 0.0d0, 0.0d0/)
  call snrt_partition_absorption(0.11d0, opacity, available, partition, ierr, unassigned, &
       inventory_scale_code=0.10d0)
  call expect(ierr == snrt_thermochemistry_err_inventory .and. &
       abs(available(1)-0.10d0) < 1.0d-15, &
       'above-tolerance unassigned absorption is rejected without mutation', failures)

  absorbed = 0.0d0
  excess = 0.0d0
  absorbed(1,7) = 0.05d0
  absorbed(2,7) = 0.002d0
  absorbed(3,7) = 0.001d0
  excess(:,7) = 200.0d0
  call snrt_thermochemistry_advance_cell(1.0d0, 0.0789474d0, 1.0d0, &
       10000.0d0, 1.0d11, 0.10d0, 0.10d0, 0.0d0, absorbed, excess, result)
  call expect(result%ierr == 0, 'native H/He photo-thermochemistry step succeeds', failures)
  simplex = result%x_helium_ii + result%x_helium_iii
  call expect(result%x_hydrogen_ii >= 0.0d0 .and. result%x_hydrogen_ii <= 1.0d0 .and. &
       result%x_helium_ii >= 0.0d0 .and. result%x_helium_iii >= 0.0d0 .and. &
       simplex <= 1.0d0 + 1.0d-12 .and. result%electron_density_cm3 >= 0.0d0, &
       'native H/He fractions remain on their physical simplex', failures)
  call expect(result%primary_hydrogen_ionizations_cm3 > 0.0d0 .and. &
       result%secondary_hydrogen_ionizations_cm3 >= 0.0d0 .and. &
       result%secondary_helium_i_ionizations_cm3 >= 0.0d0 .and. &
       result%secondary_helium_ii_ionizations_cm3 >= 0.0d0 .and. &
       result%recombination_hydrogen_cm3 >= 0.0d0 .and. &
       result%recombination_helium_ii_cm3 >= 0.0d0 .and. &
       result%recombination_helium_iii_cm3 >= 0.0d0, &
       'native primary/secondary/recombination ledgers are non-negative', failures)
  call expect(abs(result%photoelectron_energy_residual_ev_cm3) < &
       1.0d-11*max(1.0d0,result%photoelectron_energy_ev_cm3), &
       'native photoelectron energy closes into heat, ionization, and excitation', failures)
  call expect(result%heating_rate_erg_cm3_s > 0.0d0 .and. &
       result%absorbed_photon_energy_ev_cm3 >= result%photoelectron_energy_ev_cm3, &
       'only the explicitly deposited gas-heating channel is exposed to RAMSES', failures)

  absorbed = 0.0d0
  excess = 0.0d0
  absorbed(2,7) = 0.002d0
  excess(2,7) = 200.0d0
  call snrt_thermochemistry_advance_cell(1.0d0, 0.0789474d0, 1.0d0, &
       10000.0d0, 1.0d11, 1.0d0, 0.0d0, 0.0d0, absorbed, excess, result)
  call expect(result%ierr == 0 .and. result%secondary_hydrogen_ionizations_cm3 == 0.0d0 .and. &
       result%secondary_heating_energy_ev_cm3 > 0.0d0, &
       'secondary ionization unavailable in a saturated H II cell is routed to heat', failures)

  if (failures == 0) then
     write(*,'(a)') 'SNRT_NATIVE_THERMOCHEMISTRY_OK'
  else
     write(*,'(a,i0)') 'SNRT_NATIVE_THERMOCHEMISTRY_FAIL count=', failures
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

end program snrt_thermochemistry_smoke

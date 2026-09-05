! Native primordial photo-thermochemistry for the SNRT high-level feedback path.
!
! This module is deliberately independent of Python/JAX.  It loads the pinned
! Furlanetto--Stoever (2010) tables, interpolates their H II-dependent energy
! fractions, and advances a local H/He photoionization plus case-B
! recombination state.  It does not implement the full RAMSES cooling model;
! background, metal, Compton, and recombination-line cooling remain separate
! ledgers until their receivers are approved.
module snrt_thermochemistry
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: dp
  use snrt_agn_source, only: snrt_ev_to_erg
  implicit none

  private

  integer, parameter, public :: snrt_secondary_nenergy = 258
  integer, parameter, public :: snrt_secondary_nxi = 14
  integer, parameter, public :: snrt_thermochemistry_nspecies = 3
  integer, parameter, public :: snrt_thermochemistry_ngroups = 9

  real(dp), parameter, public :: snrt_hydrogen_i_threshold_ev = 13.60d0
  real(dp), parameter, public :: snrt_helium_i_threshold_ev = 24.59d0
  real(dp), parameter, public :: snrt_helium_ii_threshold_ev = 54.42d0
  real(dp), parameter, public :: snrt_nhelium_per_hydrogen = 0.0789474d0

  integer, parameter, public :: snrt_thermochemistry_ok = 0
  integer, parameter, public :: snrt_thermochemistry_err_missing = 1
  integer, parameter, public :: snrt_thermochemistry_err_open = 2
  integer, parameter, public :: snrt_thermochemistry_err_read = 3
  integer, parameter, public :: snrt_thermochemistry_err_identity = 4
  integer, parameter, public :: snrt_thermochemistry_err_table_open = 5
  integer, parameter, public :: snrt_thermochemistry_err_table_read = 6
  integer, parameter, public :: snrt_thermochemistry_err_grid = 7
  integer, parameter, public :: snrt_thermochemistry_err_values = 8
  integer, parameter, public :: snrt_thermochemistry_err_input = 9
  integer, parameter, public :: snrt_thermochemistry_err_not_loaded = 10
  integer, parameter, public :: snrt_thermochemistry_err_inventory = 11
  integer, parameter, public :: snrt_thermochemistry_err_energy = 12

  character(len=*), parameter, public :: snrt_secondary_source_id = &
       'furlanetto_stoever_2010_21cmfast'
  character(len=*), parameter, public :: snrt_secondary_upstream_commit = &
       '892f98c80cfe985ca6b399ec6b51a3aa95124b11'
  character(len=*), parameter, public :: snrt_secondary_manifest_sha256 = &
       'c74610c82414504263868704d8d5e913d7a46815c9803e30053780cd8e61d2a3'

  real(dp), parameter :: snrt_secondary_energy_minimum_ev = 10.0d0
  real(dp), parameter :: snrt_secondary_energy_maximum_ev = 9937.21d0
  real(dp), parameter :: snrt_boltzmann_erg_k = 1.380649d-16
  real(dp), parameter :: snrt_mass_fraction_hydrogen = 0.76d0
  real(dp), parameter :: snrt_mass_fraction_helium = 0.24d0
  real(dp), parameter :: snrt_energy_tolerance = 1.0d-11
  ! Transport inventories cross the CUDA FP32 boundary.  The bound is derived
  ! from the observed directional FP32 reduction margin rather than from a
  ! fixed code-density number: 256 single-precision ulps cover the measured
  ! reduction error while remaining below the 0.99995 CUDA guard band.
  real(dp), parameter, public :: snrt_inventory_fp32_ulps = 256.0d0
  real(dp), parameter, public :: snrt_inventory_host_ulps = 8.0d0

  real(dp), parameter :: snrt_xi_grid(snrt_secondary_nxi) = (/ &
       1.0d-4, 2.318d-4, 4.677d-4, 1.0d-3, 2.318d-3, 4.677d-3, &
       1.0d-2, 2.318d-2, 4.677d-2, 1.0d-1, 0.5d0, 0.9d0, 0.99d0, 0.999d0 /)
  character(len=24), parameter :: snrt_table_filename(snrt_secondary_nxi) = [ &
       character(len=24) :: &
       'log_xi_-4.0.dat', 'log_xi_-3.6.dat', 'log_xi_-3.3.dat', &
       'log_xi_-3.0.dat', 'log_xi_-2.6.dat', 'log_xi_-2.3.dat', &
       'log_xi_-2.0.dat', 'log_xi_-1.6.dat', 'log_xi_-1.3.dat', &
       'log_xi_-1.0.dat', 'xi_0.500.dat', 'xi_0.900.dat', &
       'xi_0.990.dat', 'xi_0.999.dat' ]

  real(dp), save :: snrt_secondary_energy_ev(snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_total_ionization(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_heating(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_excitation(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_hydrogen_count(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_helium_i_count(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  real(dp), save :: snrt_secondary_helium_ii_count(snrt_secondary_nxi, &
       snrt_secondary_nenergy) = 0.0d0
  logical, save, public :: snrt_secondary_tables_loaded = .false.
  character(len=128), save, public :: snrt_secondary_loaded_source_id = ''
  character(len=64), save, public :: snrt_secondary_loaded_upstream_commit = ''
  character(len=64), save, public :: snrt_secondary_loaded_manifest_sha256 = ''
  character(len=256), save, public :: snrt_thermochemistry_error_message = ''

  type, public :: snrt_thermochemistry_result
     integer :: ierr = 0
     real(dp) :: x_hydrogen_ii = 0.0d0
     real(dp) :: x_helium_ii = 0.0d0
     real(dp) :: x_helium_iii = 0.0d0
     real(dp) :: electron_density_cm3 = 0.0d0
     real(dp) :: primary_hydrogen_ionizations_cm3 = 0.0d0
     real(dp) :: primary_helium_i_ionizations_cm3 = 0.0d0
     real(dp) :: primary_helium_ii_ionizations_cm3 = 0.0d0
     real(dp) :: secondary_hydrogen_ionizations_cm3 = 0.0d0
     real(dp) :: secondary_helium_i_ionizations_cm3 = 0.0d0
     real(dp) :: secondary_helium_ii_ionizations_cm3 = 0.0d0
     real(dp) :: recombination_hydrogen_cm3 = 0.0d0
     real(dp) :: recombination_helium_ii_cm3 = 0.0d0
     real(dp) :: recombination_helium_iii_cm3 = 0.0d0
     real(dp) :: absorbed_photon_energy_ev_cm3 = 0.0d0
     real(dp) :: primary_ionization_energy_ev_cm3 = 0.0d0
     real(dp) :: photoelectron_energy_ev_cm3 = 0.0d0
     real(dp) :: secondary_heating_energy_ev_cm3 = 0.0d0
     real(dp) :: secondary_ionization_energy_ev_cm3 = 0.0d0
     real(dp) :: excitation_energy_ev_cm3 = 0.0d0
     real(dp) :: photoelectron_energy_residual_ev_cm3 = 0.0d0
     real(dp) :: heating_rate_erg_cm3_s = 0.0d0
  end type snrt_thermochemistry_result

  public :: snrt_secondary_tables_load
  public :: snrt_secondary_tables_load_from_environment
  public :: snrt_secondary_tables_reset
  public :: snrt_secondary_fractions
  public :: snrt_alpha_hydrogen_case_b
  public :: snrt_alpha_helium_ii_radiative_case_b
  public :: snrt_alpha_helium_ii_dielectronic_case_b
  public :: snrt_alpha_helium_ii_case_b
  public :: snrt_alpha_helium_iii_case_b
  public :: snrt_mean_molecular_weight
  public :: snrt_partition_absorption
  public :: snrt_inventory_tolerance
  public :: snrt_thermochemistry_advance_cell
  public :: snrt_thermochemistry_error_name

contains

  subroutine snrt_secondary_tables_load_from_environment(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: length, env_status

    filename = ''
    call get_environment_variable('SNRT_SECONDARY_TABLE_CONTRACT', filename, &
         length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0 .or. length > len(filename)) then
       call snrt_secondary_tables_reset()
       ierr = snrt_thermochemistry_err_missing
       snrt_thermochemistry_error_message = &
            'SNRT_SECONDARY_TABLE_CONTRACT is not set'
       return
    end if
    call snrt_secondary_tables_load(filename(1:length), ierr)
  end subroutine snrt_secondary_tables_load_from_environment

  subroutine snrt_secondary_tables_load(contract_filename, ierr)
    character(len=*), intent(in) :: contract_filename
    integer, intent(out) :: ierr
    integer :: unit, open_ierr, read_ierr
    integer :: contract_version, nenergy_contract, nxi_contract
    character(len=1024) :: table_directory
    character(len=128) :: source_id, upstream_commit
    character(len=128) :: manifest_sha256, license_name
    character(len=256) :: read_message
    character(len=1024) :: effective_directory

    namelist /snrt_secondary_table_contract/ contract_version, &
         table_directory, source_id, upstream_commit, manifest_sha256, &
         license_name, nenergy_contract, nxi_contract

    call snrt_secondary_tables_reset()
    ierr = snrt_thermochemistry_ok
    if (len_trim(contract_filename) == 0) then
       ierr = snrt_thermochemistry_err_missing
       snrt_thermochemistry_error_message = 'secondary contract path is empty'
       return
    end if

    contract_version = 0
    table_directory = ''
    source_id = ''
    upstream_commit = ''
    manifest_sha256 = ''
    license_name = ''
    nenergy_contract = 0
    nxi_contract = 0

    open(newunit=unit, file=trim(contract_filename), status='old', &
         action='read', form='formatted', iostat=open_ierr)
    if (open_ierr /= 0) then
       ierr = snrt_thermochemistry_err_open
       snrt_thermochemistry_error_message = &
            'secondary table contract could not be opened'
       return
    end if
    read_message = ''
    read(unit, nml=snrt_secondary_table_contract, iostat=read_ierr, &
         iomsg=read_message)
    close(unit)
    if (read_ierr /= 0) then
       ierr = snrt_thermochemistry_err_read
       snrt_thermochemistry_error_message = trim(read_message)
       return
    end if
    if (contract_version /= 1 .or. nenergy_contract /= snrt_secondary_nenergy .or. &
         nxi_contract /= snrt_secondary_nxi) then
       ierr = snrt_thermochemistry_err_identity
       snrt_thermochemistry_error_message = &
            'secondary table contract version or dimensions are unsupported'
       return
    end if
    if (trim(source_id) /= snrt_secondary_source_id .or. &
         trim(upstream_commit) /= snrt_secondary_upstream_commit .or. &
         trim(manifest_sha256) /= snrt_secondary_manifest_sha256 .or. &
         trim(license_name) /= 'MIT' .or. len_trim(table_directory) == 0) then
       ierr = snrt_thermochemistry_err_identity
       snrt_thermochemistry_error_message = &
            'secondary table source identity, license, or directory is invalid'
       return
    end if

    ! The directory is part of the provenance-bound contract.  Do not
    ! permit a process environment variable to substitute a different table
    ! directory after the contract has been validated; the shell runner is
    ! responsible for the byte-level manifest gate before a production launch.
    effective_directory = trim(table_directory)
    call snrt_secondary_tables_load_payload(effective_directory, ierr)
    if (ierr /= snrt_thermochemistry_ok) return

    snrt_secondary_loaded_source_id = trim(source_id)
    snrt_secondary_loaded_upstream_commit = trim(upstream_commit)
    snrt_secondary_loaded_manifest_sha256 = trim(manifest_sha256)
  end subroutine snrt_secondary_tables_load

  subroutine snrt_secondary_tables_load_payload(directory, ierr)
    character(len=*), intent(in) :: directory
    integer, intent(out) :: ierr
    integer :: ix, ie, header, unit, open_ierr, read_ierr
    character(len=1024) :: filename
    character(len=512) :: line
    real(dp) :: energy, fion, fheat, fexc, flya, nhi, nhei, nheii, fshull
    real(dp) :: raw_channel_sum

    ierr = snrt_thermochemistry_ok
    if (len_trim(directory) == 0) then
       ierr = snrt_thermochemistry_err_table_open
       snrt_thermochemistry_error_message = 'secondary table directory is empty'
       return
    end if

    do ix = 1, snrt_secondary_nxi
       filename = trim(directory)//'/'//trim(snrt_table_filename(ix))
       open(newunit=unit, file=trim(filename), status='old', action='read', &
            form='formatted', iostat=open_ierr)
       if (open_ierr /= 0) then
          ierr = snrt_thermochemistry_err_table_open
          snrt_thermochemistry_error_message = &
               'one or more FS2010 table files could not be opened'
          call snrt_secondary_tables_reset_payload()
          return
       end if
       do header = 1, 3
          read(unit, '(A)', iostat=read_ierr) line
          if (read_ierr /= 0) exit
       end do
       if (read_ierr /= 0) then
          close(unit)
          ierr = snrt_thermochemistry_err_table_read
          snrt_thermochemistry_error_message = &
               'FS2010 table header is truncated'
          call snrt_secondary_tables_reset_payload()
          return
       end if
       do ie = 1, snrt_secondary_nenergy
          read(unit, *, iostat=read_ierr) energy, fion, fheat, fexc, flya, &
               nhi, nhei, nheii, fshull
          if (read_ierr /= 0) exit
          if (.not. ieee_is_finite(energy) .or. .not. ieee_is_finite(fion) .or. &
               .not. ieee_is_finite(fheat) .or. .not. ieee_is_finite(fexc) .or. &
               .not. ieee_is_finite(flya) .or. .not. ieee_is_finite(nhi) .or. &
               .not. ieee_is_finite(nhei) .or. .not. ieee_is_finite(nheii) .or. &
               .not. ieee_is_finite(fshull) .or. energy <= 0.0d0 .or. &
               fion < 0.0d0 .or. fheat < 0.0d0 .or. fexc < 0.0d0 .or. &
               flya < 0.0d0 .or. nhi < 0.0d0 .or. nhei < 0.0d0 .or. &
               nheii < 0.0d0 .or. fshull < 0.0d0) then
             read_ierr = 1
             exit
          end if
          if (ix == 1) then
             snrt_secondary_energy_ev(ie) = energy
          else if (abs(energy - snrt_secondary_energy_ev(ie)) > &
               1.0d-10 * max(1.0d0, abs(energy))) then
             read_ierr = 2
             exit
          end if
          snrt_secondary_total_ionization(ix,ie) = fion
          snrt_secondary_heating(ix,ie) = fheat
          snrt_secondary_excitation(ix,ie) = fexc
          snrt_secondary_hydrogen_count(ix,ie) = nhi
          snrt_secondary_helium_i_count(ix,ie) = nhei
          snrt_secondary_helium_ii_count(ix,ie) = nheii
       end do
       close(unit)
       if (read_ierr /= 0) then
          ierr = snrt_thermochemistry_err_table_read
          snrt_thermochemistry_error_message = &
               'FS2010 table row, value, or shared-energy-grid validation failed'
          call snrt_secondary_tables_reset_payload()
          return
       end if
    end do

    if (abs(snrt_secondary_energy_ev(1) - snrt_secondary_energy_minimum_ev) > &
         1.0d-12 .or. abs(snrt_secondary_energy_ev(snrt_secondary_nenergy) - &
         snrt_secondary_energy_maximum_ev) > 1.0d-8 .or. &
         any(snrt_secondary_energy_ev(2:) <= snrt_secondary_energy_ev(:snrt_secondary_nenergy-1))) then
       ierr = snrt_thermochemistry_err_grid
       snrt_thermochemistry_error_message = &
            'FS2010 energy grid is not the pinned increasing 10--9937.21 eV grid'
       call snrt_secondary_tables_reset_payload()
       return
    end if
    do ix = 1, snrt_secondary_nxi
       do ie = 1, snrt_secondary_nenergy
          raw_channel_sum = snrt_secondary_total_ionization(ix,ie) + &
               snrt_secondary_heating(ix,ie) + snrt_secondary_excitation(ix,ie)
          if (abs(raw_channel_sum-1.0d0) > 1.0d-4) then
             ierr = snrt_thermochemistry_err_values
             snrt_thermochemistry_error_message = &
                  'FS2010 raw ionization, heating, and excitation channels do not close'
             call snrt_secondary_tables_reset_payload()
             return
          end if
       end do
    end do
    snrt_secondary_tables_loaded = .true.
  end subroutine snrt_secondary_tables_load_payload

  subroutine snrt_secondary_tables_reset()
    call snrt_secondary_tables_reset_payload()
    snrt_secondary_loaded_source_id = ''
    snrt_secondary_loaded_upstream_commit = ''
    snrt_secondary_loaded_manifest_sha256 = ''
    snrt_thermochemistry_error_message = ''
  end subroutine snrt_secondary_tables_reset

  subroutine snrt_secondary_tables_reset_payload()
    snrt_secondary_energy_ev = 0.0d0
    snrt_secondary_total_ionization = 0.0d0
    snrt_secondary_heating = 0.0d0
    snrt_secondary_excitation = 0.0d0
    snrt_secondary_hydrogen_count = 0.0d0
    snrt_secondary_helium_i_count = 0.0d0
    snrt_secondary_helium_ii_count = 0.0d0
    snrt_secondary_tables_loaded = .false.
  end subroutine snrt_secondary_tables_reset_payload

  subroutine snrt_secondary_fractions(electron_energy_ev, hydrogen_ionized_fraction, &
       heating_fraction, hydrogen_ionization_fraction, helium_i_ionization_fraction, &
       helium_ii_ionization_fraction, excitation_fraction, ierr, &
       total_ionization_fraction)
    real(dp), intent(in) :: electron_energy_ev, hydrogen_ionized_fraction
    real(dp), intent(out) :: heating_fraction, hydrogen_ionization_fraction
    real(dp), intent(out) :: helium_i_ionization_fraction
    real(dp), intent(out) :: helium_ii_ionization_fraction, excitation_fraction
    integer, intent(out) :: ierr
    real(dp), intent(out), optional :: total_ionization_fraction
    real(dp) :: energy, xh, fion, fheat, fexc, nhi, nhei, nheii
    real(dp) :: ionization_weight, total_weight
    real(dp) :: wh, whei, wheii

    heating_fraction = 0.0d0
    hydrogen_ionization_fraction = 0.0d0
    helium_i_ionization_fraction = 0.0d0
    helium_ii_ionization_fraction = 0.0d0
    excitation_fraction = 0.0d0
    if (present(total_ionization_fraction)) total_ionization_fraction = 0.0d0
    ierr = snrt_thermochemistry_ok
    if (.not. snrt_secondary_tables_loaded) then
       ierr = snrt_thermochemistry_err_not_loaded
       return
    end if
    if (.not. ieee_is_finite(electron_energy_ev) .or. &
         .not. ieee_is_finite(hydrogen_ionized_fraction) .or. &
         electron_energy_ev < 0.0d0 .or. hydrogen_ionized_fraction < 0.0d0 .or. &
         hydrogen_ionized_fraction > 1.0d0) then
       ierr = snrt_thermochemistry_err_input
       return
    end if
    if (electron_energy_ev < snrt_secondary_energy_minimum_ev) then
       heating_fraction = 1.0d0
       return
    end if

    energy = min(electron_energy_ev, snrt_secondary_energy_maximum_ev)
    xh = min(max(hydrogen_ionized_fraction, snrt_xi_grid(1)), &
         snrt_xi_grid(snrt_secondary_nxi))
    fion = max(0.0d0, snrt_bilinear(snrt_secondary_total_ionization, energy, xh))
    fheat = max(0.0d0, snrt_bilinear(snrt_secondary_heating, energy, xh))
    fexc = max(0.0d0, snrt_bilinear(snrt_secondary_excitation, energy, xh))
    nhi = max(0.0d0, snrt_bilinear(snrt_secondary_hydrogen_count, energy, xh))
    nhei = max(0.0d0, snrt_bilinear(snrt_secondary_helium_i_count, energy, xh))
    nheii = max(0.0d0, snrt_bilinear(snrt_secondary_helium_ii_count, energy, xh))
    wh = nhi * snrt_hydrogen_i_threshold_ev
    whei = nhei * snrt_helium_i_threshold_ev
    wheii = nheii * snrt_helium_ii_threshold_ev
    ionization_weight = wh + whei + wheii
    if (ionization_weight > tiny(ionization_weight)) then
       hydrogen_ionization_fraction = fion * wh / ionization_weight
       helium_i_ionization_fraction = fion * whei / ionization_weight
       helium_ii_ionization_fraction = fion * wheii / ionization_weight
    end if
    total_weight = fheat + hydrogen_ionization_fraction + &
         helium_i_ionization_fraction + helium_ii_ionization_fraction + fexc
    if (.not. ieee_is_finite(total_weight) .or. total_weight <= 0.0d0) then
       ierr = snrt_thermochemistry_err_values
       return
    end if
    heating_fraction = fheat / total_weight
    hydrogen_ionization_fraction = hydrogen_ionization_fraction / total_weight
    helium_i_ionization_fraction = helium_i_ionization_fraction / total_weight
    helium_ii_ionization_fraction = helium_ii_ionization_fraction / total_weight
    excitation_fraction = fexc / total_weight
    if (present(total_ionization_fraction)) total_ionization_fraction = fion
  end subroutine snrt_secondary_fractions

  real(dp) function snrt_bilinear(values, energy, xh) result(value)
    real(dp), intent(in) :: values(:,:), energy, xh
    integer :: ie, ix
    real(dp) :: we, wx
    real(dp) :: lower_x, upper_x, lower_e, upper_e

    if (energy <= snrt_secondary_energy_ev(1)) then
       ie = 1
       we = 0.0d0
    else if (energy >= snrt_secondary_energy_ev(snrt_secondary_nenergy)) then
       ie = snrt_secondary_nenergy - 1
       we = 1.0d0
    else
       ie = 1
       do while (ie < snrt_secondary_nenergy - 1 .and. &
            energy > snrt_secondary_energy_ev(ie+1))
          ie = ie + 1
       end do
       we = (energy - snrt_secondary_energy_ev(ie)) / &
            (snrt_secondary_energy_ev(ie+1) - snrt_secondary_energy_ev(ie))
    end if
    if (xh <= snrt_xi_grid(1)) then
       ix = 1
       wx = 0.0d0
    else if (xh >= snrt_xi_grid(snrt_secondary_nxi)) then
       ix = snrt_secondary_nxi - 1
       wx = 1.0d0
    else
       ix = 1
       do while (ix < snrt_secondary_nxi - 1 .and. xh > snrt_xi_grid(ix+1))
          ix = ix + 1
       end do
       wx = (xh - snrt_xi_grid(ix)) / (snrt_xi_grid(ix+1) - snrt_xi_grid(ix))
    end if
    lower_e = (1.0d0-we)*values(ix,ie) + we*values(ix,ie+1)
    upper_e = (1.0d0-we)*values(ix+1,ie) + we*values(ix+1,ie+1)
    value = (1.0d0-wx)*lower_e + wx*upper_e
  end function snrt_bilinear

  real(dp) function snrt_alpha_hydrogen_case_b(temperature_k) result(alpha)
    real(dp), intent(in) :: temperature_k
    real(dp) :: temperature, lambda_h

    temperature = snrt_safe_temperature(temperature_k)
    lambda_h = 315614.0d0 / temperature
    alpha = 2.753d-14 * lambda_h**1.5d0 / &
         (1.0d0 + (lambda_h/2.740d0)**0.407d0)**2.242d0
  end function snrt_alpha_hydrogen_case_b

  real(dp) function snrt_alpha_helium_ii_case_b(temperature_k) result(alpha)
    real(dp), intent(in) :: temperature_k
    alpha = snrt_alpha_helium_ii_radiative_case_b(temperature_k) + &
         snrt_alpha_helium_ii_dielectronic_case_b(temperature_k)
  end function snrt_alpha_helium_ii_case_b

  real(dp) function snrt_alpha_helium_ii_radiative_case_b(temperature_k) result(alpha)
    real(dp), intent(in) :: temperature_k
    real(dp) :: temperature, lambda_he

    temperature = snrt_safe_temperature(temperature_k)
    lambda_he = 2.0d0 * 285335.0d0 / temperature
    alpha = 1.26d-14 * lambda_he**0.75d0
  end function snrt_alpha_helium_ii_radiative_case_b

  real(dp) function snrt_alpha_helium_ii_dielectronic_case_b(temperature_k) result(alpha)
    real(dp), intent(in) :: temperature_k
    real(dp) :: temperature

    temperature = snrt_safe_temperature(temperature_k)
    alpha = 1.9d-3 / temperature**1.5d0 * exp(-4.7d5/temperature) * &
         (1.0d0 + 0.3d0*exp(-9.4d4/temperature))
  end function snrt_alpha_helium_ii_dielectronic_case_b

  real(dp) function snrt_alpha_helium_iii_case_b(temperature_k) result(alpha)
    real(dp), intent(in) :: temperature_k

    alpha = 2.0d0 * snrt_alpha_hydrogen_case_b(snrt_safe_temperature(temperature_k)/4.0d0)
  end function snrt_alpha_helium_iii_case_b

  real(dp) function snrt_mean_molecular_weight(x_hydrogen_ii, x_helium_ii, &
       x_helium_iii) result(mu)
    real(dp), intent(in) :: x_hydrogen_ii, x_helium_ii, x_helium_iii
    real(dp) :: h, heii, heiii

    h = min(max(x_hydrogen_ii,0.0d0),1.0d0)
    heii = min(max(x_helium_ii,0.0d0),1.0d0)
    heiii = min(max(x_helium_iii,0.0d0),max(0.0d0,1.0d0-heii))
    mu = 1.0d0 / (snrt_mass_fraction_hydrogen*(1.0d0+h) + &
         0.25d0*snrt_mass_fraction_helium*(1.0d0+heii+2.0d0*heiii))
  end function snrt_mean_molecular_weight

  subroutine snrt_partition_absorption(total_absorbed_code, opacity_species, &
       available_species_code, absorbed_species_code, ierr, &
       unassigned_absorption_code, inventory_scale_code)
    real(dp), intent(in) :: total_absorbed_code
    real(dp), intent(in) :: opacity_species(3)
    real(dp), intent(inout) :: available_species_code(3)
    real(dp), intent(out) :: absorbed_species_code(3)
    integer, intent(out) :: ierr
    real(dp), intent(out), optional :: unassigned_absorption_code
    ! Optional pre-partition scale keeps the tolerance tied to the cell's
    ! original inventory even after earlier spectral groups consume the local
    ! working copy.  Existing scalar callers may omit it.
    real(dp), intent(in), optional :: inventory_scale_code
    real(dp) :: requested(3), headroom(3), opacity_sum, remaining, &
         headroom_sum, addition, scale, tolerance, inventory_total, &
         eligible_inventory, target_absorbed
    integer :: species, pass

    ierr = snrt_thermochemistry_ok
    absorbed_species_code = 0.0d0
    if (present(unassigned_absorption_code)) unassigned_absorption_code = 0.0d0
    if (.not. ieee_is_finite(total_absorbed_code) .or. &
         total_absorbed_code < 0.0d0 .or. any(.not. ieee_is_finite(opacity_species)) .or. &
         any(.not. ieee_is_finite(available_species_code)) .or. &
         any(available_species_code < 0.0d0)) then
       ierr = snrt_thermochemistry_err_input
       return
    end if
    opacity_sum = sum(max(opacity_species,0.0d0))
    if (total_absorbed_code == 0.0d0) return
    if (opacity_sum <= 0.0d0) then
       ierr = snrt_thermochemistry_err_inventory
       return
    end if
    inventory_total = sum(max(available_species_code,0.0d0))
    eligible_inventory = 0.0d0
    do species = 1, 3
       if (opacity_species(species) > 0.0d0) &
            eligible_inventory = eligible_inventory + available_species_code(species)
    end do
    scale = max(total_absorbed_code, inventory_total, eligible_inventory)
    if (present(inventory_scale_code)) scale = max(scale, abs(inventory_scale_code))
    tolerance = snrt_inventory_tolerance(scale)
    target_absorbed = total_absorbed_code
    if (total_absorbed_code > inventory_total) then
       if (total_absorbed_code-inventory_total > tolerance) then
          ierr = snrt_thermochemistry_err_inventory
          return
       end if
       target_absorbed = inventory_total
       if (present(unassigned_absorption_code)) &
            unassigned_absorption_code = total_absorbed_code-inventory_total
    end if
    if (target_absorbed > eligible_inventory) then
       if (target_absorbed-eligible_inventory > tolerance) then
          ierr = snrt_thermochemistry_err_inventory
          return
       end if
       if (present(unassigned_absorption_code)) &
            unassigned_absorption_code = unassigned_absorption_code + &
            target_absorbed-eligible_inventory
       target_absorbed = eligible_inventory
    end if
    if (target_absorbed <= 0.0d0) return
    if (target_absorbed > inventory_total + tolerance) then
       ierr = snrt_thermochemistry_err_inventory
       return
    end if
    requested = target_absorbed * max(opacity_species,0.0d0) / opacity_sum
    do species = 1, 3
       absorbed_species_code(species) = min(requested(species), &
            available_species_code(species))
    end do
    remaining = target_absorbed - sum(absorbed_species_code)
    do pass = 1, 4
       if (remaining <= tolerance) exit
       headroom = 0.0d0
       do species = 1, 3
          if (opacity_species(species) > 0.0d0) headroom(species) = &
               max(available_species_code(species)-absorbed_species_code(species),0.0d0)
       end do
       headroom_sum = sum(headroom)
       if (headroom_sum <= 0.0d0) then
          ierr = snrt_thermochemistry_err_inventory
          return
       end if
       do species = 1, 3
          addition = min(headroom(species), remaining*headroom(species)/headroom_sum)
          absorbed_species_code(species) = absorbed_species_code(species) + addition
       end do
       remaining = target_absorbed - sum(absorbed_species_code)
    end do
    if (remaining > 0.0d0 .and. remaining <= tolerance) then
       do species = 1, 3
          if (opacity_species(species) <= 0.0d0) cycle
          addition = min(max(available_species_code(species)- &
               absorbed_species_code(species),0.0d0), remaining)
          absorbed_species_code(species) = absorbed_species_code(species) + addition
          remaining = target_absorbed - sum(absorbed_species_code)
          if (remaining <= 0.0d0) exit
       end do
    end if
    if (abs(remaining) > tolerance) then
       ierr = snrt_thermochemistry_err_inventory
       return
    end if
    available_species_code = max(0.0d0, available_species_code - &
         absorbed_species_code)
  end subroutine snrt_partition_absorption

  subroutine snrt_thermochemistry_advance_cell(n_hydrogen_cm3, n_helium_cm3, &
       n_h_unit_cm3, temperature_k, delta_t_s, x_hydrogen_ii, x_helium_ii, &
       x_helium_iii, absorbed_species_code, excess_energy_ev, result)
    real(dp), intent(in) :: n_hydrogen_cm3, n_helium_cm3, n_h_unit_cm3
    real(dp), intent(in) :: temperature_k, delta_t_s
    real(dp), intent(in) :: x_hydrogen_ii, x_helium_ii, x_helium_iii
    real(dp), intent(in) :: absorbed_species_code(:,:), excess_energy_ev(:,:)
    type(snrt_thermochemistry_result), intent(out) :: result
    real(dp) :: primary(3), remaining(3), secondary(3)
    real(dp) :: xh_photo, xheii_photo, xheiii_photo
    real(dp) :: heating_energy, ionization_energy, excitation_energy
    real(dp) :: input_energy, residual, electron_energy
    real(dp) :: fh, fhi, fhei, fheii, fexc
    real(dp) :: threshold(3), target_available(3), used_count
    real(dp) :: state_tolerance, inventory_tolerance
    integer :: species, group, ierr_local

    result = snrt_thermochemistry_result()
    result%x_hydrogen_ii = x_hydrogen_ii
    result%x_helium_ii = x_helium_ii
    result%x_helium_iii = x_helium_iii
    if (.not. snrt_secondary_tables_loaded) then
       result%ierr = snrt_thermochemistry_err_not_loaded
       return
    end if
    if (size(absorbed_species_code,1) /= 3 .or. &
         size(excess_energy_ev,1) /= 3 .or. &
         size(absorbed_species_code,2) /= size(excess_energy_ev,2) .or. &
         size(absorbed_species_code,2) /= snrt_thermochemistry_ngroups) then
       result%ierr = snrt_thermochemistry_err_input
       return
    end if
    state_tolerance = snrt_inventory_tolerance(1.0d0)
    if (.not. ieee_is_finite(n_hydrogen_cm3) .or. &
         .not. ieee_is_finite(n_helium_cm3) .or. &
         .not. ieee_is_finite(n_h_unit_cm3) .or. &
         .not. ieee_is_finite(temperature_k) .or. &
         .not. ieee_is_finite(delta_t_s) .or. n_hydrogen_cm3 <= 0.0d0 .or. &
         n_helium_cm3 < 0.0d0 .or. n_h_unit_cm3 <= 0.0d0 .or. &
         temperature_k <= 0.0d0 .or. delta_t_s < 0.0d0 .or. &
         .not. ieee_is_finite(x_hydrogen_ii) .or. &
         .not. ieee_is_finite(x_helium_ii) .or. &
         .not. ieee_is_finite(x_helium_iii) .or. x_hydrogen_ii < 0.0d0 .or. &
         x_hydrogen_ii > 1.0d0 .or. x_helium_ii < 0.0d0 .or. &
         x_helium_iii < 0.0d0 .or. x_helium_ii + x_helium_iii > 1.0d0 + &
         state_tolerance .or. any(.not. ieee_is_finite(absorbed_species_code)) .or. &
         any(absorbed_species_code < 0.0d0) .or. any(.not. ieee_is_finite(excess_energy_ev)) .or. &
         any(excess_energy_ev < 0.0d0)) then
       result%ierr = snrt_thermochemistry_err_input
       return
    end if

    threshold = (/ snrt_hydrogen_i_threshold_ev, snrt_helium_i_threshold_ev, &
         snrt_helium_ii_threshold_ev /)
    primary = 0.0d0
    do species = 1, 3
       primary(species) = sum(absorbed_species_code(species,:)) * n_h_unit_cm3
    end do
    remaining(1) = n_hydrogen_cm3 * (1.0d0-x_hydrogen_ii)
    remaining(2) = n_helium_cm3 * (1.0d0-x_helium_ii-x_helium_iii)
    remaining(3) = n_helium_cm3 * x_helium_ii
    inventory_tolerance = snrt_inventory_tolerance(max(maxval(abs(remaining)), &
         maxval(abs(primary))))
    if (any(primary > remaining + inventory_tolerance)) then
       result%ierr = snrt_thermochemistry_err_inventory
       return
    end if
    remaining = max(0.0d0, remaining-primary)
    result%primary_hydrogen_ionizations_cm3 = primary(1)
    result%primary_helium_i_ionizations_cm3 = primary(2)
    result%primary_helium_ii_ionizations_cm3 = primary(3)

    input_energy = 0.0d0
    heating_energy = 0.0d0
    ionization_energy = 0.0d0
    excitation_energy = 0.0d0
    secondary = 0.0d0
    do group = 1, size(absorbed_species_code,2)
       do species = 1, 3
          electron_energy = absorbed_species_code(species,group) * n_h_unit_cm3 * &
               excess_energy_ev(species,group)
          if (electron_energy <= 0.0d0) cycle
          input_energy = input_energy + electron_energy
          call snrt_secondary_fractions(excess_energy_ev(species,group), &
               x_hydrogen_ii, fh, fhi, fhei, fheii, fexc, ierr_local)
          if (ierr_local /= snrt_thermochemistry_ok) then
             result%ierr = ierr_local
             return
          end if
          heating_energy = heating_energy + electron_energy*fh
          excitation_energy = excitation_energy + electron_energy*fexc
          target_available = remaining
          target_available(1) = max(0.0d0, target_available(1))
          target_available(2) = max(0.0d0, target_available(2))
          target_available(3) = max(0.0d0, target_available(3))
          used_count = min(electron_energy*fhi/threshold(1), target_available(1))
          secondary(1) = secondary(1) + used_count
          remaining(1) = remaining(1) - used_count
          heating_energy = heating_energy + electron_energy*fhi - &
               used_count*threshold(1)
          used_count = min(electron_energy*fhei/threshold(2), target_available(2))
          secondary(2) = secondary(2) + used_count
          remaining(2) = remaining(2) - used_count
          heating_energy = heating_energy + electron_energy*fhei - &
               used_count*threshold(2)
          used_count = min(electron_energy*fheii/threshold(3), target_available(3))
          secondary(3) = secondary(3) + used_count
          remaining(3) = remaining(3) - used_count
          heating_energy = heating_energy + electron_energy*fheii - &
               used_count*threshold(3)
        end do
    end do
    ionization_energy = secondary(1)*threshold(1) + secondary(2)*threshold(2) + &
         secondary(3)*threshold(3)
    residual = input_energy - heating_energy - ionization_energy - excitation_energy
    result%photoelectron_energy_ev_cm3 = input_energy
    result%secondary_heating_energy_ev_cm3 = heating_energy
    result%secondary_ionization_energy_ev_cm3 = ionization_energy
    result%excitation_energy_ev_cm3 = excitation_energy
    result%photoelectron_energy_residual_ev_cm3 = residual
    if (abs(residual) > snrt_energy_tolerance*max(1.0d0, abs(input_energy))) then
       result%ierr = snrt_thermochemistry_err_energy
       return
    end if
    result%secondary_hydrogen_ionizations_cm3 = secondary(1)
    result%secondary_helium_i_ionizations_cm3 = secondary(2)
    result%secondary_helium_ii_ionizations_cm3 = secondary(3)
    result%absorbed_photon_energy_ev_cm3 = &
         primary(1)*threshold(1) + primary(2)*threshold(2) + &
         primary(3)*threshold(3) + input_energy
    result%primary_ionization_energy_ev_cm3 = &
         primary(1)*threshold(1) + primary(2)*threshold(2) + primary(3)*threshold(3)

    xh_photo = x_hydrogen_ii + (primary(1)+secondary(1))/n_hydrogen_cm3
    if (n_helium_cm3 > 0.0d0) then
       xheii_photo = x_helium_ii + (primary(2)+secondary(2)- &
            primary(3)-secondary(3))/n_helium_cm3
       xheiii_photo = x_helium_iii + (primary(3)+secondary(3))/n_helium_cm3
    else
       xheii_photo = 0.0d0
       xheiii_photo = 0.0d0
    end if
    if (xh_photo > 1.0d0 + state_tolerance .or. &
         xheii_photo < -state_tolerance .or. &
         xheiii_photo < -state_tolerance .or. &
         xheii_photo + xheiii_photo > 1.0d0 + state_tolerance) then
       result%ierr = snrt_thermochemistry_err_inventory
       return
    end if
    xh_photo = min(max(xh_photo,0.0d0),1.0d0)
    xheii_photo = min(max(xheii_photo,0.0d0),1.0d0)
    xheiii_photo = min(max(xheiii_photo,0.0d0),max(0.0d0,1.0d0-xheii_photo))
    call snrt_apply_case_b_recombination(n_hydrogen_cm3, n_helium_cm3, &
         temperature_k, delta_t_s, xh_photo, xheii_photo, xheiii_photo, result)
    if (result%ierr /= snrt_thermochemistry_ok) return
    if (delta_t_s > 0.0d0) result%heating_rate_erg_cm3_s = &
         heating_energy*snrt_ev_to_erg/delta_t_s
  end subroutine snrt_thermochemistry_advance_cell

  pure real(dp) function snrt_inventory_tolerance(scale) result(tolerance)
    real(dp), intent(in) :: scale
    real(dp) :: magnitude

    magnitude = max(abs(scale), tiny(1.0d0))
    tolerance = max(snrt_inventory_fp32_ulps*1.1920928955078125d-7*magnitude, &
         snrt_inventory_host_ulps*epsilon(1.0d0)*magnitude)
  end function snrt_inventory_tolerance

  subroutine snrt_apply_case_b_recombination(n_hydrogen_cm3, n_helium_cm3, &
       temperature_k, delta_t_s, x_hydrogen_ii_initial, x_helium_ii_initial, &
       x_helium_iii_initial, result)
    real(dp), intent(in) :: n_hydrogen_cm3, n_helium_cm3, temperature_k, delta_t_s
    real(dp), intent(in) :: x_hydrogen_ii_initial, x_helium_ii_initial
    real(dp), intent(in) :: x_helium_iii_initial
    type(snrt_thermochemistry_result), intent(inout) :: result
    real(dp) :: alpha_h, alpha_heii, alpha_heiii
    real(dp) :: lower, upper, midpoint, residual, implied, maximum_electron
    real(dp) :: xh, xheii, xheiii, electron_density
    real(dp) :: denominator_h, denominator_heii, denominator_heiii
    integer :: iteration

    alpha_h = snrt_alpha_hydrogen_case_b(temperature_k)
    alpha_heii = snrt_alpha_helium_ii_case_b(temperature_k)
    alpha_heiii = snrt_alpha_helium_iii_case_b(temperature_k)
    maximum_electron = n_hydrogen_cm3 + 2.0d0*n_helium_cm3
    if (delta_t_s <= 0.0d0 .or. maximum_electron <= 0.0d0) then
       result%x_hydrogen_ii = x_hydrogen_ii_initial
       result%x_helium_ii = x_helium_ii_initial
       result%x_helium_iii = x_helium_iii_initial
       result%electron_density_cm3 = n_hydrogen_cm3*x_hydrogen_ii_initial + &
            n_helium_cm3*(x_helium_ii_initial+2.0d0*x_helium_iii_initial)
       return
    end if
    lower = 0.0d0
    upper = maximum_electron
    do iteration = 1, 64
       midpoint = 0.5d0*(lower+upper)
       denominator_h = 1.0d0 + delta_t_s*alpha_h*midpoint
       denominator_heiii = 1.0d0 + delta_t_s*alpha_heiii*midpoint
       denominator_heii = 1.0d0 + delta_t_s*alpha_heii*midpoint
       xh = x_hydrogen_ii_initial/denominator_h
       xheiii = x_helium_iii_initial/denominator_heiii
       xheii = (x_helium_ii_initial + delta_t_s*alpha_heiii*midpoint*xheiii) / &
            denominator_heii
       implied = n_hydrogen_cm3*xh + n_helium_cm3*(xheii+2.0d0*xheiii)
       residual = midpoint-implied
       if (residual > 0.0d0) then
          upper = midpoint
       else
          lower = midpoint
       end if
    end do
    electron_density = 0.5d0*(lower+upper)
    denominator_h = 1.0d0 + delta_t_s*alpha_h*electron_density
    denominator_heiii = 1.0d0 + delta_t_s*alpha_heiii*electron_density
    denominator_heii = 1.0d0 + delta_t_s*alpha_heii*electron_density
    xh = x_hydrogen_ii_initial/denominator_h
    xheiii = x_helium_iii_initial/denominator_heiii
    xheii = (x_helium_ii_initial + delta_t_s*alpha_heiii*electron_density*xheiii) / &
         denominator_heii
    xh = min(max(xh,0.0d0),1.0d0)
    xheii = min(max(xheii,0.0d0),1.0d0)
    xheiii = min(max(xheiii,0.0d0),max(0.0d0,1.0d0-xheii))
    result%x_hydrogen_ii = xh
    result%x_helium_ii = xheii
    result%x_helium_iii = xheiii
    result%electron_density_cm3 = n_hydrogen_cm3*xh + &
         n_helium_cm3*(xheii+2.0d0*xheiii)
    result%recombination_hydrogen_cm3 = n_hydrogen_cm3*delta_t_s*alpha_h* &
         result%electron_density_cm3*xh
    result%recombination_helium_ii_cm3 = n_helium_cm3*delta_t_s*alpha_heii* &
         result%electron_density_cm3*xheii
    result%recombination_helium_iii_cm3 = n_helium_cm3*delta_t_s*alpha_heiii* &
         result%electron_density_cm3*xheiii
  end subroutine snrt_apply_case_b_recombination

  real(dp) function snrt_safe_temperature(temperature_k) result(temperature)
    real(dp), intent(in) :: temperature_k

    if (.not. ieee_is_finite(temperature_k)) then
       temperature = 1.0d0
    else
       temperature = max(1.0d0, temperature_k)
    end if
  end function snrt_safe_temperature

  function snrt_thermochemistry_error_name(ierr) result(name)
    integer, intent(in) :: ierr
    character(len=40) :: name

    select case (ierr)
    case (snrt_thermochemistry_ok)
       name = 'ok'
    case (snrt_thermochemistry_err_missing)
       name = 'missing_secondary_contract'
    case (snrt_thermochemistry_err_open)
       name = 'secondary_contract_open_failed'
    case (snrt_thermochemistry_err_read)
       name = 'secondary_contract_read_failed'
    case (snrt_thermochemistry_err_identity)
       name = 'secondary_source_identity_invalid'
    case (snrt_thermochemistry_err_table_open)
       name = 'secondary_table_open_failed'
    case (snrt_thermochemistry_err_table_read)
       name = 'secondary_table_read_failed'
    case (snrt_thermochemistry_err_grid)
       name = 'secondary_energy_grid_invalid'
    case (snrt_thermochemistry_err_values)
       name = 'secondary_values_invalid'
    case (snrt_thermochemistry_err_input)
       name = 'thermochemistry_input_invalid'
    case (snrt_thermochemistry_err_not_loaded)
       name = 'secondary_tables_not_loaded'
    case (snrt_thermochemistry_err_inventory)
       name = 'species_inventory_exceeded'
    case (snrt_thermochemistry_err_energy)
       name = 'photoelectron_energy_not_closed'
    case default
       name = 'unknown_thermochemistry_error'
    end select
  end function snrt_thermochemistry_error_name

end module snrt_thermochemistry

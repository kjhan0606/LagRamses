! DUST-10: native source-bound opacity and thermal contract admission.
!
! The upstream JSON tooling owns scientific sidecar construction and file
! hashing.  This runtime-facing namelist is an explicit, bounded transport
! representation.  Candidate records are readable but never runtime-enabled.
module snrt_dust_contract
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter, public :: snrt_dust_contract_dp = real64
  integer, parameter, public :: snrt_dust_contract_max_groups = 32
  integer, parameter, public :: snrt_dust_contract_max_temperature = 256
  integer, parameter, public :: snrt_dust_contract_ok = 0
  integer, parameter, public :: snrt_dust_contract_err_missing = 1
  integer, parameter, public :: snrt_dust_contract_err_open = 2
  integer, parameter, public :: snrt_dust_contract_err_read = 3
  integer, parameter, public :: snrt_dust_contract_err_version = 4
  integer, parameter, public :: snrt_dust_contract_err_identity = 5
  integer, parameter, public :: snrt_dust_contract_err_values = 6
  integer, parameter, public :: snrt_dust_contract_err_status = 7
  integer, parameter, public :: snrt_dust_contract_err_not_loaded = 8

  integer, save, public :: snrt_dust_contract_number_groups = 0
  integer, save, public :: snrt_dust_contract_number_temperature = 0
  integer, save, public :: snrt_dust_contract_version = 0
  real(real64), save, public :: snrt_dust_contract_group_edges_ev( &
       snrt_dust_contract_max_groups + 1) = 0.0d0
  real(real64), save, public :: snrt_dust_contract_absorption_per_h_cm2( &
       snrt_dust_contract_max_groups) = 0.0d0
  real(real64), save, public :: snrt_dust_contract_absorption_mean_energy_ev( &
       snrt_dust_contract_max_groups) = 0.0d0
  real(real64), save, public :: snrt_dust_contract_temperature_k( &
       snrt_dust_contract_max_temperature) = 0.0d0
  real(real64), save, public :: snrt_dust_contract_emitted_power_per_h_erg_s( &
       snrt_dust_contract_max_temperature) = 0.0d0
  real(real64), save, public :: snrt_dust_contract_mass_per_h_g = 0.0d0
  real(real64), save, public :: snrt_dust_contract_heat_capacity_per_h_erg_k = 0.0d0
  logical, save, public :: snrt_dust_contract_loaded = .false.
  logical, save, public :: snrt_dust_contract_runtime_allowed = .false.
  character(len=64), save, public :: snrt_dust_contract_opacity_status = ''
  character(len=64), save, public :: snrt_dust_contract_thermal_status = ''
  character(len=128), save, public :: snrt_dust_contract_source_id = ''
  character(len=128), save, public :: snrt_dust_contract_source_sha256 = ''
  character(len=128), save, public :: snrt_dust_contract_source_table_sha256 = ''
  character(len=128), save, public :: snrt_dust_contract_group_edges_sha256 = ''
  character(len=128), save, public :: snrt_dust_contract_approval_id = ''
  character(len=256), save, public :: snrt_dust_contract_thermal_source = ''
  character(len=256), save, public :: snrt_dust_contract_error_message = ''

  public :: snrt_dust_contract_load
  public :: snrt_dust_contract_load_from_environment
  public :: snrt_dust_contract_reset
  public :: snrt_dust_contract_error_name

contains

  subroutine snrt_dust_contract_load(filename, ierr)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr
    integer :: unit, open_ierr, read_ierr
    integer :: contract_version, ngroups_input, ntemperature_input
    character(len=64) :: opacity_status, thermal_status
    character(len=128) :: source_id, source_sha256, source_table_sha256
    character(len=128) :: group_edges_sha256, approval_id
    character(len=256) :: thermal_source, read_message
    real(real64) :: edges_input(snrt_dust_contract_max_groups + 1)
    real(real64) :: absorption_input(snrt_dust_contract_max_groups)
    real(real64) :: mean_energy_input(snrt_dust_contract_max_groups)
    real(real64) :: temperature_input(snrt_dust_contract_max_temperature)
    real(real64) :: power_input(snrt_dust_contract_max_temperature)
    real(real64) :: mass_per_h_input
    real(real64) :: heat_capacity_per_h_erg_k_input

    namelist /snrt_dust_contract/ contract_version, ngroups_input, &
         ntemperature_input, opacity_status, thermal_status, source_id, &
         source_sha256, source_table_sha256, group_edges_sha256, approval_id, &
         thermal_source, edges_input, absorption_input, mean_energy_input, &
         temperature_input, power_input, mass_per_h_input, &
         heat_capacity_per_h_erg_k_input

    call snrt_dust_contract_reset()
    ierr = snrt_dust_contract_ok
    if (len_trim(filename) == 0) then
       ierr = snrt_dust_contract_err_missing
       snrt_dust_contract_error_message = 'dust contract path is empty'
       return
    end if

    contract_version = 0
    ngroups_input = 0
    ntemperature_input = 0
    opacity_status = ''
    thermal_status = ''
    source_id = ''
    source_sha256 = ''
    source_table_sha256 = ''
    group_edges_sha256 = ''
    approval_id = ''
    thermal_source = ''
    edges_input = -1.0d0
    absorption_input = -1.0d0
    mean_energy_input = -1.0d0
    temperature_input = -1.0d0
    power_input = -1.0d0
    mass_per_h_input = -1.0d0
    heat_capacity_per_h_erg_k_input = -1.0d0

    open(newunit=unit, file=trim(filename), status='old', action='read', &
         form='formatted', iostat=open_ierr)
    if (open_ierr /= 0) then
       ierr = snrt_dust_contract_err_open
       snrt_dust_contract_error_message = 'dust contract file could not be opened'
       return
    end if
    read_message = ''
    read(unit, nml=snrt_dust_contract, iostat=read_ierr, iomsg=read_message)
    close(unit)
    if (read_ierr /= 0) then
       ierr = snrt_dust_contract_err_read
       snrt_dust_contract_error_message = trim(read_message)
       return
    end if
    if (contract_version /= 1 .and. contract_version /= 2) then
       ierr = snrt_dust_contract_err_version
       snrt_dust_contract_error_message = &
            'only snrt_dust_contract versions 1 and 2 are supported'
       return
    end if
    if (.not. known_opacity_status(opacity_status) .or. &
         .not. known_thermal_status(thermal_status)) then
       ierr = snrt_dust_contract_err_status
       snrt_dust_contract_error_message = 'dust opacity or thermal status is not recognized'
       return
    end if
    if (len_trim(source_id) == 0 .or. .not. is_sha256(source_sha256) .or. &
         .not. is_sha256(source_table_sha256) .or. &
         .not. is_sha256(group_edges_sha256) .or. len_trim(thermal_source) == 0) then
       ierr = snrt_dust_contract_err_identity
       snrt_dust_contract_error_message = 'dust source identity or hash token is incomplete'
       return
    end if
    if (ngroups_input < 1 .or. ngroups_input > snrt_dust_contract_max_groups .or. &
         ntemperature_input < 2 .or. &
         ntemperature_input > snrt_dust_contract_max_temperature) then
       ierr = snrt_dust_contract_err_values
       snrt_dust_contract_error_message = 'dust contract array bounds are invalid'
       return
    end if
    if (mass_per_h_input <= 0.0d0 .or. &
         .not. ieee_is_finite(mass_per_h_input)) then
       ierr = snrt_dust_contract_err_values
       snrt_dust_contract_error_message = 'reference dust mass per H is invalid'
       return
    end if
    if (contract_version >= 2 .and. &
         (heat_capacity_per_h_erg_k_input <= 0.0d0 .or. &
          .not. ieee_is_finite(heat_capacity_per_h_erg_k_input))) then
       ierr = snrt_dust_contract_err_values
       snrt_dust_contract_error_message = &
            'version 2 requires a positive dust heat capacity per H'
       return
    end if
    if (.not. valid_opacity_values(ngroups_input, edges_input, absorption_input, &
         mean_energy_input) .or. .not. valid_thermal_values(ntemperature_input, &
         temperature_input, power_input)) then
       ierr = snrt_dust_contract_err_values
       snrt_dust_contract_error_message = 'dust opacity or thermal table values are invalid'
       return
    end if

    snrt_dust_contract_version = contract_version
    snrt_dust_contract_number_groups = ngroups_input
    snrt_dust_contract_number_temperature = ntemperature_input
    snrt_dust_contract_group_edges_ev(1:ngroups_input + 1) = &
         edges_input(1:ngroups_input + 1)
    snrt_dust_contract_absorption_per_h_cm2(1:ngroups_input) = &
         absorption_input(1:ngroups_input)
    snrt_dust_contract_absorption_mean_energy_ev(1:ngroups_input) = &
         mean_energy_input(1:ngroups_input)
    snrt_dust_contract_temperature_k(1:ntemperature_input) = &
         temperature_input(1:ntemperature_input)
    snrt_dust_contract_emitted_power_per_h_erg_s(1:ntemperature_input) = &
         power_input(1:ntemperature_input)
    snrt_dust_contract_mass_per_h_g = mass_per_h_input
    if (contract_version >= 2) snrt_dust_contract_heat_capacity_per_h_erg_k = &
         heat_capacity_per_h_erg_k_input
    snrt_dust_contract_opacity_status = trim(opacity_status)
    snrt_dust_contract_thermal_status = trim(thermal_status)
    snrt_dust_contract_source_id = trim(source_id)
    snrt_dust_contract_source_sha256 = trim(source_sha256)
    snrt_dust_contract_source_table_sha256 = trim(source_table_sha256)
    snrt_dust_contract_group_edges_sha256 = trim(group_edges_sha256)
    snrt_dust_contract_approval_id = trim(approval_id)
    snrt_dust_contract_thermal_source = trim(thermal_source)
    snrt_dust_contract_loaded = .true.
    snrt_dust_contract_runtime_allowed = &
         contract_version >= 2 .and. &
         trim(opacity_status) == 'approved_production' .and. &
         trim(thermal_status) == 'approved_thermal_production' .and. &
         len_trim(approval_id) > 0 .and. &
         heat_capacity_per_h_erg_k_input > 0.0d0
  end subroutine snrt_dust_contract_load


  subroutine snrt_dust_contract_load_from_environment(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: length, env_status

    filename = ''
    call get_environment_variable('SNRT_DUST_CONTRACT', filename, &
         length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0 .or. length > len(filename)) then
       call snrt_dust_contract_reset()
       ierr = snrt_dust_contract_err_missing
       snrt_dust_contract_error_message = 'SNRT_DUST_CONTRACT is not set'
       return
    end if
    call snrt_dust_contract_load(filename(1:length), ierr)
  end subroutine snrt_dust_contract_load_from_environment


  subroutine snrt_dust_contract_reset()
    snrt_dust_contract_number_groups = 0
    snrt_dust_contract_number_temperature = 0
    snrt_dust_contract_version = 0
    snrt_dust_contract_group_edges_ev = 0.0d0
    snrt_dust_contract_absorption_per_h_cm2 = 0.0d0
    snrt_dust_contract_absorption_mean_energy_ev = 0.0d0
    snrt_dust_contract_temperature_k = 0.0d0
    snrt_dust_contract_emitted_power_per_h_erg_s = 0.0d0
    snrt_dust_contract_mass_per_h_g = 0.0d0
    snrt_dust_contract_heat_capacity_per_h_erg_k = 0.0d0
    snrt_dust_contract_loaded = .false.
    snrt_dust_contract_runtime_allowed = .false.
    snrt_dust_contract_opacity_status = ''
    snrt_dust_contract_thermal_status = ''
    snrt_dust_contract_source_id = ''
    snrt_dust_contract_source_sha256 = ''
    snrt_dust_contract_source_table_sha256 = ''
    snrt_dust_contract_group_edges_sha256 = ''
    snrt_dust_contract_approval_id = ''
    snrt_dust_contract_thermal_source = ''
    snrt_dust_contract_error_message = ''
  end subroutine snrt_dust_contract_reset


  function snrt_dust_contract_error_name(ierr) result(name)
    integer, intent(in) :: ierr
    character(len=32) :: name

    select case (ierr)
    case (snrt_dust_contract_ok)
       name = 'ok'
    case (snrt_dust_contract_err_missing)
       name = 'missing'
    case (snrt_dust_contract_err_open)
       name = 'open'
    case (snrt_dust_contract_err_read)
       name = 'read'
    case (snrt_dust_contract_err_version)
       name = 'version'
    case (snrt_dust_contract_err_identity)
       name = 'identity'
    case (snrt_dust_contract_err_values)
       name = 'values'
    case (snrt_dust_contract_err_status)
       name = 'status'
    case (snrt_dust_contract_err_not_loaded)
       name = 'not_loaded'
    case default
       name = 'unknown'
    end select
  end function snrt_dust_contract_error_name


  logical function valid_opacity_values(ngroups, edges, absorption, mean_energy)
    integer, intent(in) :: ngroups
    real(real64), intent(in) :: edges(:), absorption(:), mean_energy(:)
    integer :: group

    valid_opacity_values = .false.
    if (size(edges) < ngroups + 1 .or. size(absorption) < ngroups .or. &
         size(mean_energy) < ngroups) return
    if (.not. all(ieee_is_finite(edges(1:ngroups + 1))) .or. &
         .not. all(ieee_is_finite(absorption(1:ngroups))) .or. &
         .not. all(ieee_is_finite(mean_energy(1:ngroups)))) return
    if (any(edges(1:ngroups + 1) <= 0.0d0) .or. &
         any(edges(2:ngroups + 1) <= edges(1:ngroups)) .or. &
         any(absorption(1:ngroups) < 0.0d0) .or. &
         any(mean_energy(1:ngroups) <= 0.0d0)) return
    do group = 1, ngroups
       if (mean_energy(group) < edges(group) .or. &
            mean_energy(group) > edges(group + 1)) return
    end do
    valid_opacity_values = .true.
  end function valid_opacity_values


  logical function valid_thermal_values(ntemperature, temperature, power)
    integer, intent(in) :: ntemperature
    real(real64), intent(in) :: temperature(:), power(:)

    valid_thermal_values = .false.
    if (size(temperature) < ntemperature .or. size(power) < ntemperature) return
    if (.not. all(ieee_is_finite(temperature(1:ntemperature))) .or. &
         .not. all(ieee_is_finite(power(1:ntemperature)))) return
    if (any(temperature(1:ntemperature) <= 0.0d0) .or. &
         any(power(1:ntemperature) <= 0.0d0) .or. &
         any(temperature(2:ntemperature) <= temperature(1:ntemperature - 1)) .or. &
         any(power(2:ntemperature) <= power(1:ntemperature - 1))) return
    valid_thermal_values = .true.
  end function valid_thermal_values


  logical function known_opacity_status(status)
    character(len=*), intent(in) :: status
    known_opacity_status = trim(status) == 'reference_control' .or. &
         trim(status) == 'reference_scattering_control' .or. &
         trim(status) == 'candidate_source_sed_matched' .or. &
         trim(status) == 'candidate_scattering_isotropic' .or. &
         trim(status) == 'approved_production'
  end function known_opacity_status


  logical function known_thermal_status(status)
    character(len=*), intent(in) :: status
    known_thermal_status = trim(status) == 'reference_thermal_control' .or. &
         trim(status) == 'candidate_kirchhoff_equilibrium' .or. &
         trim(status) == 'approved_thermal_production'
  end function known_thermal_status


  logical function is_sha256(value)
    character(len=*), intent(in) :: value
    integer :: i, code

    is_sha256 = .false.
    if (len_trim(value) /= 64) return
    do i = 1, 64
       code = iachar(value(i:i))
       if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
            (code >= iachar('a') .and. code <= iachar('f')) .or. &
            (code >= iachar('A') .and. code <= iachar('F')))) return
    end do
    is_sha256 = .true.
  end function is_sha256

end module snrt_dust_contract

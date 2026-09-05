! Native multigroup spectral contract for the SNRT source/chemistry boundary.
!
! The photon state is indexed by group number, so changing the group table
! without changing the state checkpoint and source closure together would
! silently change the physical meaning of an existing run.  This module owns
! the canonical nine-group dimensions and loads all source-dependent moments
! from an explicit Fortran namelist.  It deliberately does not parse JSON in
! the RAMSES process; the JSON ledger remains the review/provenance sidecar.
!
! A reference-control contract may be used for wiring tests, but it requires
! the explicit SNRT_ALLOW_REFERENCE_CONTROL=1 opt-in at runtime. Candidate
! physical SEDs are readable for inspection but are not runtime-admissible.
! Only an explicitly approved production contract or the explicitly enabled
! reference control can pass the runtime gate.
module snrt_spectral_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: dp
  implicit none

  private

  integer, parameter, public :: snrt_ngroups = 9
  integer, parameter, public :: snrt_nedges = snrt_ngroups + 1

  real(dp), parameter, public :: snrt_group_edges_ev(snrt_nedges) = &
       (/ 0.01d0, 1.0d0, 5.6d0, 11.2d0, 13.6d0, 24.59d0, &
          54.42d0, 500.0d0, 2000.0d0, 10000.0d0 /)
  character(len=*), parameter, public :: snrt_group_interval_convention = &
       'left_closed_right_open_except_final_closed'
  character(len=*), parameter, public :: snrt_group_edges_sha256 = &
       'd28f78f1703730c6c0b9a7d183edfe0c5e6337979e737ce002a572b66fc53ff1'

  integer, parameter, public :: snrt_spectral_contract_ok = 0
  integer, parameter, public :: snrt_spectral_contract_err_missing = 1
  integer, parameter, public :: snrt_spectral_contract_err_open = 2
  integer, parameter, public :: snrt_spectral_contract_err_read = 3
  integer, parameter, public :: snrt_spectral_contract_err_version = 4
  integer, parameter, public :: snrt_spectral_contract_err_identity = 5
  integer, parameter, public :: snrt_spectral_contract_err_status = 6
  integer, parameter, public :: snrt_spectral_contract_err_edges = 7
  integer, parameter, public :: snrt_spectral_contract_err_values = 8
  integer, parameter, public :: snrt_spectral_contract_err_not_loaded = 9
  integer, parameter, public :: snrt_spectral_contract_err_fraction_semantics = 10

  character(len=*), parameter, public :: snrt_spectral_status_reference = &
       'reference_control'
  character(len=*), parameter, public :: snrt_spectral_status_candidate = &
       'candidate_explicit_sed'
  character(len=*), parameter, public :: snrt_spectral_status_production = &
       'approved_production'

  ! These are populated only after a complete contract passes validation.
  ! The H I array retains the historical public name used by the driver.
  real(dp), save, public :: snrt_group_mean_energy_ev(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_energy_fraction(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_cross_section_cm2(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_cross_section_hei_cm2(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_cross_section_heii_cm2(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_photoelectron_excess_energy_ev(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_photoelectron_excess_hei_ev(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_photoelectron_excess_heii_ev(snrt_ngroups) = 0.0d0
  real(dp), save, public :: snrt_group_energy_fraction_sum = 0.0d0
  real(dp), save, public :: snrt_group_unrepresented_energy_fraction = 1.0d0

  logical, save, public :: snrt_spectral_contract_loaded = .false.
  logical, save, public :: snrt_spectral_contract_runtime_allowed = .false.
  character(len=64), save, public :: snrt_spectral_contract_status = ''
  character(len=128), save, public :: snrt_spectral_contract_source_id = ''
  character(len=128), save, public :: snrt_spectral_contract_source_sha256 = ''
  character(len=128), save, public :: snrt_spectral_contract_source_commit_binding = ''
  character(len=128), save, public :: snrt_spectral_contract_approval_id = ''
  character(len=128), save, public :: snrt_spectral_contract_group_edges_sha256 = ''
  character(len=128), save, public :: snrt_spectral_contract_interval_convention = ''
  character(len=64), save, public :: snrt_spectral_contract_fraction_semantics = ''
  character(len=256), save, public :: snrt_spectral_contract_error_message = ''

  public :: snrt_spectral_contract_load
  public :: snrt_spectral_contract_load_from_environment
  public :: snrt_spectral_contract_validate_values
  public :: snrt_spectral_contract_error_name
  public :: snrt_spectral_contract_checkpoint_identity_matches

contains

  subroutine snrt_spectral_contract_load(filename, ierr)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr

    integer :: unit, open_ierr, read_ierr
    character(len=256) :: read_message
    integer :: contract_version
    character(len=64) :: contract_status
    character(len=128) :: source_id, source_sha256, source_commit_binding
    character(len=128) :: approval_id
    character(len=128) :: edges_sha256, interval_convention
    character(len=64) :: fraction_semantics
    real(dp) :: edges_input(snrt_nedges)
    real(dp) :: mean_input(snrt_ngroups), fraction_input(snrt_ngroups)
    real(dp) :: hi_input(snrt_ngroups), hei_input(snrt_ngroups)
    real(dp) :: heii_input(snrt_ngroups)
    real(dp) :: excess_hi_input(snrt_ngroups), excess_hei_input(snrt_ngroups)
    real(dp) :: excess_heii_input(snrt_ngroups)

    namelist /snrt_group_contract/ contract_version, contract_status, source_id, &
         source_sha256, source_commit_binding, approval_id, edges_sha256, &
         interval_convention, fraction_semantics, edges_input, &
         mean_input, fraction_input, hi_input, hei_input, heii_input, &
         excess_hi_input, excess_hei_input, excess_heii_input

    call snrt_spectral_contract_reset()
    ierr = snrt_spectral_contract_ok
    if (len_trim(filename) == 0) then
       ierr = snrt_spectral_contract_err_missing
       snrt_spectral_contract_error_message = 'contract path is empty'
       return
    end if

    contract_version = 0
    contract_status = ''
    source_id = ''
    source_sha256 = ''
    source_commit_binding = ''
    approval_id = ''
    edges_sha256 = ''
    interval_convention = ''
    fraction_semantics = ''
    edges_input = -1.0d0
    mean_input = -1.0d0
    fraction_input = -1.0d0
    hi_input = -1.0d0
    hei_input = -1.0d0
    heii_input = -1.0d0
    excess_hi_input = -1.0d0
    excess_hei_input = -1.0d0
    excess_heii_input = -1.0d0

    open(newunit=unit, file=trim(filename), status='old', action='read', &
         form='formatted', iostat=open_ierr)
    if (open_ierr /= 0) then
       ierr = snrt_spectral_contract_err_open
       snrt_spectral_contract_error_message = 'contract file could not be opened'
       return
    end if
    read_message = ''
    read(unit, nml=snrt_group_contract, iostat=read_ierr, iomsg=read_message)
    close(unit)
    if (read_ierr /= 0) then
       snrt_spectral_contract_error_message = trim(read_message)
       ierr = snrt_spectral_contract_err_read
       return
    end if
    if (contract_version /= 1) then
       ierr = snrt_spectral_contract_err_version
       snrt_spectral_contract_error_message = &
            'only snrt_group_contract version 1 is supported'
       return
    end if
    if (.not. snrt_spectral_status_is_known(contract_status)) then
       ierr = snrt_spectral_contract_err_status
       snrt_spectral_contract_error_message = 'contract status is not recognized'
       return
    end if
    if (.not. valid_identity(source_id, source_sha256, source_commit_binding, &
         approval_id, edges_sha256)) then
       ierr = snrt_spectral_contract_err_identity
       snrt_spectral_contract_error_message = &
            'source/hash/commit/approval identity is incomplete or malformed'
       return
    end if
    if (.not. snrt_spectral_fraction_semantics_is_known(fraction_semantics)) then
       ierr = snrt_spectral_contract_err_fraction_semantics
       snrt_spectral_contract_error_message = &
            'fraction_semantics must be intrinsic or escaped'
       return
    end if
    if (trim(edges_sha256) /= trim(snrt_group_edges_sha256) .or. &
         trim(interval_convention) /= trim(snrt_group_interval_convention)) then
       ierr = snrt_spectral_contract_err_edges
       snrt_spectral_contract_error_message = &
            'edge digest or interval convention does not match canonical ledger'
       return
    end if
    call snrt_spectral_contract_validate_values(edges_input, mean_input, fraction_input, &
         hi_input, hei_input, heii_input, excess_hi_input, excess_hei_input, &
         excess_heii_input, ierr)
    if (ierr /= snrt_spectral_contract_ok) then
       snrt_spectral_contract_error_message = &
            trim(snrt_spectral_contract_error_name(ierr))
       return
    end if

    snrt_group_mean_energy_ev = mean_input
    snrt_group_energy_fraction = fraction_input
    snrt_group_cross_section_cm2 = hi_input
    snrt_group_cross_section_hei_cm2 = hei_input
    snrt_group_cross_section_heii_cm2 = heii_input
    snrt_group_photoelectron_excess_energy_ev = excess_hi_input
    snrt_group_photoelectron_excess_hei_ev = excess_hei_input
    snrt_group_photoelectron_excess_heii_ev = excess_heii_input
    snrt_group_energy_fraction_sum = sum(fraction_input)
    snrt_group_unrepresented_energy_fraction = max(0.0d0, &
         1.0d0 - snrt_group_energy_fraction_sum)
    snrt_spectral_contract_status = trim(contract_status)
    snrt_spectral_contract_source_id = trim(source_id)
    snrt_spectral_contract_source_sha256 = trim(source_sha256)
    snrt_spectral_contract_source_commit_binding = trim(source_commit_binding)
    snrt_spectral_contract_approval_id = trim(approval_id)
    snrt_spectral_contract_group_edges_sha256 = trim(edges_sha256)
    snrt_spectral_contract_interval_convention = trim(interval_convention)
    snrt_spectral_contract_fraction_semantics = trim(fraction_semantics)
    snrt_spectral_contract_loaded = .true.
    snrt_spectral_contract_runtime_allowed = runtime_status_allowed(contract_status) .and. &
         trim(fraction_semantics) == 'escaped'
    if (runtime_status_allowed(contract_status) .and. &
         trim(fraction_semantics) /= 'escaped') then
       snrt_spectral_contract_error_message = &
            'resolved-domain SNRT requires escaped fraction semantics'
    end if
  end subroutine snrt_spectral_contract_load

  subroutine snrt_spectral_contract_load_from_environment(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: length, env_status

    filename = ''
    call get_environment_variable('SNRT_GROUP_CONTRACT', filename, &
         length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) then
       call snrt_spectral_contract_reset()
       ierr = snrt_spectral_contract_err_missing
       snrt_spectral_contract_error_message = &
            'SNRT_GROUP_CONTRACT is not set'
       return
    end if
    call snrt_spectral_contract_load(filename(1:length), ierr)
  end subroutine snrt_spectral_contract_load_from_environment

  subroutine snrt_spectral_contract_reset()
    snrt_group_mean_energy_ev = 0.0d0
    snrt_group_energy_fraction = 0.0d0
    snrt_group_cross_section_cm2 = 0.0d0
    snrt_group_cross_section_hei_cm2 = 0.0d0
    snrt_group_cross_section_heii_cm2 = 0.0d0
    snrt_group_photoelectron_excess_energy_ev = 0.0d0
    snrt_group_photoelectron_excess_hei_ev = 0.0d0
    snrt_group_photoelectron_excess_heii_ev = 0.0d0
    snrt_group_energy_fraction_sum = 0.0d0
    snrt_group_unrepresented_energy_fraction = 1.0d0
    snrt_spectral_contract_loaded = .false.
    snrt_spectral_contract_runtime_allowed = .false.
    snrt_spectral_contract_status = ''
    snrt_spectral_contract_source_id = ''
    snrt_spectral_contract_source_sha256 = ''
    snrt_spectral_contract_source_commit_binding = ''
    snrt_spectral_contract_approval_id = ''
    snrt_spectral_contract_group_edges_sha256 = ''
    snrt_spectral_contract_interval_convention = ''
    snrt_spectral_contract_fraction_semantics = ''
    snrt_spectral_contract_error_message = ''
  end subroutine snrt_spectral_contract_reset

  subroutine snrt_spectral_contract_validate_values(edges, mean_energy, fraction, hi, hei, &
       heii, excess_hi, excess_hei, excess_heii, ierr)
    real(dp), intent(in) :: edges(:), mean_energy(:), fraction(:), hi(:), hei(:)
    real(dp), intent(in) :: heii(:), excess_hi(:), excess_hei(:), excess_heii(:)
    integer, intent(out) :: ierr
    integer :: igroup
    real(dp), parameter :: relative_tolerance = 1.0d-12
    real(dp), parameter :: absolute_tolerance = 1.0d-14

    ierr = snrt_spectral_contract_ok
    if (size(edges) /= snrt_nedges .or. size(mean_energy) /= snrt_ngroups .or. &
         size(fraction) /= snrt_ngroups .or. size(hi) /= snrt_ngroups .or. &
         size(hei) /= snrt_ngroups .or. size(heii) /= snrt_ngroups .or. &
         size(excess_hi) /= snrt_ngroups .or. size(excess_hei) /= snrt_ngroups .or. &
         size(excess_heii) /= snrt_ngroups) then
       ierr = snrt_spectral_contract_err_values
       return
    end if
    if (any(.not. ieee_is_finite(edges)) .or. any(edges <= 0.0d0) .or. &
         any(edges(2:) <= edges(:snrt_ngroups))) then
       ierr = snrt_spectral_contract_err_edges
       return
    end if
    do igroup = 1, snrt_nedges
       if (abs(edges(igroup)-snrt_group_edges_ev(igroup)) > &
            absolute_tolerance + relative_tolerance * &
            max(abs(edges(igroup)),abs(snrt_group_edges_ev(igroup)))) then
          ierr = snrt_spectral_contract_err_edges
          return
       end if
    end do
    if (any(.not. ieee_is_finite(mean_energy)) .or. &
         any(.not. ieee_is_finite(fraction)) .or. &
         any(.not. ieee_is_finite(hi)) .or. any(.not. ieee_is_finite(hei)) .or. &
         any(.not. ieee_is_finite(heii)) .or. &
         any(.not. ieee_is_finite(excess_hi)) .or. &
         any(.not. ieee_is_finite(excess_hei)) .or. &
         any(.not. ieee_is_finite(excess_heii))) then
       ierr = snrt_spectral_contract_err_values
       return
    end if
    if (any(mean_energy < edges(:snrt_ngroups)) .or. &
         any(mean_energy > edges(2:))) then
       ierr = snrt_spectral_contract_err_values
       return
    end if
    if (any(fraction < 0.0d0) .or. sum(fraction) <= 0.0d0 .or. &
         sum(fraction) > 1.0d0 + 1.0d-12) then
       ierr = snrt_spectral_contract_err_values
       return
    end if
    if (any(hi < 0.0d0) .or. any(hei < 0.0d0) .or. any(heii < 0.0d0) .or. &
         any(excess_hi < 0.0d0) .or. any(excess_hei < 0.0d0) .or. &
         any(excess_heii < 0.0d0)) then
       ierr = snrt_spectral_contract_err_values
       return
    end if

    call validate_species_table(edges, hi, excess_hi, 13.6d0, ierr)
    if (ierr /= snrt_spectral_contract_ok) return
    call validate_species_table(edges, hei, excess_hei, 24.59d0, ierr)
    if (ierr /= snrt_spectral_contract_ok) return
    call validate_species_table(edges, heii, excess_heii, 54.42d0, ierr)
  end subroutine snrt_spectral_contract_validate_values

  subroutine validate_species_table(edges, cross_section, excess_energy, &
       threshold, ierr)
    real(dp), intent(in) :: edges(:), cross_section(:), excess_energy(:)
    real(dp), intent(in) :: threshold
    integer, intent(out) :: ierr
    integer :: igroup
    real(dp), parameter :: tolerance = 1.0d-12
    ! The largest threshold-adjacent primordial cross section in the
    ! supported H/He tables is about 7.4e-18 cm^2 (He I).  A 1e-17 ceiling
    ! leaves a conservative factor-of-two margin for group averaging while
    ! rejecting a one-decade decimal slip in a reference-scale value.
    real(dp), parameter :: max_cross_section_cm2 = 1.0d-17
    real(dp) :: max_excess_energy

    ierr = snrt_spectral_contract_ok
    do igroup = 1, snrt_ngroups
       if (edges(igroup+1) <= threshold + tolerance) then
          ! Cross sections are in cm^2 and can be far below 1e-12; any
          ! positive sub-threshold value is an actual opacity bug, not a
          ! floating-point noise that should be tolerated.
          if (cross_section(igroup) > 0.0d0 .or. &
               excess_energy(igroup) > 0.0d0) then
             ierr = snrt_spectral_contract_err_values
             return
          end if
       else if (edges(igroup) >= threshold - tolerance) then
          ! The group-averaged cross section is bounded conservatively above
          ! the largest primordial H/He photoionization cross section.  The
          ! excess energy is an absorber-weighted mean and cannot exceed the
          ! top of this group's interval measured from the threshold.  These
          ! upper bounds catch a decimal/transcription error that positivity
          ! checks alone would admit.
          max_excess_energy = max(0.0d0, edges(igroup+1) - threshold)
          if (cross_section(igroup) > max_cross_section_cm2 .or. &
               excess_energy(igroup) > max_excess_energy + &
               tolerance * max(1.0d0, abs(max_excess_energy))) then
             ierr = snrt_spectral_contract_err_values
             return
          end if
          if (cross_section(igroup) <= 0.0d0 .or. excess_energy(igroup) <= 0.0d0) then
             ierr = snrt_spectral_contract_err_values
             return
          end if
       end if
    end do
  end subroutine validate_species_table

  logical function valid_identity(source_id, source_sha256, source_commit_binding, &
       approval_id, edges_sha256)
    character(len=*), intent(in) :: source_id, source_sha256
    character(len=*), intent(in) :: source_commit_binding, approval_id, edges_sha256

    valid_identity = len_trim(source_id) > 0 .and. &
         len_trim(approval_id) > 0 .and. is_hex(source_sha256, 64) .and. &
         is_hex(source_commit_binding, 40) .and. is_hex(edges_sha256, 64)
  end function valid_identity

  logical function is_hex(value, required_length)
    character(len=*), intent(in) :: value
    integer, intent(in) :: required_length
    integer :: i

    is_hex = len_trim(value) == required_length
    if (.not. is_hex) return
    do i = 1, required_length
       if (index('0123456789abcdefABCDEF', value(i:i)) == 0) then
          is_hex = .false.
          return
       end if
    end do
  end function is_hex

  logical function snrt_spectral_status_is_known(status)
    character(len=*), intent(in) :: status

    snrt_spectral_status_is_known = trim(status) == snrt_spectral_status_reference .or. &
         trim(status) == snrt_spectral_status_candidate .or. &
         trim(status) == snrt_spectral_status_production
  end function snrt_spectral_status_is_known

  logical function snrt_spectral_fraction_semantics_is_known(semantics)
    character(len=*), intent(in) :: semantics

    snrt_spectral_fraction_semantics_is_known = &
         trim(semantics) == 'intrinsic' .or. trim(semantics) == 'escaped'
  end function snrt_spectral_fraction_semantics_is_known

  logical function runtime_status_allowed(status)
    character(len=*), intent(in) :: status
    character(len=16) :: opt_in
    integer :: opt_in_length, opt_in_status

    runtime_status_allowed = trim(status) == snrt_spectral_status_production
    if (trim(status) /= snrt_spectral_status_reference) return

    opt_in = ''
    opt_in_length = 0
    call get_environment_variable('SNRT_ALLOW_REFERENCE_CONTROL', opt_in, &
         length=opt_in_length, status=opt_in_status)
    ! Do not rely on Fortran short-circuit evaluation: an overlong value can
    ! report its full length while the destination is only len(opt_in).
    if (opt_in_status /= 0 .or. opt_in_length <= 0 .or. &
         opt_in_length > len(opt_in)) return
    runtime_status_allowed = trim(opt_in(1:opt_in_length)) == '1'
  end function runtime_status_allowed

  logical function snrt_spectral_contract_checkpoint_identity_matches(source_id, &
       source_sha256, source_commit_binding, approval_id, edges_sha256, &
       interval_convention, fraction_semantics, contract_status)
    character(len=*), intent(in) :: source_id, source_sha256
    character(len=*), intent(in) :: source_commit_binding, approval_id
    character(len=*), intent(in) :: edges_sha256, interval_convention
    character(len=*), intent(in) :: fraction_semantics, contract_status

    snrt_spectral_contract_checkpoint_identity_matches = &
         snrt_spectral_contract_loaded .and. &
         snrt_spectral_contract_runtime_allowed .and. &
         trim(source_id) == trim(snrt_spectral_contract_source_id) .and. &
         trim(source_sha256) == trim(snrt_spectral_contract_source_sha256) .and. &
         trim(source_commit_binding) == &
            trim(snrt_spectral_contract_source_commit_binding) .and. &
         trim(approval_id) == trim(snrt_spectral_contract_approval_id) .and. &
         trim(edges_sha256) == trim(snrt_spectral_contract_group_edges_sha256) .and. &
         trim(interval_convention) == &
            trim(snrt_spectral_contract_interval_convention) .and. &
         trim(fraction_semantics) == &
            trim(snrt_spectral_contract_fraction_semantics) .and. &
         trim(contract_status) == trim(snrt_spectral_contract_status)
  end function snrt_spectral_contract_checkpoint_identity_matches

  function snrt_spectral_contract_error_name(ierr) result(name)
    integer, intent(in) :: ierr
    character(len=40) :: name

    select case (ierr)
    case (snrt_spectral_contract_ok)
       name = 'ok'
    case (snrt_spectral_contract_err_missing)
       name = 'missing_environment_path'
    case (snrt_spectral_contract_err_open)
       name = 'file_open_failed'
    case (snrt_spectral_contract_err_read)
       name = 'namelist_read_failed'
    case (snrt_spectral_contract_err_version)
       name = 'unsupported_contract_version'
    case (snrt_spectral_contract_err_identity)
       name = 'invalid_source_identity'
    case (snrt_spectral_contract_err_status)
       name = 'unsupported_contract_status'
    case (snrt_spectral_contract_err_edges)
       name = 'canonical_edges_mismatch'
    case (snrt_spectral_contract_err_values)
       name = 'spectral_values_invalid'
    case (snrt_spectral_contract_err_not_loaded)
       name = 'contract_not_loaded'
    case (snrt_spectral_contract_err_fraction_semantics)
       name = 'invalid_fraction_semantics'
    case default
       name = 'unknown_error'
    end select
  end function snrt_spectral_contract_error_name

end module snrt_spectral_contract

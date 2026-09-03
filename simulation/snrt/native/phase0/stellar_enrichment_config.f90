! Phase 0 configuration for the stellar enrichment contract.
!
! The number and ordering of physical elements are compile-time properties.
! Runtime namelist options only enable or disable already allocated fields and
! source channels.  This keeps RAMSES nvar, MPI buffers, and restart layouts
! stable.

module stellar_enrichment_config
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer, parameter :: stellar_dp = kind(1.0d0)

  ! Change only for a separately compiled field layout.  The production
  ! Paper-IIIb layout uses all eleven physical species.
  integer, parameter :: n_stellar_elements = 11
  integer, parameter :: elem_h  = 1
  integer, parameter :: elem_he = 2
  integer, parameter :: elem_c  = 3
  integer, parameter :: elem_n  = 4
  integer, parameter :: elem_o  = 5
  integer, parameter :: elem_ne = 6
  integer, parameter :: elem_mg = 7
  integer, parameter :: elem_si = 8
  integer, parameter :: elem_s  = 9
  integer, parameter :: elem_ca = 10
  integer, parameter :: elem_fe = 11

  integer, parameter :: n_stellar_channels = 5
  integer, parameter :: channel_wind = 1
  integer, parameter :: channel_agb   = 2
  integer, parameter :: channel_snii  = 3
  integer, parameter :: channel_snia  = 4
  integer, parameter :: channel_pisn  = 5

  ! Configuration identifiers are kept in the shared contract module so the
  ! runtime cannot silently select an IMF that differs from the approved
  ! population-synthesis asset.
  integer, parameter, public :: stellar_imf_salpeter = 0
  integer, parameter, public :: stellar_imf_kroupa = 1
  integer, parameter, public :: stellar_imf_chabrier = 2
  integer, parameter, public :: stellar_imf_popiii = 3
  integer, parameter, public :: population_single_star_ssp = 0
  integer, parameter, public :: population_binary_ssp = 1
  integer, parameter, public :: enrichment_namelist_err_missing = 1005
  integer, parameter, public :: yield_basis_per_star_cumulative = 0
  integer, parameter, public :: yield_basis_per_event_cumulative = 1
  integer, parameter, public :: yield_basis_ssp_cumulative = 2
  integer, parameter, public :: yield_basis_ssp_rate = 3
  character(len=64), save, public :: stellar_fate_policy = &
       'review_only_unresolved'
  character(len=128), save, public :: stellar_fate_map_sha256 = ''
  character(len=128), save, public :: stellar_fate_approval_id = ''

  ! This is a diagnostic mirror of the JSON F-P1 partition.  It records the
  ! initial stellar mass whose terminal fate is not yet physically resolved;
  ! it is never added to returned mass or deposited as feedback.  Promotion
  ! requires the sidecar audit to prove that these bounds match the map.
  integer, parameter, public :: n_unresolved_fate_intervals = 2
  real(stellar_dp), dimension(n_unresolved_fate_intervals), save, public :: &
       unresolved_fate_mass_min = (/0.8d0, 40.0d0/)
  real(stellar_dp), dimension(n_unresolved_fate_intervals), save, public :: &
       unresolved_fate_mass_max = (/1.0d0, 120.0d0/)

  real(stellar_dp), dimension(n_stellar_channels), save, public :: &
       configured_channel_mass_min = &
       (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
  real(stellar_dp), dimension(n_stellar_channels), save, public :: &
       configured_channel_mass_max = &
       (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
  logical, dimension(n_stellar_channels), save, public :: &
       channel_owns_terminal_remnant = (/ .false., .true., .true., .false., .false. /)

  ! Runtime feedback implementation.  The channel-resolved path is the
  ! production default; legacy preserves the historical lagRamses behaviour
  ! for controlled reproduction of existing runs.
  character(len=32) :: stellar_feedback_mode = 'channel_resolved'
  integer :: default_imf_id = stellar_imf_kroupa
  integer :: population_model_id = population_single_star_ssp
  integer :: yield_source_basis_id = yield_basis_per_star_cumulative
  real(stellar_dp) :: configured_imf_mass_min = 0.08_stellar_dp
  real(stellar_dp) :: configured_imf_mass_max = 120.0_stellar_dp
  real(stellar_dp) :: configured_binary_fraction = 0.0_stellar_dp

  ! Namelist-controlled runtime switches.  The arrays remain allocated for
  ! every compile-time element regardless of these switches.
  logical :: active_element(n_stellar_elements) = .true.
  logical :: enable_wind = .true.
  logical :: enable_agb  = .true.
  logical :: enable_snii = .true.
  ! Type-Ia requires an explicit delay-time distribution.  A prompt SN table
  ! must not be interpreted as a Type-Ia history.
  logical :: enable_snia = .false.
  logical :: enable_pisn = .false.

contains

  subroutine set_enrichment_defaults()
    active_element = .true.
    enable_wind = .true.
    enable_agb  = .true.
    enable_snii = .true.
    enable_snia = .false.
    enable_pisn = .false.
    stellar_feedback_mode = 'channel_resolved'
    default_imf_id = stellar_imf_kroupa
    population_model_id = population_single_star_ssp
    yield_source_basis_id = yield_basis_per_star_cumulative
    configured_imf_mass_min = 0.08_stellar_dp
    configured_imf_mass_max = 120.0_stellar_dp
    configured_binary_fraction = 0.0_stellar_dp
    stellar_fate_policy = 'review_only_unresolved'
    stellar_fate_map_sha256 = ''
    stellar_fate_approval_id = ''
    unresolved_fate_mass_min = (/0.8d0, 40.0d0/)
    unresolved_fate_mass_max = (/1.0d0, 120.0d0/)
    configured_channel_mass_min = &
         (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
    configured_channel_mass_max = &
         (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
    channel_owns_terminal_remnant = &
         (/ .false., .true., .true., .false., .false. /)
  end subroutine set_enrichment_defaults

  subroutine read_enrichment_namelist(iunit, iostat_out)
    integer, intent(in) :: iunit
    integer, intent(out) :: iostat_out

    logical :: use_h, use_he, use_c, use_n, use_o, use_ne
    logical :: use_mg, use_si, use_s, use_ca, use_fe
    logical :: use_wind, use_agb, use_snii, use_snia, use_pisn
    character(len=32) :: feedback_mode
    character(len=32) :: population_model
    character(len=32) :: yield_source_basis
    character(len=32) :: parsed_feedback_mode
    integer :: imf_id, parsed_population_model_id, parsed_yield_source_basis_id
    real(stellar_dp) :: imf_mass_min_msun, imf_mass_max_msun, binary_fraction
    real(stellar_dp) :: channel_mass_min_msun(n_stellar_channels)
    real(stellar_dp) :: channel_mass_max_msun(n_stellar_channels)
    integer :: channel

    namelist /stellar_enrichment_params/ use_h, use_he, use_c, use_n, use_o, &
         use_ne, use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, &
         use_snii, use_snia, use_pisn, feedback_mode, imf_id, &
         population_model, yield_source_basis, imf_mass_min_msun, &
         imf_mass_max_msun, binary_fraction, channel_mass_min_msun, &
         channel_mass_max_msun

    use_h  = active_element(elem_h)
    use_he = active_element(elem_he)
    use_c  = active_element(elem_c)
    use_n  = active_element(elem_n)
    use_o  = active_element(elem_o)
    use_ne = active_element(elem_ne)
    use_mg = active_element(elem_mg)
    use_si = active_element(elem_si)
    use_s  = active_element(elem_s)
    use_ca = active_element(elem_ca)
    use_fe = active_element(elem_fe)
    use_wind = enable_wind
    use_agb  = enable_agb
    use_snii = enable_snii
    use_snia = enable_snia
    use_pisn = enable_pisn
    feedback_mode = ''
    imf_id = -1
    population_model = ''
    yield_source_basis = ''
    imf_mass_min_msun = -1.0_stellar_dp
    imf_mass_max_msun = -1.0_stellar_dp
    binary_fraction = -1.0_stellar_dp
    channel_mass_min_msun = -1.0_stellar_dp
    channel_mass_max_msun = -1.0_stellar_dp
    parsed_feedback_mode = stellar_feedback_mode
    parsed_population_model_id = population_model_id
    parsed_yield_source_basis_id = yield_source_basis_id

    read(iunit, nml=stellar_enrichment_params, iostat=iostat_out)
    if (iostat_out < 0) then
       iostat_out = enrichment_namelist_err_missing
       return
    end if
    if (iostat_out > 0) return

    call lowercase_ascii(feedback_mode)
    select case (trim(adjustl(feedback_mode)))
    case ('channel_resolved')
       parsed_feedback_mode = 'channel_resolved'
    case ('legacy')
       ! Legacy mode does not consume the population/yield-basis fields, but
       ! its element and channel switches still belong to this namelist.
       ! Commit them together only after the complete namelist read succeeds.
       call commit_runtime_switches(use_h, use_he, use_c, use_n, use_o, &
            use_ne, use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, &
            use_snii, use_snia, use_pisn)
       stellar_feedback_mode = 'legacy'
       return
    case default
       iostat_out = 1001
       return
    end select

    call lowercase_ascii(yield_source_basis)
    select case (trim(adjustl(yield_source_basis)))
    case ('per_star_cumulative')
       parsed_yield_source_basis_id = yield_basis_per_star_cumulative
    case ('per_event_cumulative')
       parsed_yield_source_basis_id = yield_basis_per_event_cumulative
    case ('ssp_cumulative_per_initial_mass')
       parsed_yield_source_basis_id = yield_basis_ssp_cumulative
    case ('ssp_rate_per_initial_mass')
       parsed_yield_source_basis_id = yield_basis_ssp_rate
    case default
       iostat_out = 1006
       return
    end select

    if (imf_id < stellar_imf_salpeter .or. imf_id > stellar_imf_popiii) then
       iostat_out = 1002
       return
    end if
    call lowercase_ascii(population_model)
    select case (trim(adjustl(population_model)))
    case ('single_star_ssp')
       parsed_population_model_id = population_single_star_ssp
    case ('binary_ssp')
       parsed_population_model_id = population_binary_ssp
    case default
       iostat_out = 1003
       return
    end select
    if (.not. ieee_is_finite(imf_mass_min_msun) .or. &
         .not. ieee_is_finite(imf_mass_max_msun) .or. &
         imf_mass_min_msun < 0.08_stellar_dp .or. &
         imf_mass_max_msun <= imf_mass_min_msun) then
       iostat_out = 1007
       return
    end if
    if (.not. ieee_is_finite(binary_fraction) .or. binary_fraction < 0.0_stellar_dp .or. &
         binary_fraction > 1.0_stellar_dp .or. &
         (parsed_population_model_id == population_single_star_ssp .and. &
          binary_fraction /= 0.0_stellar_dp) .or. &
         (parsed_population_model_id == population_binary_ssp .and. &
          binary_fraction <= 0.0_stellar_dp)) then
       iostat_out = 1008
       return
    end if
    do channel = 1, n_stellar_channels
       if (.not. ieee_is_finite(channel_mass_min_msun(channel)) .or. &
            .not. ieee_is_finite(channel_mass_max_msun(channel)) .or. &
            channel_mass_min_msun(channel) <= 0.0_stellar_dp .or. &
            channel_mass_max_msun(channel) <= channel_mass_min_msun(channel)) then
          iostat_out = 1004
          return
       end if
    end do
    do channel = channel_wind, channel_snii
       if ((channel == channel_wind .and. .not. use_wind) .or. &
            (channel == channel_agb .and. .not. use_agb) .or. &
            (channel == channel_snii .and. .not. use_snii)) cycle
       if (channel_mass_min_msun(channel) < imf_mass_min_msun .or. &
            channel_mass_max_msun(channel) > imf_mass_max_msun) then
          iostat_out = 1009
          return
       end if
    end do

    call commit_runtime_switches(use_h, use_he, use_c, use_n, use_o, use_ne, &
         use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, use_snii, &
         use_snia, use_pisn)
    stellar_feedback_mode = parsed_feedback_mode
    default_imf_id = imf_id
    population_model_id = parsed_population_model_id
    yield_source_basis_id = parsed_yield_source_basis_id
    configured_imf_mass_min = imf_mass_min_msun
    configured_imf_mass_max = imf_mass_max_msun
    configured_binary_fraction = binary_fraction
    configured_channel_mass_min = channel_mass_min_msun
    configured_channel_mass_max = channel_mass_max_msun
  end subroutine read_enrichment_namelist

  subroutine commit_runtime_switches(use_h, use_he, use_c, use_n, use_o, &
       use_ne, use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, &
       use_snii, use_snia, use_pisn)
    logical, intent(in) :: use_h, use_he, use_c, use_n, use_o, use_ne
    logical, intent(in) :: use_mg, use_si, use_s, use_ca, use_fe
    logical, intent(in) :: use_wind, use_agb, use_snii, use_snia, use_pisn

    active_element(elem_h)  = use_h
    active_element(elem_he) = use_he
    active_element(elem_c)  = use_c
    active_element(elem_n)  = use_n
    active_element(elem_o)  = use_o
    active_element(elem_ne) = use_ne
    active_element(elem_mg) = use_mg
    active_element(elem_si) = use_si
    active_element(elem_s)  = use_s
    active_element(elem_ca) = use_ca
    active_element(elem_fe) = use_fe
    enable_wind = use_wind
    enable_agb  = use_agb
    enable_snii = use_snii
    enable_snia = use_snia
    enable_pisn = use_pisn
  end subroutine commit_runtime_switches

  logical function use_channel_resolved_feedback()
    use_channel_resolved_feedback = &
         trim(stellar_feedback_mode) == 'channel_resolved'
  end function use_channel_resolved_feedback

  logical function production_source_model_supported()
    production_source_model_supported = &
         use_channel_resolved_feedback() .and. &
         population_model_id == population_single_star_ssp .and. &
         yield_source_basis_id == yield_basis_per_star_cumulative .and. &
         .not. enable_snia .and. .not. enable_pisn .and. &
         production_fate_policy_supported()
  end function production_source_model_supported

  logical function production_fate_policy_supported()
    ! F-P1 is deliberately fail-closed until a reviewed, complete terminal
    ! fate map and its immutable source package are promoted together.  The
    ! current interval map records the 0.8--1 and 40--120 Msun seams as
    ! unresolved; it is evidence, not a production admission token.
    production_fate_policy_supported = &
         trim(stellar_fate_policy) == 'approved_terminal_map_v1' .and. &
         valid_sha256(stellar_fate_map_sha256) .and. &
         len_trim(stellar_fate_approval_id) > 0
  end function production_fate_policy_supported

  logical function valid_sha256(value)
    character(len=*), intent(in) :: value
    integer :: i, code
    character :: digit

    valid_sha256 = len_trim(value) == 64
    if(.not.valid_sha256)return
    do i=1,64
       digit=value(i:i)
       code=iachar(digit)
       if(.not.((code>=iachar('0') .and. code<=iachar('9')) .or. &
            (code>=iachar('a') .and. code<=iachar('f')) .or. &
            (code>=iachar('A') .and. code<=iachar('F'))))then
          valid_sha256=.false.
          return
       endif
    end do
  end function valid_sha256

  function yield_source_basis_name() result(name)
    character(len=32) :: name

    select case (yield_source_basis_id)
    case (yield_basis_per_star_cumulative)
       name = 'per_star_cumulative'
    case (yield_basis_per_event_cumulative)
       name = 'per_event_cumulative'
    case (yield_basis_ssp_cumulative)
       name = 'ssp_cumulative_per_initial_mass'
    case (yield_basis_ssp_rate)
       name = 'ssp_rate_per_initial_mass'
    case default
       name = 'invalid'
    end select
  end function yield_source_basis_name

  subroutine lowercase_ascii(value)
    character(len=*), intent(inout) :: value
    integer :: i, code

    do i = 1, len_trim(value)
       code = iachar(value(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          value(i:i) = achar(code + iachar('a') - iachar('A'))
       end if
    end do
  end subroutine lowercase_ascii

end module stellar_enrichment_config

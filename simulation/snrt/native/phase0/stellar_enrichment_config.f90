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

  real(stellar_dp), dimension(n_stellar_channels), save, public :: &
       configured_channel_mass_min = &
       (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
  real(stellar_dp), dimension(n_stellar_channels), save, public :: &
       configured_channel_mass_max = &
       (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
  logical, dimension(n_stellar_channels), save, public :: &
       channel_owns_terminal_remnant = (/ .false., .true., .true., .false., .true. /)

  ! Runtime feedback implementation.  The channel-resolved path is the
  ! production default; legacy preserves the historical lagRamses behaviour
  ! for controlled reproduction of existing runs.
  character(len=32) :: stellar_feedback_mode = 'channel_resolved'
  integer :: default_imf_id = stellar_imf_kroupa
  integer :: population_model_id = population_single_star_ssp

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
    configured_channel_mass_min = &
         (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
    configured_channel_mass_max = &
         (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
    channel_owns_terminal_remnant = &
         (/ .false., .true., .true., .false., .true. /)
  end subroutine set_enrichment_defaults

  subroutine read_enrichment_namelist(iunit, iostat_out)
    integer, intent(in) :: iunit
    integer, intent(out) :: iostat_out

    logical :: use_h, use_he, use_c, use_n, use_o, use_ne
    logical :: use_mg, use_si, use_s, use_ca, use_fe
    logical :: use_wind, use_agb, use_snii, use_snia, use_pisn
    character(len=32) :: feedback_mode
    character(len=32) :: population_model
    integer :: imf_id
    real(stellar_dp) :: channel_mass_min_msun(n_stellar_channels)
    real(stellar_dp) :: channel_mass_max_msun(n_stellar_channels)
    integer :: channel

    namelist /stellar_enrichment_params/ use_h, use_he, use_c, use_n, use_o, &
         use_ne, use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, &
         use_snii, use_snia, use_pisn, feedback_mode, imf_id, &
         population_model, channel_mass_min_msun, channel_mass_max_msun

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
    feedback_mode = stellar_feedback_mode
    imf_id = default_imf_id
    population_model = 'single_star_ssp'
    if (population_model_id == population_binary_ssp) then
       population_model = 'binary_ssp'
    end if
    channel_mass_min_msun = configured_channel_mass_min
    channel_mass_max_msun = configured_channel_mass_max

    read(iunit, nml=stellar_enrichment_params, iostat=iostat_out)
    if (iostat_out /= 0) return

    call lowercase_ascii(feedback_mode)
    select case (trim(adjustl(feedback_mode)))
    case ('channel_resolved')
       stellar_feedback_mode = 'channel_resolved'
    case ('legacy')
       stellar_feedback_mode = 'legacy'
    case default
       iostat_out = 1001
       return
    end select

    if (imf_id < stellar_imf_salpeter .or. imf_id > stellar_imf_popiii) then
       iostat_out = 1002
       return
    end if
    call lowercase_ascii(population_model)
    select case (trim(adjustl(population_model)))
    case ('single_star_ssp')
       population_model_id = population_single_star_ssp
    case ('binary_ssp')
       population_model_id = population_binary_ssp
    case default
       iostat_out = 1003
       return
    end select
    do channel = 1, n_stellar_channels
       if (.not. ieee_is_finite(channel_mass_min_msun(channel)) .or. &
            .not. ieee_is_finite(channel_mass_max_msun(channel)) .or. &
            channel_mass_min_msun(channel) <= 0.0_stellar_dp .or. &
            channel_mass_max_msun(channel) <= channel_mass_min_msun(channel)) then
          iostat_out = 1004
          return
       end if
    end do

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
    default_imf_id = imf_id
    configured_channel_mass_min = channel_mass_min_msun
    configured_channel_mass_max = channel_mass_max_msun
  end subroutine read_enrichment_namelist

  logical function use_channel_resolved_feedback()
    use_channel_resolved_feedback = &
         trim(stellar_feedback_mode) == 'channel_resolved'
  end function use_channel_resolved_feedback

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

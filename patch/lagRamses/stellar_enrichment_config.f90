! Phase 0 configuration for the stellar enrichment contract.
!
! The number and ordering of physical elements are compile-time properties.
! Runtime namelist options only enable or disable already allocated fields and
! source channels.  This keeps RAMSES nvar, MPI buffers, and restart layouts
! stable.

module stellar_enrichment_config
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

  ! Runtime feedback implementation.  The channel-resolved path is the
  ! production default; legacy preserves the historical lagRamses behaviour
  ! for controlled reproduction of existing runs.
  character(len=32) :: stellar_feedback_mode = 'channel_resolved'

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
  end subroutine set_enrichment_defaults

  subroutine read_enrichment_namelist(iunit, iostat_out)
    integer, intent(in) :: iunit
    integer, intent(out) :: iostat_out

    logical :: use_h, use_he, use_c, use_n, use_o, use_ne
    logical :: use_mg, use_si, use_s, use_ca, use_fe
    logical :: use_wind, use_agb, use_snii, use_snia, use_pisn
    character(len=32) :: feedback_mode

    namelist /stellar_enrichment_params/ use_h, use_he, use_c, use_n, use_o, &
         use_ne, use_mg, use_si, use_s, use_ca, use_fe, use_wind, use_agb, &
         use_snii, use_snia, use_pisn, feedback_mode

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

! Shared AGN efficiency and supplied-inflow convention.
!
! This module is deliberately independent of the RAMSES sink implementation.
! Both the coarse-state writer and the SNRT source driver call this routine so
! that the coefficient used for Lbol and photons cannot silently diverge.
module snrt_agn_efficiency
  use amr_parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  private
  public :: snrt_agn_resolve_efficiency
  public :: snrt_agn_rt_requested
  public :: snrt_agn_efficiency_status_name
  public :: snrt_agn_efficiency_mode_name
  public :: snrt_agn_eff_status_ok
  public :: snrt_agn_eff_status_spin_disabled_default
  public :: snrt_agn_eff_status_spin_uninitialized
  public :: snrt_agn_eff_status_raw_nonfinite
  public :: snrt_agn_eff_status_raw_nonpositive
  public :: snrt_agn_eff_status_raw_ge_one
  public :: snrt_agn_eff_status_rate_nonfinite
  public :: snrt_agn_eff_status_rate_clipped
  public :: snrt_agn_eff_status_zero_eddington
  public :: snrt_agn_eff_status_floor_disabled
  public :: snrt_agn_eff_status_floor_nonfinite
  public :: snrt_agn_eff_status_effective_invalid
  public :: snrt_agn_eff_mode_thermal
  public :: snrt_agn_eff_mode_mad_high
  public :: snrt_agn_eff_mode_mad_quenched
  public :: snrt_agn_eff_mode_mad_floor_disabled
  public :: snrt_agn_eff_mode_invalid

  ! Status values are bit flags so a diagnostic can preserve, for example,
  ! rate clipping together with a controlled zero-Eddington state.
  integer, parameter :: snrt_agn_eff_status_ok = 0
  integer, parameter :: snrt_agn_eff_status_spin_disabled_default = 1
  integer, parameter :: snrt_agn_eff_status_spin_uninitialized = 2
  integer, parameter :: snrt_agn_eff_status_raw_nonfinite = 4
  integer, parameter :: snrt_agn_eff_status_raw_nonpositive = 8
  integer, parameter :: snrt_agn_eff_status_raw_ge_one = 16
  integer, parameter :: snrt_agn_eff_status_rate_nonfinite = 32
  integer, parameter :: snrt_agn_eff_status_rate_clipped = 64
  integer, parameter :: snrt_agn_eff_status_zero_eddington = 128
  integer, parameter :: snrt_agn_eff_status_floor_disabled = 256
  integer, parameter :: snrt_agn_eff_status_floor_nonfinite = 512
  integer, parameter :: snrt_agn_eff_status_effective_invalid = 1024

  integer, parameter :: snrt_agn_eff_mode_thermal = 0
  integer, parameter :: snrt_agn_eff_mode_mad_high = 1
  integer, parameter :: snrt_agn_eff_mode_mad_quenched = 2
  integer, parameter :: snrt_agn_eff_mode_mad_floor_disabled = 3
  integer, parameter :: snrt_agn_eff_mode_invalid = 9

contains

  logical function snrt_agn_rt_requested() result(requested)
    ! One process-lifetime latch, shared by namelist preflight, legacy
    ! dispatch and the live driver. Preflight checks MPI agreement explicitly.
    logical, save :: resolved=.false., latched=.false.
    character(len=8) :: value
    integer :: length, status
    if (.not.resolved) then
       value=''
       call get_environment_variable('SNRT_RT_ENABLE',value,length=length,status=status)
       latched=status==0 .and. length==1 .and. value(1:1)=='1'
       resolved=.true.
    end if
    requested=latched
  end function snrt_agn_rt_requested

  pure subroutine snrt_agn_resolve_efficiency(raw_efficiency, spin_bh, bondi_rate, &
       eddington_rate, mad_jet, x_floor, resolved_base_efficiency, &
       effective_efficiency, inflow_rate, edd_ratio, status, source_allowed, &
       mode)
    real(dp), intent(in) :: raw_efficiency, bondi_rate, eddington_rate, x_floor
    logical, intent(in) :: spin_bh, mad_jet
    real(dp), intent(out) :: resolved_base_efficiency, effective_efficiency
    real(dp), intent(out) :: inflow_rate, edd_ratio
    integer, intent(out) :: status, mode
    logical, intent(out) :: source_allowed
    real(dp), parameter :: default_efficiency = 0.1d0
    real(dp) :: bondi_nonnegative, eddington_nonnegative
    logical :: rates_valid, floor_valid

    resolved_base_efficiency = default_efficiency
    effective_efficiency = 0.0d0
    inflow_rate = 0.0d0
    edd_ratio = 0.0d0
    status = snrt_agn_eff_status_ok
    mode = snrt_agn_eff_mode_thermal
    source_allowed = .true.

    ! The spin-disabled branch is deliberate model behaviour.  The raw array
    ! is not consulted in this branch, so a stale/uninitialised raw value does
    ! not change the selected default; the status remains visible in output.
    if (.not. spin_bh) then
       status = ior(status, snrt_agn_eff_status_spin_disabled_default)
    else if (.not. ieee_is_finite(raw_efficiency)) then
       status = ior(status, snrt_agn_eff_status_raw_nonfinite)
       source_allowed = .false.
    else if (raw_efficiency == 0.0d0) then
       ! pm/init_sink.f90 starts eps_sink at zero.  The old accretion path
       ! reads that zero before kjhan_growspin writes the first value.  Keep a
       ! finite review value for the ledger but make the divergence explicit.
       status = ior(status, snrt_agn_eff_status_spin_uninitialized)
       source_allowed = .false.
    else if (raw_efficiency < 0.0d0) then
       status = ior(status, snrt_agn_eff_status_raw_nonpositive)
       source_allowed = .false.
    else if (raw_efficiency >= 1.0d0) then
       status = ior(status, snrt_agn_eff_status_raw_ge_one)
       resolved_base_efficiency = 0.0d0
       source_allowed = .false.
    else
       resolved_base_efficiency = raw_efficiency
    end if

    if (.not. ieee_is_finite(bondi_rate) .or. &
         .not. ieee_is_finite(eddington_rate)) then
       status = ior(status, snrt_agn_eff_status_rate_nonfinite)
       source_allowed = .false.
       rates_valid = .false.
       bondi_nonnegative = 0.0d0
       eddington_nonnegative = 0.0d0
    else
       rates_valid = .true.
       bondi_nonnegative = max(bondi_rate, 0.0d0)
       eddington_nonnegative = max(eddington_rate, 0.0d0)
       if (bondi_rate < 0.0d0 .or. eddington_rate < 0.0d0) then
          status = ior(status, snrt_agn_eff_status_rate_clipped)
          source_allowed = .false.
       end if
    end if

    ! Supplied inflow is the Bondi supply limited by the Eddington cap.  The
    ! non-negative projection is part of the shared convention; negative
    ! inputs are nevertheless marked non-promotable above.
    inflow_rate = min(bondi_nonnegative, eddington_nonnegative)
    if (rates_valid .and. eddington_rate > 0.0d0) then
       edd_ratio = bondi_rate / eddington_rate
    else
       edd_ratio = 0.0d0
       if (rates_valid) status = ior(status, snrt_agn_eff_status_zero_eddington)
    end if

    floor_valid = ieee_is_finite(x_floor)
    if (.not. floor_valid) then
       status = ior(status, snrt_agn_eff_status_floor_nonfinite)
       source_allowed = .false.
    else if (mad_jet .and. x_floor <= 0.0d0) then
       ! Match the legacy branch (no MAD reduction) while making an invalid
       ! floor visible and non-promotable in the live source path.
       status = ior(status, snrt_agn_eff_status_floor_disabled)
       source_allowed = .false.
    end if

    effective_efficiency = resolved_base_efficiency
    if (mad_jet .and. floor_valid .and. x_floor > 0.0d0) then
       if (edd_ratio < x_floor) then
          mode = snrt_agn_eff_mode_mad_quenched
          effective_efficiency = resolved_base_efficiency * &
               max(edd_ratio, 0.0d0) / x_floor
       else
          mode = snrt_agn_eff_mode_mad_high
       end if
    else if (mad_jet .and. floor_valid) then
       mode = snrt_agn_eff_mode_mad_floor_disabled
    end if

    if (.not. ieee_is_finite(resolved_base_efficiency) .or. &
         resolved_base_efficiency <= 0.0d0 .or. &
         resolved_base_efficiency >= 1.0d0 .or. &
         .not. ieee_is_finite(effective_efficiency) .or. &
         effective_efficiency < 0.0d0 .or. effective_efficiency >= 1.0d0) then
       status = ior(status, snrt_agn_eff_status_effective_invalid)
       effective_efficiency = 0.0d0
       mode = snrt_agn_eff_mode_invalid
       source_allowed = .false.
    end if
  end subroutine snrt_agn_resolve_efficiency

  function snrt_agn_efficiency_status_name(status) result(name)
    integer, intent(in) :: status
    character(len=256) :: name
    logical :: first

    name = ''
    first = .true.
    if (status == snrt_agn_eff_status_ok) then
       name = 'ok'
       return
    end if
    if (iand(status, snrt_agn_eff_status_spin_disabled_default) /= 0) &
         call append_name('spin_disabled_default')
    if (iand(status, snrt_agn_eff_status_spin_uninitialized) /= 0) &
         call append_name('spin_enabled_uninitialized')
    if (iand(status, snrt_agn_eff_status_raw_nonfinite) /= 0) &
         call append_name('raw_nonfinite')
    if (iand(status, snrt_agn_eff_status_raw_nonpositive) /= 0) &
         call append_name('raw_nonpositive')
    if (iand(status, snrt_agn_eff_status_raw_ge_one) /= 0) &
         call append_name('raw_ge_one')
    if (iand(status, snrt_agn_eff_status_rate_nonfinite) /= 0) &
         call append_name('rate_nonfinite')
    if (iand(status, snrt_agn_eff_status_rate_clipped) /= 0) &
         call append_name('rate_clipped')
    if (iand(status, snrt_agn_eff_status_zero_eddington) /= 0) &
         call append_name('zero_eddington')
    if (iand(status, snrt_agn_eff_status_floor_disabled) /= 0) &
         call append_name('floor_disabled')
    if (iand(status, snrt_agn_eff_status_floor_nonfinite) /= 0) &
         call append_name('floor_nonfinite')
    if (iand(status, snrt_agn_eff_status_effective_invalid) /= 0) &
         call append_name('effective_invalid')

  contains

    subroutine append_name(label)
      character(len=*), intent(in) :: label
      if (.not. first) name = trim(name)//'|'
      name = trim(name)//label
      first = .false.
    end subroutine append_name
  end function snrt_agn_efficiency_status_name

  function snrt_agn_efficiency_mode_name(mode) result(name)
    integer, intent(in) :: mode
    character(len=32) :: name

    select case (mode)
    case (snrt_agn_eff_mode_thermal)
       name = 'THERMAL'
    case (snrt_agn_eff_mode_mad_high)
       name = 'MAD_HIGH'
    case (snrt_agn_eff_mode_mad_quenched)
       name = 'MAD_QUENCHED'
    case (snrt_agn_eff_mode_mad_floor_disabled)
       name = 'MAD_FLOOR_DISABLED'
    case default
       name = 'INVALID'
    end select
  end function snrt_agn_efficiency_mode_name

end module snrt_agn_efficiency

program snrt_agn_efficiency_smoke
  use amr_parameters, only: dp
  use snrt_agn_efficiency
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none

  real(dp) :: base, effective, inflow, ratio
  integer :: status, mode
  logical :: allowed

  call resolve(0.1d0,.true.,2.0d0,1.0d0,.false.,0.01d0,0.1d0,0.1d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_thermal,.true.)
  call resolve(0.1d0,.true.,0.001d0,1.0d0,.true.,0.01d0,0.1d0,0.01d0, &
       0.001d0,0.001d0,snrt_agn_eff_mode_mad_quenched,.true.)
  call resolve(0.1d0,.true.,0.2d0,1.0d0,.true.,0.01d0,0.1d0,0.1d0, &
       0.2d0,0.2d0,snrt_agn_eff_mode_mad_high,.true.)
  ! The boundary is high-state by the declared strict '< X_floor' rule.
  call resolve(0.1d0,.true.,0.01d0,1.0d0,.true.,0.01d0,0.1d0,0.1d0, &
       0.01d0,0.01d0,snrt_agn_eff_mode_mad_high,.true.)
  call resolve(0.1d0,.true.,0.2d0,0.0d0,.true.,0.01d0,0.1d0,0.0d0, &
       0.0d0,0.0d0,snrt_agn_eff_mode_mad_quenched,.true.)
  call resolve(0.0d0,.false.,2.0d0,1.0d0,.false.,0.01d0,0.1d0,0.1d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_thermal,.true., &
       snrt_agn_eff_status_spin_disabled_default)
  call resolve(0.0d0,.true.,2.0d0,1.0d0,.false.,0.01d0,0.1d0,0.1d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_thermal,.false., &
       snrt_agn_eff_status_spin_uninitialized)
  call resolve(-0.1d0,.true.,2.0d0,1.0d0,.false.,0.01d0,0.1d0,0.1d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_thermal,.false., &
       snrt_agn_eff_status_raw_nonpositive)
  call resolve(1.0d0,.true.,2.0d0,1.0d0,.false.,0.01d0,0.0d0,0.0d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_invalid,.false., &
       snrt_agn_eff_status_raw_ge_one)
  call resolve(0.1d0,.true.,ieee_value(0.0d0,ieee_quiet_nan),1.0d0, &
       .false.,0.01d0,0.1d0,0.1d0,0.0d0,0.0d0,snrt_agn_eff_mode_thermal, &
       .false.,snrt_agn_eff_status_rate_nonfinite)
  call resolve(0.1d0,.true.,-1.0d0,1.0d0,.false.,0.01d0,0.1d0,0.1d0, &
       0.0d0,-1.0d0,snrt_agn_eff_mode_thermal,.false., &
       snrt_agn_eff_status_rate_clipped)
  call resolve(0.1d0,.true.,2.0d0,1.0d0,.true.,0.0d0,0.1d0,0.1d0, &
       1.0d0,2.0d0,snrt_agn_eff_mode_mad_floor_disabled,.false., &
       snrt_agn_eff_status_floor_disabled)

  write(*,'(A)') 'SNRT_AGN_EFFICIENCY_OK thermal=1 mad_low=0.01 boundary=high '&
       //'spin_init=visible invalid=fail_closed'

contains

  subroutine resolve(raw, spin, bondi, edd, mad, floor, expected_base, &
       expected_effective, expected_inflow, expected_ratio, expected_mode, &
       expected_allowed, expected_flag)
    real(dp), intent(in) :: raw, bondi, edd, floor
    logical, intent(in) :: spin, mad, expected_allowed
    real(dp), intent(in) :: expected_base, expected_effective, expected_inflow, expected_ratio
    integer, intent(in) :: expected_mode
    integer, intent(in), optional :: expected_flag
    integer :: expected_status

    call snrt_agn_resolve_efficiency(raw,spin,bondi,edd,mad,floor,base,effective, &
         inflow,ratio,status,allowed,mode)
    if (abs(base-expected_base) > 1.0d-13 .or. &
         abs(effective-expected_effective) > 1.0d-13 .or. &
         abs(inflow-expected_inflow) > 1.0d-13 .or. &
         abs(ratio-expected_ratio) > 1.0d-13 .or. &
         mode /= expected_mode .or. allowed .neqv. expected_allowed) error stop 1
    if (.not. ieee_is_finite(base) .or. .not. ieee_is_finite(effective) .or. &
         effective < 0.0d0 .or. effective >= 1.0d0) error stop 2
    expected_status = snrt_agn_eff_status_ok
    if (present(expected_flag)) expected_status = expected_flag
    if (iand(status,expected_status) /= expected_status) error stop 3
  end subroutine resolve
end program snrt_agn_efficiency_smoke

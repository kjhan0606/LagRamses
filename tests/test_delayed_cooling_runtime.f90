program test_delayed_cooling_runtime
  use iso_fortran_env, only: real64
  use delayed_cooling_runtime
  implicit none
  real(real64) :: marker, remaining
  real(real64), parameter :: rho_floor = 1.0e-30_real64

  if (delayed_cooling_model_id('first_event_fixed') /= delayed_cooling_first_event_fixed) error stop 1
  if (delayed_cooling_model_id('legacy_accumulating') /= delayed_cooling_legacy_accumulating) error stop 2
  if (delayed_cooling_model_id('unknown') /= -1) error stop 3

  ! The first SN event arms exactly t_diss after the gas mass return.
  marker = delayed_cooling_fixed_after_source(0.0_real64, 2.0_real64, &
       2.5_real64, .true., 20.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 2.5_real64, rho_floor)
  if (abs(remaining - 20.0_real64) > 1.0e-13_real64) error stop 4

  ! Let time pass before a second SN. It must retain the existing deadline,
  ! not restart a fresh t_diss interval, even after the returned gas changes
  ! the density used to encode the remaining time.
  marker = delayed_cooling_fixed_after_step(marker, 2.5_real64, 3.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 2.5_real64, rho_floor)
  if (abs(remaining - 17.0_real64) > 1.0e-13_real64) error stop 5
  marker = delayed_cooling_fixed_after_source(marker, 2.5_real64, &
       3.0_real64, .true., 20.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 3.0_real64, rho_floor)
  if (abs(remaining - 17.0_real64) > 1.0e-13_real64) error stop 6

  ! Non-SN mass return also preserves an active deadline, while never arming
  ! an inactive timer.
  marker = delayed_cooling_fixed_after_source(marker, 3.0_real64, &
       4.0_real64, .false., 20.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 4.0_real64, rho_floor)
  if (abs(remaining - 17.0_real64) > 1.0e-13_real64) error stop 7
  marker = delayed_cooling_fixed_after_step(marker, 4.0_real64, 5.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 4.0_real64, rho_floor)
  if (abs(remaining - 12.0_real64) > 1.0e-13_real64) error stop 8
  marker = delayed_cooling_fixed_after_source(0.0_real64, 3.0_real64, &
       4.0_real64, .false., 20.0_real64, rho_floor)
  if (marker /= 0.0_real64) error stop 9

  marker = delayed_cooling_fixed_after_step(marker, 4.0_real64, 12.0_real64, rho_floor)
  if (marker /= 0.0_real64) error stop 10
  marker = delayed_cooling_fixed_after_source(marker, 4.0_real64, &
       4.5_real64, .true., 20.0_real64, rho_floor)
  remaining = delayed_cooling_remaining_myr(marker, 4.5_real64, rho_floor)
  if (abs(remaining - 20.0_real64) > 1.0e-13_real64) error stop 11

  write(*,'(A)') 'delayed cooling runtime: PASS'
end program test_delayed_cooling_runtime

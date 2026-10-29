module delayed_cooling_runtime
  use iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: delayed_cooling_legacy_accumulating = 0
  integer, parameter, public :: delayed_cooling_first_event_fixed = 1

  public :: delayed_cooling_model_id
  public :: delayed_cooling_remaining_myr
  public :: delayed_cooling_fixed_after_source
  public :: delayed_cooling_fixed_after_step

contains

  pure integer function delayed_cooling_model_id(name)
    character(len=*), intent(in) :: name

    select case (trim(adjustl(name)))
    case ('legacy_accumulating')
       delayed_cooling_model_id = delayed_cooling_legacy_accumulating
    case ('first_event_fixed')
       delayed_cooling_model_id = delayed_cooling_first_event_fixed
    case default
       delayed_cooling_model_id = -1
    end select
  end function delayed_cooling_model_id

  pure real(real64) function delayed_cooling_remaining_myr(marker, rho, rho_floor)
    real(real64), intent(in) :: marker, rho, rho_floor
    real(real64) :: safe_rho

    delayed_cooling_remaining_myr = 0.0_real64
    if (marker <= 0.0_real64) return
    safe_rho = max(rho, max(rho_floor, tiny(1.0_real64)))
    delayed_cooling_remaining_myr = marker / safe_rho
  end function delayed_cooling_remaining_myr

  pure real(real64) function delayed_cooling_fixed_after_source( &
       marker, rho_before, rho_after, shock_event, t_diss_myr, rho_floor)
    real(real64), intent(in) :: marker, rho_before, rho_after, t_diss_myr, rho_floor
    logical, intent(in) :: shock_event
    real(real64) :: remaining_myr, safe_rho_after

    safe_rho_after = max(rho_after, max(rho_floor, tiny(1.0_real64)))
    remaining_myr = delayed_cooling_remaining_myr(marker, rho_before, rho_floor)
    if (remaining_myr > 0.0_real64) then
       delayed_cooling_fixed_after_source = remaining_myr * safe_rho_after
    else if (shock_event .and. t_diss_myr > 0.0_real64) then
       delayed_cooling_fixed_after_source = t_diss_myr * safe_rho_after
    else
       delayed_cooling_fixed_after_source = 0.0_real64
    end if
  end function delayed_cooling_fixed_after_source

  pure real(real64) function delayed_cooling_fixed_after_step( &
       marker, rho, dt_myr, rho_floor)
    real(real64), intent(in) :: marker, rho, dt_myr, rho_floor
    real(real64) :: remaining_myr, safe_rho

    safe_rho = max(rho, max(rho_floor, tiny(1.0_real64)))
    remaining_myr = delayed_cooling_remaining_myr(marker, safe_rho, rho_floor)
    delayed_cooling_fixed_after_step = &
         max(remaining_myr - max(dt_myr, 0.0_real64), 0.0_real64) * safe_rho
  end function delayed_cooling_fixed_after_step

end module delayed_cooling_runtime

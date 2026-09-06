! Explicit unit conversions at the RAMSES/stellar-source boundary.
!
! The yield table stores ages in yr on disk and in Gyr in memory.  RAMSES
! supplies a dimensionless code time.  The cosmological time convention used
! by the current lagRamses patch is retained here, in one auditable place:
!
!   age_gyr = code_age * scale_t / (seconds_per_gyr * aexp**2)

module stellar_native_units
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp
  implicit none

  private

  real(stellar_dp), parameter, public :: seconds_per_gyr = &
       1.0d9 * 365.25d0 * 86400.0d0
  real(stellar_dp), parameter, public :: solar_mass_cgs = 1.98847d33

  integer, parameter, public :: units_ok = 0
  integer, parameter, public :: units_err_argument = 1
  integer, parameter, public :: units_err_nonfinite = 2
  integer, parameter, public :: units_err_result = 4

  public :: code_time_to_age_gyr
  public :: code_interval_to_age_gyr
  public :: mass_code_to_msun
  public :: mass_msun_to_code
  public :: momentum_cgs_to_code
  public :: energy_erg_to_code

contains

  subroutine code_time_to_age_gyr(code_time, scale_t, aexp, age_gyr, ierr)
    real(stellar_dp), intent(in) :: code_time, scale_t, aexp
    real(stellar_dp), intent(out) :: age_gyr
    integer, intent(out) :: ierr

    call convert_code_time(code_time, scale_t, aexp, age_gyr, ierr)
  end subroutine code_time_to_age_gyr

  subroutine code_interval_to_age_gyr(code_interval, scale_t, aexp, dt_gyr, &
       ierr)
    real(stellar_dp), intent(in) :: code_interval, scale_t, aexp
    real(stellar_dp), intent(out) :: dt_gyr
    integer, intent(out) :: ierr

    call convert_code_time(code_interval, scale_t, aexp, dt_gyr, ierr)
  end subroutine code_interval_to_age_gyr

  subroutine convert_code_time(code_time, scale_t, aexp, result, ierr)
    real(stellar_dp), intent(in) :: code_time, scale_t, aexp
    real(stellar_dp), intent(out) :: result
    integer, intent(out) :: ierr
    real(stellar_dp) :: denominator

    result = 0.0_stellar_dp
    ierr = units_ok
    if (.not. ieee_is_finite(code_time) .or. .not. ieee_is_finite(scale_t) &
         .or. .not. ieee_is_finite(aexp)) then
       ierr = units_err_nonfinite
       return
    end if
    if (code_time < 0.0_stellar_dp .or. scale_t <= 0.0_stellar_dp .or. &
         aexp <= 0.0_stellar_dp) then
       ierr = units_err_argument
       return
    end if

    ! Keep aexp**2 explicit: this is the current lagRamses convention.
    denominator = seconds_per_gyr * aexp**2
    result = code_time * scale_t / denominator
    if (.not. ieee_is_finite(denominator) .or. &
         .not. ieee_is_finite(result)) then
       ierr = units_err_result
       result = 0.0_stellar_dp
    end if
  end subroutine convert_code_time

  subroutine mass_code_to_msun(mass_code, scale_mass, mass_msun, ierr)
    real(stellar_dp), intent(in) :: mass_code, scale_mass
    real(stellar_dp), intent(out) :: mass_msun
    integer, intent(out) :: ierr

    call scale_value(mass_code, scale_mass, mass_msun, ierr)
  end subroutine mass_code_to_msun

  subroutine mass_msun_to_code(mass_msun, scale_mass, mass_code, ierr)
    real(stellar_dp), intent(in) :: mass_msun, scale_mass
    real(stellar_dp), intent(out) :: mass_code
    integer, intent(out) :: ierr

    call divide_value(mass_msun, scale_mass, mass_code, ierr)
  end subroutine mass_msun_to_code

  subroutine momentum_cgs_to_code(momentum_cgs, scale_momentum, momentum_code, &
       ierr)
    real(stellar_dp), intent(in) :: momentum_cgs(3), scale_momentum
    real(stellar_dp), intent(out) :: momentum_code(3)
    integer, intent(out) :: ierr
    integer :: i

    momentum_code = 0.0_stellar_dp
    ierr = units_ok
    if (scale_momentum <= 0.0_stellar_dp .or. &
         .not. ieee_is_finite(scale_momentum)) then
       if (.not. ieee_is_finite(scale_momentum)) then
          ierr = units_err_nonfinite
       else
          ierr = units_err_argument
       end if
       return
    end if
    do i = 1, 3
       if (.not. ieee_is_finite(momentum_cgs(i))) then
          ierr = units_err_nonfinite
          return
       end if
       momentum_code(i) = momentum_cgs(i) / scale_momentum
       if (.not. ieee_is_finite(momentum_code(i))) then
          ierr = units_err_result
          momentum_code = 0.0_stellar_dp
          return
       end if
    end do
  end subroutine momentum_cgs_to_code

  subroutine energy_erg_to_code(energy_erg, scale_energy, energy_code, ierr)
    real(stellar_dp), intent(in) :: energy_erg, scale_energy
    real(stellar_dp), intent(out) :: energy_code
    integer, intent(out) :: ierr

    call divide_value(energy_erg, scale_energy, energy_code, ierr)
  end subroutine energy_erg_to_code

  subroutine scale_value(value, multiplier, result, ierr)
    real(stellar_dp), intent(in) :: value, multiplier
    real(stellar_dp), intent(out) :: result
    integer, intent(out) :: ierr

    result = 0.0_stellar_dp
    ierr = units_ok
    if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(multiplier)) then
       ierr = units_err_nonfinite
       return
    end if
    if (multiplier <= 0.0_stellar_dp) then
       ierr = units_err_argument
       return
    end if
    result = value * multiplier
    if (.not. ieee_is_finite(result)) then
       ierr = units_err_result
       result = 0.0_stellar_dp
    end if
  end subroutine scale_value

  subroutine divide_value(value, divisor, result, ierr)
    real(stellar_dp), intent(in) :: value, divisor
    real(stellar_dp), intent(out) :: result
    integer, intent(out) :: ierr

    result = 0.0_stellar_dp
    ierr = units_ok
    if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(divisor)) then
       ierr = units_err_nonfinite
       return
    end if
    if (divisor <= 0.0_stellar_dp) then
       ierr = units_err_argument
       return
    end if
    result = value / divisor
    if (.not. ieee_is_finite(result)) then
       ierr = units_err_result
       result = 0.0_stellar_dp
    end if
  end subroutine divide_value

end module stellar_native_units

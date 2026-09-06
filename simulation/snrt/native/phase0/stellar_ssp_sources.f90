! Phase 0 SSP source integrator.
!
! A table row describes one initial star.  This module integrates those rows
! over the IMF to obtain cumulative source quantities for one star particle,
! which is treated as a single-age stellar population.  It deliberately does
! not deposit into AMR cells and does not infer a Type-Ia DTD from a prompt
! SN table.

module stellar_ssp_sources
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, stellar_imf_salpeter, stellar_imf_kroupa, &
       stellar_imf_chabrier, stellar_imf_popiii, &
       yield_basis_per_star_cumulative
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, clear_cumulative
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_yield_provider, only: evaluate_channel_cumulative, provider_ok
  implicit none

  private

  ! These values match the legacy C IMF enumeration.
  integer, parameter, public :: imf_salpeter = stellar_imf_salpeter
  integer, parameter, public :: imf_kroupa = stellar_imf_kroupa
  integer, parameter, public :: imf_chabrier = stellar_imf_chabrier
  integer, parameter, public :: imf_popiii = stellar_imf_popiii

  integer, parameter, public :: ssp_source_ok = 0
  integer, parameter, public :: ssp_source_err_argument = 1
  integer, parameter, public :: ssp_source_err_imf = 2
  integer, parameter, public :: ssp_source_err_provider = 3
  integer, parameter, public :: ssp_source_err_nonfinite = 4
  integer, parameter, public :: ssp_source_err_basis = 5

  real(stellar_dp), parameter :: chabrier_high_amplitude = &
       exp(-(log10(1.0_stellar_dp)-log10(0.079_stellar_dp))**2 / &
       (2.0_stellar_dp*0.69_stellar_dp**2))

  public :: integrate_ssp_channel
  public :: calculate_imf_normalization
  public :: calculate_imf_mass_fraction
  public :: evaluate_imf

contains

  subroutine integrate_ssp_channel(table, population, channel_id, age_gyr, &
       mass_min, mass_max, n_mass_bins, state, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    integer, intent(in) :: channel_id, n_mass_bins
    real(stellar_dp), intent(in) :: age_gyr, mass_min, mass_max
    type(stellar_cumulative_t), intent(out) :: state
    integer, intent(out) :: ierr

    real(stellar_dp) :: imf_norm, log_min, log_max, dlog_mass
    real(stellar_dp) :: log_mass, mass, dm, n_stars
    real(stellar_dp) :: star_weight
    type(stellar_cumulative_t) :: star_state
    integer :: bin, provider_ierr

    call clear_cumulative(state)
    ierr = ssp_source_ok

    if (.not. ieee_is_finite(population%initial_mass) .or. &
         .not. ieee_is_finite(population%birth_metallicity) .or. &
         .not. ieee_is_finite(population%imf_mass_min) .or. &
         .not. ieee_is_finite(population%imf_mass_max) .or. &
         .not. ieee_is_finite(age_gyr) .or. &
         population%initial_mass <= 0.0_stellar_dp .or. age_gyr < 0.0_stellar_dp .or. &
         mass_min <= 0.0_stellar_dp .or. mass_max <= mass_min .or. &
         population%imf_mass_min <= 0.0_stellar_dp .or. &
         population%imf_mass_max <= population%imf_mass_min .or. &
         mass_min < population%imf_mass_min .or. &
         mass_max > population%imf_mass_max .or. &
         n_mass_bins <= 0 .or. channel_id < 1 .or. &
         channel_id > n_stellar_channels) then
       ierr = ssp_source_err_argument
       return
    end if
    if (population%yield_basis_id /= yield_basis_per_star_cumulative) then
       ierr = ssp_source_err_basis
       return
    end if

    call calculate_imf_normalization(population%imf_id, &
         population%imf_mass_min, population%imf_mass_max, imf_norm, &
         provider_ierr)
    if (provider_ierr /= 0) then
       ierr = ssp_source_err_imf
       return
    end if

    log_min = log10(mass_min)
    log_max = log10(mass_max)
    dlog_mass = (log_max - log_min) / real(n_mass_bins, stellar_dp)

    do bin = 1, n_mass_bins
       log_mass = log_min + (real(bin, stellar_dp) - 0.5_stellar_dp) * dlog_mass
       mass = 10.0_stellar_dp ** log_mass
       dm = mass * log(10.0_stellar_dp) * dlog_mass
       n_stars = population%initial_mass * imf_norm * &
            evaluate_imf(mass, population%imf_id) * dm
       if (n_stars <= 0.0_stellar_dp) cycle

       call evaluate_channel_cumulative(table, channel_id, mass, &
            population%birth_metallicity, age_gyr, star_state, provider_ierr)
       if (provider_ierr /= provider_ok) then
          ierr = ssp_source_err_provider
          return
       end if

       star_weight = n_stars
       state%ejected_mass = state%ejected_mass + &
            star_weight * star_state%ejected_mass
       state%net_yield = state%net_yield + star_weight * star_state%net_yield
       state%returned_mass = state%returned_mass + &
            star_weight * star_state%returned_mass
       state%remnant_mass = state%remnant_mass + &
            star_weight * star_state%remnant_mass
       state%energy = state%energy + star_weight * star_state%energy
       state%momentum = state%momentum + star_weight * star_state%momentum
       state%channel_returned_mass = state%channel_returned_mass + &
            star_weight * star_state%channel_returned_mass
       state%channel_energy = state%channel_energy + &
            star_weight * star_state%channel_energy
       state%channel_momentum = state%channel_momentum + &
            star_weight * star_state%channel_momentum
       state%channel_ejected_mass = state%channel_ejected_mass + &
            star_weight * star_state%channel_ejected_mass
       state%channel_net_yield = state%channel_net_yield + &
            star_weight * star_state%channel_net_yield
    end do

    if (.not. cumulative_values_finite(state)) then
       call clear_cumulative(state)
       ierr = ssp_source_err_nonfinite
       return
    end if

    ! The remnant value is retained per channel so the population ledger can
    ! enforce terminal-channel ownership.  Living mass remains a population-
    ! level quantity and is computed only after all channel states are known.
  end subroutine integrate_ssp_channel

  subroutine calculate_imf_normalization(imf_id, mass_min, mass_max, &
       normalization, ierr)
    integer, intent(in) :: imf_id
    real(stellar_dp), intent(in) :: mass_min, mass_max
    real(stellar_dp), intent(out) :: normalization
    integer, intent(out) :: ierr

    real(stellar_dp) :: integral

    ierr = 0
    normalization = 0.0_stellar_dp
    if (.not. ieee_is_finite(mass_min) .or. .not. ieee_is_finite(mass_max) .or. &
         mass_min < 0.08_stellar_dp .or. mass_max <= mass_min) then
       ierr = 1
       return
    end if
    select case (imf_id)
    case (imf_salpeter)
       integral = integrate_power_law(mass_min, mass_max, -1.35_stellar_dp, &
            1.0_stellar_dp)
    case (imf_kroupa)
       integral = 0.0_stellar_dp
       if (mass_min < 0.5_stellar_dp) integral = integral + &
            integrate_power_law(mass_min, min(mass_max,0.5_stellar_dp), &
            -0.3_stellar_dp, 2.0_stellar_dp)
       if (mass_max > 0.5_stellar_dp) integral = integral + &
            integrate_power_law(max(mass_min,0.5_stellar_dp), mass_max, &
            -1.3_stellar_dp, 1.0_stellar_dp)
    case (imf_chabrier)
       integral = 0.0_stellar_dp
       if (mass_min < 1.0_stellar_dp) integral = integral + &
            integrate_chabrier_low(mass_min, min(mass_max,1.0_stellar_dp))
       if (mass_max > 1.0_stellar_dp) integral = integral + &
            integrate_power_law(max(mass_min,1.0_stellar_dp), mass_max, &
            -1.3_stellar_dp, chabrier_high_amplitude)
    case (imf_popiii)
       integral = 0.0_stellar_dp
       if (mass_max > 10.0_stellar_dp .and. mass_min < 100.0_stellar_dp) &
            integral = integral + integrate_power_law(max(mass_min,10.0_stellar_dp), &
            min(mass_max,100.0_stellar_dp), 1.5_stellar_dp, 0.1_stellar_dp)
       if (mass_max > 100.0_stellar_dp) integral = integral + &
            100.0_stellar_dp * (mass_max-max(mass_min,100.0_stellar_dp))
    case default
       ierr = 1
       return
    end select

    if (.not. ieee_is_finite(integral) .or. integral <= 0.0_stellar_dp) then
       ierr = 1
       return
    end if
    normalization = 1.0_stellar_dp / integral
  end subroutine calculate_imf_normalization

  subroutine calculate_imf_mass_fraction(imf_id, mass_min, mass_max, &
       interval_min, interval_max, fraction, ierr)
    ! Return the initial-mass fraction in one explicitly supplied interval.
    ! The interval must already be clipped to the configured IMF support;
    ! silently clipping it here would hide a fate-map/domain mismatch.
    integer, intent(in) :: imf_id
    real(stellar_dp), intent(in) :: mass_min, mass_max
    real(stellar_dp), intent(in) :: interval_min, interval_max
    real(stellar_dp), intent(out) :: fraction
    integer, intent(out) :: ierr

    real(stellar_dp) :: normalization, integral

    fraction = 0.0_stellar_dp
    ierr = 0
    if (.not. ieee_is_finite(mass_min) .or. .not. ieee_is_finite(mass_max) .or. &
         .not. ieee_is_finite(interval_min) .or. &
         .not. ieee_is_finite(interval_max) .or. mass_min < 0.08_stellar_dp .or. &
         mass_max <= mass_min .or. interval_max <= interval_min .or. &
         interval_min < mass_min .or. interval_max > mass_max) then
       ierr = 1
       return
    end if

    call calculate_imf_normalization(imf_id, mass_min, mass_max, normalization, ierr)
    if (ierr /= 0) return

    select case (imf_id)
    case (imf_salpeter)
       integral = integrate_power_law(interval_min, interval_max, &
            -1.35_stellar_dp, 1.0_stellar_dp)
    case (imf_kroupa)
       integral = 0.0_stellar_dp
       if (interval_min < 0.5_stellar_dp) integral = integral + &
            integrate_power_law(interval_min, min(interval_max,0.5_stellar_dp), &
            -0.3_stellar_dp, 2.0_stellar_dp)
       if (interval_max > 0.5_stellar_dp) integral = integral + &
            integrate_power_law(max(interval_min,0.5_stellar_dp), interval_max, &
            -1.3_stellar_dp, 1.0_stellar_dp)
    case (imf_chabrier)
       integral = 0.0_stellar_dp
       if (interval_min < 1.0_stellar_dp) integral = integral + &
            integrate_chabrier_low(interval_min, min(interval_max,1.0_stellar_dp))
       if (interval_max > 1.0_stellar_dp) integral = integral + &
            integrate_power_law(max(interval_min,1.0_stellar_dp), interval_max, &
            -1.3_stellar_dp, chabrier_high_amplitude)
    case (imf_popiii)
       integral = 0.0_stellar_dp
       if (interval_max > 10.0_stellar_dp .and. interval_min < 100.0_stellar_dp) &
            integral = integral + integrate_power_law(max(interval_min,10.0_stellar_dp), &
            min(interval_max,100.0_stellar_dp), 1.5_stellar_dp, 0.1_stellar_dp)
       if (interval_max > 100.0_stellar_dp) integral = integral + &
            100.0_stellar_dp * (interval_max-max(interval_min,100.0_stellar_dp))
    case default
       ierr = 1
       return
    end select

    fraction = normalization * integral
    if (.not. ieee_is_finite(fraction) .or. fraction < 0.0_stellar_dp .or. &
         fraction > 1.0_stellar_dp + 1.0e-12_stellar_dp) then
       fraction = 0.0_stellar_dp
       ierr = 1
    end if
  end subroutine calculate_imf_mass_fraction

  pure real(stellar_dp) function integrate_power_law(mass_min, mass_max, &
       exponent, amplitude)
    real(stellar_dp), intent(in) :: mass_min, mass_max, exponent, amplitude

    integrate_power_law = 0.0_stellar_dp
    if (mass_max <= mass_min) return
    integrate_power_law = amplitude * &
         (mass_max**(exponent+1.0_stellar_dp) - &
          mass_min**(exponent+1.0_stellar_dp)) / &
         (exponent+1.0_stellar_dp)
  end function integrate_power_law

  pure real(stellar_dp) function integrate_chabrier_low(mass_min, mass_max)
    real(stellar_dp), intent(in) :: mass_min, mass_max
    real(stellar_dp), parameter :: log_mean = log10(0.079_stellar_dp)
    real(stellar_dp), parameter :: sigma = 0.69_stellar_dp
    real(stellar_dp), parameter :: ln10 = log(10.0_stellar_dp)
    real(stellar_dp), parameter :: sqrt_two = sqrt(2.0_stellar_dp)
    real(stellar_dp), parameter :: sqrt_pi_over_two = &
         sqrt(acos(-1.0_stellar_dp)/2.0_stellar_dp)
    real(stellar_dp) :: shifted_mean, prefactor

    integrate_chabrier_low = 0.0_stellar_dp
    if (mass_max <= mass_min) return
    shifted_mean = log_mean + ln10*sigma**2
    prefactor = ln10 * exp(ln10*log_mean + &
         0.5_stellar_dp*(ln10*sigma)**2) * sigma * sqrt_pi_over_two
    integrate_chabrier_low = prefactor * &
         (erf((log10(mass_max)-shifted_mean)/(sqrt_two*sigma)) - &
          erf((log10(mass_min)-shifted_mean)/(sqrt_two*sigma)))
  end function integrate_chabrier_low

  real(stellar_dp) function evaluate_imf(mass, imf_id)
    real(stellar_dp), intent(in) :: mass
    integer, intent(in) :: imf_id

    evaluate_imf = 0.0_stellar_dp
    if (mass < 0.08_stellar_dp) return

    select case (imf_id)
    case (imf_salpeter)
       evaluate_imf = mass ** (-2.35_stellar_dp)
    case (imf_kroupa)
       if (mass < 0.5_stellar_dp) then
          evaluate_imf = (mass / 0.5_stellar_dp) ** (-1.3_stellar_dp) * &
               0.5_stellar_dp ** (-2.3_stellar_dp)
       else
          evaluate_imf = mass ** (-2.3_stellar_dp)
       end if
    case (imf_chabrier)
       if (mass < 1.0_stellar_dp) then
          evaluate_imf = (1.0_stellar_dp / mass) * &
               exp(-((log10(mass) - log10(0.079_stellar_dp)) ** 2) / &
               (2.0_stellar_dp * 0.69_stellar_dp ** 2))
       else
          evaluate_imf = chabrier_high_amplitude * mass ** (-2.3_stellar_dp)
       end if
    case (imf_popiii)
       if (mass < 10.0_stellar_dp) then
          evaluate_imf = 0.0_stellar_dp
       else if (mass < 100.0_stellar_dp) then
          evaluate_imf = (mass / 100.0_stellar_dp) ** 0.5_stellar_dp
       else
          evaluate_imf = (mass / 100.0_stellar_dp) ** (-1.0_stellar_dp)
       end if
    end select
  end function evaluate_imf

  logical function cumulative_values_finite(state)
    type(stellar_cumulative_t), intent(in) :: state
    integer :: i, j

    cumulative_values_finite = ieee_is_finite(state%returned_mass) .and. &
         ieee_is_finite(state%remnant_mass) .and. &
         ieee_is_finite(state%living_mass) .and. &
         ieee_is_finite(state%energy)
    do i = 1, 3
       cumulative_values_finite = cumulative_values_finite .and. &
            ieee_is_finite(state%momentum(i))
    end do
    do i = 1, n_stellar_elements
       cumulative_values_finite = cumulative_values_finite .and. &
            ieee_is_finite(state%ejected_mass(i)) .and. &
            ieee_is_finite(state%net_yield(i))
    end do
    do i = 1, n_stellar_channels
       cumulative_values_finite = cumulative_values_finite .and. &
            ieee_is_finite(state%channel_returned_mass(i)) .and. &
            ieee_is_finite(state%channel_energy(i))
       do j = 1, 3
          cumulative_values_finite = cumulative_values_finite .and. &
               ieee_is_finite(state%channel_momentum(i,j))
       end do
       do j = 1, n_stellar_elements
          cumulative_values_finite = cumulative_values_finite .and. &
               ieee_is_finite(state%channel_ejected_mass(i,j)) .and. &
               ieee_is_finite(state%channel_net_yield(i,j))
       end do
    end do
  end function cumulative_values_finite

end module stellar_ssp_sources

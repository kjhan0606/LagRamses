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
       stellar_imf_chabrier, stellar_imf_popiii
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

  public :: integrate_ssp_channel

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
         .not. ieee_is_finite(age_gyr) .or. &
         population%initial_mass <= 0.0_stellar_dp .or. age_gyr < 0.0_stellar_dp .or. &
         mass_min <= 0.0_stellar_dp .or. mass_max <= mass_min .or. &
         n_mass_bins <= 0 .or. channel_id < 1 .or. &
         channel_id > n_stellar_channels) then
       ierr = ssp_source_err_argument
       return
    end if

    call calculate_imf_normalization(population%imf_id, imf_norm, provider_ierr)
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

  subroutine calculate_imf_normalization(imf_id, normalization, ierr)
    integer, intent(in) :: imf_id
    real(stellar_dp), intent(out) :: normalization
    integer, intent(out) :: ierr

    integer, parameter :: n_bins = 1024
    real(stellar_dp) :: mass_min, mass_max, dlog_mass
    real(stellar_dp) :: mass, dm, integral, log_mass
    integer :: bin

    ierr = 0
    normalization = 0.0_stellar_dp
    select case (imf_id)
    case (imf_salpeter, imf_kroupa, imf_chabrier)
       mass_min = 0.08_stellar_dp
       mass_max = 120.0_stellar_dp
    case (imf_popiii)
       mass_min = 10.0_stellar_dp
       mass_max = 300.0_stellar_dp
    case default
       ierr = 1
       return
    end select

    dlog_mass = (log10(mass_max) - log10(mass_min)) / &
         real(n_bins, stellar_dp)
    integral = 0.0_stellar_dp
    do bin = 1, n_bins
       log_mass = log10(mass_min) + &
            (real(bin, stellar_dp) - 0.5_stellar_dp) * dlog_mass
       mass = 10.0_stellar_dp ** log_mass
       dm = mass * log(10.0_stellar_dp) * dlog_mass
       integral = integral + mass * evaluate_imf(mass, imf_id) * dm
    end do

    if (integral <= 0.0_stellar_dp) then
       ierr = 1
       return
    end if
    normalization = 1.0_stellar_dp / integral
  end subroutine calculate_imf_normalization

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
          evaluate_imf = mass ** (-2.3_stellar_dp)
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

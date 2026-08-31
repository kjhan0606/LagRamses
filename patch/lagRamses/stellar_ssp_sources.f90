! Phase 0 SSP source integrator.
!
! A table row describes one initial star.  This module integrates those rows
! over the IMF to obtain cumulative source quantities for one star particle,
! which is treated as a single-age stellar population.  It deliberately does
! not deposit into AMR cells and does not infer a Type-Ia DTD from a prompt
! SN table.

module stellar_ssp_sources
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, clear_cumulative
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_yield_provider, only: evaluate_channel_cumulative, provider_ok
  implicit none

  private

  ! These values match the legacy C IMF enumeration.
  integer, parameter, public :: imf_salpeter = 0
  integer, parameter, public :: imf_kroupa = 1
  integer, parameter, public :: imf_chabrier = 2
  integer, parameter, public :: imf_popiii = 3

  integer, parameter, public :: ssp_source_ok = 0
  integer, parameter, public :: ssp_source_err_argument = 1
  integer, parameter, public :: ssp_source_err_imf = 2
  integer, parameter, public :: ssp_source_err_provider = 3

  public :: integrate_ssp_channel

contains

  subroutine integrate_ssp_channel(table, population, channel_id, age, &
       mass_min, mass_max, n_mass_bins, state, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    integer, intent(in) :: channel_id, n_mass_bins
    real(stellar_dp), intent(in) :: age, mass_min, mass_max
    type(stellar_cumulative_t), intent(out) :: state
    integer, intent(out) :: ierr

    real(stellar_dp) :: imf_norm, log_min, log_max, dlog_mass
    real(stellar_dp) :: log_mass, mass, dm, n_stars
    real(stellar_dp) :: star_weight
    type(stellar_cumulative_t) :: star_state
    integer :: bin, provider_ierr

    call clear_cumulative(state)
    ierr = ssp_source_ok

    if (population%initial_mass <= 0.0_stellar_dp .or. age < 0.0_stellar_dp .or. &
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
            population%birth_metallicity, age, star_state, provider_ierr)
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

    ! Remnant and living masses are population-level ledger quantities.  They
    ! are intentionally left unset here because a channel-specific table row
    ! must not cause the same remnant to be counted once per channel.
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

end module stellar_ssp_sources

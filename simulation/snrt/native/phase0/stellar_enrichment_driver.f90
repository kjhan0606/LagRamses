! Phase 0 top-level stellar source driver.
!
! The caller supplies the physically appropriate initial-mass interval for
! each channel.  This is intentional: AGB, SNII, SNIa, and PISN do not share
! the same progenitor or binary normalization.  The driver combines only
! already-differenced timestep sources and does not perform cell deposition.

module stellar_enrichment_driver
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_channels, &
       channel_wind, channel_agb, channel_snii, channel_snia, channel_pisn, &
       enable_wind, enable_agb, enable_snii, enable_snia, enable_pisn, &
       channel_owns_terminal_remnant
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, stellar_source_t, clear_cumulative, clear_source
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_source_increment, only: integrate_ssp_channel_increment, &
       source_increment_ok
  use stellar_ssp_sources, only: integrate_ssp_channel, ssp_source_ok
  use stellar_population_ledger, only: stellar_population_ledger_t, &
       clear_population_ledger, finalize_population_ledger, population_ledger_ok
  implicit none

  private
  integer, parameter, public :: enrichment_driver_ok = 0
  integer, parameter, public :: enrichment_driver_err_argument = 1
  integer, parameter, public :: enrichment_driver_err_channel = 2
  integer, parameter, public :: enrichment_driver_err_source = 3
  integer, parameter, public :: enrichment_driver_err_unsupported = 4
  integer, parameter, public :: enrichment_driver_err_ledger = 8

  public :: compute_stellar_source_increment
  public :: compute_stellar_cumulative

contains

  subroutine compute_stellar_source_increment(table, population, age_gyr, &
       dt_gyr, &
       channel_mass_min, channel_mass_max, n_mass_bins, source, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    real(stellar_dp), intent(in) :: age_gyr, dt_gyr
    real(stellar_dp), intent(in) :: channel_mass_min(n_stellar_channels)
    real(stellar_dp), intent(in) :: channel_mass_max(n_stellar_channels)
    integer, intent(in) :: n_mass_bins
    type(stellar_source_t), intent(out) :: source
    integer, intent(out) :: ierr

    type(stellar_source_t) :: channel_source
    integer :: channel, source_ierr

    call clear_source(source)
    ierr = enrichment_driver_ok

    if (.not. ieee_is_finite(age_gyr) .or. .not. ieee_is_finite(dt_gyr) .or. &
         age_gyr < 0.0_stellar_dp .or. dt_gyr < 0.0_stellar_dp .or. &
         n_mass_bins <= 0) then
       ierr = enrichment_driver_err_argument
       return
    end if

    do channel = 1, n_stellar_channels
       if (.not. channel_is_enabled(channel)) cycle
       if (channel == channel_snia .or. channel == channel_pisn) then
          ! SNIa needs a DTD convolution and PISN needs an explicit
          ! population/core-mass gate; neither is an IMF-only SSP channel.
          ierr = enrichment_driver_err_unsupported
          return
       end if
       if (channel_mass_min(channel) <= 0.0_stellar_dp .or. &
            channel_mass_max(channel) <= channel_mass_min(channel)) then
          ierr = enrichment_driver_err_channel
          return
       end if

       call integrate_ssp_channel_increment(table, population, channel, &
            age_gyr, dt_gyr, channel_mass_min(channel), &
            channel_mass_max(channel), &
            n_mass_bins, channel_source, source_ierr)
       if (source_ierr /= source_increment_ok) then
          ierr = enrichment_driver_err_source
          return
       end if
       call add_source(source, channel_source)
    end do
  end subroutine compute_stellar_source_increment

  subroutine compute_stellar_cumulative(table, population, age_gyr, &
       channel_mass_min, channel_mass_max, n_mass_bins, channel_states, &
       ledger, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    real(stellar_dp), intent(in) :: age_gyr
    real(stellar_dp), intent(in) :: channel_mass_min(n_stellar_channels)
    real(stellar_dp), intent(in) :: channel_mass_max(n_stellar_channels)
    integer, intent(in) :: n_mass_bins
    type(stellar_cumulative_t), intent(out) :: channel_states(n_stellar_channels)
    type(stellar_population_ledger_t), intent(out) :: ledger
    integer, intent(out) :: ierr

    logical :: channel_enabled(n_stellar_channels)
    integer :: channel, ssp_ierr, ledger_ierr

    ierr = enrichment_driver_ok
    channel_enabled = .false.
    call clear_population_ledger(ledger)
    do channel = 1, n_stellar_channels
       call clear_cumulative(channel_states(channel))
    end do

    if (.not. ieee_is_finite(age_gyr) .or. age_gyr < 0.0_stellar_dp .or. &
         n_mass_bins <= 0) then
       ierr = enrichment_driver_err_argument
       call finalize_population_ledger(population, channel_states, &
            channel_enabled, channel_owns_terminal_remnant, 1.0e-10_stellar_dp, &
            ledger, ledger_ierr)
       return
    end if

    do channel = 1, n_stellar_channels
       if (.not. channel_is_enabled(channel)) cycle
       if (channel == channel_snia .or. channel == channel_pisn) then
          ierr = enrichment_driver_err_unsupported
          return
       end if
       if (.not. ieee_is_finite(channel_mass_min(channel)) .or. &
            .not. ieee_is_finite(channel_mass_max(channel)) .or. &
            channel_mass_min(channel) <= 0.0_stellar_dp .or. &
            channel_mass_max(channel) <= channel_mass_min(channel)) then
          ierr = enrichment_driver_err_channel
          return
       end if

       call integrate_ssp_channel(table, population, channel, age_gyr, &
            channel_mass_min(channel), channel_mass_max(channel), n_mass_bins, &
            channel_states(channel), ssp_ierr)
       if (ssp_ierr /= ssp_source_ok) then
          ierr = enrichment_driver_err_source
          return
       end if
       channel_enabled(channel) = .true.
    end do

    call finalize_population_ledger(population, channel_states, channel_enabled, &
         channel_owns_terminal_remnant, 1.0e-10_stellar_dp, ledger, ledger_ierr)
    if (ledger_ierr /= population_ledger_ok) then
       ierr = enrichment_driver_err_ledger
    end if
  end subroutine compute_stellar_cumulative

  logical function channel_is_enabled(channel)
    integer, intent(in) :: channel

    select case (channel)
    case (channel_wind)
       channel_is_enabled = enable_wind
    case (channel_agb)
       channel_is_enabled = enable_agb
    case (channel_snii)
       channel_is_enabled = enable_snii
    case (channel_snia)
       channel_is_enabled = enable_snia
    case (channel_pisn)
       channel_is_enabled = enable_pisn
    case default
       channel_is_enabled = .false.
    end select
  end function channel_is_enabled

  subroutine add_source(total, component)
    type(stellar_source_t), intent(inout) :: total
    type(stellar_source_t), intent(in) :: component

    total%ejected_mass = total%ejected_mass + component%ejected_mass
    total%net_yield = total%net_yield + component%net_yield
    total%returned_mass = total%returned_mass + component%returned_mass
    total%energy = total%energy + component%energy
    total%momentum = total%momentum + component%momentum
    total%channel_returned_mass = total%channel_returned_mass + &
         component%channel_returned_mass
    total%channel_energy = total%channel_energy + component%channel_energy
    total%channel_momentum = total%channel_momentum + &
         component%channel_momentum
    total%channel_ejected_mass = total%channel_ejected_mass + &
         component%channel_ejected_mass
    total%channel_net_yield = total%channel_net_yield + &
         component%channel_net_yield
  end subroutine add_source

end module stellar_enrichment_driver

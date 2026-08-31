! Phase 0 timestep source increment.
!
! This module converts cumulative SSP quantities into the incremental source
! for one hydrodynamic timestep.  It is intentionally independent of the AMR
! deposition code; the returned stellar_source_t is the only object that the
! deposition layer needs to consume.

module stellar_source_increment
  use stellar_enrichment_config, only: stellar_dp
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t, stellar_source_t, clear_source, &
       cumulative_difference
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_ssp_sources, only: integrate_ssp_channel, ssp_source_ok
  implicit none

  private
  integer, parameter, public :: source_increment_ok = 0
  integer, parameter, public :: source_increment_err_argument = 1
  integer, parameter, public :: source_increment_err_ssp = 2

  public :: integrate_ssp_channel_increment

contains

  subroutine integrate_ssp_channel_increment(table, population, channel_id, &
       age, timestep, mass_min, mass_max, n_mass_bins, source, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    integer, intent(in) :: channel_id, n_mass_bins
    real(stellar_dp), intent(in) :: age, timestep, mass_min, mass_max
    type(stellar_source_t), intent(out) :: source
    integer, intent(out) :: ierr

    type(stellar_cumulative_t) :: previous, later
    integer :: ssp_ierr

    call clear_source(source)
    ierr = source_increment_ok

    if (age < 0.0_stellar_dp .or. timestep < 0.0_stellar_dp .or. &
         mass_min <= 0.0_stellar_dp .or. mass_max <= mass_min .or. &
         n_mass_bins <= 0) then
       ierr = source_increment_err_argument
       return
    end if

    if (timestep == 0.0_stellar_dp) return

    call integrate_ssp_channel(table, population, channel_id, age, mass_min, &
         mass_max, n_mass_bins, previous, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       ierr = source_increment_err_ssp
       return
    end if

    call integrate_ssp_channel(table, population, channel_id, &
         age + timestep, mass_min, mass_max, n_mass_bins, later, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       ierr = source_increment_err_ssp
       return
    end if

    call cumulative_difference(later, previous, source)
  end subroutine integrate_ssp_channel_increment

end module stellar_source_increment

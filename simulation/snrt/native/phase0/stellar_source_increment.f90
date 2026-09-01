! Phase 0 timestep source increment.
!
! This module converts cumulative SSP quantities into the incremental source
! for one hydrodynamic timestep.  It is intentionally independent of the AMR
! deposition code; the returned stellar_source_t is the only object that the
! deposition layer needs to consume.

module stellar_source_increment
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
  integer, parameter, public :: source_increment_err_nonfinite = 4
  integer, parameter, public :: source_increment_err_negative = 8

  public :: integrate_ssp_channel_increment

contains

  subroutine integrate_ssp_channel_increment(table, population, channel_id, &
       age_gyr, timestep_gyr, mass_min, mass_max, n_mass_bins, source, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    type(stellar_population_t), intent(in) :: population
    integer, intent(in) :: channel_id, n_mass_bins
    real(stellar_dp), intent(in) :: age_gyr, timestep_gyr, mass_min, mass_max
    type(stellar_source_t), intent(out) :: source
    integer, intent(out) :: ierr

    type(stellar_cumulative_t) :: previous, later
    integer :: ssp_ierr

    call clear_source(source)
    ierr = source_increment_ok

    if (.not. ieee_is_finite(age_gyr) .or. &
         .not. ieee_is_finite(timestep_gyr) .or. &
         .not. ieee_is_finite(mass_min) .or. .not. ieee_is_finite(mass_max) .or. &
         age_gyr < 0.0_stellar_dp .or. timestep_gyr < 0.0_stellar_dp .or. &
         mass_min <= 0.0_stellar_dp .or. mass_max <= mass_min .or. &
         n_mass_bins <= 0) then
       ierr = source_increment_err_argument
       return
    end if

    if (timestep_gyr == 0.0_stellar_dp) return

    call integrate_ssp_channel(table, population, channel_id, age_gyr, mass_min, &
         mass_max, n_mass_bins, previous, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       ierr = source_increment_err_ssp
       return
    end if

    call integrate_ssp_channel(table, population, channel_id, &
         age_gyr + timestep_gyr, mass_min, mass_max, n_mass_bins, later, ssp_ierr)
    if (ssp_ierr /= ssp_source_ok) then
       ierr = source_increment_err_ssp
       return
    end if

    call cumulative_difference(later, previous, source)
    if (.not. source_values_finite(source)) then
       call clear_source(source)
       ierr = source_increment_err_nonfinite
       return
    end if
    if (.not. source_values_nonnegative(source)) then
       call clear_source(source)
       ierr = source_increment_err_negative
       return
    end if
  end subroutine integrate_ssp_channel_increment

  logical function source_values_finite(source)
    type(stellar_source_t), intent(in) :: source
    integer :: i, j

    source_values_finite = ieee_is_finite(source%returned_mass) .and. &
         ieee_is_finite(source%energy)
    do i = 1, 3
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%momentum(i))
    end do
    do i = 1, size(source%ejected_mass)
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%ejected_mass(i)) .and. &
            ieee_is_finite(source%net_yield(i))
    end do
    do i = 1, size(source%channel_returned_mass)
       source_values_finite = source_values_finite .and. &
            ieee_is_finite(source%channel_returned_mass(i)) .and. &
            ieee_is_finite(source%channel_energy(i))
       do j = 1, 3
          source_values_finite = source_values_finite .and. &
               ieee_is_finite(source%channel_momentum(i,j))
       end do
       do j = 1, size(source%channel_ejected_mass, 2)
          source_values_finite = source_values_finite .and. &
               ieee_is_finite(source%channel_ejected_mass(i,j)) .and. &
               ieee_is_finite(source%channel_net_yield(i,j))
       end do
    end do
  end function source_values_finite

  logical function source_values_nonnegative(source)
    type(stellar_source_t), intent(in) :: source
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp
    real(stellar_dp) :: scale, tracked_mass, channel_tracked_mass
    integer :: channel, element

    source_values_nonnegative = source%returned_mass >= -tolerance .and. &
         source%energy >= -tolerance .and. &
         minval(source%ejected_mass) >= -tolerance .and. &
         minval(source%channel_returned_mass) >= -tolerance .and. &
         minval(source%channel_energy) >= -tolerance .and. &
         minval(source%channel_ejected_mass) >= -tolerance
    if (.not. source_values_nonnegative) return

    tracked_mass = sum(source%ejected_mass)
    scale = max(1.0_stellar_dp, abs(source%returned_mass), abs(tracked_mass))
    if (tracked_mass > source%returned_mass + tolerance * scale) then
       source_values_nonnegative = .false.
       return
    end if
    if (abs(sum(source%channel_returned_mass) - source%returned_mass) > &
         tolerance * scale) then
       source_values_nonnegative = .false.
       return
    end if
    do element = 1, size(source%ejected_mass)
       scale = max(1.0_stellar_dp, abs(source%ejected_mass(element)), &
            abs(sum(source%channel_ejected_mass(:,element))))
       if (abs(sum(source%channel_ejected_mass(:,element)) - &
            source%ejected_mass(element)) > tolerance * scale) then
          source_values_nonnegative = .false.
          return
       end if
    end do
    do channel = 1, size(source%channel_returned_mass)
       channel_tracked_mass = sum(source%channel_ejected_mass(channel,:))
       scale = max(1.0_stellar_dp, &
            abs(source%channel_returned_mass(channel)), &
            abs(channel_tracked_mass))
       if (channel_tracked_mass > source%channel_returned_mass(channel) + &
            tolerance * scale) then
          source_values_nonnegative = .false.
          return
       end if
    end do
  end function source_values_nonnegative

end module stellar_source_increment

! Phase 0 adapter from one interpolated table row to the common cumulative
! stellar state.  SSP integration and timestep differencing are kept outside
! this module.

module stellar_yield_provider
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels
  use stellar_enrichment_contract, only: stellar_cumulative_t, &
       clear_cumulative
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_yield_interpolation, only: interpolate_yield_row, &
       interpolation_ok
  implicit none

  private
  integer, parameter, public :: provider_ok = 0
  integer, parameter, public :: provider_err_argument = 1
  integer, parameter, public :: provider_err_interpolation = 2
  integer, parameter, public :: provider_err_nonfinite = 4

  public :: evaluate_channel_cumulative

contains

  subroutine evaluate_channel_cumulative(table, channel_id, initial_mass, &
       birth_metallicity, age_gyr, state, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: initial_mass, birth_metallicity, age_gyr
    type(stellar_cumulative_t), intent(out) :: state
    integer, intent(out) :: ierr

    real(stellar_dp) :: returned_mass, remnant_mass, energy
    real(stellar_dp) :: momentum(3)
    real(stellar_dp) :: ejected_mass(n_stellar_elements)
    real(stellar_dp) :: net_yield(n_stellar_elements)
    integer :: interpolation_ierr

    call clear_cumulative(state)
    ierr = provider_ok

    if (channel_id < 1 .or. channel_id > n_stellar_channels .or. &
         initial_mass <= 0.0_stellar_dp .or. birth_metallicity < 0.0_stellar_dp .or. &
         age_gyr < 0.0_stellar_dp .or. &
         .not. ieee_is_finite(initial_mass) .or. &
         .not. ieee_is_finite(birth_metallicity) .or. &
         .not. ieee_is_finite(age_gyr)) then
       ierr = provider_err_argument
       return
    end if

    call interpolate_yield_row(table, channel_id, initial_mass, &
         birth_metallicity, age_gyr, returned_mass, remnant_mass, energy, &
         momentum, ejected_mass, net_yield, interpolation_ierr)
    if (interpolation_ierr /= interpolation_ok) then
       ierr = provider_err_interpolation
       return
    end if

    state%ejected_mass = ejected_mass
    state%net_yield = net_yield
    state%returned_mass = returned_mass
    state%remnant_mass = remnant_mass
    state%living_mass = initial_mass - returned_mass - remnant_mass
    state%energy = energy
    state%momentum = momentum

    state%channel_returned_mass(channel_id) = returned_mass
    state%channel_energy(channel_id) = energy
    state%channel_momentum(channel_id,:) = momentum
    state%channel_ejected_mass(channel_id,:) = ejected_mass
    state%channel_net_yield(channel_id,:) = net_yield
    if (.not. ieee_is_finite(state%living_mass)) then
       call clear_cumulative(state)
       ierr = provider_err_nonfinite
    end if
  end subroutine evaluate_channel_cumulative

end module stellar_yield_provider

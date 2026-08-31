! Phase 0 contract for Fortran-native stellar mass return and enrichment.
!
! This module is deliberately independent of RAMSES internals.  It defines
! the common element ordering and the conserved quantities exchanged between
! the yield engine and the AMR deposition layer.

module stellar_enrichment_contract
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, elem_h, elem_he, elem_c, elem_n, elem_o, elem_ne, &
       elem_mg, elem_si, elem_s, elem_ca, elem_fe, channel_wind, channel_agb, &
       channel_snii, channel_snia, channel_pisn
  implicit none

  ! Metadata required to evaluate a single-age stellar population.
  type :: stellar_population_t
     real(stellar_dp) :: formation_time
     real(stellar_dp) :: initial_mass
     real(stellar_dp) :: current_mass
     real(stellar_dp) :: birth_metallicity
     real(stellar_dp) :: birth_mass_fraction(n_stellar_elements)
     integer :: imf_id
     integer :: population_id
     logical :: pisn_enabled
  end type stellar_population_t

  ! Cumulative state at one stellar population age.
  !
  ! ejected_mass is actual material returned to the gas.  net_yield is the
  ! newly produced/consumed component and may be negative.  The latter is not
  ! a gas-mass source term.
  type :: stellar_cumulative_t
     real(stellar_dp) :: ejected_mass(n_stellar_elements)
     real(stellar_dp) :: net_yield(n_stellar_elements)
     real(stellar_dp) :: returned_mass
     real(stellar_dp) :: remnant_mass
     real(stellar_dp) :: living_mass
     real(stellar_dp) :: energy
     real(stellar_dp) :: momentum(3)
     real(stellar_dp) :: channel_returned_mass(n_stellar_channels)
     real(stellar_dp) :: channel_energy(n_stellar_channels)
     real(stellar_dp) :: channel_momentum(n_stellar_channels,3)
     real(stellar_dp) :: channel_ejected_mass(n_stellar_channels, &
                                               n_stellar_elements)
     real(stellar_dp) :: channel_net_yield(n_stellar_channels, &
                                            n_stellar_elements)
  end type stellar_cumulative_t

  ! Increment to be deposited during one hydrodynamic timestep.
  type :: stellar_source_t
     real(stellar_dp) :: ejected_mass(n_stellar_elements)
     real(stellar_dp) :: net_yield(n_stellar_elements)
     real(stellar_dp) :: returned_mass
     real(stellar_dp) :: energy
     real(stellar_dp) :: momentum(3)
     real(stellar_dp) :: channel_returned_mass(n_stellar_channels)
     real(stellar_dp) :: channel_energy(n_stellar_channels)
     real(stellar_dp) :: channel_momentum(n_stellar_channels,3)
     real(stellar_dp) :: channel_ejected_mass(n_stellar_channels, &
                                               n_stellar_elements)
     real(stellar_dp) :: channel_net_yield(n_stellar_channels, &
                                            n_stellar_elements)
  end type stellar_source_t

contains

  subroutine clear_cumulative(state)
    type(stellar_cumulative_t), intent(out) :: state

    state%ejected_mass = 0.0_stellar_dp
    state%net_yield = 0.0_stellar_dp
    state%returned_mass = 0.0_stellar_dp
    state%remnant_mass = 0.0_stellar_dp
    state%living_mass = 0.0_stellar_dp
    state%energy = 0.0_stellar_dp
    state%momentum = 0.0_stellar_dp
    state%channel_returned_mass = 0.0_stellar_dp
    state%channel_energy = 0.0_stellar_dp
    state%channel_momentum = 0.0_stellar_dp
    state%channel_ejected_mass = 0.0_stellar_dp
    state%channel_net_yield = 0.0_stellar_dp
  end subroutine clear_cumulative

  subroutine clear_source(source)
    type(stellar_source_t), intent(out) :: source

    source%ejected_mass = 0.0_stellar_dp
    source%net_yield = 0.0_stellar_dp
    source%returned_mass = 0.0_stellar_dp
    source%energy = 0.0_stellar_dp
    source%momentum = 0.0_stellar_dp
    source%channel_returned_mass = 0.0_stellar_dp
    source%channel_energy = 0.0_stellar_dp
    source%channel_momentum = 0.0_stellar_dp
    source%channel_ejected_mass = 0.0_stellar_dp
    source%channel_net_yield = 0.0_stellar_dp
  end subroutine clear_source

  subroutine cumulative_difference(later, earlier, source)
    type(stellar_cumulative_t), intent(in) :: later
    type(stellar_cumulative_t), intent(in) :: earlier
    type(stellar_source_t), intent(out) :: source

    source%ejected_mass = later%ejected_mass - earlier%ejected_mass
    source%net_yield = later%net_yield - earlier%net_yield
    source%returned_mass = later%returned_mass - earlier%returned_mass
    source%energy = later%energy - earlier%energy
    source%momentum = later%momentum - earlier%momentum
    source%channel_returned_mass = later%channel_returned_mass - &
                                   earlier%channel_returned_mass
    source%channel_energy = later%channel_energy - earlier%channel_energy
    source%channel_momentum = later%channel_momentum - &
                              earlier%channel_momentum
    source%channel_ejected_mass = later%channel_ejected_mass - &
                                  earlier%channel_ejected_mass
    source%channel_net_yield = later%channel_net_yield - &
                               earlier%channel_net_yield
  end subroutine cumulative_difference

end module stellar_enrichment_contract

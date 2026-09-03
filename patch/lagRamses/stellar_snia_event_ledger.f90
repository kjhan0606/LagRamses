! Review-only SNIa event-to-source ledger.
!
! This module converts an already interval-integrated expected event count and
! a per-event source record into the common stellar source contract.  It does
! not select a DTD, binary population, event yield, or physical normalization.
! The production stellar driver remains fail-closed for channel 4 until the
! corresponding approval sidecar is present.

module stellar_snia_event_ledger
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       channel_snia
  use stellar_enrichment_contract, only: stellar_source_t, clear_source
  implicit none

  private

  integer, parameter, public :: snia_event_ledger_ok = 0
  integer, parameter, public :: snia_event_ledger_err_argument = 1
  integer, parameter, public :: snia_event_ledger_err_nonfinite = 2
  integer, parameter, public :: snia_event_ledger_err_mass = 4

  public :: build_snia_event_source

contains

  subroutine build_snia_event_source(expected_events, returned_mass_per_event, &
       ejecta_per_event, net_yield_per_event, energy_per_event, &
       momentum_per_event, source, ierr)
    real(stellar_dp), intent(in) :: expected_events
    real(stellar_dp), intent(in) :: returned_mass_per_event
    real(stellar_dp), intent(in) :: ejecta_per_event(n_stellar_elements)
    real(stellar_dp), intent(in) :: net_yield_per_event(n_stellar_elements)
    real(stellar_dp), intent(in) :: energy_per_event
    real(stellar_dp), intent(in) :: momentum_per_event(3)
    type(stellar_source_t), intent(out) :: source
    integer, intent(out) :: ierr

    real(stellar_dp) :: tracked_ejecta, scale

    call clear_source(source)
    ierr = snia_event_ledger_ok
    if (.not. ieee_is_finite(expected_events) .or. &
         .not. ieee_is_finite(returned_mass_per_event) .or. &
         .not. all(ieee_is_finite(ejecta_per_event)) .or. &
         .not. all(ieee_is_finite(net_yield_per_event)) .or. &
         .not. ieee_is_finite(energy_per_event) .or. &
         .not. all(ieee_is_finite(momentum_per_event)) .or. &
         expected_events < 0.0_stellar_dp .or. &
         returned_mass_per_event < 0.0_stellar_dp .or. &
         minval(ejecta_per_event) < 0.0_stellar_dp .or. &
         energy_per_event < 0.0_stellar_dp) then
       ierr = snia_event_ledger_err_argument
       return
    end if

    tracked_ejecta = sum(ejecta_per_event)
    scale = max(1.0_stellar_dp, returned_mass_per_event, tracked_ejecta)
    if (.not. ieee_is_finite(tracked_ejecta) .or. &
         tracked_ejecta > returned_mass_per_event + 1.0e-12_stellar_dp * scale) then
       ierr = snia_event_ledger_err_mass
       return
    end if

    source%ejected_mass = expected_events * ejecta_per_event
    source%net_yield = expected_events * net_yield_per_event
    source%returned_mass = expected_events * returned_mass_per_event
    source%energy = expected_events * energy_per_event
    source%momentum = expected_events * momentum_per_event
    source%channel_returned_mass(channel_snia) = source%returned_mass
    source%channel_energy(channel_snia) = source%energy
    source%channel_momentum(channel_snia,:) = source%momentum
    source%channel_ejected_mass(channel_snia,:) = source%ejected_mass
    source%channel_net_yield(channel_snia,:) = source%net_yield

    if (.not. all(ieee_is_finite(source%ejected_mass)) .or. &
         .not. all(ieee_is_finite(source%net_yield)) .or. &
         .not. ieee_is_finite(source%returned_mass) .or. &
         .not. ieee_is_finite(source%energy) .or. &
         .not. all(ieee_is_finite(source%momentum)) .or. &
         source%returned_mass < 0.0_stellar_dp .or. &
         source%energy < 0.0_stellar_dp) then
       call clear_source(source)
       ierr = snia_event_ledger_err_nonfinite
    end if
  end subroutine build_snia_event_source

end module stellar_snia_event_ledger

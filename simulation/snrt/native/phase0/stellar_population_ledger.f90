! Population-level mass ledger for channel-resolved stellar returns.
!
! Channel rows may share a progenitor mass range when they represent different
! evolutionary phases.  Their returned material is additive, but the terminal
! remnant is owned by exactly the designated terminal channel.  Living mass is
! therefore computed once after aggregation; it must not be copied from one
! channel's local provider state.

module stellar_population_ledger
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_cumulative_t
  implicit none

  private

  integer, parameter, public :: population_ledger_ok = 0
  integer, parameter, public :: population_ledger_err_argument = 1
  integer, parameter, public :: population_ledger_err_nonfinite = 2
  integer, parameter, public :: population_ledger_err_mass = 4
  integer, parameter, public :: population_ledger_err_owner = 8

  type, public :: stellar_population_ledger_t
     real(stellar_dp) :: initial_mass = 0.0_stellar_dp
     real(stellar_dp) :: returned_mass = 0.0_stellar_dp
     real(stellar_dp) :: tracked_ejecta_mass = 0.0_stellar_dp
     real(stellar_dp) :: untracked_ejecta_mass = 0.0_stellar_dp
     real(stellar_dp) :: remnant_mass = 0.0_stellar_dp
     real(stellar_dp) :: living_mass = 0.0_stellar_dp
     real(stellar_dp) :: channel_returned_mass(n_stellar_channels) = 0.0_stellar_dp
     real(stellar_dp) :: channel_tracked_ejecta_mass(n_stellar_channels) = &
          0.0_stellar_dp
     real(stellar_dp) :: channel_untracked_ejecta_mass(n_stellar_channels) = &
          0.0_stellar_dp
     real(stellar_dp) :: channel_remnant_mass(n_stellar_channels) = 0.0_stellar_dp
  end type stellar_population_ledger_t

  public :: clear_population_ledger
  public :: finalize_population_ledger

contains

  subroutine clear_population_ledger(ledger)
    type(stellar_population_ledger_t), intent(out) :: ledger

    ledger%initial_mass = 0.0_stellar_dp
    ledger%returned_mass = 0.0_stellar_dp
    ledger%tracked_ejecta_mass = 0.0_stellar_dp
    ledger%untracked_ejecta_mass = 0.0_stellar_dp
    ledger%remnant_mass = 0.0_stellar_dp
    ledger%living_mass = 0.0_stellar_dp
    ledger%channel_returned_mass = 0.0_stellar_dp
    ledger%channel_tracked_ejecta_mass = 0.0_stellar_dp
    ledger%channel_untracked_ejecta_mass = 0.0_stellar_dp
    ledger%channel_remnant_mass = 0.0_stellar_dp
  end subroutine clear_population_ledger

  subroutine finalize_population_ledger(population, channel_states, &
       channel_enabled, terminal_remnant_owner, tolerance, ledger, ierr)
    type(stellar_population_t), intent(in) :: population
    type(stellar_cumulative_t), intent(in) :: channel_states(n_stellar_channels)
    logical, intent(in) :: channel_enabled(n_stellar_channels)
    logical, intent(in) :: terminal_remnant_owner(n_stellar_channels)
    real(stellar_dp), intent(in) :: tolerance
    type(stellar_population_ledger_t), intent(out) :: ledger
    integer, intent(out) :: ierr

    real(stellar_dp) :: tol, ejecta_sum, channel_ejecta_sum, scale
    integer :: channel

    call clear_population_ledger(ledger)
    ierr = population_ledger_ok
    if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_stellar_dp) then
       ierr = population_ledger_err_argument
       return
    end if
    tol = max(tolerance, 1.0e-12_stellar_dp)
    ledger%initial_mass = population%initial_mass
    if (.not. ieee_is_finite(ledger%initial_mass) .or. &
         ledger%initial_mass < 0.0_stellar_dp) then
       ierr = population_ledger_err_argument
       return
    end if

    do channel = 1, n_stellar_channels
       if (.not. channel_enabled(channel)) cycle
       if (.not. cumulative_is_finite(channel_states(channel))) then
          ierr = ior(ierr, population_ledger_err_nonfinite)
          cycle
       end if
       if (channel_states(channel)%returned_mass < -tol .or. &
            channel_states(channel)%remnant_mass < -tol .or. &
            minval(channel_states(channel)%ejected_mass) < -tol) then
          ierr = ior(ierr, population_ledger_err_mass)
          cycle
       end if
       ejecta_sum = sum(channel_states(channel)%ejected_mass)
       channel_ejecta_sum = sum(channel_states(channel)% &
            channel_ejected_mass(channel,:))
       scale = max(1.0_stellar_dp, abs(ejecta_sum), &
            abs(channel_ejecta_sum), &
            abs(channel_states(channel)%returned_mass), &
            abs(channel_states(channel)%channel_returned_mass(channel)))
       if (ejecta_sum > channel_states(channel)%returned_mass + tol * scale .or. &
            channel_ejecta_sum > &
            channel_states(channel)%channel_returned_mass(channel) + &
            tol * scale) then
          ierr = ior(ierr, population_ledger_err_mass)
       end if
       if (abs(ejecta_sum - channel_ejecta_sum) > tol * scale .or. &
            abs(channel_states(channel)%returned_mass - &
            channel_states(channel)%channel_returned_mass(channel)) > tol * scale) then
          ierr = ior(ierr, population_ledger_err_mass)
       end if
       ledger%channel_returned_mass(channel) = &
            channel_states(channel)%channel_returned_mass(channel)
       ledger%channel_tracked_ejecta_mass(channel) = max(0.0_stellar_dp, &
            channel_ejecta_sum)
       ledger%channel_untracked_ejecta_mass(channel) = max(0.0_stellar_dp, &
            channel_states(channel)%channel_returned_mass(channel) - &
            channel_ejecta_sum)
       if (terminal_remnant_owner(channel)) then
          ledger%channel_remnant_mass(channel) = &
               channel_states(channel)%remnant_mass
       else if (channel_states(channel)%remnant_mass > tol) then
          ! A non-terminal channel may not smuggle a second remnant into the
          ! population ledger.  The source table must assign ownership first.
          ierr = ior(ierr, population_ledger_err_owner)
       end if
    end do

    ledger%returned_mass = sum(ledger%channel_returned_mass)
    ledger%tracked_ejecta_mass = sum(ledger%channel_tracked_ejecta_mass)
    ledger%untracked_ejecta_mass = sum(ledger%channel_untracked_ejecta_mass)
    ledger%remnant_mass = sum(ledger%channel_remnant_mass)
    ledger%living_mass = ledger%initial_mass - ledger%returned_mass - &
         ledger%remnant_mass
    if (.not. ieee_is_finite(ledger%returned_mass) .or. &
         .not. ieee_is_finite(ledger%remnant_mass) .or. &
         .not. ieee_is_finite(ledger%living_mass)) then
       ierr = ior(ierr, population_ledger_err_nonfinite)
       return
    end if
    scale = max(1.0_stellar_dp, ledger%initial_mass)
    if (ledger%returned_mass < -tol .or. ledger%remnant_mass < -tol .or. &
         ledger%tracked_ejecta_mass < -tol .or. &
         ledger%untracked_ejecta_mass < -tol .or. &
         ledger%living_mass < -tol .or. &
         abs(ledger%returned_mass - ledger%tracked_ejecta_mass - &
         ledger%untracked_ejecta_mass) > tol * scale .or. &
         abs(ledger%initial_mass - ledger%living_mass - &
         ledger%remnant_mass - ledger%returned_mass) > tol * scale) then
       ierr = ior(ierr, population_ledger_err_mass)
    end if
  end subroutine finalize_population_ledger

  logical function cumulative_is_finite(state)
    type(stellar_cumulative_t), intent(in) :: state
    integer :: i, j

    cumulative_is_finite = ieee_is_finite(state%returned_mass) .and. &
         ieee_is_finite(state%remnant_mass) .and. &
         ieee_is_finite(state%living_mass) .and. &
         ieee_is_finite(state%energy)
    do i = 1, 3
       cumulative_is_finite = cumulative_is_finite .and. &
            ieee_is_finite(state%momentum(i))
    end do
    do i = 1, n_stellar_elements
       cumulative_is_finite = cumulative_is_finite .and. &
            ieee_is_finite(state%ejected_mass(i)) .and. &
            ieee_is_finite(state%net_yield(i))
    end do
    do i = 1, size(state%channel_returned_mass)
       cumulative_is_finite = cumulative_is_finite .and. &
            ieee_is_finite(state%channel_returned_mass(i)) .and. &
            ieee_is_finite(state%channel_energy(i))
       do j = 1, 3
          cumulative_is_finite = cumulative_is_finite .and. &
               ieee_is_finite(state%channel_momentum(i,j))
       end do
       do j = 1, n_stellar_elements
          cumulative_is_finite = cumulative_is_finite .and. &
               ieee_is_finite(state%channel_ejected_mass(i,j)) .and. &
               ieee_is_finite(state%channel_net_yield(i,j))
       end do
    end do
  end function cumulative_is_finite

end module stellar_population_ledger

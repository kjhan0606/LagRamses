! Review-only physical contract for a Type-Ia event source.
!
! This module does not choose a progenitor model or populate any physical
! number.  It defines the interface that a future approved source must obey:
! an event consumes an explicitly supplied white-dwarf reservoir, and its
! returned mass, terminal remnant, energy, and signed momentum are checked
! before a source can be built.  The runtime does not call this module while
! SNIa approval is absent.

module stellar_snia_physical_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp
  implicit none

  private

  integer, parameter, public :: snia_contract_ok = 0
  integer, parameter, public :: snia_contract_err_argument = 1
  integer, parameter, public :: snia_contract_err_unapproved = 2
  integer, parameter, public :: snia_contract_err_reservoir = 4
  integer, parameter, public :: snia_contract_err_mass = 8
  integer, parameter, public :: snia_contract_err_momentum = 16

  integer, parameter, public :: snia_wd_debit_per_event = 1
  integer, parameter, public :: snia_momentum_source_frame_vector = 1
  integer, parameter, public :: snia_momentum_isotropic_zero_vector = 2

  type, public :: snia_physical_contract_t
     logical :: approved = .false.
     integer :: wd_debit_policy = 0
     integer :: momentum_policy = 0
     real(stellar_dp) :: returned_mass_per_event = -1.0_stellar_dp
     real(stellar_dp) :: terminal_remnant_per_event = -1.0_stellar_dp
     real(stellar_dp) :: wd_debit_per_event = -1.0_stellar_dp
     real(stellar_dp) :: energy_per_event = -1.0_stellar_dp
     real(stellar_dp) :: momentum_per_event(3) = 0.0_stellar_dp
  end type snia_physical_contract_t

  type, public :: snia_event_budget_t
     real(stellar_dp) :: wd_reservoir_debit = 0.0_stellar_dp
     real(stellar_dp) :: returned_mass = 0.0_stellar_dp
     real(stellar_dp) :: terminal_remnant_mass = 0.0_stellar_dp
     real(stellar_dp) :: energy = 0.0_stellar_dp
     real(stellar_dp) :: momentum(3) = 0.0_stellar_dp
  end type snia_event_budget_t

  public :: validate_snia_physical_contract
  public :: build_snia_event_budget

contains

  subroutine validate_snia_physical_contract(contract, ierr)
    type(snia_physical_contract_t), intent(in) :: contract
    integer, intent(out) :: ierr

    ierr = snia_contract_ok
    if (.not. contract%approved) then
       ierr = snia_contract_err_unapproved
       return
    end if
    if (contract%wd_debit_policy /= snia_wd_debit_per_event .or. &
         contract%momentum_policy /= snia_momentum_source_frame_vector .and. &
         contract%momentum_policy /= snia_momentum_isotropic_zero_vector) then
       ierr = snia_contract_err_argument
       return
    end if
    if (.not. ieee_is_finite(contract%returned_mass_per_event) .or. &
         .not. ieee_is_finite(contract%terminal_remnant_per_event) .or. &
         .not. ieee_is_finite(contract%wd_debit_per_event) .or. &
         .not. ieee_is_finite(contract%energy_per_event) .or. &
         .not. all(ieee_is_finite(contract%momentum_per_event)) .or. &
         contract%returned_mass_per_event < 0.0_stellar_dp .or. &
         contract%terminal_remnant_per_event < 0.0_stellar_dp .or. &
         contract%wd_debit_per_event < 0.0_stellar_dp .or. &
         contract%energy_per_event < 0.0_stellar_dp) then
       ierr = snia_contract_err_argument
       return
    end if
    if (contract%returned_mass_per_event + &
         contract%terminal_remnant_per_event > &
         contract%wd_debit_per_event + 1.0e-12_stellar_dp * &
         max(1.0_stellar_dp, contract%wd_debit_per_event)) then
       ierr = snia_contract_err_mass
       return
    end if
    if (contract%momentum_policy == snia_momentum_isotropic_zero_vector .and. &
         maxval(abs(contract%momentum_per_event)) > &
         1.0e-12_stellar_dp * max(1.0_stellar_dp, &
         maxval(abs(contract%momentum_per_event)))) then
       ierr = snia_contract_err_momentum
    end if
  end subroutine validate_snia_physical_contract

  subroutine build_snia_event_budget(contract, expected_events, &
       wd_reservoir_available, budget, ierr)
    type(snia_physical_contract_t), intent(in) :: contract
    real(stellar_dp), intent(in) :: expected_events
    real(stellar_dp), intent(in) :: wd_reservoir_available
    type(snia_event_budget_t), intent(out) :: budget
    integer, intent(out) :: ierr

    real(stellar_dp) :: scale
    integer :: contract_ierr

    budget%wd_reservoir_debit = 0.0_stellar_dp
    budget%returned_mass = 0.0_stellar_dp
    budget%terminal_remnant_mass = 0.0_stellar_dp
    budget%energy = 0.0_stellar_dp
    budget%momentum = 0.0_stellar_dp
    ierr = snia_contract_ok
    if (.not. ieee_is_finite(expected_events) .or. &
         .not. ieee_is_finite(wd_reservoir_available) .or. &
         expected_events < 0.0_stellar_dp .or. &
         wd_reservoir_available < 0.0_stellar_dp) then
       ierr = snia_contract_err_argument
       return
    end if

    call validate_snia_physical_contract(contract, contract_ierr)
    if (contract_ierr /= snia_contract_ok) then
       ierr = contract_ierr
       return
    end if

    budget%wd_reservoir_debit = expected_events * contract%wd_debit_per_event
    budget%returned_mass = expected_events * contract%returned_mass_per_event
    budget%terminal_remnant_mass = expected_events * &
         contract%terminal_remnant_per_event
    budget%energy = expected_events * contract%energy_per_event
    budget%momentum = expected_events * contract%momentum_per_event
    scale = max(1.0_stellar_dp, wd_reservoir_available, &
         budget%wd_reservoir_debit)
    if (.not. all(ieee_is_finite((/budget%wd_reservoir_debit, &
         budget%returned_mass, budget%terminal_remnant_mass, budget%energy/))) .or. &
         .not. all(ieee_is_finite(budget%momentum))) then
       budget%wd_reservoir_debit = 0.0_stellar_dp
       budget%returned_mass = 0.0_stellar_dp
       budget%terminal_remnant_mass = 0.0_stellar_dp
       budget%energy = 0.0_stellar_dp
       budget%momentum = 0.0_stellar_dp
       ierr = snia_contract_err_argument
       return
    end if
    if (budget%wd_reservoir_debit > wd_reservoir_available + &
         1.0e-12_stellar_dp * scale) then
       budget%wd_reservoir_debit = 0.0_stellar_dp
       budget%returned_mass = 0.0_stellar_dp
       budget%terminal_remnant_mass = 0.0_stellar_dp
       budget%energy = 0.0_stellar_dp
       budget%momentum = 0.0_stellar_dp
       ierr = snia_contract_err_reservoir
    end if
  end subroutine build_snia_event_budget

end module stellar_snia_physical_contract

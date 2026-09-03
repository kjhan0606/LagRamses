! Physical contract for the approved baseline Type-Ia event source.
!
! The selected HESMA n100 source is populated through the versioned F-P2
! sidecar; an event consumes an explicitly supplied white-dwarf reservoir,
! and its returned mass, terminal remnant, energy, and signed momentum are
! checked before a source can be built.  Runtime activation remains a separate
! gate until the AMR/MPI caller is connected.

module stellar_snia_physical_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  implicit none

  private

  integer, parameter, public :: snia_contract_ok = 0
  integer, parameter, public :: snia_contract_err_argument = 1
  integer, parameter, public :: snia_contract_err_unapproved = 2
  integer, parameter, public :: snia_contract_err_reservoir = 4
  integer, parameter, public :: snia_contract_err_mass = 8
  integer, parameter, public :: snia_contract_err_momentum = 16
  integer, parameter, public :: snia_contract_err_yield = 32

  integer, parameter, public :: snia_wd_debit_per_event = 1
  integer, parameter, public :: snia_shortfall_reject_event_interval = 1
  integer, parameter, public :: snia_momentum_source_frame_vector = 1
  integer, parameter, public :: snia_momentum_isotropic_zero_vector = 2
  integer, parameter, public :: snia_momentum_radial_magnitude = 3

  type, public :: snia_physical_contract_t
     logical :: approved = .false.
     integer :: wd_debit_policy = 0
     integer :: shortfall_policy = 0
     integer :: momentum_policy = 0
     real(stellar_dp) :: returned_mass_per_event = -1.0_stellar_dp
     real(stellar_dp) :: terminal_remnant_per_event = -1.0_stellar_dp
     real(stellar_dp) :: wd_debit_per_event = -1.0_stellar_dp
     real(stellar_dp) :: energy_per_event = -1.0_stellar_dp
     real(stellar_dp) :: momentum_per_event(3) = 0.0_stellar_dp
     real(stellar_dp) :: radial_momentum_per_event = -1.0_stellar_dp
     ! Tracked ejecta are actual returned material, not net yields.  Any
     ! untracked returned mass remains the explicit residual in returned_mass.
     real(stellar_dp) :: ejected_mass_per_event(n_stellar_elements) = 0.0_stellar_dp
     real(stellar_dp) :: net_yield_per_event(n_stellar_elements) = 0.0_stellar_dp
     character(len=128) :: yield_source_id = ''
     character(len=128) :: yield_source_sha256 = ''
     character(len=128) :: source_commit_binding = ''
     character(len=128) :: conversion_code_sha256 = ''
     character(len=128) :: approval_id = ''
  end type snia_physical_contract_t

  type, public :: snia_event_budget_t
     real(stellar_dp) :: wd_reservoir_debit = 0.0_stellar_dp
     real(stellar_dp) :: returned_mass = 0.0_stellar_dp
     real(stellar_dp) :: terminal_remnant_mass = 0.0_stellar_dp
     real(stellar_dp) :: energy = 0.0_stellar_dp
     real(stellar_dp) :: momentum(3) = 0.0_stellar_dp
     real(stellar_dp) :: ejected_mass(n_stellar_elements) = 0.0_stellar_dp
     real(stellar_dp) :: net_yield(n_stellar_elements) = 0.0_stellar_dp
  end type snia_event_budget_t

  public :: validate_snia_physical_contract
  public :: resolve_snia_event_momentum
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
         contract%shortfall_policy /= snia_shortfall_reject_event_interval .or. &
         (contract%momentum_policy /= snia_momentum_source_frame_vector .and. &
         contract%momentum_policy /= snia_momentum_isotropic_zero_vector .and. &
         contract%momentum_policy /= snia_momentum_radial_magnitude)) then
       ierr = snia_contract_err_argument
       return
    end if
    if (.not. ieee_is_finite(contract%returned_mass_per_event) .or. &
         .not. ieee_is_finite(contract%terminal_remnant_per_event) .or. &
         .not. ieee_is_finite(contract%wd_debit_per_event) .or. &
         .not. ieee_is_finite(contract%energy_per_event) .or. &
         .not. ieee_is_finite(contract%radial_momentum_per_event) .or. &
         .not. all(ieee_is_finite(contract%momentum_per_event)) .or. &
         .not. all(ieee_is_finite(contract%ejected_mass_per_event)) .or. &
         .not. all(ieee_is_finite(contract%net_yield_per_event)) .or. &
         contract%returned_mass_per_event < 0.0_stellar_dp .or. &
         contract%terminal_remnant_per_event < 0.0_stellar_dp .or. &
         contract%wd_debit_per_event < 0.0_stellar_dp .or. &
         contract%energy_per_event < 0.0_stellar_dp .or. &
         minval(contract%ejected_mass_per_event) < 0.0_stellar_dp .or. &
         (contract%momentum_policy == snia_momentum_radial_magnitude .and. &
         contract%radial_momentum_per_event < 0.0_stellar_dp)) then
       ierr = snia_contract_err_argument
       return
    end if
    if (abs(contract%returned_mass_per_event + &
         contract%terminal_remnant_per_event - &
         contract%wd_debit_per_event) > 1.0e-12_stellar_dp * &
         max(1.0_stellar_dp, contract%wd_debit_per_event, &
         contract%returned_mass_per_event, contract%terminal_remnant_per_event)) then
       ierr = snia_contract_err_mass
       return
    end if
    if (sum(contract%ejected_mass_per_event) > &
         contract%returned_mass_per_event + 1.0e-12_stellar_dp * &
         max(1.0_stellar_dp, contract%returned_mass_per_event)) then
       ierr = snia_contract_err_yield
       return
    end if
    if (len_trim(contract%yield_source_id) == 0 .or. &
         .not. is_hex_digest(contract%yield_source_sha256, 64) .or. &
         .not. is_hex_digest(contract%source_commit_binding, 40) .or. &
         .not. is_hex_digest(contract%conversion_code_sha256, 64) .or. &
         len_trim(contract%approval_id) == 0) then
       ierr = snia_contract_err_argument
       return
    end if
    if (contract%momentum_policy == snia_momentum_isotropic_zero_vector .and. &
         maxval(abs(contract%momentum_per_event)) > &
         1.0e-12_stellar_dp * max(1.0_stellar_dp, &
         maxval(abs(contract%momentum_per_event)))) then
       ierr = snia_contract_err_momentum
    end if
  end subroutine validate_snia_physical_contract

  subroutine resolve_snia_event_momentum(contract, direction, momentum, ierr)
    type(snia_physical_contract_t), intent(in) :: contract
    real(stellar_dp), intent(in), optional :: direction(3)
    real(stellar_dp), intent(out) :: momentum(3)
    integer, intent(out) :: ierr

    real(stellar_dp) :: direction_norm
    integer :: contract_ierr

    momentum = 0.0_stellar_dp
    call validate_snia_physical_contract(contract, contract_ierr)
    if (contract_ierr /= snia_contract_ok) then
       ierr = contract_ierr
       return
    end if
    ierr = snia_contract_ok
    select case (contract%momentum_policy)
    case (snia_momentum_source_frame_vector)
       momentum = contract%momentum_per_event
    case (snia_momentum_isotropic_zero_vector)
       momentum = 0.0_stellar_dp
    case (snia_momentum_radial_magnitude)
       if (.not. present(direction)) then
          ierr = snia_contract_err_momentum
          return
       end if
       if (.not. all(ieee_is_finite(direction))) then
          ierr = snia_contract_err_momentum
          return
       end if
       direction_norm = sqrt(sum(direction**2))
       if (.not. ieee_is_finite(direction_norm) .or. &
            abs(direction_norm - 1.0_stellar_dp) > 1.0e-12_stellar_dp) then
          ierr = snia_contract_err_momentum
          return
       end if
       momentum = contract%radial_momentum_per_event * direction
    case default
       ierr = snia_contract_err_argument
    end select
  end subroutine resolve_snia_event_momentum

    subroutine build_snia_event_budget(contract, expected_events, &
       wd_reservoir_available, budget, ierr, deposition_direction)
    type(snia_physical_contract_t), intent(in) :: contract
    real(stellar_dp), intent(in) :: expected_events
    real(stellar_dp), intent(in) :: wd_reservoir_available
    type(snia_event_budget_t), intent(out) :: budget
    integer, intent(out) :: ierr
    real(stellar_dp), intent(in), optional :: deposition_direction(3)

    real(stellar_dp) :: scale, event_momentum(3)
    integer :: contract_ierr, momentum_ierr

    budget%wd_reservoir_debit = 0.0_stellar_dp
    budget%returned_mass = 0.0_stellar_dp
    budget%terminal_remnant_mass = 0.0_stellar_dp
    budget%energy = 0.0_stellar_dp
    budget%momentum = 0.0_stellar_dp
    budget%ejected_mass = 0.0_stellar_dp
    budget%net_yield = 0.0_stellar_dp
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
    budget%ejected_mass = expected_events * contract%ejected_mass_per_event
    budget%net_yield = expected_events * contract%net_yield_per_event
    call resolve_snia_event_momentum(contract, deposition_direction, &
         event_momentum, momentum_ierr)
    if (momentum_ierr /= snia_contract_ok) then
       budget%wd_reservoir_debit = 0.0_stellar_dp
       budget%returned_mass = 0.0_stellar_dp
       budget%terminal_remnant_mass = 0.0_stellar_dp
       budget%energy = 0.0_stellar_dp
       budget%momentum = 0.0_stellar_dp
       budget%ejected_mass = 0.0_stellar_dp
       budget%net_yield = 0.0_stellar_dp
       ierr = momentum_ierr
       return
    end if
    budget%momentum = expected_events * event_momentum
    scale = max(1.0_stellar_dp, wd_reservoir_available, &
         budget%wd_reservoir_debit)
    if (.not. all(ieee_is_finite((/budget%wd_reservoir_debit, &
         budget%returned_mass, budget%terminal_remnant_mass, budget%energy/))) .or. &
         .not. all(ieee_is_finite(budget%momentum)) .or. &
         .not. all(ieee_is_finite(budget%ejected_mass)) .or. &
         .not. all(ieee_is_finite(budget%net_yield))) then
       budget%wd_reservoir_debit = 0.0_stellar_dp
       budget%returned_mass = 0.0_stellar_dp
       budget%terminal_remnant_mass = 0.0_stellar_dp
       budget%energy = 0.0_stellar_dp
       budget%momentum = 0.0_stellar_dp
       budget%ejected_mass = 0.0_stellar_dp
       budget%net_yield = 0.0_stellar_dp
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
       budget%ejected_mass = 0.0_stellar_dp
       budget%net_yield = 0.0_stellar_dp
       ierr = snia_contract_err_reservoir
    end if
  end subroutine build_snia_event_budget

  logical function is_hex_digest(value, required_length)
    character(len=*), intent(in) :: value
    integer, intent(in) :: required_length
    integer :: index, code

    is_hex_digest = len_trim(value) == required_length
    if (.not. is_hex_digest) return
    do index = 1, required_length
       code = iachar(value(index:index))
       if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
            (code >= iachar('a') .and. code <= iachar('f')) .or. &
            (code >= iachar('A') .and. code <= iachar('F')))) then
          is_hex_digest = .false.
          return
       end if
    end do
  end function is_hex_digest

end module stellar_snia_physical_contract

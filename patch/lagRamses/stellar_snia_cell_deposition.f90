! Review-only SNIa event-to-cell deposition adapter.
!
! The event budget is expressed in physical units: returned mass in Msun,
! momentum in g cm s^-1, and event energy in erg.  This adapter converts one
! already-selected target cell into conserved density increments.  It does
! not select AMR neighbours, perform MPI exchange, update a RAMSES array, or
! activate SNIa at runtime.  The source momentum is assumed to have already
! been resolved in the target-cell frame by stellar_snia_physical_contract.
!
! The only admitted energy policy in this review-only stage puts the complete
! event energy into the hydro total-energy increment.  A future split between
! thermal and non-thermal reservoirs must add an explicit receiver and
! conservation test before a fractional policy can be admitted.

module stellar_snia_cell_deposition
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_native_units, only: solar_mass_cgs
  use stellar_snia_physical_contract, only: snia_event_budget_t
  use stellar_enrichment_contract, only: generic_metal_ejecta_mass
  implicit none

  private

  integer, parameter, public :: snia_deposition_ok = 0
  integer, parameter, public :: snia_deposition_err_argument = 1
  integer, parameter, public :: snia_deposition_err_policy = 2
  integer, parameter, public :: snia_deposition_err_budget = 4
  integer, parameter, public :: snia_deposition_err_result = 8

  integer, parameter, public :: snia_thermal_all_to_total_energy = 1

  type, public :: snia_thermal_coupling_t
     logical :: approved = .false.
     integer :: mode = 0
     real(stellar_dp) :: thermal_fraction = -1.0_stellar_dp
     ! If event momentum is nonzero, its kinetic energy and bulk cross-term
     ! must be included explicitly.  A source model whose energy already
     ! includes that kinetic term must use zero-net momentum instead.
     logical :: include_event_momentum_kinetic = .false.
  end type snia_thermal_coupling_t

  type, public :: snia_cell_increment_t
     real(stellar_dp) :: mass_density = 0.0_stellar_dp
     real(stellar_dp) :: momentum_density(3) = 0.0_stellar_dp
     real(stellar_dp) :: event_energy_density = 0.0_stellar_dp
     real(stellar_dp) :: bulk_kinetic_energy_density = 0.0_stellar_dp
     real(stellar_dp) :: total_energy_density = 0.0_stellar_dp
     real(stellar_dp) :: element_mass_density(n_stellar_elements) = 0.0_stellar_dp
     real(stellar_dp) :: total_metal_density = 0.0_stellar_dp
  end type snia_cell_increment_t

  public :: validate_snia_thermal_coupling
  public :: build_snia_cell_increment

contains

  subroutine validate_snia_thermal_coupling(coupling, ierr)
    type(snia_thermal_coupling_t), intent(in) :: coupling
    integer, intent(out) :: ierr

    ierr = snia_deposition_ok
    if (.not. coupling%approved) then
       ierr = snia_deposition_err_policy
       return
    end if
    if (coupling%mode /= snia_thermal_all_to_total_energy) then
       ierr = snia_deposition_err_policy
       return
    end if
    if (.not. ieee_is_finite(coupling%thermal_fraction) .or. &
         coupling%thermal_fraction < 0.0_stellar_dp .or. &
         coupling%thermal_fraction > 1.0_stellar_dp) then
       ierr = snia_deposition_err_policy
       return
    end if
    if (abs(coupling%thermal_fraction - 1.0_stellar_dp) > &
         1.0e-12_stellar_dp) then
       ! A partial fraction has no explicitly declared non-thermal receiver.
       ierr = snia_deposition_err_policy
    end if
  end subroutine validate_snia_thermal_coupling

  subroutine build_snia_cell_increment(budget, volume_cm3, bulk_velocity_cm_s, &
       coupling, increment, ierr)
    type(snia_event_budget_t), intent(in) :: budget
    real(stellar_dp), intent(in) :: volume_cm3
    real(stellar_dp), intent(in) :: bulk_velocity_cm_s(3)
    type(snia_thermal_coupling_t), intent(in) :: coupling
    type(snia_cell_increment_t), intent(out) :: increment
    integer, intent(out) :: ierr

    real(stellar_dp) :: returned_mass_cgs
    real(stellar_dp) :: bulk_momentum(3)
    real(stellar_dp) :: bulk_kinetic_energy
    integer :: policy_ierr

    increment%mass_density = 0.0_stellar_dp
    increment%momentum_density = 0.0_stellar_dp
    increment%event_energy_density = 0.0_stellar_dp
    increment%bulk_kinetic_energy_density = 0.0_stellar_dp
    increment%total_energy_density = 0.0_stellar_dp
    increment%element_mass_density = 0.0_stellar_dp
    increment%total_metal_density = 0.0_stellar_dp
    ierr = snia_deposition_ok

    call validate_snia_thermal_coupling(coupling, policy_ierr)
    if (policy_ierr /= snia_deposition_ok) then
       ierr = policy_ierr
       return
    end if
    if (.not. ieee_is_finite(volume_cm3) .or. volume_cm3 <= 0.0_stellar_dp .or. &
         .not. all(ieee_is_finite(bulk_velocity_cm_s))) then
       ierr = snia_deposition_err_argument
       return
    end if
    if (.not. all(ieee_is_finite((/&
         budget%wd_reservoir_debit, budget%returned_mass, &
         budget%terminal_remnant_mass, budget%energy/))) .or. &
         .not. all(ieee_is_finite(budget%momentum)) .or. &
         .not. all(ieee_is_finite(budget%ejected_mass)) .or. &
         .not. all(ieee_is_finite(budget%net_yield)) .or. &
         budget%wd_reservoir_debit < 0.0_stellar_dp .or. &
         budget%returned_mass < 0.0_stellar_dp .or. &
         budget%terminal_remnant_mass < 0.0_stellar_dp .or. &
         budget%energy < 0.0_stellar_dp) then
       ierr = snia_deposition_err_budget
       return
    end if
    if (abs(budget%returned_mass + budget%terminal_remnant_mass - &
         budget%wd_reservoir_debit) > 1.0e-12_stellar_dp * &
         max(1.0_stellar_dp, budget%wd_reservoir_debit, &
         budget%returned_mass, budget%terminal_remnant_mass) .or. &
         sum(budget%ejected_mass) > budget%returned_mass + &
         1.0e-12_stellar_dp * max(1.0_stellar_dp, budget%returned_mass)) then
       ierr = snia_deposition_err_budget
       return
    end if
    if (maxval(abs(budget%momentum)) > 0.0_stellar_dp .and. &
         .not. coupling%include_event_momentum_kinetic) then
       ierr = snia_deposition_err_policy
       return
    end if

    returned_mass_cgs = budget%returned_mass * solar_mass_cgs
    bulk_momentum = returned_mass_cgs * bulk_velocity_cm_s + budget%momentum
    bulk_kinetic_energy = 0.5_stellar_dp * returned_mass_cgs * &
         sum(bulk_velocity_cm_s**2) + sum(bulk_velocity_cm_s * budget%momentum)
    if (returned_mass_cgs > 0.0_stellar_dp) then
       bulk_kinetic_energy = bulk_kinetic_energy + &
            0.5_stellar_dp * sum(budget%momentum**2) / returned_mass_cgs
    else if (maxval(abs(budget%momentum)) > 0.0_stellar_dp) then
       ierr = snia_deposition_err_budget
       return
    end if
    increment%mass_density = returned_mass_cgs / volume_cm3
    increment%momentum_density = bulk_momentum / volume_cm3
    increment%event_energy_density = coupling%thermal_fraction * &
         budget%energy / volume_cm3
    increment%bulk_kinetic_energy_density = bulk_kinetic_energy / volume_cm3
    increment%element_mass_density = budget%ejected_mass * solar_mass_cgs / volume_cm3
    increment%total_metal_density = generic_metal_ejecta_mass( &
         budget%returned_mass, budget%ejected_mass) * solar_mass_cgs / volume_cm3
    increment%total_energy_density = increment%event_energy_density + &
         increment%bulk_kinetic_energy_density

    if (.not. ieee_is_finite(increment%mass_density) .or. &
         .not. all(ieee_is_finite(increment%momentum_density)) .or. &
         .not. ieee_is_finite(increment%event_energy_density) .or. &
         .not. ieee_is_finite(increment%bulk_kinetic_energy_density) .or. &
         .not. ieee_is_finite(increment%total_energy_density) .or. &
         .not. all(ieee_is_finite(increment%element_mass_density)) .or. &
         .not. ieee_is_finite(increment%total_metal_density)) then
       increment%mass_density = 0.0_stellar_dp
       increment%momentum_density = 0.0_stellar_dp
       increment%event_energy_density = 0.0_stellar_dp
       increment%bulk_kinetic_energy_density = 0.0_stellar_dp
       increment%total_energy_density = 0.0_stellar_dp
       increment%element_mass_density = 0.0_stellar_dp
       increment%total_metal_density = 0.0_stellar_dp
       ierr = snia_deposition_err_result
    end if
  end subroutine build_snia_cell_increment

end module stellar_snia_cell_deposition

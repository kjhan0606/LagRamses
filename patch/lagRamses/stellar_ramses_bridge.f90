! Phase 0 RAMSES conserved-variable bridge.
!
! RAMSES-specific variable indices are supplied by the caller.  This keeps
! the stellar source modules independent of a particular lagRamses nvar
! layout while allowing the deposition operation to update uold directly.
! The SNIa bridge below additionally requires the code length, density, and
! velocity scales so physical event increments cannot be written into uold as
! if they were already code units.

module stellar_ramses_bridge
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_enrichment_contract, only: stellar_source_t, &
       generic_metal_ejecta_mass
  use stellar_snia_physical_contract, only: snia_event_budget_t
  use stellar_snia_cell_deposition, only: snia_thermal_coupling_t, &
       snia_cell_increment_t, &
       snia_deposition_ok, build_snia_cell_increment
  implicit none

  private
  integer, parameter, public :: ramses_bridge_ok = 0
  integer, parameter, public :: ramses_bridge_err_argument = 1
  integer, parameter, public :: ramses_bridge_err_source = 2
  integer, parameter, public :: ramses_bridge_err_index = 4
  integer, parameter, public :: ramses_bridge_err_closure = 8
  integer, parameter, public :: ramses_bridge_err_snia = 16
  integer, parameter, public :: ramses_bridge_err_result = 32

  public :: deposit_source_to_uold
  public :: deposit_snia_budget_to_uold

contains

  subroutine deposit_source_to_uold(source, nvar, n_cells, cell_volume, &
       weights, density_var, energy_var, momentum_var, element_var, uold, &
       tolerance, ierr, total_metal_var)
    type(stellar_source_t), intent(in) :: source
    integer, intent(in) :: nvar, n_cells
    real(stellar_dp), intent(in) :: cell_volume(n_cells)
    real(stellar_dp), intent(in) :: weights(n_cells)
    integer, intent(in) :: density_var, energy_var
    integer, intent(in) :: momentum_var(3)
    integer, intent(in) :: element_var(n_stellar_elements)
    real(stellar_dp), intent(inout) :: uold(nvar,n_cells)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    integer, intent(in), optional :: total_metal_var

    real(stellar_dp) :: tol, weight_sum, ejected_sum, scale
    real(stellar_dp) :: normalized_weight, cell_mass, element_mass
    real(stellar_dp) :: generic_metal_mass
    integer :: cell, element

    ierr = ramses_bridge_ok
    tol = max(tolerance, 1.0e-12_stellar_dp)

    if (nvar <= 0 .or. n_cells <= 0 .or. density_var < 1 .or. &
         density_var > nvar .or. energy_var < 1 .or. energy_var > nvar .or. &
         any(momentum_var < 1) .or. any(momentum_var > nvar)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (any(cell_volume <= 0.0_stellar_dp) .or. &
         any(weights < 0.0_stellar_dp)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (any(element_var < 0) .or. any(element_var > nvar)) then
       ierr = ramses_bridge_err_index
       return
    end if
    if (present(total_metal_var)) then
       if (total_metal_var < 0 .or. total_metal_var > nvar) then
          ierr = ramses_bridge_err_index
          return
       end if
    end if

    weight_sum = sum(weights)
    if (weight_sum <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    ejected_sum = sum(source%ejected_mass)
    scale = max(tiny(1.0_stellar_dp), abs(source%returned_mass), &
         abs(ejected_sum))
    if (source%returned_mass < -tol * scale .or. &
         minval(source%ejected_mass) < -tol * scale) then
       ierr = ramses_bridge_err_source
       return
    end if
    if (ejected_sum > source%returned_mass + tol * scale) then
       ierr = ramses_bridge_err_closure
       return
    end if
    generic_metal_mass = generic_metal_ejecta_mass(source%returned_mass, &
         source%ejected_mass)

    ! All validation is completed before uold is modified.
    do cell = 1, n_cells
       normalized_weight = weights(cell) / weight_sum
       cell_mass = normalized_weight * source%returned_mass / cell_volume(cell)
       uold(density_var,cell) = uold(density_var,cell) + cell_mass
       if (present(total_metal_var)) then
          if (total_metal_var > 0) then
             uold(total_metal_var,cell) = uold(total_metal_var,cell) + &
                  normalized_weight * generic_metal_mass / cell_volume(cell)
          end if
       end if

       do element = 1, n_stellar_elements
          if (element_var(element) == 0) cycle
          element_mass = normalized_weight * source%ejected_mass(element) / &
               cell_volume(cell)
          uold(element_var(element),cell) = &
               uold(element_var(element),cell) + element_mass
       end do

       uold(energy_var,cell) = uold(energy_var,cell) + &
            normalized_weight * source%energy / cell_volume(cell)
       uold(momentum_var(1),cell) = uold(momentum_var(1),cell) + &
            normalized_weight * source%momentum(1) / cell_volume(cell)
       uold(momentum_var(2),cell) = uold(momentum_var(2),cell) + &
            normalized_weight * source%momentum(2) / cell_volume(cell)
       uold(momentum_var(3),cell) = uold(momentum_var(3),cell) + &
            normalized_weight * source%momentum(3) / cell_volume(cell)
    end do
  end subroutine deposit_source_to_uold

  subroutine deposit_snia_budget_to_uold(budget, coupling, bulk_velocity_cm_s, &
       scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, &
       n_cells, cell_volume_code, weights, density_var, energy_var, &
       momentum_var, uold, tolerance, ierr)
    type(snia_event_budget_t), intent(in) :: budget
    type(snia_thermal_coupling_t), intent(in) :: coupling
    real(stellar_dp), intent(in) :: bulk_velocity_cm_s(3)
    real(stellar_dp), intent(in) :: scale_length_cgs, scale_density_cgs
    real(stellar_dp), intent(in) :: scale_velocity_cgs
    integer, intent(in) :: nvar, n_cells
    real(stellar_dp), intent(in) :: cell_volume_code(n_cells)
    real(stellar_dp), intent(in) :: weights(n_cells)
    integer, intent(in) :: density_var, energy_var
    integer, intent(in) :: momentum_var(3)
    real(stellar_dp), intent(inout) :: uold(nvar,n_cells)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr

    type(snia_cell_increment_t), allocatable :: increments(:)
    real(stellar_dp) :: weight_sum, normalized_weight
    real(stellar_dp) :: volume_cm3, density_scale, momentum_scale
    real(stellar_dp) :: energy_scale
    integer :: cell, allocation_status, increment_ierr

    ierr = ramses_bridge_ok
    if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_stellar_dp .or. &
         .not. all(ieee_is_finite(bulk_velocity_cm_s)) .or. &
         .not. ieee_is_finite(scale_length_cgs) .or. &
         .not. ieee_is_finite(scale_density_cgs) .or. &
         .not. ieee_is_finite(scale_velocity_cgs) .or. &
         scale_length_cgs <= 0.0_stellar_dp .or. &
         scale_density_cgs <= 0.0_stellar_dp .or. &
         scale_velocity_cgs <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (nvar <= 0 .or. n_cells <= 0 .or. density_var < 1 .or. &
         density_var > nvar .or. energy_var < 1 .or. energy_var > nvar .or. &
         any(momentum_var < 1) .or. any(momentum_var > nvar)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (energy_var == density_var .or. &
         any(momentum_var == density_var) .or. &
         any(momentum_var == energy_var)) then
       ierr = ramses_bridge_err_index
       return
    end if
    if (.not. all(ieee_is_finite(cell_volume_code)) .or. &
         .not. all(ieee_is_finite(weights)) .or. &
         any(cell_volume_code <= 0.0_stellar_dp) .or. &
         any(weights < 0.0_stellar_dp)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    weight_sum = sum(weights)
    if (.not. ieee_is_finite(weight_sum) .or. weight_sum <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if

    allocate(increments(n_cells), stat=allocation_status)
    if (allocation_status /= 0) then
       ierr = ramses_bridge_err_result
       return
    end if

    ! Build every cell increment before modifying uold.  This makes the
    ! bridge transactional when a policy or a budget fails at any cell.
    density_scale = scale_density_cgs
    momentum_scale = scale_density_cgs * scale_velocity_cgs
    energy_scale = scale_density_cgs * scale_velocity_cgs**2
    do cell = 1, n_cells
       volume_cm3 = cell_volume_code(cell) * scale_length_cgs**3
       if (.not. ieee_is_finite(volume_cm3) .or. volume_cm3 <= 0.0_stellar_dp) then
          deallocate(increments)
          ierr = ramses_bridge_err_result
          return
       end if
       call build_snia_cell_increment(budget, volume_cm3, &
            bulk_velocity_cm_s, coupling, increments(cell), increment_ierr)
       if (increment_ierr /= snia_deposition_ok) then
          deallocate(increments)
          ierr = ramses_bridge_err_snia
          return
       end if
    end do

    do cell = 1, n_cells
       if (weights(cell) == 0.0_stellar_dp) cycle
       normalized_weight = weights(cell) / weight_sum
       uold(density_var,cell) = uold(density_var,cell) + &
            normalized_weight * increments(cell)%mass_density / density_scale
       uold(energy_var,cell) = uold(energy_var,cell) + &
            normalized_weight * increments(cell)%total_energy_density / energy_scale
       uold(momentum_var(1),cell) = uold(momentum_var(1),cell) + &
            normalized_weight * increments(cell)%momentum_density(1) / momentum_scale
       uold(momentum_var(2),cell) = uold(momentum_var(2),cell) + &
            normalized_weight * increments(cell)%momentum_density(2) / momentum_scale
       uold(momentum_var(3),cell) = uold(momentum_var(3),cell) + &
            normalized_weight * increments(cell)%momentum_density(3) / momentum_scale
    end do
    deallocate(increments)
  end subroutine deposit_snia_budget_to_uold

end module stellar_ramses_bridge

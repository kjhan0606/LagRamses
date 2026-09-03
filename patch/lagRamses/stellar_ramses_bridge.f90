! Phase 0 RAMSES conserved-variable bridge.
!
! RAMSES-specific variable indices are supplied by the caller.  This keeps
! the stellar source modules independent of a particular lagRamses nvar
! layout while allowing the deposition operation to update uold directly.

module stellar_ramses_bridge
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_enrichment_contract, only: stellar_source_t, &
       generic_metal_ejecta_mass
  implicit none

  private
  integer, parameter, public :: ramses_bridge_ok = 0
  integer, parameter, public :: ramses_bridge_err_argument = 1
  integer, parameter, public :: ramses_bridge_err_source = 2
  integer, parameter, public :: ramses_bridge_err_index = 4
  integer, parameter, public :: ramses_bridge_err_closure = 8

  public :: deposit_source_to_uold

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

end module stellar_ramses_bridge

! Phase 0 local gas-cell deposition adapter.
!
! The AMR/MPI layer supplies the target-cell list and non-negative weights.
! This routine updates conserved densities only; it does not select neighbors,
! perform MPI exchange, or overwrite a metallicity fraction directly.

module stellar_cell_deposition
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_enrichment_contract, only: stellar_source_t
  implicit none

  private
  integer, parameter, public :: deposition_ok = 0
  integer, parameter, public :: deposition_err_argument = 1
  integer, parameter, public :: deposition_err_source = 2
  integer, parameter, public :: deposition_err_closure = 4

  public :: deposit_stellar_source

contains

  subroutine deposit_stellar_source(source, n_cells, cell_volume, weights, &
       gas_density, gas_element_density, gas_energy_density, &
       gas_momentum_density, tolerance, ierr)
    type(stellar_source_t), intent(in) :: source
    integer, intent(in) :: n_cells
    real(stellar_dp), intent(in) :: cell_volume(n_cells)
    real(stellar_dp), intent(in) :: weights(n_cells)
    real(stellar_dp), intent(inout) :: gas_density(n_cells)
    real(stellar_dp), intent(inout) :: &
         gas_element_density(n_stellar_elements,n_cells)
    real(stellar_dp), intent(inout), optional :: gas_energy_density(n_cells)
    real(stellar_dp), intent(inout), optional :: &
         gas_momentum_density(3,n_cells)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr

    real(stellar_dp) :: tol, weight_sum, ejected_sum, scale
    real(stellar_dp) :: normalized_weight, delta_mass, delta_element
    integer :: cell, element

    ierr = deposition_ok
    tol = max(tolerance, 1.0e-12_stellar_dp)

    if (n_cells <= 0) then
       ierr = deposition_err_argument
       return
    end if
    if (any(cell_volume <= 0.0_stellar_dp) .or. &
         any(weights < 0.0_stellar_dp)) then
       ierr = deposition_err_argument
       return
    end if

    weight_sum = sum(weights)
    if (weight_sum <= 0.0_stellar_dp) then
       ierr = deposition_err_argument
       return
    end if

    if (source%returned_mass < -tol .or. &
         minval(source%ejected_mass) < -tol) then
       ierr = deposition_err_source
       return
    end if

    ejected_sum = sum(source%ejected_mass)
    scale = max(1.0_stellar_dp, abs(source%returned_mass), &
         abs(ejected_sum))
    if (abs(ejected_sum - source%returned_mass) > tol * scale) then
       ierr = deposition_err_closure
       return
    end if

    ! All checks occur before the first update, so a rejected source leaves
    ! every gas array unchanged.
    do cell = 1, n_cells
       normalized_weight = weights(cell) / weight_sum
       delta_mass = normalized_weight * source%returned_mass
       gas_density(cell) = gas_density(cell) + &
            delta_mass / cell_volume(cell)

       do element = 1, n_stellar_elements
          delta_element = normalized_weight * source%ejected_mass(element)
          gas_element_density(element,cell) = &
               gas_element_density(element,cell) + &
               delta_element / cell_volume(cell)
       end do

       if (present(gas_energy_density)) then
          gas_energy_density(cell) = gas_energy_density(cell) + &
               normalized_weight * source%energy / cell_volume(cell)
       end if
       if (present(gas_momentum_density)) then
          gas_momentum_density(:,cell) = gas_momentum_density(:,cell) + &
               normalized_weight * source%momentum / cell_volume(cell)
       end if
    end do
  end subroutine deposit_stellar_source

end module stellar_cell_deposition

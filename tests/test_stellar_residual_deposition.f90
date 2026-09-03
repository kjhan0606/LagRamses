program test_stellar_residual_deposition
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       elem_h, elem_he, elem_c, elem_fe
  use stellar_enrichment_contract, only: stellar_source_t, clear_source
  use stellar_cell_deposition, only: deposit_stellar_source, deposition_ok, &
       deposition_err_closure
  use stellar_ramses_bridge, only: deposit_source_to_uold, ramses_bridge_ok, &
       ramses_bridge_err_closure
  use stellar_ramses_field_map, only: stellar_field_map_t, clear_field_map
  use stellar_ramses_mapped_bridge, only: deposit_source_to_uold_mapped
  implicit none

  integer, parameter :: nvar = 17
  type(stellar_source_t) :: source
  type(stellar_field_map_t) :: field_map
  real(stellar_dp) :: volume(1), weights(1), gas_density(1)
  real(stellar_dp) :: gas_metal_density(1), gas_energy(1), gas_momentum(3,1)
  real(stellar_dp) :: gas_elements(n_stellar_elements,1), uold(nvar,1)
  integer :: element_var(n_stellar_elements), momentum_var(3)
  integer :: ierr, element

  call clear_source(source)
  source%returned_mass = 10.0_stellar_dp
  source%ejected_mass(elem_h) = 4.0_stellar_dp
  source%ejected_mass(elem_he) = 3.0_stellar_dp
  source%ejected_mass(elem_c) = 1.0_stellar_dp
  source%ejected_mass(elem_fe) = 0.5_stellar_dp
  volume = 1.0_stellar_dp
  weights = 1.0_stellar_dp
  gas_density = 0.0_stellar_dp
  gas_metal_density = 0.0_stellar_dp
  gas_energy = 0.0_stellar_dp
  gas_momentum = 0.0_stellar_dp
  gas_elements = 0.0_stellar_dp

  call deposit_stellar_source(source, 1, volume, weights, gas_density, &
       gas_elements, gas_energy, gas_momentum, 1.0e-10_stellar_dp, ierr, &
       gas_metal_density)
  if (ierr /= deposition_ok) error stop 1
  call expect_close(gas_density(1), 10.0_stellar_dp, 2)
  call expect_close(sum(gas_elements(:,1)), 8.5_stellar_dp, 3)
  call expect_close(gas_metal_density(1), 3.0_stellar_dp, 4)

  momentum_var = (/2, 3, 4/)
  do element = 1, n_stellar_elements
     element_var(element) = 6 + element
  end do
  uold = 0.0_stellar_dp
  call deposit_source_to_uold(source, nvar, 1, volume, weights, 1, 5, &
       momentum_var, element_var, uold, 1.0e-10_stellar_dp, ierr, 6)
  if (ierr /= ramses_bridge_ok) error stop 5
  call expect_close(uold(1,1), 10.0_stellar_dp, 6)
  call expect_close(uold(6,1), 3.0_stellar_dp, 7)
  call expect_close(sum(uold(7:17,1)), 8.5_stellar_dp, 8)

  call clear_field_map(field_map)
  field_map%density_index = 1
  field_map%momentum_index = momentum_var
  field_map%energy_index = 5
  field_map%total_metal_index = 6
  field_map%element_index = element_var
  uold = 0.0_stellar_dp
  call deposit_source_to_uold_mapped(uold, nvar, 1, 3, volume, weights, &
       source, field_map, ierr)
  if (ierr /= 0) error stop 9
  call expect_close(uold(1,1), 10.0_stellar_dp, 10)
  call expect_close(uold(6,1), 3.0_stellar_dp, 11)
  call expect_close(sum(uold(7:17,1)), 8.5_stellar_dp, 12)

  source%ejected_mass = 0.0_stellar_dp
  source%ejected_mass(elem_h) = 10.1_stellar_dp
  call deposit_stellar_source(source, 1, volume, weights, gas_density, &
       gas_elements, gas_energy, gas_momentum, 1.0e-10_stellar_dp, ierr, &
       gas_metal_density)
  if (ierr /= deposition_err_closure) error stop 13
  call deposit_source_to_uold(source, nvar, 1, volume, weights, 1, 5, &
       momentum_var, element_var, uold, 1.0e-10_stellar_dp, ierr, 6)
  if (ierr /= ramses_bridge_err_closure) error stop 14
  call deposit_source_to_uold_mapped(uold, nvar, 1, 3, volume, weights, &
       source, field_map, ierr)
  if (ierr /= 7) error stop 15

  write(*,'(A)') 'stellar residual deposition: PASS'

contains

  subroutine expect_close(value, expected, code)
    real(stellar_dp), intent(in) :: value, expected
    integer, intent(in) :: code
    if (abs(value-expected) > 1.0e-12_stellar_dp) error stop code
  end subroutine expect_close

end program test_stellar_residual_deposition

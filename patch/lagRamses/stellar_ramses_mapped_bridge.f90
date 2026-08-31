module stellar_ramses_mapped_bridge
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, elem_c
  use stellar_enrichment_contract, only: stellar_source_t
  use stellar_ramses_field_map, only: stellar_field_map_t, validate_field_map
  implicit none
  private

  public :: deposit_source_to_uold_mapped

contains

  subroutine deposit_source_to_uold_mapped(uold, nvar, ncell, ndim, cell_volume, &
       weights, source, field_map, ierr, message)
    integer, intent(in) :: nvar, ncell, ndim
    real(stellar_dp), intent(inout) :: uold(nvar, ncell)
    real(stellar_dp), intent(in) :: cell_volume(ncell)
    real(stellar_dp), intent(in) :: weights(ncell)
    type(stellar_source_t), intent(in) :: source
    type(stellar_field_map_t), intent(in) :: field_map
    integer, intent(out) :: ierr
    character(len=*), intent(out), optional :: message
    real(stellar_dp) :: weight_sum, volume_factor
    integer :: i, j, map_ierr
    character(len=256) :: map_message

    ierr = 0
    if (present(message)) message = ''

    call validate_field_map(field_map, nvar, ndim, map_ierr, map_message)
    if (map_ierr /= 0) then
       ierr = 100 + map_ierr
       if (present(message)) message = trim(map_message)
       return
    end if

    if (ncell <= 0) then
       ierr = 1
       if (present(message)) message = 'ncell must be positive'
       return
    end if
    if (any(weights < 0.0_stellar_dp)) then
       ierr = 2
       if (present(message)) message = 'deposition weights must be nonnegative'
       return
    end if
    if (any(cell_volume <= 0.0_stellar_dp)) then
       ierr = 3
       if (present(message)) message = 'cell volumes must be positive'
       return
    end if

    weight_sum = sum(weights)
    if (weight_sum <= 0.0_stellar_dp) then
       ierr = 4
       if (present(message)) message = 'deposition weights must have a positive sum'
       return
    end if

    if (source%returned_mass < 0.0_stellar_dp) then
       ierr = 5
       if (present(message)) message = 'returned stellar mass must be nonnegative'
       return
    end if
    if (any(source%ejected_mass < 0.0_stellar_dp)) then
       ierr = 6
       if (present(message)) message = 'actual ejecta masses must be nonnegative'
       return
    end if
    if (abs(sum(source%ejected_mass) - source%returned_mass) > &
        1.0e-10_stellar_dp * max(1.0_stellar_dp, source%returned_mass)) then
       ierr = 7
       if (present(message)) message = 'element ejecta do not close to returned mass'
       return
    end if

    do i = 1, ncell
       if (weights(i) == 0.0_stellar_dp) cycle
       volume_factor = weights(i) / (weight_sum * cell_volume(i))
       uold(field_map%density_index, i) = uold(field_map%density_index, i) + &
            volume_factor * source%returned_mass
       uold(field_map%energy_index, i) = uold(field_map%energy_index, i) + &
            volume_factor * source%energy
       do j = 1, ndim
          uold(field_map%momentum_index(j), i) = &
               uold(field_map%momentum_index(j), i) + &
               volume_factor * source%momentum(j)
       end do
       if (field_map%total_metal_index /= 0) then
          uold(field_map%total_metal_index, i) = &
               uold(field_map%total_metal_index, i) + volume_factor * &
               sum(source%ejected_mass(elem_c:n_stellar_elements))
       end if
       do j = 1, n_stellar_elements
          if (field_map%element_index(j) == 0) cycle
          uold(field_map%element_index(j), i) = &
               uold(field_map%element_index(j), i) + &
               volume_factor * source%ejected_mass(j)
       end do
    end do
  end subroutine deposit_source_to_uold_mapped

end module stellar_ramses_mapped_bridge

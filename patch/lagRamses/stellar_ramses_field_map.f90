module stellar_ramses_field_map
  use stellar_enrichment_config, only: n_stellar_elements
  implicit none
  private

  public :: stellar_field_map_t
  public :: clear_field_map
  public :: validate_field_map

  type :: stellar_field_map_t
     integer :: density_index = 0
     integer :: energy_index = 0
     integer :: momentum_index(3) = 0
     integer :: total_metal_index = 0
     integer :: element_index(n_stellar_elements) = 0
     logical :: volume_is_physical = .true.
  end type stellar_field_map_t

contains

  subroutine clear_field_map(field_map)
    type(stellar_field_map_t), intent(out) :: field_map

    field_map%density_index = 0
    field_map%energy_index = 0
    field_map%momentum_index = 0
    field_map%total_metal_index = 0
    field_map%element_index = 0
    field_map%volume_is_physical = .true.
  end subroutine clear_field_map

  subroutine validate_field_map(field_map, nvar, ndim, ierr, message)
    type(stellar_field_map_t), intent(in) :: field_map
    integer, intent(in) :: nvar
    integer, intent(in) :: ndim
    integer, intent(out) :: ierr
    character(len=*), intent(out), optional :: message
    integer :: i, j, index

    ierr = 0
    if (present(message)) message = ''

    if (nvar <= 0) then
       ierr = 1
       if (present(message)) message = 'nvar must be positive'
       return
    end if
    if (ndim < 1 .or. ndim > 3) then
       ierr = 2
       if (present(message)) message = 'ndim must be in the range 1..3'
       return
    end if

    if (.not. valid_index(field_map%density_index, nvar)) then
       ierr = 3
       if (present(message)) message = 'invalid density field index'
       return
    end if
    if (.not. valid_index(field_map%energy_index, nvar)) then
       ierr = 4
       if (present(message)) message = 'invalid energy field index'
       return
    end if
    if (field_map%energy_index == field_map%density_index) then
       ierr = 5
       if (present(message)) message = 'density and energy fields overlap'
       return
    end if

    if (field_map%total_metal_index /= 0) then
       if (.not. valid_index(field_map%total_metal_index, nvar)) then
          ierr = 13
          if (present(message)) message = 'invalid total-metal field index'
          return
       end if
       if (field_map%total_metal_index == field_map%density_index .or. &
           field_map%total_metal_index == field_map%energy_index) then
          ierr = 14
          if (present(message)) message = 'total-metal field overlaps a core field'
          return
       end if
    end if

    do i = 1, ndim
       index = field_map%momentum_index(i)
       if (.not. valid_index(index, nvar)) then
          ierr = 6
          if (present(message)) write(message, '(a,i0)') &
               'invalid momentum field index for component ', i
          return
       end if
       if (index == field_map%density_index .or. &
           index == field_map%energy_index) then
          ierr = 7
          if (present(message)) write(message, '(a,i0)') &
               'momentum field overlaps a core field for component ', i
          return
       end if
       if (field_map%total_metal_index /= 0 .and. &
           index == field_map%total_metal_index) then
          ierr = 15
          if (present(message)) write(message, '(a,i0)') &
               'momentum field overlaps total-metal field for component ', i
          return
       end if
       do j = 1, i - 1
          if (index == field_map%momentum_index(j)) then
             ierr = 8
             if (present(message)) write(message, '(a,i0)') &
                  'momentum field components overlap at component ', i
             return
          end if
       end do
    end do

    do i = 1, n_stellar_elements
       index = field_map%element_index(i)
       if (index == 0) cycle
       if (.not. valid_index(index, nvar)) then
          ierr = 9
          if (present(message)) write(message, '(a,i0)') &
               'invalid element field index for element ', i
          return
       end if
       if (index == field_map%density_index .or. &
           index == field_map%energy_index) then
          ierr = 10
          if (present(message)) write(message, '(a,i0)') &
               'element field overlaps a core field for element ', i
          return
       end if
       if (field_map%total_metal_index /= 0 .and. &
           index == field_map%total_metal_index) then
          ierr = 16
          if (present(message)) write(message, '(a,i0)') &
               'element field overlaps total-metal field for element ', i
          return
       end if
       do j = 1, ndim
          if (index == field_map%momentum_index(j)) then
             ierr = 11
             if (present(message)) write(message, '(a,i0)') &
                  'element field overlaps momentum component for element ', i
             return
          end if
       end do
       do j = 1, i - 1
          if (index /= 0 .and. index == field_map%element_index(j)) then
             ierr = 12
             if (present(message)) write(message, '(a,i0)') &
                  'element field index is duplicated for element ', i
             return
          end if
       end do
    end do
  end subroutine validate_field_map

  logical function valid_index(index, nvar)
    integer, intent(in) :: index
    integer, intent(in) :: nvar

    valid_index = index >= 1 .and. index <= nvar
  end function valid_index

end module stellar_ramses_field_map

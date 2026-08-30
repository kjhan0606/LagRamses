! Persistent local S_N state, indexed by the RAMSES leaf-cell identifier.
module snrt_state
  use amr_parameters, only: MAXLEVEL, dp
  use iso_c_binding, only: c_float
#include "amr_index.h"
  implicit none

  integer, parameter, public :: snrt_ndirection = 80
  ! Four spectral groups are carried by the persistent state.  Group
  ! energies and cross sections are supplied by the NLTE/chemistry layer.
  integer, parameter, public :: snrt_ngroups = 4
  real(dp), parameter, public :: snrt_group_mean_energy_ev(snrt_ngroups) = &
       (/ 18.0d0, 35.0d0, 70.0d0, 200.0d0 /)
  real(dp), parameter, public :: snrt_group_energy_fraction(snrt_ngroups) = &
       (/ 0.45d0, 0.30d0, 0.20d0, 0.05d0 /)
  real(dp), parameter, public :: snrt_group_cross_section_cm2(snrt_ngroups) = &
       (/ 6.3d-18*(18.0d0/13.6d0)**(-3.0d0), &
          6.3d-18*(35.0d0/13.6d0)**(-3.0d0), &
          6.3d-18*(70.0d0/13.6d0)**(-3.0d0), &
          6.3d-18*(200.0d0/13.6d0)**(-3.0d0) /)
  integer, public :: snrt_nslot = 0
  real(c_float), allocatable, public :: snrt_intensity(:, :, :)
  real(c_float), allocatable, public :: snrt_neutral_fraction(:)

  integer :: snrt_capacity = 0
  integer, allocatable :: snrt_slot_of_cell(:)
  integer, allocatable :: snrt_cell_id(:)

  private
  public :: snrt_state_sync_level, snrt_state_get_slot, snrt_state_get_cell
  public :: snrt_state_checkpoint_write, snrt_state_checkpoint_read

contains

  function snrt_state_get_slot(icell) result(islot)
    integer, intent(in) :: icell
    integer :: islot

    islot = 0
    if (.not. allocated(snrt_slot_of_cell)) return
    if (icell < 1 .or. icell > size(snrt_slot_of_cell)) return
    islot = snrt_slot_of_cell(icell)
  end function snrt_state_get_slot

  function snrt_state_get_cell(islot) result(icell)
    integer, intent(in) :: islot
    integer :: icell

    icell = 0
    if (.not. allocated(snrt_cell_id)) return
    if (islot < 1 .or. islot > snrt_nslot) return
    icell = snrt_cell_id(islot)
  end function snrt_state_get_cell

  subroutine snrt_state_checkpoint_write(unit_id, ierr)
    integer, intent(in) :: unit_id
    integer, intent(out) :: ierr
    integer, parameter :: checkpoint_version = 3
    integer :: ios

    ierr = 0
    write(unit_id, iostat=ios) checkpoint_version, snrt_ndirection, &
         snrt_ngroups, snrt_nslot
    if (ios /= 0) then
       ierr = 1
       return
    end if
    if (snrt_nslot <= 0) return
    write(unit_id, iostat=ios) snrt_cell_id(1:snrt_nslot)
    if (ios /= 0) then
       ierr = 2
       return
    end if
    write(unit_id, iostat=ios) snrt_intensity(:,:,1:snrt_nslot)
    if (ios /= 0) then
       ierr = 3
       return
    end if
    write(unit_id, iostat=ios) snrt_neutral_fraction(1:snrt_nslot)
    if (ios /= 0) ierr = 4
  end subroutine snrt_state_checkpoint_write

  subroutine snrt_state_checkpoint_read(unit_id, ierr)
    integer, intent(in) :: unit_id
    integer, intent(out) :: ierr
    integer, parameter :: checkpoint_version = 3
    integer :: ios, version, ndirection_file, ngroups_file, nslot_file
    integer :: islot, icell
    integer, allocatable :: saved_cell_id(:)
    real(c_float), allocatable :: saved_intensity(:,:,:), saved_neutral(:)

    ierr = 0
    read(unit_id, iostat=ios) version, ndirection_file, ngroups_file, nslot_file
    if (ios /= 0) then
       ierr = 1
       return
    end if
    if (version /= checkpoint_version .or. ndirection_file /= snrt_ndirection .or. &
         ngroups_file /= snrt_ngroups .or. nslot_file < 0) then
       ierr = 2
       return
    end if
    if (nslot_file == 0) then
       snrt_nslot = 0
       if (allocated(snrt_slot_of_cell)) snrt_slot_of_cell = 0
       return
    end if

    allocate(saved_cell_id(nslot_file), &
         saved_intensity(snrt_ndirection,snrt_ngroups,nslot_file), &
         saved_neutral(nslot_file))
    read(unit_id, iostat=ios) saved_cell_id
    if (ios /= 0) then
       ierr = 3
       deallocate(saved_cell_id, saved_intensity, saved_neutral)
       return
    end if
    read(unit_id, iostat=ios) saved_intensity
    if (ios /= 0) then
       ierr = 4
       deallocate(saved_cell_id, saved_intensity, saved_neutral)
      return
    end if
    read(unit_id, iostat=ios) saved_neutral
    if (ios /= 0) then
       ierr = 5
       deallocate(saved_cell_id, saved_intensity, saved_neutral)
       return
    end if

    if (.not. allocated(snrt_slot_of_cell)) call snrt_state_initialize()
    if (.not. allocated(snrt_slot_of_cell)) then
       ierr = 6
       deallocate(saved_cell_id, saved_intensity, saved_neutral)
       return
    end if
    call snrt_state_grow(nslot_file)
    if (snrt_nslot > 0) snrt_intensity(:,:,1:snrt_nslot) = 0.0_c_float
    snrt_slot_of_cell = 0
    snrt_nslot = nslot_file
    snrt_cell_id(1:nslot_file) = saved_cell_id
    snrt_intensity(:,:,1:nslot_file) = saved_intensity
    snrt_neutral_fraction(1:nslot_file) = saved_neutral
    do islot = 1, nslot_file
       icell = saved_cell_id(islot)
       if (icell >= 1 .and. icell <= size(snrt_slot_of_cell)) &
            snrt_slot_of_cell(icell) = islot
    end do
    deallocate(saved_cell_id, saved_intensity, saved_neutral)
  end subroutine snrt_state_checkpoint_read

  subroutine snrt_state_sync_level(ilevel, nleaf, nnew)
    use amr_commons, only: active, ncoarse, twotondim, son
    implicit none

    integer, intent(in) :: ilevel
    integer, intent(out) :: nleaf, nnew
    integer :: i, ind, igrid, icell, islot

    call snrt_state_initialize()
    nleaf = 0
    nnew = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ICELL_OF(igrid, ind)
          if (son(icell) /= 0) cycle
          nleaf = nleaf + 1
          islot = snrt_slot_of_cell(icell)
          if (islot == 0) then
             call snrt_state_grow(snrt_nslot + 1)
             snrt_nslot = snrt_nslot + 1
             islot = snrt_nslot
             snrt_slot_of_cell(icell) = islot
             snrt_cell_id(islot) = icell
             snrt_intensity(:, :, islot) = 0.0_c_float
             nnew = nnew + 1
          endif
       enddo
    enddo
  end subroutine snrt_state_sync_level

  subroutine snrt_state_initialize()
    use amr_commons, only: ncoarse, ngridmax, twotondim
    implicit none

    integer :: cell_capacity

    if (allocated(snrt_slot_of_cell)) return
    cell_capacity = ICELL_OF(ngridmax, twotondim)
    allocate(snrt_slot_of_cell(cell_capacity))
    snrt_slot_of_cell = 0
  end subroutine snrt_state_initialize

  subroutine snrt_state_grow(required)
    implicit none

    integer, intent(in) :: required
    integer :: next_capacity
    integer, allocatable :: next_cell_id(:)
    real(c_float), allocatable :: next_intensity(:, :, :)
    real(c_float), allocatable :: next_neutral_fraction(:)

    if (required <= snrt_capacity) return
    next_capacity = max(required, max(1024, 2 * snrt_capacity))
    allocate(next_cell_id(next_capacity))
    allocate(next_intensity(snrt_ndirection, snrt_ngroups, next_capacity))
    allocate(next_neutral_fraction(next_capacity))
    next_cell_id = 0
    next_intensity = 0.0_c_float
    next_neutral_fraction = 1.0_c_float
    if (snrt_nslot > 0) then
       next_cell_id(1:snrt_nslot) = snrt_cell_id(1:snrt_nslot)
       next_intensity(:, :, 1:snrt_nslot) = snrt_intensity(:, :, 1:snrt_nslot)
       if (allocated(snrt_neutral_fraction)) then
          next_neutral_fraction(1:snrt_nslot) = &
               snrt_neutral_fraction(1:snrt_nslot)
       end if
    endif
    call move_alloc(next_cell_id, snrt_cell_id)
    call move_alloc(next_intensity, snrt_intensity)
    call move_alloc(next_neutral_fraction, snrt_neutral_fraction)
    snrt_capacity = next_capacity
  end subroutine snrt_state_grow

end module snrt_state

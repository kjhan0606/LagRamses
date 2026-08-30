module snrt_amr_topology
  ! Construct the local same-level portion of the AMR transport stencil.
  ! Faces touching coarse-fine interfaces or another MPI rank are classified
  ! and retained for the interface-state and flux-register stages.
  use amr_parameters
  use amr_commons
  use snrt_state, only: snrt_state_sync_level, snrt_state_get_slot, &
       snrt_state_get_cell, snrt_nslot, snrt_intensity
#include "amr_index.h"
  implicit none

  integer, parameter :: SNRT_FACE_LOCAL          = 0
  integer, parameter :: SNRT_FACE_PHYSICAL       = 1
  integer, parameter :: SNRT_FACE_MPI            = 2
  integer, parameter :: SNRT_FACE_COARSE_TO_FINE = 3
  integer, parameter :: SNRT_FACE_FINE_TO_COARSE = 4
  integer, parameter :: SNRT_FACE_UNMAPPED      = 5
  integer, allocatable, save :: snrt_face_kind(:,:)
  integer, allocatable, save :: snrt_face_cell(:,:)

contains

  subroutine snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
       neighbor, nleaf, n_interface_face)
    integer, intent(in) :: ilevel
    integer, allocatable, intent(out) :: leaf_cell(:), leaf_slot(:)
    integer, allocatable, intent(out) :: neighbor(:,:)
    integer, intent(out) :: nleaf, n_interface_face
    integer :: i, ind, face, igrid, icell, islot, ilocal, nnew
    integer :: neighbor_cell, neighbor_slot, local_cell
    integer :: j, direct_match
    integer :: hash_size, hash_pos, hash_probe
    integer(kind=8) :: hash_code
    integer, allocatable :: cell_hash_key(:), cell_hash_value(:)
    integer, allocatable :: raw_neighbor(:,:)
    logical, save :: report_hash_miss = .false.
    integer, dimension(1:nvector) :: igrid_one
    integer, dimension(1:nvector,0:twondim) :: igridn
    integer, dimension(1:nvector,1:twondim) :: indn

    nleaf = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ICELL_OF(igrid, ind)
          if (son(icell) == 0) nleaf = nleaf + 1
       end do
    end do

    call snrt_state_sync_level(ilevel, nleaf, nnew)
    allocate(leaf_cell(nleaf), leaf_slot(nleaf), neighbor(2*ndim,nleaf))
    neighbor = 0
    if (nleaf == 0) then
       n_interface_face = 0
       return
    end if

    ! The persistent slot map can contain cell IDs from an earlier AMR mesh.
    ! Build a per-call cell-ID map so a stale slot can never turn a current
    ! same-level neighbor into an artificial unmapped face.
    hash_size = 1
    do while (hash_size < max(2,2*nleaf))
       hash_size = 2*hash_size
    end do
    allocate(cell_hash_key(hash_size), cell_hash_value(hash_size), &
         raw_neighbor(2*ndim,nleaf))
    cell_hash_key = 0
    cell_hash_value = 0
    raw_neighbor = 0
    ilocal = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ICELL_OF(igrid, ind)
          if (son(icell) /= 0) cycle
          ilocal = ilocal + 1
          islot = snrt_state_get_slot(icell)
          leaf_cell(ilocal) = icell
          leaf_slot(ilocal) = islot
          hash_code = modulo(int(icell,8)*1000003_8 + 17_8, int(hash_size,8))
          hash_pos = 1 + int(hash_code)
          do hash_probe = 1, hash_size
             if (cell_hash_key(hash_pos) == 0) then
                cell_hash_key(hash_pos) = icell
                cell_hash_value(hash_pos) = ilocal
                exit
             end if
             hash_pos = 1 + mod(hash_pos,hash_size)
          end do
       end do
    end do

    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       ! flag_fine releases nbor_active_cache before amr_step reaches SNRT.
       ! Recompute this one-grid stencil through the canonical RAMSES helper
       ! instead of dereferencing that transient cache after it is deallocated.
       igrid_one = 0
       igrid_one(1) = igrid
       igridn = 0
       call getnborgrids_check(igrid_one, igridn, 1)
       do ind = 1, twotondim
          icell = ICELL_OF(igrid, ind)
          if (son(icell) /= 0) cycle
          islot = snrt_state_get_slot(icell)
          hash_code = modulo(int(icell,8)*1000003_8 + 17_8, int(hash_size,8))
          hash_pos = 1 + int(hash_code)
          local_cell = 0
          do hash_probe = 1, hash_size
             if (cell_hash_key(hash_pos) == icell) then
                local_cell = cell_hash_value(hash_pos)
                exit
             end if
             if (cell_hash_key(hash_pos) == 0) exit
             hash_pos = 1 + mod(hash_pos,hash_size)
          end do
          ilocal = local_cell
          if (ilocal <= 0) cycle
          call getnborcells(igridn, ind, indn, 1)
          do face = 1, 2*ndim
             raw_neighbor(face,ilocal) = indn(1,face)
             neighbor_cell = indn(1,face)
             if (neighbor_cell > 0 .and. son(neighbor_cell) == 0) then
                hash_code = modulo(int(neighbor_cell,8)*1000003_8 + 17_8, &
                     int(hash_size,8))
                hash_pos = 1 + int(hash_code)
                local_cell = 0
                do hash_probe = 1, hash_size
                   if (cell_hash_key(hash_pos) == neighbor_cell) then
                      local_cell = cell_hash_value(hash_pos)
                      exit
                   end if
                   if (cell_hash_key(hash_pos) == 0) exit
                   hash_pos = 1 + mod(hash_pos,hash_size)
                end do
                if (local_cell == 0 .and. cpu_map(neighbor_cell) == myid .and. &
                     .not. report_hash_miss) then
                   direct_match = 0
                   do j = 1, nleaf
                      if (leaf_cell(j) == neighbor_cell) then
                         direct_match = j
                         exit
                      end if
                   end do
                   write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                        ' SNRT hash miss rank=', myid, ' current=', icell, &
                        ' neighbor=', neighbor_cell, ' nleaf=', nleaf, &
                        ' hash_size=', hash_size, ' direct_match=', direct_match, &
                        ' son=', son(neighbor_cell)
                   report_hash_miss = .true.
                end if
                if (local_cell > 0) neighbor(face,ilocal) = local_cell
             end if
          end do
       end do
    end do

    call snrt_amr_classify_faces(leaf_cell, nleaf, neighbor, raw_neighbor)
    n_interface_face = count(neighbor(1:2*ndim,:) == 0)
    deallocate(cell_hash_key, cell_hash_value, raw_neighbor)
  end subroutine snrt_amr_build_same_level_neighbors

  subroutine snrt_amr_classify_faces(leaf_cell, nleaf, neighbor, raw_neighbor)
    use amr_commons
    implicit none
    integer, intent(in) :: nleaf
    integer, dimension(:), intent(in) :: leaf_cell
    integer, dimension(:,:), intent(in) :: neighbor
    integer, dimension(:,:), intent(in) :: raw_neighbor
    integer :: ilocal, face, idim, inbor, ind, igrid
    integer :: parent_cell, neighbor_grid, neighbor_cell, neighbor_ind
    integer :: bitmask, current_bit, cross_boundary, grid_face
    integer :: cell_index, nunmapped
    logical, save :: report_unmapped = .false.

    if (allocated(snrt_face_kind)) deallocate(snrt_face_kind)
    if (allocated(snrt_face_cell)) deallocate(snrt_face_cell)
    allocate(snrt_face_kind(2*ndim,nleaf))
    allocate(snrt_face_cell(2*ndim,nleaf))
    snrt_face_kind = SNRT_FACE_UNMAPPED
    snrt_face_cell = 0
    nunmapped = 0

    do ilocal=1,nleaf
       cell_index = leaf_cell(ilocal)
       ind = ICHILD_OF(cell_index)
       igrid = IGRID_OF(cell_index)

       do face=1,2*ndim
          if (neighbor(face,ilocal) > 0) then
             snrt_face_kind(face,ilocal) = SNRT_FACE_LOCAL
             cycle
          end if

          ! A raw same-level neighbor is authoritative when the canonical
          ! cell stencil found it.  It may be remote, or it may be a current
          ! coarse cell refined at the next level.
          neighbor_cell = raw_neighbor(face,ilocal)
          if (neighbor_cell > 0) then
             snrt_face_cell(face,ilocal) = neighbor_cell
             if (son(neighbor_cell) > 0) then
                snrt_face_kind(face,ilocal) = SNRT_FACE_COARSE_TO_FINE
             else if (cpu_map(neighbor_cell) /= myid) then
                snrt_face_kind(face,ilocal) = SNRT_FACE_MPI
             else
                snrt_face_kind(face,ilocal) = SNRT_FACE_UNMAPPED
             end if
             if (snrt_face_kind(face,ilocal) /= SNRT_FACE_UNMAPPED) cycle
          end if

          idim = (face-1)/2 + 1
          inbor = mod(face-1,2) + 1
          bitmask = 2**(idim-1)
          current_bit = mod((ind-1)/bitmask,2)
          cross_boundary = 0
          if ((inbor == 1 .and. current_bit == 0) .or. &
               (inbor == 2 .and. current_bit == 1)) cross_boundary = 1

          parent_cell = 0
          neighbor_grid = igrid
          if (cross_boundary == 1) then
             grid_face = 2*(idim-1) + inbor
             parent_cell = nbor(igrid,grid_face)
             if (parent_cell <= 0) then
                snrt_face_kind(face,ilocal) = SNRT_FACE_PHYSICAL
                cycle
             end if
             neighbor_grid = son(parent_cell)
          end if

          if (neighbor_grid <= 0) then
             snrt_face_cell(face,ilocal) = parent_cell
             snrt_face_kind(face,ilocal) = SNRT_FACE_FINE_TO_COARSE
             cycle
          end if

          neighbor_ind = ind
          if (cross_boundary == 1) then
             if (inbor == 1) then
                neighbor_ind = ind + bitmask
             else
                neighbor_ind = ind - bitmask
             end if
          else
             if (inbor == 1) then
                neighbor_ind = ind - bitmask
             else
                neighbor_ind = ind + bitmask
             end if
          end if
          neighbor_cell = ICELL_OF(neighbor_grid, neighbor_ind)
          snrt_face_cell(face,ilocal) = neighbor_cell
          if (son(neighbor_cell) > 0) then
             snrt_face_kind(face,ilocal) = SNRT_FACE_COARSE_TO_FINE
          else if (cpu_map(neighbor_cell) /= myid) then
             snrt_face_kind(face,ilocal) = SNRT_FACE_MPI
          else
             snrt_face_kind(face,ilocal) = SNRT_FACE_UNMAPPED
             nunmapped = nunmapped + 1
             if (.not. report_unmapped) then
                write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                     ' SNRT topology unmapped rank=', myid, ' level_cell=', &
                     cell_index, ' face=', face, ' ind=', ind, ' igrid=', igrid, &
                     ' neighbor_cell=', neighbor_cell, ' parent=', parent_cell, &
                     ' neighbor_grid=', neighbor_grid, ' son=', son(neighbor_cell), &
                     ' owner=', cpu_map(neighbor_cell), ' slot=', &
                     snrt_state_get_slot(neighbor_cell)
                report_unmapped = .true.
             end if
          end if
       end do
    end do
    if (nunmapped > 0) then
       write(*,'(A,I0,A,I0)') ' SNRT topology unmapped count rank=', myid, &
            ' count=', nunmapped
    end if
  end subroutine snrt_amr_classify_faces

  subroutine snrt_amr_exchange_interface_state(ilevel, leaf_cell, state, &
       ghost_kind, ghost_cell, ghost_face, ghost_state, ierr)
    use iso_c_binding, only: c_float
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    implicit none
    integer, intent(in) :: ilevel
    integer, intent(in) :: leaf_cell(:), ghost_kind(:), ghost_cell(:), ghost_face(:)
    real(c_float), intent(in) :: state(:,:,:)
    real(c_float), intent(out) :: ghost_state(:,:,:)
    integer, intent(out) :: ierr
    real(dp), allocatable :: field(:)
    integer :: ilocal, idir, igroup, ighost, nfield, islot, icell
    integer :: child_grid, child_cell, child, idim, inbor, bit, target_bit
    integer :: nchild
    integer :: local_need_same, local_need_coarse, local_need_fine
    integer :: global_need_same, global_need_coarse, global_need_fine, info
    real(dp) :: face_average
    logical :: need_same, need_coarse, need_fine

    ierr = 0
    ghost_state = 0.0_c_float
    if (size(state,1) < size(leaf_cell) .or. size(ghost_kind) < size(ghost_cell) .or. &
         size(ghost_face) < size(ghost_cell) .or. size(ghost_state,1) < size(ghost_cell) .or. &
         size(ghost_state,2) < size(state,2) .or. size(ghost_state,3) < size(state,3)) then
       ierr = 1
       return
    end if

    nfield = ICELL_OF(ngridmax, twotondim)
    allocate(field(nfield))
    need_same = .false.
    need_coarse = .false.
    need_fine = .false.
    do ighost = 1, size(ghost_cell)
       select case (ghost_kind(ighost))
       case (SNRT_FACE_MPI)
          need_same = .true.
       case (SNRT_FACE_FINE_TO_COARSE)
          need_coarse = ilevel > 1
       case (SNRT_FACE_COARSE_TO_FINE)
          need_fine = ilevel < nlevelmax
       case default
          ierr = 2
          deallocate(field)
          return
       end select
    end do

    ! make_virtual_fine_dp is collective.  A rank can have no local
    ! interface of a given kind even though another rank does, so exchange
    ! decisions must be global before entering the RAMSES communicator.
    local_need_same = 0
    local_need_coarse = 0
    local_need_fine = 0
    if (need_same) local_need_same = 1
    if (need_coarse) local_need_coarse = 1
    if (need_fine) local_need_fine = 1
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_need_same, global_need_same, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, info)
    call MPI_ALLREDUCE(local_need_coarse, global_need_coarse, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, info)
    call MPI_ALLREDUCE(local_need_fine, global_need_fine, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, info)
#else
    global_need_same = local_need_same
    global_need_coarse = local_need_coarse
    global_need_fine = local_need_fine
#endif
    need_same = global_need_same > 0
    need_coarse = global_need_coarse > 0
    need_fine = global_need_fine > 0
    if (.not. need_same .and. .not. need_coarse .and. .not. need_fine) return

    do igroup = 1, size(state,3)
       do idir = 1, size(state,2)
          field = 0.0d0
          do islot = 1, snrt_nslot
             icell = snrt_state_get_cell(islot)
             if (icell >= 1 .and. icell <= nfield) then
                field(icell) = real(snrt_intensity(idir,igroup,islot),dp)
             end if
          end do
          do ilocal = 1, size(leaf_cell)
             if (leaf_cell(ilocal) >= 1 .and. leaf_cell(ilocal) <= nfield) then
                field(leaf_cell(ilocal)) = real(state(ilocal,idir,igroup),dp)
             end if
          end do
          if (need_same) call make_virtual_fine_dp(field, ilevel)
          if (need_coarse) call make_virtual_fine_dp(field, ilevel-1)
          if (need_fine) call make_virtual_fine_dp(field, ilevel+1)
          do ighost = 1, size(ghost_cell)
             select case (ghost_kind(ighost))
             case (SNRT_FACE_MPI, SNRT_FACE_FINE_TO_COARSE)
                if (ghost_cell(ighost) < 1 .or. ghost_cell(ighost) > nfield) then
                   ierr = 3
                   cycle
                end if
                ghost_state(ighost,idir,igroup) = real(field(ghost_cell(ighost)),c_float)
             case (SNRT_FACE_COARSE_TO_FINE)
                child_grid = son(ghost_cell(ighost))
                if (child_grid <= 0) then
                   ierr = 4
                   cycle
                end if
                idim = (ghost_face(ighost)-1)/2 + 1
                inbor = mod(ghost_face(ighost)-1,2) + 1
                target_bit = 1
                if (inbor == 2) target_bit = 0
                face_average = 0.0d0
                nchild = 0
                do child = 1, twotondim
                   bit = mod((child-1)/(2**(idim-1)),2)
                   if (bit /= target_bit) cycle
                   child_cell = ICELL_OF(child_grid, child)
                   if (child_cell < 1 .or. child_cell > nfield) then
                      ierr = 5
                      cycle
                   end if
                   face_average = face_average + field(child_cell)
                   nchild = nchild + 1
                end do
                if (nchild > 0) then
                   ghost_state(ighost,idir,igroup) = real(face_average/dble(nchild),c_float)
                else
                   ierr = 6
                end if
             case default
                ierr = 7
             end select
          end do
       end do
    end do
    deallocate(field)
  end subroutine snrt_amr_exchange_interface_state

  subroutine snrt_amr_apply_coarse_flux_correction(ilevel, leaf_cell, state_work, &
       ghost_kind, ghost_cell, ghost_face, ghost_local, cdt_over_dx, direction_dp, ierr)
    use iso_c_binding, only: c_float
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    implicit none
    integer, intent(in) :: ilevel
    integer, intent(in) :: leaf_cell(:), ghost_kind(:), ghost_cell(:)
    integer, intent(in) :: ghost_face(:), ghost_local(:)
    real(c_float), intent(in) :: state_work(:,:,:)
    real(dp), intent(in) :: cdt_over_dx, direction_dp(:,:)
    integer, intent(out) :: ierr
    real(dp), allocatable :: correction(:)
    integer :: nfield, nface_child, ilocal, idir, igroup, ighost
    integer :: parent_cell, idim, fine_sign, coarse_sign, islot, icell
    integer :: local_has, global_has, info
    real(dp) :: mu, q_upstream, face_flux

    ierr = 0
    local_has = 0
    do ighost = 1, size(ghost_kind)
       if (ghost_kind(ighost) == SNRT_FACE_FINE_TO_COARSE) local_has = 1
    end do
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_has, global_has, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, info)
#else
    global_has = local_has
#endif
    if (global_has == 0) return
    if (ilevel <= 1 .or. size(state_work,1) < size(leaf_cell) .or. &
         size(direction_dp,1) < 1 .or. size(direction_dp,2) < 3) then
       ierr = 1
       return
    end if

    nfield = ICELL_OF(ngridmax, twotondim)
    nface_child = max(1,twotondim/2)
    allocate(correction(nfield))

    do igroup = 1, size(state_work,3)
       do idir = 1, size(state_work,2)
          correction = 0.0d0
          do ighost = 1, size(ghost_kind)
             if (ghost_kind(ighost) /= SNRT_FACE_FINE_TO_COARSE) cycle
             ilocal = ghost_local(ighost)
             parent_cell = ghost_cell(ighost)
             if (ilocal < 1 .or. ilocal > size(leaf_cell) .or. &
                  parent_cell < 1 .or. parent_cell > nfield .or. &
                  size(leaf_cell) + ighost > size(state_work,1)) then
                ierr = 2
                cycle
             end if
             idim = (ghost_face(ighost)-1)/2 + 1
             fine_sign = -1
             if (mod(ghost_face(ighost),2) == 0) fine_sign = 1
             coarse_sign = -fine_sign
             mu = direction_dp(idir,idim)
             if (mu*dble(fine_sign) >= 0.0d0) then
                q_upstream = real(state_work(ilocal,idir,igroup),dp)
             else
                q_upstream = real(state_work(size(leaf_cell)+ighost,idir,igroup),dp)
             end if
             face_flux = mu*q_upstream
             correction(parent_cell) = correction(parent_cell) - &
                  cdt_over_dx * 0.5d0/dble(nface_child) * &
                  dble(coarse_sign)*face_flux
          end do

          call make_virtual_reverse_dp(correction, ilevel-1)
          do islot = 1, snrt_nslot
             icell = snrt_state_get_cell(islot)
             if (icell >= 1 .and. icell <= nfield .and. correction(icell) /= 0.0d0) then
                snrt_intensity(idir,igroup,islot) = snrt_intensity(idir,igroup,islot) + &
                     real(correction(icell),c_float)
             end if
          end do
       end do
    end do
    deallocate(correction)
  end subroutine snrt_amr_apply_coarse_flux_correction

end module snrt_amr_topology

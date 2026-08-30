module snrt_amr_topology
  ! Construct the local same-level portion of the AMR transport stencil.
  ! Faces touching coarse-fine interfaces or another MPI rank remain zero;
  ! the later AMR flux-register stage owns those interfaces explicitly.
  use amr_parameters
  use amr_commons
  use snrt_state, only: snrt_state_sync_level, snrt_state_get_slot, snrt_nslot
  implicit none

contains

  subroutine snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
       neighbor, nleaf, n_interface_face)
    integer, intent(in) :: ilevel
    integer, allocatable, intent(out) :: leaf_cell(:), leaf_slot(:)
    integer, allocatable, intent(out) :: neighbor(:,:)
    integer, intent(out) :: nleaf, n_interface_face
    integer :: i, ind, face, igrid, icell, islot, ilocal, nnew
    integer :: neighbor_cell, neighbor_slot
    integer, allocatable :: local_of_slot(:)
    integer, dimension(1:nvector) :: igrid_one
    integer, dimension(1:nvector,0:twondim) :: igridn
    integer, dimension(1:nvector,1:twondim) :: indn

    nleaf = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ncoarse + (ind-1) * ngridmax + igrid
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

    allocate(local_of_slot(snrt_nslot))
    local_of_slot = 0
    ilocal = 0
    do i = 1, active(ilevel)%ngrid
       igrid = active(ilevel)%igrid(i)
       do ind = 1, twotondim
          icell = ncoarse + (ind-1) * ngridmax + igrid
          if (son(icell) /= 0) cycle
          ilocal = ilocal + 1
          islot = snrt_state_get_slot(icell)
          leaf_cell(ilocal) = icell
          leaf_slot(ilocal) = islot
          local_of_slot(islot) = ilocal
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
          icell = ncoarse + (ind-1) * ngridmax + igrid
          if (son(icell) /= 0) cycle
          islot = snrt_state_get_slot(icell)
          ilocal = local_of_slot(islot)
          call getnborcells(igridn, ind, indn, 1)
          do face = 1, 2*ndim
             neighbor_cell = indn(1,face)
             if (neighbor_cell > 0 .and. son(neighbor_cell) == 0) then
                neighbor_slot = snrt_state_get_slot(neighbor_cell)
                if (neighbor_slot > 0) neighbor(face,ilocal) = &
                     local_of_slot(neighbor_slot)
             end if
          end do
       end do
    end do

    n_interface_face = count(neighbor(1:2*ndim,:) == 0)
    deallocate(local_of_slot)
  end subroutine snrt_amr_build_same_level_neighbors

end module snrt_amr_topology

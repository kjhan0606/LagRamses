module snrt_transport_step
  ! AMR leaf-slot <-> GPU sparse-stencil adapter.  Only same-level local
  ! faces are populated at this stage; coarse-fine and MPI faces remain a
  ! separate conservative flux-register responsibility.
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection, snrt_ngroups, snrt_nslot, snrt_intensity
  use snrt_angular_quadrature, only: snrt_angular_init
  use snrt_amr_topology, only: snrt_amr_build_same_level_neighbors
  use snrt_cuda_sparse_transport_interface, only: snrt_cuda_upwind_sparse
  use snrt_cuda_limited_rt_step_interface, only: snrt_cuda_transport_absorb_limited
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step
  implicit none

contains

  subroutine snrt_transport_level(ilevel, cdt_over_dx, ierr, nleaf, &
       n_interface_face)
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: cdt_over_dx
    integer, intent(out) :: ierr, nleaf, n_interface_face
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:)
    real(dp) :: angular_cfl
    integer :: igroup, ilocal
    integer(c_int) :: cuda_ierr

    ierr = 0
    nleaf = 0
    n_interface_face = 0
    if (cdt_over_dx < 0.0d0) then
       ierr = 1
       return
    end if

    call snrt_angular_init(direction_dp, weight)
    angular_cfl = maxval(sum(abs(direction_dp), dim=2))
    if (cdt_over_dx * angular_cfl > 1.0d0) then
       ierr = 2
       return
    end if

    call snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
         neighbor, nleaf, n_interface_face)
    if (nleaf == 0) return

    allocate(neighbor_c(6,nleaf), packed(nleaf,snrt_ndirection))
    neighbor_c = int(neighbor, c_int)
    direction_c = real(transpose(direction_dp), c_float)

    do igroup = 1, snrt_ngroups
       do ilocal = 1, nleaf
          packed(ilocal,:) = snrt_intensity(:,igroup,leaf_slot(ilocal))
       end do
       cuda_ierr = snrt_cuda_upwind_sparse(packed, direction_c, neighbor_c, &
            int(nleaf,c_int), int(snrt_ndirection,c_int), &
            real(cdt_over_dx,c_float))
       if (cuda_ierr /= 0_c_int) then
          ierr = 100 + int(cuda_ierr)
          return
       end if
       do ilocal = 1, nleaf
          snrt_intensity(:,igroup,leaf_slot(ilocal)) = packed(ilocal,:)
       end do
    end do
  end subroutine snrt_transport_level

  subroutine snrt_transport_absorb_limited_level(ilevel, cdt_over_dx, &
       optical_depth_by_slot, neutral_hydrogen_by_slot, absorbed_by_slot, &
       ierr, nleaf, n_interface_face)
    ! optical_depth_by_slot and neutral_hydrogen_by_slot are supplied by the
    ! NLTE chemistry layer in the same n_gamma / n_H,unit convention used by
    ! snrt_intensity.  One photon group is currently enforced explicitly.
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_slot(:)
    real(c_float), intent(in) :: neutral_hydrogen_by_slot(:)
    real(c_float), intent(out) :: absorbed_by_slot(:)
    integer, intent(out) :: ierr, nleaf, n_interface_face
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:), tau(:), neutral(:), absorbed(:)
    real(c_float), allocatable :: tau_groups(:,:)
    integer :: igroup
    integer :: ilocal
    integer(c_int) :: cuda_ierr

    ierr = 0
    nleaf = 0
    n_interface_face = 0
    absorbed_by_slot = 0.0_c_float
    if (cdt_over_dx < 0.0d0) then
       ierr = 1
       return
    end if
    if (size(optical_depth_by_slot) < snrt_nslot .or. &
         size(neutral_hydrogen_by_slot) < snrt_nslot .or. &
         size(absorbed_by_slot) < snrt_nslot) then
       ierr = 2
       return
    end if

    ! Preserve the old scalar-tau entry point after the state becomes
    ! multigroup.  A caller with no group-resolved chemistry gets the same
    ! tau in every group; group-resolved callers should use the routine below.
    if (snrt_ngroups /= 1) then
       allocate(tau_groups(snrt_nslot,snrt_ngroups))
       do igroup = 1, snrt_ngroups
          tau_groups(:,igroup) = optical_depth_by_slot(1:snrt_nslot)
       end do
       call snrt_transport_absorb_multigroup_level(ilevel, cdt_over_dx, &
            tau_groups, neutral_hydrogen_by_slot, absorbed_by_slot, ierr, &
            nleaf, n_interface_face)
       deallocate(tau_groups)
       return
    end if

    call snrt_angular_init(direction_dp, weight)
    angular_cfl = maxval(sum(abs(direction_dp), dim=2))
    if (cdt_over_dx * angular_cfl > 1.0d0) then
       ierr = 4
       return
    end if

    call snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
         neighbor, nleaf, n_interface_face)
    if (nleaf == 0) return

    allocate(neighbor_c(6,nleaf), packed(nleaf,snrt_ndirection), tau(nleaf), &
         neutral(nleaf), absorbed(nleaf))
    neighbor_c = int(neighbor, c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       packed(ilocal,:) = snrt_intensity(:,1,leaf_slot(ilocal))
       tau(ilocal) = optical_depth_by_slot(leaf_slot(ilocal))
       neutral(ilocal) = neutral_hydrogen_by_slot(leaf_slot(ilocal))
    end do

    cuda_ierr = snrt_cuda_transport_absorb_limited(packed, direction_c, &
         neighbor_c, tau, neutral, absorbed, int(nleaf,c_int), &
         int(snrt_ndirection,c_int), real(cdt_over_dx,c_float))
    if (cuda_ierr /= 0_c_int) then
       ierr = 100 + int(cuda_ierr)
       return
    end if
    do ilocal = 1, nleaf
       snrt_intensity(:,1,leaf_slot(ilocal)) = packed(ilocal,:)
       absorbed_by_slot(leaf_slot(ilocal)) = absorbed(ilocal)
    end do
  end subroutine snrt_transport_absorb_limited_level

  subroutine snrt_transport_absorb_multigroup_level(ilevel, cdt_over_dx, &
       optical_depth_by_slot_group, neutral_hydrogen_by_slot, &
       absorbed_by_slot, ierr, nleaf, n_interface_face)
    ! Group-resolved adapter.  The optical-depth array is indexed as
    ! (slot, group), while the CUDA ABI receives (leaf, group) in its
    ! cell-major layout.  Absorption is capped against the neutral-H budget
    ! after summing all directions and groups.
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_slot_group(:,:)
    real(c_float), intent(in) :: neutral_hydrogen_by_slot(:)
    real(c_float), intent(out) :: absorbed_by_slot(:)
    integer, intent(out) :: ierr, nleaf, n_interface_face
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:,:), tau(:,:), neutral(:), absorbed(:)
    real(c_float), allocatable :: absorbed_group(:,:)
    integer :: ilocal, igroup
    integer(c_int) :: cuda_ierr

    ierr = 0
    nleaf = 0
    n_interface_face = 0
    absorbed_by_slot = 0.0_c_float
    if (cdt_over_dx < 0.0d0) then
       ierr = 1
       return
    end if
    if (size(optical_depth_by_slot_group,1) < snrt_nslot .or. &
         size(optical_depth_by_slot_group,2) < snrt_ngroups .or. &
         size(neutral_hydrogen_by_slot) < snrt_nslot .or. &
         size(absorbed_by_slot) < snrt_nslot) then
       ierr = 2
       return
    end if

    call snrt_angular_init(direction_dp, weight)
    angular_cfl = maxval(sum(abs(direction_dp), dim=2))
    if (cdt_over_dx * angular_cfl > 1.0d0) then
       ierr = 3
       return
    end if

    call snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
         neighbor, nleaf, n_interface_face)
    if (nleaf == 0) return

    allocate(neighbor_c(6,nleaf), &
         packed(nleaf,snrt_ndirection,snrt_ngroups), &
         tau(nleaf,snrt_ngroups), neutral(nleaf), absorbed(nleaf), &
         absorbed_group(nleaf,snrt_ngroups))
    neighbor_c = int(neighbor, c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       neutral(ilocal) = neutral_hydrogen_by_slot(leaf_slot(ilocal))
       do igroup = 1, snrt_ngroups
          packed(ilocal,:,igroup) = snrt_intensity(:,igroup,leaf_slot(ilocal))
          tau(ilocal,igroup) = optical_depth_by_slot_group(leaf_slot(ilocal),igroup)
       end do
    end do

    cuda_ierr = snrt_cuda_multigroup_rt_step(packed, direction_c, neighbor_c, &
         tau, neutral, absorbed, absorbed_group, int(nleaf,c_int), &
         int(snrt_ndirection,c_int), &
         int(snrt_ngroups,c_int), real(cdt_over_dx,c_float))
    if (cuda_ierr /= 0_c_int) then
       ierr = 100 + int(cuda_ierr)
       return
    end if
    do ilocal = 1, nleaf
       absorbed_by_slot(leaf_slot(ilocal)) = absorbed(ilocal)
       do igroup = 1, snrt_ngroups
          snrt_intensity(:,igroup,leaf_slot(ilocal)) = packed(ilocal,:,igroup)
       end do
    end do
  end subroutine snrt_transport_absorb_multigroup_level

  subroutine snrt_transport_absorb_multigroup_prepared(leaf_slot, neighbor, &
       cdt_over_dx, optical_depth_by_leaf_group, neutral_hydrogen_by_leaf, &
       absorbed_by_leaf_group, ierr)
    ! Prepared-cell ABI used by the RAMSES driver.  Keeping topology outside
    ! this routine lets the driver construct hydro/NLTE fields without a
    ! second AMR traversal.
    integer, intent(in) :: leaf_slot(:)
    integer, intent(in) :: neighbor(:,:)
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_leaf_group(:,:)
    real(c_float), intent(in) :: neutral_hydrogen_by_leaf(:)
    real(c_float), intent(out) :: absorbed_by_leaf_group(:,:)
    integer, intent(out) :: ierr
    integer :: nleaf, ilocal, igroup, nsub, isub
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:,:), tau(:,:), neutral(:)
    real(c_float), allocatable :: absorbed_total(:), absorbed_group(:,:)
    integer(c_int) :: cuda_ierr

    ierr = 0
    absorbed_by_leaf_group = 0.0_c_float
    nleaf = size(leaf_slot)
    if (cdt_over_dx < 0.0d0) then
       ierr = 1
       return
    end if
    if (size(neighbor,1) < 6 .or. size(neighbor,2) < nleaf .or. &
         size(optical_depth_by_leaf_group,1) < nleaf .or. &
         size(optical_depth_by_leaf_group,2) < snrt_ngroups .or. &
         size(neutral_hydrogen_by_leaf) < nleaf .or. &
         size(absorbed_by_leaf_group,1) < nleaf .or. &
         size(absorbed_by_leaf_group,2) < snrt_ngroups) then
       ierr = 2
       return
    end if
    if (nleaf == 0) return

    call snrt_angular_init(direction_dp, weight)
    angular_cfl = maxval(sum(abs(direction_dp), dim=2))
    nsub = max(1, ceiling(cdt_over_dx * angular_cfl))

    allocate(neighbor_c(6,nleaf), packed(nleaf,snrt_ndirection,snrt_ngroups), &
         tau(nleaf,snrt_ngroups), neutral(nleaf), absorbed_total(nleaf), &
         absorbed_group(nleaf,snrt_ngroups))
    neighbor_c = int(neighbor(:,1:nleaf), c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       neutral(ilocal) = neutral_hydrogen_by_leaf(ilocal)
       do igroup = 1, snrt_ngroups
          packed(ilocal,:,igroup) = snrt_intensity(:,igroup,leaf_slot(ilocal))
          tau(ilocal,igroup) = optical_depth_by_leaf_group(ilocal,igroup) / &
               real(nsub,c_float)
       end do
    end do

    absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) = 0.0_c_float
    do isub = 1, nsub
       absorbed_group = 0.0_c_float
       cuda_ierr = snrt_cuda_multigroup_rt_step(packed, direction_c, neighbor_c, &
            tau, neutral, absorbed_total, absorbed_group, int(nleaf,c_int), &
            int(snrt_ndirection,c_int), int(snrt_ngroups,c_int), &
            real(cdt_over_dx/real(nsub,dp),c_float))
       if (cuda_ierr /= 0_c_int) then
          ierr = 100 + int(cuda_ierr)
          return
       end if
       absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) = &
            absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) + absorbed_group
       ! The CUDA cap is a per-call neutral-H budget.  Carry the remaining
       ! budget across transport substeps so an optically thick cell cannot
       ! consume the same atoms repeatedly within one hydro step.
       neutral = max(0.0_c_float, neutral - absorbed_total)
    end do

    do ilocal = 1, nleaf
       do igroup = 1, snrt_ngroups
          snrt_intensity(:,igroup,leaf_slot(ilocal)) = packed(ilocal,:,igroup)
       end do
    end do
  end subroutine snrt_transport_absorb_multigroup_prepared

end module snrt_transport_step

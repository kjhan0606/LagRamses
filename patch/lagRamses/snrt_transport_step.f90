module snrt_transport_step
  ! AMR leaf-slot <-> GPU sparse-stencil adapter.  Same-level MPI faces use
  ! RAMSES virtual cells as read-only GPU ghosts; coarse-fine faces remain a
  ! separate conservative flux-register responsibility.
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection, snrt_ngroups, snrt_nslot, snrt_intensity
  use snrt_angular_quadrature, only: snrt_angular_init
  use snrt_amr_topology, only: snrt_amr_build_same_level_neighbors, &
       snrt_amr_exchange_interface_state, snrt_amr_apply_coarse_flux_correction, &
       snrt_face_kind, snrt_face_cell, &
       SNRT_FACE_MPI, SNRT_FACE_PHYSICAL, SNRT_FACE_COARSE_TO_FINE, &
       SNRT_FACE_FINE_TO_COARSE, SNRT_FACE_UNMAPPED
  use snrt_cuda_sparse_transport_interface, only: snrt_cuda_upwind_sparse
  use snrt_cuda_limited_rt_step_interface, only: snrt_cuda_transport_absorb_limited
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step, &
       snrt_cuda_multigroup_rt_step_owned, snrt_cuda_multigroup_rt_step_species, &
       snrt_cuda_multigroup_rt_step_species_dust
  use snrt_runtime_backend, only: snrt_runtime_species_dust_step
  implicit none

contains

  subroutine snrt_transport_collective_error(local_error, global_error)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    integer, intent(in) :: local_error
    integer, intent(out) :: global_error
    integer :: info

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_error, global_error, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, info)
    if (info /= 0) global_error = max(global_error, 999)
#else
    global_error = local_error
#endif
  end subroutine snrt_transport_collective_error

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
       optical_depth_by_slot, absorbing_atom_by_slot, absorbed_by_slot, &
       ierr, nleaf, n_interface_face)
    ! optical_depth_by_slot and absorbing_atom_by_slot are supplied by the
    ! NLTE chemistry layer in the same n_gamma / n_H,unit convention used by
    ! snrt_intensity.  The budget is the total number of absorbing H/He
    ! nuclei in code density units, not a hydrogen-only budget.
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_slot(:)
    real(c_float), intent(in) :: absorbing_atom_by_slot(:)
    real(c_float), intent(out) :: absorbed_by_slot(:)
    integer, intent(out) :: ierr, nleaf, n_interface_face
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:), tau(:), atom_budget(:), absorbed(:)
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
         size(absorbing_atom_by_slot) < snrt_nslot .or. &
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
            tau_groups, absorbing_atom_by_slot, absorbed_by_slot, ierr, &
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
         atom_budget(nleaf), absorbed(nleaf))
    neighbor_c = int(neighbor, c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       packed(ilocal,:) = snrt_intensity(:,1,leaf_slot(ilocal))
       tau(ilocal) = optical_depth_by_slot(leaf_slot(ilocal))
       atom_budget(ilocal) = absorbing_atom_by_slot(leaf_slot(ilocal))
    end do

    cuda_ierr = snrt_cuda_transport_absorb_limited(packed, direction_c, &
         neighbor_c, tau, atom_budget, absorbed, int(nleaf,c_int), &
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
       optical_depth_by_slot_group, absorbing_atom_by_slot, &
       absorbed_by_slot, ierr, nleaf, n_interface_face)
    ! Group-resolved adapter.  The optical-depth array is indexed as
    ! (slot, group), while the CUDA ABI receives (leaf, group) in its
    ! cell-major layout.  This is a compatibility adapter for scalar-budget
    ! benchmarks; the production RAMSES driver uses the species-aware
    ! prepared adapter below.
    integer, intent(in) :: ilevel
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_slot_group(:,:)
    real(c_float), intent(in) :: absorbing_atom_by_slot(:)
    real(c_float), intent(out) :: absorbed_by_slot(:)
    integer, intent(out) :: ierr, nleaf, n_interface_face
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:,:), tau(:,:), atom_budget(:), absorbed(:)
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
         size(absorbing_atom_by_slot) < snrt_nslot .or. &
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
         tau(nleaf,snrt_ngroups), atom_budget(nleaf), absorbed(nleaf), &
         absorbed_group(nleaf,snrt_ngroups))
    neighbor_c = int(neighbor, c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       atom_budget(ilocal) = absorbing_atom_by_slot(leaf_slot(ilocal))
       do igroup = 1, snrt_ngroups
          packed(ilocal,:,igroup) = snrt_intensity(:,igroup,leaf_slot(ilocal))
          tau(ilocal,igroup) = optical_depth_by_slot_group(leaf_slot(ilocal),igroup)
       end do
    end do

    cuda_ierr = snrt_cuda_multigroup_rt_step(packed, direction_c, neighbor_c, &
         tau, atom_budget, absorbed, absorbed_group, int(nleaf,c_int), &
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

  subroutine snrt_transport_absorb_multigroup_prepared_trial(leaf_slot, neighbor, &
       cdt_over_dx, optical_depth_by_leaf_group, optical_depth_by_leaf_species, &
       available_species_by_leaf, incoming_intensity, trial_intensity, &
       coarse_flux_trial, absorbed_by_leaf_group, ierr, leaf_cell, ilevel)
    ! Backward-compatible H/He-only prepared ABI.  The DUST-8 entry point is
    ! separate so existing callers and link checks retain their symbol.
    integer, intent(in) :: leaf_slot(:)
    integer, intent(in) :: neighbor(:,:)
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_leaf_group(:,:)
    real(c_float), intent(in) :: optical_depth_by_leaf_species(:,:,:)
    real(c_float), intent(in) :: available_species_by_leaf(:,:)
    real(c_float), intent(in) :: incoming_intensity(:,:,:)
    real(c_float), intent(out) :: trial_intensity(:,:,:)
    real(c_float), intent(out) :: coarse_flux_trial(:,:,:)
    real(c_float), intent(out) :: absorbed_by_leaf_group(:,:)
    integer, intent(out) :: ierr
    integer, intent(in), optional :: leaf_cell(:)
    integer, intent(in), optional :: ilevel

    call snrt_transport_absorb_multigroup_prepared_trial_core(leaf_slot, neighbor, &
         cdt_over_dx, optical_depth_by_leaf_group, optical_depth_by_leaf_species, &
         available_species_by_leaf, incoming_intensity, trial_intensity, &
         coarse_flux_trial, absorbed_by_leaf_group, ierr, .false., leaf_cell, ilevel)
  end subroutine snrt_transport_absorb_multigroup_prepared_trial

  subroutine snrt_transport_absorb_multigroup_prepared_dust_trial(leaf_slot, neighbor, &
       cdt_over_dx, optical_depth_by_leaf_group, optical_depth_by_leaf_species, &
       optical_depth_by_leaf_dust, available_species_by_leaf, incoming_intensity, &
       trial_intensity, coarse_flux_trial, raw_by_leaf_group, &
       absorbed_hhe_by_leaf_group_species, absorbed_dust_by_leaf_group, &
       returned_by_leaf_group, absorbed_by_leaf_group, ierr, leaf_cell, ilevel)
    ! DUST-8 prepared ABI.  H/He species output uses (leaf,group,species),
    ! the same contiguous layout consumed by the CUDA DUST-7 wrapper.
    integer, intent(in) :: leaf_slot(:)
    integer, intent(in) :: neighbor(:,:)
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_leaf_group(:,:)
    real(c_float), intent(in) :: optical_depth_by_leaf_species(:,:,:)
    real(c_float), intent(in) :: optical_depth_by_leaf_dust(:,:)
    real(c_float), intent(in) :: available_species_by_leaf(:,:)
    real(c_float), intent(in) :: incoming_intensity(:,:,:)
    real(c_float), intent(out) :: trial_intensity(:,:,:)
    real(c_float), intent(out) :: coarse_flux_trial(:,:,:)
    real(c_float), intent(out) :: raw_by_leaf_group(:,:)
    real(c_float), intent(out) :: absorbed_hhe_by_leaf_group_species(:,:,:)
    real(c_float), intent(out) :: absorbed_dust_by_leaf_group(:,:)
    real(c_float), intent(out) :: returned_by_leaf_group(:,:)
    real(c_float), intent(out) :: absorbed_by_leaf_group(:,:)
    integer, intent(out) :: ierr
    integer, intent(in), optional :: leaf_cell(:)
    integer, intent(in), optional :: ilevel

    call snrt_transport_absorb_multigroup_prepared_trial_core(leaf_slot, neighbor, &
         cdt_over_dx, optical_depth_by_leaf_group, optical_depth_by_leaf_species, &
         available_species_by_leaf, incoming_intensity, trial_intensity, &
         coarse_flux_trial, absorbed_by_leaf_group, ierr, .true., leaf_cell, ilevel, &
         optical_depth_by_leaf_dust, absorbed_hhe_by_leaf_group_species, &
         absorbed_dust_by_leaf_group, returned_by_leaf_group, raw_by_leaf_group)
  end subroutine snrt_transport_absorb_multigroup_prepared_dust_trial

  subroutine snrt_transport_absorb_multigroup_prepared_trial_core(leaf_slot, neighbor, &
       cdt_over_dx, optical_depth_by_leaf_group, optical_depth_by_leaf_species, &
       available_species_by_leaf, incoming_intensity, trial_intensity, &
       coarse_flux_trial, absorbed_by_leaf_group, ierr, use_dust, leaf_cell, ilevel, &
       optical_depth_by_leaf_dust, absorbed_hhe_by_leaf_group_species, &
       absorbed_dust_by_leaf_group, returned_by_leaf_group, raw_by_leaf_group)
    ! Prepared-cell ABI used by the RAMSES driver.  Keeping topology outside
    ! this routine lets the driver construct hydro/NLTE fields without a
    ! second AMR traversal.
    integer, intent(in) :: leaf_slot(:)
    integer, intent(in) :: neighbor(:,:)
    real(dp), intent(in) :: cdt_over_dx
    real(c_float), intent(in) :: optical_depth_by_leaf_group(:,:)
    ! The third dimension is (H I, He I, He II).  Keeping this array in the
    ! same (leaf, group, species) order as the Fortran caller makes its
    ! contiguous C view (species, group, leaf) explicit in the CUDA kernel.
    real(c_float), intent(in) :: optical_depth_by_leaf_species(:,:,:)
    ! The species inventory is copied into a working array and consumed by
    ! CUDA in group order.  It is therefore shared across groups and all
    ! transport substeps, not reset once per group.
    real(c_float), intent(in) :: available_species_by_leaf(:,:)
    ! Incoming and returned local-leaf fields are explicit trial arrays.  No
    ! persistent photon state is changed by this routine.
    real(c_float), intent(in) :: incoming_intensity(:,:,:)
    real(c_float), intent(out) :: trial_intensity(:,:,:)
    ! Coarse/interface corrections are returned as a full-slot trial buffer;
    ! snrt_rt_transaction commits them only after the coupled level succeeds.
    real(c_float), intent(out) :: coarse_flux_trial(:,:,:)
    real(c_float), intent(out) :: absorbed_by_leaf_group(:,:)
    integer, intent(out) :: ierr
    logical, intent(in) :: use_dust
    integer, intent(in), optional :: leaf_cell(:)
    integer, intent(in), optional :: ilevel
    real(c_float), intent(in), optional :: optical_depth_by_leaf_dust(:,:)
    real(c_float), intent(out), optional :: absorbed_hhe_by_leaf_group_species(:,:,:)
    real(c_float), intent(out), optional :: absorbed_dust_by_leaf_group(:,:)
    real(c_float), intent(out), optional :: returned_by_leaf_group(:,:)
    real(c_float), intent(out), optional :: raw_by_leaf_group(:,:)
    integer :: nleaf, ilocal, igroup, nsub, isub
    integer :: nmpi, nwork, iwork, iface, ighost
    integer :: face_kind
    integer :: local_error, global_error
    integer, allocatable :: neighbor_work(:,:), ghost_kind(:), ghost_cell(:), &
         ghost_face(:), ghost_local(:)
    integer(c_int), allocatable :: neighbor_c(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
    real(dp) :: angular_cfl
    real(c_float) :: direction_c(3,snrt_ndirection)
    real(c_float), allocatable :: packed(:,:,:), packed_work(:,:,:)
    real(c_float), allocatable :: ghost_state(:,:,:), tau(:,:), species_tau(:,:,:)
    real(c_float), allocatable :: dust_tau(:,:)
    real(c_float), allocatable :: species_budget(:,:)
    real(c_float), allocatable :: absorbed_total(:), absorbed_group(:,:)
    real(c_float), allocatable :: absorbed_hhe_group(:,:,:), absorbed_dust_group(:,:), &
         returned_group(:,:), raw_group(:,:)
    integer(c_int) :: cuda_ierr

    ierr = 0
    absorbed_by_leaf_group = 0.0_c_float
    nleaf = size(leaf_slot)
    trial_intensity = 0.0_c_float
    coarse_flux_trial = 0.0_c_float
    if (present(absorbed_hhe_by_leaf_group_species)) &
         absorbed_hhe_by_leaf_group_species = 0.0_c_float
    if (present(absorbed_dust_by_leaf_group)) &
         absorbed_dust_by_leaf_group = 0.0_c_float
    if (present(returned_by_leaf_group)) returned_by_leaf_group = 0.0_c_float
    if (present(raw_by_leaf_group)) raw_by_leaf_group = 0.0_c_float
    local_error = 0
    if (cdt_over_dx < 0.0d0) local_error = 1
    if (size(neighbor,1) < 6 .or. size(neighbor,2) < nleaf .or. &
         size(optical_depth_by_leaf_group,1) < nleaf .or. &
         size(optical_depth_by_leaf_group,2) < snrt_ngroups .or. &
         size(optical_depth_by_leaf_species,1) < nleaf .or. &
         size(optical_depth_by_leaf_species,2) < snrt_ngroups .or. &
         size(optical_depth_by_leaf_species,3) < 3 .or. &
         size(available_species_by_leaf,1) < nleaf .or. &
         size(available_species_by_leaf,2) < 3 .or. &
         size(incoming_intensity,1) < snrt_ndirection .or. &
         size(incoming_intensity,2) < snrt_ngroups .or. &
         size(incoming_intensity,3) < nleaf .or. &
         size(trial_intensity,1) < snrt_ndirection .or. &
         size(trial_intensity,2) < snrt_ngroups .or. &
         size(trial_intensity,3) < nleaf .or. &
         size(coarse_flux_trial,1) < snrt_ndirection .or. &
         size(coarse_flux_trial,2) < snrt_ngroups .or. &
         size(coarse_flux_trial,3) < size(snrt_intensity,3) .or. &
         size(absorbed_by_leaf_group,1) < nleaf .or. &
         size(absorbed_by_leaf_group,2) < snrt_ngroups) then
       local_error = max(local_error,2)
    end if
    if (use_dust) then
       if (.not. present(optical_depth_by_leaf_dust) .or. &
            .not. present(absorbed_hhe_by_leaf_group_species) .or. &
            .not. present(absorbed_dust_by_leaf_group) .or. &
            .not. present(returned_by_leaf_group) .or. .not. present(raw_by_leaf_group)) then
          local_error = max(local_error,2)
       else if (size(optical_depth_by_leaf_dust,1) < nleaf .or. &
            size(optical_depth_by_leaf_dust,2) < snrt_ngroups .or. &
            size(absorbed_hhe_by_leaf_group_species,1) < nleaf .or. &
            size(absorbed_hhe_by_leaf_group_species,2) < snrt_ngroups .or. &
            size(absorbed_hhe_by_leaf_group_species,3) < 3 .or. &
            size(absorbed_dust_by_leaf_group,1) < nleaf .or. &
            size(absorbed_dust_by_leaf_group,2) < snrt_ngroups .or. &
            size(returned_by_leaf_group,1) < nleaf .or. &
            size(returned_by_leaf_group,2) < snrt_ngroups .or. &
            size(raw_by_leaf_group,1) < nleaf .or. &
            size(raw_by_leaf_group,2) < snrt_ngroups) then
          local_error = max(local_error,2)
       end if
    end if
    if (nleaf > 0) then
       if (.not. allocated(snrt_face_kind) .or. .not. allocated(snrt_face_cell)) then
          local_error = max(local_error,3)
       else if (size(snrt_face_kind,1) < size(neighbor,1) .or. &
            size(snrt_face_kind,2) < nleaf .or. size(snrt_face_cell,1) < size(neighbor,1) .or. &
            size(snrt_face_cell,2) < nleaf) then
          local_error = max(local_error,3)
       end if
    end if

    nmpi = 0
    if (local_error == 0 .and. nleaf > 0) then
       nmpi = count(snrt_face_kind(1:size(neighbor,1),1:nleaf) == SNRT_FACE_MPI) + &
            count(snrt_face_kind(1:size(neighbor,1),1:nleaf) == SNRT_FACE_FINE_TO_COARSE)
       do ilocal = 1, nleaf
          do iface = 1, 6
             face_kind = snrt_face_kind(iface,ilocal)
             if (face_kind == SNRT_FACE_UNMAPPED) &
                  local_error = max(local_error,3)
          end do
       end do
       if (nmpi > 0 .and. (.not. present(leaf_cell) .or. .not. present(ilevel))) then
          local_error = max(local_error,4)
       end if
    end if
    call snrt_transport_collective_error(local_error, global_error)
    if (global_error /= 0) then
       ierr = global_error
       return
    end if
    nwork = nleaf + nmpi
    allocate(neighbor_work(6,nleaf), packed_work(nwork,snrt_ndirection,snrt_ngroups))
    neighbor_work = neighbor(:,1:nleaf)
    allocate(ghost_kind(nmpi), ghost_cell(nmpi), ghost_face(nmpi), ghost_local(nmpi), &
         ghost_state(nmpi,snrt_ndirection,snrt_ngroups))
    if (nmpi > 0) then
       iwork = nleaf
       ighost = 0
       do ilocal = 1, nleaf
          do iface = 1, 6
             face_kind = snrt_face_kind(iface,ilocal)
             if (face_kind == SNRT_FACE_MPI .or. &
                  face_kind == SNRT_FACE_FINE_TO_COARSE) then
                ighost = ighost + 1
                iwork = iwork + 1
                ghost_kind(ighost) = face_kind
                ghost_cell(ighost) = snrt_face_cell(iface,ilocal)
                ghost_face(ighost) = iface
                ghost_local(ighost) = ilocal
                neighbor_work(iface,ilocal) = iwork
             else if (face_kind == SNRT_FACE_COARSE_TO_FINE .or. &
                  face_kind == SNRT_FACE_PHYSICAL) then
                neighbor_work(iface,ilocal) = ilocal
             end if
          end do
       end do
    else
       do ilocal = 1, nleaf
          do iface = 1, 6
             if (snrt_face_kind(iface,ilocal) == SNRT_FACE_PHYSICAL) &
                  neighbor_work(iface,ilocal) = ilocal
             if (snrt_face_kind(iface,ilocal) == SNRT_FACE_COARSE_TO_FINE) &
                  neighbor_work(iface,ilocal) = ilocal
          end do
       end do
    end if

    call snrt_angular_init(direction_dp, weight)
    angular_cfl = maxval(sum(abs(direction_dp), dim=2))
    nsub = max(1, ceiling(cdt_over_dx * angular_cfl))

    allocate(neighbor_c(6,nleaf), packed(nleaf,snrt_ndirection,snrt_ngroups), &
         tau(nleaf,snrt_ngroups), species_tau(nleaf,snrt_ngroups,3), &
         species_budget(nleaf,3), absorbed_total(nleaf), &
         absorbed_group(nleaf,snrt_ngroups))
    if (use_dust) then
       allocate(dust_tau(nleaf,snrt_ngroups), &
            absorbed_hhe_group(nleaf,snrt_ngroups,3), &
            absorbed_dust_group(nleaf,snrt_ngroups), returned_group(nleaf,snrt_ngroups), &
            raw_group(nleaf,snrt_ngroups))
    end if
    neighbor_c = int(neighbor_work, c_int)
    direction_c = real(transpose(direction_dp), c_float)
    do ilocal = 1, nleaf
       species_budget(ilocal,1:3) = available_species_by_leaf(ilocal,1:3)
       do igroup = 1, snrt_ngroups
          packed(ilocal,:,igroup) = incoming_intensity(:,igroup,ilocal)
          tau(ilocal,igroup) = optical_depth_by_leaf_group(ilocal,igroup) / &
               real(nsub,c_float)
          species_tau(ilocal,igroup,1:3) = &
               optical_depth_by_leaf_species(ilocal,igroup,1:3) / real(nsub,c_float)
          if (use_dust) dust_tau(ilocal,igroup) = &
               optical_depth_by_leaf_dust(ilocal,igroup) / real(nsub,c_float)
       end do
    end do

    absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) = 0.0_c_float
    do isub = 1, nsub
       ! The AMR exchange is collective even when this rank has no local
       ! ghost faces; the topology routine reduces interface requirements
       ! across all ranks before entering make_virtual_fine_dp.
       call snrt_amr_exchange_interface_state(ilevel, leaf_cell, packed, &
            ghost_kind, ghost_cell, ghost_face, ghost_state, ierr)
       local_error = 0
       if (ierr /= 0) local_error = 10 + ierr
       call snrt_transport_collective_error(local_error, global_error)
       if (global_error /= 0) then
          ierr = global_error
          return
       end if
       packed_work = 0.0_c_float
       packed_work(1:nleaf,1:snrt_ndirection,1:snrt_ngroups) = packed
       do ighost = 1, nmpi
          packed_work(nleaf+ighost,1:snrt_ndirection,1:snrt_ngroups) = &
               ghost_state(ighost,1:snrt_ndirection,1:snrt_ngroups)
       end do
       call snrt_amr_apply_coarse_flux_correction(ilevel, leaf_cell, packed_work, &
            ghost_kind, ghost_cell, ghost_face, ghost_local, &
            cdt_over_dx/real(nsub,dp), direction_dp, coarse_flux_trial, ierr)
       local_error = 0
       if (ierr /= 0) local_error = 20 + ierr
       call snrt_transport_collective_error(local_error, global_error)
       if (global_error /= 0) then
          ierr = global_error
          return
       end if
       absorbed_group = 0.0_c_float
       if (use_dust) then
          absorbed_hhe_group = 0.0_c_float
          absorbed_dust_group = 0.0_c_float
          returned_group = 0.0_c_float
          raw_group = 0.0_c_float
       end if
       cuda_ierr = 0_c_int
       if (nleaf > 0) then
          if (use_dust) then
             cuda_ierr = snrt_runtime_species_dust_step(packed_work, direction_c, &
                  neighbor_c, tau, species_tau, dust_tau, species_budget, &
                  absorbed_hhe_group, absorbed_dust_group, returned_group, raw_group, &
                  absorbed_group, absorbed_total, int(nleaf,c_int), int(nwork,c_int), &
                  int(snrt_ndirection,c_int), int(snrt_ngroups,c_int), &
                  real(cdt_over_dx/real(nsub,dp),c_float))
          else
             cuda_ierr = snrt_cuda_multigroup_rt_step_species(packed_work, direction_c, &
                  neighbor_c, tau, species_tau, species_budget, absorbed_total, &
                  absorbed_group, int(nleaf,c_int), int(nwork,c_int), &
                  int(snrt_ndirection,c_int), int(snrt_ngroups,c_int), &
                  real(cdt_over_dx/real(nsub,dp),c_float))
          end if
       end if
       local_error = 0
       if (cuda_ierr /= 0_c_int) local_error = 100 + int(cuda_ierr)
       call snrt_transport_collective_error(local_error, global_error)
       if (global_error /= 0) then
          ierr = global_error
          return
       end if
       absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) = &
            absorbed_by_leaf_group(1:nleaf,1:snrt_ngroups) + absorbed_group
       if (use_dust) then
          absorbed_hhe_by_leaf_group_species(1:nleaf,1:snrt_ngroups,1:3) = &
               absorbed_hhe_by_leaf_group_species(1:nleaf,1:snrt_ngroups,1:3) + &
               absorbed_hhe_group
          absorbed_dust_by_leaf_group(1:nleaf,1:snrt_ngroups) = &
               absorbed_dust_by_leaf_group(1:nleaf,1:snrt_ngroups) + absorbed_dust_group
          returned_by_leaf_group(1:nleaf,1:snrt_ngroups) = &
               returned_by_leaf_group(1:nleaf,1:snrt_ngroups) + returned_group
          raw_by_leaf_group(1:nleaf,1:snrt_ngroups) = &
               raw_by_leaf_group(1:nleaf,1:snrt_ngroups) + raw_group
       end if
       packed = packed_work(1:nleaf,1:snrt_ndirection,1:snrt_ngroups)
    end do

    do ilocal = 1, nleaf
       do igroup = 1, snrt_ngroups
          trial_intensity(:,igroup,ilocal) = packed(ilocal,:,igroup)
       end do
    end do
    deallocate(neighbor_work, packed_work, ghost_kind, ghost_cell, ghost_face, &
         ghost_local, ghost_state, neighbor_c, packed, tau, species_tau, &
         species_budget, absorbed_total, absorbed_group)
    if (use_dust) deallocate(dust_tau, absorbed_hhe_group, absorbed_dust_group, &
         returned_group, raw_group)
  end subroutine snrt_transport_absorb_multigroup_prepared_trial_core

end module snrt_transport_step

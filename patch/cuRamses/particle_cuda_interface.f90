!==========================================================================
! ISO_C_BINDING interface to particle_cuda_kernels.cu — GPU CIC force
! gather + kick/drift for move_fine/synchro_fine (strategy-A hybrid).
!==========================================================================
module particle_cuda_interface
  use iso_c_binding
  implicit none

  integer(c_int), parameter :: PM_MODE_MOVE = 0
  integer(c_int), parameter :: PM_MODE_SYNC = 1
  integer(c_int), parameter :: PM_MODE_RHO  = 2

  interface

     subroutine cuda_pm_mesh_upload_c(f, son, phi, ncell, with_phi, &
          & ncoarse_c, ngridmax_c, hw, block_size, child_count) &
          & bind(C, name='cuda_pm_mesh_upload')
       import :: c_double, c_int, c_long_long
       real(c_double), dimension(*), intent(in) :: f, phi
       integer(c_int), dimension(*), intent(in) :: son
       integer(c_long_long), value :: ncell, ncoarse_c
       integer(c_int), value :: with_phi, ngridmax_c, hw
       integer(c_int), value :: block_size, child_count
     end subroutine cuda_pm_mesh_upload_c

     function cuda_pm_is_ready_c() result(ready) &
          & bind(C, name='cuda_pm_is_ready')
       import :: c_int
       integer(c_int) :: ready
     end function cuda_pm_is_ready_c

     function cuda_pm_flush_c(slot, ng, np, x0, nbf, px, pv, pg, dteff, &
          & params, ngridmax_c, ncoarse_c, block_size, child_count, &
          & new_v, new_x, phi_out) &
          & result(ierr) bind(C, name='cuda_pm_flush')
       import :: c_double, c_int
       integer(c_int), value :: slot, ng, np, ngridmax_c, ncoarse_c
       integer(c_int), value :: block_size, child_count
       real(c_double), dimension(*), intent(in) :: x0, px, pv, dteff, params
       integer(c_int), dimension(*), intent(in) :: nbf, pg
       real(c_double), dimension(*), intent(out) :: new_v, new_x, phi_out
       integer(c_int) :: ierr
     end function cuda_pm_flush_c

     subroutine cuda_pm_rho_begin_c(son, ncell, ncoarse_c, ngridmax_c, hw, &
          & block_size, child_count) &
          & bind(C, name='cuda_pm_rho_begin')
       import :: c_int, c_long_long
       integer(c_int), dimension(*), intent(in) :: son
       integer(c_long_long), value :: ncell, ncoarse_c
       integer(c_int), value :: ngridmax_c, hw, block_size, child_count
     end subroutine cuda_pm_rho_begin_c

     function cuda_pm_rho_is_ready_c() result(ready) &
          & bind(C, name='cuda_pm_rho_is_ready')
       import :: c_int
       integer(c_int) :: ready
     end function cuda_pm_rho_is_ready_c

     function cuda_pm_deposit_flush_c(slot, ng, np, x0, nbf, px, mass, pg, &
          & params, ngridmax_c, ncoarse_c, block_size, child_count) &
          & result(ierr) &
          & bind(C, name='cuda_pm_deposit_flush')
       import :: c_double, c_int
       integer(c_int), value :: slot, ng, np, ngridmax_c, ncoarse_c
       integer(c_int), value :: block_size, child_count
       real(c_double), dimension(*), intent(in) :: x0, px, mass, params
       integer(c_int), dimension(*), intent(in) :: nbf, pg
       integer(c_int) :: ierr
     end function cuda_pm_deposit_flush_c

     subroutine cuda_pm_rho_end_c(rho_add, phiw_add, ncell) &
          & bind(C, name='cuda_pm_rho_end')
       import :: c_double, c_long_long
       real(c_double), dimension(*), intent(out) :: rho_add, phiw_add
       integer(c_long_long), value :: ncell
     end subroutine cuda_pm_rho_end_c

     subroutine cuda_pm_report_c() bind(C, name='cuda_pm_report')
     end subroutine cuda_pm_report_c

     subroutine cuda_pm_finalize_c() bind(C, name='cuda_pm_finalize')
     end subroutine cuda_pm_finalize_c

  end interface

end module particle_cuda_interface

!==========================================================================
! Per-stream-slot superbatch staging for the particle GPU path.
! Each GPU-worker thread owns one slot; batches assembled by the grid
! loop are appended here and flushed to the GPU in large chunks.
! The original (ind_grid, ind_part) batch records are kept so a CUDA
! failure can replay every appended batch through the CPU routine.
!==========================================================================
module pm_gpu_commons
  use amr_parameters, only: dp
  implicit none

  integer, parameter :: PM_SUPER   = 16384  ! particles per flush
  integer, parameter :: PM_SUPER_G = 16384  ! grids per flush (worst case)
  integer, parameter :: PM_MAX_SLOT = 16    ! matches MAX_CUDA_STREAMS

  logical :: pm_gpu_inited = .false.
  logical :: pm_gpu_dead   = .false.  ! set on CUDA failure: no new GPU batches

  ! Packed superbatch (per slot)
  real(dp), allocatable :: pmg_x0(:,:)     ! (3*PM_SUPER_G, 0:PM_MAX_SLOT-1)
  integer,  allocatable :: pmg_nbf(:,:)    ! (27*PM_SUPER_G, slot)
  real(dp), allocatable :: pmg_px(:,:)     ! (3*PM_SUPER, slot)
  real(dp), allocatable :: pmg_pv(:,:)     ! (3*PM_SUPER, slot)
  real(dp), allocatable :: pmg_dt(:,:)     ! (PM_SUPER, slot) dteff
  integer,  allocatable :: pmg_pg(:,:)     ! (PM_SUPER, slot) 0-based grid idx
  integer,  allocatable :: pmg_idx(:,:)    ! (PM_SUPER, slot) ind_part
  integer,  allocatable :: pmg_lvl(:,:)    ! (PM_SUPER, slot) old levelp (sync)
  ! Flush outputs
  real(dp), allocatable :: pmg_nv(:,:)     ! (3*PM_SUPER, slot)
  real(dp), allocatable :: pmg_nx(:,:)     ! (3*PM_SUPER, slot)
  real(dp), allocatable :: pmg_pho(:,:)    ! (PM_SUPER, slot)
  ! Replay records (CPU fallback on CUDA failure)
  integer,  allocatable :: pmg_rgrid(:,:)  ! (PM_SUPER_G, slot) AMR igrid
  integer,  allocatable :: pmg_nb(:)       ! (slot) number of batches
  integer,  allocatable :: pmg_bg0(:,:)    ! (PM_SUPER_G+1, slot) batch grid offsets
  integer,  allocatable :: pmg_bp0(:,:)    ! (PM_SUPER+1, slot) batch particle offsets
  ! Counters
  integer,  allocatable :: pmg_ng(:)       ! (slot) grids appended
  integer,  allocatable :: pmg_np(:)       ! (slot) particles appended
  ! Host merge buffers for the GPU deposit (rho mode)
  real(dp), allocatable :: pmg_rho_add(:), pmg_phiw_add(:)

contains

  subroutine pm_gpu_alloc()
    if (pm_gpu_inited) return
    allocate(pmg_x0 (3*PM_SUPER_G, 0:PM_MAX_SLOT-1))
    allocate(pmg_nbf(27*PM_SUPER_G, 0:PM_MAX_SLOT-1))
    allocate(pmg_px (3*PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_pv (3*PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_pg (PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_dt (PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_idx(PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_lvl(PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_nv (3*PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_nx (3*PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_pho(PM_SUPER, 0:PM_MAX_SLOT-1))
    allocate(pmg_rgrid(PM_SUPER_G, 0:PM_MAX_SLOT-1))
    allocate(pmg_nb (0:PM_MAX_SLOT-1))
    allocate(pmg_bg0(PM_SUPER_G+1, 0:PM_MAX_SLOT-1))
    allocate(pmg_bp0(PM_SUPER+1, 0:PM_MAX_SLOT-1))
    allocate(pmg_ng (0:PM_MAX_SLOT-1))
    allocate(pmg_np (0:PM_MAX_SLOT-1))
    pmg_ng = 0
    pmg_np = 0
    pmg_nb = 0
    pm_gpu_inited = .true.
  end subroutine pm_gpu_alloc

  subroutine pm_gpu_reset(slot)
    integer, intent(in) :: slot
    pmg_ng(slot) = 0
    pmg_np(slot) = 0
    pmg_nb(slot) = 0
  end subroutine pm_gpu_reset

end module pm_gpu_commons

!==========================================================================
! Superbatch append/flush shared by move_fine and synchro_fine.
! append: CPU does the tree work (get3cubefather, x0) and packs one
!         nvector batch; flush: one GPU round trip + host writeback.
! A CUDA failure marks the GPU path dead and replays every appended
! batch through the unchanged CPU routine (move1/sync).
!==========================================================================
module pm_gpu_dispatch
  use amr_parameters, only: dp
  implicit none
  integer :: pm_with_phi = 0   ! set at mesh upload (move only)

contains

  !------------------------------------------------------------------
  ! Highest grid index in use, over every level and every cpu (plus
  ! non-periodic boundary grids). ngridmax is only an allocation
  ! ceiling, so the mesh upload copies [1..hw] of each oct slot
  ! instead of the whole array. Cost is one pass over the grid lists.
  !------------------------------------------------------------------
  integer function pm_grid_high_water() result(hw)
    use amr_commons
    implicit none
    integer :: ilev, icpu, ibnd, i, igrid

    hw = 0
    do ilev = 1, nlevelmax
       do icpu = 1, ncpu
          igrid = headl(icpu,ilev)
          do i = 1, numbl(icpu,ilev)
             if (igrid <= 0) exit
             if (igrid > hw) hw = igrid
             igrid = next(igrid)
          end do
       end do
    end do
    do ibnd = 1, nboundary
       do ilev = 1, nlevelmax
          igrid = headb(ibnd,ilev)
          do i = 1, numbb(ibnd,ilev)
             if (igrid <= 0) exit
             if (igrid > hw) hw = igrid
             igrid = next(igrid)
          end do
       end do
    end do
  end function pm_grid_high_water

  !------------------------------------------------------------------
  ! Particles on this level, used to decide whether the GPU can earn
  ! back its per-call mesh upload.
  !------------------------------------------------------------------
  integer function pm_level_npart(ilevel) result(np)
    use amr_commons
    use pm_commons
    implicit none
    integer, intent(in) :: ilevel
    integer :: jgrid, igrid

    np = 0
    igrid = headl(myid,ilevel)
    do jgrid = 1, numbl(myid,ilevel)
       if (igrid <= 0) exit
       np = np + numbp(igrid)
       igrid = next(igrid)
    end do
  end function pm_level_npart

  subroutine pm_gpu_append(slot, mode, ind_grid, ind_part, ind_grid_part, &
       & ng, np, ilevel)
    use amr_commons
    use pm_commons
    use pm_gpu_commons
    use poisson_commons, only: multipole
    use particle_cuda_interface, only: PM_MODE_MOVE, PM_MODE_SYNC, PM_MODE_RHO
    implicit none
    integer, intent(in) :: slot, mode, ng, np, ilevel
    integer, dimension(1:nvector), intent(in) :: ind_grid, ind_part, ind_grid_part

    integer, dimension(1:nvector) :: father_cell
    integer, dimension(1:nvector,1:threetondim) :: nbors_father_cells
    integer, dimension(1:nvector,1:twotondim)  :: nbors_father_grids
    real(dp) :: dx
    real(dp) :: mm_p, mx_p, my_p, mz_p
    integer :: i, j, k, base_g, base_p, nb, lvl

    if (pmg_np(slot)+np > PM_SUPER .or. pmg_ng(slot)+ng > PM_SUPER_G) then
       call pm_gpu_flush(slot, mode, ilevel)
    end if

    dx = 0.5D0**ilevel
    do i = 1, ng
       father_cell(i) = father(ind_grid(i))
    end do
    call get3cubefather(father_cell, nbors_father_cells, &
         & nbors_father_grids, ng, ilevel)

    base_g = pmg_ng(slot)
    do i = 1, ng
       pmg_rgrid(base_g+i, slot) = ind_grid(i)
       do k = 1, ndim
          pmg_x0(3*(base_g+i-1)+k, slot) = xg(ind_grid(i),k) - 3.0D0*dx
       end do
       do k = 1, threetondim
          pmg_nbf(27*(base_g+i-1)+k, slot) = nbors_father_cells(i,k)
       end do
    end do

    base_p = pmg_np(slot)
    mm_p = 0d0; mx_p = 0d0; my_p = 0d0; mz_p = 0d0
    do j = 1, np
       do k = 1, ndim
          pmg_px(3*(base_p+j-1)+k, slot) = xp(ind_part(j),k)
       end do
       pmg_pg(base_p+j, slot) = base_g + ind_grid_part(j) - 1
       pmg_idx(base_p+j, slot) = ind_part(j)
       if (mode == PM_MODE_SYNC) then
          do k = 1, ndim
             pmg_pv(3*(base_p+j-1)+k, slot) = vp(ind_part(j),k)
          end do
          lvl = levelp(ind_part(j))
          pmg_lvl(base_p+j, slot) = lvl
          if (lvl >= ilevel) then
             pmg_dt(base_p+j, slot) = dtnew(lvl)
          else
             pmg_dt(base_p+j, slot) = dtold(lvl)
          end if
          levelp(ind_part(j)) = ilevel
       else if (mode == PM_MODE_RHO) then
          pmg_dt(base_p+j, slot) = mp(ind_part(j))
          if (ilevel == levelmin) then
             mm_p = mm_p + mp(ind_part(j))
             mx_p = mx_p + mp(ind_part(j))*xp(ind_part(j),1)
             my_p = my_p + mp(ind_part(j))*xp(ind_part(j),2)
             mz_p = mz_p + mp(ind_part(j))*xp(ind_part(j),3)
          end if
       else
          do k = 1, ndim
             pmg_pv(3*(base_p+j-1)+k, slot) = vp(ind_part(j),k)
          end do
          pmg_dt(base_p+j, slot) = dtnew(ilevel)
       end if
    end do

    ! Multipole accumulation (cic_amr does this inside its own batch)
    if (mode == PM_MODE_RHO .and. ilevel == levelmin) then
!$omp atomic update
       multipole(1) = multipole(1) + mm_p
!$omp atomic update
       multipole(2) = multipole(2) + mx_p
!$omp atomic update
       multipole(3) = multipole(3) + my_p
!$omp atomic update
       multipole(4) = multipole(4) + mz_p
    end if

    nb = pmg_nb(slot)
    pmg_bg0(nb+1, slot) = base_g
    pmg_bp0(nb+1, slot) = base_p
    pmg_nb(slot) = nb + 1
    pmg_ng(slot) = base_g + ng
    pmg_np(slot) = base_p + np

  end subroutine pm_gpu_append

  subroutine pm_gpu_flush(slot, mode, ilevel)
    use amr_commons
    use pm_commons
    use pm_gpu_commons
    use particle_cuda_interface
    use scalar_de_commons, only: sde_phip_of_a
    use iso_c_binding
    implicit none
    integer, intent(in) :: slot, mode, ilevel

    real(c_double), dimension(1:10) :: params
    integer(c_int) :: ierr
    integer :: j, k, p, b, nb, g0, g1, p0, p1, ng_b, np_b
    integer, dimension(1:nvector) :: r_grid, r_part, r_gpart
    external :: move1, sync

    if (pmg_np(slot) == 0) return

    if (mode == PM_MODE_RHO) then
       call pm_gpu_flush_rho(slot, ilevel)
       return
    end if

    params = 0.0d0
    params(1) = boxlen/dble(icoarse_max-icoarse_min+1)
    params(2) = 0.5D0**ilevel
    params(3) = dble(icoarse_min)
    params(4) = dble(jcoarse_min)
    params(5) = dble(kcoarse_min)
    if (use_coupled_de .and. cde_friction .and. use_quintessence .and. cosmo) then
       params(6) = beta_cde*sde_phip_of_a(aexp)*hexp*0.5D0
    end if
    params(7) = dtnew(ilevel)
    if (static) params(8) = 1.0d0
    if (mode == PM_MODE_MOVE) params(9) = dble(pm_with_phi)
    params(10) = dble(mode)

    ierr = cuda_pm_flush_c(int(slot,c_int), int(pmg_ng(slot),c_int), &
         & int(pmg_np(slot),c_int), &
         & pmg_x0(:,slot), pmg_nbf(:,slot), pmg_px(:,slot), pmg_pv(:,slot), &
         & pmg_pg(:,slot), pmg_dt(:,slot), params, &
         & int(ngridmax,c_int), int(ncoarse,c_int), &
         & int(amr_block_size,c_int), int(twotondim,c_int), &
         & pmg_nv(:,slot), pmg_nx(:,slot), pmg_pho(:,slot))

    if (ierr == 0) then
       do j = 1, pmg_np(slot)
          p = pmg_idx(j, slot)
          do k = 1, ndim
             vp(p,k) = pmg_nv(3*(j-1)+k, slot)
          end do
          if (mode == PM_MODE_MOVE) then
             do k = 1, ndim
                xp(p,k) = pmg_nx(3*(j-1)+k, slot)
             end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
             if (pm_with_phi == 1) ptcl_phi(p) = pmg_pho(j, slot)
#endif
          end if
       end do
    else if (ierr == 1) then
       write(*,*) 'problem in move/sync (gpu path): particle outside 0.5..5.5'
       stop
    else
       ! CUDA failure: disable the GPU path and replay on the CPU
       pm_gpu_dead = .true.
       write(*,*) 'WARNING: particle GPU flush failed, replaying batches on CPU'
       if (mode == PM_MODE_SYNC) then
          do j = 1, pmg_np(slot)
             levelp(pmg_idx(j,slot)) = pmg_lvl(j,slot)
          end do
       end if
       nb = pmg_nb(slot)
       do b = 1, nb
          g0 = pmg_bg0(b, slot)
          p0 = pmg_bp0(b, slot)
          if (b < nb) then
             g1 = pmg_bg0(b+1, slot)
             p1 = pmg_bp0(b+1, slot)
          else
             g1 = pmg_ng(slot)
             p1 = pmg_np(slot)
          end if
          ng_b = g1 - g0
          np_b = p1 - p0
          do j = 1, ng_b
             r_grid(j) = pmg_rgrid(g0+j, slot)
          end do
          do j = 1, np_b
             r_part(j) = pmg_idx(p0+j, slot)
             r_gpart(j) = pmg_pg(p0+j, slot) - g0 + 1
          end do
          if (mode == PM_MODE_MOVE) then
             call move1(r_grid, r_part, r_gpart, ng_b, np_b, ilevel)
          else
             call sync(r_grid, r_part, r_gpart, ng_b, np_b, ilevel)
          end if
       end do
    end if

    call pm_gpu_reset(slot)

  end subroutine pm_gpu_flush

  !------------------------------------------------------------------
  ! Deposit flush (rho mode): one GPU round trip into the device
  ! rho/phiw accumulators; no host writeback here (pm_rho_merge adds
  ! the accumulated deposits after the dispatch loop).
  !------------------------------------------------------------------
  subroutine pm_gpu_flush_rho(slot, ilevel)
    use amr_commons
    use pm_commons
    use pm_gpu_commons
    use poisson_commons, only: multipole
    use particle_cuda_interface
    use iso_c_binding
    implicit none
    integer, intent(in) :: slot, ilevel

    real(c_double), dimension(1:10) :: params
    integer(c_int) :: ierr
    integer :: i, j, b, nb, g0, g1, p0, p1, ng_b, np_b, p
    integer, dimension(1:nvector) :: r_grid, r_part, r_gpart, r_cell
    real(dp), dimension(1:nvector,1:ndim) :: r_x0
    real(dp) :: dx, mm_p, mx_p, my_p, mz_p
    external :: cic_amr

    if (pmg_np(slot) == 0) return

    params = 0.0d0
    params(1) = boxlen/dble(icoarse_max-icoarse_min+1)
    params(2) = 0.5D0**ilevel
    params(3) = dble(icoarse_min)
    params(4) = dble(jcoarse_min)
    params(5) = dble(kcoarse_min)
    params(6) = (params(2)*params(1))**ndim   ! vol_loc
    if (static) params(7) = 1.0d0
    params(8) = mass_cut_refine

    ierr = cuda_pm_deposit_flush_c(int(slot,c_int), int(pmg_ng(slot),c_int), &
         & int(pmg_np(slot),c_int), &
         & pmg_x0(:,slot), pmg_nbf(:,slot), pmg_px(:,slot), pmg_dt(:,slot), &
         & pmg_pg(:,slot), params, int(ngridmax,c_int), int(ncoarse,c_int), &
         & int(amr_block_size,c_int), int(twotondim,c_int))

    if (ierr == 1) then
       write(*,*) 'problem in cic (gpu path): particle outside 0.5..5.5'
       stop
    else if (ierr < 0) then
       ! CUDA failure: disable the GPU path and replay on the CPU.
       ! (Failures happen at allocation/H2D, before any deposit lands.)
       pm_gpu_dead = .true.
       write(*,*) 'WARNING: particle GPU deposit failed, replaying batches on CPU'
       dx = 0.5D0**ilevel
       nb = pmg_nb(slot)
       do b = 1, nb
          g0 = pmg_bg0(b, slot)
          p0 = pmg_bp0(b, slot)
          if (b < nb) then
             g1 = pmg_bg0(b+1, slot)
             p1 = pmg_bp0(b+1, slot)
          else
             g1 = pmg_ng(slot)
             p1 = pmg_np(slot)
          end if
          ng_b = g1 - g0
          np_b = p1 - p0
          do j = 1, ng_b
             r_grid(j) = pmg_rgrid(g0+j, slot)
             r_cell(j) = father(r_grid(j))
             do i = 1, ndim
                r_x0(j,i) = xg(r_grid(j),i) - 3.0D0*dx
             end do
          end do
          do j = 1, np_b
             r_part(j) = pmg_idx(p0+j, slot)
             r_gpart(j) = pmg_pg(p0+j, slot) - g0 + 1
          end do
          ! cic_amr repeats the multipole sums internally; subtract the
          ! pack-time accumulation for this batch to avoid double count
          if (ilevel == levelmin) then
             mm_p=0d0; mx_p=0d0; my_p=0d0; mz_p=0d0
             do j = 1, np_b
                p = r_part(j)
                mm_p = mm_p + mp(p)
                mx_p = mx_p + mp(p)*xp(p,1)
                my_p = my_p + mp(p)*xp(p,2)
                mz_p = mz_p + mp(p)*xp(p,3)
             end do
!$omp atomic update
             multipole(1) = multipole(1) - mm_p
!$omp atomic update
             multipole(2) = multipole(2) - mx_p
!$omp atomic update
             multipole(3) = multipole(3) - my_p
!$omp atomic update
             multipole(4) = multipole(4) - mz_p
          end if
!$omp critical (cic_deposit)
          call cic_amr(r_cell, r_part, r_gpart, r_x0, ng_b, np_b, ilevel)
!$omp end critical (cic_deposit)
       end do
    end if

    call pm_gpu_reset(slot)

  end subroutine pm_gpu_flush_rho

  !------------------------------------------------------------------
  ! Add the GPU-accumulated deposits into the host rho / phi work
  ! arrays (called once after the dispatch loop, outside OMP).
  !------------------------------------------------------------------
  subroutine pm_rho_merge()
    use amr_commons
    use pm_gpu_commons
    use poisson_commons, only: rho, phi
    use particle_cuda_interface
    use iso_c_binding
    implicit none
    integer(c_long_long) :: ncell_c
    integer(i8b) :: i, ntot

    ntot = int(ncoarse,i8b) + int(twotondim,i8b)*int(ngridmax,i8b)
    if (.not.allocated(pmg_rho_add)) then
       allocate(pmg_rho_add(1:ntot))
       allocate(pmg_phiw_add(1:ntot))
    else if (size(pmg_rho_add,kind=i8b) < ntot) then
       deallocate(pmg_rho_add, pmg_phiw_add)
       allocate(pmg_rho_add(1:ntot))
       allocate(pmg_phiw_add(1:ntot))
    end if

    ncell_c = int(ntot, c_long_long)
    call cuda_pm_rho_end_c(pmg_rho_add, pmg_phiw_add, ncell_c)

!$omp parallel do private(i) schedule(static)
    do i = 1, ntot
       rho(i) = rho(i) + pmg_rho_add(i)
       phi(i) = phi(i) + pmg_phiw_add(i)
    end do
!$omp end parallel do

  end subroutine pm_rho_merge

end module pm_gpu_dispatch

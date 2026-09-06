!###########################################################################
! restore_hdf5.f90 — HDF5 parallel input for cuRAMSES
!
! Reads from output_NNNNN/data_NNNNN.h5 and restores:
!   - AMR tree (grid positions, son_flag, cpu_map)
!   - Hydro (uold conservative variables)
!   - Gravity (phi, f)
!   - Particles (xp, vp, mp, idp, levelp, tp, zp)
!   - Sinks
!
! Supports distributed variable-ncpu/cross-ordering restart to Hilbert:
!   - An xg-only pre-pass builds grid-load-balanced Hilbert bounds
!   - Recomputes cpu_map via generic cmp_cpumap
!   - Ranks read disjoint row slices and redistribute locally owned grids
!   - First coarse step forces load_balance for optimal rebalancing
! HDF5 variable-ncpu/cross-ordering restart to ksection is rejected explicitly;
! the existing same-ncpu ksection distributed path remains available.
!###########################################################################
#ifdef HDF5

!###########################################################################
! AMR restore from HDF5
!###########################################################################
subroutine restore_amr_hdf5()
  use amr_commons
  use amr_parameters, only: restart_output_index
  use hydro_commons
  use ksection
  use morton_keys
  use morton_hash
  use ramses_hdf5_io
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, iskip, idim, info
  integer :: hdf5_attr_status, hdf5_attr_status_all
  integer :: nlevelmax_file, nlevelmax_header
  integer :: ksec_kmax_file, ksec_nbinodes_file
  integer(i8b) :: ngrid_total
  integer, allocatable :: ngrid_per_cpu(:)
  integer, allocatable :: son_flag_buf(:), cpu_map_buf(:)
  real(dp), allocatable :: xg_all(:,:)
  integer(HID_T) :: grp_id, lvl_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str
  character(len=5) :: nchar
  character(len=128) :: ordering_file

  ! Grid creation variables
  integer :: igrid_new, igrid_prev_cpu, igrid_father, ind_cell
  integer(8) :: ix, iy, iz, ix_p, iy_p, iz_p
  integer :: nxny
  integer :: ngrid_all_int, icpu, grid_offset

  ! nbor rebuild variables
  integer :: j_nb, ig_nb, ih_nb, ind_oct_nb, igrid_par_nb, nbor_grid
  integer :: ggg_nb(8,6), hhh_nb(8,6)
  type(mkey_t) :: mkey
  real(dp) :: twotol

  ! Variable-ncpu specific variables
  integer :: father_cell
  ! cmp_cpumap expects x(nvector,ndim) and c(nvector) with explicit nvector
  ! stride. Undersizing these buffers makes the callee read past the caller's
  ! storage (stack OOB) and return garbage cpu ids — silent SIGSEGV downstream.
  real(dp) :: scale, dx, xx_cell(nvector,3), xc_off(3)
  integer :: c_tmp(nvector)
  integer :: nx_loc
  ! Vectorized cpu_map computation buffers (twotondim cells per grid)
  real(dp) :: xx_cells(nvector, 3)
  integer :: c_cells(nvector)
  integer :: ix_off, iy_off, iz_off
  ! Streaming chunk buffers for varcpu Phase 1 (avoid full-N allocation)
  real(dp), allocatable :: xg_chunk(:,:)
  integer, allocatable :: son_flag_chunk(:)
  integer :: chunk_size, this_chunk, j, i_file
  integer(i8b) :: i_global
  ! Tier 2 parallel restore: per-rank slice + ksection redistribution
  integer(i8b) :: my_offset_i8, my_end_i8, my_count_i8
  integer(i8b) :: chunk_off_i8, global_off_i8
  integer :: nchunk_local, nchunk_global, ichunk, nprops_pack, nrecv
  real(dp), allocatable :: sendbuf(:,:), recvbuf(:,:)
  integer,  allocatable :: dest_cpu_arr(:)
  ! Hilbert varcpu xg-only pre-pass and exact main-pass diagnostic
  integer, parameter :: varcpu_min_nbin = 65536
  integer, parameter :: varcpu_bins_per_domain = 256
  integer :: varcpu_nbin, varcpu_ibin, varcpu_ilo, varcpu_ihi, varcpu_imid
  integer :: varcpu_nvec, varcpu_j0
  integer(kind=8), allocatable :: varcpu_hist(:), varcpu_domain_load(:)
  integer(kind=8) :: varcpu_total_grids, varcpu_gmin, varcpu_gmax
  real(dp) :: varcpu_gmean, varcpu_imbalance
  real(qdp) :: varcpu_key_span, varcpu_key(1:nvector)

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore AMR: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  nxny = nx * ny

  !=====================================================
  ! Step 1: Read header (time, cosmology, nstep, etc.)
  !=====================================================
  call hdf5_open_group('/header', grp_id)
  call hdf5_read_attr_int_checked(grp_id, 'ncpu', ncpu_file, hdf5_attr_status)
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0 .or. ncpu_file <= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/ncpu status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  nlevelmax_header = 0
  call hdf5_read_attr_int_checked(grp_id, 'nlevelmax', nlevelmax_header, hdf5_attr_status)
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0 .or. nlevelmax_header < 1) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/nlevelmax status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  if(nlevelmax_header>nlevelmax) then
     if(myid==1) write(*,'(A,I4,A,I4)') &
          ' HDF5 restart level guard: checkpoint levelmax=', nlevelmax_header, &
          ' exceeds requested levelmax=', nlevelmax
     call hdf5_close_group(grp_id)
     call hdf5_close_file()
     call clean_stop
  end if
  call hdf5_restore_header_dp_checked(grp_id, 'time', t)
  call hdf5_restore_header_int_checked(grp_id, 'nstep', nstep)
  call hdf5_restore_header_int_checked(grp_id, 'nstep_coarse', nstep_coarse)
  ! Output scheduling (noutput/tout/aout) is set by the namelist.
  ! Only restore the output counters (iout/ifout) from the file.  The stored
  ! tout/aout arrays are deliberately not read: the namelist is authoritative
  ! and a corrupt noutput attribute must not control an allocation or HDF5
  ! attribute read length for values that are immediately discarded.
  block
    call hdf5_restore_header_int_checked(grp_id, 'iout', iout)
    call hdf5_restore_header_int_checked(grp_id, 'ifout', ifout)
  end block
  block
     real(dp), allocatable :: dtold_file(:), dtnew_file(:)
     allocate(dtold_file(nlevelmax_header), dtnew_file(nlevelmax_header))
     call hdf5_restore_header_dp_array_checked(grp_id, 'dtold', dtold_file, &
          nlevelmax_header)
     call hdf5_restore_header_dp_array_checked(grp_id, 'dtnew', dtnew_file, &
          nlevelmax_header)
     dtold = 0.0d0
     dtnew = 0.0d0
     dtold(1:nlevelmax_header) = dtold_file
     dtnew(1:nlevelmax_header) = dtnew_file
     deallocate(dtold_file, dtnew_file)
  end block
  call hdf5_restore_header_dp_checked(grp_id, 'const', const)
  call hdf5_restore_header_dp_checked(grp_id, 'mass_tot_0', mass_tot_0)
  call hdf5_restore_header_dp_checked(grp_id, 'rho_tot', rho_tot)
  call hdf5_restore_header_dp_checked(grp_id, 'omega_m', omega_m)
  call hdf5_restore_header_dp_checked(grp_id, 'omega_l', omega_l)
  call hdf5_restore_header_dp_checked(grp_id, 'omega_k', omega_k)
  call hdf5_restore_header_dp_checked(grp_id, 'omega_b', omega_b)
  call hdf5_restore_header_dp_checked(grp_id, 'h0', h0)
  call hdf5_restore_header_dp_checked(grp_id, 'aexp_ini', aexp_ini)
  call hdf5_restore_header_dp_checked(grp_id, 'boxlen_ini', boxlen_ini)
  call hdf5_restore_header_dp_checked(grp_id, 'aexp', aexp)
  call hdf5_restore_header_dp_checked(grp_id, 'hexp', hexp)
  call hdf5_restore_header_dp_checked(grp_id, 'aexp_old', aexp_old)
  call hdf5_restore_header_dp_checked(grp_id, 'epot_tot_int', epot_tot_int)
  call hdf5_restore_header_dp_checked(grp_id, 'epot_tot_old', epot_tot_old)
  call hdf5_restore_header_dp_checked(grp_id, 'mass_sph', mass_sph)
  call hdf5_restore_header_string_checked(grp_id, 'ordering', ordering_file)
  call hdf5_close_group(grp_id)

  iout=restart_output_index(cosmo,noutput,aout,tout,aexp,t)
  nstep_coarse_old = nstep_coarse

  if(myid==1) write(*,*) 'Restarting at t=', t, ' nstep_coarse=', nstep_coarse

  !=====================================================
  ! Step 2: Check ncpu_file vs ncpu
  !=====================================================
  if(ncpu_file /= ncpu .or. trim(ordering_file) /= trim(ordering)) then
     ! Variable-ncpu, ksection, or cross-ordering: use distributed varcpu path
     ! (the same-ncpu path replicates ALL grids on every rank, requiring
     !  ngridmax >= total_grids, which is prohibitively expensive)
     varcpu_restart = .true.
     varcpu_restart_done = .true.   ! flag in amr_commons: force load_balance on first step
     if(myid==1) then
        write(*,*) '============================================================'
        write(*,*) '============================================================'
        write(*,*) '=====  VARIABLE NCPU / CROSS-ORDERING RESTART          ====='
        write(*,*) '=====  File ncpu = ', ncpu_file, '  ordering = ', trim(ordering_file)
        write(*,*) '=====  New  ncpu = ', ncpu,      '  ordering = ', trim(ordering)
        write(*,*) '=====  Rebuilding domain decomposition                 ====='
        write(*,*) '=====  Load balance will be forced on first coarse step ====='
        write(*,*) '============================================================'
        write(*,*) '============================================================'
     end if
     ! Allocate mapping arrays
     allocate(varcpu_grid_file_idx(1:ngridmax))
     varcpu_grid_file_idx = 0
     allocate(varcpu_ngrid_file(1:nlevelmax))
     varcpu_ngrid_file = 0
  else
     varcpu_restart = .false.
  end if

  !=====================================================
  ! Step 2b: Read coarse grid (all ranks read full coarse)
  !=====================================================
  call hdf5_open_group('/coarse', grp_id)
  call hdf5_read_dataset_all_int_checked(grp_id, 'son', son(1:ncoarse), ncoarse, &
       hdf5_attr_status)
  if(hdf5_attr_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for coarse/son status=', &
          hdf5_attr_status
     call hdf5_restart_abort
  end if
  call hdf5_read_dataset_all_int_checked(grp_id, 'cpu_map', cpu_map(1:ncoarse), ncoarse, &
       hdf5_attr_status)
  if(hdf5_attr_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for coarse/cpu_map status=', &
          hdf5_attr_status
     call hdf5_restart_abort
  end if
  ! Try to read flag1 (new format); fall back to computing from son
  block
     integer :: h5err
     logical :: dset_exists
     call h5lexists_f(grp_id, 'flag1', dset_exists, h5err)
     if(dset_exists .and. h5err == 0) then
        call hdf5_read_dataset_all_int_checked(grp_id, 'flag1', flag1(1:ncoarse), ncoarse, &
             hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for coarse/flag1 status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if
     else
        do i = 1, ncoarse
           if(son(i) > 0) then
              flag1(i) = 1
           else
              flag1(i) = 0
           end if
        end do
     end if
  end block
  call hdf5_close_group(grp_id)

  ! Reset coarse son — will be set when level 1 grids are created
  son(1:ncoarse) = 0

  !=====================================================
  ! Step 3: Domain decomposition
  !=====================================================
  if(varcpu_restart) then
     !--- Variable ncpu: rebuild domain decomposition for new ncpu ---
     nx_loc = icoarse_max - icoarse_min + 1
     scale = boxlen / dble(nx_loc)

     if(ordering == 'ksection') then
        ! A variable-ncpu or cross-ordering restart cannot safely reuse a
        ! ksection tree: its dimensions and ownership are tied to the file
        ! decomposition.  The same-ncpu/same-ordering case is handled by the
        ! checked checkpoint-tree branch below.
        if(myid==1) then
           write(*,*) 'ERROR: HDF5 variable-ncpu/cross-ordering restart to ksection is unsupported.'
           write(*,*) '       Restart with Hilbert ordering or the checkpoint ncpu/ordering.'
        end if
        call hdf5_restart_abort
     else
        ! Hilbert ordering: compute the global key extent now.  The actual
        ! grid-balanced bounds are placed by the xg-only pre-pass below.
        block
           real(qdp), dimension(1:1) :: order_min_tmp, order_max_tmp
           real(dp) :: x_tmp(1:1,1:3), dx_loc
           integer(8) :: iix, iiy, iiz

           order_all_min =  huge(0.0_qdp)
           order_all_max = -huge(0.0_qdp)
           dx_loc = 0.5d0**1
           do iiz = kcoarse_min, kcoarse_max
           do iiy = jcoarse_min, jcoarse_max
           do iix = icoarse_min, icoarse_max
              x_tmp(1,1) = (dble(iix) + 0.5d0 - dble(icoarse_min)) * scale
              x_tmp(1,2) = (dble(iiy) + 0.5d0 - dble(jcoarse_min)) * scale
              x_tmp(1,3) = (dble(iiz) + 0.5d0 - dble(kcoarse_min)) * scale
              call cmp_minmaxorder(x_tmp, order_min_tmp, order_max_tmp, dx_loc, 1)
              order_all_min = min(order_all_min, order_min_tmp(1))
              order_all_max = max(order_all_max, order_max_tmp(1))
           end do
           end do
           end do
           if(.not. allocated(bound_key)) allocate(bound_key(0:ndomain))
           if(.not. allocated(bound_key2)) allocate(bound_key2(0:ndomain))
        end block
     end if
  else
     !--- Same ncpu: read ksection tree from file ---
     call hdf5_open_group('/domain', grp_id)
     if(ordering=='ksection') then
        call hdf5_read_attr_int_checked(grp_id, 'ksec_kmax', ksec_kmax_file, hdf5_attr_status)
        call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
             MPI_MAX, MPI_COMM_WORLD, info)
        if(hdf5_attr_status_all /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for domain/ksec_kmax status=', &
                hdf5_attr_status_all
           call hdf5_restart_abort
        end if
        call hdf5_read_attr_int_checked(grp_id, 'ksec_nbinodes', ksec_nbinodes_file, hdf5_attr_status)
        call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
             MPI_MAX, MPI_COMM_WORLD, info)
        if(hdf5_attr_status_all /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for domain/ksec_nbinodes status=', &
                hdf5_attr_status_all
           call hdf5_restart_abort
        end if
        if(ksec_kmax_file < 2 .or. ksec_kmax_file > ncpu .or. &
             ksec_nbinodes_file < 1 .or. ksec_nbinodes_file > 2*ncpu-1) then
           if(myid==1) write(*,*) 'ERROR: invalid HDF5 ksection dimensions kmax=', &
                ksec_kmax_file, ' nbinodes=', ksec_nbinodes_file, ' ncpu=', ncpu
           call hdf5_restart_abort
        end if
        if(.not. allocated(ksec_wall) .or. .not. allocated(ksec_next) .or. &
             .not. allocated(ksec_indx)) then
           if(myid==1) write(*,*) 'ERROR: ksection arrays are not allocated before HDF5 restore'
           call hdf5_restart_abort
        end if
        if(ksec_kmax_file /= ksec_kmax .or. ksec_nbinodes_file /= ksec_nbinodes .or. &
             size(ksec_wall,1) /= ksec_nbinodes_file .or. &
             size(ksec_next,1) /= ksec_nbinodes_file .or. &
             size(ksec_indx) /= ksec_nbinodes_file .or. &
             size(ksec_wall,2) /= max(ksec_kmax_file-1,1) .or. &
             size(ksec_next,2) /= ksec_kmax_file) then
           if(myid==1) write(*,*) 'ERROR: HDF5 ksection dimensions do not match runtime allocation', &
                ' file=', ksec_kmax_file, ksec_nbinodes_file, &
                ' runtime=', ksec_kmax, ksec_nbinodes
           call hdf5_restart_abort
        end if
        call hdf5_read_dataset_all_int_checked(grp_id, 'ksec_factor', ksec_factor, &
             nksec_levels, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ksec_factor status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if
        call hdf5_read_dataset_all_int_checked(grp_id, 'ksec_dir', ksec_dir, &
             nksec_levels, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ksec_dir status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if
        call hdf5_read_dataset_all_int_checked(grp_id, 'ksec_indx', ksec_indx, &
             ksec_nbinodes, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ksec_indx status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if
        ! Read flattened 2D arrays
        block
           integer :: nw_cols, nn_cols, nw_total, nn_total
           real(dp), allocatable :: wall_flat(:)
           integer, allocatable :: next_flat(:)
           nw_cols = max(ksec_kmax - 1, 1)
           nn_cols = ksec_kmax
           nw_total = ksec_nbinodes * nw_cols
           nn_total = ksec_nbinodes * nn_cols
           allocate(wall_flat(nw_total))
           allocate(next_flat(nn_total))
           call hdf5_read_dataset_all_dp_checked(grp_id, 'ksec_wall', wall_flat, &
                nw_total, hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ksec_wall status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
           call hdf5_read_dataset_all_int_checked(grp_id, 'ksec_next', next_flat, &
                nn_total, hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ksec_next status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
           ksec_wall(1:ksec_nbinodes, 1:nw_cols) = &
                reshape(wall_flat, (/ksec_nbinodes, nw_cols/))
           ksec_next(1:ksec_nbinodes, 1:nn_cols) = &
                reshape(next_flat, (/ksec_nbinodes, nn_cols/))
           deallocate(wall_flat, next_flat)
        end block
        ! Read CPU bounding boxes
        block
           real(dp), allocatable :: cpubox_flat(:)
           integer :: ncb
           ncb = ncpu * ndim
           allocate(cpubox_flat(ncb))
           call hdf5_read_dataset_all_dp_checked(grp_id, 'bisec_cpubox_min', cpubox_flat, ncb, &
                hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for bisec_cpubox_min status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
           bisec_cpubox_min(1:ncpu, 1:ndim) = reshape(cpubox_flat, (/ncpu, ndim/))
           call hdf5_read_dataset_all_dp_checked(grp_id, 'bisec_cpubox_max', cpubox_flat, ncb, &
                hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for bisec_cpubox_max status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
           bisec_cpubox_max(1:ncpu, 1:ndim) = reshape(cpubox_flat, (/ncpu, ndim/))
           deallocate(cpubox_flat)
        end block
        bisec_cpubox_min2 = bisec_cpubox_min
        bisec_cpubox_max2 = bisec_cpubox_max
        ! Rebuild tree navigation arrays
        call rebuild_ksec_cpuranges()
        call compute_ksec_cpu_path()
        if(myid==1) write(*,*) 'HDF5: ksection tree restored'
     end if
     call hdf5_close_group(grp_id)
  end if

  ! Set cpu_map2 = cpu_map for coarse
  cpu_map2(1:ncoarse) = cpu_map(1:ncoarse)

  !=====================================================
  ! Step 4: Read nlevelmax_file from /amr group
  !=====================================================
  call hdf5_open_group('/amr', grp_id)
  nlevelmax_file = 0
  call hdf5_read_attr_int_checked(grp_id, 'nlevelmax_file', nlevelmax_file, hdf5_attr_status)
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  call hdf5_close_group(grp_id)
  if(hdf5_attr_status_all /= 0 .or. nlevelmax_file < 1 .or. &
       nlevelmax_file > nlevelmax) then
     if(myid==1) write(*,*) 'ERROR: invalid HDF5 nlevelmax_file=', nlevelmax_file, &
          ' requested=', nlevelmax, ' status=', hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  if(myid==1) write(*,*) 'HDF5: nlevelmax_file =', nlevelmax_file

  !=====================================================
  ! Step 5: Initialize Morton hash tables
  !=====================================================
  if(.not. allocated(mort_table)) allocate(mort_table(1:nlevelmax))
  if(.not. allocated(grid_level)) then
     allocate(grid_level(1:ngridmax))
     grid_level = 0
  end if
  ! Initialize empty hash tables for all levels
  do ilevel = 1, nlevelmax
     call morton_hash_init(mort_table(ilevel), 16)
  end do

  !=====================================================
  ! Step 6: Per-level grid creation
  !=====================================================
  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()

  ! Compute scale for cell coordinate computation
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)

  if(varcpu_restart) then
     !===================================================
     ! Hilbert xg-only pre-pass.  The level discovery, global row split, and
     ! collective chunk schedule deliberately mirror Stage A below.  No AMR
     ! payload other than xg_1/xg_2/xg_3 is read here.
     !===================================================
     if(ordering /= 'ksection') then
        varcpu_nbin = max(varcpu_min_nbin, ndomain * varcpu_bins_per_domain)
        allocate(varcpu_hist(1:varcpu_nbin))
        varcpu_hist = 0_8
        varcpu_key_span = order_all_max - order_all_min

        do ilevel = 1, nlevelmax_file
           write(lvl_str, '(I0)') ilevel
           grp_name = '/amr/level_'//trim(lvl_str)

           block
              integer(HID_T) :: test_grp_id
              integer :: h5err
              logical :: grp_exists
              call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
              if(.not. grp_exists .or. h5err /= 0) cycle
              call h5gopen_f(hdf5_file_id, trim(grp_name), test_grp_id, h5err)
              if(h5err /= 0) cycle
              lvl_grp_id = test_grp_id
           end block

           allocate(ngrid_per_cpu(ncpu_file))
           call hdf5_read_dataset_all_int_checked(lvl_grp_id, 'ngrid_per_cpu', &
                ngrid_per_cpu, ncpu_file, hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ngrid_per_cpu status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
           ngrid_total = 0_i8b
           do i = 1, ncpu_file
              if(ngrid_per_cpu(i) < 0) then
                 if(myid==1) write(*,*) 'ERROR: invalid ngrid_per_cpu: negative entry at rank=', i
                 call hdf5_restart_abort
              end if
              ngrid_total = ngrid_total + int(ngrid_per_cpu(i), i8b)
           end do
           deallocate(ngrid_per_cpu)
           ngrid_all_int = int(ngrid_total)
           varcpu_ngrid_file(ilevel) = ngrid_all_int

           if(ngrid_all_int == 0) then
              call hdf5_close_group(lvl_grp_id)
              cycle
           end if

           chunk_size = min(1048576, ngrid_all_int)
           my_offset_i8 = int(myid-1, i8b) * int(ngrid_all_int, i8b) / int(ncpu, i8b)
           my_end_i8    = int(myid,   i8b) * int(ngrid_all_int, i8b) / int(ncpu, i8b)
           my_count_i8  = my_end_i8 - my_offset_i8
           nchunk_local = int((my_count_i8 + int(chunk_size,i8b) - 1_i8b) / &
                int(chunk_size,i8b))
#ifndef WITHOUTMPI
           call MPI_ALLREDUCE(nchunk_local, nchunk_global, 1, MPI_INTEGER, &
                MPI_MAX, MPI_COMM_WORLD, info)
#else
           nchunk_global = nchunk_local
#endif
           if(nchunk_global < 1) nchunk_global = 1

           allocate(xg_chunk(max(chunk_size,1), ndim))
           do ichunk = 0, nchunk_global - 1
              chunk_off_i8 = int(ichunk, i8b) * int(chunk_size, i8b)
              if(chunk_off_i8 >= my_count_i8) then
                 this_chunk = 0
                 global_off_i8 = my_offset_i8
              else
                 this_chunk = int(min(int(chunk_size,i8b), &
                      my_count_i8 - chunk_off_i8))
                 global_off_i8 = my_offset_i8 + chunk_off_i8
              end if

              do idim = 1, ndim
                 write(lvl_str, '(I0)') idim
                 call hdf5_read_dataset_collective_dp_checked(lvl_grp_id, &
                      'xg_'//trim(lvl_str), xg_chunk(:,idim), this_chunk, &
                      global_off_i8, int(ngrid_all_int, i8b), hdf5_attr_status)
                 if(hdf5_attr_status /= 0) then
                    if(myid==1) write(*,*) 'ERROR: checked AMR xg read failed status=', &
                         hdf5_attr_status
                    call hdf5_restart_abort
                 end if
              end do

              do varcpu_j0 = 1, this_chunk, nvector
                 varcpu_nvec = min(nvector, this_chunk-varcpu_j0+1)
                 do j = 1, varcpu_nvec
                    xx_cell(j,1) = (xg_chunk(varcpu_j0+j-1,1) - dble(icoarse_min)) * scale
                    xx_cell(j,2) = (xg_chunk(varcpu_j0+j-1,2) - dble(jcoarse_min)) * scale
                    xx_cell(j,3) = (xg_chunk(varcpu_j0+j-1,3) - dble(kcoarse_min)) * scale
                 end do
                 call cmp_ordering(xx_cell, varcpu_key, varcpu_nvec)
                 do j = 1, varcpu_nvec
                    if(varcpu_key_span > 0.0_qdp) then
                       varcpu_ibin = int(real(varcpu_nbin,qdp) * &
                            (varcpu_key(j)-order_all_min) / varcpu_key_span) + 1
                    else
                       varcpu_ibin = 1
                    end if
                    varcpu_ibin = min(varcpu_nbin, max(1,varcpu_ibin))
                    varcpu_hist(varcpu_ibin) = varcpu_hist(varcpu_ibin) + 1_8
                 end do
              end do
           end do
           deallocate(xg_chunk)
           call hdf5_close_group(lvl_grp_id)
        end do

#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, varcpu_hist, varcpu_nbin, MPI_INTEGER8, &
             MPI_SUM, MPI_COMM_WORLD, info)
#endif
        varcpu_total_grids = sum(varcpu_hist)
        call hilbert_bounds_from_hist(varcpu_hist, varcpu_nbin, varcpu_total_grids)
        deallocate(varcpu_hist)
        allocate(varcpu_domain_load(1:ndomain))
        varcpu_domain_load = 0_8
        if(myid==1) write(*,*) 'HDF5: grid-balanced Hilbert bound_key rebuilt for ncpu=', ncpu
     end if

     ! The decomposition must exist before either coarse ownership or streamed
     ! grid redistribution is computed.
     do iz = kcoarse_min, kcoarse_max
     do iy = jcoarse_min, jcoarse_max
     do ix = icoarse_min, icoarse_max
        ind = 1 + int(ix) + int(iy) * nx + int(iz) * nxny
        xx_cell(1,1) = (dble(ix) + 0.5d0 - dble(icoarse_min)) * scale
        xx_cell(1,2) = (dble(iy) + 0.5d0 - dble(jcoarse_min)) * scale
        xx_cell(1,3) = (dble(iz) + 0.5d0 - dble(kcoarse_min)) * scale
        call cmp_cpumap(xx_cell, c_tmp, 1)
        cpu_map(ind) = c_tmp(1)
     end do
     end do
     end do
     cpu_map2(1:ncoarse) = cpu_map(1:ncoarse)
     if(myid==1) write(*,*) 'HDF5: coarse cpu_map recomputed'
  end if

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: distributed active-only grid creation
     ! Each rank creates only its active grids (cpu_map(father)==myid)
     ! Virtual grids are created via refine_coarse/refine_fine
     ! This reduces per-rank memory by ~ncpu× vs reading ALL grids
     !===================================================
     balance = .true.
     shrink = .false.

     do ilevel = 1, nlevelmax_file
        ! Check if level group exists
        write(lvl_str, '(I0)') ilevel
        grp_name = '/amr/level_'//trim(lvl_str)

        block
           integer(HID_T) :: test_grp_id
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), test_grp_id, h5err)
           if(h5err /= 0) cycle
           lvl_grp_id = test_grp_id
        end block

        ! Read ngrid_per_cpu from file (ncpu_file entries)
        allocate(ngrid_per_cpu(ncpu_file))
        call hdf5_read_dataset_all_int_checked(lvl_grp_id, 'ngrid_per_cpu', &
             ngrid_per_cpu, ncpu_file, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ngrid_per_cpu status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if

        ! Compute total grids at this level
        ngrid_total = 0
        do i = 1, ncpu_file
           if(ngrid_per_cpu(i) < 0) then
              if(myid==1) write(*,*) 'ERROR: invalid ngrid_per_cpu: negative entry at rank=', i
              call hdf5_restart_abort
           end if
           ngrid_total = ngrid_total + ngrid_per_cpu(i)
        end do
        ngrid_all_int = int(ngrid_total)
        varcpu_ngrid_file(ilevel) = ngrid_all_int

        if(myid==1) write(*,'(A,I3,A,I10)') &
             ' HDF5 level ', ilevel, ' total grids: ', ngrid_total

        deallocate(ngrid_per_cpu)

        if(ngrid_all_int == 0) then
           call hdf5_close_group(lvl_grp_id)
           cycle
        end if

        ! Initialize hash table for this level.
        ! Cross-ordering (ksection→hilbert etc.) can deliver up to several×
        ! the average grid count to a single rank; oversize by 16× the
        ! average load to avoid the O(N) rehash cascade when imbalance is large.
        call morton_hash_init(mort_table(ilevel), &
             max(16 * (ngrid_all_int / max(ncpu, 1) + 1), 16))

        ! Subcell offsets at this level
        dx = 0.5d0**ilevel

        !========================================================
        ! Tier 2 parallel restore (exascale-ready):
        !   Stage A: each rank reads a disjoint global-index slice
        !            in parallel via collective hyperslab.
        !   Stage B: pack (xg, son_flag, i_file) + dest_cpu and
        !            redistribute via recursive k-section tree.
        !   Stage C: receiver creates grids exactly like the old
        !            single-rank chunked path.
        ! Sub-chunked to cap transient memory at chunk_size×96 B
        ! (~96 MB at 1M). All ranks participate in nchunk_global
        ! rounds (MAX across ranks) so collective ops stay aligned.
        !========================================================
        chunk_size = min(1048576, ngrid_all_int)

        ! Uniform global-index split: rank r gets [r*N/P, (r+1)*N/P)
        my_offset_i8 = int(myid-1, i8b) * int(ngrid_all_int, i8b) / int(ncpu, i8b)
        my_end_i8    = int(myid,   i8b) * int(ngrid_all_int, i8b) / int(ncpu, i8b)
        my_count_i8  = my_end_i8 - my_offset_i8

        nchunk_local = int( (my_count_i8 + int(chunk_size,i8b) - 1_i8b) &
                            / int(chunk_size,i8b) )
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nchunk_local, nchunk_global, 1, MPI_INTEGER, &
             MPI_MAX, MPI_COMM_WORLD, info)
#else
        nchunk_global = nchunk_local
#endif
        if(nchunk_global < 1) nchunk_global = 1   ! safety: at least one collective round

        ! Pack layout: xg(ndim) | son_flag(twotondim) | i_file → all as dp
        nprops_pack = ndim + twotondim + 1
        allocate(xg_chunk(max(chunk_size,1), ndim))
        allocate(son_flag_chunk(max(chunk_size,1) * twotondim))
        allocate(sendbuf(nprops_pack, max(chunk_size,1)))
        allocate(dest_cpu_arr(max(chunk_size,1)))

        igrid_prev_cpu = 0
        do ichunk = 0, nchunk_global - 1
           ! My contribution to this collective round (may be 0)
           chunk_off_i8 = int(ichunk, i8b) * int(chunk_size, i8b)
           if(chunk_off_i8 >= my_count_i8) then
              this_chunk = 0
              global_off_i8 = my_offset_i8   ! n=0 path ignores offset
           else
              this_chunk = int(min(int(chunk_size,i8b), &
                                   my_count_i8 - chunk_off_i8))
              global_off_i8 = my_offset_i8 + chunk_off_i8
           end if

           !--- Stage A: parallel hyperslab read (collective) ---
           do idim = 1, ndim
              write(lvl_str, '(I0)') idim
              call hdf5_read_dataset_collective_dp_checked(lvl_grp_id, &
                   'xg_'//trim(lvl_str), xg_chunk(:, idim), this_chunk, &
                   global_off_i8, int(ngrid_all_int, i8b), hdf5_attr_status)
              if(hdf5_attr_status /= 0) then
                 if(myid==1) write(*,*) 'ERROR: checked AMR xg read failed status=', &
                      hdf5_attr_status
                 call hdf5_restart_abort
              end if
           end do
           call hdf5_read_dataset_collective_int_checked(lvl_grp_id, 'son_flag', &
                son_flag_chunk, this_chunk * twotondim, &
                global_off_i8 * int(twotondim, i8b), &
                int(ngrid_all_int, i8b) * int(twotondim, i8b), hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: checked AMR son_flag read failed status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if

           !--- Stage B prep: dest_cpu via cmp_cpumap + pack ---
           do j = 1, this_chunk
              xx_cell(1,1) = (xg_chunk(j,1) - dble(icoarse_min)) * scale
              xx_cell(1,2) = (xg_chunk(j,2) - dble(jcoarse_min)) * scale
              xx_cell(1,3) = (xg_chunk(j,3) - dble(kcoarse_min)) * scale
              call cmp_cpumap(xx_cell, c_tmp, 1)
              dest_cpu_arr(j) = c_tmp(1)

              sendbuf(1:ndim, j) = xg_chunk(j, 1:ndim)
              do iskip = 1, twotondim
                 sendbuf(ndim + iskip, j) = &
                      dble(son_flag_chunk((j-1)*twotondim + iskip))
              end do
              ! 1-based global file row index (≤ 2^31-1 in current scale,
              ! exactly representable as dp up to 2^53)
              sendbuf(nprops_pack, j) = dble(global_off_i8 + int(j, i8b))
           end do

           ! Exact Hilbert domain loads for the binary-compatible diagnostic.
           ! This reuses xg already resident in the main pass and performs no
           ! additional HDF5 reads.
           if(ordering /= 'ksection') then
              do varcpu_j0 = 1, this_chunk, nvector
                 varcpu_nvec = min(nvector, this_chunk-varcpu_j0+1)
                 do j = 1, varcpu_nvec
                    xx_cell(j,1) = (xg_chunk(varcpu_j0+j-1,1) - dble(icoarse_min)) * scale
                    xx_cell(j,2) = (xg_chunk(varcpu_j0+j-1,2) - dble(jcoarse_min)) * scale
                    xx_cell(j,3) = (xg_chunk(varcpu_j0+j-1,3) - dble(kcoarse_min)) * scale
                 end do
                 call cmp_ordering(xx_cell, varcpu_key, varcpu_nvec)
                 do j = 1, varcpu_nvec
                    varcpu_ilo = 1
                    varcpu_ihi = ndomain
                    do while(varcpu_ilo < varcpu_ihi)
                       varcpu_imid = (varcpu_ilo + varcpu_ihi) / 2
                       if(varcpu_key(j) < bound_key(varcpu_imid)) then
                          varcpu_ihi = varcpu_imid
                       else
                          varcpu_ilo = varcpu_imid + 1
                       end if
                    end do
                    varcpu_domain_load(varcpu_ilo) = varcpu_domain_load(varcpu_ilo) + 1_8
                 end do
              end do
           end if

           !--- Stage B: recursive k-section redistribution ---
           call ksection_exchange_dp(sendbuf, this_chunk, dest_cpu_arr, &
                nprops_pack, recvbuf, nrecv)

           !--- Stage C: create received grids (already this rank's) ---
           do j = 1, nrecv
              i_file = nint(recvbuf(nprops_pack, j))

              if(ilevel == 1) then
                 twotol = 1.0d0
                 ix = int(recvbuf(1,j) * twotol, 8)
                 iy = int(recvbuf(2,j) * twotol, 8)
                 iz = int(recvbuf(3,j) * twotol, 8)
                 father_cell = 1 + int(ix) + int(iy) * nx + int(iz) * nxny
              else
                 twotol = 2.0d0**(ilevel-1)
                 ix = int(recvbuf(1,j) * twotol, 8)
                 iy = int(recvbuf(2,j) * twotol, 8)
                 iz = int(recvbuf(3,j) * twotol, 8)
                 ix_p = ix / 2_8
                 iy_p = iy / 2_8
                 iz_p = iz / 2_8
                 mkey = morton_encode(ix_p, iy_p, iz_p)
                 igrid_father = morton_hash_lookup(mort_table(ilevel-1), mkey)
                 if(igrid_father == 0) cycle  ! parent not on this rank
                 ind_cell = 1 + int(mod(ix,2_8)) + 2*int(mod(iy,2_8)) + 4*int(mod(iz,2_8))
                 father_cell = ICELL_OF(igrid_father,ind_cell)
              end if

              ! Defensive: ksection routing matches cpu_map decomposition,
              ! but skip on disagreement (e.g., float boundary cases)
              if(cpu_map(father_cell) /= myid) cycle

              igrid_new = headf
              if(igrid_new == 0) then
                 write(*,*) 'ERROR: out of free grids, myid=', myid, &
                      ' level=', ilevel, ' file_row=', i_file
                 call clean_stop
              end if
              headf = next(igrid_new)
              if(headf > 0) prev(headf) = 0
              numbf = numbf - 1

              xg(igrid_new, 1:ndim) = recvbuf(1:ndim, j)

              do iskip = 1, twotondim
                 ind_cell = ICELL_OF(igrid_new,iskip)
                 son(ind_cell) = 0
                 flag1(ind_cell) = nint(recvbuf(ndim + iskip, j))
              end do

              ! Vectorized cmp_cpumap for 8 sub-cells
              do iskip = 1, twotondim
                 iz_off = (iskip - 1) / 4
                 iy_off = (iskip - 1 - 4 * iz_off) / 2
                 ix_off = (iskip - 1 - 2 * iy_off - 4 * iz_off)
                 xx_cells(iskip, 1) = (recvbuf(1,j) + (dble(ix_off) - 0.5d0) * dx &
                      - dble(icoarse_min)) * scale
                 xx_cells(iskip, 2) = (recvbuf(2,j) + (dble(iy_off) - 0.5d0) * dx &
                      - dble(jcoarse_min)) * scale
                 xx_cells(iskip, 3) = (recvbuf(3,j) + (dble(iz_off) - 0.5d0) * dx &
                      - dble(kcoarse_min)) * scale
              end do
              call cmp_cpumap(xx_cells, c_cells, twotondim)
              do iskip = 1, twotondim
                 ind_cell = ICELL_OF(igrid_new,iskip)
                 cpu_map(ind_cell)  = c_cells(iskip)
                 cpu_map2(ind_cell) = c_cells(iskip)
              end do

              father(igrid_new) = father_cell
              son(father_cell)  = igrid_new

              mkey = grid_to_morton(igrid_new, ilevel)
              call morton_hash_insert(mort_table(ilevel), mkey, igrid_new)
              grid_level(igrid_new) = ilevel

              if(igrid_prev_cpu == 0) then
                 headl(myid, ilevel) = igrid_new
              else
                 next(igrid_prev_cpu) = igrid_new
              end if
              prev(igrid_new) = igrid_prev_cpu
              next(igrid_new) = 0
              taill(myid, ilevel) = igrid_new
              numbl(myid, ilevel) = numbl(myid, ilevel) + 1
              igrid_prev_cpu = igrid_new

              varcpu_grid_file_idx(igrid_new) = i_file
           end do

           if(allocated(recvbuf)) deallocate(recvbuf)
        end do

        deallocate(xg_chunk, son_flag_chunk, sendbuf, dest_cpu_arr)
        call hdf5_close_group(lvl_grp_id)

        ! Tier 2: ksection exchange delivers grids in non-monotonic file-idx
        ! order. restore_hydro/poisson walk the active list with an early-exit
        ! "fidx > chunk_end" check that assumes ascending order. Sort the active
        ! linked list by varcpu_grid_file_idx so subsequent restore stages see
        ! a monotonic chain.
        block
           integer :: nactive, k, ig_tmp
           integer, allocatable :: igrid_arr(:), perm(:)
           real(dp), allocatable :: fkey(:)
           nactive = numbl(myid, ilevel)
           if(nactive > 1) then
              allocate(igrid_arr(nactive), fkey(nactive), perm(nactive))
              ig_tmp = headl(myid, ilevel)
              k = 0
              do while(ig_tmp > 0)
                 k = k + 1
                 if(k > nactive) exit  ! defensive: numbl mismatch
                 igrid_arr(k) = ig_tmp
                 fkey(k) = dble(varcpu_grid_file_idx(ig_tmp))
                 ig_tmp = next(ig_tmp)
              end do
              if(k == nactive) then
                 call quick_sort_dp(fkey, perm, nactive)
                 headl(myid, ilevel) = igrid_arr(perm(1))
                 prev(igrid_arr(perm(1))) = 0
                 do k = 1, nactive - 1
                    next(igrid_arr(perm(k)))   = igrid_arr(perm(k+1))
                    prev(igrid_arr(perm(k+1))) = igrid_arr(perm(k))
                 end do
                 next(igrid_arr(perm(nactive))) = 0
                 taill(myid, ilevel) = igrid_arr(perm(nactive))
              end if
              deallocate(igrid_arr, fkey, perm)
           end if
        end block

        !--- Phase 2: Create virtual grids via RAMSES refine mechanism ---
        ! Active grids already set son(father)>0, so refine skips those cells
        ! Only cells with flag1==1 AND son==0 get new (virtual) grids
        if(ilevel == 1) then
           call flag_coarse
           call refine_coarse
        else
           call refine_fine(ilevel - 1)
        end if

        ! Build communicators for this level
        call build_comm(ilevel)

        ! Exchange flag1/cpu_map from active to virtual grids
        ! (needed for refine_fine at next level)
        call make_virtual_fine_int(flag1(1), ilevel)
        call make_virtual_fine_int(cpu_map(1), ilevel)
        ! cpu_map2 == cpu_map during restore: copy locally over virtual grids
        ! instead of a third MPI collective exchange.
        ! Per-icpu virtual grid linked lists occupy disjoint igrid slots, so
        ! the outer loop is OMP-safe (no write collisions across icpus).
!$OMP PARALLEL DO PRIVATE(icpu,igrid,iskip,ind_cell) SCHEDULE(DYNAMIC)
        do icpu = 1, ncpu
           if(icpu == myid) cycle
           igrid = headl(icpu, ilevel)
           do while(igrid > 0)
              do iskip = 1, twotondim
                 ind_cell = ICELL_OF(igrid,iskip)
                 cpu_map2(ind_cell) = cpu_map(ind_cell)
              end do
              igrid = next(igrid)
           end do
        end do
!$OMP END PARALLEL DO

        if(myid==1) write(*,'(A,I3,A,I10,A,I10)') &
             ' HDF5 level ', ilevel, ' active: ', numbl(myid, ilevel), &
             ' total: ', numbtot(1, ilevel)
     end do

     if(ordering /= 'ksection') then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(MPI_IN_PLACE, varcpu_domain_load, ndomain, MPI_INTEGER8, &
             MPI_SUM, MPI_COMM_WORLD, info)
#endif
        varcpu_total_grids = sum(varcpu_domain_load)
        varcpu_gmin = minval(varcpu_domain_load)
        varcpu_gmax = maxval(varcpu_domain_load)
        varcpu_gmean = dble(varcpu_total_grids) / dble(max(ndomain,1))
        varcpu_imbalance = 0.0d0
        if(varcpu_gmean > 0.0d0) &
             varcpu_imbalance = dble(varcpu_gmax) / varcpu_gmean
        if(myid == 1) write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES14.6,A,F10.6)') &
             ' VARCPU balance: ndomain=', ndomain, ' total_grids=', varcpu_total_grids, &
             ' per-domain min/max/mean=', varcpu_gmin, '/', varcpu_gmax, '/', &
             varcpu_gmean, ' imbalance=', varcpu_imbalance
        deallocate(varcpu_domain_load)
     end if

     balance = .false.
     shrink = .false.

  else
     !===================================================
     ! SAME NCPU PATH: original code (grids ordered by CPU)
     !===================================================
     allocate(ngrid_per_cpu(ncpu))

     do ilevel = 1, nlevelmax_file
        ! Check if level group exists
        write(lvl_str, '(I0)') ilevel
        grp_name = '/amr/level_'//trim(lvl_str)

        block
           integer(HID_T) :: test_grp_id
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), test_grp_id, h5err)
           if(h5err /= 0) cycle
           lvl_grp_id = test_grp_id
        end block

        ! Read ngrid_per_cpu (all ranks read full array)
        call hdf5_read_dataset_all_int_checked(lvl_grp_id, 'ngrid_per_cpu', &
             ngrid_per_cpu, ncpu, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ngrid_per_cpu status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if

        ! Compute total grids at this level
        ngrid_total = 0
        do i = 1, ncpu
           if(ngrid_per_cpu(i) < 0) then
              if(myid==1) write(*,*) 'ERROR: invalid ngrid_per_cpu: negative entry at rank=', i
              call hdf5_restart_abort
           end if
           ngrid_total = ngrid_total + ngrid_per_cpu(i)
        end do
        ngrid_all_int = int(ngrid_total)

        if(myid==1) write(*,'(A,I3,A,I10)') &
             ' HDF5 level ', ilevel, ' total grids: ', ngrid_total

        if(ngrid_all_int == 0) then
           call hdf5_close_group(lvl_grp_id)
           cycle
        end if

        ! ALL ranks read ALL grids' data at this level
        allocate(xg_all(ngrid_all_int, ndim))
        allocate(son_flag_buf(ngrid_all_int * twotondim))
        allocate(cpu_map_buf(ngrid_all_int * twotondim))

        ! Read ALL grid positions (every rank reads full datasets)
        do idim = 1, ndim
           write(lvl_str, '(I0)') idim
           call hdf5_read_dataset_all_dp_checked(lvl_grp_id, 'xg_'//trim(lvl_str), &
                xg_all(:, idim), ngrid_all_int, hdf5_attr_status)
           if(hdf5_attr_status /= 0) then
              if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for AMR xg status=', &
                   hdf5_attr_status
              call hdf5_restart_abort
           end if
        end do

        ! Read ALL cell data
        call hdf5_read_dataset_all_int_checked(lvl_grp_id, 'son_flag', &
             son_flag_buf, ngrid_all_int * twotondim, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for AMR son_flag status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if
        call hdf5_read_dataset_all_int_checked(lvl_grp_id, 'cpu_map', &
             cpu_map_buf, ngrid_all_int * twotondim, hdf5_attr_status)
        if(hdf5_attr_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: checked HDF5 read failed for AMR cpu_map status=', &
                hdf5_attr_status
           call hdf5_restart_abort
        end if

        ! Initialize hash table for this level
        call morton_hash_init(mort_table(ilevel), max(2 * ngrid_all_int, 16))

        !--- Create grid slots for ALL CPUs' grids ---
        ! Grids are ordered by CPU: CPU 1 first, then CPU 2, etc.
        grid_offset = 0
        do icpu = 1, ncpu
           igrid_prev_cpu = 0
           do i = 1, ngrid_per_cpu(icpu)
              ! Allocate grid from free list
              igrid_new = headf
              if(igrid_new == 0) then
                 write(*,*) 'ERROR: out of free grids, myid=', myid, &
                      ' level=', ilevel, ' cpu=', icpu
                 call clean_stop
              end if
              headf = next(igrid_new)
              if(headf > 0) prev(headf) = 0
              numbf = numbf - 1

              ! Set linked list for this CPU at this level
              if(igrid_prev_cpu == 0) then
                 headl(icpu, ilevel) = igrid_new
              else
                 next(igrid_prev_cpu) = igrid_new
              end if
              prev(igrid_new) = igrid_prev_cpu
              next(igrid_new) = 0
              taill(icpu, ilevel) = igrid_new
              numbl(icpu, ilevel) = numbl(icpu, ilevel) + 1

              ! Set xg, cpu_map, cpu_map2, son, flag1
              ind = grid_offset + i
              xg(igrid_new, 1:ndim) = xg_all(ind, 1:ndim)
              do iskip = 1, twotondim
                 ind_cell = ICELL_OF(igrid_new,iskip)
                 cpu_map(ind_cell) = cpu_map_buf((ind-1)*twotondim + iskip)
                 cpu_map2(ind_cell) = cpu_map_buf((ind-1)*twotondim + iskip)
                 son(ind_cell) = 0
                 flag1(ind_cell) = son_flag_buf((ind-1)*twotondim + iskip)
              end do

              ! Set father pointer
              if(ilevel == 1) then
                 twotol = 1.0d0
                 ix = int(xg_all(ind, 1) * twotol, 8)
                 iy = int(xg_all(ind, 2) * twotol, 8)
                 iz = int(xg_all(ind, 3) * twotol, 8)
                 father(igrid_new) = 1 + int(ix) + int(iy) * nx + int(iz) * nxny
              else
                 twotol = 2.0d0**(ilevel-1)
                 ix = int(xg_all(ind, 1) * twotol, 8)
                 iy = int(xg_all(ind, 2) * twotol, 8)
                 iz = int(xg_all(ind, 3) * twotol, 8)
                 ix_p = ix / 2_8
                 iy_p = iy / 2_8
                 iz_p = iz / 2_8
                 mkey = morton_encode(ix_p, iy_p, iz_p)
                 igrid_father = morton_hash_lookup(mort_table(ilevel-1), mkey)
                 if(igrid_father == 0) then
                    write(*,*) 'ERROR: father not found! myid=', myid, &
                         ' level=', ilevel, ' cpu=', icpu, ' grid=', i, &
                         ' ix_p=', ix_p, ' iy_p=', iy_p, ' iz_p=', iz_p
                    call clean_stop
                 end if
                 ind_cell = 1 + int(mod(ix, 2_8)) + 2 * int(mod(iy, 2_8)) + 4 * int(mod(iz, 2_8))
                 father(igrid_new) = ICELL_OF(igrid_father,ind_cell)
              end if

              ! Set son in parent cell
              son(father(igrid_new)) = igrid_new

              ! Insert into Morton hash
              mkey = grid_to_morton(igrid_new, ilevel)
              call morton_hash_insert(mort_table(ilevel), mkey, igrid_new)
              grid_level(igrid_new) = ilevel

              igrid_prev_cpu = igrid_new
           end do
           grid_offset = grid_offset + ngrid_per_cpu(icpu)
        end do

        deallocate(xg_all, son_flag_buf, cpu_map_buf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_per_cpu)
  end if

  ! Restore HDF5 error printing
  call hdf5_restore_errors()

  !=====================================================
  ! Step 7: Set ngrid_current, update used_mem
  !=====================================================
  ngrid_current = ngridmax - numbf
  used_mem = ngrid_current
  if(myid==1) write(*,*) 'HDF5: ngrid_current =', ngrid_current

  !=====================================================
  ! Step 7b: Rebuild nbor() array for son(nbor) compatibility
  !=====================================================
  if(size(nbor,1) > 1) then
     if(myid==1) write(*,*) 'HDF5: rebuilding nbor() array...'
     nbor = 0

     ! Level 1: compute nbor from coarse cell coordinates
     do icpu = 1, ncpu
        igrid = headl(icpu, 1)
        do while(igrid > 0)
           ix = int(xg(igrid, 1) * 1.0d0, 8)
           iy = int(xg(igrid, 2) * 1.0d0, 8)
           iz = int(xg(igrid, 3) * 1.0d0, 8)
           ! x-direction
           if(ix > 0) then
              nbor(igrid, 1) = 1 + int(ix-1) + int(iy)*nx + int(iz)*nxny
           else
              nbor(igrid, 1) = 1 + (nx-1) + int(iy)*nx + int(iz)*nxny
           end if
           if(ix < nx-1) then
              nbor(igrid, 2) = 1 + int(ix+1) + int(iy)*nx + int(iz)*nxny
           else
              nbor(igrid, 2) = 1 + 0 + int(iy)*nx + int(iz)*nxny
           end if
           ! y-direction
           if(iy > 0) then
              nbor(igrid, 3) = 1 + int(ix) + int(iy-1)*nx + int(iz)*nxny
           else
              nbor(igrid, 3) = 1 + int(ix) + (ny-1)*nx + int(iz)*nxny
           end if
           if(iy < ny-1) then
              nbor(igrid, 4) = 1 + int(ix) + int(iy+1)*nx + int(iz)*nxny
           else
              nbor(igrid, 4) = 1 + int(ix) + 0 + int(iz)*nxny
           end if
           ! z-direction
           if(iz > 0) then
              nbor(igrid, 5) = 1 + int(ix) + int(iy)*nx + int(iz-1)*nxny
           else
              nbor(igrid, 5) = 1 + int(ix) + int(iy)*nx + (nz-1)*nxny
           end if
           if(iz < nz-1) then
              nbor(igrid, 6) = 1 + int(ix) + int(iy)*nx + int(iz+1)*nxny
           else
              nbor(igrid, 6) = 1 + int(ix) + int(iy)*nx + 0
           end if
           igrid = next(igrid)
        end do
     end do

     ! Levels > 1: inline getnborgrids + getnborcells logic
     ! ggg: which neighbor (0=central, 1-6=direction), hhh: cell position in that grid
     ggg_nb(1:8,1)=(/1,0,1,0,1,0,1,0/); hhh_nb(1:8,1)=(/2,1,4,3,6,5,8,7/)
     ggg_nb(1:8,2)=(/0,2,0,2,0,2,0,2/); hhh_nb(1:8,2)=(/2,1,4,3,6,5,8,7/)
     ggg_nb(1:8,3)=(/3,3,0,0,3,3,0,0/); hhh_nb(1:8,3)=(/3,4,1,2,7,8,5,6/)
     ggg_nb(1:8,4)=(/0,0,4,4,0,0,4,4/); hhh_nb(1:8,4)=(/3,4,1,2,7,8,5,6/)
     ggg_nb(1:8,5)=(/5,5,5,5,0,0,0,0/); hhh_nb(1:8,5)=(/5,6,7,8,1,2,3,4/)
     ggg_nb(1:8,6)=(/0,0,0,0,6,6,6,6/); hhh_nb(1:8,6)=(/5,6,7,8,1,2,3,4/)

     do ilevel = 2, nlevelmax
        do icpu = 1, ncpu
           igrid = headl(icpu, ilevel)
           do while(igrid > 0)
              ! Extract father grid and oct position from father cell index
              ind_oct_nb = ICHILD_OF(father(igrid))
              igrid_par_nb = IGRID_OF(father(igrid))

              do j_nb = 1, twondim
                 ig_nb = ggg_nb(ind_oct_nb, j_nb)
                 ih_nb = hhh_nb(ind_oct_nb, j_nb)
                 if(ig_nb == 0) then
                    nbor_grid = igrid_par_nb  ! Central grid (parent itself)
                 else
                    if(nbor(igrid_par_nb, ig_nb) > 0) then
                       nbor_grid = son(nbor(igrid_par_nb, ig_nb))
                    else
                       nbor_grid = 0
                    end if
                 end if
                 if(nbor_grid > 0) then
                    nbor(igrid, j_nb) = ICELL_OF(nbor_grid,ih_nb)
                 else
                    nbor(igrid, j_nb) = 0
                 end if
              end do

              igrid = next(igrid)
           end do
        end do
     end do

     if(myid==1) write(*,*) 'HDF5: nbor() rebuild completed'
  end if

  !=====================================================
  ! Step 8: Compute numbtot via MPI_ALLREDUCE
  !=====================================================
  do ilevel = 1, nlevelmax
     numbtot(1:4, ilevel) = 0
     ! numbtot(1,ilevel) = total active grids across all CPUs
     call MPI_ALLREDUCE(numbl(myid, ilevel), numbtot(1, ilevel), 1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  end do

  !=====================================================
  ! Step 9: Build communicators for each level
  !=====================================================
  do ilevel = 1, nlevelmax
     call build_comm(ilevel)
  end do
  if(myid==1) write(*,*) 'HDF5: build_comm completed for all levels'

  !=====================================================
  ! Step 10: Close HDF5 file
  !=====================================================
  call hdf5_close_file()

  if(myid==1) write(*,*) 'HDF5 AMR restore done.'

end subroutine restore_amr_hdf5

!###########################################################################
! Hydro restore from HDF5
!###########################################################################
subroutine restore_hydro_hdf5()
#ifdef SNRT
  use snrt_hdf5, only: snrt_hdf5_read
#endif
  use amr_commons
  use hydro_commons
#ifdef DUST_LIVE
  use hydro_parameters, only: idust, idust_energy
#endif
  use ramses_hdf5_io
  ! varcpu variables now in amr_commons
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, ivar, info
  integer :: hdf5_attr_status, hdf5_attr_status_all
  integer :: ngrid_loc, nvar_file, fidx
#ifdef DUST_LIVE
  integer :: dust_mass_field_file, dust_energy_field_file
#endif
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: ubuf(:), ubuf_all(:), ubuf_chunk(:)
  integer(HID_T) :: grp_id, lvl_grp_id, hdr_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, var_str
  character(len=5) :: nchar
  ! Streaming chunk buffers for varcpu hydro restore
  integer :: chunk_size, this_chunk, local_idx, ngrid_file_lvl
  integer :: n_chunk_grids, igrid_w, j
  integer, allocatable :: chunk_igrids(:), chunk_fidx_offset(:)
  integer(i8b) :: i_global

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore hydro: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Read header to get ncpu_file and nvar_file
  call hdf5_open_group('/header', hdr_grp_id)
  call hdf5_read_attr_int_checked(hdr_grp_id, 'ncpu', ncpu_file, hdf5_attr_status)
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0 .or. ncpu_file <= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/ncpu status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  call hdf5_read_attr_int_checked(hdr_grp_id, 'nvar', nvar_file, hdf5_attr_status)
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/nvar status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
#ifdef DUST_LIVE
  call hdf5_read_attr_int_checked(hdr_grp_id, 'dust_mass_field_index', &
       dust_mass_field_file, hdf5_attr_status)
  if(hdf5_attr_status == 0) then
     if(dust_mass_field_file /= idust) hdf5_attr_status = 1
  end if
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: DUST_LIVE checkpoint dust mass field missing or mismatched'
     call hdf5_restart_abort
  end if
  call hdf5_read_attr_int_checked(hdr_grp_id, 'dust_energy_field_index', &
       dust_energy_field_file, hdf5_attr_status)
  if(hdf5_attr_status == 0) then
     if(dust_energy_field_file /= idust_energy .or. nvar_file /= nvar) hdf5_attr_status = 1
  end if
  call MPI_Allreduce(hdf5_attr_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, info)
  if(hdf5_attr_status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: DUST_LIVE checkpoint dust energy field or NVAR missing or mismatched'
     call hdf5_restart_abort
  end if
#endif
  call hdf5_close_group(hdr_grp_id)

  if(nvar_file /= nvar) then
     if(myid==1) write(*,*) 'WARNING: HDF5 nvar mismatch, file=', nvar_file, ' expected=', nvar
  end if

  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: streaming chunked read + scatter.
     ! Active grids on each rank are in file-index order (Phase 1 added them
     ! sequentially), so a single forward walk of headl→next suffices per
     ! chunk: advance pointer while fidx <= chunk_end, exit at first fidx > end.
     !===================================================
     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        write(lvl_str, '(I0)') ilevel
        grp_name = '/hydro/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        ngrid_file_lvl = varcpu_ngrid_file(ilevel)
        chunk_size = min(1048576, ngrid_file_lvl)
        allocate(ubuf_chunk(chunk_size * twotondim))
        ! Per-chunk grid index buffers (max possible = ngrid_loc owned at this level)
        allocate(chunk_igrids(max(ngrid_loc, 1)))
        allocate(chunk_fidx_offset(max(ngrid_loc, 1)))

        ! Outer chunk loop: walk linked list ONCE per chunk, read all ivars,
        ! parallelize scatter writes over chunk_igrids (independent igrids).
        igrid = headl(myid, ilevel)
        i_global = 0_i8b
        do while (i_global < int(ngrid_file_lvl, i8b))
           this_chunk = int(min(int(chunk_size, i8b), &
                int(ngrid_file_lvl, i8b) - i_global))

           ! Collect active grids whose fidx falls in (i_global, i_global+this_chunk]
           n_chunk_grids = 0
           do while(igrid > 0)
              fidx = varcpu_grid_file_idx(igrid)
              if(fidx > int(i_global) + this_chunk) exit
              n_chunk_grids = n_chunk_grids + 1
              chunk_igrids(n_chunk_grids) = igrid
              chunk_fidx_offset(n_chunk_grids) = fidx - int(i_global)
              igrid = next(igrid)
           end do

           do ivar = 1, min(nvar, nvar_file)
              write(var_str, '(I0)') ivar
              call hdf5_read_dataset_chunk_dp(lvl_grp_id, &
                   'uold_'//trim(var_str), ubuf_chunk, &
                   this_chunk * twotondim, i_global * int(twotondim, i8b))

!$OMP PARALLEL DO PRIVATE(j,igrid_w,local_idx,ind) SCHEDULE(STATIC)
              do j = 1, n_chunk_grids
                 igrid_w = chunk_igrids(j)
                 local_idx = chunk_fidx_offset(j)
                 do ind = 1, twotondim
                    uold(ICELL_OF(igrid_w,ind), ivar) = &
                         ubuf_chunk((local_idx-1)*twotondim + ind)
                 end do
              end do
!$OMP END PARALLEL DO
           end do
           i_global = i_global + int(this_chunk, i8b)
        end do

        deallocate(ubuf_chunk, chunk_igrids, chunk_fidx_offset)
        call hdf5_close_group(lvl_grp_id)
     end do
  else
     !===================================================
     ! SAME NCPU PATH: parallel I/O with per-rank offsets
     !===================================================
     allocate(ngrid_all(ncpu))

     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        call MPI_ALLGATHER(ngrid_loc, 1, MPI_INTEGER, ngrid_all, 1, MPI_INTEGER, &
             MPI_COMM_WORLD, info)
        ngrid_total = 0
        do i = 1, ncpu
           ngrid_total = ngrid_total + ngrid_all(i)
        end do
        if(ngrid_total == 0) cycle

        ncells_total = ngrid_total * twotondim
        offset_cells = 0
        do i = 1, myid - 1
           offset_cells = offset_cells + int(ngrid_all(i), i8b) * twotondim
        end do

        write(lvl_str, '(I0)') ilevel
        grp_name = '/hydro/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        if(ngrid_loc > 0) then
           allocate(ubuf(ngrid_loc * twotondim))
        else
           allocate(ubuf(1))
        end if

        ! Read uold: raw conservative variables
        do ivar = 1, min(nvar, nvar_file)
           write(var_str, '(I0)') ivar
           call hdf5_read_dataset_1d_dp(lvl_grp_id, 'uold_'//trim(var_str), &
                ubuf, ngrid_loc * twotondim, offset_cells)

           ! Scatter to uold array
           igrid = headl(myid, ilevel)
           do i = 1, ngrid_loc
              do ind = 1, twotondim
                 uold(ICELL_OF(igrid,ind), ivar) = ubuf((i-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        deallocate(ubuf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_all)
  end if

  call hdf5_restore_errors()
#ifdef SNRT
  call snrt_hdf5_read()
#endif
  call hdf5_close_file()

  ! Virtual exchange to populate ghost cells
  do ilevel = 1, nlevelmax
     if(numbtot(1, ilevel) == 0) cycle
     do ivar = 1, nvar
        call make_virtual_fine_dp(uold(1,ivar), ilevel)
     end do
  end do

  if(myid==1) write(*,*) 'HDF5 hydro restore done.'

end subroutine restore_hydro_hdf5

!###########################################################################
! Poisson restore from HDF5
!###########################################################################
subroutine restore_poisson_hdf5()
  use amr_commons
  use poisson_commons
  use ramses_hdf5_io
  ! varcpu variables now in amr_commons
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilevel, i, igrid, ind, idim, info, fidx
  integer :: ngrid_loc
  integer, allocatable :: ngrid_all(:)
  integer(i8b) :: ncells_total, offset_cells, ngrid_total
  real(dp), allocatable :: pbuf(:), pbuf_all(:), pbuf_chunk(:)
  integer(HID_T) :: lvl_grp_id
  character(len=200) :: h5filename
  character(len=40) :: grp_name
  character(len=10) :: lvl_str, dim_str
  character(len=5) :: nchar
  ! Streaming chunk buffers for varcpu poisson restore
  integer :: chunk_size, this_chunk, local_idx, ngrid_file_lvl
  integer :: n_chunk_grids, igrid_w, j
  integer, allocatable :: chunk_igrids(:), chunk_fidx_offset(:)
  integer(i8b) :: i_global
  logical :: scalar_dset_exists, scalar_restored_any
  integer :: h5err_scalar

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore poisson: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Suppress HDF5 error messages for missing level groups
  call hdf5_suppress_errors()
  scalar_restored_any=.false.

  if(varcpu_restart) then
     !===================================================
     ! VARIABLE NCPU PATH: streaming chunked read + scatter.
     ! Same pattern as hydro: active grids in file-index order, single
     ! forward walk per chunk for each component (phi + ndim force comps).
     !===================================================
     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        write(lvl_str, '(I0)') ilevel
        grp_name = '/gravity/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block
        scalar_dset_exists=.false.
        if(allocated(scalar_gr)) then
           call h5lexists_f(lvl_grp_id,'scalar_gr',scalar_dset_exists,h5err_scalar)
           scalar_dset_exists=scalar_dset_exists.and.h5err_scalar.eq.0
           scalar_restored_any=scalar_restored_any.or.scalar_dset_exists
        end if

        ngrid_file_lvl = varcpu_ngrid_file(ilevel)
        chunk_size = min(1048576, ngrid_file_lvl)
        allocate(pbuf_chunk(chunk_size * twotondim))
        allocate(chunk_igrids(max(ngrid_loc, 1)))
        allocate(chunk_fidx_offset(max(ngrid_loc, 1)))

        ! Single chunk loop: walk linked list ONCE per chunk, then read phi
        ! and each force component, parallel-scatter on the same igrid set.
        igrid = headl(myid, ilevel)
        i_global = 0_i8b
        do while (i_global < int(ngrid_file_lvl, i8b))
           this_chunk = int(min(int(chunk_size, i8b), &
                int(ngrid_file_lvl, i8b) - i_global))

           n_chunk_grids = 0
           do while(igrid > 0)
              fidx = varcpu_grid_file_idx(igrid)
              if(fidx > int(i_global) + this_chunk) exit
              n_chunk_grids = n_chunk_grids + 1
              chunk_igrids(n_chunk_grids) = igrid
              chunk_fidx_offset(n_chunk_grids) = fidx - int(i_global)
              igrid = next(igrid)
           end do

           ! phi
           call hdf5_read_dataset_chunk_dp(lvl_grp_id, 'phi', pbuf_chunk, &
                this_chunk * twotondim, i_global * int(twotondim, i8b))
!$OMP PARALLEL DO PRIVATE(j,igrid_w,local_idx,ind) SCHEDULE(STATIC)
           do j = 1, n_chunk_grids
              igrid_w = chunk_igrids(j)
              local_idx = chunk_fidx_offset(j)
              do ind = 1, twotondim
                 phi(ICELL_OF(igrid_w,ind)) = pbuf_chunk((local_idx-1)*twotondim + ind)
              end do
           end do
!$OMP END PARALLEL DO

           ! force components
           do idim = 1, ndim
              write(dim_str, '(I0)') idim
              call hdf5_read_dataset_chunk_dp(lvl_grp_id, &
                   'f_'//trim(dim_str), pbuf_chunk, &
                   this_chunk * twotondim, i_global * int(twotondim, i8b))
!$OMP PARALLEL DO PRIVATE(j,igrid_w,local_idx,ind) SCHEDULE(STATIC)
              do j = 1, n_chunk_grids
                 igrid_w = chunk_igrids(j)
                 local_idx = chunk_fidx_offset(j)
                 do ind = 1, twotondim
                    f(ICELL_OF(igrid_w,ind), idim) = &
                         pbuf_chunk((local_idx-1)*twotondim + ind)
                 end do
              end do
!$OMP END PARALLEL DO
           end do

           if(scalar_dset_exists) then
              call hdf5_read_dataset_chunk_dp(lvl_grp_id, 'scalar_gr', &
                   pbuf_chunk, this_chunk * twotondim, &
                   i_global * int(twotondim, i8b))
!$OMP PARALLEL DO PRIVATE(j,igrid_w,local_idx,ind) SCHEDULE(STATIC)
              do j = 1, n_chunk_grids
                 igrid_w = chunk_igrids(j)
                 local_idx = chunk_fidx_offset(j)
                 do ind = 1, twotondim
                    scalar_gr(ICELL_OF(igrid_w,ind)) = &
                         pbuf_chunk((local_idx-1)*twotondim + ind)
                    scalar_gr_old(ICELL_OF(igrid_w,ind)) = scalar_gr(ICELL_OF(igrid_w,ind))
                 end do
              end do
!$OMP END PARALLEL DO
           end if
           i_global = i_global + int(this_chunk, i8b)
        end do

        deallocate(pbuf_chunk, chunk_igrids, chunk_fidx_offset)
        call hdf5_close_group(lvl_grp_id)
     end do

     ! Free mapping arrays — no longer needed
     if(allocated(varcpu_grid_file_idx)) deallocate(varcpu_grid_file_idx)
     if(allocated(varcpu_ngrid_file)) deallocate(varcpu_ngrid_file)
  else
     !===================================================
     ! SAME NCPU PATH: parallel I/O with per-rank offsets
     !===================================================
     allocate(ngrid_all(ncpu))

     do ilevel = 1, nlevelmax
        ngrid_loc = numbl(myid, ilevel)
        call MPI_ALLGATHER(ngrid_loc, 1, MPI_INTEGER, ngrid_all, 1, MPI_INTEGER, &
             MPI_COMM_WORLD, info)
        ngrid_total = 0
        do i = 1, ncpu
           ngrid_total = ngrid_total + ngrid_all(i)
        end do
        if(ngrid_total == 0) cycle

        ncells_total = ngrid_total * twotondim
        offset_cells = 0
        do i = 1, myid - 1
           offset_cells = offset_cells + int(ngrid_all(i), i8b) * twotondim
        end do

        write(lvl_str, '(I0)') ilevel
        grp_name = '/gravity/level_'//trim(lvl_str)

        block
           integer :: h5err
           logical :: grp_exists
           call h5lexists_f(hdf5_file_id, trim(grp_name), grp_exists, h5err)
           if(.not. grp_exists .or. h5err /= 0) cycle
           call h5gopen_f(hdf5_file_id, trim(grp_name), lvl_grp_id, h5err)
           if(h5err /= 0) cycle
        end block

        if(ngrid_loc > 0) then
           allocate(pbuf(ngrid_loc * twotondim))
        else
           allocate(pbuf(1))
        end if

        ! Read phi
        call hdf5_read_dataset_1d_dp(lvl_grp_id, 'phi', &
             pbuf, ngrid_loc * twotondim, offset_cells)
        igrid = headl(myid, ilevel)
        do i = 1, ngrid_loc
           do ind = 1, twotondim
              phi(ICELL_OF(igrid,ind)) = pbuf((i-1)*twotondim + ind)
           end do
           igrid = next(igrid)
        end do

        ! Read force
        do idim = 1, ndim
           write(dim_str, '(I0)') idim
           call hdf5_read_dataset_1d_dp(lvl_grp_id, 'f_'//trim(dim_str), &
                pbuf, ngrid_loc * twotondim, offset_cells)
           igrid = headl(myid, ilevel)
           do i = 1, ngrid_loc
              do ind = 1, twotondim
                 f(ICELL_OF(igrid,ind), idim) = pbuf((i-1)*twotondim + ind)
              end do
              igrid = next(igrid)
           end do
        end do

        ! Read the scalar field when present.  Missing datasets identify
        ! legacy checkpoints and intentionally leave zero sentinels for
        ! the model-specific initialization routines.
        if(allocated(scalar_gr)) then
           scalar_dset_exists=.false.
           call h5lexists_f(lvl_grp_id,'scalar_gr',scalar_dset_exists,h5err_scalar)
           if(scalar_dset_exists .and. h5err_scalar.eq.0) then
              scalar_restored_any=.true.
              call hdf5_read_dataset_1d_dp(lvl_grp_id, 'scalar_gr', &
                   pbuf, ngrid_loc * twotondim, offset_cells)
              igrid = headl(myid, ilevel)
              do i = 1, ngrid_loc
                 do ind = 1, twotondim
                    scalar_gr(ICELL_OF(igrid,ind)) = pbuf((i-1)*twotondim + ind)
                    scalar_gr_old(ICELL_OF(igrid,ind)) = scalar_gr(ICELL_OF(igrid,ind))
                 end do
                 igrid = next(igrid)
              end do
           end if
        end if

        ! Read FDM psi if enabled
        if(use_fdm) then
           block
              logical :: dset_exists
              integer :: h5err_dset
              call h5lexists_f(lvl_grp_id, 'psi_re', dset_exists, h5err_dset)
              if(dset_exists) then
                 call hdf5_read_dataset_1d_dp(lvl_grp_id, 'psi_re', &
                      pbuf, ngrid_loc * twotondim, offset_cells)
                 igrid = headl(myid, ilevel)
                 do i = 1, ngrid_loc
                    do ind = 1, twotondim
                       psi_re(ICELL_OF(igrid,ind)) = pbuf((i-1)*twotondim + ind)
                    end do
                    igrid = next(igrid)
                 end do
                 call hdf5_read_dataset_1d_dp(lvl_grp_id, 'psi_im', &
                      pbuf, ngrid_loc * twotondim, offset_cells)
                 igrid = headl(myid, ilevel)
                 do i = 1, ngrid_loc
                    do ind = 1, twotondim
                       psi_im(ICELL_OF(igrid,ind)) = pbuf((i-1)*twotondim + ind)
                    end do
                    igrid = next(igrid)
                 end do
              end if
           end block
        end if

        deallocate(pbuf)
        call hdf5_close_group(lvl_grp_id)
     end do

     deallocate(ngrid_all)
  end if

  call hdf5_restore_errors()
  call hdf5_close_file()

  ! Virtual exchange to populate ghost cells
  do ilevel = 1, nlevelmax
     if(numbtot(1, ilevel) == 0) cycle
     call make_virtual_fine_dp(phi(1), ilevel)
     do idim = 1, ndim
        call make_virtual_fine_dp(f(1,idim), ilevel)
     end do
     if(use_fdm) then
        call make_virtual_fine_dp(psi_re(1), ilevel)
        call make_virtual_fine_dp(psi_im(1), ilevel)
     end if
     if(allocated(scalar_gr)) then
        call make_virtual_fine_dp(scalar_gr(1), ilevel)
        call make_virtual_fine_dp(scalar_gr_old(1), ilevel)
     end if
  end do

  if(myid==1 .and. scalar_restored_any) &
       & write(*,*)'Modified-gravity scalar field restored from HDF5 checkpoint'
  if(myid==1 .and. allocated(scalar_gr) .and. .not.scalar_restored_any) &
       & write(*,*)'Scalar field absent from legacy HDF5 checkpoint; solver will initialize it'
  if(myid==1) write(*,*) 'HDF5 poisson restore done.'

end subroutine restore_poisson_hdf5

!###########################################################################
! Particle restore from HDF5
!###########################################################################
subroutine restore_part_hdf5()
  use amr_commons
  use pm_commons
  use stellar_enrichment_config, only: stellar_hdf5_state_schema_version
  use ramses_hdf5_io
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: i, idim, info
  integer :: npart_loc
  integer(i8b) :: npart_total_file, offset_part, tmp_long
  integer, allocatable :: npart_per_cpu(:)
  real(dp), allocatable :: dbuf(:)
  integer(i8b), allocatable :: ibuf8(:)
  integer, allocatable :: ibuf(:)
  integer(HID_T) :: grp_id, hdr_grp_id
  character(len=200) :: h5filename
  character(len=10) :: dstr
  character(len=5) :: nchar
  logical :: dset_exists
  logical :: attr_exists
  integer :: stellar_state_schema_version
  integer :: h5err
  integer :: stellar_state_read_status
  integer :: hdf5_attr_status_all
  integer :: nindsink_file

  call title(nrestart, nchar)
  h5filename = 'output_'//trim(nchar)//'/data_'//trim(nchar)//'.h5'

  if(myid==1) write(*,*) 'HDF5 restore particles: ', trim(h5filename)

  call hdf5_open_parallel(h5filename, MPI_COMM_WORLD)

  ! Read header
  call hdf5_open_group('/header', hdr_grp_id)
  call hdf5_read_attr_int_checked(hdr_grp_id, 'ncpu', ncpu_file, &
       stellar_state_read_status)
  call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, h5err)
  if(hdf5_attr_status_all /= 0 .or. ncpu_file <= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/ncpu status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  call hdf5_close_group(hdr_grp_id)

  ! Read particle group
  call hdf5_open_group('/particles', grp_id)
  call hdf5_read_attr_int8_checked(grp_id, 'npart_total', npart_total_file, &
       stellar_state_read_status)
  call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, h5err)
  if(hdf5_attr_status_all /= 0 .or. npart_total_file < 0_i8b) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for particles/npart_total status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  block
     integer(i8b) :: nstar_tot_tmp
     call hdf5_read_attr_int8_checked(grp_id, 'nstar_tot', nstar_tot_tmp, &
          stellar_state_read_status)
     call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, h5err)
     if(hdf5_attr_status_all /= 0 .or. nstar_tot_tmp < 0_i8b) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for particles/nstar_tot status=', &
             hdf5_attr_status_all
        call hdf5_restart_abort
     end if
     nstar_tot = nstar_tot_tmp
  end block
  call hdf5_read_attr_dp_checked(grp_id, 'mstar_tot', mstar_tot, stellar_state_read_status)
  call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, h5err)
  if(hdf5_attr_status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for particles/mstar_tot status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if
  call hdf5_read_attr_dp_checked(grp_id, 'mstar_lost', mstar_lost, stellar_state_read_status)
  call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
       MPI_MAX, MPI_COMM_WORLD, h5err)
  if(hdf5_attr_status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for particles/mstar_lost status=', &
          hdf5_attr_status_all
     call hdf5_restart_abort
  end if

  ! Read npart_per_cpu
  allocate(npart_per_cpu(ncpu_file))
  call hdf5_read_dataset_all_int_checked(grp_id, 'npart_per_cpu', npart_per_cpu, &
       ncpu_file, stellar_state_read_status)
  if(stellar_state_read_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for npart_per_cpu status=', &
          stellar_state_read_status
     call hdf5_restart_abort
  end if
  tmp_long = 0_i8b
  do i = 1, ncpu_file
     if(npart_per_cpu(i) < 0) then
        if(myid==1) write(*,*) 'ERROR: invalid npart_per_cpu: negative entry at rank=', i
        call hdf5_restart_abort
     end if
     tmp_long = tmp_long + int(npart_per_cpu(i), i8b)
  end do
  if(tmp_long /= npart_total_file) then
     if(myid==1) write(*,*) 'ERROR: invalid npart_per_cpu: sum=', tmp_long, &
          ' expected=', npart_total_file
     call hdf5_restart_abort
  end if

  ! Determine local particle count and offset
  if(ncpu_file == ncpu) then
     ! Same rank: direct mapping
     npart_loc = npart_per_cpu(myid)
     offset_part = 0
     do i = 1, myid - 1
        offset_part = offset_part + npart_per_cpu(i)
     end do
  else
     ! Different rank: distribute evenly
     block
        integer(i8b) :: npp, remainder
        npp = npart_total_file / ncpu
        remainder = npart_total_file - npp * ncpu
        if(myid <= int(remainder)) then
           npart_loc = int(npp + 1)
           offset_part = int(npp + 1, i8b) * (myid - 1)
        else
           npart_loc = int(npp)
           offset_part = (npp + 1) * remainder + npp * (myid - 1 - remainder)
        end if
     end block
  end if

  if(npart_loc > npartmax) then
     if(cosmo.and.pic.and.npartmax_auto)then
        call grow_particle_bundle(npart_loc)
     else
        write(*,*) 'ERROR: npart_loc > npartmax', npart_loc, npartmax
        call clean_stop
     endif
  end if

  write(*,'(A,I0,A,I0)') 'HDF5 stellar release local count rank=', myid, &
       ' : ', npart_loc

  ! HDF5 particles are restored contiguously into slots 1:npart_loc.
  ! Mark every remaining slot as free before any full-array occupancy scan.
  if(npart_loc < npartmax) levelp(npart_loc+1:npartmax) = 0

  ! Allocate buffers
  allocate(dbuf(npart_loc))
  allocate(ibuf8(npart_loc))
  allocate(ibuf(npart_loc))

  ! Read positions
  do idim = 1, ndim
     write(dstr, '(I0)') idim
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'x_'//trim(dstr), dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for x_'//trim(dstr)// &
             ' status=', stellar_state_read_status
        call hdf5_restart_abort
     end if
     xp(1:npart_loc, idim) = dbuf(1:npart_loc)
  end do

  ! Read velocities
  do idim = 1, ndim
     write(dstr, '(I0)') idim
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'v_'//trim(dstr), dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for v_'//trim(dstr)// &
             ' status=', stellar_state_read_status
        call hdf5_restart_abort
     end if
     vp(1:npart_loc, idim) = dbuf(1:npart_loc)
  end do

  ! Read mass
  call hdf5_read_dataset_1d_dp_checked(grp_id, 'mass', dbuf, npart_loc, &
       offset_part, npart_total_file, stellar_state_read_status)
  if(stellar_state_read_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for mass status=', &
          stellar_state_read_status
     call hdf5_restart_abort
  end if
  mp(1:npart_loc) = dbuf(1:npart_loc)

  ! Read identity
  call hdf5_read_dataset_1d_int8_checked(grp_id, 'identity', ibuf8, npart_loc, &
       offset_part, npart_total_file, stellar_state_read_status)
  if(stellar_state_read_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for identity status=', &
          stellar_state_read_status
     call hdf5_restart_abort
  end if
  idp(1:npart_loc) = ibuf8(1:npart_loc)

  ! Read level
  call hdf5_read_dataset_1d_int_checked(grp_id, 'levelp', ibuf, npart_loc, &
       offset_part, npart_total_file, stellar_state_read_status)
  if(stellar_state_read_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for levelp status=', &
          stellar_state_read_status
     call hdf5_restart_abort
  end if
  levelp(1:npart_loc) = ibuf(1:npart_loc)

  ! Read particle type (stored as int4 in HDF5 for convenience, cast to int1)
  call hdf5_read_dataset_1d_int_checked(grp_id, 'ptypep', ibuf, npart_loc, &
       offset_part, npart_total_file, stellar_state_read_status)
  if(stellar_state_read_status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for ptypep status=', &
          stellar_state_read_status
     call hdf5_restart_abort
  end if
  ptypep(1:npart_loc) = int(ibuf(1:npart_loc), kind=1)

  ! Read birth epoch and metallicity if star/sink
  if(star .or. sink) then
     ! A stellar restart without the release state is scientifically unsafe:
     ! silently retaining init_part's zero values would skip already released
     ! material or replay it.  Require the version marker and every field,
     ! while keeping a clear failure mode for legacy HDF5 checkpoints.
     call h5aexists_f(grp_id, 'stellar_state_schema_version', attr_exists, h5err)
     if(h5err /= 0 .or. .not. attr_exists) then
        if(myid==1) write(*,*) 'ERROR: HDF5 stellar restart lacks ', &
             'stellar_state_schema_version; refusing unsafe restore'
        call hdf5_restart_abort
     end if
     stellar_state_schema_version = 0
     call hdf5_read_attr_int_checked(grp_id, 'stellar_state_schema_version', &
          stellar_state_schema_version, stellar_state_read_status)
     call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, h5err)
     if(hdf5_attr_status_all /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for particles/stellar_state_schema_version status=', &
             hdf5_attr_status_all
        call hdf5_restart_abort
     end if
     if(stellar_state_schema_version /= stellar_hdf5_state_schema_version) then
        if(myid==1) write(*,*) 'ERROR: unsupported HDF5 stellar state schema=', &
             stellar_state_schema_version, ' expected=', &
             stellar_hdf5_state_schema_version
        call hdf5_restart_abort
     end if

     call hdf5_read_dataset_1d_dp_checked(grp_id, 'birth_epoch', dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for birth_epoch status=', &
             stellar_state_read_status
        call hdf5_restart_abort
     end if
     tp(1:npart_loc) = dbuf(1:npart_loc)
     if(metal) then
        call hdf5_read_dataset_1d_dp_checked(grp_id, 'metallicity', dbuf, npart_loc, &
             offset_part, npart_total_file, stellar_state_read_status)
        if(stellar_state_read_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for metallicity status=', &
                stellar_state_read_status
           call hdf5_restart_abort
        end if
        zp(1:npart_loc) = dbuf(1:npart_loc)
     end if

     call h5lexists_f(grp_id, 'tpp', dset_exists, h5err)
     if(h5err /= 0 .or. .not. dset_exists) then
        if(myid==1) write(*,*) 'ERROR: HDF5 stellar restart lacks tpp'
        call hdf5_restart_abort
     end if
     call h5lexists_f(grp_id, 'mp0', dset_exists, h5err)
     if(h5err /= 0 .or. .not. dset_exists) then
        if(myid==1) write(*,*) 'ERROR: HDF5 stellar restart lacks mp0'
        call hdf5_restart_abort
     end if
     call h5lexists_f(grp_id, 'indtab', dset_exists, h5err)
     if(h5err /= 0 .or. .not. dset_exists) then
        if(myid==1) write(*,*) 'ERROR: HDF5 stellar restart lacks indtab'
        call hdf5_restart_abort
     end if
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'tpp', dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for tpp status=', &
             stellar_state_read_status
        call hdf5_restart_abort
     end if
     tpp(1:npart_loc) = dbuf(1:npart_loc)
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'mp0', dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for mp0 status=', &
             stellar_state_read_status
        call hdf5_restart_abort
     end if
     mp0(1:npart_loc) = dbuf(1:npart_loc)
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'indtab', dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for indtab status=', &
             stellar_state_read_status
        call hdf5_restart_abort
     end if
     indtab(1:npart_loc) = dbuf(1:npart_loc)
     if(myid==1 .and. npart_loc>0) then
        write(*,'(A,3(1X,ES24.16))') 'HDF5 stellar release state sample:', &
             tpp(1), mp0(1), indtab(1)
     end if
     if(npart_loc>0) then
        write(*,'(A,I0,A,I0,3(1X,ES24.16))') &
             'HDF5 stellar release local state rank=', myid, ' : ', npart_loc, &
             tpp(1), mp0(1), indtab(1)
     else
        write(*,'(A,I0,A,I0,A)') &
             'HDF5 stellar release local state rank=', myid, ' : ', npart_loc, &
             ' EMPTY'
     end if
  end if

  ! Atomic Dark Matter: dark internal energy (optional dataset)
  if(use_adm) then
     call hdf5_read_dataset_1d_dp_checked(grp_id, 'dark_energy_int', dbuf, npart_loc, &
          offset_part, npart_total_file, stellar_state_read_status)
     if(stellar_state_read_status /= 0) then
        if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for dark_energy_int status=', &
             stellar_state_read_status
        call hdf5_restart_abort
     end if
     edp(1:npart_loc) = dbuf(1:npart_loc)

     ! Non-equilibrium dark-H2 fraction (Phase 2; absent in older checkpoints,
     ! in which case xh2p keeps its primordial adm_fH2 seed from init_part).
     call h5lexists_f(grp_id, 'dark_h2_frac', dset_exists, h5err)
     if(dset_exists) then
        call hdf5_read_dataset_1d_dp_checked(grp_id, 'dark_h2_frac', dbuf, npart_loc, &
             offset_part, npart_total_file, stellar_state_read_status)
        if(stellar_state_read_status /= 0) then
           if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for dark_h2_frac status=', &
                stellar_state_read_status
           call hdf5_restart_abort
        end if
        xh2p(1:npart_loc) = dbuf(1:npart_loc)
     end if
  end if

  npart = npart_loc

  deallocate(dbuf, ibuf8, ibuf, npart_per_cpu)
  call hdf5_close_group(grp_id)

  ! Read sinks
  if(sink) then
     call hdf5_open_group('/sinks', grp_id)
     nsink = 0
     nindsink_file = 0
     call hdf5_read_attr_int_checked(grp_id, 'nsink', nsink, stellar_state_read_status)
     call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, h5err)
     if(hdf5_attr_status_all /= 0 .or. nsink < 0 .or. nsink > nsinkmax) then
        if(myid==1) write(*,*) 'ERROR: invalid HDF5 nsink=', nsink, &
             ' nsinkmax=', nsinkmax, ' status=', hdf5_attr_status_all
        call hdf5_restart_abort
     end if
     call hdf5_read_attr_int_checked(grp_id, 'nindsink', nindsink_file, stellar_state_read_status)
     call MPI_Allreduce(stellar_state_read_status, hdf5_attr_status_all, 1, MPI_INTEGER, &
          MPI_MAX, MPI_COMM_WORLD, h5err)
     if(hdf5_attr_status_all /= 0 .or. nindsink_file < 0) then
        if(myid==1) write(*,*) 'ERROR: invalid HDF5 nindsink status=', hdf5_attr_status_all
        call hdf5_restart_abort
     end if

     if(nsink > 0) then
        block
           real(dp), allocatable :: sbuf(:)
           integer, allocatable :: sibuf(:)
           integer :: ilevel
           character(len=20) :: stat_name

           allocate(sibuf(nsink))
           call hdf5_restore_sink_int_checked(grp_id, 'idsink', sibuf, nsink)
           idsink(1:nsink) = sibuf(1:nsink)
           if(nindsink_file < maxval(idsink(1:nsink))) then
              if(myid==1) write(*,*) 'ERROR: HDF5 nindsink is below idsink maximum: ', &
                   nindsink_file, maxval(idsink(1:nsink))
              call hdf5_restart_abort
           end if
           nindsink = nindsink_file
           deallocate(sibuf)

           allocate(sbuf(nsink))
           call hdf5_restore_sink_dp_checked(grp_id, 'msink', sbuf, nsink)
           msink(1:nsink) = sbuf(1:nsink)

           ! Position
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_restore_sink_dp_checked(grp_id, 'xsink_'//trim(dstr), sbuf, nsink)
              xsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! Velocity
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_restore_sink_dp_checked(grp_id, 'vsink_'//trim(dstr), sbuf, nsink)
              vsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! Time
           call hdf5_restore_sink_dp_checked(grp_id, 'tsink', sbuf, nsink)
           tsink(1:nsink) = sbuf(1:nsink)
           ! Accretion
           call hdf5_restore_sink_dp_checked(grp_id, 'dMsmbh', sbuf, nsink)
           dMsmbh(1:nsink) = sbuf(1:nsink)
           call hdf5_restore_sink_dp_checked(grp_id, 'dMBH_coarse', sbuf, nsink)
           dMBH_coarse(1:nsink) = sbuf(1:nsink)
           call hdf5_restore_sink_dp_checked(grp_id, 'dMEd_coarse', sbuf, nsink)
           dMEd_coarse(1:nsink) = sbuf(1:nsink)
           ! Esave
           call hdf5_restore_sink_dp_checked(grp_id, 'Esave', sbuf, nsink)
           Esave(1:nsink) = sbuf(1:nsink)
           ! Gas angular momentum (jsink)
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_restore_sink_dp_checked(grp_id, 'jsink_'//trim(dstr), sbuf, nsink)
              jsink(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! BH spin axis
           do idim = 1, ndim
              write(dstr, '(I0)') idim
              call hdf5_restore_sink_dp_checked(grp_id, 'bhspin_'//trim(dstr), sbuf, nsink)
              bhspin(1:nsink, idim) = sbuf(1:nsink)
           end do
           ! BH spin magnitude
           call hdf5_restore_sink_dp_checked(grp_id, 'spinmag', sbuf, nsink)
           spinmag(1:nsink) = sbuf(1:nsink)
           ! BH efficiency
           call hdf5_restore_sink_dp_checked(grp_id, 'eps_sink', sbuf, nsink)
           eps_sink(1:nsink) = sbuf(1:nsink)
           ! Sink statistics
           do idim = 1, ndim*2+1
              do ilevel = levelmin, nlevelmax
                 write(stat_name, '(I0,"_",I0)') idim, ilevel
                 call hdf5_restore_sink_dp_checked(grp_id, 'sink_stat_'//trim(stat_name), sbuf, nsink)
                 sink_stat(1:nsink, ilevel, idim) = sbuf(1:nsink)
              end do
           end do
           deallocate(sbuf)
        end block
     end if

     call hdf5_close_group(grp_id)
  end if

  call hdf5_close_file()

  if(myid==1) write(*,*) 'HDF5 particle restore done.'

end subroutine restore_part_hdf5

subroutine hdf5_restore_header_int_checked(grp_id, name, val)
  use amr_commons, only: myid
  use ramses_hdf5_io
  use hdf5, only: HID_T
  implicit none
  include 'mpif.h'
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  integer, intent(out) :: val
  integer :: status, status_all, ierr

  call hdf5_read_attr_int_checked(grp_id, name, val, status)
  call MPI_Allreduce(status, status_all, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if(status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/', trim(name), &
          ' status=', status_all
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_header_int_checked

subroutine hdf5_restore_header_dp_checked(grp_id, name, val)
  use amr_commons, only: myid
  use amr_parameters, only: dp
  use ramses_hdf5_io
  use hdf5, only: HID_T
  implicit none
  include 'mpif.h'
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  real(dp), intent(out) :: val
  integer :: status, status_all, ierr

  call hdf5_read_attr_dp_checked(grp_id, name, val, status)
  call MPI_Allreduce(status, status_all, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if(status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/', trim(name), &
          ' status=', status_all
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_header_dp_checked

subroutine hdf5_restore_header_dp_array_checked(grp_id, name, val, n)
  use amr_commons, only: myid
  use amr_parameters, only: dp
  use ramses_hdf5_io
  use hdf5, only: HID_T
  implicit none
  include 'mpif.h'
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n
  real(dp), intent(out) :: val(n)
  integer :: status, status_all, ierr

  call hdf5_read_attr_1d_dp_checked(grp_id, name, val, n, status)
  call MPI_Allreduce(status, status_all, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if(status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/', trim(name), &
          ' status=', status_all
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_header_dp_array_checked

subroutine hdf5_restore_header_string_checked(grp_id, name, val)
  use amr_commons, only: myid
  use ramses_hdf5_io
  use hdf5, only: HID_T
  implicit none
  include 'mpif.h'
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: val
  integer :: status, status_all, ierr

  call hdf5_read_attr_string_checked(grp_id, name, val, status)
  call MPI_Allreduce(status, status_all, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
  if(status_all /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for header/', trim(name), &
          ' status=', status_all
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_header_string_checked

subroutine hdf5_restore_sink_dp_checked(grp_id, name, data, n)
  use amr_commons, only: myid
  use amr_parameters, only: dp
  use ramses_hdf5_io
  implicit none
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n
  real(dp), intent(out) :: data(n)
  integer :: status

  call hdf5_read_dataset_all_dp_checked(grp_id, name, data, n, status)
  if(status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for sink dataset ', &
          trim(name), ' status=', status
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_sink_dp_checked

subroutine hdf5_restore_sink_int_checked(grp_id, name, data, n)
  use amr_commons, only: myid
  use ramses_hdf5_io
  implicit none
  integer(HID_T), intent(in) :: grp_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n
  integer, intent(out) :: data(n)
  integer :: status

  call hdf5_read_dataset_all_int_checked(grp_id, name, data, n, status)
  if(status /= 0) then
     if(myid==1) write(*,*) 'ERROR: failed checked HDF5 read for sink dataset ', &
          trim(name), ' status=', status
     call hdf5_restart_abort
  end if
end subroutine hdf5_restore_sink_int_checked

subroutine hdf5_restart_abort
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: ierr
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#else
  error stop 1
#endif
end subroutine hdf5_restart_abort

#endif

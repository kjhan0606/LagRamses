!################################################################
!################################################################
! Fuzzy Dark Matter (FDM) Schrödinger-Poisson time integrator
!
! DKD operator splitting per level:
!   1. Half Drift (kinetic): exp(-i hbar k^2 dt/4a^2) in k-space (base)
!                         or FD Laplacian (refined levels)
!   2. Full Kick (potential): exp(-i Phi dt / hbar) in real space
!   3. Half Drift (kinetic): same as step 1
!
! Base level (levelmin): FFT kinetic operator (spectral accuracy)
! Fine levels (>levelmin): Explicit FD kinetic operator
!################################################################
subroutine fdm_step(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  real(dp)::dt_loc,dt_cfl,dt_sub
  real(dp)::yw0,yw1
  integer::nsub_hjm,isub_hjm

  if(.not.use_fdm) return
  if(numbtot(1,ilevel)==0) return

  dt_loc = dtnew(ilevel)

  ! ---- Hybrid HJM: fluid solver on coarse levels ----
  if(fdm_use_hjm .and. ilevel < fdm_first_wave_level) then
     ! CFL guard: newly refined levels inherit parent dt which may violate CFL
     call fdm_hjm_local_cfl(ilevel, dt_cfl)
     nsub_hjm = 1
     if(dt_cfl > 0.0d0 .and. dt_loc > dt_cfl) then
        nsub_hjm = min(ceiling(dt_loc / dt_cfl), 10000)
        if(myid==1) then
           write(*,'(A,I2,A,I6,A,ES10.3,A,ES10.3)') &
              ' HJM_SUBCYCLE: lv=',ilevel,' nsub=',nsub_hjm, &
              ' dt_loc=',dt_loc,' dt_cfl=',dt_cfl
        end if
     end if
     dt_sub = dt_loc / dble(nsub_hjm)
     do isub_hjm = 1, nsub_hjm
        call fdm_hjm_step(ilevel, dt_sub)
     end do
     call timer('fdm-ghost','start')
     call make_virtual_fine_dp(psi_re(1), ilevel)
     call make_virtual_fine_dp(psi_im(1), ilevel)
     call timer('particles','start')
     return
  end if

  ! ---- Schrödinger wave solver (existing) ----
  if(fdm_split_order == 4) then
     yw1 = 1.0d0/(2.0d0 - 2.0d0**(1.0d0/3.0d0))
     yw0 = 1.0d0 - 2.0d0*yw1
     call fdm_drift(ilevel, 0.5d0*yw1*dt_loc)
     call timer('fdm-kick','start')
     call fdm_kick (ilevel,        yw1*dt_loc)
     call fdm_drift(ilevel, 0.5d0*(yw1+yw0)*dt_loc)
     call timer('fdm-kick','start')
     call fdm_kick (ilevel,        yw0*dt_loc)
     call fdm_drift(ilevel, 0.5d0*(yw0+yw1)*dt_loc)
     call timer('fdm-kick','start')
     call fdm_kick (ilevel,        yw1*dt_loc)
     call fdm_drift(ilevel, 0.5d0*yw1*dt_loc)
  else
     call fdm_drift(ilevel, 0.5d0*dt_loc)
     call timer('fdm-kick','start')
     call fdm_kick(ilevel, dt_loc)
     call fdm_drift(ilevel, 0.5d0*dt_loc)
  end if

  call timer('fdm-ghost','start')
  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  call timer('particles','start')

end subroutine fdm_step
!################################################################
!################################################################
! Kinetic operator (Drift): exp(-i hbar/(2a^2) nabla^2 dt)
!
! Base level: FFT (spectral)
! Fine levels: explicit FD with 7-point Laplacian
!################################################################
subroutine fdm_drift(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

#ifdef USE_FFTW
  if(ilevel == levelmin) then
     call fdm_drift_fft(ilevel, dt_half)
     return
  end if
#endif
  call fdm_drift_fd(ilevel, dt_half)

end subroutine fdm_drift
!################################################################
!################################################################
! FFT kinetic operator for base level (spectral accuracy)
! Applies exp(-i alpha k^2) to psi in Fourier space
! where alpha = hbar_code * dt / (2 a^2)
!################################################################
#ifdef USE_FFTW
subroutine fdm_drift_fft(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  use iso_c_binding
  use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  include 'fftw3.f03'

  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  ! Grid dimensions
  integer::fft_Nx, fft_Ny, fft_Nz, nx_loc
  integer(i8b)::N_total
  real(dp)::dx_fft, twopi, alpha
  real(dp)::denom, phase, cos_p, sin_p

  ! FFTW plans (cached) -- single complex-to-complex transform of psi.
  ! psi = psi_re + i*psi_im is ONE complex field; the drift multiplies its
  ! full spectrum by exp(i*phase). (The old two-R2C scheme was wrong: C2R
  ! re-imposes Hermitian symmetry, flipping exp(+i*phase)->exp(-i*phase) on
  ! the unstored half of k-space and making the result axis-dependent.)
  logical,save::fdm_fftw_init = .false.
  type(C_PTR),save::plan_fwd = C_NULL_PTR, plan_bwd = C_NULL_PTR
  integer,save::saved_Nx = 0

  ! Work arrays (cached)
  complex(C_DOUBLE_COMPLEX),allocatable,save::psi_3d(:), psi_hat(:)
  integer,allocatable,save::fft_map(:)

  ! Loop vars
  integer::igrid,ind,iskip,icell,ix,iy,iz,idx_3d,info
  integer(i8b)::idx_c
  real(dp)::cr, ci, new_re, new_im
  real(dp)::scale

  ! Grid dimensions
  fft_Nx = nx * 2**ilevel
  fft_Ny = ny * 2**ilevel
  fft_Nz = nz * 2**ilevel
  N_total = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz,i8b)
  dx_fft = 0.5d0**ilevel
  twopi = 2.0d0*acos(-1.0d0)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)

  ! Supercomoving time (dt_code = H0*dt/a^2) already absorbs the 1/a^2 of
  ! the comoving kinetic term: no aexp factor in code units.
  alpha = hbar_code * dt_half

#ifndef WITHOUTMPI
  ! Large-grid / multi-rank: distributed FFTW-MPI slab path (no full-grid
  ! ALLREDUCE). Mirrors the distributed Poisson solver; same 256^3 threshold.
  if(N_total >= 256_i8b**3 .and. ncpu > 1) then
     call fdm_drift_fft_distributed(ilevel, dt_half)
     return
  end if
#endif

  ! Initialize FFTW plans and work arrays (once)
  call timer('fdm-plan','start')
  if(.not.fdm_fftw_init .or. fft_Nx /= saved_Nx) then
     if(allocated(psi_3d))  deallocate(psi_3d, psi_hat, fft_map)
     allocate(psi_3d(1:N_total), psi_hat(1:N_total))
     allocate(fft_map(1:N_total))

     ! Build AMR cell -> flat 3D index map
     fft_map = 0
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           ix = int((xg(igrid,1) + (dble(iand(ind-1,1)) - 0.5d0)*dx_fft &
                - dble(icoarse_min)*scale) / dx_fft / scale) + 1
           iy = int((xg(igrid,2) + (dble(iand(ishft(ind-1,-1),1)) - 0.5d0)*dx_fft &
                - dble(jcoarse_min)*scale) / dx_fft / scale) + 1
           iz = int((xg(igrid,3) + (dble(iand(ishft(ind-1,-2),1)) - 0.5d0)*dx_fft &
                - dble(kcoarse_min)*scale) / dx_fft / scale) + 1
           ix = max(1, min(ix, fft_Nx))
           iy = max(1, min(iy, fft_Ny))
           iz = max(1, min(iz, fft_Nz))
           idx_3d = ix + (iy-1)*fft_Nx + (iz-1)*fft_Nx*fft_Ny
           fft_map(idx_3d) = icell
           igrid = next(igrid)
        end do
     end do

     ! Create complex-to-complex forward / backward plans (full spectrum)
     if(c_associated(plan_fwd)) call dfftw_destroy_plan(plan_fwd)
     if(c_associated(plan_bwd)) call dfftw_destroy_plan(plan_bwd)
     call dfftw_plan_dft_3d(plan_fwd, fft_Nx, fft_Ny, fft_Nz, &
          psi_3d, psi_hat, FFTW_FORWARD, FFTW_MEASURE)
     call dfftw_plan_dft_3d(plan_bwd, fft_Nx, fft_Ny, fft_Nz, &
          psi_hat, psi_3d, FFTW_BACKWARD, FFTW_MEASURE)

     saved_Nx = fft_Nx
     fdm_fftw_init = .true.
  end if

  ! Step 1: Gather complex psi from AMR to flat 3D array
  ! ALLREDUCE for small grids (each rank has subset of cells)
  call timer('fdm-pack','start')
  psi_3d = (0.0d0, 0.0d0)
  do idx_3d=1,N_total
     icell = fft_map(idx_3d)
     if(icell > 0) then
        psi_3d(idx_3d) = cmplx(psi_re(icell), psi_im(icell), C_DOUBLE_COMPLEX)
     end if
  end do
#ifndef WITHOUTMPI
  if(ncpu > 1) then
     call timer('fdm-allreduce','start')
     call MPI_ALLREDUCE(MPI_IN_PLACE, psi_3d, int(N_total), &
          MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, info)
  end if
#endif

  ! Step 2: Forward FFT (C2C) of the full complex field
  call timer('fdm-fft-fwd','start')
  call dfftw_execute_dft(plan_fwd, psi_3d, psi_hat)

  ! Step 3: Apply kinetic phase rotation exp(i*phase) to the FULL spectrum.
  ! phi(k) is even in k, so cos() of the raw array index gives the correct
  ! discrete-Laplacian eigenvalue per axis (negative-k wraparound handled by cos).
  call timer('fdm-kphase','start')
  ! psi_hat(k) -> psi_hat(k) * exp(i*phase),  phase = 0.5*alpha*denom/dx^2
  ! denom = 2*sum_d [cos(2pi k_d/N_d) - 1]  (discrete Laplacian eigenvalue).
  ! Index runs over the FULL complex spectrum (no half-storage), so the phase
  ! is applied uniformly to every mode -- correct for a complex field.
  !$omp parallel do private(iz,iy,ix,idx_c,denom,phase, &
  !$omp&  cos_p,sin_p,cr,ci,new_re,new_im) schedule(static)
  do iz=0,fft_Nz-1
     do iy=0,fft_Ny-1
        do ix=0,fft_Nx-1
           idx_c = int(ix+1,i8b) + int(iy,i8b)*int(fft_Nx,i8b) &
                + int(iz,i8b)*int(fft_Nx,i8b)*int(fft_Ny,i8b)

           denom = 2.0d0*( cos(twopi*dble(ix)/dble(fft_Nx)) &
                         + cos(twopi*dble(iy)/dble(fft_Ny)) &
                         + cos(twopi*dble(iz)/dble(fft_Nz)) - 3.0d0 )

           phase = 0.5d0 * alpha * denom / (dx_fft**2)
           cos_p = cos(phase)
           sin_p = sin(phase)

           cr = dble(psi_hat(idx_c)); ci = aimag(psi_hat(idx_c))
           new_re = cr*cos_p - ci*sin_p
           new_im = cr*sin_p + ci*cos_p
           psi_hat(idx_c) = cmplx(new_re, new_im, C_DOUBLE_COMPLEX)
        end do
     end do
  end do
  !$omp end parallel do

  ! Step 4: Inverse FFT (C2C)
  call timer('fdm-fft-bwd','start')
  call dfftw_execute_dft(plan_bwd, psi_hat, psi_3d)

  ! Step 5: Normalize (FFTW convention: unnormalized) and scatter back
  call timer('fdm-scatter','start')
  do idx_3d=1,N_total
     icell = fft_map(idx_3d)
     if(icell > 0) then
        psi_re(icell) = dble(psi_3d(idx_3d)) / dble(N_total)
        psi_im(icell) = aimag(psi_3d(idx_3d)) / dble(N_total)
     end if
  end do

end subroutine fdm_drift_fft
!################################################################
!################################################################
! Distributed FFTW-MPI kinetic drift for the base (uniform) level.
! psi = psi_re + i*psi_im is ONE complex field: sparse P2P slab scatter of
! (psi_re, psi_im) into a single complex slab buffer, full complex-to-complex
! FFT, pointwise phase rotation exp(i*phase) over the WHOLE spectrum, inverse
! C2C, sparse P2P gather back. NO full-grid MPI_ALLREDUCE; scales to large N^3.
! (A C2C transform of psi is the correct kinetic operator; the earlier
! two-R2C scheme was axis-dependent and wrong -- see fdm_drift_fft.)
!################################################################
subroutine fdm_drift_fft_distributed(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  use iso_c_binding
  use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  include 'fftw3-mpi.f03'

  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  ! Grid dimensions
  integer::fft_Nx, fft_Ny, fft_Nz
  integer(i8b)::N_total, N_complex
  real(dp)::dx_fft, dx2_fft, twopi, alpha

  ! Cached FFTW MPI state (single complex slab buffer: psi = re + i*im)
  logical,save::dist_init = .false.
  type(C_PTR),save::plan_fwd = C_NULL_PTR, plan_bwd = C_NULL_PTR
  integer,save::saved_Nx = 0, saved_Ny = 0, saved_Nz = 0
  type(C_PTR),save::p_c = C_NULL_PTR
  complex(C_DOUBLE_COMPLEX),pointer,save::cdata(:) => null()
  integer(C_INTPTR_T),save::loc_nx = 0, nx_start = 0   ! x-slab (input/output)
  integer(C_INTPTR_T),save::loc_ny = 0, ny_start = 0   ! y-slab (transposed spectrum)
  integer(C_INTPTR_T),save::alloc_local = 0
  integer,save::fftw_block = 0
  real(dp),allocatable,save::kfac(:)                   ! geometry-only k^2 factor

  ! Sparse P2P partner lists (cached)
  integer,allocatable,save::send_partners(:), recv_partners(:)
  integer,save::n_send = 0, n_recv = 0
  logical,save::partners_done = .false.
  real(dp),save::cached_xmin = -1.0d0, cached_xmax = -1.0d0
  integer,save::cached_block = 0
  real(dp),allocatable,save::cpubox_min(:), cpubox_max(:)

  ! Cached static cell->slab mapping (rebuilt only on grid/partition change).
  ! For the uniform levelmin=levelmax grid, Kx/idx_3d/dest_rank/slot/icell are
  ! all geometry-fixed; only psi values change between calls.
  integer,allocatable,save::pack_slot(:), pack_icell(:), pack_idx3d(:)
  logical,save::cache_done = .false.
  integer,save::cached_ngrid = -1, M_cache = 0

  ! Locals
  integer::igrid, ngrid_loc, igrid_amr, icell_amr
  integer::ind, iskip, ix, iy, iz, idx_3d, slab_real_size
  integer::Kx, Ky, Kz, kx_i, ky_i, kz_i
  integer::info, irank, ip, nreq, dest_rank, x_local
  integer::n_send_total, n_recv_total
  integer::m, slot
  integer::slab_x_lo, slab_x_hi
  integer(i8b)::idx_c
  real(dp)::denom, phase, cos_p, sin_p
  real(dp)::my_xmin, my_xmax, recv_pos_lo, recv_pos_hi
  real(dp)::cr, ci
  logical::overlap, local_changed, global_changed

  integer,allocatable::sendcounts(:), sdispls(:), recvcounts(:), rdispls(:)
  integer,allocatable::send_idx(:), reqs(:), mpi_stat(:,:)
  real(dp),allocatable::sendbuf(:), recvbuf(:)
  real(dp),allocatable::sendbuf2(:), recvbuf2(:)   ! reverse path: (re,im) only

  ! Grid dimensions
  fft_Nx = nx * 2**ilevel
  fft_Ny = ny * 2**ilevel
  fft_Nz = nz * 2**ilevel
  N_total = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz,i8b)
  dx_fft  = 0.5d0**ilevel
  dx2_fft = dx_fft * dx_fft
  twopi   = 2.0d0*acos(-1.0d0)
  alpha   = hbar_code * dt_half   ! supercomoving time: no aexp factor

  call timer('fdm-plan','start')

  ! --- FFTW init (idempotent; Poisson solver may have already done it) ---
  if(.not. dist_init) then
     info = fftw_init_threads()
     call fftw_plan_with_nthreads(int(omp_get_max_threads(), C_INT))
     call fftw_mpi_init()
     dist_init = .true.
  end if

  ! --- Plan / buffer creation on size change ---
  if(fft_Nx /= saved_Nx .or. fft_Ny /= saved_Ny .or. fft_Nz /= saved_Nz) then
     if(c_associated(plan_fwd)) call fftw_destroy_plan(plan_fwd)
     if(c_associated(plan_bwd)) call fftw_destroy_plan(plan_bwd)
     if(c_associated(p_c)) call fftw_free(p_c)
     saved_Nx = fft_Nx; saved_Ny = fft_Ny; saved_Nz = fft_Nz

     ! Full complex C2C with TRANSPOSED ordering: input/output are x-slab
     ! distributed (first dim Nx), but the forward FFT leaves the spectrum
     ! TRANSPOSED (Ny-distributed, logical order Ny,Nx,Nz) -- this skips the
     ! transpose-back communication, halving the FFT all-to-all. The phase is
     ! applied in that transposed layout; the backward plan (TRANSPOSED_IN)
     ! consumes it and restores the x-slab layout.
     alloc_local = fftw_mpi_local_size_3d_transposed( &
          int(fft_Nx, C_INTPTR_T), int(fft_Ny, C_INTPTR_T), &
          int(fft_Nz, C_INTPTR_T), MPI_COMM_WORLD, &
          loc_nx, nx_start, loc_ny, ny_start)

     p_c = fftw_alloc_complex(alloc_local)
     call c_f_pointer(p_c, cdata, [int(alloc_local)])

     ! In-place complex forward/backward plans (transposed out/in).
     plan_fwd = fftw_mpi_plan_dft_3d( &
          int(fft_Nx, C_INTPTR_T), int(fft_Ny, C_INTPTR_T), int(fft_Nz, C_INTPTR_T), &
          cdata, cdata, MPI_COMM_WORLD, FFTW_FORWARD, &
          ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_OUT))
     plan_bwd = fftw_mpi_plan_dft_3d( &
          int(fft_Nx, C_INTPTR_T), int(fft_Ny, C_INTPTR_T), int(fft_Nz, C_INTPTR_T), &
          cdata, cdata, MPI_COMM_WORLD, FFTW_BACKWARD, &
          ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_IN))

     fftw_block = (fft_Nx + ncpu - 1) / ncpu

     ! Geometry-only k^2 factor for the TRANSPOSED (Ny-distributed) layout:
     !   idx_c = ky_local*(Nx*Nz) + kx*Nz + kz,  ky = ny_start+ky_local
     ! phase(k) = alpha * kfac(k);  kfac = 0.5*denom/dx^2.  k=0 -> identity.
     if(allocated(kfac)) deallocate(kfac)
     N_complex = int(loc_ny,i8b) * int(fft_Nx,i8b) * int(fft_Nz,i8b)
     allocate(kfac(0:N_complex-1))
     do ky_i = 0, int(loc_ny)-1
        iy = int(ny_start) + ky_i
        do kx_i = 0, fft_Nx-1
           do kz_i = 0, fft_Nz-1
              idx_c = int(ky_i,i8b)*int(fft_Nx,i8b)*int(fft_Nz,i8b) &
                    + int(kx_i,i8b)*int(fft_Nz,i8b) + int(kz_i,i8b)
              denom = 2.0d0*( cos(twopi*dble(kx_i)/dble(fft_Nx)) &
                            + cos(twopi*dble(iy)/dble(fft_Ny)) &
                            + cos(twopi*dble(kz_i)/dble(fft_Nz)) - 3.0d0 )
              kfac(idx_c) = 0.5d0 * denom / dx2_fft
           end do
        end do
     end do
  end if

  ! ================================================================
  ! Sparse P2P partner discovery (cached; recomputed on rebalance)
  ! ================================================================
  my_xmin = 1.0d0; my_xmax = 0.0d0
  ngrid_loc = active(ilevel)%ngrid
  do igrid = 1, ngrid_loc
     igrid_amr = active(ilevel)%igrid(igrid)
     my_xmin = min(my_xmin, xg(igrid_amr,1) - dx_fft)
     my_xmax = max(my_xmax, xg(igrid_amr,1) + dx_fft)
  end do
  if(ngrid_loc == 0) then
     my_xmin = 0.5d0; my_xmax = 0.5d0
  end if

  local_changed = (.not. partners_done .or. cached_xmin /= my_xmin .or. &
       cached_xmax /= my_xmax .or. cached_block /= fftw_block)
  call MPI_ALLREDUCE(local_changed, global_changed, 1, MPI_LOGICAL, &
       MPI_LOR, MPI_COMM_WORLD, info)
  if(global_changed) then
     if(.not. allocated(cpubox_min)) allocate(cpubox_min(ncpu))
     if(.not. allocated(cpubox_max)) allocate(cpubox_max(ncpu))
     call MPI_ALLGATHER(my_xmin, 1, MPI_DOUBLE_PRECISION, &
          cpubox_min, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
     call MPI_ALLGATHER(my_xmax, 1, MPI_DOUBLE_PRECISION, &
          cpubox_max, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)

     if(allocated(send_partners)) deallocate(send_partners)
     if(allocated(recv_partners)) deallocate(recv_partners)

     ! Send partners: FFTW slab owners whose x-range overlaps my domain
     n_send = 0
     do irank = 0, ncpu-1
        slab_x_lo = irank * fftw_block
        slab_x_hi = min((irank+1)*fftw_block, fft_Nx) - 1
        recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
        recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)
        overlap = (cpubox_max(myid) > recv_pos_lo .and. cpubox_min(myid) < recv_pos_hi)
        if(recv_pos_lo < 0.0d0) overlap = overlap .or. (cpubox_max(myid) > recv_pos_lo+1.0d0)
        if(recv_pos_hi > 1.0d0) overlap = overlap .or. (cpubox_min(myid) < recv_pos_hi-1.0d0)
        if(overlap) n_send = n_send + 1
     end do
     allocate(send_partners(n_send))
     n_send = 0
     do irank = 0, ncpu-1
        slab_x_lo = irank * fftw_block
        slab_x_hi = min((irank+1)*fftw_block, fft_Nx) - 1
        recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
        recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)
        overlap = (cpubox_max(myid) > recv_pos_lo .and. cpubox_min(myid) < recv_pos_hi)
        if(recv_pos_lo < 0.0d0) overlap = overlap .or. (cpubox_max(myid) > recv_pos_lo+1.0d0)
        if(recv_pos_hi > 1.0d0) overlap = overlap .or. (cpubox_min(myid) < recv_pos_hi-1.0d0)
        if(overlap) then
           n_send = n_send + 1
           send_partners(n_send) = irank
        end if
     end do

     ! Recv partners: RAMSES ranks whose domain overlaps my FFTW slab
     slab_x_lo = (myid-1) * fftw_block
     slab_x_hi = min(myid*fftw_block, fft_Nx) - 1
     recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
     recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)
     n_recv = 0
     do irank = 1, ncpu
        overlap = (cpubox_max(irank) > recv_pos_lo .and. cpubox_min(irank) < recv_pos_hi)
        if(recv_pos_lo < 0.0d0) overlap = overlap .or. (cpubox_max(irank) > recv_pos_lo+1.0d0)
        if(recv_pos_hi > 1.0d0) overlap = overlap .or. (cpubox_min(irank) < recv_pos_hi-1.0d0)
        if(overlap) n_recv = n_recv + 1
     end do
     allocate(recv_partners(n_recv))
     n_recv = 0
     do irank = 1, ncpu
        overlap = (cpubox_max(irank) > recv_pos_lo .and. cpubox_min(irank) < recv_pos_hi)
        if(recv_pos_lo < 0.0d0) overlap = overlap .or. (cpubox_max(irank) > recv_pos_lo+1.0d0)
        if(recv_pos_hi > 1.0d0) overlap = overlap .or. (cpubox_min(irank) < recv_pos_hi-1.0d0)
        if(overlap) then
           n_recv = n_recv + 1
           recv_partners(n_recv) = irank - 1
        end if
     end do

     cached_xmin = my_xmin; cached_xmax = my_xmax
     cached_block = fftw_block; partners_done = .true.
  end if

  ! ================================================================
  ! Count cells per destination slab
  ! ================================================================
  allocate(sendcounts(0:ncpu-1), recvcounts(0:ncpu-1))
  allocate(sdispls(0:ncpu-1), rdispls(0:ncpu-1))
  sendcounts = 0; recvcounts = 0
  do igrid = 1, ngrid_loc
     igrid_amr = active(ilevel)%igrid(igrid)
     Kx = nint(xg(igrid_amr,1) * dble(fft_Nx))
     do ind = 1, twotondim
        ix = modulo(Kx - 1 + mod(ind-1,2), fft_Nx)
        dest_rank = min(ix / fftw_block, ncpu-1)
        sendcounts(dest_rank) = sendcounts(dest_rank) + 1
     end do
  end do

  ! Exchange counts (sparse P2P)
  nreq = 0
  allocate(reqs(n_send + n_recv))
  do ip = 1, n_recv
     nreq = nreq + 1
     call MPI_IRECV(recvcounts(recv_partners(ip)), 1, MPI_INTEGER, &
          recv_partners(ip), 711, MPI_COMM_WORLD, reqs(nreq), info)
  end do
  do ip = 1, n_send
     nreq = nreq + 1
     call MPI_ISEND(sendcounts(send_partners(ip)), 1, MPI_INTEGER, &
          send_partners(ip), 711, MPI_COMM_WORLD, reqs(nreq), info)
  end do
  if(nreq > 0) then
     allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
     call MPI_WAITALL(nreq, reqs, mpi_stat, info)
     deallocate(mpi_stat)
  end if
  deallocate(reqs)

  sdispls(0) = 0; rdispls(0) = 0
  do irank = 1, ncpu-1
     sdispls(irank) = sdispls(irank-1) + sendcounts(irank-1)
     rdispls(irank) = rdispls(irank-1) + recvcounts(irank-1)
  end do
  n_send_total = sdispls(ncpu-1) + sendcounts(ncpu-1)
  n_recv_total = rdispls(ncpu-1) + recvcounts(ncpu-1)

  ! ================================================================
  ! (Re)build the static cell->slab mapping cache when the grid/partition
  ! changed. Stores, in cell-instance order, the absolute sendbuf slot, the
  ! slab idx_3d, and the AMR cell -- everything that does NOT depend on psi.
  ! ================================================================
  if(.not. cache_done .or. global_changed .or. ngrid_loc /= cached_ngrid) then
     if(allocated(pack_slot)) deallocate(pack_slot, pack_icell, pack_idx3d)
     M_cache = ngrid_loc * twotondim
     allocate(pack_slot(M_cache), pack_icell(M_cache), pack_idx3d(M_cache))
     allocate(send_idx(0:ncpu-1))
     send_idx = sdispls
     m = 0
     do igrid = 1, ngrid_loc
        igrid_amr = active(ilevel)%igrid(igrid)
        Kx = nint(xg(igrid_amr,1) * dble(fft_Nx))
        Ky = nint(xg(igrid_amr,2) * dble(fft_Ny))
        Kz = nint(xg(igrid_amr,3) * dble(fft_Nz))
        do ind = 1, twotondim
           ix = modulo(Kx - 1 + mod(ind-1,2),     fft_Nx)
           iy = modulo(Ky - 1 + mod((ind-1)/2,2), fft_Ny)
           iz = modulo(Kz - 1 + (ind-1)/4,        fft_Nz)
           dest_rank = min(ix / fftw_block, ncpu-1)
           x_local = ix - dest_rank * fftw_block
           idx_3d = x_local * fft_Ny * fft_Nz + iy * fft_Nz + iz
           iskip = ncoarse + (ind-1)*ngridmax
           icell_amr = iskip + igrid_amr
           m = m + 1
           pack_slot(m)  = send_idx(dest_rank)
           pack_idx3d(m) = idx_3d
           pack_icell(m) = icell_amr
           send_idx(dest_rank) = send_idx(dest_rank) + 1
        end do
     end do
     deallocate(send_idx)
     cached_ngrid = ngrid_loc; cache_done = .true.
  end if

  ! ================================================================
  ! Pack (idx_3d, psi_re, psi_im) triples via the cached mapping.
  ! Each m writes a distinct sendbuf slot -> threadable, byte-identical.
  ! ================================================================
  call timer('fdm-pack','start')
  allocate(sendbuf(0:3*n_send_total-1))
  allocate(recvbuf(0:3*n_recv_total-1))
  !$omp parallel do private(slot) schedule(static)
  do m = 1, M_cache
     slot = pack_slot(m)
     sendbuf(3*slot)     = dble(pack_idx3d(m))
     sendbuf(3*slot + 1) = psi_re(pack_icell(m))
     sendbuf(3*slot + 2) = psi_im(pack_icell(m))
  end do
  !$omp end parallel do

  ! Counts/displs are kept in CELL units; forward ships 3 doubles/cell
  ! (idx_3d, re, im), reverse ships 2 doubles/cell (re, im).
  call timer('fdm-p2p','start')
  nreq = 0
  allocate(reqs(n_send + n_recv))
  do ip = 1, n_recv
     irank = recv_partners(ip)
     if(recvcounts(irank) > 0) then
        nreq = nreq + 1
        call MPI_IRECV(recvbuf(3*rdispls(irank)), 3*recvcounts(irank), &
             MPI_DOUBLE_PRECISION, irank, 712, MPI_COMM_WORLD, reqs(nreq), info)
     end if
  end do
  do ip = 1, n_send
     irank = send_partners(ip)
     if(sendcounts(irank) > 0) then
        nreq = nreq + 1
        call MPI_ISEND(sendbuf(3*sdispls(irank)), 3*sendcounts(irank), &
             MPI_DOUBLE_PRECISION, irank, 712, MPI_COMM_WORLD, reqs(nreq), info)
     end if
  end do
  if(nreq > 0) then
     allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
     call MPI_WAITALL(nreq, reqs, mpi_stat, info)
     deallocate(mpi_stat)
  end if
  deallocate(reqs)

  ! Unpack into the single complex slab buffer (psi = re + i*im).
  ! Zeroing is a pure fill -> safe to thread. The scatter-add below stays
  ! SERIAL: it accumulates (+=) at ghost-overlap cells, so threading would
  ! both race and reorder the float sum, breaking byte-identity vs serial.
  !$omp parallel do schedule(static)
  do idx_c = 1_i8b, int(alloc_local,i8b)
     cdata(idx_c) = (0.0d0, 0.0d0)
  end do
  !$omp end parallel do
  do ix = 0, n_recv_total-1
     idx_3d = nint(recvbuf(3*ix))
     cdata(idx_3d + 1) = cdata(idx_3d + 1) &
          + cmplx(recvbuf(3*ix + 1), recvbuf(3*ix + 2), C_DOUBLE_COMPLEX)
  end do

  ! ================================================================
  ! Forward C2C FFT, phase rotation over full spectrum, inverse C2C
  ! ================================================================
  call timer('fdm-fft-fwd','start')
  call fftw_mpi_execute_dft(plan_fwd, cdata, cdata)

  call timer('fdm-kphase','start')
  N_complex = int(loc_ny,i8b) * int(fft_Nx,i8b) * int(fft_Nz,i8b)
  !$omp parallel do private(idx_c,phase,cos_p,sin_p,cr,ci) schedule(static)
  do idx_c = 0, N_complex-1
     phase = alpha * kfac(idx_c)
     cos_p = cos(phase); sin_p = sin(phase)
     cr = dble(cdata(idx_c+1)); ci = aimag(cdata(idx_c+1))
     cdata(idx_c+1) = cmplx(cr*cos_p - ci*sin_p, cr*sin_p + ci*cos_p, C_DOUBLE_COMPLEX)
  end do
  !$omp end parallel do

  call timer('fdm-fft-bwd','start')
  call fftw_mpi_execute_dft(plan_bwd, cdata, cdata)
  cdata = cdata / dble(N_total)

  ! ================================================================
  ! Reverse scatter: slab -> original AMR owners (values only, 2 doubles).
  ! idx_3d is NOT re-sent: data returns in the same per-rank order it was
  ! sent, and the unpack walks send_idx identically -- so the receiver maps
  ! each value back to its cell by ordering alone. Cuts reverse comm by 1/3.
  ! ================================================================
  call timer('fdm-p2p-rev','start')
  allocate(sendbuf2(0:2*n_send_total-1))
  allocate(recvbuf2(0:2*n_recv_total-1))
  ! Gather: each ix writes distinct recvbuf2 slots (no accumulation) -> safe
  ! to thread and byte-identical to serial.
  !$omp parallel do private(idx_3d) schedule(static)
  do ix = 0, n_recv_total-1
     idx_3d = nint(recvbuf(3*ix))
     recvbuf2(2*ix)     = dble(cdata(idx_3d + 1))
     recvbuf2(2*ix + 1) = aimag(cdata(idx_3d + 1))
  end do
  !$omp end parallel do
  nreq = 0
  allocate(reqs(n_send + n_recv))
  do ip = 1, n_send
     irank = send_partners(ip)
     if(sendcounts(irank) > 0) then
        nreq = nreq + 1
        call MPI_IRECV(sendbuf2(2*sdispls(irank)), 2*sendcounts(irank), &
             MPI_DOUBLE_PRECISION, irank, 713, MPI_COMM_WORLD, reqs(nreq), info)
     end if
  end do
  do ip = 1, n_recv
     irank = recv_partners(ip)
     if(recvcounts(irank) > 0) then
        nreq = nreq + 1
        call MPI_ISEND(recvbuf2(2*rdispls(irank)), 2*recvcounts(irank), &
             MPI_DOUBLE_PRECISION, irank, 713, MPI_COMM_WORLD, reqs(nreq), info)
     end if
  end do
  if(nreq > 0) then
     allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
     call MPI_WAITALL(nreq, reqs, mpi_stat, info)
     deallocate(mpi_stat)
  end if
  deallocate(reqs)

  ! Unpack back into AMR cells via the cached mapping (same cell order as the
  ! pack). Each m writes a distinct icell -> threadable, byte-identical.
  call timer('fdm-scatter','start')
  !$omp parallel do private(slot) schedule(static)
  do m = 1, M_cache
     slot = pack_slot(m)
     psi_re(pack_icell(m)) = sendbuf2(2*slot)
     psi_im(pack_icell(m)) = sendbuf2(2*slot + 1)
  end do
  !$omp end parallel do

  deallocate(sendcounts, sdispls, recvcounts, rdispls)
  deallocate(sendbuf, recvbuf, sendbuf2, recvbuf2)

end subroutine fdm_drift_fft_distributed
#endif
!################################################################
!################################################################
! Finite-difference kinetic operator for AMR fine levels
! Crank-Nicolson implicit (unitary, unconditionally stable)
!
! Solves: (1 + i*alpha*L) psi^{n+1} = (1 - i*alpha*L) psi^n
! where alpha = hbar*dt/(4*a^2*dx^2), L = discrete Laplacian
!
! Implemented as fixed-point iteration (Jacobi):
!   psi^{k+1} = [(1-i*alpha*L)*psi^n - i*alpha*L_offdiag*psi^k] / (1+6*i*alpha)
! Converges in ~5-10 iterations for typical alpha.
! Falls back to explicit if alpha < 0.1 (CFL safe, explicit is faster).
!################################################################
subroutine fdm_drift_fd(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  ! Fine-level kinetic drift selector (fdm_kinetic key)
  !   0 = explicit FD with automatic sub-cycling (default, dx^2-stable)
  !   1 = Crank-Nicolson implicit (unconditionally stable; allows large dt)
  if(fdm_kinetic == 1) then
     call timer('fdm-drift-fd','start')
     call fdm_drift_fd_cn(ilevel, dt_half)
  else
     call timer('fdm-drift-fd','start')
     call fdm_drift_fd_explicit(ilevel, dt_half)
  end if

end subroutine fdm_drift_fd
!################################################################
! Crank-Nicolson (Cayley) implicit FD kinetic step, solved matrix-free
! by a complex BiCGSTAB Krylov iteration.
!
! Approximates exp(i*alpha*L) by (1 - i*g*L)^{-1}(1 + i*g*L), g=alpha/2,
! unitary (norm-preserving) and 2nd-order in dt for ANY alpha.
! L = 7-point Laplacian (sum_nbr - 6*center); alpha = hbar*dt/(2 a^2 dx^2).
!
! We solve  A x = b  for the leaf-cell unknowns (x = correction to psi^n):
!   A     = (1 + 6 i g) I - i g N        (N = sum of 6 same-level neighbours)
!   b     = (1 - 6 i g) psi^n + i g N[psi^n]
! Coarse-fine / physical boundaries are held at psi^n (Dirichlet); inter-rank
! leaf neighbours are coupled EXACTLY through make_virtual_fine_dp inside every
! matvec, so BiCGSTAB converges to the true distributed CN solution.
!
! Why BiCGSTAB instead of literal line-ADI: a tridiagonal line solve needs
! structured 1-D lines, which the octree + MPI domain decomposition does not
! provide (lines fragment at refinement and rank boundaries, and a parallel
! cyclic-tridiagonal solver would itself need outer iteration). BiCGSTAB
! reuses the existing matrix-free neighbour stencil + ghost exchange, converges
! in ~sqrt(cond(A)) iterations even for large g (where the old Jacobi sweep,
! spectral radius 6g/sqrt(1+36g^2) -> 1, stalled), and is fully parallel-exact.
!################################################################
subroutine fdm_drift_fd_cn(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  integer::igrid,ind,iskip,icell,iter,nx_loc,ntot
  integer,parameter::max_iter=300
  real(dp)::dx,scale,dx_loc,alpha,g,cntol,bnorm,rnorm
  real(dp)::omr,omi,betr,beti,alr,ali,tr1,ti1
  complex(dp)::crho,crho_old,alp,om,bet,zrv,zts,ztt
  ! x = correction; r,rhat,p,v,s,t = BiCGSTAB work vectors (real/imag split
  ! so make_virtual_fine_dp can ghost-sync them like a grid array).
  real(dp),dimension(:),allocatable::br,bi,rr,ri,rhr,rhi,pr,pi
  real(dp),dimension(:),allocatable::vr,vi,sr,si,trr,tii,xr,xi
  ! Conservative boundary bookkeeping (hybrid)
  real(dp),dimension(:),allocatable::dmb
  integer::idim,inbor,icn,icpn
  real(dp)::cfac,vratio,psb_re,psb_im,drho_w,rho_old,sfac,gre,gim
  logical::gfound

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  alpha = hbar_code * dt_half / (2.0d0 * dx_loc**2)
  g = 0.5d0 * alpha
  cntol = 1.0d-10
  ntot = ncoarse + twotondim*ngridmax

  allocate(br(ntot),bi(ntot),rr(ntot),ri(ntot),rhr(ntot),rhi(ntot))
  allocate(pr(ntot),pi(ntot),vr(ntot),vi(ntot),sr(ntot),si(ntot))
  allocate(trr(ntot),tii(ntot),xr(ntot),xi(ntot))
  br=0d0;bi=0d0;rr=0d0;ri=0d0;rhr=0d0;rhi=0d0
  pr=0d0;pi=0d0;vr=0d0;vi=0d0;sr=0d0;si=0d0
  trr=0d0;tii=0d0;xr=0d0;xi=0d0

  ! b = (1-6ig)psi + ig N[psi] = A_{-g}[psi]   (boundary psi included;
  ! hybrid: coarse-fine ghosts enter b with weight 2, see cn_matvec)
  call cn_matvec(ilevel, -g, psi_re, psi_im, br, bi, 1)
  ! r0 = b - A_{g}[psi]  (x0 = 0 correction). Store A psi in r then subtract.
  call cn_matvec(ilevel,  g, psi_re, psi_im, rr, ri, 0)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rr(icell) = br(icell) - rr(icell)
           ri(icell) = bi(icell) - ri(icell)
           rhr(icell) = rr(icell)
           rhi(icell) = ri(icell)
        end if
        igrid = next(igrid)
     end do
  end do

  call cn_cdot(ilevel, br,bi, br,bi, zts)   ! <b,b> (real part = |b|^2)
  bnorm = sqrt(max(0.0d0, dble(zts)))
  if(bnorm == 0.0d0) then
     deallocate(br,bi,rr,ri,rhr,rhi,pr,pi,vr,vi,sr,si,trr,tii,xr,xi)
     return
  end if

  crho_old=(1.0d0,0.0d0); alp=(1.0d0,0.0d0); om=(1.0d0,0.0d0)

  do iter=1,max_iter
     call cn_cdot(ilevel, rhr,rhi, rr,ri, crho)
     if(abs(crho) == 0.0d0) exit                       ! breakdown
     bet = (crho/crho_old)*(alp/om)
     betr=dble(bet); beti=aimag(bet); omr=dble(om); omi=aimag(om)
     ! p = r + bet*(p - om*v)
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              tr1 = pr(icell) - (omr*vr(icell) - omi*vi(icell))
              ti1 = pi(icell) - (omr*vi(icell) + omi*vr(icell))
              pr(icell) = rr(icell) + (betr*tr1 - beti*ti1)
              pi(icell) = ri(icell) + (betr*ti1 + beti*tr1)
           end if
           igrid = next(igrid)
        end do
     end do
     call cn_matvec(ilevel, g, pr,pi, vr,vi, 0)        ! v = A p
     call cn_cdot(ilevel, rhr,rhi, vr,vi, zrv)
     alp = crho/zrv
     alr=dble(alp); ali=aimag(alp)
     ! s = r - alp*v
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              sr(icell) = rr(icell) - (alr*vr(icell) - ali*vi(icell))
              si(icell) = ri(icell) - (alr*vi(icell) + ali*vr(icell))
           end if
           igrid = next(igrid)
        end do
     end do
     call cn_cdot(ilevel, sr,si, sr,si, zts)
     rnorm = sqrt(max(0.0d0, dble(zts)))
     if(rnorm < cntol*bnorm) then
        ! x = x + alp*p ; converged on the half-step
        do ind=1,twotondim
           iskip = ncoarse + (ind-1)*ngridmax
           igrid = headl(myid, ilevel)
           do while(igrid > 0)
              icell = igrid + iskip
              if(son(icell) == 0) then
                 xr(icell) = xr(icell) + (alr*pr(icell) - ali*pi(icell))
                 xi(icell) = xi(icell) + (alr*pi(icell) + ali*pr(icell))
              end if
              igrid = next(igrid)
           end do
        end do
        exit
     end if
     call cn_matvec(ilevel, g, sr,si, trr,tii, 0)       ! t = A s
     call cn_cdot(ilevel, trr,tii, sr,si, zts)          ! <t,s>
     call cn_cdot(ilevel, trr,tii, trr,tii, ztt)        ! <t,t>
     if(abs(ztt) == 0.0d0) exit
     om = zts/ztt
     omr=dble(om); omi=aimag(om)
     ! x = x + alp*p + om*s ;  r = s - om*t
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              xr(icell) = xr(icell) + (alr*pr(icell) - ali*pi(icell)) &
                                    + (omr*sr(icell) - omi*si(icell))
              xi(icell) = xi(icell) + (alr*pi(icell) + ali*pr(icell)) &
                                    + (omr*si(icell) + omi*sr(icell))
              rr(icell) = sr(icell) - (omr*trr(icell) - omi*tii(icell))
              ri(icell) = si(icell) - (omr*tii(icell) + omi*trr(icell))
           end if
           igrid = next(igrid)
        end do
     end do
     call cn_cdot(ilevel, rr,ri, rr,ri, zts)
     rnorm = sqrt(max(0.0d0, dble(zts)))
     if(rnorm < cntol*bnorm) exit
     if(abs(om) == 0.0d0) exit
     crho_old = crho
  end do

  ! Commit: psi^{n+1} = psi^n + x   (leaf cells)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           psi_re(icell) = psi_re(icell) + xr(icell)
           psi_im(icell) = psi_im(icell) + xi(icell)
        end if
        igrid = next(igrid)
     end do
  end do

  ! ================================================================
  ! Conservative boundary bookkeeping (hybrid): CN with Dirichlet
  ! ghosts exchanges mass through coarse-fine faces. The change of
  ! |psi_i|^2 due to a ghost face is EXACTLY (identity, no approx)
  !   drho_i = -(hbar*dt/dx^2) * Im(psibar_i^* psi_g)
  ! with psibar = (psi^n + psi^{n+1})/2 = psi^{n+1} - x/2.
  ! Credit -drho_i * (V_fine/V_parent) to the parent-level neighbour
  ! cell so total fluid+wave mass is machine-conserved.
  ! ================================================================
  if(fdm_use_hjm) then
     allocate(dmb(ntot)); dmb = 0.0d0
     cfac = hbar_code * dt_half / dx_loc**2
     vratio = 0.5d0**ndim
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              do idim=1,ndim
                 do inbor=1,2
                    call fdm_neighbor_cell(igrid, ilevel, ind, idim, inbor, icn)
                    if(icn == 0) then
                       call fdm_wave_ghost(igrid, ilevel, idim, inbor, gre, gim, gfound)
                       if(gfound) then
                          psb_re = psi_re(icell) - 0.5d0*xr(icell)
                          psb_im = psi_im(icell) - 0.5d0*xi(icell)
                          drho_w = -cfac * (psb_re*gim - psb_im*gre)
                          icpn = nbor(igrid, 2*(idim-1)+inbor)
                          if(icpn > 0) dmb(icpn) = dmb(icpn) - drho_w*vratio
                       end if
                    end if
                 end do
              end do
           end if
           igrid = next(igrid)
        end do
     end do
     ! Reduce contributions accumulated on virtual parent cells to owners
     call make_virtual_reverse_dp(dmb(1), ilevel-1)
     ! Apply to local leaf cells at the parent level
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel-1)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0 .and. dmb(icell) /= 0.0d0) then
              if(ilevel-1 < fdm_first_wave_level) then
                 ! Fluid parent: psi_re IS rho
                 psi_re(icell) = max(psi_re(icell) + dmb(icell), 1.0d-15)
              else
                 ! Wave parent: rescale amplitude, keep phase
                 rho_old = psi_re(icell)**2 + psi_im(icell)**2
                 if(rho_old > 0.0d0) then
                    sfac = sqrt(max(0.0d0, (rho_old + dmb(icell))/rho_old))
                    psi_re(icell) = psi_re(icell)*sfac
                    psi_im(icell) = psi_im(icell)*sfac
                 end if
              end if
           end if
           igrid = next(igrid)
        end do
     end do
     call make_virtual_fine_dp(psi_re(1), ilevel-1)
     call make_virtual_fine_dp(psi_im(1), ilevel-1)
     deallocate(dmb)
  end if

  deallocate(br,bi,rr,ri,rhr,rhi,pr,pi,vr,vi,sr,si,trr,tii,xr,xi)

end subroutine fdm_drift_fd_cn
!################################################################
! Matrix-free apply of the CN operator on leaf cells at ilevel:
!   out = (1 + 6 i gg) in - i gg N[in]
! N[in] = sum of 6 same-level neighbours (Neumann/zero-gradient fallback:
! an absent neighbour contributes the centre value, matching the explicit
! and old Jacobi stencils). Ghosts of `in` are synced first so inter-rank
! leaf neighbours couple exactly. `in` must be 0 on non-leaf/boundary cells
! (correction vectors satisfy this; for the b/Apsi calls `in`=psi is the
! intended Dirichlet boundary data).
!################################################################
subroutine cn_matvec(ilevel, gg, inr, ini, outr, outi, bmode)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::gg
  real(dp),dimension(*)::inr,ini
  real(dp),dimension(*)::outr,outi
  ! bmode: closure at coarse-fine faces (missing same-level neighbour).
  !   Hybrid (fdm_use_hjm):
  !     0 = correction vector: Dirichlet, contributes 0
  !     1 = b-build from psi^n: contributes 2*ghost (psi_g^n + psi_g^{n+1})
  !         with the ghost interpolated/converted from the fluid parent
  !   Legacy (no hybrid): Neumann (centre value), byte-identical to before
  integer,intent(in)::bmode
  integer::igrid,ind,iskip,icell,idim,inbor,icn
  real(dp)::nr,ni,gre,gim
  logical::gfound

  call make_virtual_fine_dp(inr(1), ilevel)
  call make_virtual_fine_dp(ini(1), ilevel)

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           nr = 0.0d0; ni = 0.0d0
           do idim=1,ndim
              do inbor=1,2
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, inbor, icn)
                 if(icn > 0) then
                    nr = nr + inr(icn); ni = ni + ini(icn)
                 else if(fdm_use_hjm) then
                    if(bmode == 1) then
                       call fdm_wave_ghost(igrid, ilevel, idim, inbor, gre, gim, gfound)
                       if(gfound) then
                          nr = nr + 2.0d0*gre; ni = ni + 2.0d0*gim
                       else
                          nr = nr + inr(icell); ni = ni + ini(icell)
                       end if
                    end if
                    ! bmode==0: Dirichlet — missing link contributes nothing
                 else
                    nr = nr + inr(icell); ni = ni + ini(icell)
                 end if
              end do
           end do
           ! (1+6i gg)(inr+i ini) - i gg (nr + i ni)
           outr(icell) = inr(icell) - 6.0d0*gg*ini(icell) + gg*ni
           outi(icell) = ini(icell) + 6.0d0*gg*inr(icell) - gg*nr
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine cn_matvec
!################################################################
! Complex inner product <a,b> = sum_leaf conj(a)*b over ilevel leaf cells,
! reduced across ranks. Returns a complex(dp) scalar.
!################################################################
subroutine cn_cdot(ilevel, ar, ai, br, bi, res)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),dimension(*)::ar,ai,br,bi
  complex(dp),intent(out)::res
  integer::igrid,ind,iskip,icell,info
  real(dp)::loc(2),glob(2)

  loc(1)=0.0d0; loc(2)=0.0d0
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           loc(1) = loc(1) + ar(icell)*br(icell) + ai(icell)*bi(icell)
           loc(2) = loc(2) + ar(icell)*bi(icell) - ai(icell)*br(icell)
        end if
        igrid = next(igrid)
     end do
  end do
  glob = loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(loc, glob, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
  res = cmplx(glob(1), glob(2), kind=dp)

end subroutine cn_cdot
!################################################################
! Copy scratch (new_re,new_im) back into psi_re,psi_im for leaf cells.
!################################################################
subroutine cn_copy_to_psi(ilevel, src_re, src_im)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),dimension(:),intent(in)::src_re,src_im
  integer::igrid,ind,iskip,icell
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           psi_re(icell) = src_re(icell)
           psi_im(icell) = src_im(icell)
        end if
        igrid = next(igrid)
     end do
  end do
end subroutine cn_copy_to_psi
!################################################################
! Maximum squared de Broglie wavenumber |grad theta|^2 over leaf cells at
! ilevel, used by the Madelung/hybrid timestep (fdm_hybrid). theta is the
! wavefunction phase; grad theta is the Madelung velocity field (up to hbar/a):
!   d(theta)/dx_i = Im(conj(psi) d psi/dx_i)/|psi|^2
!                 = (psi_re*d psi_im - psi_im*d psi_re)/(|psi|^2)
! Computed branch-cut-free from central differences of psi (no atan2). A
! density floor guards near-empty cells. Result is a global MPI max so the
! timestep is set by the highest-k (most quantum) cell anywhere on the level.
!################################################################
subroutine fdm_kdb_max_level(ilevel, k2max)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(out)::k2max

  integer::igrid,ind,iskip,icell,idim,icp,icm,nx_loc,info
  real(dp)::dx,scale,dx_loc,inv2dx,rho2,rho_floor
  real(dp)::dre,dim,dth,k2,k2loc,k2glob

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  inv2dx = 0.5d0 / dx_loc
  rho_floor = 1.0d-6            ! mean |psi|^2 = 1; ignore near-empty cells

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  k2loc = 0.0d0
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rho2 = psi_re(icell)**2 + psi_im(icell)**2
           if(rho2 > rho_floor) then
              k2 = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icp)
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icm)
                 if(icp > 0 .and. icm > 0) then
                    dre = (psi_re(icp) - psi_re(icm)) * inv2dx
                    dim = (psi_im(icp) - psi_im(icm)) * inv2dx
                    ! d(theta)/dx_i
                    dth = (psi_re(icell)*dim - psi_im(icell)*dre) / rho2
                    k2 = k2 + dth*dth
                 end if
              end do
              if(k2 > k2loc) k2loc = k2
           end if
        end if
        igrid = next(igrid)
     end do
  end do

  k2glob = k2loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(k2loc, k2glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#endif
  k2max = k2glob

end subroutine fdm_kdb_max_level
!################################################################
! Compute max Madelung phase gradient |nabla theta| at a level.
! Phase gradient = Im(psi* nabla psi)/|psi|^2 — isolates velocity
! from amplitude gradient noise.
! Cells with C1 > fdm_hjm_C1 (about to be refined to wave solver)
! are excluded — their interference-driven |nabla theta| is
! unphysical for the fluid CFL.
!################################################################
subroutine fdm_vmax_level(ilevel, dtheta_max)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(out)::dtheta_max

  integer::igrid,ind,iskip,icell,idim,icp,icm,nx_loc,info
  real(dp)::dx,scale,dx_loc,inv2dx,rho2,sqrho_c,sqrho_L,sqrho_R
  real(dp)::dre,dim,dth,dth2,dth2_loc,dth2_glob
  real(dp)::c1_dim,c1_cell
  real(dp)::S_L,S_R,re_L,re_R,im_L,im_R
  logical::is_fluid,ok_p,ok_m
  integer::n_total,n_leaf,n_lowrho,n_c1skip,n_valid
  integer::n_total_g,n_leaf_g,n_lowrho_g,n_c1skip_g,n_valid_g
  integer,parameter::NBINS_C1=8
  real(dp),parameter::c1_edges(NBINS_C1+1) = &
       (/0.0d0, 0.01d0, 0.03d0, 0.1d0, 0.3d0, 1.0d0, 3.0d0, 10.0d0, 1.0d30/)
  integer::c1_hist(NBINS_C1), c1_hist_g(NBINS_C1), ibin
  real(dp)::c1_max_loc, c1_max_glob
  real(dp)::dth_max_loc, dth_at_c1max

  n_total=0; n_leaf=0; n_lowrho=0; n_c1skip=0; n_valid=0
  c1_hist = 0; c1_max_loc = 0.0d0; dth_max_loc = 0.0d0
  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  inv2dx = 0.5d0 / dx_loc
  is_fluid = (fdm_use_hjm .and. ilevel < fdm_first_wave_level)

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  dth2_loc = 0.0d0
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        n_total = n_total + 1
        if(son(icell) == 0) then
           n_leaf = n_leaf + 1
           if(is_fluid) then
              rho2 = max(psi_re(icell),0.0d0)   ! psi_re = rho on fluid levels
           else
              rho2 = psi_re(icell)**2 + psi_im(icell)**2
           end if
           sqrho_c = sqrt(rho2)
           if(sqrho_c > 1.0d-8) then
              c1_cell = 0.0d0
              dth2 = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icp)
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icm)
                 ! A refined neighbour is not a co-level cell. Its stored value
                 ! is whatever the restriction left there, and on a fluid level
                 ! that S was rebased onto the parent's own earlier value, which
                 ! keeps it continuous in time but not with the neighbouring
                 ! leaf's HJM action. Differencing across that seam manufactures
                 ! a huge grad(S), which collapses dt_adv through dtheta_max.
                 ! Mirror the centre instead, exactly as the HJM solver does.
                 ok_p = .false.; ok_m = .false.
                 if(icp > 0) ok_p = (son(icp) == 0)
                 if(icm > 0) ok_m = (son(icm) == 0)
                 if(ok_p .or. ok_m) then
                    if(is_fluid) then
                       sqrho_R = sqrho_c; sqrho_L = sqrho_c
                       S_R = psi_im(icell); S_L = psi_im(icell)
                       if(ok_p)then
                          sqrho_R = sqrt(max(psi_re(icp),0.0d0)); S_R = psi_im(icp)
                       end if
                       if(ok_m)then
                          sqrho_L = sqrt(max(psi_re(icm),0.0d0)); S_L = psi_im(icm)
                       end if
                       c1_dim = abs(sqrho_R - 2.0d0*sqrho_c + sqrho_L) / sqrho_c
                       c1_cell = max(c1_cell, c1_dim)
                       ! dtheta = grad(S)/hbar (S unwrapped)
                       dth = (S_R - S_L) * inv2dx / hbar_code
                       dth2 = dth2 + dth*dth
                    else
                       re_R = psi_re(icell); im_R = psi_im(icell)
                       re_L = psi_re(icell); im_L = psi_im(icell)
                       if(ok_p)then
                          re_R = psi_re(icp); im_R = psi_im(icp)
                       end if
                       if(ok_m)then
                          re_L = psi_re(icm); im_L = psi_im(icm)
                       end if
                       sqrho_R = sqrt(re_R**2 + im_R**2)
                       sqrho_L = sqrt(re_L**2 + im_L**2)
                       c1_dim = abs(sqrho_R - 2.0d0*sqrho_c + sqrho_L) / sqrho_c
                       c1_cell = max(c1_cell, c1_dim)
                       dre = (re_R - re_L) * inv2dx
                       dim = (im_R - im_L) * inv2dx
                       dth = (psi_re(icell)*dim - psi_im(icell)*dre) / rho2
                       dth2 = dth2 + dth*dth
                    end if
                 end if
              end do
              do ibin=1,NBINS_C1
                 if(c1_cell < c1_edges(ibin+1)) then
                    c1_hist(ibin) = c1_hist(ibin) + 1
                    exit
                 end if
              end do
              if(c1_cell > c1_max_loc) then
                 c1_max_loc = c1_cell
                 dth_at_c1max = sqrt(dth2)
              end if
              if(c1_cell < fdm_hjm_C1) then
                 n_valid = n_valid + 1
                 if(dth2 > dth2_loc) then
                    dth2_loc = dth2
                    dth_max_loc = c1_cell
                 end if
              else
                 n_c1skip = n_c1skip + 1
              end if
           else
              n_lowrho = n_lowrho + 1
           end if
        end if
        igrid = next(igrid)
     end do
  end do

  n_total_g=n_total; n_leaf_g=n_leaf; n_lowrho_g=n_lowrho
  n_c1skip_g=n_c1skip; n_valid_g=n_valid
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(n_total,  n_total_g,  1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(n_leaf,   n_leaf_g,   1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(n_lowrho, n_lowrho_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(n_c1skip, n_c1skip_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(n_valid,  n_valid_g,  1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#endif
  c1_hist_g = c1_hist
  c1_max_glob = c1_max_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(c1_hist, c1_hist_g, NBINS_C1, MPI_INTEGER, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(c1_max_loc, c1_max_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)
#endif
  if(myid==1 .and. nstep_coarse_old < 60) then
     write(*,'(" VMAX_DBG: total=",I12," leaf=",I12," lowrho=",I12,&
          &" c1skip=",I12," valid=",I12)') &
          n_total_g, n_leaf_g, n_lowrho_g, n_c1skip_g, n_valid_g
     write(*,'(" C1_HIST: <0.01=",I10," <0.03=",I10," <0.1=",I10,&
          &" <0.3=",I10)') c1_hist_g(1),c1_hist_g(2),c1_hist_g(3),c1_hist_g(4)
     write(*,'("          <1.0=",I10," <3.0=",I10," <10=",I10,&
          &" >=10=",I10)') c1_hist_g(5),c1_hist_g(6),c1_hist_g(7),c1_hist_g(8)
     write(*,'(" C1_MAX=",1PE10.3," dth_at_c1max=",1PE10.3,&
          &" dth_max_c1=",1PE10.3)') c1_max_glob, dth_at_c1max, dth_max_loc
  end if

  dth2_glob = dth2_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dth2_loc, dth2_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#endif
  dtheta_max = sqrt(dth2_glob)

end subroutine fdm_vmax_level
!################################################################
! Compute max Madelung C1 = |sqrt(rho)_R - 2 sqrt(rho)_C + sqrt(rho)_L| / sqrt(rho)_C
! across all leaf cells at a level (for C1-adaptive quantum CFL).
!################################################################
subroutine fdm_c1max_level(ilevel, c1max)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(out)::c1max

  integer::igrid,ind,iskip,icell,idim,icp,icm,info
  real(dp)::sqrho_c,sqrho_L,sqrho_R,c1_dim,c1_cell
  real(dp)::c1_loc,c1_glob

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  c1_loc = 0.0d0
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           sqrho_c = sqrt(psi_re(icell)**2 + psi_im(icell)**2)
           if(sqrho_c > 1.0d-8) then
              c1_cell = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icm)
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icp)
                 if(icm > 0 .and. icp > 0) then
                    sqrho_L = sqrt(psi_re(icm)**2 + psi_im(icm)**2)
                    sqrho_R = sqrt(psi_re(icp)**2 + psi_im(icp)**2)
                    c1_dim = abs(sqrho_R - 2.0d0*sqrho_c + sqrho_L) / sqrho_c
                    c1_cell = max(c1_cell, c1_dim)
                 end if
              end do
              if(c1_cell > c1_loc) c1_loc = c1_cell
           end if
        end if
        igrid = next(igrid)
     end do
  end do

  c1_glob = c1_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(c1_loc, c1_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
#endif
  c1max = c1_glob

end subroutine fdm_c1max_level
!################################################################
! Explicit FD kinetic step (with automatic sub-cycling for stability)
!################################################################
subroutine fdm_drift_fd_explicit(ilevel, dt_half)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_half

  integer::igrid,ind,iskip,icell,isub,nsub
  real(dp)::dx,scale,dx_loc,alpha,alpha_sub,dt_sub
  real(dp)::lap_re,lap_im,psi_re_c,psi_im_c
  real(dp)::gre,gim
  logical::gfound
  integer::nx_loc,idim,inbor,icell_nbor

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Full alpha for the requested dt_half
  alpha = hbar_code * dt_half / (2.0d0 * dx_loc**2)

  ! Sub-cycle if alpha > 1/6 (CFL limit)
  nsub = max(1, int(6.0d0*alpha) + 1)
  dt_sub = dt_half / dble(nsub)
  alpha_sub = alpha / dble(nsub)

  do isub=1,nsub
     ! Ghost zone exchange at start of each sub-step
     call make_virtual_fine_dp(psi_re(1), ilevel)
     call make_virtual_fine_dp(psi_im(1), ilevel)

     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              psi_re_c = psi_re(icell)
              psi_im_c = psi_im(icell)

              ! 7-point Laplacian
              lap_re = -6.0d0 * psi_re_c
              lap_im = -6.0d0 * psi_im_c
              do idim=1,ndim
                 do inbor=1,2
                    call fdm_neighbor_cell(igrid, ilevel, ind, idim, inbor, icell_nbor)
                    if(icell_nbor > 0) then
                       lap_re = lap_re + psi_re(icell_nbor)
                       lap_im = lap_im + psi_im(icell_nbor)
                    else if(fdm_use_hjm) then
                       ! Dirichlet ghost from the (fluid) parent level.
                       ! NOTE: no conservative bookkeeping here (only the CN
                       ! path has it) — use fdm_kinetic=1 for hybrid runs.
                       call fdm_wave_ghost(igrid, ilevel, idim, inbor, gre, gim, gfound)
                       if(gfound) then
                          lap_re = lap_re + gre
                          lap_im = lap_im + gim
                       else
                          lap_re = lap_re + psi_re_c
                          lap_im = lap_im + psi_im_c
                       end if
                    else
                       lap_re = lap_re + psi_re_c
                       lap_im = lap_im + psi_im_c
                    end if
                 end do
              end do

              ! Explicit update: dpsi = i*(hbar/2a^2)*L*dt
              psi_re(icell) = psi_re_c - alpha_sub * lap_im
              psi_im(icell) = psi_im_c + alpha_sub * lap_re
           end if
           igrid = next(igrid)
        end do
     end do
  end do

end subroutine fdm_drift_fd_explicit
!################################################################
! Find the neighbor cell of subcell ind in grid igrid
! in direction (idim, inbor): idim=1,2,3; inbor=1(left),2(right)
!################################################################
subroutine fdm_neighbor_cell(igrid, ilevel, ind, idim, inbor, icell_nbor)
  use amr_commons
  use morton_hash
  implicit none
  integer,intent(in)::igrid, ilevel, ind, idim, inbor
  integer,intent(out)::icell_nbor

  integer::ind_bit, ind_nbor, igrid_nbor, iskip

  ind_bit = iand(ishft(ind-1, -(idim-1)), 1)  ! 0 or 1

  if(inbor == 2) then
     ! Right neighbor
     if(ind_bit == 0) then
        ! Neighbor is in same grid, flip this dimension bit
        ind_nbor = ind + ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid + iskip
     else
        ! Neighbor is in adjacent grid (Morton lookup, defrag-safe)
        igrid_nbor = morton_nbor_grid(igrid, ilevel, 2*idim)
        if(igrid_nbor == 0) then
           icell_nbor = 0; return
        end if
        ind_nbor = ind - ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid_nbor + iskip
     end if
  else
     ! Left neighbor
     if(ind_bit == 1) then
        ! Neighbor is in same grid
        ind_nbor = ind - ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid + iskip
     else
        ! Neighbor is in adjacent grid (Morton lookup, defrag-safe)
        igrid_nbor = morton_nbor_grid(igrid, ilevel, 2*idim-1)
        if(igrid_nbor == 0) then
           icell_nbor = 0; return
        end if
        ind_nbor = ind + ishft(1, idim-1)
        iskip = ncoarse + (ind_nbor-1)*ngridmax
        icell_nbor = igrid_nbor + iskip
     end if
  end if

end subroutine fdm_neighbor_cell
!################################################################
! Dirichlet ghost value for a wave-level cell whose same-level
! neighbour is missing (coarse-fine boundary, hybrid runs).
! Interpolates the father-level field to the ghost fine-cell
! position along the face normal (w=0.25 toward our own father)
! and converts fluid (rho,S) -> psi when the parent is a fluid
! level. The fluid HJ equation evolves the SAME action S as the
! Schrodinger phase (up to QP), so the converted phase stays
! coherent with the wave patch interior.
! found=.false. -> caller falls back to the legacy closure.
!################################################################
subroutine fdm_wave_ghost(igrid, ilevel, idim, inbor, gre, gim, found)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::igrid, ilevel, idim, inbor
  real(dp),intent(out)::gre, gim
  logical,intent(out)::found

  integer::icell_pf, icell_pn
  real(dp)::rho_g, S_g, amp, theta
  real(dp),parameter::w = 0.25d0

  found = .false.
  gre = 0.0d0; gim = 0.0d0

  icell_pf = father(igrid)
  icell_pn = nbor(igrid, 2*(idim-1)+inbor)
  if(icell_pf <= 0 .or. icell_pn <= 0) return

  if(fdm_use_hjm .and. ilevel-1 < fdm_first_wave_level) then
     ! Parent holds (rho, S): interpolate, then convert to psi
     rho_g = (1.0d0-w)*psi_re(icell_pn) + w*psi_re(icell_pf)
     S_g   = (1.0d0-w)*psi_im(icell_pn) + w*psi_im(icell_pf)
     amp = sqrt(max(rho_g, 0.0d0))
     theta = S_g / hbar_code
     gre = amp*cos(theta)
     gim = amp*sin(theta)
  else
     ! Parent holds psi: interpolate (Re, Im) directly
     gre = (1.0d0-w)*psi_re(icell_pn) + w*psi_re(icell_pf)
     gim = (1.0d0-w)*psi_im(icell_pn) + w*psi_im(icell_pf)
  end if
  found = .true.

end subroutine fdm_wave_ghost
!################################################################
!################################################################
! Potential operator (Kick): exp(-i Phi dt / hbar)
! Local pointwise rotation in real space
!################################################################
subroutine fdm_kick(ilevel, dt_loc)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_loc

  integer::igrid,ind,iskip,icell
  real(dp)::phase,cos_p,sin_p,re_old,im_old

  if(hbar_code <= 0.0d0) return

  ! Loop over all leaf cells at this level
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           ! Schrödinger: i*hbar*dpsi/dt = m*Phi*psi
           ! => dpsi/dt = -i*(Phi/hbar)*psi
           ! => psi(t+dt) = psi(t) * exp(-i*Phi*dt/hbar)
           ! exp(-i*theta) = cos(theta) - i*sin(theta)
           ! where theta = Phi*dt/hbar
           phase = phi(icell) * dt_loc / hbar_code

           cos_p = cos(phase)
           sin_p = sin(phase)

           re_old = psi_re(icell)
           im_old = psi_im(icell)

           ! psi_new = psi * exp(-i*phase) = (re+i*im)*(cos-i*sin)
           psi_re(icell) =  re_old*cos_p + im_old*sin_p
           psi_im(icell) = -re_old*sin_p + im_old*cos_p
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_kick
!################################################################
!################################################################
! Compute hbar_code from physical constants and code units
! Called once at initialization
!################################################################
subroutine fdm_compute_hbar()
  use amr_commons
  use fdm_commons
  implicit none

  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp)::m_axion_g

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  ! m_axion in grams: m [eV] * eV_to_erg / c^2
  m_axion_g = m_axion * eV_to_erg / c_light**2

  ! hbar_code = hbar_cgs / (m_axion * scale_l * scale_v)
  ! This is the effective Planck constant in code units
  hbar_code = hbar_cgs / (m_axion_g * scale_l * scale_v)

  if(myid==1) then
     write(*,'(A,ES10.3,A)') ' FDM: m_axion = ', m_axion, ' eV'
     write(*,'(A,ES10.3)')   '   m_axion [g] = ', m_axion_g
     write(*,'(A,ES10.3)')   '   hbar_code   = ', hbar_code
     write(*,'(A,ES10.3,A)') '   lambda_dB(100km/s) = ', &
          2.0d0*acos(-1.0d0)*hbar_code/(100.0d5/scale_v) * scale_l/3.0857d21, ' kpc'
  end if

  ! Initialize HJM hybrid solver if enabled
  call fdm_hjm_init()

end subroutine fdm_compute_hbar
!################################################################
!################################################################
! Initialize psi from the grafic initial conditions (Madelung transform)
! For fresh start (nrestart=0):
!   amplitude  |psi| = sqrt(rho),  rho = 1 + dfact*delta   (from ic_deltac)
!   phase      theta = (aexp^2/hbar_code) * phi_v,
!              where phi_v is the velocity potential, nabla phi_v = v_code,
!              solved by FFT on the periodic base grid from ic_velc{x,y,z}.
!   psi_re = |psi|*cos(theta),  psi_im = |psi|*sin(theta)
! The phase solve requires FFTW (-DUSE_FFTW); without it the wavefunction is
! initialized with zero phase (amplitude only) and a warning is emitted.
!################################################################
subroutine fdm_init_psi()
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,ind,iskip,icell,info,i1,i2,i3
  integer::ix,iy,iz,nx_loc,ng1,ng2,ng3
  integer(i8b)::N_total,idx
  real(dp)::dx,dxL,scale,dx_loc,rho_cell,theta_cell,amp
  real(dp)::dm_share
  real(dp)::xx1,xx2,xx3,p1,p2,p3,mass_loc,mass_glob
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),allocatable::rho_3d(:),theta_3d(:)
  real(kind=4),allocatable::init_plane(:,:)
  character(LEN=256)::filename
  logical::ok
  integer,parameter::tag=2207

  if(.not.use_fdm) return

#if defined(USE_FFTW) && !defined(WITHOUTMPI)
  ! Multi-rank: fully distributed init (FFTW-MPI slabs, no full-grid arrays).
  if(ncpu > 1) then
     call fdm_init_psi_distributed()
     return
  end if
#endif

  ! Base (levelmin) grafic grid dimensions
  ng1 = n1(levelmin); ng2 = n2(levelmin); ng3 = n3(levelmin)
  N_total = int(ng1,i8b)*int(ng2,i8b)*int(ng3,i8b)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen/dble(nx_loc)
  dxL = 0.5d0**levelmin
  dx_loc = dxL*scale
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  ! DM share of total matter: with hydro the gas carries Omega_b, so
  ! |psi|^2 must average (1 - Omega_b/Omega_m) for the total Poisson
  ! source to have mean 1 (mirrors the particle-mass convention).
  dm_share = 1.0d0
  if(hydro .and. omega_m > 0.0d0) dm_share = 1.0d0 - omega_b/omega_m
  if(myid==1) write(*,'(A,F8.5)') ' FDM: DM share of total matter = ', dm_share

  allocate(rho_3d(1:N_total))
  allocate(theta_3d(1:N_total))
  rho_3d = 1.0d0      ! mean density (overwritten if a density IC is found)
  theta_3d = 0.0d0    ! zero-phase default

  !-------------------------------------------------------------
  ! Amplitude: read density overdensity, rho = 1 + dfact*delta
  ! (rank 1 reads each plane, broadcasts to all ranks)
  !-------------------------------------------------------------
  allocate(init_plane(1:ng1,1:ng2))
  filename = TRIM(initfile(levelmin))//'/ic_deltac'
  INQUIRE(file=filename,exist=ok)
  if(.not.ok)then
     filename = TRIM(initfile(levelmin))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok)
  end if
  if(ok)then
     if(myid==1)write(*,*)' FDM: reading density IC '//TRIM(filename)
     if(myid==1)then
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10) ! skip header
     end if
     do i3=1,ng3
        if(myid==1)then
           read(10) ((init_plane(i1,i2),i1=1,ng1),i2=1,ng2)
        else
           init_plane = 0.0
        end if
#ifndef WITHOUTMPI
        call MPI_BCAST(init_plane,ng1*ng2,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
        do i2=1,ng2
           do i1=1,ng1
              idx = int(i1,i8b)+int(i2-1,i8b)*int(ng1,i8b) &
                   +int(i3-1,i8b)*int(ng1,i8b)*int(ng2,i8b)
              rho_3d(idx) = 1.0d0 + dfact(levelmin)*dble(init_plane(i1,i2))
           end do
        end do
     end do
     if(myid==1)close(10)
  else
     if(myid==1)write(*,*)' FDM WARNING: no density IC (ic_deltac/ic_deltab) found, using rho=1'
  end if
  deallocate(init_plane)

  !-------------------------------------------------------------
  ! Phase: solve velocity potential by FFT (rank 1), broadcast theta
  !-------------------------------------------------------------
#ifdef USE_FFTW
  if(hbar_code <= 0.0d0)then
     if(myid==1)write(*,*)' FDM WARNING: hbar_code <= 0, initializing psi with zero phase'
  else
     if(myid==1) call fdm_init_phase_fft(ng1,ng2,ng3,N_total,dx_loc,theta_3d)
#ifndef WITHOUTMPI
     call MPI_BCAST(theta_3d,int(N_total),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
#endif
  end if
#else
  if(myid==1)write(*,*)' FDM WARNING: USE_FFTW not enabled — psi initialized with zero phase (no velocity potential)'
#endif

  !-------------------------------------------------------------
  ! Scatter rho,theta onto AMR leaf cells -> psi = sqrt(rho) e^{i theta}
  !-------------------------------------------------------------
  mass_loc = 0.0d0
  do ilevel=levelmin,nlevelmax
     dx = 0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5d0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5d0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0)then
              ! Map cell centre onto the levelmin grafic grid
              p1 = xg(igrid,1)+xc(ind,1)-skip_loc(1)
              p2 = xg(igrid,2)+xc(ind,2)-skip_loc(2)
              p3 = xg(igrid,3)+xc(ind,3)-skip_loc(3)
              xx1=(p1*(dxini(levelmin)/dxL)-xoff1(levelmin))/dxini(levelmin)
              xx2=(p2*(dxini(levelmin)/dxL)-xoff2(levelmin))/dxini(levelmin)
              xx3=(p3*(dxini(levelmin)/dxL)-xoff3(levelmin))/dxini(levelmin)
              i1=max(1,min(ng1,int(xx1)+1))
              i2=max(1,min(ng2,int(xx2)+1))
              i3=max(1,min(ng3,int(xx3)+1))
              idx = int(i1,i8b)+int(i2-1,i8b)*int(ng1,i8b) &
                   +int(i3-1,i8b)*int(ng1,i8b)*int(ng2,i8b)
              rho_cell = dm_share*max(rho_3d(idx),0.0d0)
              theta_cell = theta_3d(idx)
              if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
                 ! Fluid level: store (rho, S) with S unwrapped
                 psi_re(icell) = rho_cell
                 psi_im(icell) = hbar_code*theta_cell
              else
                 amp = sqrt(rho_cell)
                 psi_re(icell) = amp*cos(theta_cell)
                 psi_im(icell) = amp*sin(theta_cell)
              end if
              mass_loc = mass_loc + rho_cell * dx**3
           end if
           igrid = next(igrid)
        end do
     end do
  end do

  deallocate(rho_3d,theta_3d)

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc, mass_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
#else
  mass_glob = mass_loc
#endif

  ! Ghost zones for psi at every level
  do ilevel=levelmin,nlevelmax
     call make_virtual_fine_dp(psi_re(1), ilevel)
     call make_virtual_fine_dp(psi_im(1), ilevel)
  end do

  if(myid==1) write(*,'(A,ES12.5)') ' FDM: initial |psi|^2 total mass = ', mass_glob

end subroutine fdm_init_psi
#if defined(USE_FFTW) && !defined(WITHOUTMPI)
!################################################################
! Fully distributed Madelung init for psi (multi-rank).
!
! Memory model: NO full N_total array is ever held on any rank.
!   * The levelmin grafic grid is decomposed into FFTW-MPI slabs along
!     the i3 (plane) axis -- rank r owns planes [start0+1 .. start0+loc_n0].
!   * Density (rho) and phase (theta) are computed only on the local slab.
!   * psi = sqrt(rho)*exp(i*theta) is formed on the slab, then delivered to
!     each rank's AMR leaf cells by a request/reply MPI_Alltoallv gather.
!
! This mirrors the slab decomposition of fdm_drift_fft_distributed but uses
! a one-shot Alltoallv gather (init runs once) instead of cached sparse P2P.
!################################################################
subroutine fdm_init_psi_distributed()
  use amr_commons
  use poisson_commons
  use fdm_commons
  use iso_c_binding
  implicit none
  include 'mpif.h'
  include 'fftw3-mpi.f03'

  integer::ng1,ng2,ng3,nx_loc
  integer::ilevel,ind,iskip,icell,igrid,info
  integer::i1,i2,i3,ix,iy,iz,zl
  integer::owner,irank,p
  integer(i8b)::N_total,gidx,locsize,lidx
  integer(C_INTPTR_T)::loc_n0,start0,alloc_local
  integer::loc_n0_i,start0_i
  real(dp)::dx,dxL,scale,dx_loc,amp,rho_c,re,im,dvol
  real(dp)::dm_share
  real(dp)::p1,p2,p3,xx1,xx2,xx3
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::mass_loc,mass_glob

  real(dp),allocatable::theta_slab(:)      ! 0:locsize-1
  real(dp),allocatable::rho_slab(:)
  real(dp),allocatable::psislab_re(:),psislab_im(:)
  real(kind=4),allocatable::init_plane(:,:)
  character(LEN=256)::filename
  logical::ok

  integer,allocatable::start0_all(:),locn0_all(:)
  integer,allocatable::scount(:),rcount(:),sdispls(:),rdispls(:),pos(:)
  integer,allocatable::scount2(:),rcount2(:),sdispls2(:),rdispls2(:)
  integer(i8b),allocatable::send_gidx(:),recv_gidx(:)
  integer,allocatable::cellmap(:)
  integer,allocatable::leaflvl(:)
  real(dp),allocatable::celldx(:)
  real(dp),allocatable::reply(:),recvreply(:)
  integer::n_leaf,n_recv,slot,s

  ! ---- Grid / unit setup (identical to serial fdm_init_psi) ----
  ! DM share: gas carries Omega_b under hydro (see serial routine)
  dm_share = 1.0d0
  if(hydro .and. omega_m > 0.0d0) dm_share = 1.0d0 - omega_b/omega_m
  if(myid==1) write(*,'(A,F8.5)') ' FDM: DM share of total matter = ', dm_share
  ng1 = n1(levelmin); ng2 = n2(levelmin); ng3 = n3(levelmin)
  N_total = int(ng1,i8b)*int(ng2,i8b)*int(ng3,i8b)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen/dble(nx_loc)
  dxL = 0.5d0**levelmin
  dx_loc = dxL*scale
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  ! ---- FFTW-MPI slab decomposition along i3 (plane) axis ----
  ! Logical FFT array has dims (ng3,ng2,ng1): first axis = i3 (distributed),
  ! last axis = i1 (contiguous, matches grafic plane storage i1-fastest).
  call fftw_mpi_init()
  alloc_local = fftw_mpi_local_size_3d( &
       int(ng3,C_INTPTR_T), int(ng2,C_INTPTR_T), int(ng1,C_INTPTR_T), &
       MPI_COMM_WORLD, loc_n0, start0)
  loc_n0_i = int(loc_n0); start0_i = int(start0)
  locsize = int(loc_n0,i8b)*int(ng2,i8b)*int(ng1,i8b)

  allocate(theta_slab(0:max(locsize,1_i8b)-1))
  allocate(rho_slab(0:max(locsize,1_i8b)-1))
  theta_slab = 0.0d0
  rho_slab   = 1.0d0

  ! ---- Phase: distributed velocity-potential FFT -> theta on slab ----
  if(hbar_code <= 0.0d0)then
     if(myid==1)write(*,*)' FDM WARNING: hbar_code <= 0, initializing psi with zero phase'
  else
     call fdm_init_phase_fft_dist(ng1,ng2,ng3,loc_n0,start0,alloc_local, &
          dx_loc,theta_slab,locsize)
  end if

  ! ---- Amplitude: read density IC plane-by-plane, keep only local slab ----
  allocate(init_plane(1:ng1,1:ng2))
  filename = TRIM(initfile(levelmin))//'/ic_deltac'
  INQUIRE(file=filename,exist=ok)
  if(.not.ok)then
     filename = TRIM(initfile(levelmin))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok)
  end if
  if(ok)then
     if(myid==1)write(*,*)' FDM: reading density IC '//TRIM(filename)
     if(myid==1)then
        open(10,file=filename,form='unformatted')
        rewind 10
        read(10) ! skip header
     end if
     do i3=1,ng3
        if(myid==1)then
           read(10) ((init_plane(i1,i2),i1=1,ng1),i2=1,ng2)
        else
           init_plane = 0.0
        end if
        call MPI_BCAST(init_plane,ng1*ng2,MPI_REAL,0,MPI_COMM_WORLD,info)
        if(i3-1>=start0_i .and. i3-1<start0_i+loc_n0_i)then
           zl = i3-1-start0_i
           do i2=1,ng2
              do i1=1,ng1
                 lidx = int(zl,i8b)*int(ng2,i8b)*int(ng1,i8b) &
                      + int(i2-1,i8b)*int(ng1,i8b) + int(i1-1,i8b)
                 rho_slab(lidx) = 1.0d0 + dfact(levelmin)*dble(init_plane(i1,i2))
              end do
           end do
        end if
     end do
     if(myid==1)close(10)
  else
     if(myid==1)write(*,*)' FDM WARNING: no density IC (ic_deltac/ic_deltab) found, using rho=1'
  end if
  deallocate(init_plane)

  ! ---- Form the slab field ----
  ! Hybrid: ship (rho, S=hbar*theta) so fluid levels can store the
  ! UNWRAPPED action directly (atan2 of psi would wrap S and break
  ! HJM at large boxes where |dS| per cell > pi*hbar).
  ! Pure wave: ship psi = sqrt(rho)*exp(i theta) as before.
  allocate(psislab_re(0:max(locsize,1_i8b)-1))
  allocate(psislab_im(0:max(locsize,1_i8b)-1))
  if(fdm_use_hjm)then
     do lidx=0,locsize-1
        psislab_re(lidx) = dm_share*max(rho_slab(lidx),0.0d0)
        psislab_im(lidx) = hbar_code*theta_slab(lidx)
     end do
  else
     do lidx=0,locsize-1
        rho_c = dm_share*max(rho_slab(lidx),0.0d0)
        amp = sqrt(rho_c)
        psislab_re(lidx) = amp*cos(theta_slab(lidx))
        psislab_im(lidx) = amp*sin(theta_slab(lidx))
     end do
  end if
  deallocate(theta_slab,rho_slab)

  ! ---- Slab ownership table (which rank owns each i3 plane) ----
  allocate(start0_all(ncpu),locn0_all(ncpu))
  call MPI_ALLGATHER(start0_i,1,MPI_INTEGER,start0_all,1,MPI_INTEGER,MPI_COMM_WORLD,info)
  call MPI_ALLGATHER(loc_n0_i,1,MPI_INTEGER,locn0_all,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! ================================================================
  ! Gather slab psi -> local AMR leaf cells (request/reply Alltoallv)
  ! ================================================================
  allocate(scount(0:ncpu-1),rcount(0:ncpu-1))
  allocate(sdispls(0:ncpu-1),rdispls(0:ncpu-1),pos(0:ncpu-1))
  scount = 0

  ! Pass 1: count requests per owning rank
  n_leaf = 0
  do ilevel=levelmin,nlevelmax
     dx = 0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4; iy=(ind-1-4*iz)/2; ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5d0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5d0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0)then
              p1 = xg(igrid,1)+xc(ind,1)-skip_loc(1)
              p2 = xg(igrid,2)+xc(ind,2)-skip_loc(2)
              p3 = xg(igrid,3)+xc(ind,3)-skip_loc(3)
              xx1=(p1*(dxini(levelmin)/dxL)-xoff1(levelmin))/dxini(levelmin)
              xx2=(p2*(dxini(levelmin)/dxL)-xoff2(levelmin))/dxini(levelmin)
              xx3=(p3*(dxini(levelmin)/dxL)-xoff3(levelmin))/dxini(levelmin)
              i3=max(1,min(ng3,int(xx3)+1))
              owner = slab_owner(i3,start0_all,locn0_all,ncpu)
              scount(owner) = scount(owner) + 1
              n_leaf = n_leaf + 1
           end if
           igrid = next(igrid)
        end do
     end do
  end do

  sdispls(0) = 0
  do irank=1,ncpu-1
     sdispls(irank) = sdispls(irank-1) + scount(irank-1)
  end do

  allocate(send_gidx(max(n_leaf,1)))
  allocate(cellmap(max(n_leaf,1)))
  allocate(celldx(max(n_leaf,1)))
  allocate(leaflvl(max(n_leaf,1)))
  pos = sdispls

  ! Pass 2: fill request buffers (bucketed by owner), record cell back-map
  do ilevel=levelmin,nlevelmax
     dx = 0.5d0**ilevel
     dvol = dx**3
     do ind=1,twotondim
        iz=(ind-1)/4; iy=(ind-1-4*iz)/2; ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5d0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5d0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0)then
              p1 = xg(igrid,1)+xc(ind,1)-skip_loc(1)
              p2 = xg(igrid,2)+xc(ind,2)-skip_loc(2)
              p3 = xg(igrid,3)+xc(ind,3)-skip_loc(3)
              xx1=(p1*(dxini(levelmin)/dxL)-xoff1(levelmin))/dxini(levelmin)
              xx2=(p2*(dxini(levelmin)/dxL)-xoff2(levelmin))/dxini(levelmin)
              xx3=(p3*(dxini(levelmin)/dxL)-xoff3(levelmin))/dxini(levelmin)
              i1=max(1,min(ng1,int(xx1)+1))
              i2=max(1,min(ng2,int(xx2)+1))
              i3=max(1,min(ng3,int(xx3)+1))
              gidx = int(i1-1,i8b)+int(i2-1,i8b)*int(ng1,i8b) &
                   +int(i3-1,i8b)*int(ng1,i8b)*int(ng2,i8b)
              owner = slab_owner(i3,start0_all,locn0_all,ncpu)
              slot = pos(owner) + 1            ! 1-based slot
              pos(owner) = pos(owner) + 1
              send_gidx(slot) = gidx
              cellmap(slot)   = icell
              celldx(slot)    = dvol
              leaflvl(slot)   = ilevel
           end if
           igrid = next(igrid)
        end do
     end do
  end do

  ! Exchange request counts and build recv layout
  call MPI_ALLTOALL(scount,1,MPI_INTEGER,rcount,1,MPI_INTEGER,MPI_COMM_WORLD,info)
  rdispls(0) = 0
  do irank=1,ncpu-1
     rdispls(irank) = rdispls(irank-1) + rcount(irank-1)
  end do
  n_recv = rdispls(ncpu-1) + rcount(ncpu-1)

  allocate(recv_gidx(max(n_recv,1)))
  call MPI_ALLTOALLV(send_gidx,scount,sdispls,MPI_INTEGER8, &
       recv_gidx,rcount,rdispls,MPI_INTEGER8,MPI_COMM_WORLD,info)

  ! Owner side: look up psi for each requested grafic index
  allocate(reply(0:max(2*n_recv,1)-1))
  do p=1,n_recv
     gidx = recv_gidx(p)
     i3 = int(gidx/(int(ng1,i8b)*int(ng2,i8b)))
     i2 = int( mod(gidx,int(ng1,i8b)*int(ng2,i8b)) / int(ng1,i8b) )
     i1 = int( mod(gidx,int(ng1,i8b)) )
     zl = i3 - start0_i
     lidx = int(zl,i8b)*int(ng2,i8b)*int(ng1,i8b) &
          + int(i2,i8b)*int(ng1,i8b) + int(i1,i8b)
     reply(2*(p-1))   = psislab_re(lidx)
     reply(2*(p-1)+1) = psislab_im(lidx)
  end do

  ! Reply path: counts/displs are doubled (2 reals per request)
  allocate(scount2(0:ncpu-1),rcount2(0:ncpu-1))
  allocate(sdispls2(0:ncpu-1),rdispls2(0:ncpu-1))
  do irank=0,ncpu-1
     scount2(irank)  = 2*scount(irank)
     rcount2(irank)  = 2*rcount(irank)
     sdispls2(irank) = 2*sdispls(irank)
     rdispls2(irank) = 2*rdispls(irank)
  end do
  allocate(recvreply(0:max(2*n_leaf,1)-1))
  ! owner -> requester: send rcount2 (from rdispls2), recv scount2 (into sdispls2)
  call MPI_ALLTOALLV(reply,rcount2,rdispls2,MPI_DOUBLE_PRECISION, &
       recvreply,scount2,sdispls2,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)

  ! Scatter replies back into AMR cells (same order as send_gidx)
  mass_loc = 0.0d0
  if(fdm_use_hjm)then
     ! reply carries (rho, S): fluid levels store as-is (S unwrapped);
     ! wave levels convert to psi.
     do s=1,n_leaf
        icell = cellmap(s)
        rho_c = recvreply(2*(s-1))
        if(leaflvl(s) < fdm_first_wave_level)then
           psi_re(icell) = rho_c
           psi_im(icell) = recvreply(2*(s-1)+1)
        else
           amp = sqrt(max(rho_c,0.0d0))
           psi_re(icell) = amp*cos(recvreply(2*(s-1)+1)/hbar_code)
           psi_im(icell) = amp*sin(recvreply(2*(s-1)+1)/hbar_code)
        end if
        mass_loc = mass_loc + rho_c*celldx(s)
     end do
  else
     do s=1,n_leaf
        icell = cellmap(s)
        re = recvreply(2*(s-1))
        im = recvreply(2*(s-1)+1)
        psi_re(icell) = re
        psi_im(icell) = im
        mass_loc = mass_loc + (re*re+im*im)*celldx(s)
     end do
  end if

  call MPI_ALLREDUCE(mass_loc,mass_glob,1,MPI_DOUBLE_PRECISION, &
       MPI_SUM,MPI_COMM_WORLD,info)

  ! Ghost zones for psi at every level
  do ilevel=levelmin,nlevelmax
     call make_virtual_fine_dp(psi_re(1), ilevel)
     call make_virtual_fine_dp(psi_im(1), ilevel)
  end do

  if(myid==1) write(*,'(A,ES12.5)') ' FDM: initial |psi|^2 total mass = ', mass_glob

  deallocate(psislab_re,psislab_im)
  deallocate(start0_all,locn0_all)
  deallocate(scount,rcount,sdispls,rdispls,pos)
  deallocate(scount2,rcount2,sdispls2,rdispls2)
  deallocate(send_gidx,recv_gidx,cellmap,leaflvl,celldx,reply,recvreply)

contains
  integer function slab_owner(i3,s0,l0,np)
    integer,intent(in)::i3,np
    integer,intent(in)::s0(np),l0(np)
    integer::r
    slab_owner = np-1
    do r=1,np
       if(i3-1>=s0(r) .and. i3-1<s0(r)+l0(r))then
          slab_owner = r-1
          return
       end if
    end do
  end function slab_owner
end subroutine fdm_init_psi_distributed
!################################################################
! Distributed velocity-potential FFT (FFTW-MPI slabs along i3).
! Fills theta_slab(0:locsize-1) = (aexp^2/hbar_code)*phi_v on the local slab.
! Logical FFT dims (ng3,ng2,ng1): i1 contiguous (last), i3 distributed (first).
!   phi_v(k) = (i*dx_loc/denom)*[Sx*Vx + Sy*Vy + Sz*Vz]
!   Sx=sin(2pi K1/ng1) (x-vel, last axis), Sz=sin(2pi K3/ng3) (z-vel, first axis)
!   denom = 2(cosK1+cosK2+cosK3-3),  k=0 -> 0
!################################################################
subroutine fdm_init_phase_fft_dist(ng1,ng2,ng3,loc_n0,start0,alloc_local, &
     dx_loc,theta_slab,locsize)
  use amr_commons
  use fdm_commons
  use iso_c_binding
  implicit none
  include 'mpif.h'
  include 'fftw3-mpi.f03'
  integer,intent(in)::ng1,ng2,ng3
  integer(C_INTPTR_T),intent(in)::loc_n0,start0,alloc_local
  real(dp),intent(in)::dx_loc
  integer(i8b),intent(in)::locsize
  real(dp),intent(out)::theta_slab(0:locsize-1)

  integer(i8b)::N_total,lidx
  integer::i1,i2,i3,ivar,zl,info,al,b,c,start0_i,loc_n0_i
  integer::K1,K3
  real(dp)::twopi,Sx,Sy,Sz,denom,Xr,Xi,fac,vscale
  real(kind=4),allocatable::plane(:,:)
  type(C_PTR)::pvx,pvy,pvz,plan_fwd,plan_bwd
  complex(C_DOUBLE_COMPLEX),pointer::cvx(:),cvy(:),cvz(:),cv(:)
  character(LEN=256)::filename
  logical::ok

  twopi = 2.0d0*acos(-1.0d0)
  N_total = int(ng1,i8b)*int(ng2,i8b)*int(ng3,i8b)
  start0_i = int(start0); loc_n0_i = int(loc_n0)
  vscale = dfact(levelmin)*vfact(1)*dx_loc/dxini(levelmin)/vfact(levelmin)

  pvx = fftw_alloc_complex(alloc_local); call c_f_pointer(pvx,cvx,[int(alloc_local)])
  pvy = fftw_alloc_complex(alloc_local); call c_f_pointer(pvy,cvy,[int(alloc_local)])
  pvz = fftw_alloc_complex(alloc_local); call c_f_pointer(pvz,cvz,[int(alloc_local)])

  plan_fwd = fftw_mpi_plan_dft_3d( &
       int(ng3,C_INTPTR_T),int(ng2,C_INTPTR_T),int(ng1,C_INTPTR_T), &
       cvx,cvx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
  plan_bwd = fftw_mpi_plan_dft_3d( &
       int(ng3,C_INTPTR_T),int(ng2,C_INTPTR_T),int(ng1,C_INTPTR_T), &
       cvx,cvx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE)

  allocate(plane(1:ng1,1:ng2))

  ! ---- Read each velocity component into its slab buffer ----
  do ivar=1,3
     if(ivar==1)then; cv=>cvx; filename=TRIM(initfile(levelmin))//'/ic_velcx'; end if
     if(ivar==2)then; cv=>cvy; filename=TRIM(initfile(levelmin))//'/ic_velcy'; end if
     if(ivar==3)then; cv=>cvz; filename=TRIM(initfile(levelmin))//'/ic_velcz'; end if
     cv = (0.0d0,0.0d0)
     INQUIRE(file=filename,exist=ok)
     if(.not.ok)then
        if(myid==1)write(*,*)' FDM WARNING: velocity IC missing '//TRIM(filename)//' -> treated as 0'
        cycle
     end if
     if(myid==1)write(*,*)' FDM: reading velocity IC '//TRIM(filename)
     if(myid==1)then
        open(11,file=filename,form='unformatted')
        rewind 11
        read(11) ! skip header
     end if
     do i3=1,ng3
        if(myid==1)then
           read(11) ((plane(i1,i2),i1=1,ng1),i2=1,ng2)
        else
           plane = 0.0
        end if
        call MPI_BCAST(plane,ng1*ng2,MPI_REAL,0,MPI_COMM_WORLD,info)
        if(i3-1>=start0_i .and. i3-1<start0_i+loc_n0_i)then
           zl = i3-1-start0_i
           do i2=1,ng2
              do i1=1,ng1
                 lidx = int(zl,i8b)*int(ng2,i8b)*int(ng1,i8b) &
                      + int(i2-1,i8b)*int(ng1,i8b) + int(i1-1,i8b)
                 cv(lidx+1) = cmplx(vscale*dble(plane(i1,i2)),0.0d0,C_DOUBLE_COMPLEX)
              end do
           end do
        end if
     end do
     if(myid==1)close(11)
  end do

  ! ---- Forward FFTs (in-place, new-array execute) ----
  call fftw_mpi_execute_dft(plan_fwd,cvx,cvx)
  call fftw_mpi_execute_dft(plan_fwd,cvy,cvy)
  call fftw_mpi_execute_dft(plan_fwd,cvz,cvz)

  ! ---- Combine: phih = (i*dx/denom)(Sx Vx + Sy Vy + Sz Vz) into cvx ----
  ! Local spectrum (non-transposed): first axis K3 = start0+al (z-vel),
  ! middle K2 = b (y-vel), last K1 = c (x-vel).
  do al=0,loc_n0_i-1
     K3 = start0_i + al
     Sz = sin(twopi*dble(K3)/dble(ng3))
     do b=0,ng2-1
        Sy = sin(twopi*dble(b)/dble(ng2))
        do c=0,ng1-1
           K1 = c
           Sx = sin(twopi*dble(K1)/dble(ng1))
           lidx = int(al,i8b)*int(ng2,i8b)*int(ng1,i8b) &
                + int(b,i8b)*int(ng1,i8b) + int(c,i8b)
           denom = 2.0d0*( cos(twopi*dble(c)/dble(ng1)) &
                         + cos(twopi*dble(b)/dble(ng2)) &
                         + cos(twopi*dble(K3)/dble(ng3)) - 3.0d0 )
           if(abs(denom) < 1.0d-30)then
              cvx(lidx+1) = (0.0d0,0.0d0)
           else
              Xr = Sx*dble(cvx(lidx+1)) + Sy*dble(cvy(lidx+1)) + Sz*dble(cvz(lidx+1))
              Xi = Sx*aimag(cvx(lidx+1)) + Sy*aimag(cvy(lidx+1)) + Sz*aimag(cvz(lidx+1))
              fac = dx_loc/denom
              cvx(lidx+1) = cmplx(-fac*Xi, fac*Xr, C_DOUBLE_COMPLEX)
           end if
        end do
     end do
  end do

  ! ---- Inverse FFT, normalize, theta = phi_v/hbar (v_code = hbar*grad theta) ----
  call fftw_mpi_execute_dft(plan_bwd,cvx,cvx)
  do lidx=0,locsize-1
     theta_slab(lidx) = dble(cvx(lidx+1))/(hbar_code*dble(N_total))
  end do

  deallocate(plane)
  call fftw_destroy_plan(plan_fwd)
  call fftw_destroy_plan(plan_bwd)
  call fftw_free(pvx); call fftw_free(pvy); call fftw_free(pvz)

end subroutine fdm_init_phase_fft_dist
#endif
#ifdef USE_FFTW
!################################################################
! Solve the velocity potential phi_v (nabla phi_v = v_code) on the
! periodic base grid by FFT, return theta = (aexp^2/hbar_code)*phi_v.
! Serial: executed on rank 1 only over the full ng1*ng2*ng3 grid.
!   phi_v(k) = (i*dx_loc/denom) * [Sx*Vx + Sy*Vy + Sz*Vz]
!   Sd = sin(2pi Kd/Nd),  denom = 2(cosKx+cosKy+cosKz-3),  k=0 -> 0
!################################################################
subroutine fdm_init_phase_fft(ng1,ng2,ng3,N_total,dx_loc,theta_3d)
  use amr_commons
  use fdm_commons
  use iso_c_binding
  implicit none
  include 'fftw3.f03'
  integer,intent(in)::ng1,ng2,ng3
  integer(i8b),intent(in)::N_total
  real(dp),intent(in)::dx_loc
  real(dp),intent(inout)::theta_3d(1:N_total)

  integer(i8b)::N_complex,idx,idx_c
  integer::i1,i2,i3,ivar,Kx,Ky,Kz
  real(dp)::twopi,Sx,Sy,Sz,denom,Xr,Xi,fac,vscale,fmin,fmax
  real(kind=4),allocatable::plane(:,:)
  real(dp),allocatable::vx(:),vy(:),vz(:),phi_v(:)
  complex(C_DOUBLE_COMPLEX),allocatable::vxh(:),vyh(:),vzh(:),phih(:)
  type(C_PTR)::plan_x,plan_y,plan_z,plan_bwd
  character(LEN=256)::filename
  logical::ok

  twopi = 2.0d0*acos(-1.0d0)
  N_complex = int(ng1/2+1,i8b)*int(ng2,i8b)*int(ng3,i8b)

  allocate(vx(1:N_total),vy(1:N_total),vz(1:N_total),phi_v(1:N_total))
  allocate(vxh(1:N_complex),vyh(1:N_complex),vzh(1:N_complex),phih(1:N_complex))
  allocate(plane(1:ng1,1:ng2))
  vx=0.0d0; vy=0.0d0; vz=0.0d0

  ! Velocity rescale to code units (identical to init_flow_fine, level=levelmin)
  vscale = dfact(levelmin)*vfact(1)*dx_loc/dxini(levelmin)/vfact(levelmin)

  do ivar=1,3
     if(ivar==1) filename = TRIM(initfile(levelmin))//'/ic_velcx'
     if(ivar==2) filename = TRIM(initfile(levelmin))//'/ic_velcy'
     if(ivar==3) filename = TRIM(initfile(levelmin))//'/ic_velcz'
     INQUIRE(file=filename,exist=ok)
     if(.not.ok)then
        write(*,*)' FDM WARNING: velocity IC missing '//TRIM(filename)//' -> treated as 0'
        cycle
     end if
     write(*,*)' FDM: reading velocity IC '//TRIM(filename)
     open(11,file=filename,form='unformatted')
     rewind 11
     read(11) ! skip header
     do i3=1,ng3
        read(11) ((plane(i1,i2),i1=1,ng1),i2=1,ng2)
        do i2=1,ng2
           do i1=1,ng1
              idx = int(i1,i8b)+int(i2-1,i8b)*int(ng1,i8b) &
                   +int(i3-1,i8b)*int(ng1,i8b)*int(ng2,i8b)
              if(ivar==1) vx(idx) = vscale*dble(plane(i1,i2))
              if(ivar==2) vy(idx) = vscale*dble(plane(i1,i2))
              if(ivar==3) vz(idx) = vscale*dble(plane(i1,i2))
           end do
        end do
     end do
     close(11)
  end do

  ! Forward R2C transforms
  call dfftw_plan_dft_r2c_3d(plan_x, ng1, ng2, ng3, vx, vxh, FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan_x, vx, vxh)
  call dfftw_destroy_plan(plan_x)
  call dfftw_plan_dft_r2c_3d(plan_y, ng1, ng2, ng3, vy, vyh, FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan_y, vy, vyh)
  call dfftw_destroy_plan(plan_y)
  call dfftw_plan_dft_r2c_3d(plan_z, ng1, ng2, ng3, vz, vzh, FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan_z, vz, vzh)
  call dfftw_destroy_plan(plan_z)

  ! phi_v(k) = (i*dx_loc/denom) * (Sx*Vx + Sy*Vy + Sz*Vz)
  do Kz=0,ng3-1
     do Ky=0,ng2-1
        do Kx=0,ng1/2
           idx_c = int(Kx+1,i8b)+int(Ky,i8b)*int(ng1/2+1,i8b) &
                +int(Kz,i8b)*int(ng1/2+1,i8b)*int(ng2,i8b)
           denom = 2.0d0*( cos(twopi*dble(Kx)/dble(ng1)) &
                         + cos(twopi*dble(Ky)/dble(ng2)) &
                         + cos(twopi*dble(Kz)/dble(ng3)) - 3.0d0 )
           if(abs(denom) < 1.0d-30)then
              phih(idx_c) = cmplx(0.0d0,0.0d0,C_DOUBLE_COMPLEX)
           else
              Sx = sin(twopi*dble(Kx)/dble(ng1))
              Sy = sin(twopi*dble(Ky)/dble(ng2))
              Sz = sin(twopi*dble(Kz)/dble(ng3))
              Xr = Sx*dble(vxh(idx_c)) + Sy*dble(vyh(idx_c)) + Sz*dble(vzh(idx_c))
              Xi = Sx*aimag(vxh(idx_c)) + Sy*aimag(vyh(idx_c)) + Sz*aimag(vzh(idx_c))
              fac = dx_loc/denom
              ! multiply (Xr+iXi) by i*fac -> (-fac*Xi) + i(fac*Xr)
              phih(idx_c) = cmplx(-fac*Xi, fac*Xr, C_DOUBLE_COMPLEX)
           end if
        end do
     end do
  end do

  ! Inverse C2R, normalize, theta = phi_v/hbar_code (v_code = hbar*grad theta)
  call dfftw_plan_dft_c2r_3d(plan_bwd, ng1, ng2, ng3, phih, phi_v, FFTW_ESTIMATE)
  call dfftw_execute_dft_c2r(plan_bwd, phih, phi_v)
  call dfftw_destroy_plan(plan_bwd)

  fmin=1.0d30; fmax=-1.0d30
  do idx=1,N_total
     phi_v(idx) = phi_v(idx)/dble(N_total)
     theta_3d(idx) = phi_v(idx)/hbar_code
     if(theta_3d(idx)<fmin) fmin=theta_3d(idx)
     if(theta_3d(idx)>fmax) fmax=theta_3d(idx)
  end do
  write(*,'(A,2ES12.4)') ' FDM: initial phase theta range [min,max] = ', fmin, fmax

  deallocate(vx,vy,vz,phi_v,vxh,vyh,vzh,phih,plane)

end subroutine fdm_init_phase_fft
#endif
!################################################################
!################################################################
! Phase 4: AMR prolongation for psi (mass-conserving)
! Called when new refined grids are created at ilevel from ilevel-1
! Trilinear interpolation + renormalization to preserve |psi|^2 integral
!################################################################
subroutine fdm_prolong(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,icell_coarse,ind_child,iskip_child,icell_child
  integer::igrid_p,ind_p,ival,idim,icL,icR,ix,iy,iz
  real(dp)::psi_re_c,psi_im_c
  real(dp)::grad_re(3),grad_im(3)
  real(dp)::sx(3)
  real(dp)::rho_c,S_c,rho_L,rho_R,S_L,S_R
  real(dp)::grad_rho(3),grad_S(3)
  real(dp)::rho_child,S_child,amp_child,theta_child,dS
  real(dp)::pi_h,twopi_h
  real(dp)::rho_parent,sum_rho_child,renorm
  logical::use_rhoS,child_is_wave

  if(.not.use_fdm) return
  if(ilevel <= levelmin) return

  ! use_rhoS: the PARENT level (ilevel-1) is a fluid level holding (rho,S)
  use_rhoS = (fdm_use_hjm .and. ilevel-1 < fdm_first_wave_level)
  child_is_wave = (.not.fdm_use_hjm) .or. (ilevel >= fdm_first_wave_level)
  pi_h = 3.14159265358979d0 * hbar_code
  twopi_h = 2.0d0 * pi_h

  call make_virtual_fine_dp(psi_re(1), ilevel-1)
  call make_virtual_fine_dp(psi_im(1), ilevel-1)

  igrid = headl(myid, ilevel)
  do while(igrid > 0)
     icell_coarse = father(igrid)
     if(icell_coarse <= 0) then
        igrid = next(igrid)
        cycle
     end if

     psi_re_c = psi_re(icell_coarse)
     psi_im_c = psi_im(icell_coarse)

     ival = icell_coarse - ncoarse
     ind_p = (ival - 1) / ngridmax + 1
     igrid_p = ival - (ind_p - 1) * ngridmax

     if(use_rhoS) then
        ! Parent is a fluid level: psi_re=rho, psi_im=S (unwrapped) — read
        ! directly, no atan2. Plain S differences are exact.
        rho_c = psi_re_c
        S_c = psi_im_c

        do idim=1,ndim
           call fdm_neighbor_cell(igrid_p, ilevel-1, ind_p, idim, 1, icL)
           call fdm_neighbor_cell(igrid_p, ilevel-1, ind_p, idim, 2, icR)

           if(icL > 0) then
              rho_L = psi_re(icL); S_L = psi_im(icL)
           else
              rho_L = rho_c; S_L = S_c
           end if
           if(icR > 0) then
              rho_R = psi_re(icR); S_R = psi_im(icR)
           else
              rho_R = rho_c; S_R = S_c
           end if

           if(icL > 0 .and. icR > 0) then
              grad_rho(idim) = (rho_R - rho_L) * 0.5d0
              grad_S(idim) = (S_R - S_L) * 0.5d0
           else if(icR > 0) then
              grad_rho(idim) = rho_R - rho_c
              grad_S(idim) = S_R - S_c
           else if(icL > 0) then
              grad_rho(idim) = rho_c - rho_L
              grad_S(idim) = S_c - S_L
           else
              grad_rho(idim) = 0.0d0
              grad_S(idim) = 0.0d0
           end if
        end do

        do ind_child=1,twotondim
           iskip_child = ncoarse + (ind_child-1)*ngridmax
           icell_child = igrid + iskip_child
           ix = ibits(ind_child-1, 0, 1)
           iy = ibits(ind_child-1, 1, 1)
           iz = ibits(ind_child-1, 2, 1)
           sx(1) = dble(2*ix - 1)
           sx(2) = dble(2*iy - 1)
           sx(3) = dble(2*iz - 1)
           rho_child = rho_c + 0.25d0 * &
                (sx(1)*grad_rho(1) + sx(2)*grad_rho(2) + sx(3)*grad_rho(3))
           S_child = S_c + 0.25d0 * &
                (sx(1)*grad_S(1) + sx(2)*grad_S(2) + sx(3)*grad_S(3))
           if(rho_child < 1.0d-15) rho_child = 1.0d-15
           if(child_is_wave) then
              ! Fluid -> wave boundary: convert to psi
              amp_child = sqrt(rho_child)
              theta_child = S_child / hbar_code
              psi_re(icell_child) = amp_child * cos(theta_child)
              psi_im(icell_child) = amp_child * sin(theta_child)
           else
              ! Fluid -> fluid: store (rho, S) directly
              psi_re(icell_child) = rho_child
              psi_im(icell_child) = S_child
           end if
        end do
     else
        ! Standard (Re, Im) interpolation for wave levels
        do idim=1,ndim
           call fdm_neighbor_cell(igrid_p, ilevel-1, ind_p, idim, 1, icL)
           call fdm_neighbor_cell(igrid_p, ilevel-1, ind_p, idim, 2, icR)
           if(icL > 0 .and. icR > 0) then
              grad_re(idim) = (psi_re(icR) - psi_re(icL)) * 0.5d0
              grad_im(idim) = (psi_im(icR) - psi_im(icL)) * 0.5d0
           else if(icR > 0) then
              grad_re(idim) = psi_re(icR) - psi_re_c
              grad_im(idim) = psi_im(icR) - psi_im_c
           else if(icL > 0) then
              grad_re(idim) = psi_re_c - psi_re(icL)
              grad_im(idim) = psi_im_c - psi_im(icL)
           else
              grad_re(idim) = 0.0d0
              grad_im(idim) = 0.0d0
           end if
        end do

        do ind_child=1,twotondim
           iskip_child = ncoarse + (ind_child-1)*ngridmax
           icell_child = igrid + iskip_child
           ix = ibits(ind_child-1, 0, 1)
           iy = ibits(ind_child-1, 1, 1)
           iz = ibits(ind_child-1, 2, 1)
           sx(1) = dble(2*ix - 1)
           sx(2) = dble(2*iy - 1)
           sx(3) = dble(2*iz - 1)
           psi_re(icell_child) = psi_re_c + 0.25d0 * &
                (sx(1)*grad_re(1) + sx(2)*grad_re(2) + sx(3)*grad_re(3))
           psi_im(icell_child) = psi_im_c + 0.25d0 * &
                (sx(1)*grad_im(1) + sx(2)*grad_im(2) + sx(3)*grad_im(3))
        end do
        ! Renormalize to conserve mass: Σ|ψ_child|² = 8|ψ_parent|²
        rho_parent = psi_re_c**2 + psi_im_c**2
        sum_rho_child = 0.0d0
        do ind_child=1,twotondim
           iskip_child = ncoarse + (ind_child-1)*ngridmax
           icell_child = igrid + iskip_child
           sum_rho_child = sum_rho_child + &
                psi_re(icell_child)**2 + psi_im(icell_child)**2
        end do
        if(sum_rho_child > 0.0d0) then
           renorm = sqrt(8.0d0 * rho_parent / sum_rho_child)
           do ind_child=1,twotondim
              iskip_child = ncoarse + (ind_child-1)*ngridmax
              icell_child = igrid + iskip_child
              psi_re(icell_child) = psi_re(icell_child) * renorm
              psi_im(icell_child) = psi_im(icell_child) * renorm
           end do
        end if
     end if

     igrid = next(igrid)
  end do

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

end subroutine fdm_prolong
!################################################################
!################################################################
! Lightweight CFL check for HJM fluid levels.
! Computes max |nabla theta| and returns the advection CFL dt.
! Used to prevent CFL violation when newly refined levels inherit
! the parent's (too large) dt.
!################################################################
subroutine fdm_hjm_local_cfl(ilevel, dt_cfl)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  real(dp),intent(out)::dt_cfl

  integer::igrid,ind,iskip,icell,idim,icp,icm,nx_loc,info
  real(dp)::dx,scale,dx_loc,inv2dx,rho2
  real(dp)::dre,dim_val,dth,dth2,dth2_loc,dth2_glob
  logical::is_fluid

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale
  inv2dx = 0.5d0 / dx_loc
  is_fluid = (fdm_use_hjm .and. ilevel < fdm_first_wave_level)

  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  dth2_loc = 0.0d0
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           if(is_fluid) then
              rho2 = max(psi_re(icell),0.0d0)
           else
              rho2 = psi_re(icell)**2 + psi_im(icell)**2
           end if
           if(rho2 > 1.0d-16) then
              dth2 = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icp)
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icm)
                 if(icp > 0 .and. icm > 0) then
                    if(is_fluid) then
                       dth = (psi_im(icp) - psi_im(icm)) * inv2dx / hbar_code
                    else
                       dre = (psi_re(icp) - psi_re(icm)) * inv2dx
                       dim_val = (psi_im(icp) - psi_im(icm)) * inv2dx
                       dth = (psi_re(icell)*dim_val - psi_im(icell)*dre) / rho2
                    end if
                    dth2 = dth2 + dth*dth
                 end if
              end do
              if(dth2 > dth2_loc) dth2_loc = dth2
           end if
        end if
        igrid = next(igrid)
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dth2_loc, dth2_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)
#else
  dth2_glob = dth2_loc
#endif

  if(dth2_glob > 0.0d0) then
     dt_cfl = fdm_courant * dx_loc / (hbar_code * sqrt(dth2_glob))
  else
     dt_cfl = 1.0d30
  end if

end subroutine fdm_hjm_local_cfl
!################################################################
!################################################################
! Phase 4: AMR restriction for psi (mass-conserving)
! Average children to parent cell when derefinement occurs
!################################################################
subroutine fdm_restrict(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,igrid_child
  integer::ind_child,iskip_child,icell_child
  real(dp)::sum_re,sum_im,rho_avg,S_avg,S_old,twopi_h,amp2,sfac
  logical::parent_fluid,child_fluid

  if(.not.use_fdm) return
  if(ilevel < levelmin) return

  parent_fluid = (fdm_use_hjm .and. ilevel   < fdm_first_wave_level)
  child_fluid  = (fdm_use_hjm .and. ilevel+1 < fdm_first_wave_level)
  twopi_h = 2.0d0 * 3.14159265358979d0 * hbar_code

  ! For each cell at ilevel that has children (son /= 0),
  ! average the children back to parent.
  ! fluid->fluid: plain (rho,S) average.  wave->wave: plain psi average.
  ! wave child -> fluid parent: rho = <|psi|^2>, S = hbar*atan2(<psi>)
  ! rebased onto the parent's previous unwrapped S by continuity.
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        igrid_child = son(icell)
        if(igrid_child > 0) then
           sum_re = 0.0d0
           sum_im = 0.0d0
           rho_avg = 0.0d0
           do ind_child=1,twotondim
              iskip_child = ncoarse + (ind_child-1)*ngridmax
              icell_child = igrid_child + iskip_child
              sum_re = sum_re + psi_re(icell_child)
              sum_im = sum_im + psi_im(icell_child)
              if(.not.child_fluid) then
                 rho_avg = rho_avg + psi_re(icell_child)**2 + psi_im(icell_child)**2
              end if
           end do
           if(parent_fluid .and. .not.child_fluid) then
              rho_avg = rho_avg / dble(twotondim)
              S_avg = hbar_code * atan2(sum_im, sum_re)
              S_old = psi_im(icell)
              S_avg = S_avg + twopi_h * nint((S_old - S_avg)/twopi_h)
              psi_re(icell) = rho_avg
              psi_im(icell) = S_avg
           else if(.not.child_fluid) then
              ! wave->wave: plain complex averaging DESTROYS mass under
              ! interference (|<psi>|^2 <= <|psi|^2>, Cauchy-Schwarz).
              ! Conserve mass: amplitude from <|psi|^2>, phase from <psi>.
              rho_avg = rho_avg / dble(twotondim)
              amp2 = (sum_re**2 + sum_im**2)
              if(amp2 > 0.0d0) then
                 sfac = sqrt(rho_avg / amp2)
                 psi_re(icell) = sum_re * sfac
                 psi_im(icell) = sum_im * sfac
              else
                 ! perfect phase cancellation: keep mass, arbitrary phase
                 psi_re(icell) = sqrt(rho_avg)
                 psi_im(icell) = 0.0d0
              end if
           else
              ! fluid->fluid: (rho, S) average is linear and conservative
              psi_re(icell) = sum_re / dble(twotondim)
              psi_im(icell) = sum_im / dble(twotondim)
           end if
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_restrict
!################################################################
!################################################################
! Phase 5: de Broglie wavelength refinement criterion
! Flag cells where lambda_dB < fdm_nrefine_dB * dx for refinement
!################################################################
subroutine fdm_refine_flag(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,idim,inbor,icell_nbor
  real(dp)::dx,scale,dx_loc,grad2,psi2,lambda_dB
  real(dp)::dpsi_re,dpsi_im
  integer::nx_loc

  if(.not.use_fdm) return
  if(hbar_code <= 0.0d0) return
  ! Without HJM: only flag at levelmin (prolongated levels have parent
  ! gradients causing runaway cascade).
  ! With HJM: flag at wave levels (>= fdm_first_wave_level) where psi
  ! is actually evolved by the Schrodinger solver; skip fluid levels
  ! (handled by fdm_madelung_refine_flag).
  if(fdm_use_hjm) then
     if(ilevel < fdm_first_wave_level) return
  else
     if(ilevel > levelmin) return
  end if

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) > 0) then
           flag1(icell) = 1
        else
           psi2 = psi_re(icell)**2 + psi_im(icell)**2
           if(psi2 > 1.0d-30 .and. psi2 >= fdm_refine_rho_min) then
              grad2 = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icell_nbor)
                 if(icell_nbor > 0) then
                    dpsi_re = psi_re(icell_nbor) - psi_re(icell)
                    dpsi_im = psi_im(icell_nbor) - psi_im(icell)
                 else
                    dpsi_re = 0.0d0; dpsi_im = 0.0d0
                 end if
                 grad2 = grad2 + dpsi_re**2 + dpsi_im**2
              end do
              grad2 = grad2 / (dx_loc**2)
              if(grad2 > 0.0d0) then
                 lambda_dB = 2.0d0*acos(-1.0d0) * sqrt(psi2) / sqrt(grad2)
                 if(lambda_dB / dx_loc < dble(fdm_nrefine_dB)) then
                    flag1(icell) = 1
                 end if
              end if
           end if
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_refine_flag
!################################################################
!################################################################
! Phase 7: FDM diagnostics — total mass and min lambda_dB
! Called from adaptive_loop
!################################################################
subroutine fdm_diagnostics()
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,igrid,ind,iskip,icell,info
  real(dp)::mass_loc,mass_glob,rho_max_loc,rho_max_glob
  real(dp)::psi2,vol

  if(.not.use_fdm) return

  mass_loc = 0.0d0
  rho_max_loc = 0.0d0

  do ilevel=levelmin,nlevelmax
     vol = (0.5d0**ilevel * boxlen/dble(icoarse_max-icoarse_min+1))**3
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              if(fdm_use_hjm .and. ilevel < fdm_first_wave_level) then
                 psi2 = psi_re(icell)          ! fluid level: psi_re = rho
              else
                 psi2 = psi_re(icell)**2 + psi_im(icell)**2
              end if
              mass_loc = mass_loc + psi2 * vol
              if(psi2 > rho_max_loc) rho_max_loc = psi2
           end if
           igrid = next(igrid)
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass_loc, mass_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(rho_max_loc, rho_max_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_MAX, MPI_COMM_WORLD, info)
#else
  mass_glob = mass_loc
  rho_max_glob = rho_max_loc
#endif

  if(myid==1) then
     write(*,'(A,ES12.5,A,ES10.3)') &
          ' FDM: M_tot=', mass_glob, '  rho_max=', rho_max_glob
  end if

end subroutine fdm_diagnostics
!################################################################
!################################################################
! Per-operation mass check for debugging mass conservation.
! Prints total mass and per-level breakdown.
! label: short tag identifying the checkpoint location.
!################################################################
subroutine fdm_mass_check(label)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=*),intent(in)::label
  integer::ilevel,igrid,ind,iskip,icell,info
  real(dp)::mass8_loc,mass8_glob,mass9_loc,mass9_glob
  real(dp)::psi2,vol8,vol9,total

  mass8_loc = 0.0d0
  mass9_loc = 0.0d0
  vol8 = (0.5d0**levelmin * boxlen/dble(icoarse_max-icoarse_min+1))**3
  vol9 = vol8 / 8.0d0

  do ilevel=levelmin,min(levelmin+1,nlevelmax)
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        igrid = headl(myid, ilevel)
        do while(igrid > 0)
           icell = igrid + iskip
           if(son(icell) == 0) then
              if(fdm_use_hjm .and. ilevel < fdm_first_wave_level) then
                 psi2 = psi_re(icell)          ! fluid level: psi_re = rho
              else
                 psi2 = psi_re(icell)**2 + psi_im(icell)**2
              end if
              if(ilevel == levelmin) then
                 mass8_loc = mass8_loc + psi2 * vol8
              else
                 mass9_loc = mass9_loc + psi2 * vol9
              end if
           end if
           igrid = next(igrid)
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mass8_loc, mass8_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(mass9_loc, mass9_glob, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_COMM_WORLD, info)
#else
  mass8_glob = mass8_loc
  mass9_glob = mass9_loc
#endif

  total = mass8_glob + mass9_glob
  if(myid==1) then
     write(*,'(A,A16,A,ES16.9,A,ES12.5,A,ES12.5)') &
          ' MCHK:', label, ' tot=', total, &
          ' L8=', mass8_glob, ' L9=', mass9_glob
  end if

end subroutine fdm_mass_check

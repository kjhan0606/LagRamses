! -----------------------------------------------------------------------
! Power spectrum measurement module for cuRAMSES / lagRamses
! Computes P(k) from the dark-matter density field using NGP + FFTW3.
!   * CDM mode : NGP deposit of particle positions.
!   * FDM mode : deposit of the wavefunction density |psi|^2 from the
!                base-level AMR cells (no DM particles exist under FDM).
! Called from dump_all when dump_pk=.true.
! -----------------------------------------------------------------------
module power_spectrum_mod
  use amr_commons, only: dp, i8b
  implicit none
  private
  public :: compute_power_spectrum
contains

subroutine compute_power_spectrum(ilevel, filedir, nchar)
#ifdef USE_FFTW
  use amr_commons
  use pm_commons, only: xp, levelp, npartmax
  use poisson_commons, only: psi_re, psi_im
  use iso_c_binding
  use omp_lib
  implicit none
#ifndef WITHOUTMPI
  include "mpif.h"
#endif
  include 'fftw3.f03'

  integer, intent(in) :: ilevel
  character(LEN=*), intent(in) :: filedir, nchar

  ! Grid dimensions
  integer :: fft_Nx, fft_Ny, fft_Nz
  integer(i8b) :: N_total, N_complex
  integer(i8b), parameter :: MAX_REPLICATED_FFT_CELLS = 134217728_i8b

  ! FFTW arrays
  real(C_DOUBLE), allocatable :: rhs_3d(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: cdata(:)
  type(C_PTR) :: plan_r2c

  ! Local variables
  integer :: ipart
  integer :: ix, iy, iz
  integer :: kx_i, ky_i, kz_i
  integer :: info, ierr
  integer(i8b) :: idx_3d, idx_c
  real(dp) :: twopi, mean_count

  ! FDM field deposit
  integer :: ind, iskip, igrid_loc, icell, nx_loc
  real(dp) :: dx_lev
  real(dp), dimension(1:3) :: skip_loc
  real(dp), dimension(1:twotondim,1:3) :: xc

  ! Density arrays (for MPI reduction)
  real(dp), allocatable :: rhs_local(:)

  ! P(k) binning
  integer, parameter :: NKBIN = 50
  real(dp) :: k_fund, k_nyq, k_min_bin, k_max_bin
  real(dp) :: dk_log, k_val, pk_val, power_re, power_im, power_sq
  real(dp) :: pk_sum(NKBIN), pk_cnt(NKBIN)
  real(dp) :: k_bin_center(NKBIN), k_bin_edge(NKBIN+1)
  integer :: ibin
  real(dp) :: V_box, norm_fft
  real(dp) :: shot_noise

  ! xi(r) — two-point correlation function via Hankel transform
  real(dp) :: r_min_bin, r_max_bin, dr_log
  real(dp) :: r_bin_center(NKBIN), r_bin_edge(NKBIN+1)
  real(dp) :: xi_val(NKBIN), xi_r
  real(dp) :: k_damp, gauss_w, sinc_kr, dk_i
  integer :: jbin

  ! Physical units
  real(dp) :: boxlen_phys

  ! Output
  character(LEN=256) :: pk_filename
  integer :: ilun

  ! Particle / cell count
  integer(i8b) :: npart_loc, npart_all

  ! ================================================================
  ! Step 0: Grid dimensions (uniform grid at levelmin)
  ! ================================================================
  fft_Nx = nx * 2**ilevel
  fft_Ny = ny * 2**ilevel
  fft_Nz = nz * 2**ilevel
  N_total = int(fft_Nx, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b)
  twopi   = 2.0d0 * acos(-1.0d0)

  if(myid==1) write(*,'(A,I3,A,I5,A,I5,A,I5)') &
       ' Power spectrum: level=', ilevel, &
       ' grid=', fft_Nx, 'x', fft_Ny, 'x', fft_Nz

  ! The current estimator replicates the real-space density arrays on every
  ! MPI rank.  Refuse grids larger than 512^3 instead of exhausting node
  ! memory; production-scale spectra must be measured by distributed
  ! post-processing.
  if(N_total > MAX_REPLICATED_FFT_CELLS)then
     if(myid==1)then
        write(*,'(A,I0,A)') ' WARNING: dump_pk skipped for ', N_total, &
             ' cells; replicated in-situ FFT is limited to 512^3.'
        write(*,'(A)') '          Use distributed snapshot post-processing.'
     end if
     return
  end if

  ! ================================================================
  ! Step 1: Density deposit
  !   CDM: NGP of current particle positions (avoids stale rho array).
  !   FDM: |psi|^2 from base-level cells (no DM particles exist).
  ! ================================================================
  allocate(rhs_local(0:N_total-1))
  rhs_local = 0.0d0
  npart_loc = 0

  if(use_fdm)then
     ! Map base-level (ilevel) cells onto the uniform FFT grid. At ilevel
     ! the base grid fully tiles the box, so each cell hits one FFT node.
     nx_loc = icoarse_max - icoarse_min + 1
     skip_loc = 0.0d0
     if(ndim>0) skip_loc(1) = dble(icoarse_min)
     if(ndim>1) skip_loc(2) = dble(jcoarse_min)
     if(ndim>2) skip_loc(3) = dble(kcoarse_min)
     dx_lev = 0.5d0**ilevel
     do ind=1,twotondim
        iz = (ind-1)/4
        iy = (ind-1-4*iz)/2
        ix = (ind-1-2*iy-4*iz)
        xc(ind,1) = (dble(ix)-0.5d0)*dx_lev
        xc(ind,2) = (dble(iy)-0.5d0)*dx_lev
        xc(ind,3) = (dble(iz)-0.5d0)*dx_lev
     end do
     do ind=1,twotondim
        iskip = ncoarse + (ind-1)*ngridmax
        do igrid_loc = 1, active(ilevel)%ngrid
           icell = active(ilevel)%igrid(igrid_loc) + iskip
           ! Normalized cell-center position in [0,1)
           ix = int(((xg(active(ilevel)%igrid(igrid_loc),1)+xc(ind,1)-skip_loc(1)) &
                & / dble(nx_loc)) * dble(fft_Nx))
           iy = int(((xg(active(ilevel)%igrid(igrid_loc),2)+xc(ind,2)-skip_loc(2)) &
                & / dble(nx_loc)) * dble(fft_Ny))
           iz = int(((xg(active(ilevel)%igrid(igrid_loc),3)+xc(ind,3)-skip_loc(3)) &
                & / dble(nx_loc)) * dble(fft_Nz))
           ix = modulo(ix, fft_Nx)
           iy = modulo(iy, fft_Ny)
           iz = modulo(iz, fft_Nz)
           idx_3d = int(ix, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b) &
                  + int(iy, i8b) * int(fft_Nz, i8b) + int(iz, i8b)
           if(fdm_use_hjm .and. ilevel < fdm_first_wave_level)then
              rhs_local(idx_3d) = rhs_local(idx_3d) + psi_re(icell)
           else
              rhs_local(idx_3d) = rhs_local(idx_3d) &
                   + psi_re(icell)**2 + psi_im(icell)**2
           end if
           npart_loc = npart_loc + 1
        end do
     end do
  else
     do ipart = 1, npartmax
        if(levelp(ipart) > 0) then
           npart_loc = npart_loc + 1
           ix = int(xp(ipart, 1) * dble(fft_Nx))
           iy = int(xp(ipart, 2) * dble(fft_Ny))
           iz = int(xp(ipart, 3) * dble(fft_Nz))
           ix = modulo(ix, fft_Nx)
           iy = modulo(iy, fft_Ny)
           iz = modulo(iz, fft_Nz)
           idx_3d = int(ix, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b) &
                  + int(iy, i8b) * int(fft_Nz, i8b) + int(iz, i8b)
           rhs_local(idx_3d) = rhs_local(idx_3d) + 1.0d0
        end if
     end do
  end if

  ! ================================================================
  ! Step 2: MPI_ALLREDUCE to get global density field + count
  ! ================================================================
  allocate(rhs_3d(0:N_total-1))
  rhs_3d = 0.0d0

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_local, rhs_3d, int(N_total), &
       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(npart_loc, npart_all, 1, MPI_INTEGER8, MPI_SUM, &
       MPI_COMM_WORLD, info)
#else
  rhs_3d = rhs_local
  npart_all = npart_loc
#endif
  deallocate(rhs_local)

  ! Convert to overdensity: delta = rho(x)/rho_mean - 1
  ! CDM: rho_mean = <particle count per cell>.
  ! FDM: rho_mean = < |psi|^2 > averaged over all cells.
  if(use_fdm)then
     mean_count = sum(rhs_3d) / dble(N_total)
  else
     mean_count = dble(npart_all) / dble(N_total)
  end if
  if(mean_count > 0.0d0)then
     rhs_3d = rhs_3d / mean_count - 1.0d0
  else
     rhs_3d = 0.0d0
  end if

  ! Diagnostic output
  if(myid == 1) then
     if(use_fdm)then
        write(*,'(A,I15)')    '   P(k) N_cell        = ', npart_all
     else
        write(*,'(A,I15)')    '   P(k) N_part        = ', npart_all
     end if
     write(*,'(A,ES15.8)') '   P(k) boxlen_ini    = ', boxlen_ini
     write(*,'(A,ES15.8)') '   P(k) mean_count    = ', mean_count
     write(*,'(A,ES15.8)') '   P(k) delta mean    = ', sum(rhs_3d) / dble(N_total)
     write(*,'(A,ES15.8)') '   P(k) delta rms     = ', &
          sqrt(sum(rhs_3d**2) / dble(N_total))
  end if

  ! ================================================================
  ! Step 3: FFT on rank 0 only (small grid, N<=256^3)
  ! ================================================================
  N_complex = int(fft_Nx, i8b) * int(fft_Ny, i8b) * int(fft_Nz/2+1, i8b)
  allocate(cdata(0:N_complex-1))

  if(myid == 1) then
     ! Initialize FFTW threads
     ierr = fftw_init_threads()
     call fftw_plan_with_nthreads(int(omp_get_max_threads(), C_INT))

     ! Create R2C plan
     plan_r2c = fftw_plan_dft_r2c_3d( &
          int(fft_Nx, C_INT), int(fft_Ny, C_INT), int(fft_Nz, C_INT), &
          rhs_3d, cdata, FFTW_ESTIMATE)

     ! Execute FFT
     call fftw_execute_dft_r2c(plan_r2c, rhs_3d, cdata)
     call fftw_destroy_plan(plan_r2c)
  end if

  deallocate(rhs_3d)

  ! ================================================================
  ! Step 4: P(k) binning (rank 0 only)
  ! ================================================================
  ! Physical box size in Mpc/h
  boxlen_phys = boxlen_ini  ! already in Mpc/h (set in init_cosmo)

  ! Fundamental and Nyquist wavenumbers
  k_fund = twopi / boxlen_phys  ! h/Mpc
  k_nyq  = k_fund * dble(fft_Nx) / 2.0d0

  ! Logarithmic bins from k_fund to k_nyq
  k_min_bin = k_fund
  k_max_bin = k_nyq
  dk_log = (log10(k_max_bin) - log10(k_min_bin)) / dble(NKBIN)

  do ibin = 1, NKBIN+1
     k_bin_edge(ibin) = 10.0d0**(log10(k_min_bin) + dble(ibin-1) * dk_log)
  end do
  do ibin = 1, NKBIN
     k_bin_center(ibin) = sqrt(k_bin_edge(ibin) * k_bin_edge(ibin+1))
  end do

  pk_sum = 0.0d0
  pk_cnt = 0.0d0

  ! Volume and normalization
  V_box = boxlen_phys**3  ! (Mpc/h)^3
  norm_fft = V_box / dble(N_total)**2

  if(myid == 1) then

     do kx_i = 0, fft_Nx - 1
        do ky_i = 0, fft_Ny - 1
           do kz_i = 0, fft_Nz/2
              ! Skip DC mode
              if(kx_i == 0 .and. ky_i == 0 .and. kz_i == 0) cycle

              idx_c = int(kx_i, i8b) * int(fft_Ny, i8b) * int(fft_Nz/2+1, i8b) &
                    + int(ky_i, i8b) * int(fft_Nz/2+1, i8b) + int(kz_i, i8b)

              ! Physical k in h/Mpc
              ix = kx_i; if(ix > fft_Nx/2) ix = ix - fft_Nx
              iy = ky_i; if(iy > fft_Ny/2) iy = iy - fft_Ny
              iz = kz_i  ! R2C: kz only goes to Nz/2

              k_val = k_fund * sqrt(dble(ix)**2 + dble(iy)**2 + dble(iz)**2)

              ! Find bin
              if(k_val < k_min_bin .or. k_val >= k_max_bin) cycle
              ibin = int((log10(k_val) - log10(k_min_bin)) / dk_log) + 1
              if(ibin < 1 .or. ibin > NKBIN) cycle

              ! Power: |delta_hat(k)|^2
              power_re = dble(cdata(idx_c))
              power_im = aimag(cdata(idx_c))
              power_sq = power_re**2 + power_im**2

              ! Weight: modes at kz=0 and kz=Nz/2 appear once,
              ! others have a conjugate partner (counted by R2C)
              if(kz_i == 0 .or. kz_i == fft_Nz/2) then
                 pk_sum(ibin) = pk_sum(ibin) + power_sq * norm_fft
                 pk_cnt(ibin) = pk_cnt(ibin) + 1.0d0
              else
                 pk_sum(ibin) = pk_sum(ibin) + 2.0d0 * power_sq * norm_fft
                 pk_cnt(ibin) = pk_cnt(ibin) + 2.0d0
              end if
           end do
        end do
     end do
  end if

  deallocate(cdata)

  ! ================================================================
  ! Step 5: Shot noise correction
  !   FDM density is a field: no particle Poisson shot noise.
  ! ================================================================
  if(use_fdm)then
     shot_noise = 0.0d0
  else if(npart_all > 0) then
     shot_noise = V_box / dble(npart_all)  ! (Mpc/h)^3
  else
     shot_noise = 0.0d0
  end if

  ! ================================================================
  ! Step 6: Two-point correlation xi(r) via Hankel transform
  ! xi(r) = 1/(2pi^2) * int [P(k)-Pshot] * W(k) * k^2 * sin(kr)/(kr) dk
  ! W(k) = exp(-(k/k_damp)^2)  Gaussian to suppress ringing
  ! ================================================================
  r_min_bin = boxlen_phys / dble(fft_Nx)   ! grid cell size [Mpc/h]
  r_max_bin = boxlen_phys / 2.0d0          ! half box [Mpc/h]
  dr_log = (log10(r_max_bin) - log10(r_min_bin)) / dble(NKBIN)

  do ibin = 1, NKBIN+1
     r_bin_edge(ibin) = 10.0d0**(log10(r_min_bin) + dble(ibin-1) * dr_log)
  end do
  do ibin = 1, NKBIN
     r_bin_center(ibin) = sqrt(r_bin_edge(ibin) * r_bin_edge(ibin+1))
  end do

  k_damp = k_nyq  ! Gaussian rolls off at Nyquist
  xi_val = 0.0d0

  if(myid == 1) then
     do jbin = 1, NKBIN
        xi_r = r_bin_center(jbin)
        do ibin = 1, NKBIN
           if(pk_cnt(ibin) < 0.5d0) cycle
           pk_val = pk_sum(ibin) / pk_cnt(ibin) - shot_noise
           dk_i   = k_bin_edge(ibin+1) - k_bin_edge(ibin)
           gauss_w = exp(-(k_bin_center(ibin) / k_damp)**2)
           sinc_kr = sin(k_bin_center(ibin) * xi_r) &
                   / (k_bin_center(ibin) * xi_r)
           xi_val(jbin) = xi_val(jbin) &
                + pk_val * gauss_w * k_bin_center(ibin)**2 * sinc_kr * dk_i
        end do
        ! 1/(2*pi^2) = 2/twopi^2
        xi_val(jbin) = xi_val(jbin) * 2.0d0 / twopi**2
     end do
  end if

  ! ================================================================
  ! Step 7: Write P(k) + xi(r) output file (rank 0 only)
  ! ================================================================
  if(myid == 1) then
     pk_filename = TRIM(filedir)//'pk_'//TRIM(nchar)//'.dat'
     ilun = 42

     open(unit=ilun, file=pk_filename, form='formatted', status='replace')
     write(ilun, '(A,ES15.8)') '# Power spectrum at a_exp = ', aexp
     write(ilun, '(A,ES15.8)') '# boxlen (Mpc/h) = ', boxlen_phys
     write(ilun, '(A,I10)')    '# N_grid = ', fft_Nx
     if(use_fdm)then
        write(ilun, '(A,I15)') '# N_cell (FDM |psi|^2) = ', npart_all
     else
        write(ilun, '(A,I15)') '# N_part = ', npart_all
     end if
     write(ilun, '(A,ES15.8)') '# shot_noise (Mpc/h)^3 = ', shot_noise
     write(ilun, '(A,ES15.8)') '# k_damp (h/Mpc) = ', k_damp
     write(ilun, '(A)') '# k [h/Mpc]    P(k) [(Mpc/h)^3]    r [Mpc/h]    xi(r)    Nmodes'

     do ibin = 1, NKBIN
        if(pk_cnt(ibin) > 0.5d0) then
           pk_val = pk_sum(ibin) / pk_cnt(ibin)
           write(ilun, '(5ES16.8)') k_bin_center(ibin), pk_val, &
                r_bin_center(ibin), xi_val(ibin), pk_cnt(ibin)
        end if
     end do
     close(ilun)

     write(*,'(A,A)') '   P(k) written to ', trim(pk_filename)
  end if

#else
  ! No FFTW — stub
  use amr_commons, only: myid
  implicit none
  integer, intent(in) :: ilevel
  character(LEN=*), intent(in) :: filedir, nchar
  if(myid==1) write(*,*) 'WARNING: dump_pk=T but not compiled with USE_FFTW'
#endif

end subroutine compute_power_spectrum

end module power_spectrum_mod

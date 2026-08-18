! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ fine and coarse MG levels
!   * communicator building
!   * MPI routines
!   * helper functions
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------

#ifdef FDMDEBUG
module mg_omp_profile_m
   use omp_lib, only: omp_get_wtime
   implicit none

   integer, parameter :: MGP_NPHASE = 14
   integer, parameter :: MGP_TOTAL=1, MGP_BUILD=2, MGP_MASK=3
   integer, parameter :: MGP_SCAN=4, MGP_BCRHS=5, MGP_FINE_SMOOTH=6
   integer, parameter :: MGP_FINE_RESID=7, MGP_FINE_COMM=8
   integer, parameter :: MGP_COARSE_SMOOTH=9, MGP_COARSE_RESID=10
   integer, parameter :: MGP_MG_FORWARD=11, MGP_MG_REVERSE=12
   integer, parameter :: MGP_RESTRICT=13, MGP_INTERP=14

   real(kind=8), save :: mgp_wall(MGP_NPHASE)=0d0
   real(kind=8), save :: mgp_cpu(MGP_NPHASE)=0d0
   integer(kind=8), save :: mgp_calls(MGP_NPHASE)=0_8
   integer(kind=8), save :: mgp_work(MGP_NPHASE)=0_8
   integer, save :: mgp_last_report=-1

contains

   subroutine mgp_start(wall_start,cpu_start)
      real(kind=8), intent(out) :: wall_start,cpu_start
      wall_start=omp_get_wtime()
      call cpu_time(cpu_start)
   end subroutine mgp_start

   subroutine mgp_stop(iphase,wall_start,cpu_start,nwork)
      integer, intent(in) :: iphase
      real(kind=8), intent(in) :: wall_start,cpu_start
      integer(kind=8), intent(in), optional :: nwork
      real(kind=8) :: wall_stop,cpu_stop

      wall_stop=omp_get_wtime()
      call cpu_time(cpu_stop)
      mgp_wall(iphase)=mgp_wall(iphase)+wall_stop-wall_start
      mgp_cpu(iphase)=mgp_cpu(iphase)+cpu_stop-cpu_start
      mgp_calls(iphase)=mgp_calls(iphase)+1_8
      if(present(nwork))mgp_work(iphase)=mgp_work(iphase)+nwork
   end subroutine mgp_stop

   function mgp_label(iphase) result(label)
      integer, intent(in) :: iphase
      character(len=24) :: label

      select case(iphase)
      case(MGP_TOTAL);         label='multigrid total'
      case(MGP_BUILD);         label='build parent comms'
      case(MGP_MASK);          label='make fine mask'
      case(MGP_SCAN);          label='coarse scan flags'
      case(MGP_BCRHS);         label='fine BC RHS'
      case(MGP_FINE_SMOOTH);   label='fine smoother'
      case(MGP_FINE_RESID);    label='fine residual'
      case(MGP_FINE_COMM);     label='fine ghost exchange'
      case(MGP_COARSE_SMOOTH); label='coarse smoother'
      case(MGP_COARSE_RESID);  label='coarse residual'
      case(MGP_MG_FORWARD);    label='MG forward exchange'
      case(MGP_MG_REVERSE);    label='MG reverse exchange'
      case(MGP_RESTRICT);      label='restriction'
      case(MGP_INTERP);        label='interpolation'
      case default;            label='other/setup'
      end select
   end function mgp_label

   subroutine mgp_report(step,myid,ncpu,nthreads)
      integer, intent(in) :: step,myid,ncpu,nthreads
      real(kind=8) :: wall_sum(MGP_NPHASE),wall_max(MGP_NPHASE)
      real(kind=8) :: cpu_sum(MGP_NPHASE),wall_avg,eff,balance
      integer(kind=8) :: calls_sum(MGP_NPHASE),work_sum(MGP_NPHASE)
      integer :: iphase,info
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      if(step<5 .or. mgp_last_report>=0)return
#ifndef WITHOUTMPI
      call MPI_REDUCE(mgp_wall,wall_sum,MGP_NPHASE,MPI_DOUBLE_PRECISION, &
           MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(mgp_wall,wall_max,MGP_NPHASE,MPI_DOUBLE_PRECISION, &
           MPI_MAX,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(mgp_cpu,cpu_sum,MGP_NPHASE,MPI_DOUBLE_PRECISION, &
           MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(mgp_calls,calls_sum,MGP_NPHASE,MPI_INTEGER8, &
           MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(mgp_work,work_sum,MGP_NPHASE,MPI_INTEGER8, &
           MPI_SUM,0,MPI_COMM_WORLD,info)
#else
      wall_sum=mgp_wall
      wall_max=mgp_wall
      cpu_sum=mgp_cpu
      calls_sum=mgp_calls
      work_sum=mgp_work
#endif
      if(myid==1)then
         write(*,'(A,I0,A,I0,A)') &
              ' MG_OMP_PROFILE step=',step,' threads=',nthreads, &
              ' (cumulative; CPU/(threads*wall))'
         write(*,'(A)') &
              ' MGPROF phase                     wall-avg wall-max  omp-eff balance calls/rank work/rank'
         do iphase=1,MGP_NPHASE
            if(calls_sum(iphase)==0_8)cycle
            wall_avg=wall_sum(iphase)/dble(ncpu)
            eff=0d0
            if(wall_sum(iphase)>0d0)eff=cpu_sum(iphase)/ &
                 (dble(max(1,nthreads))*wall_sum(iphase))
            balance=0d0
            if(wall_max(iphase)>0d0)balance=wall_avg/wall_max(iphase)
            write(*,'(A,A24,2F10.4,2F9.3,2F12.1)') ' MGPROF ', &
                 mgp_label(iphase),wall_avg,wall_max(iphase),eff,balance, &
                 dble(calls_sum(iphase))/dble(ncpu), &
                 dble(work_sum(iphase))/dble(ncpu)
         end do
      end if
      mgp_last_report=step
   end subroutine mgp_report

end module mg_omp_profile_m

module fftw_omp_profile_m
   use omp_lib, only: omp_get_wtime
   implicit none

   integer, parameter :: FFTP_NPHASE=12, FFTP_NPATH=2
   integer, parameter :: FFTP_TOTAL=1, FFTP_SETUP=2, FFTP_PARTNER=3
   integer, parameter :: FFTP_COUNT=4, FFTP_COUNTCOMM=5, FFTP_PACK=6
   integer, parameter :: FFTP_FWDCOMM=7, FFTP_R2C=8, FFTP_GREEN=9
   integer, parameter :: FFTP_C2R=10, FFTP_REVCOMM=11, FFTP_SCATTER=12

   real(kind=8), save :: fftp_wall(FFTP_NPHASE,FFTP_NPATH)=0d0
   real(kind=8), save :: fftp_cpu(FFTP_NPHASE,FFTP_NPATH)=0d0
   integer(kind=8), save :: fftp_calls(FFTP_NPHASE,FFTP_NPATH)=0_8
   integer(kind=8), save :: fftp_work(FFTP_NPHASE,FFTP_NPATH)=0_8
   integer, save :: fftp_last_report=-1

contains

   subroutine fftp_start(wall_start,cpu_start)
      real(kind=8), intent(out) :: wall_start,cpu_start
      wall_start=omp_get_wtime()
      call cpu_time(cpu_start)
   end subroutine fftp_start

   subroutine fftp_stop(ipath,iphase,wall_start,cpu_start,nwork)
      integer, intent(in) :: ipath,iphase
      real(kind=8), intent(in) :: wall_start,cpu_start
      integer(kind=8), intent(in), optional :: nwork
      real(kind=8) :: wall_stop,cpu_stop

      wall_stop=omp_get_wtime()
      call cpu_time(cpu_stop)
      fftp_wall(iphase,ipath)=fftp_wall(iphase,ipath)+wall_stop-wall_start
      fftp_cpu(iphase,ipath)=fftp_cpu(iphase,ipath)+cpu_stop-cpu_start
      fftp_calls(iphase,ipath)=fftp_calls(iphase,ipath)+1_8
      if(present(nwork))fftp_work(iphase,ipath)=fftp_work(iphase,ipath)+nwork
   end subroutine fftp_stop

   function fftp_label(iphase) result(label)
      integer, intent(in) :: iphase
      character(len=20) :: label

      select case(iphase)
      case(FFTP_TOTAL);     label='total'
      case(FFTP_SETUP);     label='setup/plan'
      case(FFTP_PARTNER);   label='partner discovery'
      case(FFTP_COUNT);     label='cell count'
      case(FFTP_COUNTCOMM); label='count exchange'
      case(FFTP_PACK);      label='gather/pack'
      case(FFTP_FWDCOMM);   label='forward exchange'
      case(FFTP_R2C);       label='forward FFT'
      case(FFTP_GREEN);     label='spectral kernel'
      case(FFTP_C2R);       label='inverse FFT'
      case(FFTP_REVCOMM);   label='reverse exchange'
      case(FFTP_SCATTER);   label='normalize/scatter'
      case default;         label='other'
      end select
   end function fftp_label

   subroutine fftp_report(step,myid,ncpu,nthreads)
      integer, intent(in) :: step,myid,ncpu,nthreads
      real(kind=8) :: wall_sum(FFTP_NPHASE,FFTP_NPATH)
      real(kind=8) :: wall_max(FFTP_NPHASE,FFTP_NPATH)
      real(kind=8) :: cpu_sum(FFTP_NPHASE,FFTP_NPATH)
      real(kind=8) :: wall_avg,eff,balance
      integer(kind=8) :: calls_sum(FFTP_NPHASE,FFTP_NPATH)
      integer(kind=8) :: work_sum(FFTP_NPHASE,FFTP_NPATH)
      integer :: iphase,ipath,info
      character(len=11) :: path_label
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      if(step<3 .or. fftp_last_report>=0)return
#ifndef WITHOUTMPI
      call MPI_REDUCE(fftp_wall,wall_sum,FFTP_NPHASE*FFTP_NPATH, &
           MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(fftp_wall,wall_max,FFTP_NPHASE*FFTP_NPATH, &
           MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(fftp_cpu,cpu_sum,FFTP_NPHASE*FFTP_NPATH, &
           MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(fftp_calls,calls_sum,FFTP_NPHASE*FFTP_NPATH, &
           MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,info)
      call MPI_REDUCE(fftp_work,work_sum,FFTP_NPHASE*FFTP_NPATH, &
           MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,info)
#else
      wall_sum=fftp_wall
      wall_max=fftp_wall
      cpu_sum=fftp_cpu
      calls_sum=fftp_calls
      work_sum=fftp_work
#endif
      if(myid==1)then
         write(*,'(A,I0,A,I0,A)') &
              ' FFTW_OMP_PROFILE step=',step,' threads=',nthreads, &
              ' (cumulative; CPU/(threads*wall))'
         write(*,'(A)') &
              ' FFTPROF path        phase                 wall-avg wall-max  omp-eff balance calls/rank work/rank'
         do ipath=1,FFTP_NPATH
            path_label=merge('distributed','replicated ',ipath==2)
            do iphase=1,FFTP_NPHASE
               if(calls_sum(iphase,ipath)==0_8)cycle
               wall_avg=wall_sum(iphase,ipath)/dble(ncpu)
               eff=0d0
               if(wall_sum(iphase,ipath)>0d0)eff=cpu_sum(iphase,ipath)/ &
                    (dble(max(1,nthreads))*wall_sum(iphase,ipath))
               balance=0d0
               if(wall_max(iphase,ipath)>0d0)balance=wall_avg/wall_max(iphase,ipath)
               write(*,'(A,A11,1X,A20,2F10.4,2F9.3,2F12.1)') &
                    ' FFTPROF ',path_label,fftp_label(iphase), &
                    wall_avg,wall_max(iphase,ipath),eff,balance, &
                    dble(calls_sum(iphase,ipath))/dble(ncpu), &
                    dble(work_sum(iphase,ipath))/dble(ncpu)
            end do
         end do
      end if
      fftp_last_report=step
   end subroutine fftp_report

end module fftw_omp_profile_m
#endif


! ------------------------------------------------------------------------
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid_fine(ilevel,icount)
   use amr_commons
   use poisson_commons
   use poisson_parameters
#ifdef FDMDEBUG
   use mg_omp_profile_m
   use omp_lib, only: omp_get_max_threads
#endif
#ifdef HYDRO_CUDA
   use poisson_cuda_interface
   use iso_c_binding
   use cuda_commons, only: cuda_pool_is_initialized_c
#endif

   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel,icount

   interface
      subroutine cmp_residual_mg_fine(ilevel, norm2)
         use amr_commons, only: dp
         integer, intent(in) :: ilevel
         real(dp), intent(out), optional :: norm2
      end subroutine
      subroutine cmp_residual_norm2_fine(ilevel, norm2)
         use amr_commons, only: dp
         integer, intent(in) :: ilevel
         real(dp), intent(out) :: norm2
      end subroutine
   end interface

   real(dp), parameter :: SAFE_FACTOR = 0.5

   integer  :: ifine, i, iter, info, icpu
   real(kind=8) :: res_norm2, i_res_norm2, i_res_norm2_tot, res_norm2_tot
   real(kind=8) :: debug_norm2, debug_norm2_tot
   real(kind=8) :: err, last_err

   logical :: allmasked, allmasked_tot, use_restored_phi, mg_failed
#ifdef FDMDEBUG
   real(kind=8) :: mgp_wall_start,mgp_cpu_start,mgp_phase_wall,mgp_phase_cpu
#endif

   ! FFT direct solve variables (shared by FFTW3 and cuFFT)
   logical :: is_uniform_fft
   integer(i8b) :: expected_grids

   ! GPU Poisson MG variables
   logical :: use_mg_gpu, use_ri_gpu
#ifdef HYDRO_CUDA
   integer(c_long_long) :: ncell_tot_c
   integer(i8b) :: ncell_tot
   integer :: safe_int, ri_flag
   real(kind=8) :: dx_mg, dx2_mg, oneoverdx2_mg, dx2_norm_mg
   real(kind=8) :: gpu_norm2, dummy_norm2
#endif

   if(gravity_type>0)return
   if(numbtot(1,ilevel)==0)return

#ifdef FDMDEBUG
   call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

   if(verbose) print '(A,I2)','Entering fine multigrid at level ',ilevel

   ! ---------------------------------------------------------------------
   ! FFTW3 CPU direct-solve FAST PATH (fully-uniform periodic level).
   ! The FFT solver needs only rho (already built by rho_fine) and writes
   ! phi directly via scatter, so we bypass the ENTIRE multigrid hierarchy
   ! setup (mask, parent comms, neighbor grids, BC-rhs) and its teardown.
   ! That per-call build+teardown -- not the FFT itself -- dominates the
   ! base-level Poisson cost. numbtot is global (identical on all ranks),
   ! so the uniform test needs no MPI consensus. Refined levels fall through
   ! to the multigrid V-cycle below.
   ! ---------------------------------------------------------------------
#ifdef USE_FFTW
   if(use_fftw) then
      expected_grids = int(nx,i8b)*int(ny,i8b)*int(nz,i8b)*8_i8b**(ilevel-1)
      if(numbtot(1,ilevel) == expected_grids) then
         call fftw_poisson_solve_uniform(ilevel, icount)
         call make_virtual_fine_dp_mg_profile(phi(1), ilevel) ! ghost update for force_fine
         if(myid==1) write(*,'(A,I5,A)') &
              '   ==> Level=',ilevel,' FFT direct solve DONE (fast path)'
         return
      end if
   end if
#endif

   ! ---------------------------------------------------------------------
   ! Prepare first guess, mask and BCs at finest level
   ! ---------------------------------------------------------------------

   use_restored_phi=.false.
   if(allocated(phi_restart_available))then
      use_restored_phi=phi_restart_available(ilevel)
   endif
   if(use_restored_phi)then
      phi_restart_available(ilevel)=.false.
      if(myid==1)write(*,'(A,I3)') &
           ' Poisson warm start from restored phi at level ',ilevel
   else if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)         ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)       ! Fill with simple initial guess
   endif
   call make_virtual_fine_dp_mg_profile(phi(1),ilevel) ! Update boundaries
   call make_boundary_phi(ilevel)           ! Update physical boundaries

   call make_fine_mask  (ilevel)            ! Fill the fine mask
   call make_virtual_fine_dp_mg_profile(f(:,3),ilevel) ! Communicate mask
   call make_boundary_mask(ilevel)          ! Set mask to -1 in phys bounds

   ! The flat neighbor cache is consumed only by the CUDA MG upload path.
   ! CPU MG performs its own Morton lookup and must not pay this serial setup
   ! or retain the otherwise-unused (twondim+1)*ngrid integer array.
#ifdef HYDRO_CUDA
   call precompute_nbor_grid_fine(ilevel)
#endif

   call make_fine_bc_rhs(ilevel,icount)            ! Fill BC-modified RHS

   ! ---------------------------------------------------------------------
   ! Build communicators up
   ! ---------------------------------------------------------------------

   ! @ finer level
   call build_parent_comms_mg(active(ilevel),ilevel)
   ! @ coarser levels
   do ifine=(ilevel-1),2,-1
      call build_parent_comms_mg(active_mg(myid,ifine),ifine)
   end do

   ! ---------------------------------------------------------------------
   ! Restrict mask up, then set scan flag
   ! ---------------------------------------------------------------------
   ! @ finer level

   if(ilevel>1) then
      ! Restrict and communicate mask
#ifdef FDMDEBUG
      call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
      call restrict_mask_fine_reverse(ilevel)
#ifdef FDMDEBUG
      call mgp_stop(MGP_RESTRICT,mgp_phase_wall,mgp_phase_cpu, &
           int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
      call make_reverse_mg_dp(4,ilevel-1)
      call make_virtual_mg_dp(4,ilevel-1)

      ! Convert volume fraction to mask value
      do icpu=1,ncpu
         if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
         active_mg(icpu,ilevel-1)%u(:,4)=2d0*active_mg(icpu,ilevel-1)%u(:,4)-1d0
      end do

      ! Check active mask state
      if(active_mg(myid,ilevel-1)%ngrid>0) then
         allmasked=(maxval(active_mg(myid,ilevel-1)%u(:,4))<=0d0)
      else
         allmasked=.true.
      end if

      ! Allreduce on mask state
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(allmasked, allmasked_tot, 1, MPI_LOGICAL, &
           & MPI_LAND, MPI_COMM_WORLD, info)
      allmasked=allmasked_tot
#endif
   else
      allmasked=.true.
   endif

   ! @ coarser levels
   ! Restrict mask and compute levelmin_mg in the process
   if (.not. allmasked) then
      levelmin_mg=1
      do ifine=(ilevel-1),2,-1

         ! Restrict and communicate mask
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call restrict_mask_coarse_reverse(ifine)
#ifdef FDMDEBUG
         call mgp_stop(MGP_RESTRICT,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifine)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
         call make_reverse_mg_dp(4,ifine-1)
         call make_virtual_mg_dp(4,ifine-1)

         ! Convert volume fraction to mask value
         do icpu=1,ncpu
            if(active_mg(icpu,ifine-1)%ngrid==0) cycle
            active_mg(icpu,ifine-1)%u(:,4)=2d0*active_mg(icpu,ifine-1)%u(:,4)-1d0
         end do

         ! Check active mask state
         if(active_mg(myid,ifine-1)%ngrid>0) then
            allmasked=(maxval(active_mg(myid,ifine-1)%u(:,4))<=0d0)
         else
            allmasked=.true.
         end if

         ! Allreduce on mask state
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(allmasked,allmasked_tot,1,MPI_LOGICAL, &
                 & MPI_LAND,MPI_COMM_WORLD,info)
         allmasked=allmasked_tot
#endif

         if(allmasked) then ! Coarser level is fully masked: stop here
            levelmin_mg=ifine
            exit
         end if
      end do
   else
      levelmin_mg=ilevel
   end if
   if(nboundary>0)levelmin_mg=max(levelmin_mg,2)

   ! Update flag with scan flag
   call set_scan_flag_fine(ilevel)
   do ifine=levelmin_mg,ilevel-1
      call set_scan_flag_coarse_omp(ifine)
   end do

   ! Precompute neighbor grids for coarse MG levels
   if(ilevel>1 .and. levelmin_mg < ilevel) then
      call precompute_nbor_grid_coarse(levelmin_mg, ilevel-1)
   end if

   ! ---------------------------------------------------------------------
   ! Initiate solve at fine level
   ! ---------------------------------------------------------------------

   ! -----------------------------------------------------------------
   ! FFT direct solve for fully uniform level (periodic BC)
   ! Priority: FFTW3 (CPU) > cuFFT (GPU) > MG V-cycle
   ! -----------------------------------------------------------------
   use_mg_gpu = .false.
   use_ri_gpu = .false.
   is_uniform_fft = .false.

#ifdef USE_FFTW
   ! FFTW3 CPU direct solve (highest priority if enabled)
   if(use_fftw) then
      expected_grids = int(nx,i8b)*int(ny,i8b)*int(nz,i8b) &
                     * 8_i8b**(ilevel-1)
      if(numbtot(1,ilevel) == expected_grids) is_uniform_fft = .true.
   end if
#endif

#ifdef HYDRO_CUDA
   ! cuFFT GPU direct solve (if FFTW not handling it)
   if(.not. is_uniform_fft .and. gpu_fft .and. cuda_pool_is_initialized_c() /= 0) then
      expected_grids = int(nx,i8b)*int(ny,i8b)*int(nz,i8b) &
                     * 8_i8b**(ilevel-1)
      if(numbtot(1,ilevel) == expected_grids) is_uniform_fft = .true.
   end if
   ! Force consensus across MPI ranks (MIN: if any rank can't use FFT direct,
   ! all fall back to MG — otherwise ranks diverge into incompatible collectives).
#ifndef WITHOUTMPI
   ri_flag = 0; if(is_uniform_fft) ri_flag = 1
   call MPI_ALLREDUCE(MPI_IN_PLACE, ri_flag, 1, &
        MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, info)
   is_uniform_fft = (ri_flag == 1)
#endif

   ! GPU MG setup: controlled by gpu_poisson namelist parameter
   ! All ranks must enter this block together (contains MPI collectives)
   if(gpu_poisson .and. .not. is_uniform_fft .and. cuda_pool_is_initialized_c() /= 0) then
      ncell_tot = int(ncoarse,i8b) + int(twotondim,i8b)*int(ngridmax,i8b)
      ncell_tot_c = int(ncell_tot, c_long_long)
      dx_mg  = 0.5d0**ilevel
      dx2_mg = dx_mg*dx_mg
      oneoverdx2_mg = 1.0d0/dx2_mg
      dx2_norm_mg = dx_mg**(ndim)
      ! Upload only if this rank has active grids
      if(active(ilevel)%ngrid > 0) then
         call cuda_mg_upload_c( &
              phi, f, flag2, ncell_tot_c, &
              nbor_grid_fine, active(ilevel)%igrid, &
              int(active(ilevel)%ngrid, c_int))
         use_mg_gpu = (cuda_mg_is_ready_c() /= 0)
      else
         use_mg_gpu = .false.
      end if
      ! Synchronize use_mg_gpu: all ranks must agree (V-cycle has MPI comms)
#ifndef WITHOUTMPI
      ri_flag = 0; if(use_mg_gpu) ri_flag = 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, ri_flag, 1, &
           MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, info)
      use_mg_gpu = (ri_flag == 1)
#endif
      if(use_mg_gpu) call build_mg_halo_indices(ilevel)
      ! Setup GPU restrict/interp if MG GPU is ready and coarse levels exist
      if(use_mg_gpu .and. ilevel > 1) then
         call precompute_mg_gpu_restrict_interp(ilevel)
         use_ri_gpu = (cuda_mg_ri_is_ready_c() /= 0)
      end if
      ! Synchronize use_ri_gpu across all ranks
#ifndef WITHOUTMPI
      ri_flag = 0; if(use_ri_gpu) ri_flag = 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, ri_flag, 1, &
           MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, info)
      use_ri_gpu = (ri_flag == 1)
#endif
      if(myid==1) write(*,'(A,I3,A,L1,A,L1,A,I12,A,I12,A,I12,A,I10)') &
           ' MG GPU: level=',ilevel,' ready=',use_mg_gpu, &
           ' ri=',use_ri_gpu, &
           ' ncell=',ncell_tot,' nco=',ncoarse,' ngm=',ngridmax, &
           ' ngrid=',active(ilevel)%ngrid
   end if
#endif

   if(is_uniform_fft) then
#ifdef USE_FFTW
      if(use_fftw) then
         call fftw_poisson_solve_uniform(ilevel, icount)
      else
#endif
#ifdef HYDRO_CUDA
         call fft_poisson_solve_uniform(ilevel, icount)
#endif
#ifdef USE_FFTW
      end if
#endif
      ! phi is already scattered by the FFT solver
      call make_virtual_fine_dp_mg_profile(phi(1), ilevel)
#ifdef HYDRO_CUDA
      ! Cleanup GPU state
      if(use_mg_gpu) then
         call cuda_mg_halo_free_c()
         call cuda_mg_free_c()
      end if
      if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
      if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
      if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
      if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)
      mg_halo_n_emit = 0
      mg_halo_n_recv = 0
      if(use_ri_gpu) call cuda_mg_ri_free_c()
      if(allocated(mg_ri_flat_offset)) deallocate(mg_ri_flat_offset)
      if(allocated(mg_ri_coarse_rhs))  deallocate(mg_ri_coarse_rhs)
      if(allocated(mg_ri_coarse_phi))  deallocate(mg_ri_coarse_phi)
      mg_ri_total_coarse = 0
#endif
      ! Free precomputed neighbors
      if(allocated(nbor_grid_fine)) deallocate(nbor_grid_fine)
      if(ilevel>1 .and. levelmin_mg < ilevel) then
         call cleanup_nbor_grid_coarse(levelmin_mg, ilevel-1)
      end if
      ! Cleanup MG levels
      do ifine=1,ilevel-1
         call cleanup_mg_level(ifine)
      end do
      if(myid==1) write(*,'(A,I5,A)') &
           '   ==> Level=',ilevel,' FFT direct solve DONE'
#ifdef HYDRO_CUDA
      ! Print cuFFT phase timers immediately after first solve (reviewer C11).
      ! Doing it here (not at coarse-step end) ensures we capture data even if
      ! a later phase of the step crashes before adaptive_loop's print.
      if(myid==1 .and. .not. use_fftw) call cuda_fft_print_timers_c(int(myid, c_int))
#endif
      return
   end if

   iter = 0
   err = 1.0d0
   main_iteration_loop: do
      iter=iter+1

      ! Pre-smoothing
      do i=1,ngs_fine
#ifdef HYDRO_CUDA
         if(use_mg_gpu) then
            safe_int = 0
            if(safe_mode(ilevel)) safe_int = 1
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 0, safe_int)
            if(.not.mg_merged_rb) call make_virtual_fine_dp_gpu(ilevel)
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 1, safe_int)
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
#ifdef FDMDEBUG
            call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
            call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
#ifdef FDMDEBUG
            call mgp_stop(MGP_FINE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
                 int(active(ilevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
            if(.not.mg_merged_rb) call make_virtual_fine_dp_mg_profile(phi(1),ilevel) ! Communicate phi (Red)
#ifdef FDMDEBUG
            call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
            call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
#ifdef FDMDEBUG
            call mgp_stop(MGP_FINE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
                 int(active(ilevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
            call make_virtual_fine_dp_mg_profile(phi(1),ilevel) ! Communicate phi (Black)
#ifdef HYDRO_CUDA
         end if
#endif
      end do

      ! Compute residual and restrict into upper level RHS
#ifdef HYDRO_CUDA
      if(use_mg_gpu) then
         gpu_norm2 = 0.0d0
         if(iter==1) then
            call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), &
                 oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
                 gpu_norm2, 1)
            i_res_norm2 = gpu_norm2
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                    & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
            i_res_norm2=i_res_norm2_tot
#endif
         else
            call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), &
                 oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
                 dummy_norm2, 0)
         end if
         if(.not. use_ri_gpu) call cuda_mg_download_f1_c(f, ncell_tot_c)
      else
#endif
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call cmp_residual_mg_fine(ilevel)
#ifdef FDMDEBUG
         call mgp_stop(MGP_FINE_RESID,mgp_phase_wall,mgp_phase_cpu, &
              int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
#ifdef HYDRO_CUDA
      end if
#endif
      call make_virtual_fine_dp_mg_profile(f(1,1),ilevel) ! communicate residual
      ! Compute norm AFTER communication (SRC-compatible ordering)
      if(iter==1) then
         call cmp_residual_norm2_fine(ilevel, i_res_norm2)
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                 & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         i_res_norm2=i_res_norm2_tot
#endif
      end if

      ! Restrict residual into upper level RHS
#ifdef HYDRO_CUDA
      if(use_ri_gpu) then
         ! GPU restriction: d_mg_f1 → d_coarse_rhs_flat → host → active_mg
         call cuda_mg_restrict_execute_c(int(active(ilevel)%ngrid,c_int), &
              int(ngridmax,c_int), int(ncoarse,c_int))
         call cuda_mg_restrict_download_c(mg_ri_coarse_rhs, &
              int(mg_ri_total_coarse,c_int))
         call scatter_coarse_rhs_from_flat(ilevel)
      else
#endif
         ! First clear the rhs in coarser reception comms
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,2)=0.0d0
         end do
         ! Restrict
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call restrict_residual_fine_reverse(ilevel)
#ifdef FDMDEBUG
         call mgp_stop(MGP_RESTRICT,mgp_phase_wall,mgp_phase_cpu, &
              int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
#ifdef HYDRO_CUDA
      end if
#endif
      call make_reverse_mg_dp(2,ilevel-1) ! communicate rhs

      if(ilevel>1) then
         ! Reset correction at upper level before solve
         do icpu=1,ncpu
            if(active_mg(icpu,ilevel-1)%ngrid==0) cycle
            active_mg(icpu,ilevel-1)%u(:,1)=0.0d0
         end do

         ! Multigrid-solve the upper level
         call recursive_multigrid_coarse(ilevel-1, safe_mode(ilevel))

         ! Interpolate coarse solution and correct fine solution
#ifdef HYDRO_CUDA
         if(use_ri_gpu) then
            ! GPU interpolation: active_mg → flat → GPU → correct d_mg_phi
            call gather_coarse_phi_to_flat(ilevel)
            call cuda_mg_interp_upload_c(mg_ri_coarse_phi, &
                 int(mg_ri_total_coarse,c_int))
            call cuda_mg_interp_execute_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int))
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
#ifdef HYDRO_CUDA
            if(use_mg_gpu) then
               call cuda_mg_download_phi_c(phi, ncell_tot_c)
            end if
#endif
#ifdef FDMDEBUG
            call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
            call interpolate_and_correct_fine(ilevel)
#ifdef FDMDEBUG
            call mgp_stop(MGP_INTERP,mgp_phase_wall,mgp_phase_cpu, &
                 int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
            call make_virtual_fine_dp_mg_profile(phi(1),ilevel)
#ifdef HYDRO_CUDA
            if(use_mg_gpu) then
               call cuda_mg_upload_phi_c(phi, ncell_tot_c)
            end if
         end if
#endif
      end if

      ! Post-smoothing
      do i=1,ngs_fine
#ifdef HYDRO_CUDA
         if(use_mg_gpu) then
            safe_int = 0
            if(safe_mode(ilevel)) safe_int = 1
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 0, safe_int)
            if(.not.mg_merged_rb) call make_virtual_fine_dp_gpu(ilevel)
            call cuda_mg_gauss_seidel_c(int(active(ilevel)%ngrid,c_int), &
                 int(ngridmax,c_int), int(ncoarse,c_int), dx2_mg, 1, safe_int)
            call make_virtual_fine_dp_gpu(ilevel)
         else
#endif
#ifdef FDMDEBUG
            call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
            call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
#ifdef FDMDEBUG
            call mgp_stop(MGP_FINE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
                 int(active(ilevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
            if(.not.mg_merged_rb) call make_virtual_fine_dp_mg_profile(phi(1),ilevel) ! Communicate phi (Red)
#ifdef FDMDEBUG
            call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
            call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
#ifdef FDMDEBUG
            call mgp_stop(MGP_FINE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
                 int(active(ilevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
            call make_virtual_fine_dp_mg_profile(phi(1),ilevel) ! Communicate phi (Black)
#ifdef HYDRO_CUDA
         end if
#endif
      end do

      ! Update fine residual (fused with norm computation)
#ifdef HYDRO_CUDA
      if(use_mg_gpu) then
         gpu_norm2 = 0.0d0
         call cuda_mg_residual_c(int(active(ilevel)%ngrid,c_int), &
              int(ngridmax,c_int), int(ncoarse,c_int), &
              oneoverdx2_mg, dble(twondim), dx2_norm_mg, &
              gpu_norm2, 1)
         res_norm2 = gpu_norm2
      else
#endif
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call cmp_residual_mg_fine(ilevel)
#ifdef FDMDEBUG
         call mgp_stop(MGP_FINE_RESID,mgp_phase_wall,mgp_phase_cpu, &
              int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
#ifdef HYDRO_CUDA
      end if
#endif
      call make_virtual_fine_dp_mg_profile(f(1,1),ilevel) ! communicate residual
      ! Compute norm AFTER communication (SRC-compatible ordering)
      call cmp_residual_norm2_fine(ilevel, res_norm2)
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1, &
              & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      res_norm2=res_norm2_tot
#endif

      last_err = err
      err = sqrt(res_norm2/(i_res_norm2+1d-20*rho_tot**2))

      ! Verbosity
      if(verbose .or. use_mg_gpu) then
         if(myid==1) print '(A,I5,A,1pE10.3,A,I3,A,L1)', &
              '   ==> Step=',iter,' Error=',err, &
              ' level=',ilevel,' gpu=',use_mg_gpu
      end if

      ! Converged?
      ! Also exit when absolute residual is below the rho_tot floor
      ! (happens when fine-level density is injected from coarse and phi
      !  is already the solution — i_res_norm2 ≈ 0 inflates relative err)
      if(err<epsilon .or. res_norm2<1d-20*rho_tot**2 &
           .or. iter>=max(1,maxiter_fine)) exit

      ! Not converged, check error and possibly enable safe mode for the level
      if(err > last_err*SAFE_FACTOR .and. (.not. safe_mode(ilevel))) then
         if(verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
         safe_mode(ilevel) = .true.
      end if

   end do main_iteration_loop

   ! Download final phi from GPU to CPU before cleanup
#ifdef HYDRO_CUDA
   if(use_mg_gpu) then
      call cuda_mg_download_phi_c(phi, ncell_tot_c)
   end if
#endif

   ! Cleanup GPU MG state
#ifdef HYDRO_CUDA
   if(use_ri_gpu) then
      call cuda_mg_ri_free_c()
   end if
   if(allocated(mg_ri_flat_offset))  deallocate(mg_ri_flat_offset)
   if(allocated(mg_ri_coarse_rhs))   deallocate(mg_ri_coarse_rhs)
   if(allocated(mg_ri_coarse_phi))   deallocate(mg_ri_coarse_phi)
   mg_ri_total_coarse = 0
   if(use_mg_gpu) then
      call cuda_mg_halo_free_c()
      call cuda_mg_free_c()
   end if
   if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
   if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
   if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
   if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)
   mg_halo_n_emit = 0
   mg_halo_n_recv = 0
#endif

   ! Free pre-computed neighbor grids
   if(allocated(nbor_grid_fine)) deallocate(nbor_grid_fine)

   ! Free pre-computed coarse neighbor cache
   if(ilevel>1 .and. levelmin_mg < ilevel) then
      call cleanup_nbor_grid_coarse(levelmin_mg, ilevel-1)
   end if

   if(myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel, ' Step=', &
            iter,' Error=',err
   mg_failed=iter>=max(1,maxiter_fine) .and. err>=epsilon &
        .and. res_norm2>=1d-20*rho_tot**2
   if(myid==1 .and. mg_failed) &
      print *,'WARN: Fine multigrid Poisson failed to converge...'
   if(mg_failed .and. abort_on_mg_nonconvergence)then
      if(myid==1) print *,'FATAL: abort_on_mg_nonconvergence is enabled'
#ifndef WITHOUTMPI
      call MPI_ABORT(MPI_COMM_WORLD,914,info)
#else
      stop 914
#endif
   end if

   ! ---------------------------------------------------------------------
   ! Cleanup MG levels after solve complete
   ! ---------------------------------------------------------------------
   do ifine=1,ilevel-1
      call cleanup_mg_level(ifine)
   end do

#ifdef FDMDEBUG
   call mgp_stop(MGP_TOTAL,mgp_wall_start,mgp_cpu_start, &
        int(active(ilevel)%ngrid,kind=8)*int(twotondim,kind=8))
   call mgp_report(nstep_coarse,myid,ncpu,omp_get_max_threads())
#endif

end subroutine multigrid_fine


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Build flat halo cell index arrays for GPU MG phi exchange
! Enumerates emission and reception cell indices in the same order
! as make_virtual_fine_dp packing.
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine build_mg_halo_indices(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
  use amr_index, only: icell_of
   implicit none

   integer, intent(in) :: ilevel

   integer :: icpu, i, j, idx

   ! Count total emission and reception cells
   mg_halo_n_emit = 0
   mg_halo_n_recv = 0
   do icpu = 1, ncpu
      mg_halo_n_emit = mg_halo_n_emit + emission(icpu,ilevel)%ngrid * twotondim
      mg_halo_n_recv = mg_halo_n_recv + reception(icpu,ilevel)%ngrid * twotondim
   end do

   ! Allocate flat arrays
   if(allocated(mg_halo_emit_cells)) deallocate(mg_halo_emit_cells)
   if(allocated(mg_halo_recv_cells)) deallocate(mg_halo_recv_cells)
   if(allocated(mg_halo_emit_buf))   deallocate(mg_halo_emit_buf)
   if(allocated(mg_halo_recv_buf))   deallocate(mg_halo_recv_buf)

   allocate(mg_halo_emit_cells(1:max(mg_halo_n_emit,1)))
   allocate(mg_halo_recv_cells(1:max(mg_halo_n_recv,1)))
   allocate(mg_halo_emit_buf(1:max(mg_halo_n_emit,1)))
   allocate(mg_halo_recv_buf(1:max(mg_halo_n_recv,1)))

   ! Build emission cell indices
   ! Order: for each CPU, for each child cell j, for each grid i
   ! This matches the packing order in make_virtual_fine_dp
   idx = 0
   do icpu = 1, ncpu
      if(emission(icpu,ilevel)%ngrid > 0) then
         do j = 1, twotondim
            do i = 1, emission(icpu,ilevel)%ngrid
               idx = idx + 1
               mg_halo_emit_cells(idx) = icell_of(emission(icpu,ilevel)%igrid(i),j)
            end do
         end do
      end if
   end do

   ! Build reception cell indices (same order)
   idx = 0
   do icpu = 1, ncpu
      if(reception(icpu,ilevel)%ngrid > 0) then
         do j = 1, twotondim
            do i = 1, reception(icpu,ilevel)%ngrid
               idx = idx + 1
               mg_halo_recv_cells(idx) = icell_of(reception(icpu,ilevel)%igrid(i),j)
            end do
         end do
      end if
   end do

   if(myid==1) write(*,'(A,I3,A,I8,A,I8)') &
        ' MG halo indices: level=',ilevel, &
        ' n_emit=',mg_halo_n_emit,' n_recv=',mg_halo_n_recv

   ! Upload cell indices to GPU
   call cuda_mg_halo_setup_c( &
        mg_halo_emit_cells, int(mg_halo_n_emit, c_int), &
        mg_halo_recv_cells, int(mg_halo_n_recv, c_int))

end subroutine build_mg_halo_indices
#endif


subroutine make_virtual_fine_dp_mg_profile(xx,ilevel)
   use amr_commons
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
   implicit none

   integer, intent(in) :: ilevel
   real(dp), dimension(1:ncoarse+ngridmax*twotondim) :: xx
#ifdef FDMDEBUG
   real(kind=8) :: wall_start,cpu_start
   integer(kind=8) :: nwork_total
   integer :: icpu
#endif

   if(numbtot(1,ilevel)==0)return
#ifdef FDMDEBUG
   nwork_total=0_8
   do icpu=1,ncpu
      nwork_total=nwork_total+int(emission(icpu,ilevel)%ngrid,kind=8) &
           +int(reception(icpu,ilevel)%ngrid,kind=8)
   end do
   nwork_total=nwork_total*int(twotondim,kind=8)
   call mgp_start(wall_start,cpu_start)
#endif

   call make_virtual_fine_dp(xx,ilevel)
#ifdef FDMDEBUG
   call mgp_stop(MGP_FINE_COMM,wall_start,cpu_start,nwork_total)
#endif

end subroutine make_virtual_fine_dp_mg_profile

! ------------------------------------------------------------------------
! GPU phi exchange for MG smoothing
! Full D2H → MPI exchange → full H2D
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine make_virtual_fine_dp_gpu(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
   implicit none

   integer, intent(in) :: ilevel
   integer :: i

   if(mg_halo_n_emit == 0 .and. mg_halo_n_recv == 0) return

   ! Step 1: GPU gather — emission cell values from GPU phi → flat host buffer
   if(mg_halo_n_emit > 0) then
      call cuda_mg_halo_gather_c(mg_halo_emit_buf, int(mg_halo_n_emit, c_int))
      ! Step 2: Scatter flat buffer → host phi at emission positions
      do i = 1, mg_halo_n_emit
         phi(mg_halo_emit_cells(i)) = mg_halo_emit_buf(i)
      end do
   end if

   ! Step 3: MPI exchange (packs phi at emission cells, unpacks to reception cells)
   call make_virtual_fine_dp_mg_profile(phi(1), ilevel)

   ! Step 4: Gather host phi at reception positions → flat buffer
   if(mg_halo_n_recv > 0) then
      do i = 1, mg_halo_n_recv
         mg_halo_recv_buf(i) = phi(mg_halo_recv_cells(i))
      end do
      ! Step 5: GPU scatter — flat host buffer → GPU phi at reception positions
      call cuda_mg_halo_scatter_c(mg_halo_recv_buf, int(mg_halo_n_recv, c_int))
   end if

end subroutine make_virtual_fine_dp_gpu
#endif


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Precompute GPU restrict/interp mapping arrays for V-cycle
! Builds restrict_target and interp_nbor_flat, then uploads to GPU.
! Also allocates flat host buffers for coarse data transfer.
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine precompute_mg_gpu_restrict_interp(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_cuda_interface
   use iso_c_binding
  use amr_index, only: ichild_of, igrid_of
   implicit none

   integer, intent(in) :: ilevel

   integer :: icoarselevel, ngrid_fine
   integer :: igrid_f_mg, igrid_f_amr, icell_c_amr, ind_c_cell
   integer :: igrid_c_amr, igrid_c_mg, cpu_amr, icell_c_mg
   integer :: i, j, istart, nbatch, icpu, ngrid_c

   integer, allocatable :: h_restrict_target(:)
   integer, allocatable :: h_interp_nbor_flat(:,:)

   real(dp) :: a, b, c, d
   real(dp) :: bbb(8)
   integer :: ccc(8,8), ccc_gpu(8,8)

   ! Work arrays for get3cubefather
   integer :: ind_cell_father_loc(1:nvector)
   integer :: nbors_father_cells_loc(1:nvector, 1:threetondim)
   integer :: nbors_father_grids_loc(1:nvector, 1:twotondim)

   icoarselevel = ilevel - 1
   ngrid_fine = active(ilevel)%ngrid

   if(ngrid_fine == 0 .or. icoarselevel < 1) return

   ! Compute flat_offset for all CPUs
   if(allocated(mg_ri_flat_offset)) deallocate(mg_ri_flat_offset)
   allocate(mg_ri_flat_offset(1:ncpu))
   mg_ri_flat_offset(1) = 0
   do icpu = 2, ncpu
      mg_ri_flat_offset(icpu) = mg_ri_flat_offset(icpu-1) + &
           active_mg(icpu-1,icoarselevel)%ngrid * twotondim
   end do
   mg_ri_total_coarse = mg_ri_flat_offset(ncpu) + &
        active_mg(ncpu,icoarselevel)%ngrid * twotondim

   if(mg_ri_total_coarse == 0) return

   ! Allocate host flat buffers for coarse data
   if(allocated(mg_ri_coarse_rhs)) deallocate(mg_ri_coarse_rhs)
   if(allocated(mg_ri_coarse_phi)) deallocate(mg_ri_coarse_phi)
   allocate(mg_ri_coarse_rhs(1:mg_ri_total_coarse))
   allocate(mg_ri_coarse_phi(1:mg_ri_total_coarse))

   ! Set interpolation coefficients (from interpolate_and_correct_fine)
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.0D0*a
   bbb(:) = (/a, b, b, c, b, c, c, d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Transpose ccc for GPU C row-major layout:
   ! Fortran ccc_gpu(f,a) stored at offset (a-1)*8+(f-1) matches C d_ccc[a-1][f-1]
   do i = 1, 8
      do j = 1, 8
         ccc_gpu(j, i) = ccc(i, j)
      end do
   end do

   ! === Compute restrict_target ===
   allocate(h_restrict_target(1:ngrid_fine))

   do igrid_f_mg = 1, ngrid_fine
      igrid_f_amr = active(ilevel)%igrid(igrid_f_mg)
      icell_c_amr = father(igrid_f_amr)
      ind_c_cell  = ichild_of(icell_c_amr)
      igrid_c_amr = igrid_of(icell_c_amr)
      cpu_amr     = cpu_map(father(igrid_c_amr))
      igrid_c_mg  = lookup_mg(igrid_c_amr)

      if(igrid_c_mg <= 0) then
         h_restrict_target(igrid_f_mg) = -1
      else
         ngrid_c = active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg = (ind_c_cell - 1) * ngrid_c + igrid_c_mg
         ! Check coarse mask
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4) <= 0d0) then
            h_restrict_target(igrid_f_mg) = -1
         else
            h_restrict_target(igrid_f_mg) = mg_ri_flat_offset(cpu_amr) + &
                 icell_c_mg - 1  ! 0-based for C
         end if
      end if
   end do

   ! === Compute interp_nbor_flat ===
   allocate(h_interp_nbor_flat(1:27, 1:ngrid_fine))

   do istart = 1, ngrid_fine, nvector
      nbatch = min(nvector, ngrid_fine - istart + 1)

      ! Get father cells
      do i = 1, nbatch
         ind_cell_father_loc(i) = father(active(ilevel)%igrid(istart+i-1))
      end do

      ! Get 3x3x3 neighbor cells
      call get3cubefather(ind_cell_father_loc, nbors_father_cells_loc, &
           nbors_father_grids_loc, nbatch, ilevel)

      ! Convert to flat indices
      do j = 1, threetondim  ! 27
         do i = 1, nbatch
            igrid_f_mg  = istart + i - 1
            icell_c_amr = nbors_father_cells_loc(i, j)
            ind_c_cell  = ichild_of(icell_c_amr)
            igrid_c_amr = igrid_of(icell_c_amr)
            cpu_amr     = cpu_map(father(igrid_c_amr))
            igrid_c_mg  = lookup_mg(igrid_c_amr)

            if(igrid_c_mg <= 0) then
               h_interp_nbor_flat(j, igrid_f_mg) = -1
            else
               ngrid_c = active_mg(cpu_amr,icoarselevel)%ngrid
               icell_c_mg = (ind_c_cell - 1) * ngrid_c + igrid_c_mg
               h_interp_nbor_flat(j, igrid_f_mg) = &
                    mg_ri_flat_offset(cpu_amr) + icell_c_mg - 1  ! 0-based
            end if
         end do
      end do
   end do

   ! Upload to GPU
   call cuda_mg_ri_setup_c(h_restrict_target, h_interp_nbor_flat, &
        int(ngrid_fine, c_int), int(mg_ri_total_coarse, c_int), &
        bbb, ccc_gpu)

   deallocate(h_restrict_target, h_interp_nbor_flat)

   if(myid==1) write(*,'(A,I3,A,I10,A,I10)') &
        ' MG RI precompute: level=',ilevel, &
        ' ngrid_fine=',ngrid_fine,' total_coarse=',mg_ri_total_coarse

end subroutine precompute_mg_gpu_restrict_interp
#endif


! ------------------------------------------------------------------------
! Scatter coarse RHS from flat array into active_mg structure
! Called after GPU restrict download
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine scatter_coarse_rhs_from_flat(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel
   integer :: icoarselevel, icpu, ngrid_c, ncells_c, offset

   icoarselevel = ilevel - 1

   ! Scatter flat data into active_mg%u(:,2)
   ! The flat buffer replaces (not accumulates) since it was zeroed on GPU
   do icpu = 1, ncpu
      ngrid_c = active_mg(icpu,icoarselevel)%ngrid
      if(ngrid_c == 0) cycle
      ncells_c = ngrid_c * twotondim
      offset = mg_ri_flat_offset(icpu)
      active_mg(icpu,icoarselevel)%u(1:ncells_c,2) = &
           mg_ri_coarse_rhs(offset+1:offset+ncells_c)
   end do

end subroutine scatter_coarse_rhs_from_flat
#endif


! ------------------------------------------------------------------------
! Gather coarse phi from active_mg structure into flat array
! Called before GPU interp upload
! ------------------------------------------------------------------------
#ifdef HYDRO_CUDA
subroutine gather_coarse_phi_to_flat(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel
   integer :: icoarselevel, icpu, ngrid_c, ncells_c, offset

   icoarselevel = ilevel - 1

   ! Gather active_mg%u(:,1) into flat array
   do icpu = 1, ncpu
      ngrid_c = active_mg(icpu,icoarselevel)%ngrid
      if(ngrid_c == 0) cycle
      ncells_c = ngrid_c * twotondim
      offset = mg_ri_flat_offset(icpu)
      mg_ri_coarse_phi(offset+1:offset+ncells_c) = &
           active_mg(icpu,icoarselevel)%u(1:ncells_c,1)
   end do

end subroutine gather_coarse_phi_to_flat
#endif


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------
recursive subroutine recursive_multigrid_coarse(ifinelevel, safe)
   use amr_commons
   use poisson_commons
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ifinelevel
   logical, intent(in) :: safe

   real(dp) :: debug_norm2, debug_norm2_tot
   integer :: i, icpu, info, icycle, ncycle
#ifdef FDMDEBUG
   real(kind=8) :: mgp_phase_wall,mgp_phase_cpu
#endif

   if(ifinelevel<=levelmin_mg) then
      ! Solve 'directly'
      do i=1,2*ngs_coarse
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Red)
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Black)
      end do
      return
   end if

   if(safe) then
      ncycle=ncycles_coarse_safe
   else
      ncycle=1
   endif

   do icycle=1,ncycle

      ! Pre-smoothing
      do i=1,ngs_coarse
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Red)
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Black)
      end do

      ! Compute residual and restrict into upper level RHS
#ifdef FDMDEBUG
      call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
      call cmp_residual_mg_coarse(ifinelevel)
#ifdef FDMDEBUG
      call mgp_stop(MGP_COARSE_RESID,mgp_phase_wall,mgp_phase_cpu, &
           int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
      call make_virtual_mg_dp(3,ifinelevel)  ! Communicate residual

      ! First clear the rhs in coarser reception comms
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,2)=0.0d0
      end do

      ! Restrict and do reverse-comm
#ifdef FDMDEBUG
      call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
      call restrict_residual_coarse_reverse(ifinelevel)
#ifdef FDMDEBUG
      call mgp_stop(MGP_RESTRICT,mgp_phase_wall,mgp_phase_cpu, &
           int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
      call make_reverse_mg_dp(2,ifinelevel-1) ! communicate rhs

      ! Reset correction from upper level before solve
      do icpu=1,ncpu
         if(active_mg(icpu,ifinelevel-1)%ngrid==0) cycle
         active_mg(icpu,ifinelevel-1)%u(:,1)=0.0d0
      end do

      ! Multigrid-solve the upper level
      call recursive_multigrid_coarse(ifinelevel-1, safe)

      ! Interpolate coarse solution and correct back into fine solution
#ifdef FDMDEBUG
      call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
      call interpolate_and_correct_coarse(ifinelevel)
#ifdef FDMDEBUG
      call mgp_stop(MGP_INTERP,mgp_phase_wall,mgp_phase_cpu, &
           int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim,kind=8))
#endif
      call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution

      ! Post-smoothing
      do i=1,ngs_coarse
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Red)
#ifdef FDMDEBUG
         call mgp_start(mgp_phase_wall,mgp_phase_cpu)
#endif
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
#ifdef FDMDEBUG
         call mgp_stop(MGP_COARSE_SMOOTH,mgp_phase_wall,mgp_phase_cpu, &
              int(active_mg(myid,ifinelevel)%ngrid,kind=8)*int(twotondim/2,kind=8))
#endif
         call make_virtual_mg_dp(1,ifinelevel)  ! Communicate solution (Black)
      end do

   end do

end subroutine recursive_multigrid_coarse

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid communicator building
! ------------------------------------------------------------------------
subroutine build_parent_comms_mg(active_f_comm, ifinelevel)
   use amr_commons
   use poisson_commons
   use ksection
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
   integer, dimension (MPI_STATUS_SIZE, ncpu) :: statuses
   integer :: ntotal_ksec, nrecv_ksec, idx_ksec
   real(dp), allocatable :: sbuf_ksec(:,:), rbuf_ksec(:,:)
   integer, allocatable :: dcpu_ksec(:)
#endif

   integer, intent(in) :: ifinelevel
   type(communicator), intent(in) :: active_f_comm

   integer :: icoarselevel
   integer :: ngrids, cur_grid, cur_cpu, cur_cell, newgrids
   integer :: i, nbatch, ind, icpu, istart, info

   integer :: nact_tot, nreq_tot, nreq_tot2
   integer, dimension(1:ncpu) :: nreq, nreq2

   ! Per-thread work arrays for OpenMP (Stages 1 & 4)
   integer, dimension(:,:), target, allocatable :: P_icf_bp
   integer, dimension(:,:,:), target, allocatable :: P_nfg_bp, P_nfc_bp
   integer, dimension(:), pointer :: ind_cell_father
   integer, dimension(:,:), pointer :: nbors_father_grids
   integer, dimension(:,:), pointer :: nbors_father_cells
   common /omp_build_parent_comms/ ind_cell_father, nbors_father_cells, nbors_father_grids
!$omp threadprivate(/omp_build_parent_comms/)
   integer :: mythread, nthreads
   common /openmpthreads/ mythread, nthreads
!$omp threadprivate(/openmpthreads/)

   type(communicator), dimension(1:ncpu) :: comm_send, comm_receive
   type(communicator), dimension(1:ncpu) :: comm_send2, comm_receive2

   integer, dimension(1:ncpu) :: indx, recvbuf, recvbuf2
   integer, dimension(1:ncpu) :: reqsend, reqrecv
   integer :: countrecv, countsend
   integer :: tag = 777
#ifdef FDMDEBUG
   real(kind=8) :: mgp_wall_start,mgp_cpu_start
#endif


   icoarselevel=ifinelevel-1
#ifdef FDMDEBUG
   call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

   nact_tot=0
   nreq_tot=0; nreq=0
   indx=0; recvbuf=0

   ! Setup per-thread work arrays for OpenMP
!$omp parallel
   mythread = omp_get_thread_num()
   nthreads = omp_get_num_threads()
!$omp end parallel
   allocate(P_icf_bp(1:nvector, 0:nthreads-1))
   allocate(P_nfg_bp(1:nvector, 1:twotondim, 0:nthreads-1))
   allocate(P_nfc_bp(1:nvector, 1:threetondim, 0:nthreads-1))
!$omp parallel
   ind_cell_father => P_icf_bp(:, mythread)
   nbors_father_grids => P_nfg_bp(:, :, mythread)
   nbors_father_cells => P_nfc_bp(:, :, mythread)
!$omp end parallel

   ! ---------------------------------------------------------------------
   ! STAGE 1 : Coarse grid MG activation for local grids (1st pass)
   ! ---------------------------------------------------------------------

   ! Loop over the AMR active communicator first
   ngrids = active_f_comm%ngrid
!$omp parallel do private(istart,nbatch,i,ind,cur_grid,cur_cpu) schedule(dynamic,128)
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_f_comm%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells and grids
      call get3cubefather(ind_cell_father,nbors_father_cells, &
         & nbors_father_grids,nbatch,ifinelevel)

      ! Now process the twotondim father grids
      do ind=1,twotondim
         do i=1,nbatch
            cur_grid = nbors_father_grids(i,ind)
            if(lookup_mg(cur_grid)>0) cycle ! Grid already active (pre-check)

            cur_cpu=cpu_map(father(cur_grid))
            if(cur_cpu==0) cycle

            !$omp critical(stage1_update)
            if(lookup_mg(cur_grid)<=0) then  ! Definitive check under lock
               if(cur_cpu==myid) then
                  ! Stack grid for local activation
                  nact_tot=nact_tot+1
                  flag2(nact_tot)=cur_grid
                  lookup_mg(cur_grid)=nact_tot
               else
                  ! Stack grid for remote activation
                  nreq_tot=nreq_tot+1
                  nreq(cur_cpu)=nreq(cur_cpu)+1
                  flag2(ngridmax+nreq_tot)=cur_grid
                  lookup_mg(cur_grid)=abs(lookup_mg(cur_grid))
               end if
            end if
            !$omp end critical(stage1_update)
         end do
      end do
   end do


   ! ---------------------------------------------------------------------
   ! STAGE 2 : Coarse grid MG activation request
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   if(ordering=='ksection') then
      ntotal_ksec = 0
      do icpu = 1, ncpu
         if(nreq(icpu) > 0) ntotal_ksec = ntotal_ksec + 1
      end do
      allocate(sbuf_ksec(1:2, 1:max(ntotal_ksec,1)))
      allocate(dcpu_ksec(1:max(ntotal_ksec,1)))
      idx_ksec = 0
      do icpu = 1, ncpu
         if(nreq(icpu) > 0) then
            idx_ksec = idx_ksec + 1
            dcpu_ksec(idx_ksec) = icpu
            sbuf_ksec(1, idx_ksec) = dble(myid)
            sbuf_ksec(2, idx_ksec) = dble(nreq(icpu))
         end if
      end do
      call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
      recvbuf = 0
      do idx_ksec = 1, nrecv_ksec
         recvbuf(nint(rbuf_ksec(1, idx_ksec))) = nint(rbuf_ksec(2, idx_ksec))
      end do
      deallocate(sbuf_ksec, dcpu_ksec, rbuf_ksec)
   else
      call MPI_ALLTOALL(nreq, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, &
         & MPI_COMM_WORLD, info)
   end if

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive(icpu)%ngrid=recvbuf(icpu)
      if(recvbuf(icpu)>0) allocate(comm_receive(icpu)%igrid(1:recvbuf(icpu)))
   end do

   ! Receive to-be-activated grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (activation request) communicators
   do icpu=1,ncpu
      comm_send(icpu)%ngrid=nreq(icpu)
      if(nreq(icpu)>0) allocate(comm_send(icpu)%igrid(1:nreq(icpu)))
   end do
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      comm_send(cur_cpu)%igrid(nreq(cur_cpu))=lookup_mg(cur_grid) ! Remote
   end do

   ! Send to-be-activated grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Activate requested grids
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         do i=1,ngrids
            cur_grid=comm_receive(icpu)%igrid(i) ! Local AMR index
            if(lookup_mg(cur_grid)>0) cycle      ! Already active: cycle
            ! Activate grid
            nact_tot=nact_tot+1
            flag2(nact_tot)=cur_grid
            lookup_mg(cur_grid)=nact_tot
         end do
      end if
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)
#endif

   ! ---------------------------------------------------------------------
   ! STAGE 3 : Coarse grid MG active comm gathering
   ! ---------------------------------------------------------------------

   active_mg(myid,icoarselevel)%ngrid=nact_tot
   if(nact_tot>0) then
      allocate( active_mg(myid,icoarselevel)%igrid(1:nact_tot) )
      allocate( active_mg(myid,icoarselevel)%u(1:nact_tot*twotondim,1:4) )
      allocate( active_mg(myid,icoarselevel)%f(1:nact_tot*twotondim,1:1) )
      active_mg(myid,icoarselevel)%igrid=0
      active_mg(myid,icoarselevel)%u=0.0d0
      active_mg(myid,icoarselevel)%f=0
   end if
   do i=1,nact_tot
      active_mg(myid,icoarselevel)%igrid(i)=flag2(i)
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 4 : Screen active grid neighbors for new reception grids
   ! ---------------------------------------------------------------------
   ngrids = active_mg(myid,icoarselevel)%ngrid
   nreq2 = 0
   nreq_tot2 = 0
!$omp parallel do private(istart,nbatch,i,ind,cur_cell,cur_cpu,cur_grid) schedule(dynamic,128)
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_mg(myid,icoarselevel)%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells
      call get3cubefather(ind_cell_father,nbors_father_cells,nbors_father_grids,nbatch,icoarselevel)

      ! Now process the father grids
      do ind=1,threetondim
         do i=1,nbatch
            cur_cell = nbors_father_cells(i,ind)
            cur_cpu  = cpu_map(cur_cell)
            if(cur_cpu==0) cycle
            cur_grid = son(cur_cell)
            if(cur_cpu/=myid) then
               ! Neighbor cell is not managed by current CPU
               if (cur_grid==0) cycle              ! No grid there
               if (lookup_mg(cur_grid)>0) cycle    ! Already selected (pre-check)
               !$omp critical(stage4_update)
               if(lookup_mg(cur_grid)<=0) then  ! Definitive check under lock
                  nreq_tot2=nreq_tot2+1
                  nreq2(cur_cpu)=nreq2(cur_cpu)+1
                  flag2(ngridmax+nreq_tot+nreq_tot2)=cur_grid
                  lookup_mg(cur_grid)=abs(lookup_mg(cur_grid))
               end if
               !$omp end critical(stage4_update)
            end if
         end do
      end do
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 5 : Share new reception grid requests, build emission comms
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   recvbuf2=0
   if(ordering=='ksection') then
      ntotal_ksec = 0
      do icpu = 1, ncpu
         if(nreq2(icpu) > 0) ntotal_ksec = ntotal_ksec + 1
      end do
      allocate(sbuf_ksec(1:2, 1:max(ntotal_ksec,1)))
      allocate(dcpu_ksec(1:max(ntotal_ksec,1)))
      idx_ksec = 0
      do icpu = 1, ncpu
         if(nreq2(icpu) > 0) then
            idx_ksec = idx_ksec + 1
            dcpu_ksec(idx_ksec) = icpu
            sbuf_ksec(1, idx_ksec) = dble(myid)
            sbuf_ksec(2, idx_ksec) = dble(nreq2(icpu))
         end if
      end do
      call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
      do idx_ksec = 1, nrecv_ksec
         recvbuf2(nint(rbuf_ksec(1, idx_ksec))) = nint(rbuf_ksec(2, idx_ksec))
      end do
      deallocate(sbuf_ksec, dcpu_ksec, rbuf_ksec)
   else
      call MPI_ALLTOALL(nreq2, 1, MPI_INTEGER, recvbuf2, 1, MPI_INTEGER, &
         & MPI_COMM_WORLD, info)
   end if

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive2(icpu)%ngrid=recvbuf2(icpu)
      if(recvbuf2(icpu)>0) allocate(comm_receive2(icpu)%igrid(1:recvbuf2(icpu)))
   end do

   ! Receive potential reception grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive2(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive2(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (reception request) communicators
   do icpu=1,ncpu
      comm_send2(icpu)%ngrid=nreq2(icpu)
      if(nreq2(icpu)>0) allocate(comm_send2(icpu)%igrid(1:nreq2(icpu)))
   end do
   nreq2=0
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))=lookup_mg(cur_grid) ! Remote AMR index
      ! Restore negative lookup_mg
      lookup_mg(cur_grid)=-abs(lookup_mg(cur_grid))
   end do

   ! Send reception request grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send2(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Compute local MG indices of inbound grids, alloc & fill emission comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      newgrids=0
      do i=1,recvbuf2(icpu)
         ! MAP AMR -> MG INDICES IN PLACE
         comm_receive2(icpu)%igrid(i)=lookup_mg(comm_receive2(icpu)%igrid(i))
         if(comm_receive2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate emission communicators
      ngrids=recvbuf(icpu)+newgrids
      emission_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(emission_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(emission_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4) )
         allocate(emission_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         emission_mg(icpu,icoarselevel)%igrid=0
         emission_mg(icpu,icoarselevel)%u=0.0d0
         emission_mg(icpu,icoarselevel)%f=0
      end if
      ! First part: activation request emission grids
      do i=1,recvbuf(icpu)
         emission_mg(icpu,icoarselevel)%igrid(i)=lookup_mg(comm_receive(icpu)%igrid(i))
      end do
      ! Second part: new emission grids
      cur_grid=recvbuf(icpu)
      do i=1,recvbuf2(icpu)
         if(comm_receive2(icpu)%igrid(i)>0) then
            cur_grid=cur_grid+1
            emission_mg(icpu,icoarselevel)%igrid(cur_grid)=comm_receive2(icpu)%igrid(i)
         end if
      end do
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! ---------------------------------------------------------------------
   ! STAGE 6 : Reply with local MG grid status and build reception comms
   ! ---------------------------------------------------------------------
   ! Receive MG mappings from other CPUs back into comm_send2
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=nreq2(icpu)
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Send local MG mappings to other CPUs from comm_receive
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=recvbuf2(icpu)
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_receive2(icpu)%igrid,ngrids,MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for full completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Count remotely active MG grids, and allocate and fill reception comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ! Count requested grids which are MG-active remotely
      newgrids=0
      do i=1,nreq2(icpu)
         if(comm_send2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate and fill reception communicators on the fly
      ngrids=nreq(icpu)+newgrids
      active_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(active_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(active_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4))
         allocate(active_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         active_mg(icpu,icoarselevel)%igrid=0
         active_mg(icpu,icoarselevel)%u=0.0d0
         active_mg(icpu,icoarselevel)%f=0
      end if
   end do
   
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      ! Add to reception comm
      active_mg(cur_cpu,icoarselevel)%igrid(nreq(cur_cpu))=cur_grid
      ! Backup lookup_mg into flag2
      flag2(cur_grid)=lookup_mg(cur_grid)
      ! Update lookup_mg
      lookup_mg(cur_grid)=nreq(cur_cpu)
   end do

   nreq2=0; indx=nreq
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      if(comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))>0) then
         indx(cur_cpu)=indx(cur_cpu)+1
         ! Add to reception comm
         active_mg(cur_cpu,icoarselevel)%igrid(indx(cur_cpu))=cur_grid
         ! Backup lookup_mg
         flag2(cur_grid)=-lookup_mg(cur_grid)
         ! Update lookup_mg
         lookup_mg(cur_grid)=indx(cur_cpu)
      end if
   end do

   ! Wait for full completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! Cleanup
   do icpu=1,ncpu
      if(comm_send (icpu)%ngrid>0) deallocate(comm_send (icpu)%igrid)
      if(comm_send2(icpu)%ngrid>0) deallocate(comm_send2(icpu)%igrid)
      if(comm_receive (icpu)%ngrid>0) deallocate(comm_receive (icpu)%igrid)
      if(comm_receive2(icpu)%ngrid>0) deallocate(comm_receive2(icpu)%igrid)
   end do
#endif

   deallocate(P_icf_bp, P_nfg_bp, P_nfc_bp)
#ifdef FDMDEBUG
   call mgp_stop(MGP_BUILD,mgp_wall_start,mgp_cpu_start, &
        int(active_f_comm%ngrid,kind=8))
#endif

end subroutine build_parent_comms_mg


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid level cleanup
! ------------------------------------------------------------------------
subroutine cleanup_mg_level(ilevel)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: igrid, icpu, cur_grid, cur_cpu

   ! ---------------------------------------------------------------------
   ! Cleanup lookup table
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      do igrid=1,active_mg(icpu,ilevel)%ngrid
         cur_grid=active_mg(icpu,ilevel)%igrid(igrid)
         cur_cpu=cpu_map(father(cur_grid))
         if(cur_cpu==myid) then
            lookup_mg(cur_grid)=0
         else
            lookup_mg(cur_grid)=-mod(flag2(cur_grid),ngridmax)
         end if
      end do
   end do

   ! ---------------------------------------------------------------------
   ! Deallocate communicators
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0)then
         deallocate(active_mg(icpu,ilevel)%igrid)
         deallocate(active_mg(icpu,ilevel)%u)
         deallocate(active_mg(icpu,ilevel)%f)
      endif
      active_mg(icpu,ilevel)%ngrid=0
      if(emission_mg(icpu,ilevel)%ngrid>0)then
         deallocate(emission_mg(icpu,ilevel)%igrid)
         deallocate(emission_mg(icpu,ilevel)%u)
         deallocate(emission_mg(icpu,ilevel)%f)
      endif
      emission_mg(icpu,ilevel)%ngrid=0
   end do

end subroutine cleanup_mg_level

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------
subroutine make_fine_mask(ilevel)

   use amr_commons
   use pm_commons
   use poisson_commons
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
  use amr_index, only: icell_of
   implicit none
   integer, intent(in) :: ilevel

   integer  :: ngrid
   integer  :: ind, igrid_mg, icpu, ibound
   integer  :: igrid_amr, icell_amr
#ifdef FDMDEBUG
   real(kind=8) :: mgp_wall_start,mgp_cpu_start
   integer(kind=8) :: mgp_nwork
#endif

   ngrid=active(ilevel)%ngrid
#ifdef FDMDEBUG
   mgp_nwork=int(ngrid,kind=8)*int(twotondim,kind=8)
   call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif
!$omp parallel do if(ngrid*twotondim>=4096) collapse(2) &
!$omp private(igrid_amr,icell_amr) schedule(static)
   do ind=1,twotondim
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = icell_of(igrid_amr,ind)
         ! Init mask to 1.0 on active cells :
         f(icell_amr,3) = 1.0d0
      end do
   end do
!$omp end parallel do

   do icpu=1,ncpu
      ngrid=reception(icpu,ilevel)%ngrid
      do ind=1,twotondim
         do igrid_mg=1,ngrid
            igrid_amr = reception(icpu,ilevel)%igrid(igrid_mg)
            icell_amr = icell_of(igrid_amr,ind)
            ! Init mask to 1.0 on virtual cells :
            f(icell_amr,3) = 1.0d0
         end do
      end do
   end do

   do ibound=1,nboundary
      ngrid=boundary(ibound,ilevel)%ngrid
      do ind=1,twotondim
         do igrid_mg=1,ngrid
            igrid_amr = boundary(ibound,ilevel)%igrid(igrid_mg)
            icell_amr = icell_of(igrid_amr,ind)
            ! Init mask to -1.0 on boundary cells :
            f(icell_amr,3) = -1.0d0
         end do
      end do
   end do
#ifdef FDMDEBUG
   call mgp_stop(MGP_MASK,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif

end subroutine make_fine_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Coarse-level scan flag setting
! ------------------------------------------------------------------------
subroutine set_scan_flag_coarse_omp(ilevel)
   use amr_commons
   use poisson_commons
   use morton_hash
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
   implicit none

   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr, cpu_nbor_amr
   integer :: icell_nbor_amr
   integer :: iskip_mg, icell_mg, igrid_nbor_mg, icell_nbor_mg
   integer, dimension(1:3,1:2,1:8) :: iii, jjj
#ifdef FDMDEBUG
   real(kind=8) :: mgp_wall_start,mgp_cpu_start
#endif

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active_mg(myid,ilevel)%ngrid
   if(ngrid==0)return
#ifdef FDMDEBUG
   call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif
   ! Avoid fork/join overhead on the tiny bottom levels where one serial
   ! vector sweep is cheaper than starting an OMP team.
   if(ngrid*twotondim<4096)then
      call set_scan_flag_coarse(ilevel)
#ifdef FDMDEBUG
      call mgp_stop(MGP_SCAN,mgp_wall_start,mgp_cpu_start, &
           int(ngrid,kind=8)*int(twotondim,kind=8))
#endif
      return
   end if

   ! Every (cell-within-grid, MG-grid) pair owns a distinct flag entry.
   ! Morton/hash and communicator data are read-only during this phase.
!$omp parallel do collapse(2) default(shared) &
!$omp private(iskip_mg,icell_mg,igrid_amr,scan_flag,inbor,idim,igshift, &
!$omp         igrid_nbor_amr,cpu_nbor_amr,icell_nbor_amr, &
!$omp         igrid_nbor_mg,icell_nbor_mg) schedule(static)
   do ind=1,twotondim
      do igrid_mg=1,ngrid
         iskip_mg=(ind-1)*ngrid
         igrid_amr=active_mg(myid,ilevel)%igrid(igrid_mg)
         icell_mg=iskip_mg+igrid_mg

         if(active_mg(myid,ilevel)%u(icell_mg,4)==1d0)then
            scan_flag=0
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift=iii(idim,inbor,ind)
                  if(igshift==0)then
                     igrid_nbor_amr=igrid_amr
                     cpu_nbor_amr=myid
                  else
                     igrid_nbor_amr=morton_nbor_grid(igrid_amr,ilevel,igshift)
                     icell_nbor_amr=morton_nbor_cell(igrid_amr,ilevel,igshift)
                     cpu_nbor_amr=cpu_map(icell_nbor_amr)
                  end if

                  if(igrid_nbor_amr==0)then
                     scan_flag=1
                     exit scan_flag_loop
                  end if
                  igrid_nbor_mg=lookup_mg(igrid_nbor_amr)
                  if(igrid_nbor_mg<=0)then
                     scan_flag=1
                     exit scan_flag_loop
                  end if
                  icell_nbor_mg=igrid_nbor_mg+(jjj(idim,inbor,ind)-1)* &
                       active_mg(cpu_nbor_amr,ilevel)%ngrid
                  if(active_mg(cpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0d0)then
                     scan_flag=1
                     exit scan_flag_loop
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         active_mg(myid,ilevel)%f(icell_mg,1)=scan_flag
      end do
   end do
!$omp end parallel do
#ifdef FDMDEBUG
   call mgp_stop(MGP_SCAN,mgp_wall_start,mgp_cpu_start, &
        int(ngrid,kind=8)*int(twotondim,kind=8))
#endif

end subroutine set_scan_flag_coarse_omp

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------
subroutine make_fine_bc_rhs(ilevel,icount)

   use amr_commons
   use pm_commons
   use poisson_commons
   use morton_hash
   use dark_energy_commons, only: de_table_loaded, get_de_ratio, f_de_val
   use scalar_de_commons, only: sde_dmcorr_of_a, horndeski_mu_of_a
#ifdef FDMDEBUG
   use mg_omp_profile_m
#endif
  use amr_index, only: icell_of
   implicit none
   integer, intent(in) :: ilevel,icount

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   real(dp) :: dx, oneoverdx2

   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr
   integer  :: ifathercell_nbor_amr

   ! Thread-private variables for OpenMP
   real(dp) :: phi_b, nb_mask, nb_phi, w
   real(dp), dimension(1:nvector,1:twotondim) :: phi_int
   integer,  dimension(1:nvector) :: ind_cell

   integer  :: nx_loc
   real(dp) :: scale, fourpi
   real(dp) :: omega_de_a_mg, omega_cb_mg
#ifdef FDMDEBUG
   real(kind=8) :: mgp_wall_start,mgp_cpu_start

   call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

   ! Set constants
   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) then
      if(use_neutrino .and. omega_nu > 0.0d0) then
         fourpi = 1.5D0*(omega_m - omega_nu)*aexp*scale
      else
         fourpi = 1.5D0*omega_m*aexp*scale
      end if
      ! DE perturbation boost (table-based, cs2_de~0 only)
      ! For cs2_de~0, R_DE is k-independent so MG boost is valid
      if(de_perturb .and. de_table_loaded .and. cs2_de < 1.0d-4) then
         omega_de_a_mg = omega_l * f_de_val(aexp) * aexp**3
         if(use_neutrino .and. omega_nu > 0.0d0) then
            omega_cb_mg = omega_m - omega_nu
         else
            omega_cb_mg = omega_m
         end if
         fourpi = fourpi * (1.0d0 + (omega_de_a_mg / omega_cb_mg) * get_de_ratio(1.0d-3, aexp))
      end if
      ! Coupled quintessence: DM mass evolution rho_dm*a^3/rho_dm0
      if(use_coupled_de .and. cde_vary_mass .and. use_quintessence) then
         fourpi = fourpi * sde_dmcorr_of_a(aexp)
      end if
      ! Horndeski quasi-static mu(a) (scale-independent limit)
      if(use_horndeski) fourpi = fourpi * horndeski_mu_of_a(aexp)
   end if

   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
   do ind=1,twotondim

      ! Loop over active grids — OpenMP parallelized
      ! Each igrid_mg writes to unique icell_amr, so no race condition
!$omp parallel do private(igrid_mg,igrid_amr,icell_amr,idim,inbor, &
!$omp    igshift,igrid_nbor_amr,icell_nbor_amr,ifathercell_nbor_amr, &
!$omp    nb_mask,nb_phi,w,phi_b,phi_int,ind_cell) &
!$omp schedule(dynamic,1024)
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = icell_of(igrid_amr,ind)

         ! Init BC-modified RHS to rho - rho_tot :
         f(icell_amr,2) = fourpi*(rho(icell_amr) - rho_tot)

         if(f(icell_amr,3)<=0.0) cycle ! Do not process masked cells

         ! Separate directions
         do idim=1,ndim
            ! Loop over the 2 neighbors
            do inbor=1,2
               ! Get neighbor grid shift
               igshift = iii(idim,inbor,ind)

               ! Get neighbor grid using precomputed array
               if(igshift==0) then
                  igrid_nbor_amr = igrid_amr
               else
                  igrid_nbor_amr = morton_nbor_grid(igrid_amr,ilevel,igshift)
               end if

               if(igrid_nbor_amr==0) then
                  ! No neighbor (rare boundary case): interp. phi
                  nb_mask = -1.0d0
                  ! Only call morton_nbor_cell for this rare case
                  ifathercell_nbor_amr = morton_nbor_cell(igrid_amr,ilevel,igshift)
                  ind_cell(1)=ifathercell_nbor_amr
                  call interpol_phi(ind_cell,phi_int,1,ilevel,icount)
                  nb_phi = phi_int(1,jjj(idim,inbor,ind))
               else
                  ! Fetch neighbor cell id
                  icell_nbor_amr = igrid_nbor_amr + (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  ! Check neighbor cell mask
                  nb_mask = f(icell_nbor_amr,3)
                  if(nb_mask>0) cycle ! Neighbor cell is active too: cycle
                  nb_phi  = phi(icell_nbor_amr)
               end if
               ! phi(#) interpolated with mask:
               w = nb_mask/(nb_mask-f(icell_amr,3)) ! Linear parameter
               phi_b = ((1.0d0-w)*nb_phi + w*phi(icell_amr))

               ! Increment correction for current cell
               f(icell_amr,2) = f(icell_amr,2) - 2.0d0*oneoverdx2*phi_b
            end do
         end do
      end do
!$omp end parallel do
   end do
#ifdef FDMDEBUG
   call mgp_stop(MGP_BCRHS,mgp_wall_start,mgp_cpu_start, &
        int(ngrid,kind=8)*int(twotondim,kind=8))
#endif

end subroutine make_fine_bc_rhs


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! MPI routines for MG communication for CPU boundaries,
! Those are the MG versions of the make_virtual_* AMR routines
! ------------------------------------------------------------------------

subroutine make_virtual_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
#ifdef FDMDEBUG
  use mg_omp_profile_m
#endif

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifdef FDMDEBUG
  real(kind=8) :: mgp_wall_start,mgp_cpu_start
  integer(kind=8) :: mgp_nwork

  mgp_nwork=0_8
  do icpu=1,ncpu
     mgp_nwork=mgp_nwork+int(emission_mg(icpu,ilevel)%ngrid,kind=8)
  end do
  mgp_nwork=mgp_nwork*int(twotondim,kind=8)
  call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_virtual_mg_dp_ksec(ivar,ilevel)
#ifdef FDMDEBUG
     call mgp_stop(MGP_MG_FORWARD,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%u(i+step,1)=active_mg(myid,ilevel)%u(icell,ivar)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

#ifdef FDMDEBUG
  call mgp_stop(MGP_MG_FORWARD,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_dp

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
#ifdef FDMDEBUG
  use mg_omp_profile_m
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  integer::icpu,i,j,ncache,iskip,step,icell
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifdef FDMDEBUG
  real(kind=8) :: mgp_wall_start,mgp_cpu_start
  integer(kind=8) :: mgp_nwork

  mgp_nwork=0_8
  do icpu=1,ncpu
     mgp_nwork=mgp_nwork+int(emission_mg(icpu,ilevel)%ngrid,kind=8)
  end do
  mgp_nwork=mgp_nwork*int(twotondim,kind=8)
  call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_virtual_mg_int_ksec(ilevel)
#ifdef FDMDEBUG
     call mgp_stop(MGP_MG_FORWARD,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%f(1,1),ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%f(i+step,1)=active_mg(myid,ilevel)%f(icell,1)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

#ifdef FDMDEBUG
  call mgp_stop(MGP_MG_FORWARD,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_int

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
#ifdef FDMDEBUG
  use mg_omp_profile_m
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifdef FDMDEBUG
  real(kind=8) :: mgp_wall_start,mgp_cpu_start
  integer(kind=8) :: mgp_nwork

  mgp_nwork=0_8
  do icpu=1,ncpu
     mgp_nwork=mgp_nwork+int(emission_mg(icpu,ilevel)%ngrid,kind=8)
  end do
  mgp_nwork=mgp_nwork*int(twotondim,kind=8)
  call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_reverse_mg_dp_ksec(ivar,ilevel)
#ifdef FDMDEBUG
     call mgp_stop(MGP_MG_REVERSE,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%u(icell,ivar)=active_mg(myid,ilevel)%u(icell,ivar)+ &
                   & emission_mg(icpu,ilevel)%u(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

#ifdef FDMDEBUG
  call mgp_stop(MGP_MG_REVERSE,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_dp

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
#ifdef FDMDEBUG
  use mg_omp_profile_m
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifdef FDMDEBUG
  real(kind=8) :: mgp_wall_start,mgp_cpu_start
  integer(kind=8) :: mgp_nwork

  mgp_nwork=0_8
  do icpu=1,ncpu
     mgp_nwork=mgp_nwork+int(emission_mg(icpu,ilevel)%ngrid,kind=8)
  end do
  mgp_nwork=mgp_nwork*int(twotondim,kind=8)
  call mgp_start(mgp_wall_start,mgp_cpu_start)
#endif

#ifndef WITHOUTMPI
  if(ordering=='ksection') then
     call make_reverse_mg_int_ksec(ilevel)
#ifdef FDMDEBUG
     call mgp_stop(MGP_MG_REVERSE,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif
     return
  end if

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%f(icell,1)=active_mg(myid,ilevel)%f(icell,1)+&
                 & emission_mg(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

#ifdef FDMDEBUG
  call mgp_stop(MGP_MG_REVERSE,mgp_wall_start,mgp_cpu_start,mgp_nwork)
#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_int

! ########################################################################
! ########################################################################
! Ksection-based MG communication routines
! ########################################################################
! ########################################################################

subroutine make_virtual_mg_dp_ksec(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ivar,ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based forward MG ghost zone exchange for double precision.
  ! Packs emission_mg data with metadata, exchanges via ksection tree,
  ! then scatters to active_mg(sender) grids on the receiver.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(myid,ilevel)%ngrid
           sendbuf(j, idx) = active_mg(myid,ilevel)%u( &
                & emission_mg(icpu,ilevel)%igrid(i) + step, ivar)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data to active_mg(sender) on this CPU
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step = (j-1)*active_mg(sender,ilevel)%ngrid
        active_mg(sender,ilevel)%u(ridx + step, ivar) = recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_mg_dp_ksec

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_int_ksec(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based forward MG ghost zone exchange for integer arrays.
  ! Converts int to dp, exchanges via ksection tree, converts back.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total emission items
  ntotal = 0
  do icpu = 1, ncpu
     ntotal = ntotal + emission_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     do i = 1, emission_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(myid,ilevel)%ngrid
           sendbuf(j, idx) = dble(active_mg(myid,ilevel)%f( &
                & emission_mg(icpu,ilevel)%igrid(i) + step, 1))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Scatter received data (dp->int) to active_mg(sender)%f
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step = (j-1)*active_mg(sender,ilevel)%ngrid
        active_mg(sender,ilevel)%f(ridx + step, 1) = nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_virtual_mg_int_ksec

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp_ksec(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ivar,ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based reverse MG ghost zone exchange for double precision.
  ! Packs active_mg(icpu) data, exchanges via ksection tree,
  ! then accumulates (+=) into local active_mg(myid) using emission_mg.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step,icell
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total items from remote active_mg
  ntotal = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     ntotal = ntotal + active_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from active_mg(icpu) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     do i = 1, active_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(icpu,ilevel)%ngrid
           sendbuf(j, idx) = active_mg(icpu,ilevel)%u(i + step, ivar)
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data into local active_mg(myid) using emission_mg
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step  = (j-1)*active_mg(myid,ilevel)%ngrid
        icell = emission_mg(sender,ilevel)%igrid(ridx) + step
        active_mg(myid,ilevel)%u(icell, ivar) = &
             & active_mg(myid,ilevel)%u(icell, ivar) + recvbuf(j, i)
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_reverse_mg_dp_ksec

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_int_ksec(ilevel)
  use amr_commons
  use poisson_commons
  use ksection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  ! -------------------------------------------------------------------
  ! Ksection-based reverse MG ghost zone exchange for integer arrays.
  ! Packs active_mg(icpu)%f, exchanges via ksection tree,
  ! then accumulates (+=) into local active_mg(myid)%f using emission_mg.
  ! -------------------------------------------------------------------
  integer::icpu,i,j,idx,ntotal,nrecv,nprops_ksec,sender,ridx,step,icell
  real(dp),allocatable::sendbuf(:,:),recvbuf(:,:)
  integer,allocatable::dest_cpu(:)

#ifndef WITHOUTMPI
  nprops_ksec = twotondim + 2

  ! Count total items from remote active_mg
  ntotal = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     ntotal = ntotal + active_mg(icpu,ilevel)%ngrid
  end do

  ! Pack sendbuf from active_mg(icpu)%f (int->dp) + dest_cpu
  allocate(sendbuf(1:nprops_ksec, 1:max(ntotal,1)))
  allocate(dest_cpu(1:max(ntotal,1)))
  idx = 0
  do icpu = 1, ncpu
     if(icpu == myid) cycle
     do i = 1, active_mg(icpu,ilevel)%ngrid
        idx = idx + 1
        dest_cpu(idx) = icpu
        do j = 1, twotondim
           step = (j-1)*active_mg(icpu,ilevel)%ngrid
           sendbuf(j, idx) = dble(active_mg(icpu,ilevel)%f(i + step, 1))
        end do
        sendbuf(twotondim+1, idx) = dble(myid)
        sendbuf(twotondim+2, idx) = dble(i)
     end do
  end do

  ! Exchange via ksection tree
  call ksection_exchange_dp(sendbuf, ntotal, dest_cpu, nprops_ksec, &
       & recvbuf, nrecv)

  ! Accumulate received data (dp->int) into active_mg(myid)%f
  do i = 1, nrecv
     sender = nint(recvbuf(twotondim+1, i))
     ridx   = nint(recvbuf(twotondim+2, i))
     do j = 1, twotondim
        step  = (j-1)*active_mg(myid,ilevel)%ngrid
        icell = emission_mg(sender,ilevel)%igrid(ridx) + step
        active_mg(myid,ilevel)%f(icell, 1) = &
             & active_mg(myid,ilevel)%f(icell, 1) + nint(recvbuf(j, i))
     end do
  end do

  deallocate(sendbuf, dest_cpu, recvbuf)
#endif

end subroutine make_reverse_mg_int_ksec

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine dump_mg_levels(ilevel,idout)
   use amr_commons
   use poisson_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'  
#endif
   integer, intent(in) :: idout, ilevel

   character(len=24)  :: cfile
   character(len=5)   :: ccpu='00000'
   character(len=5)   :: cout='00000'

   integer :: i, ngrids, igrid, icpu, idim
   
   integer,parameter::tag=1119
   integer::dummy_io,info2

   write(ccpu,'(I5.5)') myid
   write(cout,'(I5.5)') idout
   cfile='multigrid_'//cout//'.out'//ccpu
   
   ! Wait for the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
         call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
              & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
   endif
#endif
   
   open(unit=10,file=cfile,status='unknown',form='formatted')

   write(10,'(I1)') ndim
   write(10,'(I1)') myid
   write(10,'(I1)') ncpu
   write(10,'(I2)') ilevel

   ! Loop over levels
   do i=1,ilevel-1
      ! Active grids
      ngrids=active_mg(myid,i)%ngrid
      write(10,*) ngrids
      do igrid=1,ngrids
         do idim=1,ndim
            write(10,*) xg(active_mg(myid,i)%igrid(igrid),idim)
         end do
      end do

      ! Reception grids
      do icpu=1,ncpu
         if(icpu==myid)cycle
         ngrids=active_mg(icpu,i)%ngrid
         write(10,*) ngrids
         do igrid=1,ngrids
            do idim=1,ndim
               write(10,*) xg(active_mg(icpu,i)%igrid(igrid),idim)
            end do
         end do
      end do

   end do

   close(10)

        ! Send the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
         dummy_io=1
         call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
              & MPI_COMM_WORLD,info2)
      end if
   endif
#endif

end subroutine dump_mg_levels


#ifdef HYDRO_CUDA
! ########################################################################
! cuFFT direct Poisson solver for fully uniform levels (periodic BC)
! Replaces MG V-cycle iteration with single FFT solve: O(N log N) vs O(N * niter)
! ########################################################################

subroutine fft_poisson_solve_uniform(ilevel, icount)
   use amr_commons
   use poisson_commons
   use poisson_parameters
   use poisson_cuda_interface
   use neutrino_commons, only: nu_table_loaded, get_nu_ratio
   use dark_energy_commons, only: de_table_loaded, get_de_ratio, &
        f_de_val, compute_de_kspace_params, de_helmholtz_on
   use scalar_de_commons, only: sde_dmcorr_of_a, horndeski_mu_of_a
   use iso_c_binding

  use amr_index, only: icell_of
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif

   integer, intent(in) :: ilevel, icount

   ! Grid dimensions
   integer :: fft_Nx, fft_Ny, fft_Nz
   integer(i8b) :: N_total
   real(dp) :: dx_fft, dx2_fft

   ! Persistent arrays (saved across calls to avoid reallocation)
   integer,  allocatable, save :: fft_map(:)
   real(dp), allocatable, save :: rhs_3d(:), rhs_local(:)
   integer(i8b), save :: saved_N_total = 0

   ! Neutrino/DE correction
   real(dp) :: fourpi_fft, scale_fft
   integer  :: nx_loc_fft
   real(dp) :: omega_cb_loc, nu_factor, R_nu_val, k_phys, twopi_fft
   real(dp) :: de_factor, R_DE_val, omega_de_a
   real(dp) :: hs_factor, mu_a_hs
   real(dp) :: de_kappa2, de_alpha, k_tilde_sq
   integer  :: kx_i, ky_i, kz_i
   integer(i8b) :: N_complex_fft, idx_c

   ! Loop variables
   integer :: igrid, ngrid_loc, ind
   integer :: icell_amr, igrid_amr
   integer :: ix, iy, iz, idx_3d
   integer :: Kx, Ky, Kz
   integer :: info

   ! cuFFTMp dispatch
   logical :: use_fftmp
#ifdef USE_CUFFTMP
   integer :: local_Nx_mp, x_start_mp, local_Ny_mp, y_start_mp
   integer :: slab_size_mp
   real(dp), allocatable, save :: rhs_slab(:), phi_slab(:)
   integer(i8b), save :: saved_slab_size = 0

   ! MPI_ALLTOALLV variables for slab exchange
   integer, allocatable :: sendcounts(:), sdispls(:)
   integer, allocatable :: recvcounts(:), rdispls(:)
   real(dp), allocatable :: sendbuf(:), recvbuf(:)
   integer :: dest_rank, x_local, base_Nx, rem_Nx
   integer :: n_send_total, n_recv_total
   integer :: irank
   ! Scratch for counting
   integer, allocatable :: send_idx(:)
#endif

   ! Grid dimensions for this level (cell grid, not oct grid)
   fft_Nx = nx * 2**ilevel
   fft_Ny = ny * 2**ilevel
   fft_Nz = nz * 2**ilevel
   N_total = int(fft_Nx, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b)

   dx_fft  = 0.5d0**ilevel
   dx2_fft = dx_fft * dx_fft

   ! ------------------------------------------------------------------
   ! Decide: cuFFTMp (distributed) vs cuFFT (ALLREDUCE)
   ! ------------------------------------------------------------------
   use_fftmp = .false.
#ifdef USE_CUFFTMP
   if(N_total > 256_i8b**3 .and. ncpu > 1) then
      ! Only try if not permanently failed
      if(cuda_fftmp_is_ready_c() >= 0) use_fftmp = .true.
   end if
#endif

   if(myid==1) write(*,'(A,I3,A,I5,A,I5,A,I5,A,I15,A,L1)') &
        ' FFT Poisson: level=', ilevel, &
        ' grid=', fft_Nx, 'x', fft_Ny, 'x', fft_Nz, &
        ' N=', N_total, ' distributed=', use_fftmp

#ifdef USE_CUFFTMP
   if(use_fftmp) then
      ! ================================================================
      ! cuFFTMp distributed path
      ! ================================================================

      ! Setup cuFFTMp plan (once per grid size change)
      call cuda_fftmp_poisson_setup_c(int(MPI_COMM_WORLD, c_int), &
           int(fft_Nx, c_int), int(fft_Ny, c_int), int(fft_Nz, c_int), dx2_fft)

      ! Check if setup succeeded; if not, fall back to ALLREDUCE+cuFFT
      if(cuda_fftmp_is_ready_c() /= 1) then
         if(myid==1) write(*,'(A)') &
              '   cuFFTMp setup FAILED — falling back to ALLREDUCE + cuFFT'
         use_fftmp = .false.
         goto 100  ! jump to ALLREDUCE+cuFFT path
      end if

      ! Get local slab sizes from cuFFTMp
      call cuda_fftmp_get_local_sizes_c(local_Nx_mp, x_start_mp, &
           local_Ny_mp, y_start_mp)

      slab_size_mp = local_Nx_mp * fft_Ny * fft_Nz

      if(myid==1) write(*,'(A,I5,A,I5,A,I5,A,I5)') &
           '   cuFFTMp: local_Nx=', local_Nx_mp, ' x_start=', x_start_mp, &
           ' local_Ny=', local_Ny_mp, ' y_start=', y_start_mp

      ! Allocate slab arrays
      if(int(slab_size_mp,i8b) /= saved_slab_size) then
         if(allocated(rhs_slab)) deallocate(rhs_slab)
         if(allocated(phi_slab)) deallocate(phi_slab)
         allocate(rhs_slab(0:max(slab_size_mp,1)-1))
         allocate(phi_slab(0:max(slab_size_mp,1)-1))
         saved_slab_size = int(slab_size_mp, i8b)
      end if

      ! Compute X-slab decomposition parameters
      base_Nx = fft_Nx / ncpu
      rem_Nx  = mod(fft_Nx, ncpu)

      ! ------------------------------------------------------------------
      ! Step 1: Gather local RHS and prepare for ALLTOALLV exchange
      ! Each rank sends its cells to the rank that owns the X-slab
      ! ------------------------------------------------------------------
      allocate(sendcounts(0:ncpu-1))
      allocate(sdispls(0:ncpu-1))
      allocate(recvcounts(0:ncpu-1))
      allocate(rdispls(0:ncpu-1))

      ! Count how many cells this rank sends to each X-slab owner
      sendcounts = 0
      ngrid_loc = active(ilevel)%ngrid
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))

         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            ix = modulo(ix, fft_Nx)

            ! Which rank owns this X-plane?
            if(rem_Nx > 0 .and. ix < rem_Nx * (base_Nx + 1)) then
               dest_rank = ix / (base_Nx + 1)
            else
               dest_rank = rem_Nx + (ix - rem_Nx*(base_Nx+1)) / max(base_Nx,1)
            end if
            dest_rank = min(dest_rank, ncpu-1)
            sendcounts(dest_rank) = sendcounts(dest_rank) + 1
         end do
      end do

      ! Exchange counts
      call MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER, &
           recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, info)

      ! Compute displacements
      sdispls(0) = 0
      rdispls(0) = 0
      do irank = 1, ncpu-1
         sdispls(irank) = sdispls(irank-1) + sendcounts(irank-1)
         rdispls(irank) = rdispls(irank-1) + recvcounts(irank-1)
      end do
      n_send_total = sdispls(ncpu-1) + sendcounts(ncpu-1)
      n_recv_total = rdispls(ncpu-1) + recvcounts(ncpu-1)

      ! Pack send buffer: (idx_in_slab, rhs_value) interleaved
      ! We'll send 2 doubles per cell: slab index + value
      allocate(sendbuf(0:2*n_send_total-1))
      allocate(recvbuf(0:2*n_recv_total-1))
      allocate(send_idx(0:ncpu-1))
      send_idx = sdispls

      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
         Ky = nint(xg(igrid_amr, 2) * dble(fft_Ny))
         Kz = nint(xg(igrid_amr, 3) * dble(fft_Nz))

         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            iy = Ky - 1 + mod((ind-1)/2, 2)
            iz = Kz - 1 + (ind-1)/4
            ix = modulo(ix, fft_Nx)
            iy = modulo(iy, fft_Ny)
            iz = modulo(iz, fft_Nz)

            if(rem_Nx > 0 .and. ix < rem_Nx * (base_Nx + 1)) then
               dest_rank = ix / (base_Nx + 1)
            else
               dest_rank = rem_Nx + (ix - rem_Nx*(base_Nx+1)) / max(base_Nx,1)
            end if
            dest_rank = min(dest_rank, ncpu-1)

            ! Local X index within destination slab
            if(dest_rank < rem_Nx) then
               x_local = ix - dest_rank * (base_Nx + 1)
            else
               x_local = ix - rem_Nx*(base_Nx+1) - (dest_rank-rem_Nx)*base_Nx
            end if

            ! Slab-local row-major index
            idx_3d = x_local * fft_Ny * fft_Nz + iy * fft_Nz + iz

            icell_amr = icell_of(igrid_amr,ind)

            sendbuf(2*send_idx(dest_rank))     = dble(idx_3d)
            sendbuf(2*send_idx(dest_rank) + 1) = f(icell_amr, 2)
            send_idx(dest_rank) = send_idx(dest_rank) + 1
         end do
      end do

      ! Exchange (2 doubles per cell)
      sendcounts = sendcounts * 2
      recvcounts = recvcounts * 2
      sdispls = sdispls * 2
      rdispls = rdispls * 2

      call MPI_ALLTOALLV(sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
           recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
           MPI_COMM_WORLD, info)

      ! Unpack into local slab
      rhs_slab = 0.0d0
      do ix = 0, n_recv_total - 1
         idx_3d = nint(recvbuf(2*ix))
         rhs_slab(idx_3d) = rhs_slab(idx_3d) + recvbuf(2*ix + 1)
      end do

      ! ------------------------------------------------------------------
      ! Step 2: cuFFTMp distributed solve
      ! ------------------------------------------------------------------
      call cuda_fftmp_poisson_solve_c(rhs_slab, phi_slab, &
           int(local_Nx_mp, c_int), int(fft_Ny, c_int), int(fft_Nz, c_int), &
           int(y_start_mp, c_int), int(local_Ny_mp, c_int))

      ! ------------------------------------------------------------------
      ! Step 3: Send solved phi back to original cell owners
      ! Reuse ALLTOALLV: slab owner sends phi to cell owner
      ! ------------------------------------------------------------------
      ! recvbuf already has (idx_in_slab, _) entries.
      ! Replace rhs values with phi values from phi_slab
      do ix = 0, n_recv_total - 1
         idx_3d = nint(recvbuf(2*ix))
         recvbuf(2*ix + 1) = phi_slab(idx_3d)
      end do

      ! Reverse ALLTOALLV: send phi back
      ! Now recv→send, send→recv
      call MPI_ALLTOALLV(recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
           sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
           MPI_COMM_WORLD, info)

      ! Unpack phi to RAMSES cells
      send_idx = sdispls / 2  ! Reset to element indices (not byte offsets)

      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))

         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            ix = modulo(ix, fft_Nx)

            if(rem_Nx > 0 .and. ix < rem_Nx * (base_Nx + 1)) then
               dest_rank = ix / (base_Nx + 1)
            else
               dest_rank = rem_Nx + (ix - rem_Nx*(base_Nx+1)) / max(base_Nx,1)
            end if
            dest_rank = min(dest_rank, ncpu-1)

            icell_amr = icell_of(igrid_amr,ind)

            ! phi value is in sendbuf at the same position we packed it
            phi(icell_amr) = sendbuf(2*send_idx(dest_rank) + 1)
            send_idx(dest_rank) = send_idx(dest_rank) + 1
         end do
      end do

      deallocate(sendcounts, sdispls, recvcounts, rdispls)
      deallocate(sendbuf, recvbuf, send_idx)

      return
   end if
#endif

   ! ================================================================
   ! cuFFT ALLREDUCE path (default, for grids ≤256³ or single rank)
   ! Also used as fallback when cuFFTMp fails
   ! ================================================================
100 continue

   ! ------------------------------------------------------------------
   ! Persistent allocation (only allocate on first call or size change)
   ! ------------------------------------------------------------------
   if(N_total /= saved_N_total) then
      if(allocated(fft_map))   deallocate(fft_map)
      if(allocated(rhs_local)) deallocate(rhs_local)
      if(allocated(rhs_3d))    deallocate(rhs_3d)
      allocate(fft_map(0:N_total-1))
      allocate(rhs_local(0:N_total-1))
      allocate(rhs_3d(0:N_total-1))
      saved_N_total = N_total
   end if

   ! ------------------------------------------------------------------
   ! Step 1: Build fft_map and gather local RHS
   ! fft_map(idx_3d) = icell_amr (1-based), idx_3d is 0-based C index
   ! ------------------------------------------------------------------
   fft_map   = 0
   rhs_local = 0.0d0

   ! Compute fourpi for FFT RHS (NO DE boost — corrections go in Green's function)
   nx_loc_fft = icoarse_max - icoarse_min + 1
   scale_fft  = boxlen / dble(nx_loc_fft)
   if(cosmo) then
      if(use_neutrino .and. omega_nu > 0.0d0) then
         fourpi_fft = 1.5D0 * (omega_m - omega_nu) * aexp * scale_fft
      else
         fourpi_fft = 1.5D0 * omega_m * aexp * scale_fft
      end if
      ! Coupled quintessence: DM mass evolution rho_dm*a^3/rho_dm0
      if(use_coupled_de .and. cde_vary_mass .and. use_quintessence) &
           & fourpi_fft = fourpi_fft * sde_dmcorr_of_a(aexp)
      ! Horndeski scale-independent mu(a); k-dependent case (hs_mass>0)
      ! is applied in the Green's-function correction instead
      if(use_horndeski .and. hs_mass == 0.0d0) &
           & fourpi_fft = fourpi_fft * horndeski_mu_of_a(aexp)
   else
      fourpi_fft = 4.D0 * ACOS(-1.0D0) * scale_fft
   end if

   ngrid_loc = active(ilevel)%ngrid
   do igrid = 1, ngrid_loc
      igrid_amr = active(ilevel)%igrid(igrid)

      ! Grid center coordinates -> integer cell coords
      Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
      Ky = nint(xg(igrid_amr, 2) * dble(fft_Ny))
      Kz = nint(xg(igrid_amr, 3) * dble(fft_Nz))

      do ind = 1, twotondim
         ! Child cell offset within oct (0-based)
         ix = Kx - 1 + mod(ind-1, 2)
         iy = Ky - 1 + mod((ind-1)/2, 2)
         iz = Kz - 1 + (ind-1)/4

         ! Periodic wrapping
         ix = modulo(ix, fft_Nx)
         iy = modulo(iy, fft_Ny)
         iz = modulo(iz, fft_Nz)

         ! C row-major index: idx = ix * Ny * Nz + iy * Nz + iz
         idx_3d = ix * fft_Ny * fft_Nz + iy * fft_Nz + iz

         ! RAMSES cell index
         icell_amr = icell_of(igrid_amr,ind)

         fft_map(idx_3d) = icell_amr
         ! Use un-boosted fourpi (corrections in Green's function)
         rhs_local(idx_3d) = fourpi_fft * (rho(icell_amr) - rho_tot)
      end do
   end do

   ! ------------------------------------------------------------------
   ! Step 2: MPI_Allreduce to get global RHS (all ranks contribute local parts)
   ! ------------------------------------------------------------------
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(rhs_local, rhs_3d, int(N_total), &
        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
   rhs_3d = rhs_local
#endif

   ! ------------------------------------------------------------------
   ! Step 3: cuFFT setup (plans + Green's function, only on grid change)
   ! ------------------------------------------------------------------
   call cuda_fft_poisson_setup_c(fft_map, &
        int(fft_Nx, c_int), int(fft_Ny, c_int), int(fft_Nz, c_int), dx2_fft)

   ! ------------------------------------------------------------------
   ! Step 3b: Compute and upload neutrino/DE correction factors to GPU
   ! ------------------------------------------------------------------
   if((use_neutrino .and. omega_nu > 0.0d0) .or. de_perturb &
        & .or. (use_horndeski .and. hs_mass > 0.0d0)) then
      mu_a_hs = horndeski_mu_of_a(aexp)
      twopi_fft = 2.0d0 * ACOS(-1.0D0)
      N_complex_fft = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz/2+1,i8b)
      if(use_neutrino .and. omega_nu > 0.0d0) omega_cb_loc = omega_m - omega_nu
      if(.not. allocated(rhs_local)) allocate(rhs_local(0:N_complex_fft-1))
      ! Reuse rhs_local as scratch (it's already allocated >= N_complex)

      ! Precompute DE params
      if(de_perturb .and. .not. de_table_loaded .and. de_helmholtz_on()) then
         call compute_de_kspace_params(aexp, de_kappa2, de_alpha)
      end if
      if(de_perturb .and. de_table_loaded) then
         omega_de_a = omega_l * f_de_val(aexp) * aexp**3
         if(use_neutrino .and. omega_nu > 0.0d0) then
            omega_cb_loc = omega_m - omega_nu
         else
            omega_cb_loc = omega_m
         end if
      end if

      do idx_c = 0, N_complex_fft - 1
         kx_i = int(idx_c / (int(fft_Ny,i8b)*int(fft_Nz/2+1,i8b)))
         ky_i = int(mod(idx_c / int(fft_Nz/2+1,i8b), int(fft_Ny,i8b)))
         kz_i = int(mod(idx_c, int(fft_Nz/2+1,i8b)))
         if(kx_i > fft_Nx/2) kx_i = kx_i - fft_Nx
         if(ky_i > fft_Ny/2) ky_i = ky_i - fft_Ny

         nu_factor = 1.0d0
         de_factor = 1.0d0

         ! Neutrino correction
         if(use_neutrino .and. omega_nu > 0.0d0) then
            k_phys = twopi_fft * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
            if(k_phys > 0.0d0) then
               R_nu_val = get_nu_ratio(k_phys, aexp)
               nu_factor = 1.0d0 + (omega_nu / omega_cb_loc) * R_nu_val
            end if
         end if

         ! DE perturbation correction
         if(de_perturb) then
            if(de_table_loaded) then
               k_phys = twopi_fft * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
               if(k_phys > 0.0d0) then
                  R_DE_val = get_de_ratio(k_phys, aexp)
                  de_factor = 1.0d0 + (omega_de_a / omega_cb_loc) * R_DE_val
               end if
            else if(de_helmholtz_on()) then
               k_tilde_sq = twopi_fft * twopi_fft * dble(kx_i**2 + ky_i**2 + kz_i**2)
               de_factor = 1.0d0 + de_alpha / (k_tilde_sq + de_kappa2)
            end if
         end if

         ! Horndeski k-dependent mu(a,k)
         hs_factor = 1.0d0
         if(use_horndeski .and. hs_mass > 0.0d0) then
            k_phys = twopi_fft * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
            hs_factor = 1.0d0 + (mu_a_hs - 1.0d0) * k_phys*k_phys / &
                 & (k_phys*k_phys + (aexp*hs_mass)**2)
         end if

         rhs_local(idx_c) = nu_factor * de_factor * hs_factor
      end do

      call cuda_fft_set_correction_c(rhs_local, int(N_complex_fft, c_int))
   else
      ! Clear correction (use standard Green's function)
      call cuda_fft_set_correction_c(rhs_local, int(0, c_int))
   end if

   ! ------------------------------------------------------------------
   ! Step 4: cuFFT solve and scatter phi to CPU
   ! Uses cuda_fft_poisson_solve_host_c (no d_mg_phi dependency)
   ! ------------------------------------------------------------------
   call cuda_fft_poisson_solve_host_c(rhs_3d, rhs_3d, int(N_total, c_int))

   ! Scatter from 3D array to RAMSES phi via fft_map
   do idx_3d = 0, int(N_total-1)
      icell_amr = fft_map(idx_3d)
      if(icell_amr > 0) phi(icell_amr) = rhs_3d(idx_3d)
   end do

end subroutine fft_poisson_solve_uniform
#endif


#ifdef USE_FFTW
! ########################################################################
! FFTW3 CPU direct Poisson solver for fully uniform levels (periodic BC)
! Two paths:
!   Small grid (N <= 256^3): MPI_ALLREDUCE + local FFTW3+OMP
!   Large grid (N >  256^3): MPI_ALLTOALLV + FFTW3 MPI slab decomposition
! ########################################################################

subroutine fftw_poisson_solve_uniform(ilevel, icount)
   use amr_commons, only: dp
   implicit none
   integer,intent(in)::ilevel,icount

   call fftw_uniform_solve_engine(ilevel,icount,0,0d0,0d0,1d0)
end subroutine fftw_poisson_solve_uniform

subroutine fftw_scalar_solve_uniform(ilevel,m2,step_frac,relax)
   use amr_commons, only: dp
   implicit none
   integer,intent(in)::ilevel
   real(dp),intent(in)::m2,step_frac,relax

   call fftw_uniform_solve_engine(ilevel,0,1,m2,step_frac,relax)
end subroutine fftw_scalar_solve_uniform

subroutine fftw_uniform_solve_engine(ilevel,icount,solve_mode,m2,step_frac,relax)
   use amr_commons
   use poisson_commons
   use poisson_parameters
   use neutrino_commons, only: nu_table_loaded, read_neutrino_table, get_nu_ratio
   use dark_energy_commons, only: compute_de_kspace_params, de_table_loaded, &
        read_de_table, get_de_ratio, f_de_val, de_helmholtz_on
   use scalar_de_commons, only: sde_dmcorr_of_a, horndeski_mu_of_a
   use iso_c_binding
   use omp_lib
#ifdef FDMDEBUG
   use fftw_omp_profile_m
#endif
  use amr_index, only: icell_of
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif
   include 'fftw3-mpi.f03'

   integer, intent(in) :: ilevel,icount,solve_mode
   real(dp),intent(in)::m2,step_frac,relax

   ! === Grid dimensions ===
   integer :: fft_Nx, fft_Ny, fft_Nz
   integer(i8b) :: N_total
   real(dp) :: dx_fft, dx2_fft
   logical :: use_distributed

   ! === FFTW cached state (persist across calls) ===
   logical, save :: fftw_initialized = .false.
   type(C_PTR), save :: plan_r2c = C_NULL_PTR, plan_c2r = C_NULL_PTR
   integer, save :: saved_Nx = 0, saved_Ny = 0, saved_Nz = 0
   logical, save :: saved_distributed = .false.

   ! Small-grid arrays (save)
   integer,  allocatable, save :: fft_map(:)
   real(dp), allocatable, save :: rhs_local(:), rhs_3d(:)
   real(C_DOUBLE), allocatable, save :: green(:)
   complex(C_DOUBLE_COMPLEX), allocatable, save :: cdata(:)
   integer(i8b), save :: saved_N_total = 0
   integer(i8b) :: N_complex

   ! Large-grid FFTW MPI (save)
   type(C_PTR), save :: p_rdata = C_NULL_PTR, p_cdata = C_NULL_PTR
   real(C_DOUBLE), pointer, save :: rdata_fftw(:) => null()
   complex(C_DOUBLE_COMPLEX), pointer, save :: cdata_fftw(:) => null()
   integer(C_INTPTR_T), save :: local_nx_fftw = 0, nx_start_fftw = 0
   integer(C_INTPTR_T), save :: local_ny_fftw = 0, ny_start_fftw = 0  ! transposed-out (Ny-distributed)
   integer(C_INTPTR_T), save :: alloc_local = 0
   integer, save :: fftw_block = 0  ! FFTW MPI block size for slab distribution

   ! Sparse P2P partner lists (precomputed, cached)
   integer, allocatable, save :: fftw_send_partners(:)  ! 0-based MPI ranks
   integer, allocatable, save :: fftw_recv_partners(:)  ! 0-based MPI ranks
   integer, save :: n_fftw_send = 0, n_fftw_recv = 0
   logical, save :: fftw_partners_computed = .false.
   real(dp), save :: cached_box_xmin = -1.0d0
   real(dp), save :: cached_box_xmax = -1.0d0
   integer, save :: cached_fftw_block = 0
   real(dp), allocatable, save :: fftw_cpubox_min(:), fftw_cpubox_max(:)

   ! Loop variables
   integer :: igrid, ngrid_loc, ind
   integer :: icell_amr, igrid_amr
   integer :: ix, iy, iz, idx_3d
   integer :: Kx, Ky, Kz
   integer :: info, ierr
   integer :: kx_i, ky_i, kz_i
   real(dp) :: denom, twopi, scalar_green, uold, du
   integer(i8b) :: idx_c

   ! Neutrino linear response variables
   real(dp) :: omega_cb_loc, nu_factor, R_nu_val, k_phys

   ! DE perturbation variables
   real(dp) :: de_kappa2, de_alpha, k_tilde_sq, de_factor
   real(dp) :: hs_factor, mu_a_hs
   real(dp) :: R_DE_val, omega_de_a

   ! FFT RHS fourpi (without DE boost)
   real(dp) :: fourpi_fft, scale_fft
   integer  :: nx_loc_fft

   ! MPI ALLTOALLV variables for large-grid path
   integer, allocatable :: sendcounts(:), sdispls(:)
   integer, allocatable :: recvcounts(:), rdispls(:)
   real(dp), allocatable :: sendbuf(:), recvbuf(:)
   integer, allocatable :: thread_counts(:,:), thread_offsets(:,:)
   integer :: dest_rank, x_local, base_Nx, rem_Nx
   integer :: n_send_total, n_recv_total, irank
   integer :: local_Nz_half, slab_real_size
   integer :: nthreads_fft, tid, slot

   ! Sparse P2P variables
   integer, allocatable :: reqs(:)
   integer, allocatable :: mpi_stat(:,:)
   integer :: nreq, ip
   integer :: slab_x_lo, slab_x_hi
   real(dp) :: recv_pos_lo, recv_pos_hi
   real(dp) :: my_xmin, my_xmax
   logical :: overlap
   logical :: local_changed, global_changed
#ifdef FDMDEBUG
   integer :: fftp_path
   real(kind=8) :: fftp_total_wall,fftp_total_cpu
   real(kind=8) :: fftp_phase_wall,fftp_phase_cpu
#endif

   ! ================================================================
   ! Step 0: Grid dimensions
   ! ================================================================
   fft_Nx = nx * 2**ilevel
   fft_Ny = ny * 2**ilevel
   fft_Nz = nz * 2**ilevel
   N_total = int(fft_Nx, i8b) * int(fft_Ny, i8b) * int(fft_Nz, i8b)
   dx_fft  = 0.5d0**ilevel
   dx2_fft = dx_fft * dx_fft
   twopi   = 2.0d0 * acos(-1.0d0)
   use_distributed = (N_total >= 256_i8b**3 .and. ncpu > 1)
#ifdef FDMDEBUG
   fftp_path=merge(2,1,use_distributed)
   call fftp_start(fftp_total_wall,fftp_total_cpu)
   call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

   if(myid==1) write(*,'(A,A,A,I3,A,I5,A,I5,A,I5,A,I15,A,L1)') &
        ' FFTW3 ',merge('scalar ','Poisson',solve_mode==1),': level=',ilevel, &
        ' grid=', fft_Nx, 'x', fft_Ny, 'x', fft_Nz, &
        ' N=', N_total, ' distributed=', use_distributed

   ! ================================================================
   ! Step 0b: Load neutrino table (once)
   ! ================================================================
   if(use_neutrino .and. .not. nu_table_loaded) then
      call read_neutrino_table(trim(neutrino_table))
      if(myid==1) write(*,'(A)') '   Neutrino transfer function table loaded'
   end if

   ! Load DE table (once)
   if(de_perturb .and. len_trim(de_table) > 0 .and. .not. de_table_loaded) then
      call read_de_table(trim(de_table))
      if(myid==1) write(*,'(A)') '   DE transfer function table loaded'
   end if

   ! ================================================================
   ! Step 1: FFTW initialization (once)
   ! ================================================================
   if(.not. fftw_initialized) then
      ierr = fftw_init_threads()
      call fftw_plan_with_nthreads(int(omp_get_max_threads(), C_INT))
      call fftw_mpi_init()
      fftw_initialized = .true.
      if(myid==1) write(*,'(A,I3,A)') &
           '   FFTW3 initialized with ', omp_get_max_threads(), ' threads'
   end if

   ! ================================================================
   ! Step 2: Plan creation (only on grid size change)
   ! ================================================================
   if(fft_Nx /= saved_Nx .or. fft_Ny /= saved_Ny .or. fft_Nz /= saved_Nz &
      .or. (use_distributed .neqv. saved_distributed)) then
      ! Destroy old plans
      if(c_associated(plan_r2c)) call fftw_destroy_plan(plan_r2c)
      if(c_associated(plan_c2r)) call fftw_destroy_plan(plan_c2r)
      plan_r2c = C_NULL_PTR
      plan_c2r = C_NULL_PTR
      saved_Nx = fft_Nx
      saved_Ny = fft_Ny
      saved_Nz = fft_Nz
      saved_distributed = use_distributed

      if(use_distributed) then
         ! --- Large-grid: FFTW3 MPI plans (in-place R2C) ---
         ! Free old allocation (single block for in-place)
         if(c_associated(p_rdata)) call fftw_free(p_rdata)
         p_rdata = C_NULL_PTR
         p_cdata = C_NULL_PTR

         ! Query local sizes for TRANSPOSED layout:
         !   real input  : distributed over n0=Nx (local_nx_fftw / nx_start_fftw)
         !   complex out : distributed over n1=Ny (local_ny_fftw / ny_start_fftw)
         ! Transposed ordering skips the second global all-to-all transpose,
         ! roughly halving FFT communication (critical at large N^3 / many ranks).
         ! The Green's function is pointwise in k-space, so the transposed
         ! complex layout is harmless as long as k-indices are computed for it.
         alloc_local = fftw_mpi_local_size_3d_transposed( &
              int(fft_Nx, C_INTPTR_T), &
              int(fft_Ny, C_INTPTR_T), &
              int(fft_Nz/2+1, C_INTPTR_T), &
              MPI_COMM_WORLD, local_nx_fftw, nx_start_fftw, &
              local_ny_fftw, ny_start_fftw)

         ! In-place R2C: single allocation, real & complex share memory
         ! Real layout: local_nx * Ny * 2*(Nz/2+1) (padded last dim)
         ! Complex layout: local_nx * Ny * (Nz/2+1)
         p_rdata = fftw_alloc_complex(alloc_local)
         p_cdata = p_rdata  ! Same memory for in-place
         slab_real_size = 2 * int(alloc_local)
         call c_f_pointer(p_rdata, rdata_fftw, [slab_real_size])
         call c_f_pointer(p_cdata, cdata_fftw, [int(alloc_local)])

         ! In-place plans (FFTW detects same base address) + transposed I/O:
         ! forward leaves output Ny-distributed (TRANSPOSED_OUT); inverse takes
         ! the same Ny-distributed layout as input (TRANSPOSED_IN). This omits
         ! the transpose-back, halving the FFT all-to-all communication.
         plan_r2c = fftw_mpi_plan_dft_r2c_3d( &
              int(fft_Nx, C_INTPTR_T), &
              int(fft_Ny, C_INTPTR_T), &
              int(fft_Nz, C_INTPTR_T), &
              rdata_fftw, cdata_fftw, &
              MPI_COMM_WORLD, ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_OUT))

         plan_c2r = fftw_mpi_plan_dft_c2r_3d( &
              int(fft_Nx, C_INTPTR_T), &
              int(fft_Ny, C_INTPTR_T), &
              int(fft_Nz, C_INTPTR_T), &
              cdata_fftw, rdata_fftw, &
              MPI_COMM_WORLD, ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_IN))

         ! FFTW block size: ceil(Nx/ncpu)
         fftw_block = (fft_Nx + ncpu - 1) / ncpu

         if(myid==1) write(*,'(A,I6,A,I6,A,I6)') &
              '   FFTW3 MPI: local_nx=', int(local_nx_fftw), &
              ' nx_start=', int(nx_start_fftw), &
              ' block=', fftw_block

         ! Precompute Green's function for TRANSPOSED complex layout.
         ! After TRANSPOSED_OUT the complex array is Ny-distributed: the
         ! first two dims are swapped, so storage order is
         !   idx_c = ky_local*(Nx*(Nz/2+1)) + kx*(Nz/2+1) + kz
         ! with global ky = ny_start_fftw + ky_local, kx = 0..Nx-1.
         if(allocated(green)) deallocate(green)
         N_complex = int(local_ny_fftw,i8b) * int(fft_Nx,i8b) * int(fft_Nz/2+1,i8b)
         allocate(green(0:N_complex-1))

         do ky_i = 0, int(local_ny_fftw)-1
            do kx_i = 0, fft_Nx-1
               do kz_i = 0, fft_Nz/2
                  idx_c = int(ky_i,i8b)*int(fft_Nx,i8b)*int(fft_Nz/2+1,i8b) &
                        + int(kx_i,i8b)*int(fft_Nz/2+1,i8b) + int(kz_i,i8b)
                  iy = int(ny_start_fftw) + ky_i  ! global ky
                  denom = 2.0d0*(cos(twopi*dble(kx_i)/dble(fft_Nx)) &
                               + cos(twopi*dble(iy)/dble(fft_Ny)) &
                               + cos(twopi*dble(kz_i)/dble(fft_Nz)) - 3.0d0)
                  if(kx_i==0 .and. iy==0 .and. kz_i==0) then
                     green(idx_c) = 0.0d0
                  else
                     green(idx_c) = dx2_fft / denom
                  end if
               end do
            end do
         end do

      else
         ! --- Small-grid: local FFTW3+OMP plans ---
         if(allocated(green)) deallocate(green)
         if(allocated(cdata)) deallocate(cdata)

         N_complex = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz/2+1,i8b)
         allocate(cdata(0:N_complex-1))
         allocate(green(0:N_complex-1))

         ! Plans: fftw3.f03 uses bind(C), no dimension reversal needed
         ! Our data is C row-major (idx = ix*Ny*Nz + iy*Nz + iz)
         ! n0=Nx (slowest), n1=Ny, n2=Nz (fastest, gets R2C reduction)
         ! We need allocated arrays for plan creation
         if(.not. allocated(rhs_3d)) then
            allocate(rhs_3d(0:N_total-1))
         else if(size(rhs_3d) < N_total) then
            deallocate(rhs_3d)
            allocate(rhs_3d(0:N_total-1))
         end if

         plan_r2c = fftw_plan_dft_r2c_3d( &
              int(fft_Nx, C_INT), int(fft_Ny, C_INT), int(fft_Nz, C_INT), &
              rhs_3d, cdata, FFTW_ESTIMATE)

         plan_c2r = fftw_plan_dft_c2r_3d( &
              int(fft_Nx, C_INT), int(fft_Ny, C_INT), int(fft_Nz, C_INT), &
              cdata, rhs_3d, FFTW_ESTIMATE)

         ! Precompute Green's function
         do kx_i = 0, fft_Nx-1
            do ky_i = 0, fft_Ny-1
               do kz_i = 0, fft_Nz/2
                  idx_c = int(kx_i,i8b)*int(fft_Ny,i8b)*int(fft_Nz/2+1,i8b) &
                        + int(ky_i,i8b)*int(fft_Nz/2+1,i8b) + int(kz_i,i8b)
                  denom = 2.0d0*(cos(twopi*dble(kx_i)/dble(fft_Nx)) &
                               + cos(twopi*dble(ky_i)/dble(fft_Ny)) &
                               + cos(twopi*dble(kz_i)/dble(fft_Nz)) - 3.0d0)
                  if(kx_i==0 .and. ky_i==0 .and. kz_i==0) then
                     green(idx_c) = 0.0d0
                  else
                     green(idx_c) = dx2_fft / denom
                  end if
               end do
            end do
         end do
      end if
   end if  ! plan creation

   ! alloc_local is cached with the distributed FFTW plans, whereas
   ! slab_real_size is a local scratch scalar.  Reconstruct it on every
   ! solve so the OpenMP zero/normalization loops cover the complete
   ! in-place FFTW allocation after the first call as well.
   if(use_distributed) slab_real_size = 2 * int(alloc_local)


   ! Compute FFT RHS fourpi (NO DE boost — corrections go in Green's function)
   nx_loc_fft = icoarse_max - icoarse_min + 1
   scale_fft  = boxlen / dble(nx_loc_fft)
   if(cosmo) then
      if(use_neutrino .and. omega_nu > 0.0d0) then
         fourpi_fft = 1.5D0 * (omega_m - omega_nu) * aexp * scale_fft
      else
         fourpi_fft = 1.5D0 * omega_m * aexp * scale_fft
      end if
      ! Coupled quintessence: DM mass evolution rho_dm*a^3/rho_dm0
      if(use_coupled_de .and. cde_vary_mass .and. use_quintessence) &
           & fourpi_fft = fourpi_fft * sde_dmcorr_of_a(aexp)
      ! Horndeski scale-independent mu(a); k-dependent case (hs_mass>0)
      ! is applied in the Green's-function correction instead
      if(use_horndeski .and. hs_mass == 0.0d0) &
           & fourpi_fft = fourpi_fft * horndeski_mu_of_a(aexp)
   else
      fourpi_fft = 4.D0 * ACOS(-1.0D0) * scale_fft
   end if
#ifdef FDMDEBUG
   call fftp_stop(fftp_path,FFTP_SETUP,fftp_phase_wall,fftp_phase_cpu,N_total)
#endif

   if(use_distributed) then
      ! ============================================================
      ! LARGE GRID PATH: Sparse P2P + FFTW3 MPI slab decomposition
      ! Uses ISEND/IRECV with precomputed partner lists instead of
      ! ALLTOALLV — O(n_partners) messages instead of O(ncpu)
      ! ============================================================

      ! --- Precompute P2P partner lists (cached, recomputed on rebalance) ---
      ! Compute my spatial x-bounding box from active grid positions
#ifdef FDMDEBUG
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif
      my_xmin = 1.0d0
      my_xmax = 0.0d0
      ngrid_loc = active(ilevel)%ngrid
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         my_xmin = min(my_xmin, xg(igrid_amr, 1) - dx_fft)
         my_xmax = max(my_xmax, xg(igrid_amr, 1) + dx_fft)
      end do
      if(ngrid_loc == 0) then
         my_xmin = 0.5d0; my_xmax = 0.5d0
      end if

      ! Check if ANY rank's bounding box changed (must be collective
      ! because MPI_ALLGATHER inside requires all ranks to participate)
      local_changed = (.not. fftw_partners_computed .or. &
         cached_box_xmin /= my_xmin .or. &
         cached_box_xmax /= my_xmax .or. &
         cached_fftw_block /= fftw_block)
      call MPI_ALLREDUCE(local_changed, global_changed, 1, MPI_LOGICAL, &
           MPI_LOR, MPI_COMM_WORLD, info)
      if(global_changed) then

         ! Gather all ranks' bounding boxes
         if(.not. allocated(fftw_cpubox_min)) allocate(fftw_cpubox_min(ncpu))
         if(.not. allocated(fftw_cpubox_max)) allocate(fftw_cpubox_max(ncpu))
         call MPI_ALLGATHER(my_xmin, 1, MPI_DOUBLE_PRECISION, &
              fftw_cpubox_min, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
         call MPI_ALLGATHER(my_xmax, 1, MPI_DOUBLE_PRECISION, &
              fftw_cpubox_max, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)

         if(allocated(fftw_send_partners)) deallocate(fftw_send_partners)
         if(allocated(fftw_recv_partners)) deallocate(fftw_recv_partners)

         ! Send partners: FFTW slab ranks whose position range overlaps my domain
         n_fftw_send = 0
         do irank = 0, ncpu - 1
            slab_x_lo = irank * fftw_block
            slab_x_hi = min((irank + 1) * fftw_block, fft_Nx) - 1
            recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
            recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)
            overlap = (fftw_cpubox_max(myid) > recv_pos_lo .and. &
                       fftw_cpubox_min(myid) < recv_pos_hi)
            if(recv_pos_lo < 0.0d0) &
               overlap = overlap .or. (fftw_cpubox_max(myid) > recv_pos_lo + 1.0d0)
            if(recv_pos_hi > 1.0d0) &
               overlap = overlap .or. (fftw_cpubox_min(myid) < recv_pos_hi - 1.0d0)
            if(overlap) n_fftw_send = n_fftw_send + 1
         end do

         allocate(fftw_send_partners(n_fftw_send))
         n_fftw_send = 0
         do irank = 0, ncpu - 1
            slab_x_lo = irank * fftw_block
            slab_x_hi = min((irank + 1) * fftw_block, fft_Nx) - 1
            recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
            recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)
            overlap = (fftw_cpubox_max(myid) > recv_pos_lo .and. &
                       fftw_cpubox_min(myid) < recv_pos_hi)
            if(recv_pos_lo < 0.0d0) &
               overlap = overlap .or. (fftw_cpubox_max(myid) > recv_pos_lo + 1.0d0)
            if(recv_pos_hi > 1.0d0) &
               overlap = overlap .or. (fftw_cpubox_min(myid) < recv_pos_hi - 1.0d0)
            if(overlap) then
               n_fftw_send = n_fftw_send + 1
               fftw_send_partners(n_fftw_send) = irank
            end if
         end do

         ! Recv partners: RAMSES ranks whose domain overlaps my FFTW slab
         slab_x_lo = (myid - 1) * fftw_block
         slab_x_hi = min(myid * fftw_block, fft_Nx) - 1
         recv_pos_lo = (dble(slab_x_lo) - 0.5d0) / dble(fft_Nx)
         recv_pos_hi = (dble(slab_x_hi) + 1.5d0) / dble(fft_Nx)

         n_fftw_recv = 0
         do irank = 1, ncpu
            overlap = (fftw_cpubox_max(irank) > recv_pos_lo .and. &
                       fftw_cpubox_min(irank) < recv_pos_hi)
            if(recv_pos_lo < 0.0d0) &
               overlap = overlap .or. (fftw_cpubox_max(irank) > recv_pos_lo + 1.0d0)
            if(recv_pos_hi > 1.0d0) &
               overlap = overlap .or. (fftw_cpubox_min(irank) < recv_pos_hi - 1.0d0)
            if(overlap) n_fftw_recv = n_fftw_recv + 1
         end do

         allocate(fftw_recv_partners(n_fftw_recv))
         n_fftw_recv = 0
         do irank = 1, ncpu
            overlap = (fftw_cpubox_max(irank) > recv_pos_lo .and. &
                       fftw_cpubox_min(irank) < recv_pos_hi)
            if(recv_pos_lo < 0.0d0) &
               overlap = overlap .or. (fftw_cpubox_max(irank) > recv_pos_lo + 1.0d0)
            if(recv_pos_hi > 1.0d0) &
               overlap = overlap .or. (fftw_cpubox_min(irank) < recv_pos_hi - 1.0d0)
            if(overlap) then
               n_fftw_recv = n_fftw_recv + 1
               fftw_recv_partners(n_fftw_recv) = irank - 1  ! 0-based MPI rank
            end if
         end do

         cached_box_xmin = my_xmin
         cached_box_xmax = my_xmax
         cached_fftw_block = fftw_block
         fftw_partners_computed = .true.

         if(myid==1) write(*,'(A,I4,A,I4,A)') &
              '   FFTW3 sparse P2P: n_send=', n_fftw_send, &
              ', n_recv=', n_fftw_recv, ' partners'
      end if
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_PARTNER,fftp_phase_wall,fftp_phase_cpu, &
           int(ngrid_loc,kind=8))
#endif

      ! --- Step 1: Compute sendcounts ---
#ifdef FDMDEBUG
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif
      allocate(sendcounts(0:ncpu-1))
      allocate(recvcounts(0:ncpu-1))
      allocate(sdispls(0:ncpu-1))
      allocate(rdispls(0:ncpu-1))

      sendcounts = 0
      recvcounts = 0
      ngrid_loc = active(ilevel)%ngrid
      nthreads_fft = omp_get_max_threads()
      allocate(thread_counts(0:ncpu-1,0:nthreads_fft-1))
      allocate(thread_offsets(0:ncpu-1,0:nthreads_fft-1))
      thread_counts = 0
!$omp parallel num_threads(nthreads_fft) default(shared) &
!$omp private(tid,igrid,igrid_amr,Kx,ind,ix,dest_rank)
      tid = omp_get_thread_num()
!$omp do schedule(static)
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
         ! Half of a grid's children have each of the two x offsets.  Count
         ! those two destinations directly instead of repeating the same
         ! rank lookup for every y/z child.
         ix = modulo(Kx - 1, fft_Nx)
         dest_rank = min(ix / fftw_block, ncpu-1)
         thread_counts(dest_rank,tid) = thread_counts(dest_rank,tid) + twotondim/2
         ix = modulo(Kx, fft_Nx)
         dest_rank = min(ix / fftw_block, ncpu-1)
         thread_counts(dest_rank,tid) = thread_counts(dest_rank,tid) + twotondim/2
      end do
!$omp end do
!$omp end parallel
      do tid = 0, nthreads_fft-1
         sendcounts = sendcounts + thread_counts(:,tid)
      end do
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_COUNT,fftp_phase_wall,fftp_phase_cpu, &
           int(ngrid_loc,kind=8)*int(twotondim,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 2: Exchange counts via sparse P2P ---
      nreq = 0
      allocate(reqs(n_fftw_send + n_fftw_recv))

      do ip = 1, n_fftw_recv
         nreq = nreq + 1
         call MPI_IRECV(recvcounts(fftw_recv_partners(ip)), 1, MPI_INTEGER, &
              fftw_recv_partners(ip), 701, MPI_COMM_WORLD, reqs(nreq), info)
      end do
      do ip = 1, n_fftw_send
         nreq = nreq + 1
         call MPI_ISEND(sendcounts(fftw_send_partners(ip)), 1, MPI_INTEGER, &
              fftw_send_partners(ip), 701, MPI_COMM_WORLD, reqs(nreq), info)
      end do

      if(nreq > 0) then
         allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
         call MPI_WAITALL(nreq, reqs, mpi_stat, info)
         deallocate(mpi_stat)
      end if
      deallocate(reqs)

      ! Compute displacements
      sdispls(0) = 0
      rdispls(0) = 0
      do irank = 1, ncpu-1
         sdispls(irank) = sdispls(irank-1) + sendcounts(irank-1)
         rdispls(irank) = rdispls(irank-1) + recvcounts(irank-1)
      end do
      n_send_total = sdispls(ncpu-1) + sendcounts(ncpu-1)
      n_recv_total = rdispls(ncpu-1) + recvcounts(ncpu-1)
      do irank = 0, ncpu-1
         thread_offsets(irank,0) = sdispls(irank)
         do tid = 1, nthreads_fft-1
            thread_offsets(irank,tid) = thread_offsets(irank,tid-1) + &
                 thread_counts(irank,tid-1)
         end do
      end do
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_COUNTCOMM,fftp_phase_wall,fftp_phase_cpu, &
           int(n_fftw_send+n_fftw_recv,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 3: Pack send buffer: (slab_idx, rhs_value) pairs ---
      ! RHS uses un-boosted fourpi (corrections applied in Green's function)
      allocate(sendbuf(0:2*n_send_total-1))
      allocate(recvbuf(0:2*n_recv_total-1))

!$omp parallel num_threads(nthreads_fft) default(shared) &
!$omp private(tid,igrid,igrid_amr,Kx,Ky,Kz,ind,ix,iy,iz,dest_rank, &
!$omp x_local,idx_3d,icell_amr,slot)
      tid = omp_get_thread_num()
!$omp do schedule(static)
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
         Ky = nint(xg(igrid_amr, 2) * dble(fft_Ny))
         Kz = nint(xg(igrid_amr, 3) * dble(fft_Nz))
         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            iy = Ky - 1 + mod((ind-1)/2, 2)
            iz = Kz - 1 + (ind-1)/4
            ix = modulo(ix, fft_Nx)
            iy = modulo(iy, fft_Ny)
            iz = modulo(iz, fft_Nz)

            dest_rank = min(ix / fftw_block, ncpu-1)

            ! Local X index within destination slab
            x_local = ix - dest_rank * fftw_block

            ! Slab-local row-major index (in-place R2C: padded last dim)
            idx_3d = x_local * fft_Ny * 2*(fft_Nz/2+1) + iy * 2*(fft_Nz/2+1) + iz

            icell_amr = icell_of(igrid_amr,ind)

            slot = thread_offsets(dest_rank,tid)
            sendbuf(2*slot) = dble(idx_3d)
            if(solve_mode==0) then
               sendbuf(2*slot+1)=fourpi_fft*(rho(icell_amr)-rho_tot)
            else
               sendbuf(2*slot+1)=scalar_gr_old(icell_amr)
            end if
            thread_offsets(dest_rank,tid) = slot + 1
         end do
      end do
!$omp end do
!$omp end parallel
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_PACK,fftp_phase_wall,fftp_phase_cpu, &
           int(n_send_total,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 4: Forward data exchange via sparse P2P ---
      sendcounts = sendcounts * 2
      recvcounts = recvcounts * 2
      sdispls = sdispls * 2
      rdispls = rdispls * 2

      nreq = 0
      allocate(reqs(n_fftw_send + n_fftw_recv))

      do ip = 1, n_fftw_recv
         irank = fftw_recv_partners(ip)
         if(recvcounts(irank) > 0) then
            nreq = nreq + 1
            call MPI_IRECV(recvbuf(rdispls(irank)), recvcounts(irank), &
                 MPI_DOUBLE_PRECISION, irank, 702, MPI_COMM_WORLD, reqs(nreq), info)
         end if
      end do
      do ip = 1, n_fftw_send
         irank = fftw_send_partners(ip)
         if(sendcounts(irank) > 0) then
            nreq = nreq + 1
            call MPI_ISEND(sendbuf(sdispls(irank)), sendcounts(irank), &
                 MPI_DOUBLE_PRECISION, irank, 702, MPI_COMM_WORLD, reqs(nreq), info)
         end if
      end do

      if(nreq > 0) then
         allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
         call MPI_WAITALL(nreq, reqs, mpi_stat, info)
         deallocate(mpi_stat)
      end if
      deallocate(reqs)

      ! Unpack into FFTW real data array
!$omp parallel do schedule(static)
      do ix = 1, slab_real_size
         rdata_fftw(ix) = 0.0d0
      end do
!$omp end parallel do
      ! Preserve the serial accumulation order: reception/boundary entries
      ! can share a slab index even on a fully uniform solve.
      do ix = 0, n_recv_total - 1
         idx_3d = nint(recvbuf(2*ix))
         rdata_fftw(idx_3d + 1) = rdata_fftw(idx_3d + 1) + recvbuf(2*ix + 1)
      end do
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_FWDCOMM,fftp_phase_wall,fftp_phase_cpu, &
           int(n_recv_total,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 5: Forward R2C FFT ---
      call fftw_mpi_execute_dft_r2c(plan_r2c, rdata_fftw, cdata_fftw)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_R2C,fftp_phase_wall,fftp_phase_cpu,N_total)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Precompute DE kspace params if needed (quasi-static fallback) ---
      if(de_perturb .and. .not. de_table_loaded .and. de_helmholtz_on()) then
         call compute_de_kspace_params(aexp, de_kappa2, de_alpha)
      end if
      ! Precompute DE table params
      if(de_perturb .and. de_table_loaded) then
         omega_de_a = omega_l * f_de_val(aexp) * aexp**3
         if(use_neutrino .and. omega_nu > 0.0d0) then
            omega_cb_loc = omega_m - omega_nu
         else
            omega_cb_loc = omega_m
         end if
      end if

      ! --- Step 6: Apply Green's function (+ neutrino/DE corrections) ---
      ! TRANSPOSED complex layout (Ny-distributed):
      !   idx_c = ky_local*(Nx*(Nz/2+1)) + kx*(Nz/2+1) + kz, ky = ny_start+ky_local
      N_complex = int(local_ny_fftw,i8b) * int(fft_Nx,i8b) * int(fft_Nz/2+1,i8b)
      if(solve_mode==1) then
!$omp parallel do private(ky_i,kx_i,kz_i,denom,scalar_green) schedule(static)
         do idx_c=0,N_complex-1
            ky_i=int(ny_start_fftw)+int(idx_c/(int(fft_Nx,i8b)*int(fft_Nz/2+1,i8b)))
            kx_i=int(mod(idx_c/int(fft_Nz/2+1,i8b),int(fft_Nx,i8b)))
            kz_i=int(mod(idx_c,int(fft_Nz/2+1,i8b)))
            denom=2d0*(cos(twopi*dble(kx_i)/dble(fft_Nx)) &
                 +cos(twopi*dble(ky_i)/dble(fft_Ny)) &
                 +cos(twopi*dble(kz_i)/dble(fft_Nz))-3d0)-m2*dx2_fft
            if(abs(denom)>1d-30) then
               scalar_green=dx2_fft/denom
            else
               scalar_green=0d0
            end if
            cdata_fftw(idx_c+1)=cdata_fftw(idx_c+1)*scalar_green
         end do
!$omp end parallel do
      else if((use_neutrino .and. omega_nu > 0.0d0) .or. de_perturb &
           & .or. (use_horndeski .and. hs_mass > 0.0d0)) then
         mu_a_hs = horndeski_mu_of_a(aexp)
         if(use_neutrino .and. omega_nu > 0.0d0) omega_cb_loc = omega_m - omega_nu
!$omp parallel do private(ky_i,kx_i,kz_i,nu_factor,de_factor,k_phys, &
!$omp R_nu_val,R_DE_val,k_tilde_sq,hs_factor) schedule(static)
         do idx_c = 0, N_complex-1
            ! Compute integer wave numbers for this mode (transposed order)
            ky_i = int(ny_start_fftw) + int(idx_c / (int(fft_Nx,i8b)*int(fft_Nz/2+1,i8b)))
            kx_i = int(mod(idx_c / int(fft_Nz/2+1,i8b), int(fft_Nx,i8b)))
            kz_i = int(mod(idx_c, int(fft_Nz/2+1,i8b)))
            ! Fold to negative frequencies
            if(kx_i > fft_Nx/2) kx_i = kx_i - fft_Nx
            if(ky_i > fft_Ny/2) ky_i = ky_i - fft_Ny

            ! Combined correction factors (default 1.0 = no change)
            nu_factor = 1.0d0
            de_factor = 1.0d0

            ! Neutrino correction
            if(use_neutrino .and. omega_nu > 0.0d0) then
               k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
               if(k_phys > 0.0d0) then
                  R_nu_val = get_nu_ratio(k_phys, aexp)
                  nu_factor = 1.0d0 + (omega_nu / omega_cb_loc) * R_nu_val
               end if
            end if

            ! DE perturbation correction
            if(de_perturb) then
               if(de_table_loaded) then
                  ! Table-based linear response: 1 + (rho_DE/rho_cb) * R_DE(k,a)
                  k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
                  if(k_phys > 0.0d0) then
                     R_DE_val = get_de_ratio(k_phys, aexp)
                     de_factor = 1.0d0 + (omega_de_a / omega_cb_loc) * R_DE_val
                  end if
               else if(de_helmholtz_on()) then
                  ! Fallback: kappa2/alpha quasi-static method
                  k_tilde_sq = twopi * twopi * dble(kx_i**2 + ky_i**2 + kz_i**2)
                  de_factor = 1.0d0 + de_alpha / (k_tilde_sq + de_kappa2)
               end if
            end if

            ! Horndeski k-dependent mu(a,k)
            hs_factor = 1.0d0
            if(use_horndeski .and. hs_mass > 0.0d0) then
               k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
               hs_factor = 1.0d0 + (mu_a_hs - 1.0d0) * k_phys*k_phys / &
                    & (k_phys*k_phys + (aexp*hs_mass)**2)
            end if

            cdata_fftw(idx_c + 1) = cdata_fftw(idx_c + 1) * green(idx_c) * nu_factor * de_factor * hs_factor
         end do
!$omp end parallel do
      else
!$omp parallel do schedule(static)
         do idx_c = 0, N_complex-1
            cdata_fftw(idx_c + 1) = cdata_fftw(idx_c + 1) * green(idx_c)
         end do
!$omp end parallel do
      end if
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_GREEN,fftp_phase_wall,fftp_phase_cpu,N_complex)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 7: Inverse C2R FFT ---
      call fftw_mpi_execute_dft_c2r(plan_c2r, cdata_fftw, rdata_fftw)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_C2R,fftp_phase_wall,fftp_phase_cpu,N_total)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 8: Normalize ---
!$omp parallel do schedule(static)
      do ix = 1, slab_real_size
         rdata_fftw(ix) = rdata_fftw(ix) / dble(N_total)
      end do
!$omp end parallel do

      ! --- Step 9: Reverse transfer via sparse P2P ---
      ! Replace rhs values in recvbuf with phi from rdata_fftw
!$omp parallel do private(idx_3d) schedule(static)
      do ix = 0, n_recv_total - 1
         idx_3d = nint(recvbuf(2*ix))
         recvbuf(2*ix + 1) = rdata_fftw(idx_3d + 1)
      end do
!$omp end parallel do

      nreq = 0
      allocate(reqs(n_fftw_send + n_fftw_recv))

      ! Receive phi from FFTW slab owners (my send partners)
      do ip = 1, n_fftw_send
         irank = fftw_send_partners(ip)
         if(sendcounts(irank) > 0) then
            nreq = nreq + 1
            call MPI_IRECV(sendbuf(sdispls(irank)), sendcounts(irank), &
                 MPI_DOUBLE_PRECISION, irank, 703, MPI_COMM_WORLD, reqs(nreq), info)
         end if
      end do
      ! Send phi to original RAMSES senders (my recv partners)
      do ip = 1, n_fftw_recv
         irank = fftw_recv_partners(ip)
         if(recvcounts(irank) > 0) then
            nreq = nreq + 1
            call MPI_ISEND(recvbuf(rdispls(irank)), recvcounts(irank), &
                 MPI_DOUBLE_PRECISION, irank, 703, MPI_COMM_WORLD, reqs(nreq), info)
         end if
      end do

      if(nreq > 0) then
         allocate(mpi_stat(MPI_STATUS_SIZE, nreq))
         call MPI_WAITALL(nreq, reqs, mpi_stat, info)
         deallocate(mpi_stat)
      end if
      deallocate(reqs)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_REVCOMM,fftp_phase_wall,fftp_phase_cpu, &
           int(n_recv_total,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! Unpack phi to RAMSES cells
      do irank = 0, ncpu-1
         thread_offsets(irank,0) = sdispls(irank) / 2
         do tid = 1, nthreads_fft-1
            thread_offsets(irank,tid) = thread_offsets(irank,tid-1) + &
                 thread_counts(irank,tid-1)
         end do
      end do

!$omp parallel num_threads(nthreads_fft) default(shared) &
!$omp private(tid,igrid,igrid_amr,Kx,ind,ix,dest_rank, &
!$omp icell_amr,slot,uold,du)
      tid = omp_get_thread_num()
!$omp do schedule(static)
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            ix = modulo(ix, fft_Nx)
            dest_rank = min(ix / fftw_block, ncpu-1)

            icell_amr = icell_of(igrid_amr,ind)
            slot = thread_offsets(dest_rank,tid)

            if(solve_mode==0) then
               phi(icell_amr)=sendbuf(2*slot+1)
            else
               uold=scalar_gr(icell_amr)
               du=relax*sendbuf(2*slot+1)
               if(step_frac>0d0 .and. abs(uold)>0d0) &
                    du=max(-step_frac*abs(uold),min(step_frac*abs(uold),du))
               scalar_gr(icell_amr)=uold+du
            end if
            thread_offsets(dest_rank,tid) = slot + 1
         end do
      end do
!$omp end do
!$omp end parallel

      deallocate(sendcounts, sdispls, recvcounts, rdispls)
      deallocate(sendbuf, recvbuf, thread_counts, thread_offsets)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_SCATTER,fftp_phase_wall,fftp_phase_cpu, &
           int(ngrid_loc,kind=8)*int(twotondim,kind=8))
#endif

   else
      ! ============================================================
      ! SMALL GRID PATH: ALLREDUCE + local FFTW3+OMP
      ! ============================================================

      ! Persistent allocation (only on size change)
      if(N_total /= saved_N_total) then
         if(allocated(fft_map))   deallocate(fft_map)
         if(allocated(rhs_local)) deallocate(rhs_local)
         if(allocated(rhs_3d))    deallocate(rhs_3d)
         allocate(fft_map(0:N_total-1))
         allocate(rhs_local(0:N_total-1))
         allocate(rhs_3d(0:N_total-1))
         saved_N_total = N_total
      end if

      ! --- Step 1: Build fft_map and gather local RHS ---
#ifdef FDMDEBUG
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif
      fft_map   = 0
      rhs_local = 0.0d0

      ngrid_loc = active(ilevel)%ngrid
      do igrid = 1, ngrid_loc
         igrid_amr = active(ilevel)%igrid(igrid)
         Kx = nint(xg(igrid_amr, 1) * dble(fft_Nx))
         Ky = nint(xg(igrid_amr, 2) * dble(fft_Ny))
         Kz = nint(xg(igrid_amr, 3) * dble(fft_Nz))
         do ind = 1, twotondim
            ix = Kx - 1 + mod(ind-1, 2)
            iy = Ky - 1 + mod((ind-1)/2, 2)
            iz = Kz - 1 + (ind-1)/4
            ix = modulo(ix, fft_Nx)
            iy = modulo(iy, fft_Ny)
            iz = modulo(iz, fft_Nz)

            ! C row-major index
            idx_3d = ix * fft_Ny * fft_Nz + iy * fft_Nz + iz

            icell_amr = icell_of(igrid_amr,ind)

            fft_map(idx_3d) = icell_amr
            if(solve_mode==0) then
               ! Use un-boosted fourpi (corrections in Green's function)
               rhs_local(idx_3d)=fourpi_fft*(rho(icell_amr)-rho_tot)
            else
               rhs_local(idx_3d)=scalar_gr_old(icell_amr)
            end if
         end do
      end do
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_PACK,fftp_phase_wall,fftp_phase_cpu, &
           int(ngrid_loc,kind=8)*int(twotondim,kind=8))
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 2: MPI_ALLREDUCE ---
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(rhs_local, rhs_3d, int(N_total), &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
      rhs_3d = rhs_local
#endif
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_FWDCOMM,fftp_phase_wall,fftp_phase_cpu,N_total)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 3: Forward R2C FFT ---
      call fftw_execute_dft_r2c(plan_r2c, rhs_3d, cdata)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_R2C,fftp_phase_wall,fftp_phase_cpu,N_total)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Precompute DE kspace params if needed (quasi-static fallback) ---
      if(de_perturb .and. .not. de_table_loaded .and. de_helmholtz_on()) then
         call compute_de_kspace_params(aexp, de_kappa2, de_alpha)
      end if
      ! Precompute DE table params
      if(de_perturb .and. de_table_loaded) then
         omega_de_a = omega_l * f_de_val(aexp) * aexp**3
         if(use_neutrino .and. omega_nu > 0.0d0) then
            omega_cb_loc = omega_m - omega_nu
         else
            omega_cb_loc = omega_m
         end if
      end if

      ! --- Step 4: Apply Green's function (+ neutrino/DE corrections) ---
      N_complex = int(fft_Nx,i8b) * int(fft_Ny,i8b) * int(fft_Nz/2+1,i8b)
      if(solve_mode==1) then
         do idx_c=0,N_complex-1
            kx_i=int(idx_c/(int(fft_Ny,i8b)*int(fft_Nz/2+1,i8b)))
            ky_i=int(mod(idx_c/int(fft_Nz/2+1,i8b),int(fft_Ny,i8b)))
            kz_i=int(mod(idx_c,int(fft_Nz/2+1,i8b)))
            denom=2d0*(cos(twopi*dble(kx_i)/dble(fft_Nx)) &
                 +cos(twopi*dble(ky_i)/dble(fft_Ny)) &
                 +cos(twopi*dble(kz_i)/dble(fft_Nz))-3d0)-m2*dx2_fft
            if(abs(denom)>1d-30) then
               scalar_green=dx2_fft/denom
            else
               scalar_green=0d0
            end if
            cdata(idx_c)=cdata(idx_c)*scalar_green
         end do
      else if((use_neutrino .and. omega_nu > 0.0d0) .or. de_perturb &
           & .or. (use_horndeski .and. hs_mass > 0.0d0)) then
         mu_a_hs = horndeski_mu_of_a(aexp)
         if(use_neutrino .and. omega_nu > 0.0d0) omega_cb_loc = omega_m - omega_nu
         do idx_c = 0, N_complex-1
            kx_i = int(idx_c / (int(fft_Ny,i8b)*int(fft_Nz/2+1,i8b)))
            ky_i = int(mod(idx_c / int(fft_Nz/2+1,i8b), int(fft_Ny,i8b)))
            kz_i = int(mod(idx_c, int(fft_Nz/2+1,i8b)))
            if(kx_i > fft_Nx/2) kx_i = kx_i - fft_Nx
            if(ky_i > fft_Ny/2) ky_i = ky_i - fft_Ny

            ! Combined correction factors
            nu_factor = 1.0d0
            de_factor = 1.0d0

            ! Neutrino correction
            if(use_neutrino .and. omega_nu > 0.0d0) then
               k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
               if(k_phys > 0.0d0) then
                  R_nu_val = get_nu_ratio(k_phys, aexp)
                  nu_factor = 1.0d0 + (omega_nu / omega_cb_loc) * R_nu_val
               end if
            end if

            ! DE perturbation correction
            if(de_perturb) then
               if(de_table_loaded) then
                  ! Table-based linear response
                  k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
                  if(k_phys > 0.0d0) then
                     R_DE_val = get_de_ratio(k_phys, aexp)
                     de_factor = 1.0d0 + (omega_de_a / omega_cb_loc) * R_DE_val
                  end if
               else if(de_helmholtz_on()) then
                  ! Fallback: kappa2/alpha quasi-static
                  k_tilde_sq = twopi * twopi * dble(kx_i**2 + ky_i**2 + kz_i**2)
                  de_factor = 1.0d0 + de_alpha / (k_tilde_sq + de_kappa2)
               end if
            end if

            ! Horndeski k-dependent mu(a,k)
            hs_factor = 1.0d0
            if(use_horndeski .and. hs_mass > 0.0d0) then
               k_phys = twopi * sqrt(dble(kx_i**2 + ky_i**2 + kz_i**2)) / boxlen_ini
               hs_factor = 1.0d0 + (mu_a_hs - 1.0d0) * k_phys*k_phys / &
                    & (k_phys*k_phys + (aexp*hs_mass)**2)
            end if

            cdata(idx_c) = cdata(idx_c) * green(idx_c) * nu_factor * de_factor * hs_factor
         end do
      else
         do idx_c = 0, N_complex-1
            cdata(idx_c) = cdata(idx_c) * green(idx_c)
         end do
      end if
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_GREEN,fftp_phase_wall,fftp_phase_cpu,N_complex)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 5: Inverse C2R FFT ---
      call fftw_execute_dft_c2r(plan_c2r, cdata, rhs_3d)
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_C2R,fftp_phase_wall,fftp_phase_cpu,N_total)
      call fftp_start(fftp_phase_wall,fftp_phase_cpu)
#endif

      ! --- Step 6: Normalize ---
      rhs_3d = rhs_3d / dble(N_total)

      ! --- Step 7: Scatter to phi via fft_map ---
      do idx_3d = 0, int(N_total-1)
         icell_amr = fft_map(idx_3d)
         if(icell_amr>0) then
            if(solve_mode==0) then
               phi(icell_amr)=rhs_3d(idx_3d)
            else
               uold=scalar_gr(icell_amr)
               du=relax*rhs_3d(idx_3d)
               if(step_frac>0d0 .and. abs(uold)>0d0) &
                    du=max(-step_frac*abs(uold),min(step_frac*abs(uold),du))
               scalar_gr(icell_amr)=uold+du
            end if
         end if
      end do
#ifdef FDMDEBUG
      call fftp_stop(fftp_path,FFTP_SCATTER,fftp_phase_wall,fftp_phase_cpu,N_total)
#endif

   end if  ! use_distributed

#ifdef FDMDEBUG
   call fftp_stop(fftp_path,FFTP_TOTAL,fftp_total_wall,fftp_total_cpu,N_total)
   call fftp_report(nstep_coarse,myid,ncpu,omp_get_max_threads())
#endif

end subroutine fftw_uniform_solve_engine
#endif

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use pbh_commons, only: use_pbh, pbh_mark_level
  use omp_lib, only: omp_get_wtime,omp_get_max_threads
#ifdef HYDRO_CUDA
  use cuda_commons, only: cuda_pool_is_initialized_c
  use poisson_cuda_interface, only: cuda_mg_release_arrays_c
  use hydro_cuda_interface, only: cuda_mesh_free_c
#endif
#ifdef SNRT
  use snrt_ramses_driver, only: snrt_ramses_diagnose_level, &
       snrt_ramses_advance_level
#endif
#ifdef SNRT_LEDGER_DIAGNOSTIC
  use snrt_agn_ledger, only: snrt_agn_ledger_diagnose
#endif
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar,mpi_err,jgrid,igrid_lb,ind_lb
  integer(kind=8)::nleaf_lb
  logical::ok_defrag,output_now_all,lb_timing_sample
  logical::phi_topology_changed,phi_topology_changed_all
!!$  integer::i,idim,ivar
!!$  logical::ok_defrag
  logical,save::first_step=.true.
  integer:: info

  real(kind=4):: real_mem, real_mem_tot
  real(kind=8):: t_lb_level_start, t_lb_level_end
  real(kind=8):: t_lb_cpu_level_start,t_lb_cpu_level_end
  real(kind=8):: t_lb_cpu_child_start,t_lb_cpu_child_time

  ! Particle sub-timers
  integer(kind=8) :: pt_t1, pt_t2, pt_rate
  real(dp), save :: pt_mktree=0, pt_killtree=0, pt_synchro=0, pt_move=0, pt_merge=0

  ! Sink sub-timers
  integer(kind=8) :: sk_t1, sk_t2
  real(dp), save :: sk_agn_fb=0, sk_create_sink=0, sk_grow=0, sk_bondi_hoyle=0

  ! Re-entrant-safe direct cooling timings.  The legacy phase timer is a
  ! single global state and is not reliable across recursive amr_step calls.
  real(dp), save :: cool_fine_wall=0d0, cool_table_wall=0d0
  real(dp) :: cool_t1

  ! SNRT is called immediately after cooling without a legacy phase switch.
  real(dp), save :: snrt_advance_wall=0d0, snrt_diagnose_wall=0d0

#ifdef HYDRO_CUDA
  ! GPU auto-tuning framework
  ! Phase 0: first call → force CPU, record time
  ! Phase 1: second call → force GPU, record time
  ! Phase 2+: use faster path, keep booking times
  !
  ! AGN feedback (sink particles)
  integer, save :: sk_auto_phase = 0
  real(dp), save :: sk_cpu_ref = 0d0, sk_gpu_ref = 0d0
  logical, save :: sk_use_gpu = .false.
  logical, save :: sk_auto_init = .false.
  real(dp) :: sk_dt_agn
  !
  ! Hydro Godunov
  integer, save :: hy_auto_phase = 0
  real(dp), save :: hy_cpu_ref = 0d0, hy_gpu_ref = 0d0
  logical, save :: hy_use_gpu = .false.
  logical, save :: hy_auto_init = .false.
  integer(kind=8) :: hy_t1, hy_t2
  real(dp) :: hy_dt
  !
  ! Poisson MG (gpu_poisson only; gpu_fft excluded from auto-tuning)
  integer, save :: mg_auto_phase = 0
  real(dp), save :: mg_cpu_ref = 0d0, mg_gpu_ref = 0d0
  logical, save :: mg_use_gpu = .false.
  logical, save :: mg_auto_init = .false.
  logical, save :: mg_orig_poisson = .false.
  integer(kind=8) :: mg_t1, mg_t2
  real(dp) :: mg_dt
#endif

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel, levelmin

#ifdef HYDRO_CUDA
  !---------------------------------------------------
  ! GPU auto-tuning: set flags at start of coarse step
  !---------------------------------------------------
  if(ilevel==levelmin) then
     ! Hydro auto-tuning: init + flag setting
     if(gpu_auto_tune) then
        if(.not. hy_auto_init .and. gpu_hydro) then
           hy_auto_init = .true.
           hy_auto_phase = 0
        endif
        if(hy_auto_init) then
           if(hy_auto_phase == 0) then
              gpu_hydro = .false.
           else if(hy_auto_phase == 1) then
              gpu_hydro = .true.
           else
              gpu_hydro = hy_use_gpu
           endif
        endif
     endif
     ! Poisson auto-tuning: init is done near multigrid_fine call
     ! (mg_auto_init set there), flags set at levelmin entry
  endif
#endif

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
                               call timer('refine','start')
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
                 end do
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
#ifdef FDMDEBUG
           if(use_fdm .and. i>=levelmin) call fdm_mass_check('pre-refine',i)
#endif
           call refine_fine(i)
           if(poisson .and. i<nlevelmax .and. &
                allocated(phi_checkpoint_level_valid))then
              phi_topology_changed=ncreate>0 .or. nkill>0
#ifndef WITHOUTMPI
              call MPI_ALLREDUCE(phi_topology_changed,phi_topology_changed_all, &
                   1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_err)
#else
              phi_topology_changed_all=phi_topology_changed
#endif
              if(phi_topology_changed_all)then
                 phi_checkpoint_level_valid(i+1:nlevelmax)=.false.
                 if(allocated(phi_restart_available)) &
                      phi_restart_available(i+1:nlevelmax)=.false.
              end if
           end if
#ifdef FDMDEBUG
           if(use_fdm .and. i>=levelmin) call fdm_mass_check('post-refine',i)
#endif
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
                               call timer('loadbalance','start')
  ok_defrag=.false.
  ! Variable-ncpu restart: must run even when levelmin==nlevelmax
  if(ilevel==levelmin .and. varcpu_restart_done)then
     if(myid==1) write(*,*) 'Forcing load_balance after variable-ncpu restart'
     call load_balance
     ! Do not defragment the cross-ordering/variable-ncpu intermediate
     ! hierarchy here.  load_balance has already rebuilt valid sparse-grid
     ! communicators, whereas defrag assumes that every remote parent/son
     ! link is locally complete.  Normal evolution or the next output can
     ! compact the now-consistent hierarchy later.
     do i=nlevelmax,1,-1
        if(hydro)then
           do ivar=1,nvar
              call make_virtual_fine_dp(uold(1,ivar),i)
           end do
        end if
        if(poisson .and. allocated(f))then
           call make_virtual_fine_dp(phi(1),i)
           do idim=1,ndim
              call make_virtual_fine_dp(f(1,idim),i)
           end do
           if(allocated(scalar_gr))then
              call make_virtual_fine_dp(scalar_gr(1),i)
              call make_virtual_fine_dp(scalar_gr_old(1),i)
           end if
        end if
     end do
     if(use_fdm) then
        call morton_hash_rebuild()
        call restore_psi_postlb()
     end if
     if(allocated(varcpu_nactive)) deallocate(varcpu_nactive)
     if(allocated(varcpu_my_files)) deallocate(varcpu_my_files)
     if(allocated(varcpu_ngrid_file)) deallocate(varcpu_ngrid_file)
     varcpu_restart_done=.false.
     first_step=.false.
     call diag_check_nan('post_varcpu_load_balance')
     if(myid==1) then
        write(*,*) 'Variable-ncpu restart block done, entering time step'
        call flush(6)
     end if
  end if
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(lb_force_remap)then
           if(myid==1) write(*,'(A,I0)') &
                ' Runtime load balance: forcing remap at coarse step ',nstep_coarse
           call load_balance
           call defrag
           lb_force_remap=.false.
           first_step=.false.
           ok_defrag=.true.
        else if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        else if(remap_thresh>0d0)then
           ! Auto remap: check weight inhomogeneity every coarse step
           ! Skip first step on restart (already balanced before dump)
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              call check_load_imbalance(ok_defrag)
           end if
        end if
     endif
  end if

  ! Sparse measurement starts after any coarse-level remap.  It is outside
  ! the hydro block so DMO, SIDM, and modified-gravity runs are all sampled.
  lb_timing_sample=(.not.memory_balance).and.lb_timing_interval>0
  if(lb_timing_sample) lb_timing_sample= &
       mod(nstep_coarse,lb_timing_interval)==0
  if(lb_timing_sample)then
     t_lb_level_start=omp_get_wtime()
     call cpu_time(t_lb_cpu_level_start)
  endif
  t_lb_cpu_child_time=0d0

  !-----------------
  ! Particle leakage
  !-----------------
                               call timer('particles','start')
  call system_clock(pt_t1, pt_rate)
  if(pic)call make_tree_fine(ilevel)
  call system_clock(pt_t2)
  pt_mktree = pt_mktree + dble(pt_t2-pt_t1)/dble(pt_rate)
  if(ilevel==levelmin) call diag_check_nan('post_maketree')
  
  !------------------------
  ! Output results to files
  !------------------------
   if(ilevel==levelmin) then
   endif


  if(ilevel==levelmin)then
     ! check if any of the processes received a signal for output
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(output_now,output_now_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_err)




     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now_all.EQV..true.)then
                               call timer('io','start')
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)

        ! Dump lightcone
!jhshin1
!        if(lightcone) call output_cone()
!jhshin2
        if (output_now_all.EQV..true.) then
          output_now=.false.
        endif

     endif

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then ! ifort returns error for next statement if looking
                             ! beyond what is allocated as an array amovout/tmovout
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
                               call timer('io','start')
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then

     !jaehyun

        if(aexp_old2 .le. 0.0001) aexp_old2=aexp_old
	 
	 if(spherical_region) then
        do i=1,5
	   if(hydro)   call output_sphere_hydro(i)
	   call output_sphere_part(i)
	   if(sink)   call output_sphere_sink(i)
        enddo
     endif

	 if(lightcone) then
        do i=1,2 
           if(hydro) call output_cone_hydro(i)
           call output_cone_part(i)
           if(sink) call output_cone_sink(i)
        enddo
	 endif
	 aexp_old2=aexp

     !jaehyun

                               call timer('feedback','start')
        !----------------------------------------------------
        ! Kinetic feedback
        !----------------------------------------------------
     if(hydro.and.star.and.f_w>0.)call kinetic_feedback
     if(hydro.and.star.and.f_w>0.)call diag_check_eint('kinetic_fb',0)
     call diag_check_nan('post_kinfb')

     call timer('sinks','start')
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning: set gpu_sink for this step ---
     if(.not. sk_auto_init .and. gpu_sink) then
        sk_auto_init = .true.
        sk_auto_phase = 0
     endif
     if(sk_auto_init) then
        if(sk_auto_phase == 0) then
           gpu_sink = .false.  ! Phase 0: force CPU path
        else if(sk_auto_phase == 1) then
           gpu_sink = .true.   ! Phase 1: force GPU path
        else
           gpu_sink = sk_use_gpu  ! Phase 2+: use decided path
        endif
     endif
#endif
     call system_clock(sk_t1)
     if(sink .and. sink_AGN)call AGN_feedback
     if(sink .and. sink_AGN)call diag_check_eint('AGN_fb',0)
     call diag_check_nan('post_agnfb')
     call system_clock(sk_t2)
     sk_agn_fb = sk_agn_fb + dble(sk_t2-sk_t1)/dble(pt_rate)
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning: record time and update decision ---
     if(sk_auto_init .and. sink .and. sink_AGN) then
        sk_dt_agn = dble(sk_t2-sk_t1)/dble(pt_rate)
        if(sk_auto_phase == 0) then
           ! Phase 0 done: recorded CPU time
           sk_cpu_ref = sk_dt_agn
           sk_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] AGN_feedback CPU: ', sk_cpu_ref, ' s'
        else if(sk_auto_phase == 1) then
           ! Phase 1 done: recorded GPU time, decide
           sk_gpu_ref = sk_dt_agn
           sk_use_gpu = (sk_gpu_ref < sk_cpu_ref)
           sk_auto_phase = 2
           gpu_sink = sk_use_gpu
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] AGN_feedback GPU: ', sk_gpu_ref, &
                   ' s  (CPU was ', sk_cpu_ref, ' s)'
              if(sk_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: keep booking, check for switch
           if(sk_use_gpu) then
              sk_gpu_ref = sk_dt_agn
           else
              sk_cpu_ref = sk_dt_agn
              ! Switch to GPU if CPU became slower
              if(sk_cpu_ref > sk_gpu_ref .and. sk_gpu_ref > 0d0) then
                 sk_use_gpu = .true.
                 gpu_sink = .true.
                 if(myid==1) write(*,'(A,F8.3,A,F8.3,A)') &
                      ' [GPU auto-tune] Switching to GPU: CPU=', &
                      sk_cpu_ref, ' s > GPU=', sk_gpu_ref, ' s'
              endif
           endif
        endif
     endif
#endif
     !-----------------------------------------------------
     ! Create sink particles and associated cloud particles
     !-----------------------------------------------------
     call system_clock(sk_t1)
     if(sink)call create_sink
     if(sink)call diag_check_eint('create_sink',0)
     call system_clock(sk_t2)
     sk_create_sink = sk_create_sink + dble(sk_t2-sk_t1)/dble(pt_rate)

     !-----------------------------------------------------
     ! Enforce eEOS floor after sink/AGN operations
     ! Prevents negative internal energy AND extreme velocity
     ! from AGN jet momentum injection (energy-conservative)
     !-----------------------------------------------------
     if(hydro .and. eeos_poly_coeff > 0d0)then
        call enforce_eeos_after_sink
     endif

  endif

  if(ilevel==levelmin) call diag_check_nan('pre_rho')

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('poisson','start')
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     if(timer_report_interval>0) call timer('rho','start')
     call rho_fine(ilevel,icount)
     if(ilevel==levelmin) call diag_check_nan('post_rho')
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
                               call timer('particles','start')
     call system_clock(pt_t1)
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call virtual_tree_fine(ilevel)
     call system_clock(pt_t2)
     pt_killtree = pt_killtree + dble(pt_t2-pt_t1)/dble(pt_rate)
  end if
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
     call timer('poisson','start')

     if(ilevel==levelmin) call diag_check_nan('pre_sync1')

     ! Remove gravity source term with half time step and old force
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

     if(ilevel==levelmin) call diag_check_nan('post_sync1')
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     ! Compute gravitational potential
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
!jhshin1
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning for Poisson MG (NOT fft) ---
     ! gpu_fft is excluded from auto-tuning because cuFFT direct solve
     ! and MG V-cycle produce different potential scales, so switching
     ! between them mid-run causes energy conservation failure.
     if(gpu_auto_tune) then
        if(.not. mg_auto_init .and. gpu_poisson) then
           mg_auto_init = .true.
           mg_auto_phase = 0
           mg_orig_poisson = gpu_poisson
        endif
        if(mg_auto_init .and. ilevel==levelmin) then
           if(mg_auto_phase == 0) then
              gpu_poisson = .false.
           else if(mg_auto_phase == 1) then
              gpu_poisson = mg_orig_poisson
           else
              if(.not. mg_use_gpu) then
                 gpu_poisson = .false.
              endif
           endif
        endif
     endif
#endif
#ifdef HYDRO_CUDA
     call system_clock(mg_t1)
#endif
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call timer('poisson - cg', 'start')
           call phi_fine_cg(ilevel,icount)
        else
           call timer('poisson - mg AMR', 'start')
           call multigrid_fine(ilevel,icount)
        end if
     else
#ifdef USE_FFTW
        if(use_fftw) then
           call timer('poisson-fftw3 base','start')
        else
#endif
#ifdef HYDRO_CUDA
        if(gpu_fft .and. cuda_pool_is_initialized_c()/=0) then
           call timer('poisson-cuFFT base','start')
        else
#endif
           call timer('poisson - mg base', 'start')
#ifdef HYDRO_CUDA
        end if
#endif
#ifdef USE_FFTW
        end if
#endif
        call multigrid_fine(levelmin,icount)
     end if
     if(allocated(phi_checkpoint_level_valid)) &
          phi_checkpoint_level_valid(ilevel)=.true.
     ! FFT direct solves bypass the generic initial-guess block, so consume
     ! the one-shot restored-phi flag here for every solver backend.
     if(allocated(phi_restart_available)) &
          phi_restart_available(ilevel)=.false.
#ifdef HYDRO_CUDA
     call system_clock(mg_t2)
     if(mg_auto_init) then
        mg_dt = dble(mg_t2-mg_t1)/dble(pt_rate)
        if(mg_auto_phase == 0) then
           mg_cpu_ref = mg_cpu_ref + mg_dt
        else if(mg_auto_phase == 1) then
           mg_gpu_ref = mg_gpu_ref + mg_dt
        endif
     endif
#endif
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################
      call timer('poisson', 'start')
     !when there is no old potential...
!jhshin2      
     
     if (nstep==0)call save_phi_old(ilevel)
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     ! Compute gravitational acceleration
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call force_fine(ilevel,icount)
     if(ilevel==levelmin) call diag_check_nan('post_force')
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     !-------------------------------------------------
     ! MOND: phantom density + Poisson re-solve
     !-------------------------------------------------
     if(use_mond .and. mond_type >= 1) then
        call timer('mond-phantom','start')

        if(mond_type == 1) then
           ! Phase 1: QUMOND single-pass
           call compute_mond_phantom_density(ilevel, .false.)
           call make_virtual_fine_dp(rho(1),ilevel)

           if(ilevel>levelmin)then
              if(ilevel .ge. cg_levelmin) then
                 call phi_fine_cg(ilevel,icount)
              else
                 call multigrid_fine(ilevel,icount)
              end if
           else
              call multigrid_fine(levelmin,icount)
           end if

           call force_fine(ilevel,icount)

        else if(mond_type == 2) then
           ! Phase 2: AQUAL iterative
           call aqual_iterate(ilevel, icount)
        end if

        call timer('poisson','start')
     end if

     !-------------------------------------------------
     ! f(R) Hu-Sawicki gravity
     !-------------------------------------------------
     if(use_fR) then
        call timer('fR-solve','start')
        call fR_solve_level(ilevel, icount)
        call timer('poisson','start')
     end if

     !-------------------------------------------------
     ! nDGP gravity
     !-------------------------------------------------
     if(use_nDGP) then
        call timer('nDGP-solve','start')
        call nDGP_solve_level(ilevel, icount)
        call timer('poisson','start')
     end if

     !-------------------------------------------------
     ! Symmetron gravity
     !-------------------------------------------------
     if(use_symmetron) then
        call timer('symmetron','start')
        call symmetron_solve_level(ilevel, icount)
        call timer('poisson','start')
     end if

     !-------------------------------------------------
     ! Dilaton gravity
     !-------------------------------------------------
     if(use_dilaton) then
        call timer('dilaton','start')
        call dilaton_solve_level(ilevel, icount)
        call timer('poisson','start')
     end if

     !-------------------------------------------------
     ! Galileon (cubic) gravity
     !-------------------------------------------------
     if(use_galileon) then
        call timer('galileon','start')
        call galileon_solve_level(ilevel, icount)
        call timer('poisson','start')
     end if

     ! ADM hydro-particle-mesh closure. It augments the cell force before
     ! the ordinary KDK particle kicks.
     if(use_adm .and. adm_hpm) call adm_hpm_force_fine(ilevel)

     ! Synchronize remaining particles for gravity
     if(pic)then
                               call timer('particles','start')
        call system_clock(pt_t1)
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
        call synchro_fine(ilevel)
        call system_clock(pt_t2)
        pt_synchro = pt_synchro + dble(pt_t2-pt_t1)/dble(pt_rate)

        ! FDM Schrödinger-Poisson step (after Poisson solve, using latest Phi)
        if(use_fdm) call fdm_step(ilevel)

        ! SIDM scattering (after velocity sync, before timestep)
        if(sidm) call sidm_scatter(ilevel)

        ! DM-baryon drag force (after DM-DM scatter)
        if(sidm_baryon) call sidm_baryon_drag(ilevel)

        ! Atomic Dark Matter: dark-sector cooling
        if(use_adm) call dark_cooling_fine(ilevel)

     end if
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

     if(hydro)then
                               call timer('poisson','start')

!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        if(ilevel==levelmin) call diag_check_nan('post_sync2')
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################


        ! Density threshold and/or Bondi accretion onto sink particle
                               call timer('sinks','start')
        call system_clock(sk_t1)
        if(sink)then
           if(bondi)then
              call grow_bondi(ilevel)
           else
              call grow_jeans(ilevel)
           endif
           call diag_check_eint('grow_bondi',ilevel)
        endif
        call system_clock(sk_t2)
        sk_grow = sk_grow + dble(sk_t2-sk_t1)/dble(pt_rate)
!############################################################
!       call getmem(real_mem)
!       call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
!       if(myid==1) then
!           call writemem(real_mem_tot)
!       endif
!############################################################

        ! Update boundaries
                               call timer('hydro - ghostzones','start')
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(simple_boundary)call make_boundary_hydro(ilevel)

     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
                               call timer('radiative transfer','start')
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
                               call timer('radiative transfer','start')
  if(rt)call rt_set_unew(ilevel)
#endif

  ! Record this level step's starting scale factor for the PBH update
  ! (aexp still holds the step-start value here; after the recursion
  ! update_time has advanced it to the step end)
  if(pic.and.use_pbh)call pbh_mark_level(ilevel,nlevelmax,aexp)

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(lb_timing_sample)then
     call cpu_time(t_lb_cpu_child_start)
  endif
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
        ! Restrict FDM psi/(rho,S) from the completed finer level so
        ! non-leaf cells stay fresh (derefinement resumes from restricted
        ! values; wave->fluid parents get rho=<|psi|^2>, unwrapped S).
        if(use_fdm) then
#ifdef FDMDEBUG
           call fdm_mass_check('pre-restrict',ilevel)
#endif
           call fdm_restrict(ilevel)
#ifdef FDMDEBUG
           call fdm_mass_check('post-restrict',ilevel)
#endif
        end if
     else
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
     end if
  else
     call update_time(ilevel)
  end if
  if(lb_timing_sample)then
     call cpu_time(t_lb_cpu_level_end)
     t_lb_cpu_child_time=t_lb_cpu_level_end-t_lb_cpu_child_start
  endif

  ! Thermal feedback from stars (also call if no feedback, for bookkeeping)
  if(hydro.and.star) then
                               call timer('feedback','start')
     call thermal_feedback(ilevel)
  endif

  ! Evaporating-PBH update: exact mass loss + local heating into unew
  ! (must stay inside the set_unew/set_uold window, see paper appendix A)
  if(pic.and.use_pbh) then
                               call timer('pbh evap','start')
     call pbh_evap_fine(ilevel)
  endif

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
#ifdef HYDRO_CUDA
     ! --- GPU auto-tuning for hydro: flags set at levelmin entry ---
     ! (actual flag setting is done at start of amr_step(levelmin))
#endif
#ifdef HYDRO_CUDA
     call system_clock(hy_t1)
     ! Release MG Poisson GPU arrays before hydro mesh allocation
     if(gpu_hydro) call cuda_mg_release_arrays_c()
#endif
     call godunov_fine(ilevel)
#ifdef HYDRO_CUDA
     ! Free hydro mesh from GPU after godunov_fine
     if(gpu_hydro) call cuda_mesh_free_c()
     call system_clock(hy_t2)
     if(hy_auto_init) then
        hy_dt = dble(hy_t2-hy_t1)/dble(pt_rate)
        if(hy_auto_phase == 0) then
           hy_cpu_ref = hy_cpu_ref + hy_dt
        else if(hy_auto_phase == 1) then
           hy_gpu_ref = hy_gpu_ref + hy_dt
        endif
     endif
#endif

     ! Reverse update boundaries
                               call timer('hydro - rev ghostzones','start')
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
     end do
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
!    call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
                               call timer('poisson','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
     call diag_check_eint('godunov+sync',ilevel)

     ! Restriction operator
                               call timer('hydro upload fine','start')
     call upload_fine(ilevel)

  endif
 

  !---------------------
  ! Do RT/Chemistry step
  !---------------------
#ifdef RT
  if(rt .and. rt_advect) then  
                               call timer('radiative transfer','start')
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
                               call timer('cooling','start')
     if(neq_chem.or.cooling.or.T2_star>0.0)then
        cool_t1=omp_get_wtime()
        call cooling_fine(ilevel)
        cool_fine_wall=cool_fine_wall+omp_get_wtime()-cool_t1
     endif
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
                               call timer('radiative transfer','start')
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(aexp)
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
                               call timer('cooling','start')
     if(cosmo)then
        cool_t1=omp_get_wtime()
        call update_coolrates_tables(dble(aexp))
        cool_table_wall=cool_table_wall+omp_get_wtime()-cool_t1
     endif
                               call timer('radiative transfer','start')
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
                               call timer('cooling','start')
  if(neq_chem.or.cooling.or.T2_star>0.0)then
     cool_t1=omp_get_wtime()
     call cooling_fine(ilevel)
     cool_fine_wall=cool_fine_wall+omp_get_wtime()-cool_t1
  endif
  call diag_check_eint('cooling',ilevel)
#endif
#ifdef SNRT
  cool_t1=omp_get_wtime()
  call snrt_ramses_advance_level(ilevel)
  snrt_advance_wall=snrt_advance_wall+omp_get_wtime()-cool_t1
#ifndef SNRT_LEDGER_ONLY
  cool_t1=omp_get_wtime()
  call snrt_ramses_diagnose_level(ilevel)
  snrt_diagnose_wall=snrt_diagnose_wall+omp_get_wtime()-cool_t1
#endif
#endif
#ifdef SNRT_LEDGER_DIAGNOSTIC
  call snrt_agn_ledger_diagnose(ilevel)
#endif
  ! SGS turbulence source terms (production, dissipation, PdV coupling)
  if(use_sgs)call sgs_fine(ilevel)

  !---------------
  ! Move particles
  !---------------
  if(pic)then
                               call timer('particles','start')
     call system_clock(pt_t1)
     call move_fine(ilevel) ! Only remaining particles
     call system_clock(pt_t2)
     pt_move = pt_move + dble(pt_t2-pt_t1)/dble(pt_rate)
  end if

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
                               call timer('feedback','start')
  if(hydro.and.star)call star_formation(ilevel)

  ! Compute Bondi-Hoyle accretion parameters
                               call timer('sinks','start')
  call system_clock(sk_t1)
  if(sink.and.bondi)call bondi_hoyle(ilevel)
  call system_clock(sk_t2)
  sk_bondi_hoyle = sk_bondi_hoyle + dble(sk_t2-sk_t1)/dble(pt_rate)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
                               call timer('hydro - ghostzones','start')
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
                               call timer('hydro - diffusion','start')
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  ! Keep the parent FDM ghosts current for the next coarse-step's
  ! new-grid-only prolongation.  This must happen before flag_fine: the
  ! refinement map is state carried across the step boundary, so no unrelated
  ! communication belongs between flag_fine and the following refine_fine.
  if(use_fdm .and. ilevel<nlevelmax)then
     if(timer_report_interval>0) call timer('fdm-psi-ghost','start')
     call make_virtual_fine_dp(psi_re(1),ilevel)
     call make_virtual_fine_dp(psi_im(1),ilevel)
     if(timer_report_interval>0) call timer('flag','start')
  end if
  if(.not.static) call flag_fine(ilevel,icount)

  ! Accumulate exclusive rank x level work only on sparse sample steps.
  ! Count actual leaf cells, not all children of every AMR grid.
  if(lb_timing_sample)then
     t_lb_level_end=omp_get_wtime()
     call cpu_time(t_lb_cpu_level_end)
     level_time_loc(ilevel)=level_time_loc(ilevel)+ &
          max(0d0,t_lb_cpu_level_end-t_lb_cpu_level_start- &
          t_lb_cpu_child_time)
     nleaf_lb=0_8
     igrid_lb=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        do ind_lb=1,twotondim
           if(cpu_map(ICELL_OF(igrid_lb,ind_lb))==myid.and. &
                son(ICELL_OF(igrid_lb,ind_lb))==0)nleaf_lb=nleaf_lb+1_8
        end do
        igrid_lb=next(igrid_lb)
     end do
     level_ncells_loc(ilevel)=level_ncells_loc(ilevel)+nleaf_lb
     if(ilevel==levelmin)call update_work_timing_ema( &
          max(0d0,t_lb_level_end-t_lb_level_start))
  end if

  !----------------------------
  ! Merge finer level particles
  !----------------------------
                               call timer('particles','start')
  call system_clock(pt_t1)
  if(pic)call merge_tree_fine(ilevel)
  call system_clock(pt_t2)
  pt_merge = pt_merge + dble(pt_t2-pt_t1)/dble(pt_rate)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
                               call timer('aton','start')
     call rad_step(dtnew(ilevel))
  endif
#endif

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

#ifdef HYDRO_CUDA
  ! --- GPU auto-tuning: end-of-coarse-step decision for hydro & poisson ---
  if(ilevel==levelmin) then
     ! Hydro auto-tuning decision
     if(hy_auto_init) then
        if(hy_auto_phase == 0) then
           hy_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] Hydro CPU: ', hy_cpu_ref, ' s'
        else if(hy_auto_phase == 1) then
           hy_use_gpu = (hy_gpu_ref < hy_cpu_ref)
           hy_auto_phase = 2
           gpu_hydro = hy_use_gpu
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] Hydro GPU: ', hy_gpu_ref, &
                   ' s  (CPU was ', hy_cpu_ref, ' s)'
              if(hy_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Hydro decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Hydro decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: monitor accumulated time per coarse step
           if(hy_use_gpu) then
              hy_gpu_ref = 0d0  ! will re-accumulate next step
           else
              hy_cpu_ref = 0d0
           endif
        endif
     endif
     ! Poisson MG auto-tuning decision (gpu_fft excluded)
     if(mg_auto_init) then
        if(mg_auto_phase == 0) then
           mg_auto_phase = 1
           if(myid==1) write(*,'(A,F8.3,A)') &
                ' [GPU auto-tune] Poisson MG CPU: ', mg_cpu_ref, ' s'
        else if(mg_auto_phase == 1) then
           mg_use_gpu = (mg_gpu_ref < mg_cpu_ref)
           mg_auto_phase = 2
           if(.not. mg_use_gpu) then
              gpu_poisson = .false.
           endif
           if(myid==1) then
              write(*,'(A,F8.3,A,F8.3,A)') &
                   ' [GPU auto-tune] Poisson MG GPU: ', mg_gpu_ref, &
                   ' s  (CPU was ', mg_cpu_ref, ' s)'
              if(mg_use_gpu) then
                 write(*,'(A)') ' [GPU auto-tune] Poisson MG decision: GPU (faster)'
              else
                 write(*,'(A)') ' [GPU auto-tune] Poisson MG decision: CPU (faster)'
              endif
           endif
        else
           ! Phase 2+: monitor
           if(mg_use_gpu) then
              mg_gpu_ref = 0d0
           else
              mg_cpu_ref = 0d0
           endif
        endif
     endif
  endif
#endif

  ! Print particle & sink sub-timers every coarse step, then reset
  if(ilevel==levelmin .and. myid==1) then
     write(*,'(A)') ' === Direct cooling timings ==='
     write(*,'(A,F10.3,A)') '   cooling_fine : ', cool_fine_wall, ' s'
     write(*,'(A,F10.3,A)') '   cooling_table: ', cool_table_wall, ' s'
     write(*,'(A)') ' === Direct SNRT timings ==='
     write(*,'(A,F10.3,A)') '   snrt_advance : ', snrt_advance_wall, ' s'
     write(*,'(A,F10.3,A)') '   snrt_diagnose: ', snrt_diagnose_wall, ' s'
     write(*,'(A)') ' === Particle sub-timers ==='
     write(*,'(A,F10.3,A)') '   make_tree  : ', pt_mktree, ' s'
     write(*,'(A,F10.3,A)') '   kill+virt  : ', pt_killtree, ' s'
     write(*,'(A,F10.3,A)') '   synchro    : ', pt_synchro, ' s'
     write(*,'(A,F10.3,A)') '   move       : ', pt_move, ' s'
     write(*,'(A,F10.3,A)') '   merge      : ', pt_merge, ' s'
     write(*,'(A,F10.3,A)') '   TOTAL      : ', &
          pt_mktree+pt_killtree+pt_synchro+pt_move+pt_merge, ' s'
     write(*,'(A)') ' === Sink sub-timers ==='
     write(*,'(A,F10.3,A)') '   AGN_feedback : ', sk_agn_fb, ' s'
     write(*,'(A,F10.3,A)') '   create_sink  : ', sk_create_sink, ' s'
     write(*,'(A,F10.3,A)') '   grow_bondi   : ', sk_grow, ' s'
     write(*,'(A,F10.3,A)') '   bondi_hoyle  : ', sk_bondi_hoyle, ' s'
     write(*,'(A,F10.3,A)') '   TOTAL        : ', &
          sk_agn_fb+sk_create_sink+sk_grow+sk_bondi_hoyle, ' s'
     ! Reset for next coarse step
     pt_mktree=0; pt_killtree=0; pt_synchro=0; pt_move=0; pt_merge=0
     sk_agn_fb=0; sk_create_sink=0; sk_grow=0; sk_bondi_hoyle=0
     cool_fine_wall=0d0; cool_table_wall=0d0
     snrt_advance_wall=0d0; snrt_diagnose_wall=0d0
   end if

  ! Periodically emit and reset the existing mutually-exclusive phase timers.
  ! nstep_coarse is incremented by adaptive_loop after this routine returns,
  ! hence +1 identifies the coarse step that has just completed.  Starting
  ! timer-report first closes the currently active phase so the table does not
  ! lose the tail of this step.  The second start charges report overhead to
  ! the following interval instead of hiding it in a physics phase.
  if(ilevel==levelmin .and. timer_report_interval>0)then
     if(mod(nstep_coarse+1,timer_report_interval)==0)then
        call timer('timer-report','start')
        if(myid==1) write(*,'(A,I0,A,I0)') &
             ' PERF_TIMER_REPORT step=',nstep_coarse+1, &
             ' interval=',timer_report_interval
        call finalize_timer
        call fdm_reset_timer_running('timer-report')
     end if
  end if

999 format(' Entering amr_step',i1,' for level',i2, '  for a levelmin ',i3)

end subroutine amr_step

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

#ifdef RT
subroutine rt_step(ilevel)
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, t, dtnew, myid
  use rt_parameters, only: rt_isDiffuseUVsrc
  use rt_cooling_module, only: update_UVrates
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer, intent(in) :: ilevel

!--------------------------------------------------------------------------
!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries. 
!--------------------------------------------------------------------------

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=t ; t=t-t_left
  
  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2.0**(ilevel-levelmin))
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     ! If (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel    
     if (i_substep > 1) call rt_set_unew(ilevel)

     if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
     
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)
  
  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')
  
end subroutine rt_step
#endif
!###########################################################
!###########################################################
subroutine update_work_timing_ema(local_wall_time)
  ! Build the level cost model from process CPU time rather than elapsed wall
  ! time.  Process CPU time excludes ranks sleeping in MPI and counts work by
  ! all OpenMP threads, so a lightly loaded rank waiting in a collective is not
  ! misclassified as an intrinsically slow rank.
  use amr_commons
  use omp_lib, only: omp_get_max_threads
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),intent(in)::local_wall_time
  integer::ilevel,info,iref
  integer::nthreads
  real(dp)::alpha,ref_cpc,scale_min,scale_max,target_scale
  real(dp)::local_step_cpu,global_cpu_sum,global_cpu_max,global_cpu_min
  real(dp)::global_wall_sum,global_wall_max,cpu_mean,compute_imbalance
  real(dp)::omp_active,nrank,denom,cpc_average,cpc_slope,cpc_sample
  real(dp),dimension(1:MAXLEVEL,1:4)::stats_loc,stats_glob
  integer(kind=8),dimension(1:3)::fdm_work_loc,fdm_work_min,fdm_work_max

  alpha=max(0d0,min(1d0,lb_timing_ema_alpha))
  nrank=dble(max(1,ncpu))
  stats_loc=0d0
  stats_loc(:,1)=dble(level_ncells_loc)
  stats_loc(:,2)=level_time_loc
  stats_loc(:,3)=stats_loc(:,1)*stats_loc(:,1)
  stats_loc(:,4)=stats_loc(:,1)*stats_loc(:,2)
  local_step_cpu=sum(level_time_loc(levelmin:nlevelmax))
  fdm_work_loc(1)=sum(lb_cn_matvec_loc(levelmin:nlevelmax))
  fdm_work_loc(2)=sum(lb_cn_iter_loc(levelmin:nlevelmax))
  fdm_work_loc(3)=sum(lb_hjm_subcycle_loc(levelmin:nlevelmax))
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(stats_loc,stats_glob,4*MAXLEVEL, &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(local_step_cpu,global_cpu_sum,1, &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(local_step_cpu,global_cpu_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(local_step_cpu,global_cpu_min,1, &
       MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(local_wall_time,global_wall_sum,1, &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(local_wall_time,global_wall_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(fdm_work_loc,fdm_work_min,3,MPI_INTEGER8, &
       MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(fdm_work_loc,fdm_work_max,3,MPI_INTEGER8, &
       MPI_MAX,MPI_COMM_WORLD,info)
#else
  stats_glob=stats_loc
  global_cpu_sum=local_step_cpu
  global_cpu_max=local_step_cpu
  global_cpu_min=local_step_cpu
  global_wall_sum=local_wall_time
  global_wall_max=local_wall_time
  fdm_work_min=fdm_work_loc
  fdm_work_max=fdm_work_loc
#endif

  ! Fit T_cpu(level)=intercept+slope*N_leaf across ranks.  The intercept
  ! removes rank-independent setup/collective CPU overhead.  If the leaf-count
  ! variance is too small for a stable slope, use the global CPU/leaf average.
  do ilevel=levelmin,nlevelmax
     if(stats_glob(ilevel,1)>0d0.and.stats_glob(ilevel,2)>0d0)then
        cpc_average=stats_glob(ilevel,2)/stats_glob(ilevel,1)
        cpc_sample=cpc_average
        denom=nrank*stats_glob(ilevel,3)-stats_glob(ilevel,1)**2
        if(denom>1d-12*max(1d0,nrank*stats_glob(ilevel,3)))then
           cpc_slope=(nrank*stats_glob(ilevel,4)- &
                stats_glob(ilevel,1)*stats_glob(ilevel,2))/denom
           if(cpc_slope>0d0)cpc_sample=max(0.1d0*cpc_average, &
                min(10d0*cpc_average,cpc_slope))
        endif
        if(level_cell_time_ema(ilevel)<=0d0)then
           level_cell_time_ema(ilevel)=cpc_sample
        else
           level_cell_time_ema(ilevel)=(1d0-alpha)* &
                level_cell_time_ema(ilevel)+alpha*cpc_sample
        end if
     end if
  end do
  ! The old rank multiplier was derived from wall-time/leaf and could amplify
  ! collective wait imbalance.  Keep it neutral for all decomposition backends.
  level_rank_scale_ema=1d0

  iref=0
  do ilevel=levelmin,nlevelmax
     if(level_cell_time_ema(ilevel)>0d0)then
        iref=ilevel
        exit
     end if
  end do
  if(iref>0)then
     ref_cpc=level_cell_time_ema(iref)
     do ilevel=levelmin,nlevelmax
        if(level_cell_time_ema(ilevel)>0d0)then
           target_scale=max(0.25d0,min(8d0, &
                level_cell_time_ema(ilevel)/ref_cpc))
           level_mesh_scale_ema(ilevel)=max(0.25d0,min(8d0, &
                1d0+time_balance_alpha*(target_scale-1d0)))
        endif
     end do
  end if

  cpu_mean=global_cpu_sum/nrank
  compute_imbalance=0d0
  if(cpu_mean>0d0)compute_imbalance=max(0d0,global_cpu_max/cpu_mean-1d0)
  if(time_balance_alpha>0d0)then
     if(.not.lb_imbalance_ema_valid)then
        lb_imbalance_ema=compute_imbalance
        lb_imbalance_ema_valid=.true.
     else
        lb_imbalance_ema=(1d0-alpha)*lb_imbalance_ema+ &
             alpha*compute_imbalance
     end if
  endif

  ! Remap economics needs elapsed time, whereas the cost model needs CPU work.
  if(global_wall_max>0d0)then
     if(lb_step_time_ema<=0d0)then
        lb_step_time_ema=global_wall_max
     else
        lb_step_time_ema=(1d0-alpha)*lb_step_time_ema+ &
             alpha*global_wall_max
     end if
  end if

  if(myid==1.and.iref>0)then
     nthreads=max(1,omp_get_max_threads())
     omp_active=0d0
     if(global_wall_sum>0d0)omp_active=max(0d0,min(1d0, &
          global_cpu_sum/(dble(nthreads)*global_wall_sum)))
     scale_min=minval(level_mesh_scale_ema(levelmin:nlevelmax))
     scale_max=maxval(level_mesh_scale_ema(levelmin:nlevelmax))
     write(*,'(A,I8,A,F8.2,A,F6.2,A,F6.2,A,2F7.3)') &
          ' LB process timing step=',nstep_coarse,' wall-max=', &
          global_wall_max,' s, compute-balance=', &
          100d0/(1d0+compute_imbalance),'%, OMP-active=', &
          100d0*omp_active,'%, level-scale=',scale_min,scale_max
     write(*,'(A,3F10.3,A)') ' LB rank process-CPU min/mean/max=', &
          global_cpu_min,cpu_mean,global_cpu_max,' s'
     if(maxval(fdm_work_max)>0_8)write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
          ' LB FDM work ranges: CN-matvec=',fdm_work_min(1),'..', &
          fdm_work_max(1),', CN-iter=',fdm_work_min(2),'..', &
          fdm_work_max(2),', HJM-subcycle=',fdm_work_min(3),'..',fdm_work_max(3)
  end if

  level_time_loc=0d0
  level_ncells_loc=0_8
  lb_cn_matvec_loc=0_8
  lb_cn_iter_loc=0_8
  lb_hjm_subcycle_loc=0_8
end subroutine update_work_timing_ema
!###########################################################
!###########################################################
subroutine check_load_imbalance(did_remap)
  ! Check weight inhomogeneity and trigger load_balance if
  ! max/avg - 1 > remap_thresh.
  ! Called every coarse step when nremap==0 and remap_thresh>0.
  !
  ! The trigger and both decomposition algorithms use domain_leaf_cost, so
  ! remapping is judged with the same cost that the new domain cut balances.
  use amr_commons
  use pm_commons, only: count_particles_by_leaf
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical,intent(out)::did_remap

  real(dp)::my_cost,max_cost,sum_cost,imbalance,effective_imbalance
  real(dp)::predicted_saving,required_saving
  real(dp)::alpha
  integer::ilevel,info,jgrid,igrid,ind,npair_cell,steps_since
  integer,dimension(1:twotondim)::npart_leaf,ndm_leaf
  integer(kind=8)::my_cost_i8
  integer(kind=8)::cell_cost_i8
  integer(kind=8),dimension(1:MAXLEVEL)::niter_cost
  logical::worth_remap
#ifndef WITHOUTMPI
  real(dp)::buf(2),gbuf(2)
#endif

  did_remap=.false.

  niter_cost=1_8
  if((.not.memory_balance).and.cost_weighting)then
     niter_cost(levelmin)=1_8
     do ilevel=levelmin+1,nlevelmax
        if(niter_cost(ilevel-1)>huge(niter_cost(ilevel))/ &
             int(max(1,nsubcycle(ilevel-1)),kind=8))then
           if(myid==1)write(*,*)'check_load_imbalance: subcycle cost overflow'
           stop
        end if
        niter_cost(ilevel)=int(max(1,nsubcycle(ilevel-1)),kind=8)* &
             niter_cost(ilevel-1)
     end do
  end if

  my_cost_i8=0_8

  ! Coarse leaf cells are part of the domain ordering too.
  do ind=1,ncoarse
     if(cpu_map(ind)==myid.and.son(ind)==0)then
        cell_cost_i8=domain_leaf_cost(0,0,1_8, &
             level_mesh_scale_ema(levelmin))
        my_cost_i8=my_cost_i8+cell_cost_i8
     end if
  end do

  ! Sum the same actual-leaf particle and SIDM-pair proxy used by the cut.
  do ilevel=1,nlevelmax
     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        npart_leaf=0
        ndm_leaf=0
        if(pic)call count_particles_by_leaf(igrid,npart_leaf,ndm_leaf)
        do ind=1,twotondim
           if(cpu_map(ICELL_OF(igrid,ind))/=myid.or. &
                son(ICELL_OF(igrid,ind))/=0)cycle
           npair_cell=domain_sidm_pair_count(ndm_leaf(ind))
           cell_cost_i8=domain_leaf_cost(npart_leaf(ind),npair_cell, &
                niter_cost(max(levelmin,ilevel)), &
                level_mesh_scale_ema(max(levelmin,ilevel)))
           my_cost_i8=my_cost_i8+cell_cost_i8
        end do
        igrid=next(igrid)
     end do
  end do
  my_cost=dble(my_cost_i8)

#ifndef WITHOUTMPI
  buf(1)=my_cost
  buf(2)=my_cost
  call MPI_ALLREDUCE(buf(1),max_cost,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(buf(2),sum_cost,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  max_cost=my_cost
  sum_cost=my_cost
#endif

  if(sum_cost>0d0)then
     imbalance=(max_cost-sum_cost/dble(ncpu))/(sum_cost/dble(ncpu))
  else
     imbalance=0d0
  end if

  ! Timed mode uses the process-CPU imbalance measured by
  ! update_work_timing_ema.  Replacing it here with the model imbalance would
  ! make the adaptive trigger circular: the same proxy would choose the cut
  ! and then declare that cut balanced.  Analytic work and memory modes retain
  ! the inexpensive model-based trigger.
  if((.not.memory_balance).and.time_balance_alpha>0d0)then
     if(lb_imbalance_ema_valid)then
        effective_imbalance=lb_imbalance_ema
     else
        ! No process-CPU sample exists yet (notably just after a remap).
        ! Use the current proxy once, but do not seed the measured EMA with it.
        effective_imbalance=imbalance
     endif
  else
     alpha=max(0d0,min(1d0,lb_timing_ema_alpha))
     if(.not.lb_imbalance_ema_valid)then
        lb_imbalance_ema=imbalance
        lb_imbalance_ema_valid=.true.
     else
        lb_imbalance_ema=(1d0-alpha)*lb_imbalance_ema+alpha*imbalance
     end if
     effective_imbalance=lb_imbalance_ema
  end if

  steps_since=nstep_coarse-lb_last_remap_step
  worth_remap=work_remap_is_economic(effective_imbalance, &
       lb_step_time_ema,lb_remap_time_ema,steps_since, &
       predicted_saving,required_saving)

  if(worth_remap)then
     if(myid==1) write(*,'(A,F6.2,A,F6.2,A,F8.1,A,F8.1,A)') &
          ' Load imbalance EMA ',effective_imbalance*100d0, &
          '% > threshold ',remap_thresh*100d0, &
          '%, predicted/required saving=',predicted_saving,'/', &
          required_saving,' s -> rebalancing'
     call load_balance
     call defrag
     did_remap=.true.
  else if(verbose.and.myid==1.and.effective_imbalance>remap_thresh)then
     write(*,'(A,F6.2,A,I0,A,F8.1,A,F8.1,A)') &
          ' Load imbalance EMA ',effective_imbalance*100d0, &
          '%; steps since remap=',steps_since, &
          ', predicted/required=',predicted_saving,'/',required_saving, &
          ' s -> keep current domains'
  end if

end subroutine check_load_imbalance

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine load_balance
  use amr_commons
  use pm_commons
  use hydro_commons, ONLY: nvar, uold
#ifdef RT
  use rt_hydro_commons, ONLY: nrtvar, rtuold
#endif
  use poisson_commons, ONLY: phi, f, scalar_gr, scalar_gr_old, psi_re, psi_im
  use bisection
  use ksection
  use iso_c_binding, only: c_int, c_size_t
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::igrid,ncache,ilevel,i,ind,jlevel,info
  integer::idim,ivar,icpu,jcpu,kcpu
  integer::nxny,ix,iy,iz
  integer::required_grid_capacity
  integer(i8b),dimension(nlevelmax,3)::comm_buffin,comm_buffout
  integer,allocatable::numbp_save(:)
  integer::nsave,isave
  integer,save::lb_chain_depth=0
  logical::continue_bounded_remap
#ifndef WITHOUTMPI
  integer::countsend,countrecv
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
  integer,dimension(ncpu)::reqsend,reqrecv
  real(dp)::t_lb_start,t_lb_end
  real(dp)::remap_elapsed,remap_elapsed_max
  real(dp)::t0,t1,t2,t3,t4,t5,t6
  real(dp)::t_tree_collapse,t_numbp_done,t_stats_done,t_map_done,t_particle_done
  real(dp)::t_tree_stage0
  real(dp),dimension(1:nlevelmax)::t_tree_kill_loc,t_tree_virtual_loc,t_tree_merge_loc
  real(dp),dimension(1:nlevelmax)::t_tree_kill_max,t_tree_virtual_max,t_tree_merge_max
  real(dp)::te_flag,te_refine,te_bcomm,te_virt,te_phys
  real(dp)::ts_flag,ts_refine,ts_bcomm
  real(dp)::tt0,tt1
#endif

  if(ncpu==1)return

  ! [RESIZABLE] Production load balancing reserves grid headroom using the
  ! same capacity on every rank.  Do this before particle-tree changes or
  ! nonblocking communication so the collective has a rank-symmetric safe
  ! point.  The extra grid preserves the existing strict free-slot guard at
  ! an exact headroom boundary.
  if(ngridmax_auto .and. lb_grid_headroom>0d0 .and. &
       lb_grid_headroom<1d0)then
     required_grid_capacity=ceiling((dble(used_mem)+1d0)/ &
          lb_grid_headroom)
     call ensure_grid_capacity_collective(required_grid_capacity, &
          'load_balance_headroom')
  endif

  lb_chain_depth=lb_chain_depth+1

#ifndef WITHOUTMPI
  t_lb_start = MPI_WTIME()
  t0 = t_lb_start
  t_tree_collapse=t0
  t_numbp_done=t0
  if(myid==1)write(*,*)'Load balancing AMR grid...'

  ! Put all particle in main tree trunk
  if(pic.and.(.not.init))then
     call make_tree_fine(levelmin)
     do ilevel=levelmin-1,1,-1
        call merge_tree_fine(ilevel)
     end do
  endif
  t_tree_collapse=MPI_WTIME()

  ! The fast memory-balance path starts with the exact particle totals stored
  ! on levelmin grids and propagates them deterministically down the AMR tree.
  ! A levelmin grid can be a reception grid on the rank owning one of its
  ! refined descendants, so only that level needs numbp synchronized.  Exact
  ! linked-list counting visits active grids only and skips this communication.
  !
  ! After make_tree_fine/merge_tree_fine, numbp is correct for active grids but
  ! stale for reception grids. We temporarily overwrite reception grids'
  ! numbp with the remote active grid's value for the fast cost function.
  ! We save the original numbp before overwriting, and restore afterwards.
  ! Setting numbp=0 would break the particle tree because merge_tree_fine
  ! can attach particles to reception grids at coarser levels.
  if(memory_balance .and. memory_balance_fast_particles .and. &
       (.not.use_cpubox_decomp) .and. pic .and. &
       (.not.init))then
     ! Count total reception grids for save buffer
     nsave=0
     do ilevel=levelmin,levelmin
        do icpu=1,ncpu
           nsave=nsave+reception(icpu,ilevel)%ngrid
        end do
     end do
     if(nsave>0) allocate(numbp_save(1:nsave))
     isave=0
     do ilevel=levelmin,levelmin
        ! Post receives into reception%f(:,1) — contiguous buffer
        countrecv=0
        do icpu=1,ncpu
           ncache=reception(icpu,ilevel)%ngrid
           if(ncache>0) then
              ! Save original numbp before overwriting
              do i=1,ncache
                 isave=isave+1
                 numbp_save(isave)=numbp(reception(icpu,ilevel)%igrid(i))
              end do
              countrecv=countrecv+1
              call MPI_IRECV(reception(icpu,ilevel)%f(1,1),ncache, &
                   & MPI_INTEGER,icpu-1,199,MPI_COMM_WORLD,reqrecv(countrecv),info)
           end if
        end do
        ! Pack numbp for emission grids and send
        countsend=0
        do icpu=1,ncpu
           ncache=emission(icpu,ilevel)%ngrid
           if(ncache>0) then
              do i=1,ncache
                 emission(icpu,ilevel)%f(i,1)=numbp(emission(icpu,ilevel)%igrid(i))
              end do
              countsend=countsend+1
              call MPI_ISEND(emission(icpu,ilevel)%f(1,1),ncache, &
                   & MPI_INTEGER,icpu-1,199,MPI_COMM_WORLD,reqsend(countsend),info)
           end if
        end do
        ! Wait for all communication
        call MPI_WAITALL(countrecv,reqrecv,statuses,info)
        call MPI_WAITALL(countsend,reqsend,statuses,info)
        ! Scatter received remote numbp values into numbp array
        do icpu=1,ncpu
           ncache=reception(icpu,ilevel)%ngrid
           do i=1,ncache
              numbp(reception(icpu,ilevel)%igrid(i))=reception(icpu,ilevel)%f(i,1)
           end do
        end do
     end do
  end if

  t_numbp_done=MPI_WTIME()
  t1 = t_numbp_done

  balance=.true.

  if(verbose)then
     write(*,*)'Input mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  !-------------------------------------------
  ! Compute new cpu map using chosen ordering
  !-------------------------------------------
  call cmp_new_cpu_map

  ! Restore original numbp for reception grids (undo the fast-path sync).
  ! merge_tree_fine can attach particles to reception grids, so we must
  ! restore the original values rather than setting to 0.
  if(memory_balance .and. memory_balance_fast_particles .and. &
       (.not.use_cpubox_decomp) .and. pic .and. &
       (.not.init))then
     isave=0
     do ilevel=levelmin,levelmin
        do icpu=1,ncpu
           ncache=reception(icpu,ilevel)%ngrid
           do i=1,ncache
              isave=isave+1
              numbp(reception(icpu,ilevel)%igrid(i))=numbp_save(isave)
           end do
        end do
     end do
     if(allocated(numbp_save)) deallocate(numbp_save)
  end if

  ! No ownership boundary moved when the preflight found too little working
  ! space.  Avoid an otherwise pointless expand/shrink cycle: even a no-op
  ! remap temporarily consumes ghost slots and is exactly what the guard is
  ! intended to defer.  cmp_new_cpu_map is called only after the particle
  ! tree has been collapsed to the level-1 trunk above.  A normal remap
  ! scatters the particles again before returning (see the matching block
  ! below); the no-op path must do the same with the still-current cpu map.
  ! Otherwise real and virtual particles remain mixed in the trunk and a
  ! later sink-cloud rebuild can spend an unbounded time walking inconsistent
  ! particle lists.
  if(lb_remap_fraction<=0d0)then
     if(pic.and.(.not.init))then
        if(myid==1) then
           write(*,*) 'Bounded remap: restoring particle tree before no-op return'
           call flush(6)
        end if
        do ilevel=1,nlevelmax-1
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        end do
        call virtual_tree_fine(nlevelmax)
        do ilevel=nlevelmax-1,levelmin,-1
           call merge_tree_fine(ilevel)
        end do
        if(myid==1) then
           write(*,*) 'Bounded remap: particle tree restore complete'
           call flush(6)
        end if
     end if
     balance=.false.
     if(myid==1) write(*,*) 'Bounded remap: no safe progress; keeping current map'
     lb_chain_depth=lb_chain_depth-1
     return
  endif

  t2 = MPI_WTIME()

  !------------------------------------------------------
  ! Expand boundaries to account for new mesh partition
  !------------------------------------------------------
  te_flag=0d0; te_refine=0d0; te_bcomm=0d0; te_virt=0d0; te_phys=0d0

  tt0 = MPI_WTIME()
  call flag_coarse
  call refine_coarse
  tt1 = MPI_WTIME(); te_refine = te_refine + (tt1-tt0)

  tt0 = tt1
  call build_comm(1)
  tt1 = MPI_WTIME(); te_bcomm = te_bcomm + (tt1-tt0)

  tt0 = tt1
  call make_virtual_fine_int_pair(cpu_map(1),cpu_map2(1),1)
  tt1 = MPI_WTIME(); te_virt = te_virt + (tt1-tt0)

  do i=1,nlevelmax-1
     tt0 = MPI_WTIME()
     call flag_fine(i,2)
     tt1 = MPI_WTIME(); te_flag = te_flag + (tt1-tt0)

     tt0 = tt1
     call refine_fine(i)
     tt1 = MPI_WTIME(); te_refine = te_refine + (tt1-tt0)

     tt0 = tt1
     call build_comm(i+1)
     tt1 = MPI_WTIME(); te_bcomm = te_bcomm + (tt1-tt0)

     tt0 = tt1
     call make_virtual_fine_int_pair(cpu_map(1),cpu_map2(1),i+1)
     tt1 = MPI_WTIME(); te_virt = te_virt + (tt1-tt0)
  end do

  !--------------------------------------
  ! Update physical boundary conditions
  !--------------------------------------
  tt0 = MPI_WTIME()
  do ilevel=nlevelmax,1,-1
     if(hydro)then
#ifdef SOLVERmhd
        call make_virtual_fine_dp_bulk(uold,nvar+3,ilevel)
#else
        call make_virtual_fine_dp_bulk(uold,nvar,ilevel)
#endif
        if(simple_boundary)then
           call make_boundary_hydro(ilevel)
        end if
     end if
#ifdef RT
     if(rt)then
        call make_virtual_fine_dp_bulk(rtuold,nrtvar,ilevel)
        if(simple_boundary)then
           call rt_make_boundary_hydro(ilevel)
        end if
     endif
#endif
     if(poisson)then
        call make_virtual_fine_dp(phi(1),ilevel)
        call make_virtual_fine_dp_bulk(f,ndim,ilevel)
        if(allocated(scalar_gr))then
           call make_virtual_fine_dp(scalar_gr(1),ilevel)
           call make_virtual_fine_dp(scalar_gr_old(1),ilevel)
        end if
     end if
     if(use_fdm)then
        call make_virtual_fine_dp(psi_re(1),ilevel)
        call make_virtual_fine_dp(psi_im(1),ilevel)
     end if
  end do
  tt1 = MPI_WTIME(); te_phys = tt1 - tt0

  t3 = MPI_WTIME()

  !--------------------------------------
  ! Rearrange octs between cpus
  !--------------------------------------
  do ilevel=1,nlevelmax
     do icpu=1,ncpu
        if(icpu==myid)then
           ncache=active(ilevel)%ngrid
        else
           ncache=reception(icpu,ilevel)%ngrid
        end if
        ! Disconnect from old linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
              igrid=reception(icpu,ilevel)%igrid(i)
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu)then
              if(prev(igrid).ne.0) then
                 if(next(igrid).ne.0)then
                    next(prev(igrid))=next(igrid)
                    prev(next(igrid))=prev(igrid)
                 else
                    next(prev(igrid))=0
                    taill(kcpu,ilevel)=prev(igrid)
                 end if
              else
                 if(next(igrid).ne.0)then
                    prev(next(igrid))=0
                    headl(kcpu,ilevel)=next(igrid)
                 else
                    headl(kcpu,ilevel)=0
                    taill(kcpu,ilevel)=0
                 end if
              end if
              numbl(kcpu,ilevel)=numbl(kcpu,ilevel)-1 
           end if
        end do        
        ! Connect to new linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
              igrid=reception(icpu,ilevel)%igrid(i)
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu)then
              if(numbl(jcpu,ilevel)>0)then
                 next(igrid)=0
                 prev(igrid)=taill(jcpu,ilevel)
                 next(taill(jcpu,ilevel))=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=numbl(jcpu,ilevel)+1
              else
                 next(igrid)=0
                 prev(igrid)=0
                 headl(jcpu,ilevel)=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=1
              end if
           end if
        end do
     end do
  end do
  t4 = MPI_WTIME()

  !--------------------------------------
  ! Compute new grid number statistics
  !--------------------------------------
  do ilevel=1,nlevelmax
     comm_buffin(ilevel,1)=numbl(myid,ilevel)
     comm_buffin(ilevel,2)=numbl(myid,ilevel)
     comm_buffin(ilevel,3)=numbl(myid,ilevel)
  end do
#ifndef LONGINT
  call MPI_ALLREDUCE(comm_buffin(1,1),comm_buffout(1,1),nlevelmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,2),comm_buffout(1,2),nlevelmax,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,3),comm_buffout(1,3),nlevelmax,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(comm_buffin(1,1),comm_buffout(1,1),nlevelmax,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,2),comm_buffout(1,2),nlevelmax,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,3),comm_buffout(1,3),nlevelmax,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
#endif
  call MPI_ALLREDUCE(used_mem        ,used_mem_tot     ,1        ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  do ilevel=1,nlevelmax
     numbtot(1,ilevel)=comm_buffout(ilevel,1)
     numbtot(2,ilevel)=comm_buffout(ilevel,2)
     numbtot(3,ilevel)=comm_buffout(ilevel,3)
     numbtot(4,ilevel)=numbtot(1,ilevel)/ncpu
  end do
  t_stats_done=MPI_WTIME()

  !--------------------------------------
  ! Set old cpu map to new cpu map
  !--------------------------------------
  if(.not.use_cpubox_decomp) then
     bound_key=bound_key2
  else
     bisec_cpubox_min=bisec_cpubox_min2
     bisec_cpubox_max=bisec_cpubox_max2
  end if

  nxny=nx*ny
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(SHARED) PRIVATE(ix,iy,iz,ind) &
  !$OMP SCHEDULE(STATIC)
  do iz=kcoarse_min,kcoarse_max
  do iy=jcoarse_min,jcoarse_max
  do ix=icoarse_min,icoarse_max
     ind=1+ix+iy*nx+iz*nxny
     cpu_map(ind)=cpu_map2(ind)
  end do
  end do
  end do
  !$OMP END PARALLEL DO
  do ilevel=1,nlevelmax
     ! Build new communicators
     call build_comm(ilevel)
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ind) SCHEDULE(STATIC)
     do i=1,active(ilevel)%ngrid
        do ind=1,twotondim
           cpu_map(ICELL_OF(active(ilevel)%igrid(i),ind))= &
                cpu_map2(ICELL_OF(active(ilevel)%igrid(i),ind))
        end do
     end do
     !$OMP END PARALLEL DO
     call make_virtual_fine_int(cpu_map(1),ilevel)
  end do
  t_map_done=MPI_WTIME()

  t_tree_kill_loc=0d0
  t_tree_virtual_loc=0d0
  t_tree_merge_loc=0d0
  t_tree_kill_max=0d0
  t_tree_virtual_max=0d0
  t_tree_merge_max=0d0
  if(pic.and.(.not.init))then
     ! Sort particles down to nlevelmax
     do ilevel=1,nlevelmax-1
        t_tree_stage0=MPI_WTIME()
        call kill_tree_fine(ilevel)
        t_tree_kill_loc(ilevel)=MPI_WTIME()-t_tree_stage0
        t_tree_stage0=MPI_WTIME()
        call virtual_tree_fine(ilevel)
        t_tree_virtual_loc(ilevel)=MPI_WTIME()-t_tree_stage0
     end do
     t_tree_stage0=MPI_WTIME()
     call virtual_tree_fine(nlevelmax)
     t_tree_virtual_loc(nlevelmax)=MPI_WTIME()-t_tree_stage0
     do ilevel=nlevelmax-1,levelmin,-1
        t_tree_stage0=MPI_WTIME()
        call merge_tree_fine(ilevel)
        t_tree_merge_loc(ilevel)=MPI_WTIME()-t_tree_stage0
     end do
     call MPI_ALLREDUCE(t_tree_kill_loc,t_tree_kill_max,nlevelmax, &
          MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(t_tree_virtual_loc,t_tree_virtual_max,nlevelmax, &
          MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(t_tree_merge_loc,t_tree_merge_max,nlevelmax, &
          MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  end if

  t_particle_done=MPI_WTIME()
  t5 = t_particle_done

  !--------------------------------------------
  ! Shrink boundaries around new mesh partition
  !--------------------------------------------
  ts_flag=0d0; ts_refine=0d0; ts_bcomm=0d0
  shrink=.true.
  do i=nlevelmax-1,1,-1
     tt0 = MPI_WTIME()
     call flag_fine(i,2)
     tt1 = MPI_WTIME(); ts_flag = ts_flag + (tt1-tt0)

     tt0 = tt1
     call refine_fine(i)
     tt1 = MPI_WTIME(); ts_refine = ts_refine + (tt1-tt0)

     tt0 = tt1
     call build_comm(i+1)
     tt1 = MPI_WTIME(); ts_bcomm = ts_bcomm + (tt1-tt0)
  end do
  tt0 = MPI_WTIME()
  call flag_coarse
  call refine_coarse
  tt1 = MPI_WTIME(); ts_refine = ts_refine + (tt1-tt0)
  tt0 = tt1
  call build_comm(1)
  tt1 = MPI_WTIME(); ts_bcomm = ts_bcomm + (tt1-tt0)
  shrink=.false.

  balance=.false.

  t6 = MPI_WTIME()
  t_lb_end = t6
  remap_elapsed=t_lb_end-t_lb_start
  call MPI_ALLREDUCE(remap_elapsed,remap_elapsed_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  if(lb_remap_time_ema<=0d0)then
     lb_remap_time_ema=remap_elapsed_max
  else
     lb_remap_time_ema=(1d0-lb_timing_ema_alpha)*lb_remap_time_ema+ &
          lb_timing_ema_alpha*remap_elapsed_max
  end if
  lb_last_remap_step=nstep_coarse
  lb_imbalance_ema_valid=.false.
  if(myid==1) then
     write(*,'(A,F8.3,A)') ' load_balance total:         ', t_lb_end - t_lb_start, ' s'
     write(*,'(A,F8.3,A)') '   particle_tree_collapse:   ', t_tree_collapse - t0, ' s'
     write(*,'(A,F8.3,A)') '   fast_child_numbp_sync:    ', t_numbp_done - t_tree_collapse, ' s'
     write(*,'(A,F8.3,A)') '   cmp_new_cpu_map:          ', t2 - t1, ' s'
     write(*,'(A,F8.3,A)') '   expand_pass:              ', t3 - t2, ' s'
     write(*,'(A,F8.3,A)') '     flag_fine:              ', te_flag, ' s'
     write(*,'(A,F8.3,A)') '     refine:                 ', te_refine, ' s'
     write(*,'(A,F8.3,A)') '     build_comm:             ', te_bcomm, ' s'
     write(*,'(A,F8.3,A)') '     virtual_int_pair:       ', te_virt, ' s'
     write(*,'(A,F8.3,A)') '     phys_boundary:          ', te_phys, ' s'
     write(*,'(A,F8.3,A)') '   grid_migration:           ', t4 - t3, ' s'
     write(*,'(A,F8.3,A)') '   grid_stats_allreduce:     ', t_stats_done - t4, ' s'
     write(*,'(A,F8.3,A)') '   cpumap_owner_update:      ', t_map_done - t_stats_done, ' s'
     write(*,'(A,F8.3,A)') '   particle_tree_rebuild:    ', t_particle_done - t_map_done, ' s'
     if(pic.and.(.not.init))then
        do ilevel=1,nlevelmax
           if(t_tree_kill_max(ilevel)+t_tree_virtual_max(ilevel)+ &
                t_tree_merge_max(ilevel)>0d0)then
              write(*,'(A,I2,A,3(F8.3,A))') '     level ',ilevel, &
                   ': kill=',t_tree_kill_max(ilevel),' s virtual=', &
                   t_tree_virtual_max(ilevel),' s merge=', &
                   t_tree_merge_max(ilevel),' s'
           end if
        end do
     end if
     write(*,'(A,F8.3,A)') '   shrink_pass:              ', t6 - t5, ' s'
     write(*,'(A,F8.3,A)') '     flag_fine:              ', ts_flag, ' s'
     write(*,'(A,F8.3,A)') '     refine:                 ', ts_refine, ' s'
     write(*,'(A,F8.3,A)') '     build_comm:             ', ts_bcomm, ' s'
  end if

  ! Release grow-only buffers after rebalancing (comm patterns changed)
  call ksection_trim_buffers()

  ! Particle migration updates the free list, but not every remap path updates
  ! the cached local count.  Output and restart buffers are sized from npart,
  ! so synchronize it before returning from the load-balance safe point.
  if(pic.and.particle_free_list_ready)then
     npart=npartmax-numbp_free
     call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN, &
          MPI_COMM_WORLD,info)
  endif

  ! Return freed heap pages to OS (reduces RSS after bulk dealloc/realloc)
  call fortran_malloc_trim()

  ! A bounded remap deliberately moves only part of a large ownership change.
  ! The standard shrink pass above has now returned the old virtual grids to
  ! the free list, so another ordinary remap can safely advance the boundary.
  continue_bounded_remap = lb_remap_fraction>0d0 .and. &
       lb_remap_fraction<0.999999d0 .and. lb_chain_depth<16
  if(continue_bounded_remap)then
     if(myid==1) write(*,'(A,I0,A)') ' Bounded remap: starting round ', &
          lb_chain_depth+1,' after releasing old grid slots'
     call load_balance
  else if(lb_remap_fraction>0d0 .and. lb_remap_fraction<0.999999d0)then
     if(myid==1) write(*,*) ' Bounded remap: round limit reached; remaining motion deferred'
  end if

  if(verbose)then
     write(*,*)'Output mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if
#endif

  lb_chain_depth=lb_chain_depth-1

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_new_cpu_map
  use amr_commons
  use pm_commons
  use bisection
  use ksection
  use omp_lib, only: omp_get_wtime
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !---------------------------------------------------
  ! This routine computes the new cpu map using 
  ! the choosen ordering to balance load across cpus.
  !---------------------------------------------------
  integer::igrid,ncell,ncell_loc,ncache,ngrid
  integer::ncode,bit_length,ilevel,i,ind,idim
  integer::nx_loc,ny_loc,nz_loc,nfar
  integer::info,icpu,jcpu,isub,idom,jdom
  integer::nxny,ix,iy,iz,ilo,ihi,imid
  integer::ind_long
  integer::isink,igrid_sink,ind_sink,icell_sink,isubcell_sink
  integer::npair_cell,jgrid,parent_cell,parent_grid,parent_ind
  integer,dimension(1:nvector)::ind_grid,ind_cell
  integer,dimension(1:nvector,1:twotondim)::npart_leaf,ndm_leaf
  integer,allocatable::sink_per_grid(:)
  integer,allocatable::grid_particle_budget(:)
  logical,dimension(1:twotondim)::all_children

  real(dp)::dx,scale,weight
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(kind=8)::incost_tot,local_cost,cell_cost
  real(kind=8),dimension(0:ndomain)::incost_new,incost_old
  integer(kind=8),dimension(1:overload)::npart_sub
  integer(kind=8)::wflag
  integer(kind=8)::nraised_loc,nraised,ntot_grids_loc,ntot_grids
  integer,dimension(1:overload)::ncell_sub
  real(kind=8),dimension(1:ndomain)::cost_loc,cost_old,cost_new
  real(qdp),dimension(0:ndomain)::bound_key_loc
  real(qdp),dimension(0:ndomain)::lb_bound_key_target
  real(kind=8),dimension(0:ndomain)::bigdbl,bigtmp
  integer,dimension(1:nvector)::dom
  real(qdp),dimension(1:nvector)::order_min,order_max
  integer(kind=8),dimension(1:MAXLEVEL)::niter_cost

  real(dp),dimension(1:1,1:ndim) :: xx_tmp
  integer,dimension(1:1) :: c_tmp
  ! Vector scratch for ksection cpu-map (Pass 2)
  integer,dimension(1:nvector) :: c_tmp_v

  ! OMP variables for parallelized remap
  integer::batch_size, my_base, my_idx
  integer,dimension(1:overload)::ncell_sub_t
  integer(kind=8),dimension(1:overload)::npart_sub_t

  real(kind=8) :: floor_w,min_weight_loc,min_weight_global
  real(kind=8) :: grid_cap,guard_denom,predicted_maxcount,cost_imbalance
  real(kind=8) :: grid_avail
  integer(kind=8)::ngrid_own_loc,ngrid_own_max,ngrid_ext_loc,ngrid_ext_max
  integer::ilev_g,icpu_g,ibnd_g
  ! Every grid but the root is the child of exactly one refined cell, so an AMR
  ! tree of G grids carries 8G-(G-1) = 7G+1 leaf cells.  flag1 holds one entry
  ! per leaf cell, so this is the constant that converts a grid budget into the
  ! cell budget the guard works in.  Using twotondim here overstates capacity by
  ! 8/7 and helped sink job 399652 on 2026-08-03.
  real(kind=8),parameter :: LB_LEAF_PER_GRID=7d0
  integer :: guard_iter,floor_flag
  integer,parameter :: LB_GRID_GUARD_MAXITER=3
  real(kind=8),parameter :: LB_REMAP_SLOT_SHARE=0.5d0
  logical :: guard_applied
  integer,dimension(1:ncpu)::lb_nsend,lb_nrecv
  integer::lb_incoming,lb_free_slots,lb_allowed,lb_target_cpu,lb_limit_iter
  real(kind=8)::lb_fraction_local,lb_fraction_global
  logical::lb_boundary_limited
  logical::use_fast_particle_balance,need_dm_count
  integer(kind=8)::fast_grid_loc,fast_fallback_loc,exact_particle_visits_loc
  integer(kind=8)::particle_assigned_loc,particle_expected_loc
  integer(kind=8)::particle_physical_loc,particle_physical_tot
  integer(kind=8)::particle_conservation_error,particle_conservation_tol
  integer(kind=8)::fast_grid_tot,fast_fallback_tot,exact_particle_visits_tot
  integer(kind=8)::particle_assigned_tot,particle_expected_tot
  real(dp)::tpart_work_loc,tkey_work_loc,tcost_work_loc,tt_work
  real(dp)::tpart_work_max,tkey_work_max,tcost_work_max
#ifndef WITHOUTMPI
  real(dp)::tcmp_start,tcmp_key,tcmp_sort,tcmp_bound,tcmp_map,tcmp_virtual
#endif

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  lb_remap_fraction=1d0
  lb_boundary_limited=.false.
  lb_limit_iter=0
  lb_bound_key_target=0.0_qdp
  use_fast_particle_balance=memory_balance.and. &
       memory_balance_fast_particles.and.(.not.use_cpubox_decomp)
  need_dm_count=(.not.memory_balance).and.sidm.and.work_weight_sidm_pair>0
  fast_grid_loc=0_8
  fast_fallback_loc=0_8
  exact_particle_visits_loc=0_8
  particle_assigned_loc=0_8
  particle_expected_loc=0_8
  particle_physical_loc=int(max(0,npart),kind=8)
  tpart_work_loc=0d0
  tkey_work_loc=0d0
  tcost_work_loc=0d0
  all_children=.true.
  
  ! Compute AMR subcycle work factors.  Memory balance deliberately ignores
  ! these factors because allocated memory does not grow with update count.
  niter_cost=1_8
  if((.not.memory_balance).and.cost_weighting)then
     niter_cost(levelmin)=1_8
     do ilevel=levelmin+1,nlevelmax
        if(niter_cost(ilevel-1)>huge(niter_cost(ilevel))/ &
             int(max(1,nsubcycle(ilevel-1)),kind=8))then
           if(myid==1)write(*,*)'load_balance: AMR subcycle cost overflow'
           stop
        end if
        niter_cost(ilevel)=int(max(1,nsubcycle(ilevel-1)),kind=8)* &
             niter_cost(ilevel-1)
     end do
  endif

  if(verbose) print *,"Entering cmp_new_cpu_map"

#ifndef WITHOUTMPI
  tcmp_start=MPI_WTIME()
  tcmp_key=tcmp_start
  tcmp_sort=tcmp_start
  tcmp_bound=tcmp_start
  tcmp_map=tcmp_start
  tcmp_virtual=tcmp_start
#endif

  if(.not.use_cpubox_decomp) then      ! begin if not bisection/ksection

  ! Build per-grid sink particle count for cost weighting
  if(memory_balance .and. sink .and. nsink > 0 .and. mem_weight_sink > 0) then
     allocate(sink_per_grid(1:ngridmax))
     sink_per_grid = 0
     do isink = 1, nsink
        ix = int(xsink(isink,1) / scale)
        iy = int(xsink(isink,2) / scale)
        iz = int(xsink(isink,3) / scale)
        ix = min(max(ix, 0), nx-1)
        iy = min(max(iy, 0), ny-1)
        iz = min(max(iz, 0), nz-1)
        icell_sink = 1 + ix + iy*nx + iz*nxny
        do while(son(icell_sink) > 0)
           igrid_sink = son(icell_sink)
           ind_sink = 1
           if(xsink(isink,1) >= xg(igrid_sink,1)) ind_sink = ind_sink + 1
           if(xsink(isink,2) >= xg(igrid_sink,2)) ind_sink = ind_sink + 2
           if(ndim > 2) then
              if(xsink(isink,3) >= xg(igrid_sink,3)) ind_sink = ind_sink + 4
           end if
           icell_sink = ICELL_OF(igrid_sink,ind_sink)
        end do
        if(icell_sink > ncoarse) then
           isubcell_sink = ICHILD_OF(icell_sink)
           igrid_sink = IGRID_OF(icell_sink)
           sink_per_grid(igrid_sink) = sink_per_grid(igrid_sink) + 1
        end if
     end do
  end if

  ! During load_balance all particles are in the levelmin trunk.  Propagate
  ! each exact levelmin total down its AMR branches by the same deterministic
  ! quotient/remainder split used for leaf costs.  Refined-cell shares are
  ! carried by the child grid; unrefined-cell shares become leaf costs.  This
  ! conserves the global particle total without touching a linked list and
  ! handles the initial zoom state where finer-grid numbp values are zero.
  if(use_fast_particle_balance.and.pic)then
     allocate(grid_particle_budget(1:ngridmax))
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
     do i=1,ngridmax
        grid_particle_budget(i)=0
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,igrid) SCHEDULE(STATIC) &
     !$OMP REDUCTION(+:particle_expected_loc)
     do i=1,active(levelmin)%ngrid
        igrid=active(levelmin)%igrid(i)
        grid_particle_budget(igrid)=max(0,numbp(igrid))
        particle_expected_loc=particle_expected_loc+ &
             int(grid_particle_budget(igrid),kind=8)
     end do
     !$OMP END PARALLEL DO
     do icpu=1,ncpu
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,igrid) SCHEDULE(STATIC)
        do i=1,reception(icpu,levelmin)%ngrid
           igrid=reception(icpu,levelmin)%igrid(i)
           grid_particle_budget(igrid)=max(0,numbp(igrid))
        end do
        !$OMP END PARALLEL DO
     end do
     do ilevel=levelmin+1,nlevelmax
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(jgrid,igrid,parent_cell,parent_grid,parent_ind) &
        !$OMP SCHEDULE(STATIC)
        do jgrid=1,active(ilevel)%ngrid
           igrid=active(ilevel)%igrid(jgrid)
           parent_cell=father(igrid)
           if(parent_cell>ncoarse)then
              parent_grid=IGRID_OF(parent_cell)
              parent_ind=ICHILD_OF(parent_cell)
              grid_particle_budget(igrid)= &
                   grid_particle_budget(parent_grid)/twotondim
              if(parent_ind<=mod(grid_particle_budget(parent_grid), &
                   twotondim))grid_particle_budget(igrid)= &
                   grid_particle_budget(igrid)+1
           end if
        end do
        !$OMP END PARALLEL DO
        do icpu=1,ncpu
           !$OMP PARALLEL DO DEFAULT(SHARED) &
           !$OMP PRIVATE(jgrid,igrid,parent_cell,parent_grid,parent_ind) &
           !$OMP SCHEDULE(STATIC)
           do jgrid=1,reception(icpu,ilevel)%ngrid
              igrid=reception(icpu,ilevel)%igrid(jgrid)
              parent_cell=father(igrid)
              if(parent_cell>ncoarse)then
                 parent_grid=IGRID_OF(parent_cell)
                 parent_ind=ICHILD_OF(parent_cell)
                 grid_particle_budget(igrid)= &
                      grid_particle_budget(parent_grid)/twotondim
                 if(parent_ind<=mod(grid_particle_budget(parent_grid), &
                      twotondim))grid_particle_budget(igrid)= &
                      grid_particle_budget(igrid)+1
              end if
           end do
           !$OMP END PARALLEL DO
        end do
     end do
  end if

  !----------------------------------------
  ! Compute cell ordering and cost
  ! for leaf cells with cpu map = myid.
  ! Store cost in flag1 and MAXIMUM
  ! ordering key in hilbert_key of kind=16
  !----------------------------------------
  ncell=0
  npart_sub=0
  ncell_sub=0
  ncell_loc=1
  dx=1.0*scale
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     if(cpu_map(ind)==myid.and.son(ind)==0)then
        xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
        xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
        xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
        call cmp_minmaxorder(xx,order_min,order_max,dx,ncell_loc)
        if(overload>1)then
           call cmp_dommap(xx,dom,ncell_loc)
        else
           dom(1)=1
        end if
        ncell=ncell+1
        isub=(dom(1)-1)/ncpu+1
        ncell_sub(isub)=ncell_sub(isub)+1
        wflag=domain_leaf_cost(0,0,1_8,level_mesh_scale_ema(levelmin))
        if(wflag>huge(flag1(ncell)))then
           write(*,*)'load_balance: coarse leaf cost exceeds flag1 range: ',wflag
           stop
        end if
        flag1(ncell)=int(wflag)
        hilbert_key(ncell)=order_max(1)
     end if
  end do
  end do
  end do
  ! Loop over levels (OMP parallelized on igrid loop)
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(igrid,ngrid,ind,idim,ilevel,ncell_loc,batch_size,my_base,my_idx, &
  !$OMP         ind_grid,ind_cell,xx,order_min,order_max,dom,isub,wflag, &
  !$OMP         ncell_sub_t,npart_sub_t,npart_leaf,ndm_leaf,npair_cell,i, &
  !$OMP         tt_work) &
  !$OMP REDUCTION(+:fast_grid_loc,fast_fallback_loc,exact_particle_visits_loc, &
  !$OMP             particle_assigned_loc,tpart_work_loc,tkey_work_loc,tcost_work_loc)
  ncell_sub_t=0
  npart_sub_t=0
  do ilevel=1,nlevelmax
     ! Cell size and cell center offset
     !$OMP SINGLE
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do
     !$OMP END SINGLE
     ! Only active(myid) grids contribute leaf cells with cpu_map==myid.
     ! Reception(icpu,ilevel)%igrid level-ilevel cells have cpu_map=icpu
     ! by build_comm invariant, so the filter never matches for icpu/=myid.
     ncache=active(ilevel)%ngrid
     ! Loop over grids by vector sweeps (OMP workshared)
     !$OMP DO SCHEDULE(DYNAMIC,4)
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        npart_leaf(1:ngrid,:)=0
        ndm_leaf(1:ngrid,:)=0
        tt_work=omp_get_wtime()
        if(pic)then
           if(use_fast_particle_balance)then
              do i=1,ngrid
                 call distribute_particle_total_by_leaf( &
                      int(grid_particle_budget(ind_grid(i)),kind=8), &
                      all_children,npart_leaf(i,:))
                 fast_grid_loc=fast_grid_loc+1_8
              end do
           else
              do i=1,ngrid
                 if(numbp(ind_grid(i))>0)then
                    call count_particles_by_leaf(ind_grid(i), &
                         npart_leaf(i,:),ndm_leaf(i,:), &
                         count_dm=need_dm_count)
                    exact_particle_visits_loc=exact_particle_visits_loc+ &
                         int(numbp(ind_grid(i)),kind=8)
                 end if
              end do
           end if
        end if
        tpart_work_loc=tpart_work_loc+omp_get_wtime()-tt_work
        ! Loop over cells
        do ind=1,twotondim
           tt_work=omp_get_wtime()
           do i=1,ngrid
              ind_cell(i)=ICELL_OF(ind_grid(i),ind)
           end do
           do idim=1,ndim
           ncell_loc=0
           do i=1,ngrid
           if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
              ncell_loc=ncell_loc+1
              xx(ncell_loc,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
           end if
           end do
           end do
           if(ncell_loc>0)then
              call cmp_minmaxorder(xx,order_min,order_max,dx*scale,ncell_loc)
              if(overload>1)then
                 call cmp_dommap(xx,dom,ncell_loc)
              else
                 dom(1:ncell_loc)=1
              end if
           end if
           ! Reserve batch of indices atomically
           batch_size=ncell_loc
           if(batch_size>0)then
              !$OMP ATOMIC CAPTURE
              my_base=ncell
              ncell=ncell+batch_size
              !$OMP END ATOMIC
           end if
           tkey_work_loc=tkey_work_loc+omp_get_wtime()-tt_work
           tt_work=omp_get_wtime()
           ncell_loc=0
           do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                 ncell_loc=ncell_loc+1
                 my_idx=my_base+ncell_loc
                 isub=(dom(ncell_loc)-1)/ncpu+1
                 ncell_sub_t(isub)=ncell_sub_t(isub)+1
                 particle_assigned_loc=particle_assigned_loc+ &
                      int(npart_leaf(i,ind),kind=8)
                 npair_cell=domain_sidm_pair_count(ndm_leaf(i,ind))
                 wflag=domain_leaf_cost(npart_leaf(i,ind),npair_cell, &
                      niter_cost(ilevel),level_mesh_scale_ema(ilevel))
                 if(allocated(sink_per_grid))then
                    wflag=wflag+int(sink_per_grid(ind_grid(i)),kind=8)* &
                         int(mem_weight_sink,kind=8)/int(twotondim,kind=8)
                 endif
                 if(wflag>huge(flag1(my_idx)))then
                    write(*,*)'load_balance: leaf cost exceeds flag1 range: ',wflag
                    stop
                 endif
                 flag1(my_idx)=int(wflag)
                 npart_sub_t(isub)=npart_sub_t(isub)+flag1(my_idx)
                 hilbert_key(my_idx)=order_max(ncell_loc)
              end if
           end do
           tcost_work_loc=tcost_work_loc+omp_get_wtime()-tt_work
        end do
        ! End loop over cells
     end do
     !$OMP END DO
     ! End loop over grids
  end do
  ! End loop over levels
  ! Merge thread-local accumulators
  !$OMP CRITICAL(merge_remap_cost)
  do isub=1,overload
     ncell_sub(isub)=ncell_sub(isub)+ncell_sub_t(isub)
     npart_sub(isub)=npart_sub(isub)+npart_sub_t(isub)
  end do
  !$OMP END CRITICAL(merge_remap_cost)
  !$OMP END PARALLEL

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(fast_grid_loc,fast_grid_tot,1,MPI_INTEGER8,MPI_SUM, &
       MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(fast_fallback_loc,fast_fallback_tot,1,MPI_INTEGER8, &
       MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(exact_particle_visits_loc,exact_particle_visits_tot,1, &
       MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(particle_assigned_loc,particle_assigned_tot,1, &
       MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(particle_expected_loc,particle_expected_tot,1, &
       MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(particle_physical_loc,particle_physical_tot,1, &
       MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tpart_work_loc,tpart_work_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tkey_work_loc,tkey_work_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tcost_work_loc,tcost_work_max,1, &
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#else
  fast_grid_tot=fast_grid_loc
  fast_fallback_tot=fast_fallback_loc
  exact_particle_visits_tot=exact_particle_visits_loc
  particle_assigned_tot=particle_assigned_loc
  particle_expected_tot=particle_expected_loc
  particle_physical_tot=particle_physical_loc
  tpart_work_max=tpart_work_loc
  tkey_work_max=tkey_work_loc
  tcost_work_max=tcost_work_loc
#endif
  if(use_fast_particle_balance.and.pic)then
     particle_conservation_error=abs(particle_assigned_tot- &
          particle_expected_tot)
     ! A rank-boundary branch can exist only as a reception oct on the rank
     ! carrying its parent remainder.  Permit at most one 2^ndim remainder per
     ! rank; anything larger indicates a broken propagation and remains fatal.
     particle_conservation_tol=max(1_8,int(ncpu,kind=8)* &
          int(twotondim,kind=8))
     if(particle_conservation_error>particle_conservation_tol)then
        if(myid==1)write(*,'(A,I0,A,I0,A,I0,A,I0)') &
             ' ERROR fast particle balance conservation: assigned=', &
             particle_assigned_tot,' tree=',particle_expected_tot, &
             ' delta=',particle_conservation_error,' tolerance=', &
             particle_conservation_tol
        call clean_stop
     else if(myid==1.and.particle_conservation_error>0_8)then
        write(*,'(A,I0,A,I0)') &
             ' Fast particle balance boundary remainder delta=', &
             particle_conservation_error,' tolerance=',particle_conservation_tol
     end if
  end if
  if(myid==1.and.pic)then
     write(*,'(A,3(F10.3,A))') &
          ' cmp key/cost max-rank thread-work: particles=',tpart_work_max, &
          ' s keys=',tkey_work_max,' s costs=',tcost_work_max,' s'
     if(use_fast_particle_balance)then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
             ' cmp fast particle grids=',fast_grid_tot, &
             ' fallback=',fast_fallback_tot,' assigned=', &
             particle_assigned_tot,' tree=',particle_expected_tot, &
             ' physical=',particle_physical_tot
     else
        write(*,'(A,I0)') &
             ' cmp exact particle linked-list visits=',exact_particle_visits_tot
     end if
  end if

#ifndef WITHOUTMPI
  tcmp_key=MPI_WTIME()
#endif

  ! Clean up sink cost array
  if(allocated(sink_per_grid)) deallocate(sink_per_grid)
  if(allocated(grid_particle_budget)) deallocate(grid_particle_budget)

  ! Reset time-based load balancing accumulators
  level_time_loc = 0d0
  level_ncells_loc = 0

  !------------------------------------------------
  ! Sort ordering key and store new index in flag2
  !------------------------------------------------
  if (ncell>0) call quick_sort_omp(hilbert_key(1),flag2(1),ncell)

#ifndef WITHOUTMPI
  tcmp_sort=MPI_WTIME()
#endif

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc = 0 ! Compute local and global cost
  do isub=1,overload
     cost_loc(myid+(isub-1)*ncpu) = dble(npart_sub(isub))
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ndomain,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  incost_tot = 0D0
  incost_old(0) = 0D0
  do idom = 1,ndomain
     incost_tot = incost_tot + cost_old(idom)
     incost_old(idom) = incost_tot
  end do

  ! A cost-only partition can put too many cheap cells in one domain.  Impose
  ! a lower bound on every positive cost so that the equal domain budget
  ! cannot represent more than grid_cap entries.  flag1 carries one entry per
  ! leaf cell, not per grid, so the per-rank capacity ngridmax is converted to
  ! cells by twotondim before the cap is formed.  Zero-cost entries are not
  ! AMR cells and retain their historical zero cost.
  guard_applied = .false.
  if(lb_grid_headroom > 0d0 .and. ngridmax > 0 .and. incost_tot > 0d0) then

     ! Grid slots already spent on ghost copies of neighbouring domains and on
     ! physical boundaries.  A rank must fit its own grids into whatever is
     ! left.  The ghost layer follows the surface of the partition rather than
     ! its volume, so no cell budget can predict it; measure what the previous
     ! balance actually produced and subtract that, which lets each balance
     ! correct the one before it.  Job 399652 died precisely here: the guard
     ! predicted 6.8M grids against a 9.3M cap and the rank still ran out,
     ! because the ghost layer alone needed another 40 per cent on top.
     ngrid_own_loc = 0_8
     ngrid_ext_loc = 0_8
     do ilev_g = 1,nlevelmax
        ngrid_own_loc = ngrid_own_loc+int(numbl(myid,ilev_g),kind=8)
        do icpu_g = 1,ncpu
           if(icpu_g /= myid) &
                & ngrid_ext_loc = ngrid_ext_loc+int(numbl(icpu_g,ilev_g),kind=8)
        end do
        do ibnd_g = 1,nboundary
           ngrid_ext_loc = ngrid_ext_loc+int(numbb(ibnd_g,ilev_g),kind=8)
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(ngrid_own_loc,ngrid_own_max,1, &
          & MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ngrid_ext_loc,ngrid_ext_max,1, &
          & MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,info)
#else
     ngrid_own_max = ngrid_own_loc
     ngrid_ext_max = ngrid_ext_loc
#endif

     grid_avail = lb_grid_headroom*dble(ngridmax)-dble(ngrid_ext_max)
     ! Never let a huge ghost measurement drive the budget to nothing; the
     ! fallback below is a better answer than an unsatisfiable target.
     grid_avail = max(grid_avail,0.1d0*lb_grid_headroom*dble(ngridmax))
     grid_cap = LB_LEAF_PER_GRID*grid_avail/dble(overload)
     guard_denom = dble(ndomain)*grid_cap

     if(myid==1)then
        write(*,'(A,I0,A,I0,A,I0,A,F6.1,A)') &
             ' LB grid usage: own=',ngrid_own_max,' ghost+bnd=',ngrid_ext_max, &
             ' ngridmax=',ngridmax,' ghost share=', &
             1d2*dble(ngrid_ext_max)/dble(max(ngrid_own_max+ngrid_ext_max,1_8)),'%'
     end if

     min_weight_loc = huge(1d0)
     ntot_grids_loc = 0_8
     do i=1,ncell
        if(flag1(i) > 0) then
           min_weight_loc = min(min_weight_loc,dble(flag1(i)))
           ntot_grids_loc = ntot_grids_loc+1_8
        end if
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(min_weight_loc,min_weight_global,1, &
          & MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(ntot_grids_loc,ntot_grids,1, &
          & MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#else
     min_weight_global = min_weight_loc
     ntot_grids = ntot_grids_loc
#endif

     floor_w = incost_tot/guard_denom
     if(ntot_grids > 0_8 .and. floor_w > min_weight_global) then
        do guard_iter=1,LB_GRID_GUARD_MAXITER
           floor_w = incost_tot/guard_denom
           if(floor_w <= min_weight_global) exit

           ! Close the fixed-point iteration conservatively on the last pass.
           ! For current total T, K positive entries and denominator D, choosing
           ! floor >= T/(D-K) guarantees T_new <= D*floor.
           if(guard_iter == LB_GRID_GUARD_MAXITER) then
              if(guard_denom <= dble(ntot_grids)) then
                 if(myid==1) write(*,*) &
                      ' LB grid guard: occupancy exceeds the headroom;', &
                      ' falling back to pure count balancing'
                 ! The headroom target is out of reach, but an equal-count
                 ! split is still the best partition for the count limit and
                 ! remains feasible whenever occupancy stays under ngridmax.
                 ! Give every positive entry the same cost and stop iterating.
                 do i=1,ncell
                    if(flag1(i) > 0) flag1(i) = 1
                 end do
                 guard_applied = .true.
                 min_weight_global = 1d0
                 exit
              end if
              floor_w = max(floor_w,incost_tot/ &
                   & (guard_denom-dble(ntot_grids)))
           end if

           if(floor_w > dble(huge(flag1(1)))) then
              if(myid==1) write(*,*) &
                   ' wrong type for flag1 --> change to integer kind=8: floor=',floor_w
#ifndef WITHOUTMPI
              call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif
              stop
           end if
           floor_flag = ceiling(floor_w)

           nraised_loc = 0_8
           do i=1,ncell
              if(flag1(i) > 0 .and. flag1(i) < floor_flag) then
                 flag1(i) = floor_flag
                 nraised_loc = nraised_loc+1_8
              end if
           end do
#ifndef WITHOUTMPI
           call MPI_ALLREDUCE(nraised_loc,nraised,1, &
                & MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#else
           nraised = nraised_loc
#endif
           guard_applied = guard_applied .or. (nraised > 0_8)
           min_weight_global = max(min_weight_global,dble(floor_flag))

           ! Rebuild every local subdomain sum because flag1 changed.
           npart_sub = 0_8
           ncell_loc = 0
           do isub=1,overload
              do i=1,ncell_sub(isub)
                 npart_sub(isub) = npart_sub(isub) + &
                      & int(flag1(flag2(ncell_loc+i)),kind=8)
              end do
              ncell_loc = ncell_loc+ncell_sub(isub)
           end do
           cost_loc = 0d0
           do isub=1,overload
              cost_loc(myid+(isub-1)*ncpu) = dble(npart_sub(isub))
           end do
#ifndef WITHOUTMPI
           call MPI_ALLREDUCE(cost_loc,cost_old,ndomain, &
                & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
           cost_old = cost_loc
#endif
           incost_tot = 0d0
           incost_old(0) = 0d0
           do idom=1,ndomain
              incost_tot = incost_tot+cost_old(idom)
              incost_old(idom) = incost_tot
           end do
        end do

        predicted_maxcount = (incost_tot/dble(ndomain))/min_weight_global
        cost_imbalance = maxval(cost_old)/(incost_tot/dble(ndomain))
        if(myid==1 .and. guard_applied) then
           write(*,'(A,ES12.4,A,I0,A,I0,A,F14.0,A,F14.0,A,F10.4)') &
                ' LB grid guard: floor=',floor_w,' raised=',nraised,'/',ntot_grids, &
                ' predicted max cells=',predicted_maxcount,' / cap=',grid_cap, &
                ' cost imbalance=',cost_imbalance
           write(*,'(A,F14.0,A,F14.0,A,I0)') &
                '                in grids: predicted own=', &
                predicted_maxcount/LB_LEAF_PER_GRID,' + measured ghost=', &
                dble(ngrid_ext_max),' vs ngridmax=',ngridmax
        end if
     else if(myid==1 .and. verbose) then
        write(*,'(A)') ' LB grid guard: not needed'
     end if
  end if

  incost_new(0) = 0D0
  do idom = 1,ndomain
     cost_new(idom) = incost_tot/dble(ndomain) ! Exact load balancing
     incost_new(idom) = incost_new(idom-1) + cost_new(idom)
  end do

  !-----------------------------
  ! Compute new cpu boundaries
  !-----------------------------
  bound_key_loc=0.0d0; bound_key2=0.0d0
  ncell_loc=0
  do isub=1,overload
     if(ncell_sub(isub)>0)then
        ! First cpu on local domain
        idom=0
        do while(incost_new(idom)<incost_old(myid-1+(isub-1)*ncpu))
           idom=idom+1
           if (idom > ndomain) exit 
        end do
        ! Compute Hilbert key at boundaries
        i=idom
        local_cost=incost_old(myid-1+(isub-1)*ncpu)
        do ind_long=1,ncell_sub(isub)
           cell_cost=dble(flag1(flag2(ind_long+ncell_loc)))
           local_cost=local_cost+cell_cost
           if (i > ndomain) exit
           if(incost_new(i)<local_cost)then
              bound_key_loc(i)=hilbert_key(ind_long+ncell_loc)
              i=i+1
           endif
        end do
     end if
     ncell_loc=ncell_loc+ncell_sub(isub)
  end do
#ifndef WITHOUTMPI
#ifdef QUADHILBERT
  bigdbl= real(bound_key_loc,kind=8)
  bigtmp= 0.0d0
  call MPI_ALLREDUCE(bigdbl,bigtmp,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ! if call to mpi_sum with mpi_type=mpi_real16 is supported by mpi_allreduce we can do: 
  !call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_REAL16,MPI_SUM,MPI_COMM_WORLD,info)
  bound_key2         = real(bigtmp,kind=qdp)
#else
  call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#endif
  bound_key2(0)      =order_all_min
  bound_key2(ndomain)=order_all_max

  else if(ordering=='bisection') then
     ! update the bisection
     call build_bisection(update=.true.)
  else  ! ordering=='ksection'
     ! update the ksection
     call build_ksection(update=.true.)
  end if   ! end if not bisection/ksection

#ifndef WITHOUTMPI
  tcmp_bound=MPI_WTIME()
#endif

  ! Reset time-based load balancing accumulators (ksection/bisection path)
  if(use_cpubox_decomp) then
     level_time_loc = 0d0
     level_ncells_loc = 0
  end if

  ! Free on-demand histogram arrays (no longer needed after build_bisection/ksection)
  if(allocated(bisec_ind_cell))  deallocate(bisec_ind_cell)
  if(allocated(bisec_cell_level))deallocate(bisec_cell_level)
  if(allocated(bisec_cell_coord))deallocate(bisec_cell_coord)
  if(allocated(bisec_cell_cost)) deallocate(bisec_cell_cost)

  !----------------------------------------
  ! Compute new cpu map
  !----------------------------------------
210 continue
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
  do i=1,size(cpu_map2)
     cpu_map2(i)=0
  end do
  !$OMP END PARALLEL DO
  ncell_loc=1
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
     xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
     xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
     cpu_map2(ind)=ncpu ! default value                                                               

     if(.not.use_cpubox_decomp) then
        call cmp_ordering(xx,order_max,ncell_loc)
        cpu_map2(ind)=ncpu ! default value
        ilo=1
        ihi=ndomain
        do while(ilo<ihi)
           imid=(ilo+ihi)/2
           if(order_max(1)<bound_key2(imid))then
              ihi=imid
           else
              ilo=imid+1
           end if
        end do
        cpu_map2(ind)=mod(ilo-1,ncpu)+1
     else if(ordering=='bisection') then
        xx_tmp(1,:) = xx(1,:)
        call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
        cpu_map2(ind) = c_tmp(1)
     else  ! ksection
        xx_tmp(1,:) = xx(1,:)
        call cmp_ksection_cpumap(xx_tmp,c_tmp,1)
        cpu_map2(ind) = c_tmp(1)
     end if
  end do
  end do
  end do
  ! Loop over levels (OMP parallelized on igrid loop)
  do ilevel=1,nlevelmax
     ! Cell size and cell center offset
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do
     ncache=active(ilevel)%ngrid
     ! Loop over grids by vector sweeps (OMP parallelized)
     !$OMP PARALLEL DO DEFAULT(SHARED) &
     !$OMP PRIVATE(igrid,ngrid,ind_grid,ind_cell,xx,order_max,i,ind,idim,idom, &
     !$OMP         ilo,ihi,imid,xx_tmp,c_tmp,c_tmp_v) &
     !$OMP SCHEDULE(STATIC)
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           do i=1,ngrid
              ind_cell(i)=ICELL_OF(ind_grid(i),ind)
           end do
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
              end do
           end do
           if(.not.use_cpubox_decomp) then
              if(ngrid>0)call cmp_ordering(xx,order_max,ngrid)
              do i=1,ngrid
                 ilo=1
                 ihi=ndomain
                 do while(ilo<ihi)
                    imid=(ilo+ihi)/2
                    if(order_max(i)<bound_key2(imid))then
                       ihi=imid
                    else
                       ilo=imid+1
                    end if
                 end do
                 cpu_map2(ind_cell(i))=mod(ilo-1,ncpu)+1
              end do
           else if(ordering=='bisection') then
              do i=1,ngrid
                 xx_tmp(1,:) = xx(i,:)
                 call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
                 cpu_map2(ind_cell(i)) = c_tmp(1)
              end do
           else  ! ksection
              if(ngrid>0) call cmp_ksection_cpumap(xx,c_tmp_v,ngrid)
              do i=1,ngrid
                 cpu_map2(ind_cell(i)) = c_tmp_v(i)
              end do
           endif
        end do
        ! End loop over cells
     end do
     !$OMP END PARALLEL DO
     ! End loop over grids
  end do
  ! End loop over levels

  ! Bound a large ordered-domain move by the grid slots that are free now.
  ! This leaves the existing refine/build_comm/Morton machinery untouched:
  ! each pass performs a normal remap, shrinks the old boundary, and only then
  ! advances farther toward the requested boundary in the next pass.
  !
  ! Count complete source grids rather than particles or leaf cells.  The
  ! incoming count is a conservative upper bound because some incoming grids
  ! may already exist locally as virtual grids.  Half of the headroom is kept
  ! for the new ghost surface and physical boundaries.
  ! Interpolating bound_key is a Hilbert-only operation.  Other orderings use
  ! their own domain descriptors and must retain their normal remap path.
  if(trim(ordering)=='hilbert' .and. lb_remap_fraction>0d0)then
     lb_nsend=0
     do ilevel=1,nlevelmax
        do i=1,active(ilevel)%ngrid
           igrid=active(ilevel)%igrid(i)
           ! Before the expand pass, a small number of coarse/root grids can
           ! still have no parent cell.  The normal remap assigns those
           ! fathers during refine_coarse; they must not index cpu_map2(0)
           ! in this preflight estimate.  Omitting them is harmless because
           ! half of the available slots is already reserved as margin.
           if(father(igrid)<=0) cycle
           lb_target_cpu=cpu_map2(father(igrid))
           if(lb_target_cpu<1.or.lb_target_cpu>ncpu) cycle
           if(lb_target_cpu/=myid) lb_nsend(lb_target_cpu)= &
                lb_nsend(lb_target_cpu)+1
        end do
     end do
#ifndef WITHOUTMPI
     call MPI_ALLTOALL(lb_nsend,1,MPI_INTEGER,lb_nrecv,1,MPI_INTEGER, &
          MPI_COMM_WORLD,info)
#else
     lb_nrecv=lb_nsend
#endif
     lb_incoming=sum(lb_nrecv)-lb_nrecv(myid)
     lb_free_slots=max(0,int(lb_grid_headroom*dble(ngridmax))-used_mem)
     lb_allowed=max(0,int(LB_REMAP_SLOT_SHARE*dble(lb_free_slots)))
     if(lb_incoming>0)then
        lb_fraction_local=min(1d0,dble(lb_allowed)/dble(lb_incoming))
     else
        lb_fraction_local=1d0
     end if
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(lb_fraction_local,lb_fraction_global,1, &
          MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#else
     lb_fraction_global=lb_fraction_local
#endif
     if(lb_fraction_global<0.999999d0)then
        if(.not.lb_boundary_limited)then
           lb_bound_key_target=bound_key2
           lb_boundary_limited=.true.
        endif
        ! Hilbert occupancy is not linear in boundary-key distance.  Recount
        ! the actual candidate after every reduction rather than assuming a
        ! fraction of the full key motion moves the same fraction of grids.
        lb_remap_fraction=lb_remap_fraction*lb_fraction_global
        if(lb_remap_fraction<0.01d0 .or. lb_limit_iter>=12)then
           lb_remap_fraction=0d0
           bound_key2=bound_key
           if(myid==1) write(*,'(A,I0,A,I0,A)') &
                ' Bounded remap deferred: incoming=',lb_incoming, &
                ' allowed=',lb_allowed,' (insufficient free slots)'
        else
           lb_limit_iter=lb_limit_iter+1
           do idom=1,ndomain-1
              bound_key2(idom)=bound_key(idom)+ &
                   real(lb_remap_fraction,kind=qdp)* &
                   (lb_bound_key_target(idom)-bound_key(idom))
           end do
           bound_key2(0)=order_all_min
           bound_key2(ndomain)=order_all_max
           if(myid==1) write(*,'(A,I0,A,F7.3,A,I0,A,I0,A,I0)') &
                ' Bounded remap candidate ',lb_limit_iter, &
                ' fraction=',lb_remap_fraction, &
                ' incoming=',lb_incoming,' allowed=',lb_allowed, &
                ' free=',lb_free_slots
        end if
        goto 210
     else if(lb_boundary_limited .and. myid==1)then
        write(*,'(A,F7.3,A,I0,A,I0,A,I0)') &
             ' Bounded remap accepted fraction=',lb_remap_fraction, &
             ' actual incoming=',lb_incoming,' allowed=',lb_allowed, &
             ' free=',lb_free_slots
     end if
  end if

#ifndef WITHOUTMPI
  tcmp_map=MPI_WTIME()
#endif

  ! Update virtual boundaries for new cpu map
  call make_virtual_coarse_int(cpu_map2(1))
  do ilevel=1,nlevelmax
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do 

#ifndef WITHOUTMPI
  tcmp_virtual=MPI_WTIME()
  if(myid==1)then
     write(*,'(A,A,A,5(F9.3,A))') ' cmp_new_cpu_map stages [',trim(ordering), &
          ']: key/cost=',tcmp_key-tcmp_start,' s sort=',tcmp_sort-tcmp_key, &
          ' s boundary=',tcmp_bound-tcmp_sort,' s cpumap=',tcmp_map-tcmp_bound, &
          ' s virtual=',tcmp_virtual-tcmp_map,' s'
  end if
#endif

end subroutine cmp_new_cpu_map
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap(x,c,nn)
  use amr_parameters
  use amr_commons
  use bisection
  use ksection
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ilo,ihi,imid
  real(qdp),dimension(1:nvector)::order

  if(.not.use_cpubox_decomp) then
     call cmp_ordering(x,order,nn)
     do i=1,nn
        ilo=1
        ihi=ndomain
        do while(ilo<ihi)
           imid=(ilo+ihi)/2
           if(order(i)<bound_key(imid))then
              ihi=imid
           else
              ilo=imid+1
           end if
        end do
        c(i)=ilo
     end do
     do i=1,nn
        c(i)=MOD(c(i)-1,ncpu)+1
     end do
  else if(ordering=='bisection') then
     call cmp_bisection_cpumap(x,c,nn)
  else  ! ksection
     call cmp_ksection_cpumap(x,c,nn)
  end if

end subroutine cmp_cpumap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_dommap(x,c,nn)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ilo,ihi,imid
  real(qdp),dimension(1:nvector)::order

  call cmp_ordering(x,order,nn)
  do i=1,nn
     ilo=1
     ihi=ndomain
     do while(ilo<ihi)
        imid=(ilo+ihi)/2
        if(order(i)<bound_key(imid))then
           ihi=imid
        else
           ilo=imid+1
        end if
     end do
     c(i)=ilo
  end do
  
end subroutine cmp_dommap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ordering(x,order,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer ::nn
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order
  !--------------------------------------------------------
  ! This routine computes the index key of the input cell
  ! according to its position in space and for the chosen
  ! ordering. Position x are in user units.
  !-----------------------------------------------------
  integer,dimension(1:nvector)::ix,iy,iz
  integer::i,ncode,bit_length,nx_loc
  integer::temp,info
  real(kind=8)::scale,bscale,xx,yy,zz,xc,yc,zc

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  if(ordering=='planar')then
     ! Planar domain decomposition
     do i=1,nn
        order(i)=x(i,1)
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     xc=boxlen/2.
     yc=boxlen/2.
     zc=boxlen/2.
     do i=1,nn
        xx=x(i,1)-xc+1d-10
        yy=x(i,2)-yc
#if NDIM>2
        zz=x(i,3)
#endif
        if(xx>0.)then
           order(i)=atan(yy/xx)+acos(-1.)/2.
        else
           order(i)=atan(yy/xx)+acos(-1.)*3./2.
        endif
#if NDIM>2
        if(zz.gt.zc)order(i)=order(i)+2.*acos(-1.)
#endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     ncode=nx_loc*int(bscale)
     bscale=bscale/scale
     
     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

     do i=1,nn
        ix(i)=int(x(i,1)*bscale)
#if NDIM>1           
        iy(i)=int(x(i,2)*bscale)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscale)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order,bit_length,nn)
     end if

  end if

end subroutine cmp_ordering
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_minmaxorder(x,order_min,order_max,dx,nn)
  use amr_parameters
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer ::nn
  integer ::temp,info
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order_min,order_max
  !-----------------------------------------------------
  ! This routine computes the minimum and maximum index
  ! key contained in the input cell and for the chosen 
  ! ordering.
  !-----------------------------------------------------
  integer,dimension(1:nvector)::ix,iy,iz
  integer::i,ncode,bit_length,nxny,nx_loc

  real(dp)::theta1,theta2,theta3,theta4,dxx,dxmin  
  real(kind=8)::scale,bscaleloc,bscale,xx,yy,zz,xc,yc,zc
  real(qdp)::dkey,oneqdp=1.0

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dxmin=scale/dble(2**nlevelmax)

  if(ordering=='planar')then
     ! Planar domain decomposition
     dxx=0.5d0*dx
     do i=1,nn
        order_min(i)=x(i,1)-dxx
        order_max(i)=x(i,1)+dxx
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     dxx=0.5d0*dx
     xc=boxlen/2.
     yc=boxlen/2.
     zc=boxlen/2.
     do i=1,nn
        if(dx==boxlen)then
           order_min(i)=0.
           order_max(i)=4.*acos(-1.)
        else
           ! x- y-
           yy=x(i,2)-yc-dxx
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta1=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta1=atan(yy/xx)+acos(-1.)*3./2.
           endif
           ! x+ y-
           xx=x(i,1)-xc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta2=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta2=atan(yy/xx)+acos(-1.)*3./2.
           endif
           
           ! x+ y+
           yy=x(i,2)-yc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta3=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta3=atan(yy/xx)+acos(-1.)*3./2.
           endif
           ! x- y+
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta4=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta4=atan(yy/xx)+acos(-1.)*3./2.
           endif
           order_min(i)=min(theta1,theta2,theta3,theta4)
           order_max(i)=max(theta1,theta2,theta3,theta4)
#if NDIM>2
           zz=x(i,3)
           if(zz.gt.zc)then
              order_min(i)=order_min(i)+2.*acos(-1.)
              order_max(i)=order_max(i)+2.*acos(-1.)
           endif
#endif
        endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     bscaleloc=2**nlevelmax*dxmin/dx
     ncode=nx_loc*int(bscaleloc)
     bscaleloc=bscaleloc/scale
     bscale   =bscale   /scale
     
     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

     do i=1,nn
        ix(i)=int(x(i,1)*bscaleloc)
#if NDIM>1           
        iy(i)=int(x(i,2)*bscaleloc)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscaleloc)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order_min,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order_min,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order_min,bit_length,nn)
     end if

     dkey=(real(bscale,kind=qdp)/real(bscaleloc,kind=qdp))**ndim
     do i=1,nn
        order_max(i)=(order_min(i)+oneqdp)*dkey
        order_min(i)=(order_min(i))*dkey
     end do

  end if

end subroutine cmp_minmaxorder
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine defrag
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ncache,ngrid2,igridmax,i,igrid,ibound,ilevel
  integer::igrid1,igrid2,ind1,ind2,icell1,icell2
  integer::ind,idim,ivar,istart
  real(dp),allocatable::defrag_dp(:)
  integer,allocatable::defrag_map(:)
#ifndef WITHOUTMPI
  real(dp)::t_defrag_start,t_defrag_end
#endif

#ifndef WITHOUTMPI
  t_defrag_start = MPI_WTIME()
#endif
  if(verbose)write(*,*)'Defragmenting main memory...'

  ngrid2=0
  igridmax=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              igridmax=max(igridmax,igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do

  ! Allocate local scratch for old→new grid index mapping (replaces cpu_map2 in defrag)
  allocate(defrag_map(1:igridmax))
  defrag_map=0

  ngrid2=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_map(igrid)=ngrid2+i
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              icell1=father(igrid)
              if(icell1>ncoarse)then
                 ind1=ICHILD_OF(icell1)
                 igrid1=IGRID_OF(icell1)
                 igrid2=defrag_map(igrid1)
                 icell2=ICELL_OF(igrid2,ind1)
              else
                 icell2=icell1
              end if
              flag2(ngrid2+i)=icell2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     father(igrid)=flag2(igrid)
  end do

  ! nbor defrag remapping removed — computed from Morton keys

  ! Allocate local scratch for defrag (replaces hilbert_key usage)
  allocate(defrag_dp(1:igridmax))

  do idim=1,ndim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=xg(igrid,idim)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     xg(igrid,idim)=defrag_dp(igrid)
  end do
  end do

  if(pic)then

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=headp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     headp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=tailp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     tailp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=numbp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     numbp(igrid)=flag2(igrid)
  end do

  endif

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              igrid1=son(ICELL_OF(igrid,ind))
              if(igrid1>0)then
                 igrid2=defrag_map(igrid1)
              else
                 igrid2=0
              end if
              flag2(ngrid2+i)=igrid2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     son(ICELL_OF(igrid,ind))=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=cpu_map(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     cpu_map(ICELL_OF(igrid,ind))=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=flag1(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     flag1(ICELL_OF(igrid,ind))=flag2(igrid)
  end do
  end do

  if(hydro)then

#ifdef SOLVERmhd
  do ivar=1,nvar+3
#else
  do ivar=1,nvar
#endif
  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=uold(ICELL_OF(igrid,ind),ivar)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     uold(ICELL_OF(igrid,ind),ivar)=defrag_dp(igrid)
  end do
  end do
  end do

  end if

#ifdef RT
  if(rt)then

  do ivar=1,nrtvar
  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=rtuold(ICELL_OF(igrid,ind),ivar)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     rtuold(ICELL_OF(igrid,ind),ivar)=defrag_dp(igrid)
  end do
  end do
  end do

  end if
#endif

  if(poisson)then

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=phi(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     phi(ICELL_OF(igrid,ind))=defrag_dp(igrid)
  end do
  end do

  do idim=1,ndim
  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=f(ICELL_OF(igrid,ind),idim)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     f(ICELL_OF(igrid,ind),idim)=defrag_dp(igrid)
  end do
  end do
  end do

  ! Modified-gravity scalar fields follow the AMR octs through
  ! defragmentation just like phi and f.  This is required both after a
  ! variable-ncpu checkpoint restore and during ordinary load balancing.
  if(allocated(scalar_gr))then
  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=scalar_gr(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     scalar_gr(ICELL_OF(igrid,ind))=defrag_dp(igrid)
  end do
  end do

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=scalar_gr_old(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     scalar_gr_old(ICELL_OF(igrid,ind))=defrag_dp(igrid)
  end do
  end do
  end if

  end if

  if(use_fdm)then

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=psi_re(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     psi_re(ICELL_OF(igrid,ind))=defrag_dp(igrid)
  end do
  end do

  do ind=1,twotondim
  ngrid2=0
  do igrid=1,igridmax
     defrag_dp(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              defrag_dp(ngrid2+i)=psi_im(ICELL_OF(igrid,ind))
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     psi_im(ICELL_OF(igrid,ind))=defrag_dp(igrid)
  end do
  end do

  end if

  deallocate(defrag_dp)
  deallocate(defrag_map)

  ngrid2=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           if(ibound<=ncpu)then
              headl(ibound,ilevel)=ngrid2+1
              taill(ibound,ilevel)=ngrid2+ncache
           else
              headb(ibound-ncpu,ilevel)=ngrid2+1
              tailb(ibound-ncpu,ilevel)=ngrid2+ncache
           end if
           prev(ngrid2+1)=0
           do i=2,ncache
              prev(ngrid2+i)=ngrid2+i-1
           end do
           do i=1,ncache-1
              next(ngrid2+i)=ngrid2+i+1
           end do
           next(ngrid2+ncache)=0
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  headf=ngrid2+1
  tailf=ngridmax
  numbf=ngridmax-ngrid2
  prev(headf)=0
  next(tailf)=0
  do i=ngrid2+2,ngridmax
     prev(i)=i-1
  end do
  do i=ngrid2+1,ngridmax-1
     next(i)=i+1
  end do

  do i=1,nlevelmax
     call build_comm(i)
  end do

  ngrid_current=ngrid2
  ! Cached cell-index maps (notably the distributed FDM FFT pack map) must
  ! never survive grid renumbering, even when local counts and bounds match.
  amr_mesh_epoch=amr_mesh_epoch+1_i8b

#ifndef WITHOUTMPI
  t_defrag_end = MPI_WTIME()
  if(myid==1) write(*,'(A,F8.3,A)') ' defrag total: ', t_defrag_end - t_defrag_start, ' s'
#endif

end subroutine defrag
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
! Wrapper to call glibc malloc_trim(0) from Fortran.
! Returns freed heap pages to OS, reducing RSS after bulk dealloc cycles.
!#########################################################################
subroutine fortran_malloc_trim()
  use iso_c_binding, only: c_int, c_size_t
  implicit none
  interface
     integer(c_int) function malloc_trim(pad) bind(C, name='malloc_trim')
       import :: c_int, c_size_t
       integer(c_size_t), value :: pad
     end function malloc_trim
  end interface
  integer(c_int) :: rc
  rc = malloc_trim(0_c_size_t)
end subroutine fortran_malloc_trim

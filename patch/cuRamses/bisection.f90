! KNOWN BUGS:
!  * modifying flag1 in bisection causes load balancing to malfunction
! POSSIBLE IMPROVEMENTS
!  * cell cost is fixed at 1 for now

module bisection
   use amr_parameters
   use amr_commons
   use pm_commons, only: numbp, xsink, count_particles_by_leaf
   use pm_parameters, only: nsink

   implicit none

contains

   subroutine cmp_bisection_cpumap(x,c,nn)
      ! This routine takes a ndim x nvector array as input, representing nvector space points,
      ! and returns the nvector-array of integers corresponding to the matching CPU id in
      ! the domain decomposition.

      ! Array of input coordinates
      real(dp), intent(in), dimension(:,:) :: x
      ! Array of cpu ids to output
      integer, intent(out), dimension(:)   :: c
      integer, intent(in) :: nn
      
      integer :: p, dir, id, cur, half

!      if(verbose) print *, 'entering cmp_bisection_cpumap'

      ! Loop on input points
      do p=1,nn
         ! Begin splitting along the first coordinate
         dir=1
         ! Go down binary tree starting from root
         cur=bisec_root
         do while(bisec_next(cur,1)>0) ! Keep exploring tree downwards till cur is a leaf
            ! Choose relevant half
            half=1; if(x(p,dir)>bisec_wall(cur))half=2
            ! Next node in the tree, in the matching branch
            cur=bisec_next(cur,half)
            ! Next direction
            dir=dir+1; if(dir>ndim)dir=1
         end do
         ! cur should be a leaf by now
         ! Save point cpu id into the output array
         c(p)=bisec_indx(cur)
      end do
 !     if(verbose) print *, 'done with cmp_bisection_cpumap'
   end subroutine cmp_bisection_cpumap


   ! MAIN BISECTION CREATION/UPDATING ROUTINE
   subroutine build_bisection(update)
#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      logical, intent(in) :: update

      ! Tree-wide variables (needed between levels)
      integer,  dimension(1:nbinodes) :: tmp_imin, tmp_imax
      integer(i8b),  dimension(1:nbinodes) :: tmp_load
      real(dp), dimension(1:nbinodes,1:ndim) :: tmp_bxmin, tmp_bxmax

      ! Level-wide variables (needed within one level)
      logical,  dimension(1:nbileafnodes) :: skip
      real(dp), dimension(1:nbileafnodes) :: u_limit, l_limit
      real(dp), dimension(1:nbileafnodes) :: last_wall, best_wall, best_score, walls

      integer(i8b),dimension(1:nbileafnodes) :: load1, myload, totload
      integer,  dimension(1:nbileafnodes) :: lncpu1, lncpu2

      logical :: all_skip, start_bisec
      integer(i8b) :: load2, mytmp, tottmp
      real(dp) :: scale, nload1, nload2, score
      real(dp) :: mean, var, stdev
      real(dp), dimension(1:ndim) :: xmin, xmax
  
      integer :: nc, dir, i, lvl, ierr, iter
      integer :: lncpu, cpuid, child1, child2
      integer :: cur_levelstart, cur_cell

      scale=boxlen/dble(icoarse_max-icoarse_min+1)

      if(verbose) print *,'entering build_bisection with update = ',update

      ! TREE INIT
      bisec_root=1; bisec_next=0

      tmp_imin=0; tmp_imax=0
      tmp_imin(1)=1; tmp_imax(1)=ncpu
      tmp_bxmin(1,:)=0.0; tmp_bxmax(1,:)=scale

      l_limit(1)=0.0; u_limit(1)=1.0
      cur_levelstart=1; nc=0; dir=1; cur_cell=0

      if(update) then
         ! init histogram
         call init_bisection_histogram

         ! compute cell coordinates for current direction
         call compute_bisec_cell_coords(dir)

         ! build top-level histogram in x dir to get total load for current cpu
         call build_bisection_histogram(0,dir,1)

         ! compute total load for comp. box
         mytmp=bisec_hist(1,bisec_nres)   ! total load for current cpu
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(mytmp,tottmp,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
         tottmp=mytmp
#endif
         tmp_load(1)=tottmp
      end if

      ! Loop through levels
      lvl=0
      level_loop: do
         
         ! Number of cells for the current level
         nc  = 2**lvl

         ! Rebuild histograms if needed, level 0 is already done by init
         if(update .and. lvl>0 .and. lvl<nbilevelmax) then
            call compute_bisec_cell_coords(dir)
            call build_bisection_histogram(lvl,dir,nc)
         end if

         ! WALL-FINDING
         start_bisec=.true.
         skip=.false.

         iter=0
         dichotomy_loop: do ! main dichotomy loop, for levels 0..nbilevelmax-1
            if(lvl==nbilevelmax) exit    ! No need to bisect anything at very last level
            iter=iter+1
            
            all_skip=.true. ! init at true

            ! This loop sets the skip flag in some special cases
            do i=1,nc
               cur_cell = cur_levelstart + (i-1)      ! cell id
               if (tmp_imax(cur_cell)==0) then
                  ! this is an empty slot left by a leaf cell on an upper level
                  skip(i) = .true.
                  cycle
               end if

               ! calc number of left and right cpus
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if (lncpu==1) then
                  ! leaf cell, consider it done
                  skip(i) = .true.
                  cycle
               end if

               lncpu1(i) = lncpu/2
               lncpu2(i) = lncpu - lncpu1(i)

               all_skip = all_skip .and. skip(i)
            end do
            ! Check if dichotomy is over
            if (all_skip) exit

            ! skip flag set, start the init
            ! treat new tree creation separately
            build_from_scratch: if (.not. update) then
            do i=1,nc
               cur_cell = cur_levelstart + (i-1)
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if(skip(i)) cycle

               bisec_wall(cur_cell) = round_to_bisec_res( ( tmp_bxmin(cur_cell,dir)*lncpu2(i) &
                                                         + tmp_bxmax(cur_cell,dir)*lncpu1(i) )/lncpu )
               if(bisec_wall(cur_cell)==tmp_bxmin(cur_cell,dir) .or. bisec_wall(cur_cell)==tmp_bxmax(cur_cell,dir)) then
                  if(myid==1) print *,"Problem in bisection tree creation : insufficient resolution"
#ifndef WITHOUTMPI
                  call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
                  stop
               end if
            end do
            ! don't stay in the dichotomy loop
            exit dichotomy_loop
            end if build_from_scratch


            ! tree update, two cases : 1. init, 2. bisection
            if(start_bisec) then
               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
                  if(skip(i)) cycle
                  ! check whether wall position is compatible with bounding box
                  if( bisec_wall(cur_cell)<=tmp_bxmin(cur_cell,dir) .or. bisec_wall(cur_cell)>=tmp_bxmax(cur_cell,dir) ) then
                     bisec_wall(cur_cell) = round_to_bisec_res( 0.5 * (tmp_bxmin(cur_cell,dir) + &
                                                                          tmp_bxmax(cur_cell,dir) ) )
                  end if

                  ! get local load for left subcell knowing current wall pos
                  myload(i) = bisec_hist( i , floor(bisec_wall(cur_cell)/bisec_res)+1 )
               end do
               ! sum the local left loads from every cpu into load1
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(myload,load1,nbileafnodes,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               load1=myload
#endif

               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle
                  
                  ! init best wall and init score to a terrible value
                  best_wall(i) = bisec_wall(cur_cell)
#if NPRE==4
                  best_score(i) = huge(1.0e0)
#endif
#if NPRE==8
                  best_score(i) = huge(1.0d0)
#endif

                  ! setup domain limits
                  l_limit(i) = tmp_bxmin(cur_cell, dir)
                  u_limit(i) = tmp_bxmax(cur_cell, dir)
               end do
               start_bisec=.false.
            else
               ! not starting new dichotomic stuff
               do i=1,nc
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle

                  ! retrieve differential load
                  myload(i) = abs( bisec_hist( i , floor(max(bisec_wall(cur_cell),last_wall(i))/bisec_res)+1 ) &
                              - bisec_hist( i , floor(min(bisec_wall(cur_cell),last_wall(i))/bisec_res)+1 ) )
               end do

               ! sum up all differential loads into totload
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(myload,totload,nbileafnodes,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               totload=myload
#endif

               do i=1,nc
                  ! transfer differential load
                  cur_cell = cur_levelstart + (i-1)
                  if(skip(i)) cycle
                  
                  if(bisec_wall(cur_cell)>last_wall(i)) then
                     load1(i) = load1(i) + totload(i)
                  else
                     load1(i) = load1(i) - totload(i)
                  end if

               end do
            end if

            do i=1,nc
               cur_cell = cur_levelstart + (i-1)
               lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
               if(skip(i)) cycle

               ! compute imbalance
               load2 = tmp_load(cur_cell)-load1(i)
               nload1 = dble(load1(i))/dble(lncpu1(i)); nload2 = dble(load2)/dble(lncpu2(i))  ! normalized loads
               score = abs(nload1-nload2)/(nload1+nload2)

               ! tolerance met ?
               if(score < bisec_tol) then
                  skip(i)=.true.
                  cycle
               end if

               ! tolerance is not met, proceed with bisection.
               ! if wall is the best one yet, store it
               if (score < best_score(i)) then
                   best_score(i) = score
                   best_wall(i)  = bisec_wall(cur_cell)
               end if
               ! compute new wall position for next bisection step
               if(nload1>nload2) then
                  ! move left
                  u_limit(i) = bisec_wall(cur_cell)
               else
                  ! move right
                  l_limit(i) = bisec_wall(cur_cell)
               end if

               ! wall pos for next step
               last_wall(i) = bisec_wall(cur_cell)
               bisec_wall(cur_cell) = round_to_bisec_res( 0.5 * (u_limit(i) + l_limit(i)) )
               ! check if we're at resolution limit for next bisection step
               if( abs(bisec_wall(cur_cell)-u_limit(i))<0.5*bisec_res &
                       .or. abs(bisec_wall(cur_cell)-l_limit(i))<0.5*bisec_res )  then
                  ! restore best wall, mark cell done and loop
                  bisec_wall(cur_cell) = best_wall(i)
                  skip(i) = .true.
                  cycle
               end if
            end do

         end do dichotomy_loop

         ! CHILDREN CREATION AND LEAF PROCESSING
         ! this is done at every level, including the very last one (for leaf processing)
         walls=0.0
         children_and_leaves: do i=1,nc
            cur_cell = cur_levelstart + (i-1)
            if (tmp_imax(cur_cell)==0) cycle

            ! Is current cell a leaf?
            if (tmp_imin(cur_cell)==tmp_imax(cur_cell)) then
               cpuid=tmp_imin(cur_cell)
               ! save cpu id
               bisec_indx(cur_cell)=cpuid
               ! save cpu bound box.
               ! for a first computation (update=false) this goes into bisec_cpubox_min
               ! for update=true, this goes into bisec_cpubox_min2, as the old boxes are still
               ! needed for load balancing (virtual_boundaries.f90)
               if(update) then
                  bisec_cpubox_min2(cpuid,:)=tmp_bxmin(cur_cell,:)
                  bisec_cpubox_max2(cpuid,:)=tmp_bxmax(cur_cell,:)
               else
                  bisec_cpubox_min(cpuid,:)=tmp_bxmin(cur_cell,:)
                  bisec_cpubox_max(cpuid,:)=tmp_bxmax(cur_cell,:)
               end if
               ! save cpu workload
               bisec_cpu_load(cpuid)=tmp_load(cur_cell)
               ! make sure has no child
               bisec_next(cur_cell,:)=0
               ! skip to next cpu
               cycle
            end if

            ! node ids of the two children
            child1 = cur_levelstart + nc + 2*i-2
            child2 = cur_levelstart + nc + 2*i-1
            ! create links in the tree
            bisec_next(cur_cell,1)=child1
            bisec_next(cur_cell,2)=child2
            ! store bounding boxes
            tmp_bxmin(child1,:)=tmp_bxmin(cur_cell,:)
            tmp_bxmax(child1,:)=tmp_bxmax(cur_cell,:)
            tmp_bxmin(child2,:)=tmp_bxmin(cur_cell,:)
            tmp_bxmax(child2,:)=tmp_bxmax(cur_cell,:)
            
            tmp_bxmax(child1,dir)=bisec_wall(cur_cell)
            tmp_bxmin(child2,dir)=bisec_wall(cur_cell)
            ! store index ranges
            lncpu = tmp_imax(cur_cell)-tmp_imin(cur_cell)+1
            tmp_imin(child1)=tmp_imin(cur_cell); tmp_imax(child1)=tmp_imin(cur_cell)+lncpu1(i)-1
            tmp_imin(child2)=tmp_imin(cur_cell)+lncpu1(i); tmp_imax(child2)=tmp_imax(cur_cell)
            ! store load of subcells
            tmp_load(child1)=load1(i)
            tmp_load(child2)=tmp_load(cur_cell)-load1(i)
            ! store walls for histogram updating
            walls(i)=bisec_wall(cur_cell)
         end do children_and_leaves


         if(lvl==nbilevelmax) exit    ! terminate level loop here if at bottom level

         ! split and sort for next histogram computation
         if(update) call splitsort_bisection_histogram(lvl,dir,walls)

         ! LEVEL CHANGE

         ! Advance level start cursor
         cur_levelstart = cur_levelstart + nc
         
         ! Alternate direction
         dir=dir+1; if(dir>ndim) dir=1

         ! Next level
         lvl = lvl + 1

      end do level_loop


      mean  = sum(dble(bisec_cpu_load))/ncpu
      var   = sum(dble(bisec_cpu_load)*dble(bisec_cpu_load))/ncpu
      stdev = sqrt(var-mean*mean)

      if(verbose .and. update) then
      print *,"Load balancing report ..."
      print *,"   Average CPU load     :", mean
      print *,"   Standard deviation   :", stdev
      print *,"   Balancing accuracy   :", stdev/mean
      print *,"   Requested tolerance  :", bisec_tol
      end if

      if(verbose) print *,'done with build_bisection'

   end subroutine build_bisection

   
   function round_to_bisec_res(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y

      y = dble(nint(x/bisec_res))*bisec_res
   end function round_to_bisec_res
    


   subroutine splitsort_bisection_histogram(lev,dir,walls)
      ! rearranges ind_to_part and updates ind_min and ind_max according to a split of domains in
      ! direction dir at wall positions walls
      ! this effectively creates two subdomains ready for histogram building
      ! Uses pre-computed bisec_cell_coord cache for coordinate lookups.

      real(dp), intent(in), dimension(:) :: walls
      integer, intent(in) :: dir, lev

      integer :: i, nc, lmost, rmost, tmp, tmp_level
      real(dp) :: tmp_coord
      integer(i8b) :: tmp_cost

      if(verbose) print *,'entering splitsort_bisection_histogram'

      new_hist_bounds=0

      nc=2**lev
      do i=1,nc
         lmost=bisec_hist_bounds(i); rmost=bisec_hist_bounds(i+1)-1       ! note the -1
         do while(rmost-lmost>0)
            if(bisec_cell_coord(lmost)<walls(i)) then    ! NOTE : strict inequality
               lmost=lmost+1
            else
               ! swap lmost and rmost: ind_cell, coord, cost
               tmp=bisec_ind_cell(lmost); bisec_ind_cell(lmost)=bisec_ind_cell(rmost); bisec_ind_cell(rmost)=tmp
               tmp_coord=bisec_cell_coord(lmost); bisec_cell_coord(lmost)=bisec_cell_coord(rmost); bisec_cell_coord(rmost)=tmp_coord
               tmp_cost=bisec_cell_cost(lmost); bisec_cell_cost(lmost)=bisec_cell_cost(rmost); bisec_cell_cost(rmost)=tmp_cost
               tmp_level=bisec_cell_level(lmost); bisec_cell_level(lmost)=bisec_cell_level(rmost); bisec_cell_level(rmost)=tmp_level
               rmost=rmost-1
            end if
         end do
         ! NOTE : rmost==lmost by now
         ! new histogram splitting bounds
         new_hist_bounds(2*i-1) = bisec_hist_bounds(i)
         new_hist_bounds(2*i  ) = lmost
         new_hist_bounds(2*i+1) = bisec_hist_bounds(i+1)
      end do
      bisec_hist_bounds=new_hist_bounds
   end subroutine splitsort_bisection_histogram


   subroutine init_bisection_histogram()
      ! This sets up bisec_ind_cell, bisec_hist_bounds ready to start a level-0 hist build
      ! by looping over all AMR cells

      use omp_lib, only: omp_get_max_threads, omp_get_thread_num

#ifndef WITHOUTMPI
      include 'mpif.h'
#endif

      integer::igrid,ncache,ngrid,ierr
      integer::ilevel,i,ind,idim

      integer::nc,ibcell,p

      integer::nx_loc
      integer::icpu,ncell,ncell_loc,ncell_max
      integer::nxny,ix,iy,iz,iskip
      integer::icell_tmp,igrid_tmp,isubcell_tmp
      integer::isink
      integer::npair_cell
      integer,allocatable::sink_per_grid(:),sink_coarse(:)
      integer(kind=8),dimension(1:MAXLEVEL)::niter_cost
      integer,dimension(1:twotondim)::npart_leaf,ndm_leaf
      real(dp),dimension(0:MAXLEVEL)::rank_scale
      integer::guard_iter
      integer::nthreads,tid,slot,cell_count
      integer,allocatable::thread_count(:,:),thread_offset(:,:)
      integer(kind=8)::cell_cost_tot_loc,cell_cost_tot
      integer(kind=8)::ntot_loc,ntot,nraised_loc,nraised
      integer(kind=8)::min_cost_loc,min_cost_global,floor_cost
      integer,parameter::LB_CELL_GUARD_MAXITER=3

      integer,dimension(1:nvector),save::ind_grid,ind_cell

      real(dp)::dx,scale,floor_c,cell_cap,guard_denom
      real(dp)::predicted_max_cells
      real(dp)::t_compact0,t_compact1,t_compact2,t_compact3
      integer(kind=8)::scratch_bytes,omp_scratch_bytes
      real(dp),dimension(1:twotondim,1:3)::xc
      logical::guard_applied
      
      if(verbose) print *,'entering init_bisection_histogram'

      ! Local constants
      nxny=nx*ny
      nx_loc=icoarse_max-icoarse_min+1
      scale=boxlen/dble(nx_loc)
      ncell=ncoarse+twotondim*ngridmax

      ! Use the same AMR subcycle factors as the Hilbert decomposition.
      ! domain_leaf_cost ignores them in memory-balance mode.
      niter_cost=1_8
      if((.not.memory_balance).and.cost_weighting)then
         niter_cost(levelmin)=1_8
         do ilevel=levelmin+1,nlevelmax
            if(niter_cost(ilevel-1)>huge(niter_cost(ilevel))/ &
                 int(max(1,nsubcycle(ilevel-1)),kind=8))then
               if(myid==1)write(*,*)'bisection: AMR subcycle cost overflow'
               stop
            end if
            niter_cost(ilevel)=int(max(1,nsubcycle(ilevel-1)),kind=8)* &
                 niter_cost(ilevel-1)
         end do
      end if

      rank_scale=1d0
      if((.not.memory_balance).and.time_balance_alpha>0d0)then
         do ilevel=levelmin,nlevelmax
            rank_scale(ilevel)=1d0+time_balance_alpha* &
                 (level_rank_scale_ema(ilevel)-1d0)
            rank_scale(ilevel)=max(0.5d0,min(2d0,rank_scale(ilevel)))
         end do
      end if

#ifndef WITHOUTMPI
      t_compact0=MPI_WTIME()
#endif

      ! --- Pass 1: exact per-thread leaf counts.  Reception grids cannot own
      ! myid cells by the build_comm invariant, so only active grids are needed.
      ! Keeping counts per (thread,level) lets Pass 2 fill compact arrays using
      ! prefix offsets, without an atomic increment in the hot loop.
      nthreads=max(1,omp_get_max_threads())
      allocate(thread_count(0:nthreads-1,0:nlevelmax))
      allocate(thread_offset(0:nthreads-1,0:nlevelmax))
      thread_count=0
      thread_offset=0

      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
      !$OMP PRIVATE(tid,cell_count,ind)
      tid=omp_get_thread_num()
      cell_count=0
      !$OMP DO SCHEDULE(STATIC)
      do ind=1,ncoarse
         if(cpu_map(ind)==myid .and. son(ind)==0) cell_count=cell_count+1
      end do
      !$OMP END DO
      thread_count(tid,0)=cell_count
      !$OMP END PARALLEL

      do ilevel=1,nlevelmax
         ncache=active(ilevel)%ngrid
         !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
         !$OMP PRIVATE(tid,cell_count,igrid,igrid_tmp,ind,iskip,icell_tmp)
         tid=omp_get_thread_num()
         cell_count=0
         !$OMP DO SCHEDULE(STATIC)
         do igrid=1,ncache
            igrid_tmp=active(ilevel)%igrid(igrid)
            do ind=1,twotondim
               iskip=ncoarse+(ind-1)*ngridmax
               icell_tmp=igrid_tmp+iskip
               if(cpu_map(icell_tmp)==myid .and. son(icell_tmp)==0) &
                    cell_count=cell_count+1
            end do
         end do
         !$OMP END DO
         thread_count(tid,ilevel)=cell_count
         !$OMP END PARALLEL
      end do

      ncell=0
      do ilevel=0,nlevelmax
         do tid=0,nthreads-1
            thread_offset(tid,ilevel)=ncell
            ncell=ncell+thread_count(tid,ilevel)
         end do
      end do
      ncell_max=ncell

      ! On-demand compact allocation to the exact leaf count.  Use a one-item
      ! allocation for a genuinely empty mesh so lower bounds remain valid.
      if(.not. allocated(bisec_ind_cell)) &
           allocate(bisec_ind_cell(1:max(ncell_max,1)))
      if(.not. allocated(bisec_cell_level)) &
           allocate(bisec_cell_level(1:max(ncell_max,1)))
      if(.not. allocated(bisec_cell_cost)) &
           allocate(bisec_cell_cost(1:max(ncell_max,1)))
      if(.not. allocated(bisec_cell_coord)) &
           allocate(bisec_cell_coord(1:max(ncell_max,1)))
      bisec_cell_cost=0_8
#ifndef WITHOUTMPI
      t_compact1=MPI_WTIME()
#endif

      ! --- Pass 2: deterministic thread-private compact fill. ---
      !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
      !$OMP PRIVATE(tid,slot,ind)
      tid=omp_get_thread_num()
      slot=thread_offset(tid,0)
      !$OMP DO SCHEDULE(STATIC)
      do ind=1,ncoarse
         if(cpu_map(ind)==myid.and.son(ind)==0)then
            slot=slot+1
            flag1(slot)=ind
            bisec_ind_cell(slot)=ind
            bisec_cell_level(slot)=0
            bisec_cell_cost(slot)=domain_leaf_cost(0,0,1_8, &
                 level_mesh_scale_ema(levelmin))
            bisec_cell_cost(slot)=max(1_8,nint( &
                 dble(bisec_cell_cost(slot))*rank_scale(levelmin),kind=8))
         end if
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      do ilevel=1,nlevelmax
         ncache=active(ilevel)%ngrid
         !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
         !$OMP PRIVATE(tid,slot,igrid,igrid_tmp,ind,iskip,icell_tmp, &
         !$OMP         npart_leaf,ndm_leaf,npair_cell)
         tid=omp_get_thread_num()
         slot=thread_offset(tid,ilevel)
         !$OMP DO SCHEDULE(STATIC)
         do igrid=1,ncache
            igrid_tmp=active(ilevel)%igrid(igrid)
            npart_leaf=0
            ndm_leaf=0
            if(pic) call count_particles_by_leaf(igrid_tmp,npart_leaf,ndm_leaf)
            do ind=1,twotondim
               iskip=ncoarse+(ind-1)*ngridmax
               icell_tmp=igrid_tmp+iskip
               if(cpu_map(icell_tmp)==myid.and.son(icell_tmp)==0)then
                  slot=slot+1
                  bisec_cell_level(slot)=ilevel
                  flag1(slot)=icell_tmp
                  bisec_ind_cell(slot)=icell_tmp
                  npair_cell=domain_sidm_pair_count(ndm_leaf(ind))
                  bisec_cell_cost(slot)=domain_leaf_cost( &
                       npart_leaf(ind),npair_cell,niter_cost(ilevel), &
                       level_mesh_scale_ema(ilevel))
                  bisec_cell_cost(slot)=max(1_8,nint( &
                       dble(bisec_cell_cost(slot))*rank_scale(ilevel),kind=8))
               end if
            end do
         end do
         !$OMP END DO
         !$OMP END PARALLEL
      end do
      ncell=ncell_max
      deallocate(thread_count,thread_offset)
#ifndef WITHOUTMPI
      t_compact2=MPI_WTIME()
#endif
      ! Ok, bisec_ind_cell is good, init the bound arrays for a level-0 histogram
      bisec_hist_bounds=0
      ! only one region (region id=1) at level 0
      bisec_hist_bounds(1)=1; bisec_hist_bounds(2)=ncell+1

      ! Store local cell count for cache arrays
      bisec_ncells_loc = ncell

      ! Add sink particle cost (computational weight near sinks)
      if (memory_balance .and. sink .and. nsink > 0 .and. mem_weight_sink > 0) then
         allocate(sink_per_grid(1:ngridmax))
         allocate(sink_coarse(1:ncoarse))
         sink_per_grid = 0
         sink_coarse = 0

         ! Map each sink to its leaf cell via AMR tree descent
         do isink = 1, nsink
            ix = int(xsink(isink,1) / scale)
            iy = int(xsink(isink,2) / scale)
            iz = int(xsink(isink,3) / scale)
            ix = min(max(ix, 0), nx-1)
            iy = min(max(iy, 0), ny-1)
            iz = min(max(iz, 0), nz-1)
            icell_tmp = 1 + ix + iy*nx + iz*nxny

            ! Descend AMR tree to leaf cell
            do while(son(icell_tmp) > 0)
               igrid_tmp = son(icell_tmp)
               ind = 1
               if(xsink(isink,1) >= xg(igrid_tmp,1)) ind = ind + 1
               if(xsink(isink,2) >= xg(igrid_tmp,2)) ind = ind + 2
               if(ndim > 2) then
                  if(xsink(isink,3) >= xg(igrid_tmp,3)) ind = ind + 4
               end if
               icell_tmp = igrid_tmp + ncoarse + (ind-1)*ngridmax
            end do

            ! Accumulate in grid or coarse cell
            if(icell_tmp > ncoarse) then
               isubcell_tmp = ((icell_tmp - ncoarse) / ngridmax) + 1
               igrid_tmp = icell_tmp - ncoarse - ngridmax * (isubcell_tmp - 1)
               sink_per_grid(igrid_tmp) = sink_per_grid(igrid_tmp) + 1
            else
               sink_coarse(icell_tmp) = sink_coarse(icell_tmp) + 1
            end if
         end do

         ! Add sink cost to bisec_cell_cost
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP PRIVATE(i,icell_tmp,isubcell_tmp,igrid_tmp) SCHEDULE(STATIC)
         do i = 1, bisec_ncells_loc
            icell_tmp = bisec_ind_cell(i)
            if(icell_tmp > ncoarse) then
               isubcell_tmp = ((icell_tmp - ncoarse) / ngridmax) + 1
               igrid_tmp = icell_tmp - ncoarse - ngridmax * (isubcell_tmp - 1)
               bisec_cell_cost(i) = bisec_cell_cost(i) + &
                    sink_per_grid(igrid_tmp) * mem_weight_sink / twotondim
            else
               bisec_cell_cost(i) = bisec_cell_cost(i) + &
                    sink_coarse(icell_tmp) * mem_weight_sink
            end if
         end do
         !$OMP END PARALLEL DO

         deallocate(sink_per_grid, sink_coarse)
      end if

      ! Shared bisection/ksection safety guard.  The cached entries are leaf
      ! cells: one AMR grid contributes twotondim cells and its grid cost is
      ! divided by twotondim above.  Thus the per-domain grid capacity becomes
      ! a cell capacity of lb_grid_headroom*ngridmax*twotondim.  With the
      ! current bounded memory-cost model this is normally a no-op; keeping the
      ! guard here protects both decompositions if their cost model is extended.
      guard_applied = .false.
      if(lb_grid_headroom > 0d0 .and. ngridmax > 0) then
         cell_cost_tot_loc = 0_8
         ntot_loc = 0_8
         min_cost_loc = huge(min_cost_loc)
         do i=1,bisec_ncells_loc
            if(bisec_cell_cost(i) > 0) then
               cell_cost_tot_loc = cell_cost_tot_loc+ &
                    int(bisec_cell_cost(i),kind=8)
               ntot_loc = ntot_loc+1_8
               min_cost_loc = min(min_cost_loc, &
                    int(bisec_cell_cost(i),kind=8))
            end if
         end do
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(cell_cost_tot_loc,cell_cost_tot,1, &
              MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(ntot_loc,ntot,1, &
              MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE(min_cost_loc,min_cost_global,1, &
              MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
#else
         cell_cost_tot = cell_cost_tot_loc
         ntot = ntot_loc
         min_cost_global = min_cost_loc
#endif

         ! Both split trees terminate in exactly ncpu domains.
         cell_cap = lb_grid_headroom*dble(ngridmax)*dble(twotondim)
         guard_denom = dble(ncpu)*cell_cap
         floor_c = dble(cell_cost_tot)/guard_denom
         if(ntot > 0_8 .and. cell_cost_tot > 0_8 .and. &
              floor_c > dble(min_cost_global)) then
            do guard_iter=1,LB_CELL_GUARD_MAXITER
               floor_c = dble(cell_cost_tot)/guard_denom
               if(floor_c <= dble(min_cost_global)) exit

               ! On the final pass choose floor >= T/(D-K).  Then raising at
               ! most K positive entries gives T_new <= D*floor.  D<=K means
               ! the requested occupancy cap cannot contain the current cells.
               if(guard_iter == LB_CELL_GUARD_MAXITER) then
                  if(guard_denom <= dble(ntot)) then
                     if(myid==1) write(*,*) &
                          ' LB cell guard (bisection/ksection): occupancy exceeds', &
                          ' the headroom; falling back to pure count balancing'
                     ! Headroom is unreachable, but an equal-count split is
                     ! still the best partition for the count limit and stays
                     ! feasible while occupancy is under the cell capacity.
                     do i=1,bisec_ncells_loc
                        if(bisec_cell_cost(i) > 0_i8b) bisec_cell_cost(i) = 1_i8b
                     end do
                     guard_applied = .true.
                     min_cost_global = 1_i8b
                     exit
                  end if
                  floor_c = max(floor_c,dble(cell_cost_tot)/ &
                       (guard_denom-dble(ntot)))
               end if

               if(floor_c > dble(huge(0_i8b))) then
                  if(myid==1) write(*,*) &
                       ' LB cell guard (bisection/ksection): integer cost overflow, floor=',floor_c
#ifndef WITHOUTMPI
                  call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
                  stop
               end if
               floor_cost = ceiling(floor_c,kind=8)

               nraised_loc = 0_8
               do i=1,bisec_ncells_loc
                  if(bisec_cell_cost(i) > 0 .and. &
                       bisec_cell_cost(i) < floor_cost) then
                     bisec_cell_cost(i) = floor_cost
                     nraised_loc = nraised_loc+1_8
                  end if
               end do
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(nraised_loc,nraised,1, &
                    MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               nraised = nraised_loc
#endif
               guard_applied = guard_applied .or. (nraised > 0_8)
               min_cost_global = max(min_cost_global,floor_cost)

               cell_cost_tot_loc = 0_8
               do i=1,bisec_ncells_loc
                  if(bisec_cell_cost(i) > 0) cell_cost_tot_loc = &
                       cell_cost_tot_loc+int(bisec_cell_cost(i),kind=8)
               end do
#ifndef WITHOUTMPI
               call MPI_ALLREDUCE(cell_cost_tot_loc,cell_cost_tot,1, &
                    MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
               cell_cost_tot = cell_cost_tot_loc
#endif
            end do

            predicted_max_cells = (dble(cell_cost_tot)/dble(ncpu))/ &
                 dble(min_cost_global)
            if(myid==1 .and. guard_applied) then
               write(*,'(A,ES12.4,A,I0,A,I0,A,F12.2,A,F12.2)') &
                    ' LB cell guard (bisection/ksection): floor=',floor_c, &
                    ' raised=',nraised,'/',ntot, &
                    ' predicted max cells=',predicted_max_cells, &
                    ' / cap=',cell_cap
            end if
         else if(myid==1 .and. verbose) then
            write(*,'(A)') ' LB cell guard (bisection/ksection): not needed'
         end if
      end if

#ifndef WITHOUTMPI
      t_compact3=MPI_WTIME()
      scratch_bytes=int(ncell_max,kind=8)*24_8
      omp_scratch_bytes=max( &
           int(2*nthreads*(nlevelmax+1),kind=8)*int(storage_size(0)/8,kind=8), &
           int(nthreads,kind=8)*int(bisec_nres,kind=8)*8_8)
      if(myid==1)then
         write(*,'(A,I0,A,F10.2,A,F8.3,A,I0,A)') &
              ' LB compact scratch: leaves/rank=',ncell_max,' compact=', &
              dble(scratch_bytes)/1048576d0,' MiB omp_tmp=', &
              dble(omp_scratch_bytes)/1048576d0,' MiB threads=',nthreads, &
              ' atomic=0'
         write(*,'(A,3(F9.3,A))') ' LB compact stages: count=', &
              t_compact1-t_compact0,' s fill=',t_compact2-t_compact1, &
              ' s cost/guard=',t_compact3-t_compact2,' s'
      end if
#endif

   end subroutine



   subroutine build_bisection_histogram(lev,dir,nc_in)
      ! build the nc_in histograms of level lev, for a split in direction dir
      ! nc_in = number of domains at this level (2**lev for bisection, product of factors for ksection)
      ! Uses pre-computed bisec_cell_coord and bisec_cell_cost cache arrays.

      use omp_lib, only: omp_get_max_threads, omp_get_thread_num

      integer, intent(in) :: lev, dir, nc_in
      integer :: nc
      integer :: i, ibicell,ithread,nthreads
      integer :: cell_slot
      integer(i8b),allocatable::thread_hist(:,:)

      nc=nc_in

      ! reset histograms before adding up loads
      bisec_hist=0

      nthreads=max(1,omp_get_max_threads())
      if(nc<nthreads .and. bisec_ncells_loc>0)then
         ! At the root (nc=1) a parallel-over-domain loop is serial.  Give each
         ! thread a small private 1-D histogram, then reduce it without atomics.
         ! Storage is O(nthreads*bisec_nres), independent of the leaf count.
         allocate(thread_hist(1:bisec_nres,0:nthreads-1))
         do ibicell=1,nc
            thread_hist=0_i8b
            !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
            !$OMP PRIVATE(ithread,i,cell_slot)
            ithread=omp_get_thread_num()
            !$OMP DO SCHEDULE(STATIC)
            do i=bisec_hist_bounds(ibicell),bisec_hist_bounds(ibicell+1)-1
               cell_slot=floor(bisec_cell_coord(i)/bisec_res)+1
               thread_hist(cell_slot,ithread)=thread_hist(cell_slot,ithread)+ &
                    bisec_cell_cost(i)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
            do ithread=0,nthreads-1
               do i=1,bisec_nres
                  bisec_hist(ibicell,i)=bisec_hist(ibicell,i)+ &
                       thread_hist(i,ithread)
               end do
            end do
            do i=2,bisec_nres
               bisec_hist(ibicell,i)=bisec_hist(ibicell,i)+bisec_hist(ibicell,i-1)
            end do
         end do
         deallocate(thread_hist)
      else
         ! Once there are at least as many domains as threads, rows are
         ! independent and can be built without atomics or private copies.
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP PRIVATE(ibicell,i,cell_slot) SCHEDULE(DYNAMIC,1)
         do ibicell=1,nc
            do i=bisec_hist_bounds(ibicell),bisec_hist_bounds(ibicell+1)-1
               cell_slot=floor(bisec_cell_coord(i)/bisec_res)+1
               bisec_hist(ibicell,cell_slot)=bisec_hist(ibicell,cell_slot)+ &
                    bisec_cell_cost(i)
            end do
            do i=2,bisec_nres
               bisec_hist(ibicell,i)=bisec_hist(ibicell,i)+bisec_hist(ibicell,i-1)
            end do
         end do
         !$OMP END PARALLEL DO
      end if

   end subroutine build_bisection_histogram

   subroutine compute_bisec_cell_coords(dir)
      ! Compute bisec_cell_coord(1:bisec_ncells_loc) for the given split direction.
      ! Must be called before build_bisection_histogram and splitsort at each level.
      integer, intent(in) :: dir

      integer :: i, ix, iy, iz, icell, igrid, isubcell, nx_loc, nxny
      integer, dimension(1:3) :: iarray, icoarse_array
      real(dp) :: dx, scale, subcell_c

      nxny = nx * ny
      icoarse_array = (/ icoarse_min, jcoarse_min, kcoarse_min /)
      nx_loc = icoarse_max - icoarse_min + 1
      scale = boxlen / dble(nx_loc)

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP PRIVATE(i,ix,iy,iz,icell,igrid,isubcell,iarray,dx,subcell_c) SCHEDULE(STATIC)
      do i = 1, bisec_ncells_loc
         icell = bisec_ind_cell(i)
         dx = 0.5d0**bisec_cell_level(i)

         if (icell <= ncoarse) then
            iz = (icell - 1) / nxny
            iy = ((icell - 1) - iz * nxny) / nx
            ix = ((icell - 1) - iy * nx - iz * nxny)
            iarray = (/ ix, iy, iz /)
            bisec_cell_coord(i) = scale * (dble(iarray(dir)) - 0.5d0) * dx
         else
            isubcell = ((icell - ncoarse) / ngridmax) + 1
            igrid = icell - ncoarse - ngridmax * (isubcell - 1)
            iz = (isubcell - 1) / 4
            iy = (isubcell - 1 - 4 * iz) / 2
            ix = (isubcell - 1 - 2 * iy - 4 * iz)
            iarray = (/ ix, iy, iz /)
            subcell_c = (dble(iarray(dir)) - 0.5d0) * dx - dble(icoarse_array(dir))
            bisec_cell_coord(i) = scale * (xg(igrid, dir) + subcell_c)
         end if
      end do
      !$OMP END PARALLEL DO

   end subroutine compute_bisec_cell_coords

end module bisection

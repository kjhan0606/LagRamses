subroutine init_hydro
  use amr_commons
  use hydro_commons
#ifdef RT      
  use rt_parameters,only: convert_birth_times
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar,irad
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering init_hydro'
  
  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar))
  allocate(unew(1:ncell,1:nvar))
  ! uold/unew: Active cells initialized by restart reader or init_flow_fine.
  ! Free-list cells get zero from mmap(MAP_ANONYMOUS) lazy page allocation.
  ! Skip full-array zeroing to avoid paging in 18 GB at startup.
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
#ifdef HDF5
     if(informat == 'hdf5') then
        call restore_hydro_hdf5()
        if(verbose)write(*,*)'HDF5 HYDRO backup files read completed'
        call sgs_init_restart
        return
     end if
#endif
     if(varcpu_restart) then
        call restore_hydro_binary_varcpu()
#ifndef WITHOUTMPI
        call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        call diag_check_nan('post_hydro_restore')
        if(verbose)write(*,*)'Binary varcpu HYDRO backup files read completed'
        call sgs_init_restart
        return
     end if
     ilun=ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     endif



     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token                                                                                                                                                                    
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(.not.(neq_chem.or.rt) .and. nvar2.ne.nvar)then
        write(*,*)'File hydro.tmp is not compatible'
        write(*,*)'Found   =',nvar2
        write(*,*)'Expected=',nvar
        call clean_stop
     end if
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading first ',nvar2, &
                  'variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar)then ! Not OK to drop variables 
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        call clean_stop
     end if
#endif
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax

                 ! Read density and velocities --> density and momenta
                 do ivar=1,ndim+1
                    read(ilun)xx
                    if(ivar==1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,1)=xx(i)
                       end do
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 end do

#if NENER>0
                 ! Read non-thermal pressures --> non-thermal energies
                 do ivar=ndim+3,ndim+2+nener
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)/(gamma_rad(ivar-ndim-2)-1d0)
                    end do
                 end do
#endif
                 ! Read thermal pressure --> total fluid energy
                 read(ilun)xx
                 do i=1,ncache
                    xx(i)=xx(i)/(gamma-1d0)
                    if (uold(ind_grid(i)+iskip,1)>0.)then
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,2)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#if NDIM>1
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,3)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NDIM>2
                    xx(i)=xx(i)+0.5d0*uold(ind_grid(i)+iskip,4)**2/max(uold(ind_grid(i)+iskip,1),smallr)
#endif
#if NENER>0
                    do irad=1,nener
                       xx(i)=xx(i)+uold(ind_grid(i)+iskip,ndim+2+irad)
                    end do
#endif
                 else
                    xx(i)=0.
                 end if
                    uold(ind_grid(i)+iskip,ndim+2)=xx(i)
                 end do
#if NVAR>NDIM+2+NENER
                 ! Read passive scalars
                 do ivar=ndim+3+nener,min(nvar,nvar2)
                    read(ilun)xx
                    do i=1,ncache
                       uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                 end do
#endif
              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)

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



#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     call diag_check_nan('post_hydro_restore_std')
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro

subroutine restore_hydro_binary_varcpu
  implicit none
  call restore_hydro_binary_varcpu_streaming
end subroutine restore_hydro_binary_varcpu
!################################################################
! Kept temporarily as a reference for binary-format compatibility.  The
! production entry point above uses the bounded field-streaming implementation
! below instead of constructing one 91-double record per grid.
subroutine restore_hydro_binary_varcpu_legacy
  !--------------------------------------------------------------
  ! Chunked distributed I/O version: reads hydro files in chunks,
  ! exchanges via ksection_exchange_dp (O(log_k ncpu) memory),
  ! and uses Morton hash lookup for position→igrid mapping.
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use ksection
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, k, ind, ivar, irad, info, ifile
  integer :: ncache, ilun, nvar2, nvar_send, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, nprops, base
  integer :: nlocal, nrecv
  real(dp) :: gamma2, eval, rho_val, twotol
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  ! Chunking variables
  integer :: nchunk, ichunk, chunk_file_lo, chunk_file_hi
  integer :: chunk_ngrids(1:MAXLEVEL), chunk_lvl_offset(1:MAXLEVEL)
  integer :: n_from_file, xg_base

  ! Per-level chunk hydro buffers
  type hydro_level_t
     real(dp), allocatable :: udata(:,:)  ! (ngrids, nval_per_grid)
  end type
  type(hydro_level_t) :: chunk_hvl(1:MAXLEVEL)

  ! Ksection exchange buffers
  real(dp), allocatable :: sendbuf_2d(:,:), recvbuf_2d(:,:)
  integer, allocatable :: dest(:)

  ! Morton hash lookup
  integer(8) :: ixm, iym, izm
  type(mkey_t) :: mkey
  real(dp) :: xg_recv(3), xx_pos(1:nvector, 1:ndim), scale
  integer :: c_tmp(1:nvector), nx_loc, nxny

  if(myid==1) write(*,*) 'Binary varcpu hydro restore (chunked ksection): ncpu_file=', ncpu_file

  ilun = 99
  call title(nrestart, nchar)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  nxny = nx * ny

  ! Get nvar2 from file 00001
  if(myid == 1) then
     if(IOGROUPSIZEREP>0) then
        call title(1, ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     end if
     call title(1, ncharcpu)
     fileloc=TRIM(fileloc)//TRIM(ncharcpu)
     open(unit=ilun, file=fileloc, form='unformatted')
     read(ilun)  ! ncpu
     read(ilun) nvar2
     close(ilun)
  end if
  call MPI_BCAST(nvar2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)

  nvar_send = min(nvar, nvar2)
  nval_per_grid = twotondim * nvar_send
  nprops = ndim + nval_per_grid  ! xg(1:3) prepended to data

  ! Compute chunk boundaries on global file indices
  if(varcpu_chunk_nfile <= 0) then
     nchunk = 1
  else
     nchunk = (ncpu_file + varcpu_chunk_nfile - 1) / varcpu_chunk_nfile
  end if

  ! Main chunked loop
  do ichunk = 1, nchunk
     if(varcpu_chunk_nfile <= 0) then
        chunk_file_lo = 1
        chunk_file_hi = ncpu_file
     else
        chunk_file_lo = (ichunk - 1) * varcpu_chunk_nfile + 1
        chunk_file_hi = min(ichunk * varcpu_chunk_nfile, ncpu_file)
     end if

     ! Compute per-level grid counts for this chunk (my files only)
     chunk_ngrids = 0
     do j = 1, varcpu_nfiles_local
        if(varcpu_my_files(j) < chunk_file_lo .or. &
           varcpu_my_files(j) > chunk_file_hi) cycle
        do ilevel = 1, nlevelmax
           chunk_ngrids(ilevel) = chunk_ngrids(ilevel) + &
                varcpu_nactive(varcpu_my_files(j), ilevel)
        end do
     end do

     ! Allocate chunk-level hydro buffers
     do ilevel = 1, nlevelmax
        if(chunk_ngrids(ilevel) > 0) &
             allocate(chunk_hvl(ilevel)%udata(chunk_ngrids(ilevel), nval_per_grid))
     end do

     ! Read assigned hydro files in this chunk
     chunk_lvl_offset = 0
     do ifile = 1, varcpu_nfiles_local
        icpu_file = varcpu_my_files(ifile)
        if(icpu_file < chunk_file_lo .or. icpu_file > chunk_file_hi) cycle

        if(IOGROUPSIZEREP>0) then
           call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
           fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
        else
           fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
        end if
        call title(icpu_file, ncharcpu)
        fileloc=TRIM(fileloc)//TRIM(ncharcpu)

        open(unit=ilun, file=fileloc, form='unformatted')
        read(ilun) ncpu2
        read(ilun) nvar2
        read(ilun) ndim2
        read(ilun) nlevelmax2
        read(ilun) nboundary2
        read(ilun) gamma2

        do ilevel = 1, nlevelmax2
           do ibound = 1, nboundary2 + ncpu2
              read(ilun)  ! ilevel2
              read(ilun) ncache

              if(ncache > 0) then
                 if(ibound == icpu_file) then
                    allocate(xx(1:ncache))
                    ind = chunk_lvl_offset(ilevel)
                    do iskip = 1, twotondim
                       do ivar = 1, nvar_send
                          read(ilun) xx
                          do i = 1, ncache
                             chunk_hvl(ilevel)%udata(ind+i, (iskip-1)*nvar_send+ivar) = xx(i)
                          end do
                       end do
                       do ivar = nvar_send+1, nvar2
                          read(ilun)  ! skip extra variables
                       end do
                    end do
                    chunk_lvl_offset(ilevel) = ind + ncache
                    deallocate(xx)
                 else
                    do i = 1, twotondim * nvar2
                       read(ilun)
                    end do
                 end if
              end if
           end do
        end do
        close(ilun)
     end do

     ! Exchange and scatter level by level
     do ilevel = 1, nlevelmax
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        nlocal = chunk_ngrids(ilevel)
        twotol = 2.0d0**(ilevel-1)

        ! Pack sendbuf: xg(1:ndim) + udata(1:nval_per_grid)
        allocate(sendbuf_2d(1:nprops, 1:max(nlocal,1)))
        allocate(dest(max(nlocal,1)))
        k = 0
        do j = 1, varcpu_nfiles_local
           if(varcpu_my_files(j) < chunk_file_lo .or. &
              varcpu_my_files(j) > chunk_file_hi) cycle
           n_from_file = varcpu_nactive(varcpu_my_files(j), ilevel)
           xg_base = varcpu_file_start(j-1, ilevel)
           do i = 1, n_from_file
              k = k + 1
              ! Grid position from varcpu_lvl (saved during AMR restore)
              sendbuf_2d(1, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 1)
              sendbuf_2d(2, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 2)
              sendbuf_2d(3, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 3)
              ! Hydro data
              sendbuf_2d(ndim+1:nprops, k) = chunk_hvl(ilevel)%udata(k, 1:nval_per_grid)
              ! Determine owner CPU
              xx_pos(1,1) = (sendbuf_2d(1, k) - dble(icoarse_min)) * scale
              xx_pos(1,2) = (sendbuf_2d(2, k) - dble(jcoarse_min)) * scale
              xx_pos(1,3) = (sendbuf_2d(3, k) - dble(kcoarse_min)) * scale
              if(ordering == 'ksection') then
                 call cmp_ksection_cpumap(xx_pos, c_tmp, 1)
              else
                 call cmp_cpumap(xx_pos, c_tmp, 1)
              end if
              dest(k) = c_tmp(1)
           end do
        end do

        ! Ksection hierarchical exchange
        call ksection_exchange_dp(sendbuf_2d, nlocal, dest, nprops, recvbuf_2d, nrecv)
        deallocate(sendbuf_2d, dest)

        ! Scatter to local grids with primitive → conservative conversion
        do i = 1, nrecv
           xg_recv(1:ndim) = recvbuf_2d(1:ndim, i)

           ! Morton hash lookup → igrid
           ixm = int(xg_recv(1) * twotol, 8)
           iym = int(xg_recv(2) * twotol, 8)
           izm = int(xg_recv(3) * twotol, 8)
           mkey = morton_encode(ixm, iym, izm)
           igrid = morton_hash_lookup(mort_table(ilevel), mkey)
           if(igrid == 0) cycle

           do iskip = 1, twotondim
              icell = igrid + ncoarse + (iskip-1)*ngridmax
              base = ndim + (iskip-1)*nvar_send

              ! Density (file record 1)
              rho_val = recvbuf_2d(base + 1, i)
              uold(icell, 1) = rho_val

              ! Velocities → momenta (file records 2..ndim+1)
              do ivar = 2, ndim+1
                 uold(icell, ivar) = recvbuf_2d(base + ivar, i) * max(rho_val, smallr)
              end do

#if NENER>0
              ! Non-thermal pressures → energies
              do irad = 1, nener
                 uold(icell, ndim+2+irad) = recvbuf_2d(base + ndim+1+irad, i) / &
                      (gamma_rad(irad) - 1d0)
              end do
#endif

              ! Thermal pressure → total energy
              eval = recvbuf_2d(base + ndim+2+nener, i) / (gamma - 1d0)
              if(rho_val > 0d0) then
                 do ivar = 2, ndim+1
                    eval = eval + 0.5d0*uold(icell,ivar)**2 / max(rho_val, smallr)
                 end do
#if NENER>0
                 do irad = 1, nener
                    eval = eval + uold(icell, ndim+2+irad)
                 end do
#endif
              else
                 eval = 0d0
              end if
              uold(icell, ndim+2) = eval

#if NVAR>NDIM+2+NENER
              ! Passive scalars
              do ivar = ndim+3+nener, nvar_send
                 uold(icell, ivar) = recvbuf_2d(base + ivar, i) * max(rho_val, smallr)
              end do
#endif
           end do
        end do
        deallocate(recvbuf_2d)

        ! Free chunk level buffer
        if(allocated(chunk_hvl(ilevel)%udata)) deallocate(chunk_hvl(ilevel)%udata)
     end do
  end do  ! ichunk

  if(myid==1) write(*,*) 'Binary varcpu hydro restore done.'

end subroutine restore_hydro_binary_varcpu_legacy
!################################################################
subroutine restore_hydro_binary_varcpu_streaming
  ! Stream one scalar field at a time from each old-CPU file.  Direct P2P uses
  ! known peer counts, so k-section ranks never aggregate 91-double records.
  use amr_commons
  use hydro_commons
  use dynamic_exchange, only: EXCHANGE_SPARSE_P2P, exchange_dp_sorted
  use morton_keys
  use morton_hash
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: ilun,info,nvar2,ncpu2,ndim2,nlevelmax2,nboundary2
  integer :: icpu_file,nslot,islot,ilevel,ilevel_file,ibound,ncache
  integer :: i,j,k,idim,ivar,iskip,irad,nvec,nlocal,nrecv
  integer :: nvar_send,xg_base,icpu,icell,igrid,tag
  integer :: nlocal_global,nrecv_global,nmissing_local,nmissing_global
  integer :: ns(1:ncpu),nr(1:ncpu),sd(1:ncpu),rd(1:ncpu),offsets(1:ncpu)
  integer,allocatable :: target_cpu(:),sort_index(:),recv_grid(:)
  real(dp),allocatable :: file_value(:),send_value(:),recv_value(:)
  real(dp),allocatable :: recv_xg(:,:)
  real(dp) :: gamma2,rho_val,eval,twotol,scale
  real(dp) :: xpos(1:nvector,1:ndim)
  integer :: ctmp(1:nvector)
  integer(8) :: ixm,iym,izm
  type(mkey_t) :: mkey
  character(len=80) :: fileloc
  character(len=5) :: nchar,ncharcpu
  logical :: have_file

  if(myid==1) write(*,*) &
       'Binary varcpu hydro restore (field-streamed large-count P2P): ncpu_file=',ncpu_file
  ilun=99
  call title(nrestart,nchar)
  scale=boxlen/dble(icoarse_max-icoarse_min+1)

  ! Read and broadcast the common header.  Reader ranks reopen their assigned
  ! file below and advance through it exactly once.
  if(myid==1)then
     if(IOGROUPSIZEREP>0)then
        call title(1,ncharcpu)
        fileloc='output_'//trim(nchar)//'/group_'//trim(ncharcpu)// &
             '/hydro_'//trim(nchar)//'.out'
     else
        fileloc='output_'//trim(nchar)//'/hydro_'//trim(nchar)//'.out'
     endif
     call title(1,ncharcpu)
     fileloc=trim(fileloc)//trim(ncharcpu)
     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2; read(ilun)nvar2; read(ilun)ndim2
     read(ilun)nlevelmax2; read(ilun)nboundary2; read(ilun)gamma2
     close(ilun)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(ncpu2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(nvar2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(ndim2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(nlevelmax2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(nboundary2,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST(gamma2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(varcpu_nfiles_local,nslot,1,MPI_INTEGER,MPI_MAX, &
       MPI_COMM_WORLD,info)
#else
  nslot=varcpu_nfiles_local
#endif
  nvar_send=min(nvar,nvar2)

  ! One slot contains at most one old file per reader.  For 32 -> 128 all old
  ! files therefore stream in parallel in slot one.
  do islot=1,nslot
     have_file=islot<=varcpu_nfiles_local
     icpu_file=0
     if(have_file)then
        icpu_file=varcpu_my_files(islot)
        if(IOGROUPSIZEREP>0)then
           call title(((icpu_file-1)/IOGROUPSIZEREP)+1,ncharcpu)
           fileloc='output_'//trim(nchar)//'/group_'//trim(ncharcpu)// &
                '/hydro_'//trim(nchar)//'.out'
        else
           fileloc='output_'//trim(nchar)//'/hydro_'//trim(nchar)//'.out'
        endif
        call title(icpu_file,ncharcpu)
        fileloc=trim(fileloc)//trim(ncharcpu)
        open(unit=ilun,file=fileloc,form='unformatted')
        read(ilun)ncpu2; read(ilun)nvar2; read(ilun)ndim2
        read(ilun)nlevelmax2; read(ilun)nboundary2; read(ilun)gamma2
     endif

     do ilevel=1,nlevelmax2
        nlocal=0
        ! Reach this file's active-domain record.  Ranks arrive at the first
        ! coordinate exchange together even though their old-domain index is
        ! different.
        if(have_file)then
           do ibound=1,icpu_file-1
              read(ilun)ilevel_file; read(ilun)ncache
              call skip_hydro_records(ilun,ncache,nvar2)
           enddo
           read(ilun)ilevel_file; read(ilun)nlocal
           if(ilevel_file/=ilevel)then
              write(*,*)'FATAL streamed hydro level mismatch',myid,ilevel_file,ilevel
              call clean_stop
           endif
           if(nlocal/=varcpu_nactive(icpu_file,ilevel))then
              write(*,*)'FATAL streamed hydro grid-count mismatch',myid, &
                   icpu_file,ilevel,nlocal,varcpu_nactive(icpu_file,ilevel)
              call clean_stop
           endif
        endif

        allocate(target_cpu(max(nlocal,1)),sort_index(max(nlocal,1)))
        xg_base=0
        if(have_file)xg_base=varcpu_file_start(islot-1,ilevel)
        do i=1,nlocal,nvector
           nvec=min(nvector,nlocal-i+1)
           do k=1,nvec
              xpos(k,1)=(varcpu_lvl(ilevel)%xg(xg_base+i+k-1,1)- &
                   dble(icoarse_min))*scale
#if NDIM>1
              xpos(k,2)=(varcpu_lvl(ilevel)%xg(xg_base+i+k-1,2)- &
                   dble(jcoarse_min))*scale
#endif
#if NDIM>2
              xpos(k,3)=(varcpu_lvl(ilevel)%xg(xg_base+i+k-1,3)- &
                   dble(kcoarse_min))*scale
#endif
           enddo
           call cmp_cpumap(xpos,ctmp,nvec)
           target_cpu(i:i+nvec-1)=ctmp(1:nvec)
        enddo

        ns=0
        do i=1,nlocal
           ns(target_cpu(i))=ns(target_cpu(i))+1
        enddo
#ifndef WITHOUTMPI
        call MPI_ALLTOALL(ns,1,MPI_INTEGER,nr,1,MPI_INTEGER, &
             MPI_COMM_WORLD,info)
#else
        nr=ns
#endif
        sd(1)=0; rd(1)=0
        do icpu=2,ncpu
           sd(icpu)=sd(icpu-1)+ns(icpu-1)
           rd(icpu)=rd(icpu-1)+nr(icpu-1)
        enddo
        nrecv=sum(nr)
        offsets=sd
        do i=1,nlocal
           icpu=target_cpu(i)
           offsets(icpu)=offsets(icpu)+1
           sort_index(offsets(icpu))=i
        enddo

        allocate(send_value(max(nlocal,1)),recv_value(max(nrecv,1)))
        allocate(recv_xg(max(nrecv,1),ndim),recv_grid(max(nrecv,1)))
        do idim=1,ndim
           do i=1,nlocal
              send_value(i)=varcpu_lvl(ilevel)%xg(xg_base+sort_index(i),idim)
           enddo
           call exchange_dp_sorted(send_value,recv_value,ns,nr,sd,rd, &
                EXCHANGE_SPARSE_P2P,900+idim)
           if(nrecv>0)recv_xg(1:nrecv,idim)=recv_value(1:nrecv)
        enddo

        twotol=2.0d0**(ilevel-1)
        do i=1,nrecv
           ixm=int(recv_xg(i,1)*twotol,8)
           iym=int(recv_xg(i,2)*twotol,8)
           izm=int(recv_xg(i,3)*twotol,8)
           mkey=morton_encode(ixm,iym,izm)
           recv_grid(i)=morton_hash_lookup(mort_table(ilevel),mkey)
        enddo
        nmissing_local=count(recv_grid(1:nrecv)==0)
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nmissing_local,nmissing_global,1,MPI_INTEGER, &
             MPI_SUM,MPI_COMM_WORLD,info)
#else
        nmissing_global=nmissing_local
#endif
        if(nmissing_global/=0)then
           if(myid==1)write(*,*)'FATAL streamed hydro Morton misses=', &
                nmissing_global,' slot=',islot,' level=',ilevel
           call clean_stop
        endif
        deallocate(recv_xg,target_cpu)

        ! The file is primitive.  Density and momenta arrive before pressure,
        ! so total energy can be assembled immediately without retaining a
        ! multi-field record.
        allocate(file_value(max(nlocal,1)))
        do iskip=1,twotondim
           do ivar=1,nvar2
              if(have_file.and.nlocal>0)read(ilun)file_value(1:nlocal)
              if(ivar<=nvar_send)then
                 do i=1,nlocal
                    send_value(i)=file_value(sort_index(i))
                 enddo
                 tag=920+(iskip-1)*nvar2+ivar
                 call exchange_dp_sorted(send_value,recv_value,ns,nr,sd,rd, &
                      EXCHANGE_SPARSE_P2P,tag)
                 do i=1,nrecv
                    igrid=recv_grid(i)
                    if(igrid==0)cycle
                    icell=igrid+ncoarse+(iskip-1)*ngridmax
                    if(ivar==1)then
                       uold(icell,1)=recv_value(i)
                    else if(ivar>=2.and.ivar<=ndim+1)then
                       uold(icell,ivar)=recv_value(i)*max(uold(icell,1),smallr)
#if NENER>0
                    else if(ivar>=ndim+2.and.ivar<=ndim+1+nener)then
                       irad=ivar-(ndim+1)
                       uold(icell,ndim+2+irad)=recv_value(i)/ &
                            (gamma_rad(irad)-1d0)
#endif
                    else if(ivar==ndim+2+nener)then
                       rho_val=uold(icell,1)
                       eval=recv_value(i)/(gamma-1d0)
                       if(rho_val>0d0)then
                          do j=2,ndim+1
                             eval=eval+0.5d0*uold(icell,j)**2/max(rho_val,smallr)
                          enddo
#if NENER>0
                          do irad=1,nener
                             eval=eval+uold(icell,ndim+2+irad)
                          enddo
#endif
                       else
                          eval=0d0
                       endif
                       uold(icell,ndim+2)=eval
#if NVAR>NDIM+2+NENER
                    else if(ivar>=ndim+3+nener)then
                       uold(icell,ivar)=recv_value(i)*max(uold(icell,1),smallr)
#endif
                    endif
                 enddo
              endif
              if(myid==1.and.(mod(ivar,8)==0.or.ivar==nvar2)) &
                   write(*,'(A,I0,A,I0,A,I0,A,I0)') &
                   ' Hydro stream slot ',islot,' level ',ilevel,' field ', &
                   (iskip-1)*nvar2+ivar,'/',twotondim*nvar2
           enddo
        enddo
        deallocate(file_value,send_value,recv_value,recv_grid,sort_index)

        if(have_file)then
           do ibound=icpu_file+1,nboundary2+ncpu2
              read(ilun)ilevel_file; read(ilun)ncache
              call skip_hydro_records(ilun,ncache,nvar2)
           enddo
        endif
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nlocal,nlocal_global,1,MPI_INTEGER,MPI_SUM, &
             MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(nrecv,nrecv_global,1,MPI_INTEGER,MPI_SUM, &
             MPI_COMM_WORLD,info)
#else
        nlocal_global=nlocal; nrecv_global=nrecv
#endif
        if(myid==1)write(*,'(A,I0,A,I0,A,I0,A,I0)') &
             ' Hydro stream completed slot ',islot,' level ',ilevel, &
             ' sent grids=',nlocal_global,' received grids=',nrecv_global
     enddo
     if(have_file)close(ilun)
  enddo
  if(myid==1)write(*,*)'Binary varcpu hydro restore done.'

contains
  subroutine skip_hydro_records(iunit,ngrid_file,nvar_file)
    integer,intent(in)::iunit,ngrid_file,nvar_file
    integer::irec
    if(ngrid_file<=0)return
    do irec=1,twotondim*nvar_file
       read(iunit)
    enddo
  end subroutine skip_hydro_records
end subroutine restore_hydro_binary_varcpu_streaming

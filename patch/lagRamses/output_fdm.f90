subroutine backup_psi(filename)
  use amr_commons
  use poisson_commons, only: psi_re, psi_im
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=80)::filename

  integer::i,ncache,ind,ilevel,igrid,ilun,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1131
  integer::dummy_io,info2

  if(.not.use_fdm)return
  if(verbose)write(*,*)'Entering backup_psi'

  ilun=ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZEOUT>0) then
     if (mod(myid-1,IOGROUPSIZEOUT)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              ! Write Re(psi)
              do i=1,ncache
                 xdp(i)=psi_re(ICELL_OF(ind_grid(i),ind))
              end do
              write(ilun)xdp
              ! Write Im(psi)
              do i=1,ncache
                 xdp(i)=psi_im(ICELL_OF(ind_grid(i),ind))
              end do
              write(ilun)xdp
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do
  close(ilun)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZEOUT>0) then
     if(mod(myid,IOGROUPSIZEOUT)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

end subroutine backup_psi
!#########################################################################
!#########################################################################
!#########################################################################
subroutine restore_psi
  !--------------------------------------------------------------
  ! Binary-mode restart reader for the FDM wavefunction psi.
  ! Mirrors the standard binary read in init_hydro/init_poisson.
  ! psi_re/psi_im are already allocated by init_poisson under use_fdm.
  ! Reads the per-CPU fdm_<nrestart>.out file written by backup_psi.
  !--------------------------------------------------------------
  use amr_commons
  use poisson_commons, only: psi_re, psi_im
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncache,igrid,i,ilevel,ind,ilun,ibound,istart,info
  integer::ncpu2,ndim2,nlevelmax2,nboundary2,ilevel2,numbl2
  integer,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1132
  integer::dummy_io,info2

  if(.not.use_fdm)return
  if(nrestart<=0)return
  if(verbose)write(*,*)'Entering restore_psi'

  ilun=ncpu+myid+10
  call title(nrestart,nchar)
  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/fdm_'//TRIM(nchar)//'.out'
  else
     fileloc='output_'//TRIM(nchar)//'/fdm_'//TRIM(nchar)//'.out'
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
  read(ilun)ndim2
  read(ilun)nlevelmax2
  read(ilun)nboundary2
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
           write(*,*)'File fdm.tmp is not compatible'
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
              ! Read Re(psi)
              read(ilun)xx
              do i=1,ncache
                 psi_re(ICELL_OF(ind_grid(i),ind))=xx(i)
              end do
              ! Read Im(psi)
              read(ilun)xx
              do i=1,ncache
                 psi_im(ICELL_OF(ind_grid(i),ind))=xx(i)
              end do
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
  if(debug)write(*,*)'fdm.tmp read for processor ',myid
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  if(verbose)write(*,*)'FDM (psi) backup files read completed'

end subroutine restore_psi
!#########################################################################
!#########################################################################
!#########################################################################
subroutine restore_psi_binary_varcpu
  use amr_commons
  use poisson_commons, only: psi_re, psi_im
  use ksection
  use morton_keys
  use morton_hash
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_file, ilevel, ibound, i, j, k, ind, info, ifile
  integer :: ncache, ilun, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, nprops, base
  integer :: nlocal, nrecv
  real(dp), allocatable :: xx(:)
  character(LEN=80) :: fileloc
  character(LEN=5) :: nchar, ncharcpu

  integer :: nchunk, ichunk, chunk_file_lo, chunk_file_hi
  integer :: chunk_ngrids(1:MAXLEVEL), chunk_lvl_offset(1:MAXLEVEL)
  integer :: n_from_file, xg_base

  type psi_level_t
     real(dp), allocatable :: pdata(:,:)
  end type
  type(psi_level_t) :: chunk_pvl(1:MAXLEVEL)

  real(dp), allocatable :: sendbuf_2d(:,:), recvbuf_2d(:,:)
  integer, allocatable :: dest(:)

  integer(8) :: ixm, iym, izm
  type(mkey_t) :: mkey
  real(dp) :: xg_recv(3), xx_pos(1:nvector, 1:ndim), scale, twotol
  integer :: c_tmp(1:nvector), nx_loc

  integer :: nread_total, nrecv_total, nskip_total, info2
  integer :: nread_g, nrecv_g, nskip_g
  real(dp) :: psi2_read, psi2_recv, psi2_read_g, psi2_recv_g

  if(.not.use_fdm) return
  if(myid==1) write(*,*) 'Binary varcpu FDM psi restore (chunked ksection): ncpu_file=', ncpu_file
  if(myid==1) write(*,*) '  varcpu_nfiles_local=', varcpu_nfiles_local

  nread_total = 0; nrecv_total = 0; nskip_total = 0
  psi2_read = 0.0d0; psi2_recv = 0.0d0

  ilun = 99
  call title(nrestart, nchar)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  nval_per_grid = twotondim * 2
  nprops = ndim + nval_per_grid

  if(varcpu_chunk_nfile <= 0) then
     nchunk = 1
  else
     nchunk = (ncpu_file + varcpu_chunk_nfile - 1) / varcpu_chunk_nfile
  end if

  do ichunk = 1, nchunk
     if(varcpu_chunk_nfile <= 0) then
        chunk_file_lo = 1
        chunk_file_hi = ncpu_file
     else
        chunk_file_lo = (ichunk - 1) * varcpu_chunk_nfile + 1
        chunk_file_hi = min(ichunk * varcpu_chunk_nfile, ncpu_file)
     end if

     chunk_ngrids = 0
     do j = 1, varcpu_nfiles_local
        if(varcpu_my_files(j) < chunk_file_lo .or. &
           varcpu_my_files(j) > chunk_file_hi) cycle
        do ilevel = 1, nlevelmax
           chunk_ngrids(ilevel) = chunk_ngrids(ilevel) + &
                varcpu_nactive(varcpu_my_files(j), ilevel)
        end do
     end do

     do ilevel = 1, nlevelmax
        if(chunk_ngrids(ilevel) > 0) &
             allocate(chunk_pvl(ilevel)%pdata(chunk_ngrids(ilevel), nval_per_grid))
     end do

     chunk_lvl_offset = 0
     do ifile = 1, varcpu_nfiles_local
        icpu_file = varcpu_my_files(ifile)
        if(icpu_file < chunk_file_lo .or. icpu_file > chunk_file_hi) cycle

        if(IOGROUPSIZEREP>0) then
           call title(((icpu_file-1)/IOGROUPSIZEREP)+1, ncharcpu)
           fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/fdm_'//TRIM(nchar)//'.out'
        else
           fileloc='output_'//TRIM(nchar)//'/fdm_'//TRIM(nchar)//'.out'
        end if
        call title(icpu_file, ncharcpu)
        fileloc=TRIM(fileloc)//TRIM(ncharcpu)

        open(unit=ilun, file=fileloc, form='unformatted')
        read(ilun) ncpu2
        read(ilun) ndim2
        read(ilun) nlevelmax2
        read(ilun) nboundary2

        do ilevel = 1, nlevelmax2
           do ibound = 1, nboundary2 + ncpu2
              read(ilun)  ! ilevel2
              read(ilun) ncache

              if(ncache > 0) then
                 if(ibound == icpu_file) then
                    allocate(xx(1:ncache))
                    ind = chunk_lvl_offset(ilevel)
                    do iskip = 1, twotondim
                       read(ilun) xx
                       do i = 1, ncache
                          chunk_pvl(ilevel)%pdata(ind+i, (iskip-1)*2+1) = xx(i)
                       end do
                       read(ilun) xx
                       do i = 1, ncache
                          chunk_pvl(ilevel)%pdata(ind+i, (iskip-1)*2+2) = xx(i)
                       end do
                    end do
                    chunk_lvl_offset(ilevel) = ind + ncache
                    nread_total = nread_total + ncache
                    do i = 1, ncache
                       do iskip = 1, twotondim
                          psi2_read = psi2_read + &
                               chunk_pvl(ilevel)%pdata(ind+i,(iskip-1)*2+1)**2 + &
                               chunk_pvl(ilevel)%pdata(ind+i,(iskip-1)*2+2)**2
                       end do
                    end do
                    deallocate(xx)
                 else
                    do i = 1, twotondim * 2
                       read(ilun)
                    end do
                 end if
              end if
           end do
        end do
        close(ilun)
     end do

     do ilevel = 1, nlevelmax
        if(varcpu_ngrid_file(ilevel) == 0) cycle

        nlocal = chunk_ngrids(ilevel)
        twotol = 2.0d0**(ilevel-1)

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
              sendbuf_2d(1, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 1)
              sendbuf_2d(2, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 2)
              sendbuf_2d(3, k) = varcpu_lvl(ilevel)%xg(xg_base + i, 3)
              sendbuf_2d(ndim+1:nprops, k) = chunk_pvl(ilevel)%pdata(k, 1:nval_per_grid)
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

        call ksection_exchange_dp(sendbuf_2d, nlocal, dest, nprops, recvbuf_2d, nrecv)
        deallocate(sendbuf_2d, dest)

        do i = 1, nrecv
           xg_recv(1:ndim) = recvbuf_2d(1:ndim, i)
           ixm = int(xg_recv(1) * twotol, 8)
           iym = int(xg_recv(2) * twotol, 8)
           izm = int(xg_recv(3) * twotol, 8)
           mkey = morton_encode(ixm, iym, izm)
           igrid = morton_hash_lookup(mort_table(ilevel), mkey)
           if(igrid == 0) then
              nskip_total = nskip_total + 1
              cycle
           end if
           nrecv_total = nrecv_total + 1

           do iskip = 1, twotondim
              icell = ICELL_OF(igrid,iskip)
              base = ndim + (iskip-1)*2
              psi_re(icell) = recvbuf_2d(base + 1, i)
              psi_im(icell) = recvbuf_2d(base + 2, i)
              psi2_recv = psi2_recv + recvbuf_2d(base+1,i)**2 + recvbuf_2d(base+2,i)**2
           end do
        end do
        deallocate(recvbuf_2d)

        if(allocated(chunk_pvl(ilevel)%pdata)) deallocate(chunk_pvl(ilevel)%pdata)
     end do
  end do

  nread_g=nread_total; nrecv_g=nrecv_total; nskip_g=nskip_total
  psi2_read_g=psi2_read; psi2_recv_g=psi2_recv
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nread_total, nread_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(nrecv_total, nrecv_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(nskip_total, nskip_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(psi2_read, psi2_read_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(psi2_recv, psi2_recv_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info2)
#endif
  if(myid==1) then
     write(*,'(A,I12,A,I12,A,I12)') &
          ' PSI_DIAG: nread=', nread_g, ' nrecv=', nrecv_g, ' nskip=', nskip_g
     write(*,'(A,ES12.5,A,ES12.5)') &
          ' PSI_DIAG: sum|psi|2_read=', psi2_read_g, ' sum|psi|2_recv=', psi2_recv_g
  end if
  if(myid==1) write(*,*) 'Binary varcpu FDM psi restore done.'

end subroutine restore_psi_binary_varcpu

subroutine restore_psi_postlb
  use amr_commons
  use poisson_commons, only: psi_re, psi_im
  use ksection
  use morton_keys
  use morton_hash
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: icpu_f, ilevel, ibound, i, j, k, ind, info, ifile
  integer :: ncache, ilun_amr, ilun_fdm, ndim2, nlevelmax2, nboundary2, ncpu2
  integer :: igrid, iskip, icell, nval_per_grid, nprops, base
  integer :: nlocal, nrecv, nrec_skip
  real(dp), allocatable :: xx(:), xxg(:)
  integer, allocatable :: iig(:)
  character(LEN=80) :: fileloc_amr, fileloc_fdm
  character(LEN=5) :: nchar, ncharcpu

  integer, allocatable :: chunk_ngrids(:), chunk_lvl_offset(:)
  integer :: nlevel_data

  type postlb_level_t
     real(dp), allocatable :: xg_file(:,:)
     real(dp), allocatable :: pdata(:,:)
  end type
  type(postlb_level_t), allocatable :: clvl(:)

  real(dp), allocatable :: sendbuf_2d(:,:), recvbuf_2d(:,:)
  integer, allocatable :: dest(:)
  integer, allocatable :: sendcnt(:), recvcnt(:), sdispls(:), rdispls(:)
  real(dp), allocatable :: sorted_buf(:,:)
  integer :: idest, ntot_recv

  integer(8) :: ixm, iym, izm
  type(mkey_t) :: mkey
  real(dp) :: xg_recv(3), xx_pos(1:nvector, 1:ndim), scale, twotol
  integer :: c_tmp(1:nvector), nx_loc

  integer :: nread_total, nrecv_total, nskip_total, info2
  integer :: nread_g, nrecv_g, nskip_g

  integer :: ncpu_amr, ndim_amr, nlevelmax_amr, nboundary_amr
  integer, allocatable :: numbl_file(:,:), numbb_file(:,:)

  ! Repair pass variables
  real(dp), allocatable :: repair_buf(:,:), repair_all(:,:)
  integer :: nrepair, nrepair_g, nrepair_matched, nrepair_global
  integer, allocatable :: repair_counts(:), repair_displs(:)
  integer :: repair_total, nrepair_matched_g

  if(.not.use_fdm) return
  if(myid==1) write(*,*) 'Post-LB FDM psi restore: ncpu_file=', ncpu_file

  nread_total = 0; nrecv_total = 0; nskip_total = 0


  ilun_amr = 98
  ilun_fdm = 99
  call title(nrestart, nchar)
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  nval_per_grid = twotondim * 2
  nprops = ndim + nval_per_grid

  allocate(chunk_ngrids(1:nlevelmax))
  allocate(chunk_lvl_offset(1:nlevelmax))
  allocate(clvl(1:nlevelmax))
  chunk_ngrids = 0
  do j = 1, varcpu_nfiles_local
     do ilevel = 1, nlevelmax
        chunk_ngrids(ilevel) = chunk_ngrids(ilevel) + &
             varcpu_nactive(varcpu_my_files(j), ilevel)
     end do
  end do

  nlevel_data = 0
  do ilevel = nlevelmax, 1, -1
     if(chunk_ngrids(ilevel) > 0) then
        nlevel_data = ilevel
        exit
     end if
  end do
  if(myid==1) then
     write(*,'(A,I3,A,I12)') &
          ' PSI_POSTLB: nlevel_data=', nlevel_data, &
          '  chunk_ngrids(max)=', chunk_ngrids(max(nlevel_data,1))
     call flush(6)
  end if

  do ilevel = 1, nlevel_data
     if(chunk_ngrids(ilevel) > 0) then
        allocate(clvl(ilevel)%xg_file(chunk_ngrids(ilevel), ndim))
        allocate(clvl(ilevel)%pdata(chunk_ngrids(ilevel), nval_per_grid))
     end if
  end do

  chunk_lvl_offset = 0
  do ifile = 1, varcpu_nfiles_local
     icpu_f = varcpu_my_files(ifile)

     if(IOGROUPSIZEREP>0) then
        call title(((icpu_f-1)/IOGROUPSIZEREP)+1, ncharcpu)
        fileloc_amr='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/amr_'//TRIM(nchar)//'.out'
        fileloc_fdm='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/fdm_'//TRIM(nchar)//'.out'
     else
        fileloc_amr='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
        fileloc_fdm='output_'//TRIM(nchar)//'/fdm_'//TRIM(nchar)//'.out'
     end if
     call title(icpu_f, ncharcpu)
     fileloc_amr=TRIM(fileloc_amr)//TRIM(ncharcpu)
     fileloc_fdm=TRIM(fileloc_fdm)//TRIM(ncharcpu)

     open(unit=ilun_fdm, file=fileloc_fdm, form='unformatted')
     read(ilun_fdm) ncpu2
     read(ilun_fdm) ndim2
     read(ilun_fdm) nlevelmax2
     read(ilun_fdm) nboundary2

     open(unit=ilun_amr, file=fileloc_amr, form='unformatted')
     do i = 1, 8; read(ilun_amr); end do
     block
       integer :: nout_tmp, dum1, dum2
       read(ilun_amr) nout_tmp, dum1, dum2
       read(ilun_amr); read(ilun_amr)
     end block
     do i = 1, 8; read(ilun_amr); end do
     read(ilun_amr); read(ilun_amr)
     ncpu_amr = ncpu2; nlevelmax_amr = nlevelmax2; nboundary_amr = nboundary2
     allocate(numbl_file(1:ncpu_amr, 1:nlevelmax_amr))
     read(ilun_amr) numbl_file
     read(ilun_amr)
     if(nboundary_amr > 0) then
        allocate(numbb_file(1:nboundary_amr, 1:nlevelmax_amr))
        do i = 1, 3; read(ilun_amr); end do
     end if
     read(ilun_amr)
     block
       character(LEN=128) :: ord_tmp
       read(ilun_amr) ord_tmp
       if(trim(ord_tmp) == 'ksection') then
          do i = 1, 10; read(ilun_amr); end do
       else if(trim(ord_tmp) == 'bisection') then
          do i = 1, 5; read(ilun_amr); end do
       else
          read(ilun_amr)
       end if
     end block
     do i = 1, 3; read(ilun_amr); end do

     do ilevel = 1, nlevelmax2
        do ibound = 1, nboundary2 + ncpu2
           read(ilun_fdm)
           read(ilun_fdm) ncache

           if(ibound <= ncpu_amr) then
              if(numbl_file(ibound, ilevel) /= ncache .and. ibound == icpu_f) then
              end if
           end if

           if(ncache > 0) then
              if(ibound == icpu_f) then
                 ind = chunk_lvl_offset(ilevel)
                 allocate(xxg(1:ncache))
                 read(ilun_amr); read(ilun_amr); read(ilun_amr)
                 do j = 1, ndim
                    read(ilun_amr) xxg
                    clvl(ilevel)%xg_file(ind+1:ind+ncache, j) = xxg(1:ncache)
                 end do
                 read(ilun_amr)
                 do i = 1, twondim; read(ilun_amr); end do
                 do i = 1, 3*twotondim; read(ilun_amr); end do
                 deallocate(xxg)

                 allocate(xx(1:ncache))
                 do iskip = 1, twotondim
                    read(ilun_fdm) xx
                    do i = 1, ncache
                       clvl(ilevel)%pdata(ind+i, (iskip-1)*2+1) = xx(i)
                    end do
                    read(ilun_fdm) xx
                    do i = 1, ncache
                       clvl(ilevel)%pdata(ind+i, (iskip-1)*2+2) = xx(i)
                    end do
                 end do
                 deallocate(xx)
                 chunk_lvl_offset(ilevel) = ind + ncache
                 nread_total = nread_total + ncache
              else
                 nrec_skip = 4 + ndim + twondim + 3*twotondim
                 do i = 1, nrec_skip; read(ilun_amr); end do
                 do i = 1, twotondim * 2; read(ilun_fdm); end do
              end if
           end if
        end do
     end do
     close(ilun_amr)
     close(ilun_fdm)
     deallocate(numbl_file)
     if(allocated(numbb_file)) deallocate(numbb_file)
  end do

  ! Cap exchange at highest level with actual grids in post-LB mesh
  do ilevel = nlevel_data, 1, -1
     if(numbtot(1, ilevel) > 0) exit
  end do
  if(ilevel < nlevel_data) then
     if(myid==1) then
        write(*,'(A,I3,A,I3)') &
             ' PSI_POSTLB: capping exchange at level', ilevel, &
             ' (file has', nlevel_data, ')'
        call flush(6)
     end if
     do j = ilevel+1, nlevel_data
        if(allocated(clvl(j)%xg_file)) deallocate(clvl(j)%xg_file)
        if(allocated(clvl(j)%pdata)) deallocate(clvl(j)%pdata)
     end do
     nlevel_data = ilevel
  end if
#ifndef WITHOUTMPI
  ! All ranks must call the same collectives at each level.
  ! nlevel_data is rank-local (depends on file assignment), so synchronize.
  call MPI_ALLREDUCE(MPI_IN_PLACE, nlevel_data, 1, &
       MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, info)
#endif
  if(myid==1) then
     write(*,*) 'PSI_POSTLB: file reading done, starting exchange'
     call flush(6)
  end if
  do ilevel = 1, nlevel_data
     nlocal = chunk_ngrids(ilevel)

     if(myid==1) then
        write(*,'(A,I3,A,I12)') ' PSI_POSTLB: exchange lv=', ilevel, ' n=', nlocal
        call flush(6)
     end if
     twotol = 2.0d0**(ilevel-1)
     allocate(sendbuf_2d(1:nprops, 1:max(nlocal,1)))
     allocate(dest(max(nlocal,1)))
     do i = 1, nlocal
        sendbuf_2d(1, i) = clvl(ilevel)%xg_file(i, 1)
        sendbuf_2d(2, i) = clvl(ilevel)%xg_file(i, 2)
        sendbuf_2d(3, i) = clvl(ilevel)%xg_file(i, 3)
        sendbuf_2d(ndim+1:nprops, i) = clvl(ilevel)%pdata(i, 1:nval_per_grid)
        xx_pos(1,1) = (sendbuf_2d(1, i) - dble(icoarse_min)) * scale
        xx_pos(1,2) = (sendbuf_2d(2, i) - dble(jcoarse_min)) * scale
        xx_pos(1,3) = (sendbuf_2d(3, i) - dble(kcoarse_min)) * scale
        if(ordering == 'ksection') then
           call cmp_ksection_cpumap(xx_pos, c_tmp, 1)
        else
           call cmp_cpumap(xx_pos, c_tmp, 1)
        end if
        dest(i) = c_tmp(1)
     end do
     if(allocated(clvl(ilevel)%xg_file)) deallocate(clvl(ilevel)%xg_file)
     if(allocated(clvl(ilevel)%pdata))   deallocate(clvl(ilevel)%pdata)

     if(myid==1) then
        write(*,'(A,I3,A)') ' PSI_POSTLB: lv=', ilevel, ' dest done, sorting'
        call flush(6)
     end if

     ! --- MPI_Alltoallv exchange ---
     allocate(sendcnt(0:ncpu-1), recvcnt(0:ncpu-1))
     allocate(sdispls(0:ncpu-1), rdispls(0:ncpu-1))
     sendcnt = 0
     do i = 1, nlocal
        idest = dest(i) - 1
        sendcnt(idest) = sendcnt(idest) + 1
     end do
     allocate(sorted_buf(1:nprops, 1:max(nlocal,1)))
     sdispls(0) = 0
     do i = 1, ncpu-1
        sdispls(i) = sdispls(i-1) + sendcnt(i-1)
     end do
     do i = 1, nlocal
        idest = dest(i) - 1
        sdispls(idest) = sdispls(idest) + 1
        sorted_buf(1:nprops, sdispls(idest)) = sendbuf_2d(1:nprops, i)
     end do
     sdispls(0) = 0
     do i = 1, ncpu-1
        sdispls(i) = sdispls(i-1) + sendcnt(i-1)
     end do
     deallocate(sendbuf_2d, dest)

     call MPI_ALLTOALL(sendcnt, 1, MPI_INTEGER, &
                       recvcnt, 1, MPI_INTEGER, MPI_COMM_WORLD, info)
     ntot_recv = 0
     rdispls(0) = 0
     do i = 0, ncpu-1
        ntot_recv = ntot_recv + recvcnt(i)
        if(i > 0) rdispls(i) = rdispls(i-1) + recvcnt(i-1)
     end do

     do i = 0, ncpu-1
        sendcnt(i) = sendcnt(i) * nprops
        recvcnt(i) = recvcnt(i) * nprops
        sdispls(i) = sdispls(i) * nprops
        rdispls(i) = rdispls(i) * nprops
     end do
     allocate(recvbuf_2d(1:nprops, 1:max(ntot_recv,1)))
     call MPI_ALLTOALLV(sorted_buf, sendcnt, sdispls, MPI_DOUBLE_PRECISION, &
                        recvbuf_2d, recvcnt, rdispls, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, info)
     deallocate(sorted_buf, sendcnt, recvcnt, sdispls, rdispls)
     nrecv = ntot_recv

     if(myid==1) then
        write(*,'(A,I3,A,I12)') ' PSI_POSTLB: lv=', ilevel, ' exchange done, nrecv=', nrecv
        call flush(6)
     end if

     ! --- Scatter to local grids + collect misrouted items ---
     allocate(repair_buf(1:nprops, 1:max(nrecv/6, 1)))
     nrepair = 0
     do i = 1, nrecv
        xg_recv(1:ndim) = recvbuf_2d(1:ndim, i)
        ixm = int(xg_recv(1) * twotol, 8)
        iym = int(xg_recv(2) * twotol, 8)
        izm = int(xg_recv(3) * twotol, 8)
        mkey = morton_encode(ixm, iym, izm)
        igrid = morton_hash_lookup(mort_table(ilevel), mkey)
        if(igrid <= 0 .or. igrid > ngridmax) then
           nskip_total = nskip_total + 1
           nrepair = nrepair + 1
           if(nrepair > size(repair_buf, 2)) then
              call repair_buf_grow(repair_buf, nprops, nrepair)
           end if
           repair_buf(1:nprops, nrepair) = recvbuf_2d(1:nprops, i)
           cycle
        end if
        nrecv_total = nrecv_total + 1
        do iskip = 1, twotondim
           icell = ICELL_OF(igrid,iskip)
           base = ndim + (iskip-1)*2
           psi_re(icell) = recvbuf_2d(base + 1, i)
           psi_im(icell) = recvbuf_2d(base + 2, i)
        end do
     end do
     deallocate(recvbuf_2d)

     ! --- Repair pass: Allgatherv misrouted items ---
#ifndef WITHOUTMPI
     ! Sync nrepair across all ranks so conditional is uniform
     call MPI_ALLREDUCE(nrepair, nrepair_global, 1, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
     if(nrepair_global > 0) then
        allocate(repair_counts(0:ncpu-1), repair_displs(0:ncpu-1))
        call MPI_ALLGATHER(nrepair, 1, MPI_INTEGER, &
             repair_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, info)
        repair_total = sum(repair_counts)

        repair_displs(0) = 0
        do i = 1, ncpu-1
           repair_displs(i) = repair_displs(i-1) + repair_counts(i-1)
        end do
        allocate(repair_all(1:nprops, 1:repair_total))
        ! Scale for dp Allgatherv
        repair_counts = repair_counts * nprops
        repair_displs = repair_displs * nprops
        call MPI_ALLGATHERV(repair_buf, nrepair*nprops, &
             MPI_DOUBLE_PRECISION, repair_all, repair_counts, &
             repair_displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
        ! Restore counts (unscale)
        repair_counts = repair_counts / nprops
        repair_displs = repair_displs / nprops
        repair_total = sum(repair_counts)

        if(myid==1) then
           write(*,'(A,I3,A,I10)') &
                ' PSI_POSTLB: lv=', ilevel, ' repair pass, ntotal=', repair_total
           call flush(6)
        end if

        nrepair_matched = 0
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(xg_recv, ixm, iym, izm, mkey, igrid, iskip, icell, base) &
        !$OMP REDUCTION(+:nrepair_matched) &
        !$OMP SCHEDULE(STATIC, 4096)
        do i = 1, repair_total
           xg_recv(1:ndim) = repair_all(1:ndim, i)
           ixm = int(xg_recv(1) * twotol, 8)
           iym = int(xg_recv(2) * twotol, 8)
           izm = int(xg_recv(3) * twotol, 8)
           mkey = morton_encode(ixm, iym, izm)
           igrid = morton_hash_lookup(mort_table(ilevel), mkey)
           if(igrid <= 0 .or. igrid > ngridmax) cycle
           nrepair_matched = nrepair_matched + 1
           do iskip = 1, twotondim
              icell = ICELL_OF(igrid,iskip)
              base = ndim + (iskip-1)*2
              psi_re(icell) = repair_all(base + 1, i)
              psi_im(icell) = repair_all(base + 2, i)
           end do
        end do
        !$OMP END PARALLEL DO

        nrecv_total = nrecv_total + nrepair_matched
        deallocate(repair_all)

        call MPI_ALLREDUCE(nrepair_matched, nrepair_matched_g, 1, &
             MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
        if(myid==1) then
           write(*,'(A,I3,A,I10)') &
                ' PSI_POSTLB: lv=', ilevel, ' repair matched=', nrepair_matched_g
           call flush(6)
        end if
        deallocate(repair_counts, repair_displs)
     end if
#endif
     deallocate(repair_buf)
  end do

  nread_g=nread_total; nrecv_g=nrecv_total; nskip_g=nskip_total
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nread_total, nread_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(nrecv_total, nrecv_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
  call MPI_ALLREDUCE(nskip_total, nskip_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info2)
#endif
  if(myid==1) then
     write(*,'(A,I12,A,I12,A,I12)') &
          ' PSI_POSTLB: nread=', nread_g, ' nrecv=', nrecv_g, ' nskip=', nskip_g
     call flush(6)
  end if
  deallocate(chunk_ngrids, chunk_lvl_offset, clvl)
  if(myid==1) then
     write(*,*) 'Post-LB FDM psi restore done.'
     call flush(6)
  end if

contains
  subroutine repair_buf_grow(buf, np, needed)
    real(dp), allocatable, intent(inout) :: buf(:,:)
    integer, intent(in) :: np, needed
    real(dp), allocatable :: tmp(:,:)
    integer :: old_cap, new_cap
    old_cap = size(buf, 2)
    new_cap = max(2 * old_cap, needed)
    allocate(tmp(1:np, 1:new_cap))
    tmp(1:np, 1:old_cap) = buf(1:np, 1:old_cap)
    call move_alloc(tmp, buf)
  end subroutine

end subroutine restore_psi_postlb

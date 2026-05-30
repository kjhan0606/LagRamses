subroutine backup_psi(filename)
  use amr_commons
  use poisson_commons, only: psi_re, psi_im
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=80)::filename

  integer::i,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound
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
              iskip=ncoarse+(ind-1)*ngridmax
              ! Write Re(psi)
              do i=1,ncache
                 xdp(i)=psi_re(ind_grid(i)+iskip)
              end do
              write(ilun)xdp
              ! Write Im(psi)
              do i=1,ncache
                 xdp(i)=psi_im(ind_grid(i)+iskip)
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ncache,iskip,igrid,i,ilevel,ind,ilun,ibound,istart,info
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
              iskip=ncoarse+(ind-1)*ngridmax
              ! Read Re(psi)
              read(ilun)xx
              do i=1,ncache
                 psi_re(ind_grid(i)+iskip)=xx(i)
              end do
              ! Read Im(psi)
              read(ilun)xx
              do i=1,ncache
                 psi_im(ind_grid(i)+iskip)=xx(i)
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

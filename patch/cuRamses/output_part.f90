subroutine backup_part(filename)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif 
  character(LEN=80)::filename

  integer::i,idim,ilun,ipart,ielt
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  integer(kind=1),allocatable,dimension(:)::l1
  logical,allocatable,dimension(:)::nb
  integer,parameter::tag=1122
  integer::dummy_io,info2
  integer :: unit_info, ivar, idesc
  character(len=80) :: descriptor

  if(verbose)write(*,*)'Entering backup_part'
  
  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZEOUT>0) then
     if (mod(myid-1,IOGROUPSIZEOUT)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif


  ilun=2*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)npart
  write(ilun)localseed
  write(ilun)nstar_tot   
  write(ilun)mstar_tot   
  write(ilun)mstar_lost
  write(ilun)nsink
  ! Write position
  allocate(xdp(1:npart))
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=xp(i,idim)
         end if
     end do
     write(ilun)xdp
  end do
  ! Write velocity
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=vp(i,idim)
        end if
     end do
     write(ilun)xdp
  end do
  ! Write mass
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        xdp(ipart)=mp(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii8(ipart)=idp(i)
     end if
  end do
  write(ilun)ii8
  deallocate(ii8)
  ! Write level
  allocate(ll(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ll(ipart)=levelp(i)
     end if
  end do
  write(ilun)ll
  deallocate(ll)

  ! Write particle type
  allocate(l1(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        l1(ipart)=ptypep(i)
     end if
  end do
  write(ilun)l1
  deallocate(l1)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  ipart=0
  do i=1, npartmax
     if(levelp(i)>0) then
        ipart=ipart+1
        xdp(ipart)=ptcl_phi(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
#endif

  ! Write birth epoch
  if(star.or.sink)then
     allocate(xdp(1:npart))
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=tp(i)
        end if
     end do
     write(ilun)xdp
     ! Write metallicity
     if(metal)then
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=zp(i)
           end if
        end do
        write(ilun)xdp
     end if
     
     ! Write birth epoch (proper time)                                                                                   
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=tpp(i)
        end if
     end do
     write(ilun)xdp
     ! Write mass (initial)                                                                                              
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=mp0(i)
        end if
     end do
     write(ilun)xdp
     ! Write indtab (checkpoint in yield table)                                                                                                         
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=indtab(i)
        end if
     end do
     write(ilun)xdp
     ! Write metal abundances                                                                                            
!     do ielt=1,nelt
!        ipart=0
!        do i=1,npartmax
!           if(levelp(i)>0)then
!              ipart=ipart+1
!              xdp(ipart)=cep(i,ielt)
!           end if
!        end do
!        write(ilun)xdp
!     enddo
     
     deallocate(xdp)
  end if
  close(ilun)

  ! Keep the native particle stream self-describing.  The descriptor records
  ! the exact payload order used above, including the stellar release cursor;
  ! restart/evidence tooling must derive binary record positions from this
  ! file rather than duplicating an implicit Python or source-code order.
  if(myid==1) then
     idesc=index(trim(filename),'part_')
     if(idesc>1) then
        descriptor=trim(filename(1:idesc-1))//'part_file_descriptor.txt'
        open(newunit=unit_info,file=trim(descriptor),status='replace',form='formatted')
        write(unit_info,'("# version: ",i2)') 1
        write(unit_info,'("# ",a,", ",a,", ",a)') 'ivar','variable_name','variable_type'
        ivar=1
        do idim=1,ndim
           write(unit_info,'(I2,", ",a,", d")') ivar, 'position_'//char(iachar('x')+idim-1)
           ivar=ivar+1
        end do
        do idim=1,ndim
           write(unit_info,'(I2,", ",a,", d")') ivar, 'velocity_'//char(iachar('x')+idim-1)
           ivar=ivar+1
        end do
        write(unit_info,'(I2,", mass, d")') ivar
        ivar=ivar+1
        write(unit_info,'(I2,", identity, q")') ivar
        ivar=ivar+1
        write(unit_info,'(I2,", levelp, i")') ivar
        ivar=ivar+1
        write(unit_info,'(I2,", ptypep, b")') ivar
        ivar=ivar+1
#ifdef OUTPUT_PARTICLE_POTENTIAL
        write(unit_info,'(I2,", potential, d")') ivar
        ivar=ivar+1
#endif
        if(star.or.sink)then
           write(unit_info,'(I2,", birth_epoch, d")') ivar
           ivar=ivar+1
           if(metal)then
              write(unit_info,'(I2,", metallicity, d")') ivar
              ivar=ivar+1
           end if
           write(unit_info,'(I2,", tpp, d")') ivar
           ivar=ivar+1
           write(unit_info,'(I2,", mp0, d")') ivar
           ivar=ivar+1
           write(unit_info,'(I2,", indtab, d")') ivar
        end if
        close(unit_info)
     end if
  end if



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
  

end subroutine backup_part

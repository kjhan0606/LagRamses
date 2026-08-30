subroutine init_sink
  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::idim,ilevel
  integer::i,isink
  integer::ilun,nx_loc
  integer::nsinkold
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp)::sm2,dmf
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical,allocatable,dimension(:)::nb
  logical::eof,ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  integer,parameter::tag=1112,tag2=1113
  integer::dummy_io,info2


  allocate(total_volume(1:nsinkmax))
  allocate(wdens(1:nsinkmax))
  allocate(wvol(1:nsinkmax))
  allocate(wmom(1:nsinkmax,1:ndim))
  allocate(wc2(1:nsinkmax))
  allocate(wdens_new(1:nsinkmax))
  allocate(wvol_new(1:nsinkmax))
  allocate(wmom_new(1:nsinkmax,1:ndim))
  allocate(wc2_new(1:nsinkmax))
  allocate(msink(1:nsinkmax))
  allocate(msink_new(1:nsinkmax))
  allocate(msink_all(1:nsinkmax))
  allocate(idsink(1:nsinkmax))
  ! Important to set nindsink
  idsink=0
  allocate(idsink_new(1:nsinkmax))
  allocate(idsink_all(1:nsinkmax))
  allocate(tsink(1:nsinkmax))
  allocate(tsink_new(1:nsinkmax))
  allocate(tsink_all(1:nsinkmax))
  allocate(vsink(1:nsinkmax,1:ndim))
  allocate(xsink(1:nsinkmax,1:ndim))
  allocate(vsink_new(1:nsinkmax,1:ndim))
  allocate(vsink_all(1:nsinkmax,1:ndim))
  allocate(xsink_new(1:nsinkmax,1:ndim))
  allocate(xsink_all(1:nsinkmax,1:ndim))
  allocate(dMBHoverdt(1:nsinkmax))
  allocate(dMEdoverdt(1:nsinkmax))
  allocate(r2sink(1:nsinkmax))
  allocate(r2k(1:nsinkmax))
  allocate(v2sink(1:nsinkmax))
  allocate(c2sink(1:nsinkmax))
  allocate(v2sink_new(1:nsinkmax))
  allocate(c2sink_new(1:nsinkmax))
  allocate(v2sink_all(1:nsinkmax))
  allocate(c2sink_all(1:nsinkmax))
  allocate(weighted_density(1:nsinkmax,1:nlevelmax))
  allocate(weighted_volume (1:nsinkmax,1:nlevelmax))
  allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
  allocate(weighted_c2 (1:nsinkmax,1:nlevelmax))
  allocate(oksink_new(1:nsinkmax))
  allocate(oksink_all(1:nsinkmax))
  allocate(jsink(1:nsinkmax,1:ndim))
  allocate(jsink_new(1:nsinkmax,1:ndim))
  allocate(jsink_all(1:nsinkmax,1:ndim))
  allocate(dMBH_coarse    (1:nsinkmax))
  allocate(dMEd_coarse    (1:nsinkmax))
  allocate(dMsmbh         (1:nsinkmax))
  allocate(Esave          (1:nsinkmax))
  allocate(dMBH_coarse_new(1:nsinkmax))
  allocate(dMEd_coarse_new(1:nsinkmax))
  allocate(dMsmbh_new     (1:nsinkmax))
  allocate(Esave_new      (1:nsinkmax))
  allocate(dMBH_coarse_all(1:nsinkmax))
  allocate(dMEd_coarse_all(1:nsinkmax))
  allocate(dMsmbh_all     (1:nsinkmax))
  allocate(Esave_all      (1:nsinkmax))
  allocate(sink_stat      (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(sink_stat_all  (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(v_avgptr(1:nsinkmax))
  allocate(c_avgptr(1:nsinkmax))
  allocate(d_avgptr(1:nsinkmax))
  allocate(spinmag(1:nsinkmax),bhspin(1:nsinkmax,1:ndim))
  allocate(spinmag_new(1:nsinkmax),bhspin_new(1:nsinkmax,1:ndim))
  allocate(spinmag_all(1:nsinkmax),bhspin_all(1:nsinkmax,1:ndim))
  allocate(eps_sink(1:nsinkmax))
  eps_sink=0.057190958d0

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(nrestart>0)then
     ilun=4*ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/sink_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/sink_'//TRIM(nchar)//'.out'
     endif


!     call title(myid,nchar)
     fileloc=TRIM(fileloc)!//TRIM(nchar)

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
     rewind(ilun)
     read(ilun)nsink
     read(ilun)nindsink

     if(nsink>0)then
        allocate(xdp(1:nsink))
        allocate(isp(1:nsink))
        read(ilun)isp
        idsink(1:nsink)=isp
        ! Important for the indexation of sinks
        nindsink=MAXVAL(idsink)
        deallocate(isp)
        read(ilun)xdp
        msink(1:nsink)=xdp
        do idim=1,ndim
           read(ilun)xdp
           xsink(1:nsink,idim)=xdp
        end do
        do idim=1,ndim
           read(ilun)xdp
           vsink(1:nsink,idim)=xdp
        end do
        read(ilun)xdp
        tsink(1:nsink)=xdp
        read(ilun)xdp
        dMsmbh(1:nsink)=xdp
        read(ilun)xdp
        dMBH_coarse(1:nsink)=xdp
        read(ilun)xdp
        dMEd_coarse(1:nsink)=xdp
        read(ilun)xdp
        Esave(1:nsink)=xdp
        do idim=1,ndim
           read(ilun)xdp
           jsink(1:nsink,idim)=xdp
        end do
        do idim=1,ndim
           read(ilun)xdp
           bhspin(1:nsink,idim)=xdp
        end do
        read(ilun)xdp
        spinmag(1:nsink)=xdp
        read(ilun)xdp
        eps_sink(1:nsink)=xdp
        do idim=1,ndim*2+1
           do ilevel=levelmin,nlevelmax
              read(ilun)xdp
              sink_stat(1:nsink,ilevel,idim)=xdp
           enddo
        enddo
        deallocate(xdp)
     end if
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

  end if

  if (nrestart>0)then
     nsinkold=nsink  
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink_restart'
     else
        filename='ic_sink_restart'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink_restart: ',filename
     if (.not. ic_sink)then
        filename='ic_sink_restart'
        INQUIRE(FILE=filename, EXIST=ic_sink)
     end if
  else
     nsink=0
     nindsink=0
     nsinkold=0
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink'
     else
        filename='ic_sink'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink: ',filename
     if (.not. ic_sink)then
        filename='ic_sink'
        INQUIRE(FILE=filename, EXIST=ic_sink)
     end if
  end if

  if (ic_sink)then

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag2,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(10,file=filename,form='formatted')
     do
        read(10,*,end=103)mm1,xx1,xx2,xx3,vv1,vv2,vv3,ll1,ll2,ll3,sm2,dmf
        if(nsink>=nsinkmax)then
           if(myid==1)write(*,*)'init_sink: ic_sink holds more than nsinkmax sinks'
           call clean_stop
        end if
        nsink=nsink+1
        nindsink=nindsink+1
        idsink(nsink)=nindsink
        msink(nsink)=mm1
        xsink(nsink,1)=xx1+boxlen/2
        xsink(nsink,2)=xx2+boxlen/2
        xsink(nsink,3)=xx3+boxlen/2
        vsink(nsink,1)=vv1
        vsink(nsink,2)=vv2
        vsink(nsink,3)=vv3
        jsink(nsink,1)=ll1
        jsink(nsink,2)=ll2
        jsink(nsink,3)=ll3
        tsink(nsink)=t
        dMsmbh(nsink)=0d0
        dMBH_coarse(nsink)=0d0
        dMEd_coarse(nsink)=0d0
        Esave(nsink)=0d0
        bhspin(nsink,1:ndim)=0d0
        spinmag(nsink)=0d0
        sink_stat(nsink,levelmin:nlevelmax,1:ndim*2+1)=0d0
     end do
103  continue
     close(10)
     if(myid==1)write(*,*)'init_sink: loaded ',nsink,' sinks from ',TRIM(filename)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag2, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

  end if

end subroutine init_sink

!-------------------------------------------------------
! Allocate sink arrays only (no file I/O).
! Called before restore_part_hdf5 so that sink arrays
! exist before HDF5 restore writes into them.
!-------------------------------------------------------
subroutine init_sink_alloc
  use amr_commons
  use pm_commons
  implicit none

  allocate(total_volume(1:nsinkmax))
  allocate(wdens(1:nsinkmax))
  allocate(wvol(1:nsinkmax))
  allocate(wmom(1:nsinkmax,1:ndim))
  allocate(wc2(1:nsinkmax))
  allocate(wdens_new(1:nsinkmax))
  allocate(wvol_new(1:nsinkmax))
  allocate(wmom_new(1:nsinkmax,1:ndim))
  allocate(wc2_new(1:nsinkmax))
  allocate(msink(1:nsinkmax))
  allocate(msink_new(1:nsinkmax))
  allocate(msink_all(1:nsinkmax))
  allocate(idsink(1:nsinkmax))
  idsink=0
  allocate(idsink_new(1:nsinkmax))
  allocate(idsink_all(1:nsinkmax))
  allocate(tsink(1:nsinkmax))
  allocate(tsink_new(1:nsinkmax))
  allocate(tsink_all(1:nsinkmax))
  allocate(vsink(1:nsinkmax,1:ndim))
  allocate(xsink(1:nsinkmax,1:ndim))
  allocate(vsink_new(1:nsinkmax,1:ndim))
  allocate(vsink_all(1:nsinkmax,1:ndim))
  allocate(xsink_new(1:nsinkmax,1:ndim))
  allocate(xsink_all(1:nsinkmax,1:ndim))
  allocate(dMBHoverdt(1:nsinkmax))
  allocate(dMEdoverdt(1:nsinkmax))
  allocate(r2sink(1:nsinkmax))
  allocate(r2k(1:nsinkmax))
  allocate(v2sink(1:nsinkmax))
  allocate(c2sink(1:nsinkmax))
  allocate(v2sink_new(1:nsinkmax))
  allocate(c2sink_new(1:nsinkmax))
  allocate(v2sink_all(1:nsinkmax))
  allocate(c2sink_all(1:nsinkmax))
  allocate(weighted_density(1:nsinkmax,1:nlevelmax))
  allocate(weighted_volume (1:nsinkmax,1:nlevelmax))
  allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
  allocate(weighted_c2 (1:nsinkmax,1:nlevelmax))
  allocate(oksink_new(1:nsinkmax))
  allocate(oksink_all(1:nsinkmax))
  allocate(jsink(1:nsinkmax,1:ndim))
  allocate(jsink_new(1:nsinkmax,1:ndim))
  allocate(jsink_all(1:nsinkmax,1:ndim))
  allocate(dMBH_coarse    (1:nsinkmax))
  allocate(dMEd_coarse    (1:nsinkmax))
  allocate(dMsmbh         (1:nsinkmax))
  allocate(Esave          (1:nsinkmax))
  allocate(dMBH_coarse_new(1:nsinkmax))
  allocate(dMEd_coarse_new(1:nsinkmax))
  allocate(dMsmbh_new     (1:nsinkmax))
  allocate(Esave_new      (1:nsinkmax))
  allocate(dMBH_coarse_all(1:nsinkmax))
  allocate(dMEd_coarse_all(1:nsinkmax))
  allocate(dMsmbh_all     (1:nsinkmax))
  allocate(Esave_all      (1:nsinkmax))
  allocate(sink_stat      (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(sink_stat_all  (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(v_avgptr(1:nsinkmax))
  allocate(c_avgptr(1:nsinkmax))
  allocate(d_avgptr(1:nsinkmax))
  allocate(spinmag(1:nsinkmax),bhspin(1:nsinkmax,1:ndim))
  allocate(spinmag_new(1:nsinkmax),bhspin_new(1:nsinkmax,1:ndim))
  allocate(spinmag_all(1:nsinkmax),bhspin_all(1:nsinkmax,1:ndim))
  allocate(eps_sink(1:nsinkmax))
  eps_sink=0.057190958d0

end subroutine init_sink_alloc

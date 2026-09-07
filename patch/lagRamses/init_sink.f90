! Project-local override of patch/cuRamses/init_sink.f90: initialize the
! accepted AGN reservoirs (persisted in HDF5; legacy binary layout unchanged).
subroutine init_sink
  use amr_commons
  use pm_commons
  use clfind_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
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
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical,allocatable,dimension(:)::nb
  logical::eof,ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  integer,parameter::tag=1112,tag2=1113
  integer::dummy_io,info2
  integer::ic_status,ic_extra_status,ic_line,ic_owner(nvector),ic_part
  real(dp)::ic_row(12),ic_position(nvector,ndim)
  character(len=2048)::ic_record
  character(len=64)::ic_extra


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
  allocate(agn_pending_erg(1:nsinkmax))
  agn_pending_erg=0d0
  allocate(agn_mechanical_pending(4,nsinkmax))
  agn_mechanical_pending=0d0
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

  ! Read the established 12-column sink IC layout in code units:
  ! mass, centered x/y/z, vx/vy/vz, gas angular momentum x/y/z,
  ! upstream-only SMBH mass and drag fraction. This patch has one BH mass;
  ! reject nonzero upstream-only fields rather than silently discarding them.
  ! The old override only INQUIREd this file and never created its sinks.
#ifndef WITHOUTMPI
  call MPI_BCAST(ic_sink,1,MPI_LOGICAL,0,MPI_COMM_WORLD,info2)
#endif
  if(ic_sink)then
     if(ndim/=3)then
        if(myid==1)write(*,*)'ERROR: ic_sink loading requires NDIM=3'
        call clean_stop
     endif
     ic_status=0
     if(myid==1)open(newunit=ilun,file=trim(filename),status='old',action='read',iostat=ic_status)
#ifndef WITHOUTMPI
     call MPI_BCAST(ic_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,info2)
#endif
     if(ic_status/=0)then
        if(myid==1)write(*,*)'ERROR: cannot open sink IC ',trim(filename)
        call clean_stop
     endif
     ic_line=0
     do
        ic_status=0
        if(myid==1)then
           do
              read(ilun,'(A)',iostat=ic_status)ic_record
              ic_line=ic_line+1
              if(ic_status/=0)exit
              ic_record=adjustl(ic_record)
              if(len_trim(ic_record)==0)cycle
              if(ic_record(1:1)=='!'.or.ic_record(1:1)=='#')cycle
              exit
           enddo
           if(ic_status<0)then
              ic_status=-1
           else if(ic_status==0)then
              ic_row=ieee_value(0d0,ieee_quiet_nan)
              read(ic_record,*,iostat=ic_status)ic_row
              if(ic_status/=0)ic_status=1 ! An incomplete row is not file EOF.
              if(ic_status==0)then
                 read(ic_record,*,iostat=ic_extra_status)ic_row,ic_extra
                 if(ic_extra_status==0)ic_status=1
                 if(.not.all(ieee_is_finite(ic_row)))ic_status=1
                 if(ic_row(1)<=0d0.or.any(ic_row(11:12)/=0d0))ic_status=1
                 if(any(ic_row(2:4)<-boxlen/2).or.any(ic_row(2:4)>=boxlen/2))ic_status=1
              endif
           endif
        endif
#ifndef WITHOUTMPI
        call MPI_BCAST(ic_status,1,MPI_INTEGER,0,MPI_COMM_WORLD,info2)
#endif
        if(ic_status==-1)exit
        if(ic_status/=0.or.nsink>=nsinkmax.or.nindsink==huge(nindsink))then
           if(myid==1)write(*,*)'ERROR: invalid or over-capacity sink IC at line ',ic_line
           call clean_stop
        endif
#ifndef WITHOUTMPI
        call MPI_BCAST(ic_row,12,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info2)
#endif
        nsink=nsink+1
        nindsink=nindsink+1
        idsink(nsink)=nindsink
        msink(nsink)=ic_row(1)
        xsink(nsink,:)=ic_row(2:4)+boxlen/2
        vsink(nsink,:)=ic_row(5:7)
        jsink(nsink,:)=ic_row(8:10)
        tsink(nsink)=t
        dMsmbh(nsink)=0d0; dMBH_coarse(nsink)=0d0; dMEd_coarse(nsink)=0d0
        Esave(nsink)=0d0; agn_pending_erg(nsink)=0d0
        agn_mechanical_pending(:,nsink)=0d0
        bhspin(nsink,:)=0d0; spinmag(nsink)=0d0
        sink_stat(nsink,:,:)=0d0

        ! init_part calls us immediately before init_tree: append one
        ! canonical particle on its owner, not just the replicated BH table.
        ! Normal create_sink subsequently builds its accretion cloud.
        ic_position=0d0
        ic_position(1,:)=xsink(nsink,:)
        ic_owner(1)=1
#ifndef WITHOUTMPI
        call cmp_cpumap(ic_position,ic_owner,1)
#endif
        if(ic_owner(1)==myid)then
           ic_part=npart+1
           if(ic_part>npartmax)then
              if(.not.npartmax_auto)then
                 write(*,*)'ERROR: insufficient particle capacity for sink IC'
                 call clean_stop
              endif
              call grow_particle_bundle(ic_part)
           endif
           npart=ic_part
           xp(ic_part,:)=xsink(nsink,:); vp(ic_part,:)=vsink(nsink,:)
           mp(ic_part)=msink(nsink); idp(ic_part)=-int(nsink,i8b)
           ptypep(ic_part)=PTYPE_SINK; levelp(ic_part)=levelmin
           tp(ic_part)=0d0
           if(allocated(zp))zp(ic_part)=0d0
           if(allocated(tpp))tpp(ic_part)=0d0
           if(allocated(mp0))mp0(ic_part)=0d0
           if(allocated(indtab))indtab(ic_part)=0d0
        endif
     enddo
     if(myid==1)then
        close(ilun)
        write(*,'(A,I0,A,A)')' Sink IC loaded: added=',nsink-nsinkold,' file=',trim(filename)
     endif
  endif

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
  allocate(agn_pending_erg(1:nsinkmax))
  agn_pending_erg=0d0
  allocate(agn_mechanical_pending(4,nsinkmax))
  agn_mechanical_pending=0d0
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

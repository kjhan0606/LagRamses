subroutine init_sink
  use amr_commons
  use pm_commons
  use clfind_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: iostat_end
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::idim,ilevel
  integer::ic_unit,ic_status
  integer::sink_stat_marker,source_ncpu
  integer,parameter::sink_stat_global_marker=20261003
  real(dp)::ic_values(12)
  real(dp)::seed_pos(1:nvector,1:ndim)
  integer::seed_cpu(1:nvector)
  integer::i,isink
  integer::ilun,nx_loc
  integer::nsinkold
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical,allocatable,dimension(:)::nb
  logical::eof,ic_sink=.false.,sink_file_exists
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
  weighted_density=0d0
  weighted_volume=0d0
  weighted_momentum=0d0
  weighted_c2=0d0
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
     inquire(file=trim(fileloc),exist=sink_file_exists)
     if(.not.sink_file_exists.and.IOGROUPSIZEREP>0)then
        ! backup_sink writes the shared SMBH state only from rank 1,
        ! including when the particle shards use multiple I/O groups.
        fileloc='output_'//trim(nchar)//'/group_00001/sink_'//trim(nchar)//'.out'
        inquire(file=trim(fileloc),exist=sink_file_exists)
     endif
     if(.not.sink_file_exists)then
        if(myid==1)write(*,*)'ERROR: missing SMBH checkpoint ',trim(fileloc)
        call clean_stop
     endif

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
        read(ilun,iostat=ic_status)sink_stat_marker
        if(ic_status==0)then
           if(sink_stat_marker/=sink_stat_global_marker)then
              if(myid==1)write(*,*)'ERROR: unsupported SMBH sink-stat checkpoint format'
              call clean_stop
           endif
           ! The checkpoint stores the global sum once.  Assign it to one
           ! rank so the first post-restart MPI reduction recovers it exactly.
           if(myid>1)sink_stat(1:nsink,levelmin:nlevelmax,:)=0d0
        else if(ic_status==iostat_end)then
           ! Older multi-rank files contain only rank 1's local statistic.
           ! Its missing rank contributions cannot be reconstructed safely.
           source_ncpu=ncpu
           if(varcpu_restart.and.ncpu_file>0)source_ncpu=ncpu_file
           if(source_ncpu>1)then
              if(myid==1)write(*,*)'ERROR: legacy MPI sink checkpoint lacks global sink statistics'
              call clean_stop
           endif
           ! A one-rank legacy source has a complete statistic, even when
           ! restored with more ranks. Seed it once before MPI reduction.
           if(myid>1)sink_stat(1:nsink,levelmin:nlevelmax,:)=0d0
        else
           if(myid==1)write(*,*)'ERROR: cannot read SMBH sink-stat format marker'
           call clean_stop
        endif
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

  ! The cuRamses initializer previously discovered ic_sink but never read it.
  ! The VPATH-selected SMBH implementation accepts the established twelve
  ! column sink seed format (mass, position, velocity, gas angular momentum,
  ! BH mass, dark-envelope mass). Every rank reads the same small seed file.
  if(smbh .and. nrestart==0 .and. ic_sink) then
     open(newunit=ic_unit,file=trim(filename),status='old',action='read', &
          iostat=ic_status)
     if(ic_status/=0) then
        if(myid==1) write(*,*) 'ERROR: cannot open SMBH sink seeds ',trim(filename)
        call clean_stop
     endif
     do
        read(ic_unit,*,iostat=ic_status) ic_values
        if(ic_status==iostat_end) exit
        if(ic_status/=0) then
           if(myid==1) write(*,*) 'ERROR: invalid SMBH sink seed row in ',trim(filename)
           call clean_stop
        endif
        if(nsink>=nsinkmax .or. .not.all(ieee_is_finite(ic_values))) then
           if(myid==1) write(*,*) 'ERROR: SMBH sink seed capacity or finite-value failure'
           call clean_stop
        endif
        if(ic_values(1)<=0d0 .or. any(abs(ic_values(2:4))>=boxlen/2d0) .or. &
             & abs(ic_values(11)-ic_values(1))>1d-12*ic_values(1) .or. ic_values(12)/=0d0) then
           if(myid==1) write(*,*) 'ERROR: unsupported SMBH seed mass, position, or dark-envelope mass'
           call clean_stop
        endif
        nsink=nsink+1
        nindsink=nindsink+1
        idsink(nsink)=nindsink
        msink(nsink)=ic_values(1)
        xsink(nsink,1:ndim)=ic_values(2:1+ndim)+boxlen/2d0
        vsink(nsink,1:ndim)=ic_values(5:4+ndim)
        jsink(nsink,1:ndim)=ic_values(8:7+ndim)
        tsink(nsink)=t
        dMsmbh(nsink)=0d0
        dMBH_coarse(nsink)=0d0
        dMEd_coarse(nsink)=0d0
        Esave(nsink)=0d0
        bhspin(nsink,1:ndim)=0d0
        bhspin(nsink,ndim)=1d0
        spinmag(nsink)=0d0
        sink_stat(nsink,levelmin:nlevelmax,:)=0d0
        ! init_tree sees particles, not the replicated sink arrays. Seed one
        ! canonical sink particle on its owning MPI rank before init_tree.
        seed_pos(1,1:ndim)=xsink(nsink,1:ndim)
        call cmp_cpumap(seed_pos,seed_cpu,1)
        if(seed_cpu(1)==myid) then
           if(npart>=npartmax) then
              write(*,*) 'ERROR: particle capacity exhausted by SMBH seed'
              call clean_stop
           endif
           npart=npart+1
           xp(npart,1:ndim)=xsink(nsink,1:ndim)
           vp(npart,1:ndim)=vsink(nsink,1:ndim)
           mp(npart)=msink(nsink)
           idp(npart)=-int(nsink,i8b)
           ptypep(npart)=PTYPE_SINK
           levelp(npart)=levelmin
           tp(npart)=t
           tpp(npart)=0d0
           mp0(npart)=0d0
           indtab(npart)=0d0
#ifdef OUTPUT_PARTICLE_POTENTIAL
           ptcl_phi(npart)=0d0
#endif
        endif
     enddo
     close(ic_unit,iostat=ic_status)
     if(ic_status/=0) then
        if(myid==1) write(*,*) 'ERROR: cannot close SMBH sink seed file'
        call clean_stop
     endif
     if(myid==1) write(*,'(A,I0,A,A)') 'Loaded ',nsink,' SMBH seeds from ',trim(filename)
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
  weighted_density=0d0
  weighted_volume=0d0
  weighted_momentum=0d0
  weighted_c2=0d0
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

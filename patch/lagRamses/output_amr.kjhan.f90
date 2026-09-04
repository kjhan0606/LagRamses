!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine dump_all
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons, only: phi_checkpoint_level_valid
  use cooling_module
  use power_spectrum_mod, only: compute_power_spectrum
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character,dimension(:),allocatable::nml_bytes
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::filename,filedir,filedirini,filecmd
  integer::i,itest,info,ierr,ilevel
  integer(kind=8)::nml_size
  logical::marker_exists,phi_marker_valid,phi_marker_valid_all

  if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
  if(nstep_coarse==0.and.nrestart>0)return
  if(verbose)write(*,*)'Entering dump_all'

  call write_screen
  call title(ifout,nchar)
  ifout=ifout+1
  if(t>=tout(iout).or.aexp>=aout(iout))iout=iout+1
  output_done=.true.

  ! The HDF5 normal-output path stores the AMR/particle/Poisson payload in
  ! data_<output>.h5 and jumps over backup_psi.  A pure-FDM outer ledger is
  ! only consumable when the complete per-rank wave shards are present, so
  ! fail closed instead of advertising an available field that was not
  ! written.  Operators must select the original binary output format until
  ! an equivalent HDF5 FDM-field writer is implemented.
  if(use_fdm .and. fdm_outer_ledger .and. trim(outformat)=='hdf5')then
     if(myid==1)then
        write(*,'(A)') 'ERROR: fdm_outer_ledger requires outformat=original; HDF5 FDM field shards are unavailable'
        call flush(6)
     endif
     call clean_stop
  endif
  
  if(IOGROUPSIZEREP>0)call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)

  if(ndim>1)then
     if(IOGROUPSIZEREP>0) then
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif

     filecmd='mkdir -p '//TRIM(filedir)
     
     if (.not.withoutmkdir) then
#ifdef NOSYSTEM
!jhshin1
        if(IOGROUPSIZEREP>0) then
           call PXFMKDIR(TRIM(filedirini),LEN(TRIM(filedirini)),O'775',info)
        endif
        call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'775',info)
!jhshin2
#else
        call system(filecmd)
#endif
     endif
     
#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     ! Remove markers from a pre-existing directory before any component is
     ! rewritten.  A failed replacement dump must never inherit validity or
     ! completion state from an older snapshot with the same output number.
     if(myid==1)then
        filename='output_'//TRIM(nchar)//'/POISSON_PHI_VALID'
        inquire(file=TRIM(filename),exist=marker_exists)
        if(marker_exists)then
           open(unit=11,file=TRIM(filename),status='old')
           close(11,status='delete')
        end if
        filename='output_'//TRIM(nchar)//'/COMPLETE'
        inquire(file=TRIM(filename),exist=marker_exists)
        if(marker_exists)then
           open(unit=11,file=TRIM(filename),status='old')
           close(11,status='delete')
        end if
     end if
#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'Start backup header'
     ! Output header: must be called by each process !
     filename=TRIM(filedir)//'header_'//TRIM(nchar)//'.txt'
     call output_header(filename)
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif     
     if(myid==1.and.print_when_io) write(*,*)'End backup header'

     if(myid==1.and.print_when_io) write(*,*)'Start backup info etc.'
     ! Only master process
     if(myid==1)then
        filename=TRIM(filedir)//'info_'//TRIM(nchar)//'.txt'
        call output_info(filename)
        filename=TRIM(filedir)//'makefile.txt'
        call output_makefile(filename)
       ! filename=TRIM(filedir)//'patches.txt'
       ! call output_patch(filename)
        if(hydro)then
           filename=TRIM(filedir)//'hydro_file_descriptor.txt'
           call file_descriptor_hydro(filename)
        end if
        if(cooling .and. .not. neq_chem)then
           filename=TRIM(filedir)//'cooling_'//TRIM(nchar)//'.out'
           call output_cool(filename)
        end if
        if(sink)then
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.info'
           call output_sink(filename)
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.csv'
           call output_sink_csv(filename)
        endif
        ! Copy namelist file to output directory
        filename=TRIM(filedir)//'namelist.txt'
        OPEN(UNIT=10, FILE=namelist_file, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', &
             & ACTION='READ', IOSTAT=IERR)
        if(IERR/=0)then
           write(*,*)'Cannot open namelist for provenance copy: ',TRIM(namelist_file)
           call clean_stop
        endif
        INQUIRE(UNIT=10, SIZE=nml_size, IOSTAT=IERR)
        if(IERR/=0 .or. nml_size<=0)then
           write(*,*)'Cannot size namelist for provenance copy: ',TRIM(namelist_file)
           close(10)
           call clean_stop
        endif
        allocate(nml_bytes(nml_size))
        READ(UNIT=10, POS=1, IOSTAT=IERR)nml_bytes
        if(IERR/=0)then
           write(*,*)'Cannot read namelist for provenance copy: ',TRIM(namelist_file)
           deallocate(nml_bytes)
           close(10)
           call clean_stop
        endif
        OPEN(UNIT=11, FILE=filename, ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE', &
             & ACTION='WRITE', IOSTAT=IERR)
        if(IERR/=0)then
           write(*,*)'Cannot create namelist provenance copy: ',TRIM(filename)
           close(10)
           call clean_stop
        endif
        WRITE(UNIT=11, POS=1, IOSTAT=IERR)nml_bytes
        if(IERR/=0)then
           write(*,*)'Cannot write namelist provenance copy: ',TRIM(filename)
           deallocate(nml_bytes)
           close(10)
           close(11)
           call clean_stop
        endif
        deallocate(nml_bytes)
        CLOSE(10)
        CLOSE(11)
        ! Copy compilation details to output directory
        filename=TRIM(filedir)//'compilation.txt'
        OPEN(UNIT=11, FILE=filename, FORM='formatted')
        write(11,'(" compile date = ",A)')TRIM(builddate)
        write(11,'(" patch dir    = ",A)')TRIM(patchdir)
        write(11,'(" remote repo  = ",A)')TRIM(gitrepo)
        write(11,'(" local branch = ",A)')TRIM(gitbranch)
        write(11,'(" last commit  = ",A)')TRIM(githash)
        CLOSE(11)
        filename=TRIM(filedir)//'dm_run_provenance_'//TRIM(nchar)//'.txt'
        call output_dm_run_provenance(filename)
     endif
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End backup info etc.'

#ifdef HDF5
     if(outformat == 'hdf5') then
        call dump_all_hdf5(filedir, nchar)
        goto 998  ! skip binary output
     end if
#endif

     if(myid==1.and.print_when_io) write(*,*)'Start backup amr'
     filename=TRIM(filedir)//'amr_'//TRIM(nchar)//'.out'
     call backup_amr(filename)
#ifndef WITHOUTMPI
     if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1.and.print_when_io) write(*,*)'End backup amr'

     if(hydro)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup hydro'
        filename=TRIM(filedir)//'hydro_'//TRIM(nchar)//'.out'
        call backup_hydro(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup hydro'
     end if
     
#ifdef RT
     if(rt.or.neq_chem)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup rt'
        filename=TRIM(filedir)//'rt_'//TRIM(nchar)//'.out'
        call rt_backup_hydro(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup rt'
     endif
#endif
    
     if(pic)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup part'
        filename=TRIM(filedir)//'part_'//TRIM(nchar)//'.out'
        call backup_part(filename)
        if(sink)then
           filename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.out'
           call backup_sink(filename)
        end if
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup part'
     end if
     
     if(poisson)then
        if(myid==1) write(*,*)'Start backup poisson'
        if(myid==1) call flush(6)
        filename=TRIM(filedir)//'grav_'//TRIM(nchar)//'.out'
        call backup_poisson(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1) write(*,*)'End backup poisson'
        if(myid==1) call flush(6)
     end if

     if(use_fdm)then
        ! Compute the compact diagnostic before backup_psi.  The diagnostic
        ! refreshes same-level ghost cells; backup_psi uses the output-token
        ! point-to-point protocol.  Keeping those communication phases in
        ! this order avoids carrying requests from the ghost refresh into the
        ! token protocol (which otherwise can fail in MPI_Waitall after a
        ! successful shard write).
        if(fdm_outer_ledger)then
           if(myid==1) write(*,*)'Start FDM outer-wave provenance'
           call output_fdm_outer_wave_provenance(nchar)
#ifndef WITHOUTMPI
           if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
           if(myid==1) write(*,*)'End FDM outer-wave provenance'
           if(myid==1) call flush(6)
        end if
        if(myid==1) write(*,*)'Start backup fdm (psi)'
        if(myid==1) call flush(6)
        filename=TRIM(filedir)//'fdm_'//TRIM(nchar)//'.out'
        call backup_psi(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1) write(*,*)'End backup fdm (psi)'
        if(myid==1) call flush(6)
     end if
#ifdef ATON
     if(aton)then
        if(myid==1.and.print_when_io) write(*,*)'Start backup rad'
        filename=TRIM(filedir)//'rad_'//TRIM(nchar)//'.out'
        call backup_radiation(filename)
        filename=TRIM(filedir)//'radgpu_'//TRIM(nchar)//'.out'
        call store_radiation(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup rad'
     end if
#endif
     if (gadget_output) then
        if(myid==1.and.print_when_io) write(*,*)'Start backup gadget format'
        filename=TRIM(filedir)//'gsnapshot_'//TRIM(nchar)
        call savegadget(filename)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End backup gadget format'
     end if
998  continue

     ! Power spectrum measurement at base level
     if(dump_pk .and. poisson) then
        if(myid==1.and.print_when_io) write(*,*)'Start power spectrum'
        call compute_power_spectrum(levelmin, filedir, nchar)
#ifndef WITHOUTMPI
        if(synchro_when_io) call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
        if(myid==1.and.print_when_io) write(*,*)'End power spectrum'
     end if

     ! Advertise phi only when every active AMR level has completed a solve
     ! since its most recent topology change.  The marker is separate from
     ! the gravity payload so legacy and interrupted outputs fail closed.
     phi_marker_valid=poisson .and. allocated(phi_checkpoint_level_valid)
     if(phi_marker_valid)then
        do ilevel=levelmin,nlevelmax
           if(numbtot(1,ilevel)>0)then
              phi_marker_valid=phi_marker_valid .and. &
                   phi_checkpoint_level_valid(ilevel)
           end if
        end do
     end if
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(phi_marker_valid,phi_marker_valid_all,1,MPI_LOGICAL, &
          MPI_LAND,MPI_COMM_WORLD,info)
#else
     phi_marker_valid_all=phi_marker_valid
#endif
     if(myid==1 .and. phi_marker_valid_all)then
        filename='output_'//TRIM(nchar)//'/POISSON_PHI_VALID'
        open(unit=11,file=TRIM(filename),form='formatted',status='replace', &
             iostat=ierr)
        if(ierr==0)then
           write(11,'(A)')'LAGRAMSES_POISSON_PHI_VALID_V1'
           write(11,*)nstep_coarse,nlevelmax,t,aexp
           close(11)
           write(*,*)'Poisson checkpoint marker: valid warm-start phi'
        else
           write(*,*)'WARNING: could not write Poisson phi validity marker'
        end if
     else if(myid==1 .and. poisson)then
        write(*,*)'Poisson checkpoint marker omitted: predictor required on restart'
     end if

     ! Completion marker, written only once every rank has flushed every
     ! component of this dump.  A backup interrupted part way through leaves the
     ! output directory in place with its tail files truncated or missing, and a
     ! restart that trusts the highest numbered directory then reads garbage: on
     ! 2026-08-03 a filesystem stall during the poisson write cost two runs six
     ! hours each, because every retry restarted from the same broken snapshot.
     ! The marker replaces guesswork about completeness with a fact.  It sits
     ! after label 998 so that the HDF5 path is covered as well.
#ifndef WITHOUTMPI
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(myid==1)then
        filename=TRIM(filedir)//'resolved_physics_inventory_'//TRIM(nchar)//'.txt'
        call output_resolved_physics_inventory(filename,nchar,filedir,phi_marker_valid_all)
     endif

     ! The resolved-physics inventory is deliberately written before COMPLETE.
     ! It indexes raw files only; consumers must still require COMPLETE before
     ! reading any listed snapshot or treating its diagnostic as durable.
     if(myid==1)then
        filename='output_'//TRIM(nchar)//'/COMPLETE'
        open(unit=11,file=TRIM(filename),form='formatted')
        write(11,'(A)')TRIM(nchar)
        close(11)
     endif

  end if

end subroutine dump_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_dm_run_provenance(filename)
  ! Record the active DM realization at each normal output.  This is a
  ! run-level sidecar; it deliberately does not alter sink capture/merging.
  use amr_commons
  implicit none

  character(LEN=*)::filename
  character(LEN=16)::dm_model
  character(LEN=512)::iomsg
  integer::ilun,ios
  real(dp)::sidm_pmax_output

  if(use_fdm .and. sidm) then
     write(*,'(A)') 'ERROR: cannot write DM provenance for simultaneous FDM and SIDM'
     call clean_stop
  endif
  if(use_fdm) then
     dm_model='fdm'
  else if(sidm) then
     dm_model='sidm'
  else if(pic) then
     dm_model='cdm'
  else
     dm_model='none'
  endif

  iomsg=''
  open(newunit=ilun,file=TRIM(filename),status='replace',action='write', &
       & form='formatted',iostat=ios,iomsg=iomsg)
  if(ios /= 0) call dm_run_provenance_fatal('open',filename,ios,iomsg)

  write(ilun,'(A)',iostat=ios,iomsg=iomsg) '# dm_run_provenance_v1'
  if(ios /= 0) call dm_run_provenance_fatal('write schema',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'dark_matter_model = ',trim(dm_model)
  if(ios /= 0) call dm_run_provenance_fatal('write model',filename,ios,iomsg)
  write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'pic_enabled = ',pic
  if(ios /= 0) call dm_run_provenance_fatal('write PIC flag',filename,ios,iomsg)
  write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'sidm_enabled = ',sidm
  if(ios /= 0) call dm_run_provenance_fatal('write SIDM flag',filename,ios,iomsg)
  write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'fdm_enabled = ',use_fdm
  if(ios /= 0) call dm_run_provenance_fatal('write FDM flag',filename,ios,iomsg)
  write(ilun,'(A,I0)',iostat=ios,iomsg=iomsg) 'nstep_coarse = ',nstep_coarse
  if(ios /= 0) call dm_run_provenance_fatal('write coarse step',filename,ios,iomsg)
  write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'time_code = ',t
  if(ios /= 0) call dm_run_provenance_fatal('write time',filename,ios,iomsg)
  write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'aexp = ',aexp
  if(ios /= 0) call dm_run_provenance_fatal('write scale factor',filename,ios,iomsg)
  if(len_trim(githash)>0) then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'build_git_hash = ',trim(githash)
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'build_git_hash = unknown'
  endif
  if(ios /= 0) call dm_run_provenance_fatal('write git hash',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'namelist_copy = namelist.txt'
  if(ios /= 0) call dm_run_provenance_fatal('write namelist link',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'compilation_copy = compilation.txt'
  if(ios /= 0) call dm_run_provenance_fatal('write compilation link',filename,ios,iomsg)
  write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'smbh_capture_ledger_enabled = ', &
       & smbh_capture_ledger
  if(ios /= 0) call dm_run_provenance_fatal('write capture-ledger flag',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'smbh_capture_ledger_file = ', &
       & trim(smbh_capture_ledger_file)
  if(ios /= 0) call dm_run_provenance_fatal('write capture-ledger path',filename,ios,iomsg)
  if(rmerge<0d0) then
     write(*,'(A,1X,ES24.16)') 'ERROR: rmerge must be non-negative for DM provenance',rmerge
     call clean_stop
  endif
  write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'smbh_merge_radius_cells = ',rmerge
  if(ios /= 0) call dm_run_provenance_fatal('write sink merge radius',filename,ios,iomsg)
  if(rmerge==0d0) then
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'smbh_compaction_mode = no_finite_radius_rmerge_zero'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'smbh_compaction_mode = enabled'
  endif
  if(ios /= 0) call dm_run_provenance_fatal('write sink compaction mode',filename,ios,iomsg)

  if(len_trim(model_zoom_manifest_sha256)==64 .and. len_trim(model_zoom_case_id)>0 .and. &
       & len_trim(model_zoom_capture_event_sha256)==64 .and. &
       & len_trim(model_zoom_initial_conditions_sha256)==64 .and. &
       & len_trim(model_zoom_baryon_configuration_sha256)==64 .and. &
       & len_trim(model_zoom_sink_initial_conditions_sha256)==64) then
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'model_zoom_execution_identity_status = available'
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom identity status',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_manifest_sha256 = ', &
          & trim(model_zoom_manifest_sha256)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom manifest identity',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_case_id = ', &
          & trim(model_zoom_case_id)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom case identity',filename,ios,iomsg)
     write(ilun,'(A,I0)',iostat=ios,iomsg=iomsg) 'model_zoom_levelmax = ',nlevelmax
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom level',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_capture_event_sha256 = ', &
          & trim(model_zoom_capture_event_sha256)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom capture identity',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_initial_conditions_sha256 = ', &
          & trim(model_zoom_initial_conditions_sha256)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom IC identity',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_baryon_configuration_sha256 = ', &
          & trim(model_zoom_baryon_configuration_sha256)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom baryon identity',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'model_zoom_sink_initial_conditions_sha256 = ', &
          & trim(model_zoom_sink_initial_conditions_sha256)
     if(ios /= 0) call dm_run_provenance_fatal('write model zoom sink identity',filename,ios,iomsg)
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'model_zoom_execution_identity_status = unavailable'
     if(ios /= 0) call dm_run_provenance_fatal('write unavailable model zoom identity',filename,ios,iomsg)
  endif

  select case(trim(dm_model))
  case('cdm')
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'dm_transport = collisionless_nbody'
     if(ios /= 0) call dm_run_provenance_fatal('write CDM transport',filename,ios,iomsg)
     if(len_trim(cdm_zoom_plan_manifest_sha256)==64 .and. &
          & len_trim(cdm_zoom_capture_event_sha256)==64 .and. &
          & len_trim(cdm_zoom_host_orbit_initial_conditions_sha256)==64 .and. &
          & len_trim(cdm_zoom_initial_conditions_sha256)==64 .and. &
          & len_trim(cdm_zoom_sink_initial_conditions_sha256)==64) then
        write(ilun,'(A)',iostat=ios,iomsg=iomsg) &
             & 'cdm_zoom_execution_identity_status = available'
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom identity status',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'cdm_zoom_plan_manifest_sha256 = ', &
             & trim(cdm_zoom_plan_manifest_sha256)
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom plan identity',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'cdm_zoom_capture_event_sha256 = ', &
             & trim(cdm_zoom_capture_event_sha256)
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom capture identity',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) &
             & 'cdm_zoom_host_orbit_initial_conditions_sha256 = ', &
             & trim(cdm_zoom_host_orbit_initial_conditions_sha256)
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom host-orbit identity',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) &
             & 'cdm_zoom_initial_conditions_sha256 = ', &
             & trim(cdm_zoom_initial_conditions_sha256)
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom IC identity',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) &
             & 'cdm_zoom_sink_initial_conditions_sha256 = ', &
             & trim(cdm_zoom_sink_initial_conditions_sha256)
        if(ios /= 0) call dm_run_provenance_fatal('write CDM zoom sink identity',filename,ios,iomsg)
     else
        write(ilun,'(A)',iostat=ios,iomsg=iomsg) &
             & 'cdm_zoom_execution_identity_status = unavailable'
        if(ios /= 0) call dm_run_provenance_fatal('write unavailable CDM zoom identity',filename,ios,iomsg)
     endif
  case('sidm')
     sidm_pmax_output=maxval(sidm_Pmax(1:nlevelmax))
     write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'sidm_cross_section_cm2_g = ', &
          & sidm_cross_section
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM cross section',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'sidm_type = ',trim(sidm_type)
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM type',filename,ios,iomsg)
     write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'sidm_v0_km_s = ',sidm_v0
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM velocity scale',filename,ios,iomsg)
     write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'sidm_power = ',sidm_power
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM power',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'sidm_angular = ',trim(sidm_angular)
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM angular model',filename,ios,iomsg)
     write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'sidm_inelastic = ',sidm_inelastic
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM inelastic flag',filename,ios,iomsg)
     write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'sidm_max_scatter_probability = ', &
          & sidm_pmax_output
     if(ios /= 0) call dm_run_provenance_fatal('write SIDM probability',filename,ios,iomsg)
  case('fdm')
     write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'm_axion_ev = ',m_axion
     if(ios /= 0) call dm_run_provenance_fatal('write FDM axion mass',filename,ios,iomsg)
     write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'fdm_use_hjm = ',fdm_use_hjm
     if(ios /= 0) call dm_run_provenance_fatal('write FDM HJM flag',filename,ios,iomsg)
     write(ilun,'(A,I0)',iostat=ios,iomsg=iomsg) 'fdm_first_wave_level = ',fdm_first_wave_level
     if(ios /= 0) call dm_run_provenance_fatal('write FDM wave level',filename,ios,iomsg)
     write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'fdm_outer_ledger_enabled = ', &
          & fdm_outer_ledger
     if(ios /= 0) call dm_run_provenance_fatal('write FDM ledger flag',filename,ios,iomsg)
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_force_accounting = resolved_wave_only'
     if(ios /= 0) call dm_run_provenance_fatal('write FDM force accounting',filename,ios,iomsg)
  end select

  flush(unit=ilun,iostat=ios,iomsg=iomsg)
  if(ios /= 0) call dm_run_provenance_fatal('flush',filename,ios,iomsg)
  close(ilun,iostat=ios,iomsg=iomsg)
  if(ios /= 0) call dm_run_provenance_fatal('close',filename,ios,iomsg)
end subroutine output_dm_run_provenance
!#########################################################################
subroutine dm_run_provenance_fatal(operation,filename,status,message)
  use amr_commons
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  character(LEN=*),intent(in)::operation,filename,message
  integer,intent(in)::status

  write(error_unit,'(A,1X,A,1X,A,1X,I0,1X,A)') 'DM provenance I/O failure:', &
       & trim(operation),trim(filename),status,trim(message)
  call flush(error_unit)
  call clean_stop
end subroutine dm_run_provenance_fatal
!#########################################################################
subroutine output_resolved_physics_inventory(filename,output_char,output_directory, &
     & phi_checkpoint_valid)
  ! Record what raw output evidence exists for a model-specific postprocessor.
  ! This routine does not measure a profile, decompose a force, or modify any
  ! CDM/SIDM/FDM dynamics.  In particular, a missing force/conservation/SIDM
  ! scatter ledger is emitted as unavailable rather than represented by zero.
  use amr_commons
  implicit none

  character(LEN=*),intent(in)::filename,output_char,output_directory
  logical,intent(in)::phi_checkpoint_valid
  character(LEN=16)::dm_model
  character(LEN=40)::stars_status,gas_status,dm_status
  character(LEN=512)::iomsg
  integer::ilun,ios

  if(use_fdm .and. sidm) then
     write(*,'(A)') 'ERROR: cannot write resolved inventory for simultaneous FDM and SIDM'
     call clean_stop
  endif
  if(use_fdm) then
     dm_model='fdm'
     dm_status='available'
  else if(sidm) then
     dm_model='sidm'
     dm_status='available'
  else if(pic) then
     dm_model='cdm'
     dm_status='available'
  else
     dm_model='none'
     dm_status='absent'
  endif
  if(pic) then
     ! A particle dump alone cannot distinguish stars from collisionless DM.
     ! Do not infer an available stellar force channel without classification.
     stars_status='requires_particle_classification'
  else
     stars_status='absent'
  endif
  if(hydro) then
     gas_status='available'
  else
     gas_status='absent'
  endif

  iomsg=''
  open(newunit=ilun,file=TRIM(filename),status='replace',action='write', &
       & form='formatted',iostat=ios,iomsg=iomsg)
  if(ios/=0) call resolved_physics_inventory_fatal('open',filename,ios,iomsg)

  write(ilun,'(A)',iostat=ios,iomsg=iomsg) '# lagramses_resolved_physics_inventory_v2'
  if(ios/=0) call resolved_physics_inventory_fatal('write schema',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) &
       & '# Raw file availability only; this record is not a force or delay measurement.'
  if(ios/=0) call resolved_physics_inventory_fatal('write scope',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'output_number = ',trim(output_char)
  if(ios/=0) call resolved_physics_inventory_fatal('write output number',filename,ios,iomsg)
  write(ilun,'(A,I0)',iostat=ios,iomsg=iomsg) 'nstep_coarse = ',nstep_coarse
  if(ios/=0) call resolved_physics_inventory_fatal('write coarse step',filename,ios,iomsg)
  write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'time_code = ',t
  if(ios/=0) call resolved_physics_inventory_fatal('write time',filename,ios,iomsg)
  write(ilun,'(A,ES24.16)',iostat=ios,iomsg=iomsg) 'aexp = ',aexp
  if(ios/=0) call resolved_physics_inventory_fatal('write scale factor',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'dark_matter_model = ',trim(dm_model)
  if(ios/=0) call resolved_physics_inventory_fatal('write model',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'raw_snapshot_directory = ', &
       & trim(output_directory)
  if(ios/=0) call resolved_physics_inventory_fatal('write directory',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'completion_marker = COMPLETE'
  if(ios/=0) call resolved_physics_inventory_fatal('write completion marker',filename,ios,iomsg)
  write(ilun,'(A,L1)',iostat=ios,iomsg=iomsg) 'star_formation_enabled = ',star
  if(ios/=0) call resolved_physics_inventory_fatal('write star flag',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'stars_channel_status = ',trim(stars_status)
  if(ios/=0) call resolved_physics_inventory_fatal('write stars status',filename,ios,iomsg)
  if(pic)then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'stars_particle_snapshot_prefix = part_', &
          & trim(output_char)//'.out'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'stars_particle_snapshot_prefix = none'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write stars snapshot',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'gas_channel_status = ',trim(gas_status)
  if(ios/=0) call resolved_physics_inventory_fatal('write gas status',filename,ios,iomsg)
  if(hydro)then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'gas_snapshot_prefix = hydro_', &
          & trim(output_char)//'.out'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'gas_snapshot_prefix = none'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write gas snapshot',filename,ios,iomsg)
  write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'dark_matter_channel_status = ',trim(dm_status)
  if(ios/=0) call resolved_physics_inventory_fatal('write DM status',filename,ios,iomsg)
  if(pic)then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'particle_snapshot_prefix = part_', &
          & trim(output_char)//'.out'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'particle_snapshot_prefix = none'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write particle snapshot',filename,ios,iomsg)
  if(poisson)then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'potential_snapshot_prefix = grav_', &
          & trim(output_char)//'.out'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'potential_snapshot_prefix = none'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write potential snapshot',filename,ios,iomsg)
  if(.not.poisson)then
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'potential_checkpoint_status = absent'
  else if(phi_checkpoint_valid)then
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'potential_checkpoint_status = validated'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'potential_checkpoint_status = unvalidated'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write potential status',filename,ios,iomsg)
  if(sink)then
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'sink_info_file = sink_', &
          & trim(output_char)//'.info'
  else
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'sink_info_file = none'
  endif
  if(ios/=0) call resolved_physics_inventory_fatal('write sink info',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'force_source_ledger_status = unavailable'
  if(ios/=0) call resolved_physics_inventory_fatal('write force ledger status',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'force_source_ledger_reason = no_source_decomposition_in_normal_output'
  if(ios/=0) call resolved_physics_inventory_fatal('write force ledger reason',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'force_source_ledger_path = none'
  if(ios/=0) call resolved_physics_inventory_fatal('write force ledger path',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'force_source_ledger_sha256 = none'
  if(ios/=0) call resolved_physics_inventory_fatal('write force ledger SHA-256',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'conservation_ledger_status = unavailable'
  if(ios/=0) call resolved_physics_inventory_fatal('write conservation ledger status',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'conservation_ledger_reason = no_time_series_in_normal_output'
  if(ios/=0) call resolved_physics_inventory_fatal('write conservation ledger reason',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'conservation_ledger_path = none'
  if(ios/=0) call resolved_physics_inventory_fatal('write conservation ledger path',filename,ios,iomsg)
  write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'conservation_ledger_sha256 = none'
  if(ios/=0) call resolved_physics_inventory_fatal('write conservation ledger SHA-256',filename,ios,iomsg)

  select case(trim(dm_model))
  case('sidm')
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'sidm_scattering_ledger_status = unavailable'
     if(ios/=0) call resolved_physics_inventory_fatal('write SIDM ledger status',filename,ios,iomsg)
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) &
          & 'sidm_scattering_ledger_reason = no_cumulative_scatter_counter_in_normal_output'
     if(ios/=0) call resolved_physics_inventory_fatal('write SIDM ledger reason',filename,ios,iomsg)
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'sidm_scattering_ledger_path = none'
     if(ios/=0) call resolved_physics_inventory_fatal('write SIDM ledger path',filename,ios,iomsg)
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'sidm_scattering_ledger_sha256 = none'
     if(ios/=0) call resolved_physics_inventory_fatal('write SIDM ledger SHA-256',filename,ios,iomsg)
  case('fdm')
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_field_snapshot_status = available'
     if(ios/=0) call resolved_physics_inventory_fatal('write FDM field status',filename,ios,iomsg)
     write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'fdm_field_snapshot_prefix = fdm_', &
          & trim(output_char)//'.out'
     if(ios/=0) call resolved_physics_inventory_fatal('write FDM field prefix',filename,ios,iomsg)
     if(fdm_outer_ledger)then
        write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_wave_provenance_status = available'
        if(ios/=0) call resolved_physics_inventory_fatal('write FDM provenance status',filename,ios,iomsg)
        write(ilun,'(A,A)',iostat=ios,iomsg=iomsg) 'fdm_wave_provenance_path = output_', &
             & trim(output_char)//'/fdm_outer_wave_provenance_'//trim(output_char)//'.txt'
     else
        write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_wave_provenance_status = unavailable'
        if(ios/=0) call resolved_physics_inventory_fatal('write FDM provenance status',filename,ios,iomsg)
        write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_wave_provenance_path = none'
     endif
     if(ios/=0) call resolved_physics_inventory_fatal('write FDM provenance path',filename,ios,iomsg)
     write(ilun,'(A)',iostat=ios,iomsg=iomsg) 'fdm_force_accounting = resolved_wave_only'
     if(ios/=0) call resolved_physics_inventory_fatal('write FDM force accounting',filename,ios,iomsg)
  end select

  flush(unit=ilun,iostat=ios,iomsg=iomsg)
  if(ios/=0) call resolved_physics_inventory_fatal('flush',filename,ios,iomsg)
  close(ilun,iostat=ios,iomsg=iomsg)
  if(ios/=0) call resolved_physics_inventory_fatal('close',filename,ios,iomsg)
end subroutine output_resolved_physics_inventory
!#########################################################################
subroutine resolved_physics_inventory_fatal(operation,filename,status,message)
  use amr_commons
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  character(LEN=*),intent(in)::operation,filename,message
  integer,intent(in)::status

  write(error_unit,'(A,1X,A,1X,A,1X,I0,1X,A)') 'Resolved inventory I/O failure:', &
       & trim(operation),trim(filename),status,trim(message)
  call flush(error_unit)
  call clean_stop
end subroutine resolved_physics_inventory_fatal
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine backup_amr(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  use morton_hash
  use amr_index, only: icell_legacy
#include "amr_index.h"
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun
  integer::ilevel,ibound,ncache,istart,i,igrid,idim,ind
  integer,allocatable,dimension(:)::ind_grid,iig
  real(dp),allocatable,dimension(:)::xdp
  real(sp),allocatable,dimension(:)::xsp
  real(dp),dimension(1:3)::skip_loc
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale
  integer,parameter::tag=1120
  integer::dummy_io,info2

  if(verbose)write(*,*)'Entering backup_amr'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !-----------------------------------
  ! Output amr grid in file
  !-----------------------------------  
  ilun=myid+10
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
  ! Write grid variables
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)nx,ny,nz
  write(ilun)nlevelmax
  write(ilun)ngridmax
  write(ilun)nboundary
  write(ilun)ngrid_current
  write(ilun)boxlen
  ! Write time variables
  write(ilun)noutput,iout,ifout
  write(ilun)tout(1:noutput)
  write(ilun)aout(1:noutput)
  write(ilun)t
  write(ilun)dtold(1:nlevelmax)
  write(ilun)dtnew(1:nlevelmax)
  write(ilun)nstep,nstep_coarse
  write(ilun)const,mass_tot_0,rho_tot
  write(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
  write(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
  write(ilun)mass_sph
  ! Write levels variables
  write(ilun)headl(1:ncpu,1:nlevelmax)
  write(ilun)taill(1:ncpu,1:nlevelmax)
  write(ilun)numbl(1:ncpu,1:nlevelmax)
  write(ilun)numbtot(1:10,1:nlevelmax)
  ! Read boundary linked list
  if(simple_boundary)then
     write(ilun)headb(1:nboundary,1:nlevelmax)
     write(ilun)tailb(1:nboundary,1:nlevelmax)
     write(ilun)numbb(1:nboundary,1:nlevelmax)
  end if
  ! Write free memory
  write(ilun)headf,tailf,numbf,used_mem,used_mem_tot
  ! Write cpu boundaries
  write(ilun)ordering
  if(ordering=='bisection') then
     write(ilun)bisec_wall(1:nbinodes)
     write(ilun)bisec_next(1:nbinodes,1:2)
     write(ilun)bisec_indx(1:nbinodes)
     write(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
     write(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
  else if(ordering=='ksection') then
     write(ilun)nksec_levels
     write(ilun)ksec_kmax
     write(ilun)ksec_nbinodes
     write(ilun)ksec_factor(1:nksec_levels)
     write(ilun)ksec_dir(1:nksec_levels)
     write(ilun)ksec_wall(1:ksec_nbinodes,1:ksec_kmax-1)
     write(ilun)ksec_next(1:ksec_nbinodes,1:ksec_kmax)
     write(ilun)ksec_indx(1:ksec_nbinodes)
     write(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
     write(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
  else
     write(ilun)bound_key(0:ndomain)
  endif

  ! Write coarse level
  write(ilun)son(1:ncoarse)
  write(ilun)flag1(1:ncoarse)
  write(ilun)cpu_map(1:ncoarse)
  ! Write fine levels
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
           allocate(ind_grid(1:ncache),xdp(1:ncache),iig(1:ncache))
           ! Write grid index
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           write(ilun)ind_grid
           ! Write next index
           do i=1,ncache
              iig(i)=next(ind_grid(i))
           end do
           write(ilun)iig
           ! Write prev index
           do i=1,ncache
              iig(i)=prev(ind_grid(i))
           end do
           write(ilun)iig
           ! Write grid center
           do idim=1,ndim
              do i=1,ncache
                 xdp(i)=xg(ind_grid(i),idim)
              end do
              write(ilun)xdp
           end do
           ! Write father index
           do i=1,ncache
              iig(i)=icell_legacy(father(ind_grid(i)))
           end do
           write(ilun)iig
           ! Write nbor index
           do ind=1,twondim
              do i=1,ncache
                 iig(i)=icell_legacy(morton_nbor_cell(ind_grid(i),ilevel,ind))
              end do
              write(ilun)iig
           end do
           ! Write son index
           do ind=1,twotondim
              do i=1,ncache
                 iig(i)=son(ICELL_OF(ind_grid(i),ind))
              end do
              write(ilun)iig
           end do
           ! Write cpu map
           do ind=1,twotondim
              do i=1,ncache
                 iig(i)=cpu_map(ICELL_OF(ind_grid(i),ind))
              end do
              write(ilun)iig
           end do
           ! Write refinement map
           do ind=1,twotondim
              do i=1,ncache
                 iig(i)=flag1(ICELL_OF(ind_grid(i),ind))
              end do
              write(ilun)iig
           end do
           deallocate(xdp,iig,ind_grid)
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

end subroutine backup_amr
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_info(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
#ifdef PHASE0_STELLAR_ENRICHMENT
  use stellar_enrichment_config, only: stellar_feedback_mode, default_imf_id, &
       population_model_id, configured_channel_mass_min, &
       configured_channel_mass_max, n_stellar_channels, active_element, &
       n_stellar_elements, enable_wind, enable_agb, enable_snii, enable_snia, &
       enable_pisn, yield_source_basis_name, configured_imf_mass_min, &
       configured_imf_mass_max, configured_binary_fraction, stellar_fate_policy, &
       stellar_fate_map_sha256, stellar_fate_approval_id
  use stellar_ramses_runtime, only: phase0_get_runtime_identity
#endif
  implicit none
  character(LEN=80)::filename

  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
#ifdef PHASE0_STELLAR_ENRICHMENT
  integer::stellar_channel,stellar_table_rows,stellar_element
  character(LEN=1024)::stellar_table_path
  logical::stellar_table_loaded
#endif
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  character(LEN=80)::fileloc
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_info'

  ilun=11

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  
  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')scale
  write(ilun,'("time        =",E23.15)')t
  write(ilun,'("aexp        =",E23.15)')aexp
  write(ilun,'("H0          =",E23.15)')h0
  write(ilun,'("omega_m     =",E23.15)')omega_m
  write(ilun,'("omega_l     =",E23.15)')omega_l
  write(ilun,'("omega_k     =",E23.15)')omega_k
  write(ilun,'("omega_b     =",E23.15)')omega_b
  write(ilun,'("unit_l      =",E23.15)')scale_l
  write(ilun,'("unit_d      =",E23.15)')scale_d
  write(ilun,'("unit_t      =",E23.15)')scale_t
#ifdef PHASE0_STELLAR_ENRICHMENT
  write(ilun,'("feedback_mode=",A)')trim(stellar_feedback_mode)
  write(ilun,'("stellar_imf_id=",I11)')default_imf_id
  write(ilun,'("stellar_population_model_id=",I11)')population_model_id
  write(ilun,'("stellar_yield_source_basis=",A)')trim(yield_source_basis_name())
  write(ilun,'("stellar_imf_mass_support=",2(1X,E23.15))') &
       configured_imf_mass_min,configured_imf_mass_max
  write(ilun,'("stellar_binary_fraction=",E23.15)')configured_binary_fraction
  write(ilun,'("stellar_terminal_fate_policy=",A)')trim(stellar_fate_policy)
  write(ilun,'("stellar_fate_map_sha256=",A)')trim(stellar_fate_map_sha256)
  write(ilun,'("stellar_fate_approval_id=",A)')trim(stellar_fate_approval_id)
  write(ilun,'("stellar_channel_enabled=",5(1X,L1))')enable_wind,enable_agb, &
       enable_snii,enable_snia,enable_pisn
  write(ilun,'("stellar_active_elements=",11(1X,L1))') &
       (active_element(stellar_element),stellar_element=1,n_stellar_elements)
  do stellar_channel=1,n_stellar_channels
     write(ilun,'("stellar_channel_mass_window=",I2,2(1X,E23.15))') &
          stellar_channel,configured_channel_mass_min(stellar_channel), &
          configured_channel_mass_max(stellar_channel)
  end do
  call phase0_get_runtime_identity(stellar_table_path,stellar_table_rows, &
       stellar_table_loaded)
  write(ilun,'("phase0_yield_table_loaded=",L1)')stellar_table_loaded
  write(ilun,'("phase0_yield_table_rows=",I11)')stellar_table_rows
  if(stellar_table_loaded .and. len_trim(stellar_table_path)>0)then
     write(ilun,'("phase0_yield_table=",A)')trim(stellar_table_path)
  else
     write(ilun,'("phase0_yield_table=<not_loaded>")')
  endif
#endif
  write(ilun,*)
  
  ! Write ordering information
  write(ilun,'("ordering type=",A80)')ordering
  if(ordering=='bisection') then
     do icpu=1,ncpu
        ! write 2*ndim floats for cpu bound box
        write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
        ! write 1 float for cpu load
        write(ilun,'(E23.15)')dble(bisec_cpu_load(icpu))
     end do
  else if(ordering=='ksection') then
     do icpu=1,ncpu
        write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
     end do
  else
     write(ilun,'("   DOMAIN   ind_min                 ind_max")')
     do idom=1,ndomain
        write(ilun,'(I8,1X,E23.15,1X,E23.15)')idom,bound_key(idom-1),bound_key(idom)
     end do
  endif

  close(ilun)

end subroutine output_info
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine output_header(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::info,ilun,ielt
  integer(i8b)::tmp_long,npart_tot,tmp_long2
  character(LEN=80)::fileloc
  character(LEN=150)::header_string

  if(verbose)write(*,*)'Entering output_header'

  ! Compute total number of particles
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(npart,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  tmp_long=npart
  call MPI_ALLREDUCE(tmp_long,npart_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
  tmp_long2=nstar_tot
#endif
#endif
#ifdef WITHOUTMPI
  npart_tot=npart
#endif

  if(myid==1)then

     ilun=myid+10

     ! Open file
     fileloc=TRIM(filename)
     open(unit=ilun,file=fileloc,form='formatted')
     
     ! Write header information
     write(ilun,*)'Total number of particles'
     write(ilun,*)npart_tot
     write(ilun,*)'Total number of dark matter particles'
     write(ilun,*)nDM
     write(ilun,*)'Total number of star particles'
#ifndef LONGINT
     write(ilun,*)nstar_tot
#else
     write(ilun,*)tmp_long2
#endif
     write(ilun,*)'Total number of sink particles'
     write(ilun,*)nsink
     write(ilun,*)'Total number of cloud particles (including sinks)'
#ifndef LONGINT
     write(ilun,*)npart_tot-nstar_tot-nDM
#else
     write(ilun,*)npart_tot-tmp_long2-nDM
#endif
     ! Keep track of what particle fields are present
     write(ilun,*)'Particle fields'
     write(ilun,'(a)',advance='no')'pos vel mass iord level '
#ifdef OUTPUT_PARTICLE_POTENTIAL
     write(ilun,'(a)',advance='no')'phi '
#endif
     if(star.or.sink) then
        write(ilun,'(a)',advance='no')'tform metal propertime initialmass'
     endif
     
     if(sf_birth_properties) then
        write(ilun,*)''
        write(ilun,*)'Star fields'
        header_string = 'ID level birth_epoch birth_proper_time M0 x y z vx vy vz rho T/mu Met '
        if(sf_virial) header_string = trim(header_string)//'sf_eff '
        do ielt=1,nelt
           header_string = trim(header_string)//' '//elem_list(ielt)
        enddo
        write(ilun,'(a)') trim(header_string)
     endif

     close(ilun)

  endif

end subroutine output_header
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine savegadget(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename
  TYPE (gadgetheadertype) :: header
  real,allocatable,dimension(:,:)::pos, vel
  integer(i8b),allocatable,dimension(:)::ids
  integer::i, idim, ipart
  real:: gadgetvfact
  integer::info
  integer(i8b)::npart_tot, npart_loc
  real, parameter:: RHOcrit = 2.7755d11

#ifndef WITHOUTMPI
  npart_loc=npart
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_loc,npart_tot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#else
  npart_tot=npart
#endif

  allocate(pos(ndim, npart), vel(ndim, npart), ids(npart))
  gadgetvfact = 100.0 * boxlen_ini / aexp / SQRT(aexp)

  header%npart = 0
  header%npart(2) = npart
  header%mass = 0
  header%mass(2) = omega_m*RHOcrit*(boxlen_ini)**3/npart_tot/1.d10
  header%time = aexp
  header%redshift = 1.d0/aexp-1.d0
  header%flag_sfr = 0
  header%nparttotal = 0
#ifndef LONGINT
  header%nparttotal(2) = npart_tot
#else
  header%nparttotal(2) = MOD(npart_tot,4294967296_8)
#endif
  header%flag_cooling = 0
  header%numfiles = ncpu
  header%boxsize = boxlen_ini
  header%omega0 = omega_m
  header%omegalambda = omega_l
  header%hubbleparam = h0/100.0
  header%flag_stellarage = 0
  header%flag_metals = 0
  header%totalhighword = 0
#ifndef LONGINT
  header%totalhighword(2) = 0
#else
  header%totalhighword(2) = npart_tot/4294967296_8
#endif
  header%flag_entropy_instead_u = 0
  header%flag_doubleprecision = 0
  header%flag_ic_info = 0
  header%lpt_scalingfactor = 0
  header%unused = ' '

  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           if (ipart .gt. npart) then
                write(*,*) myid, "Ipart=",ipart, "exceeds", npart
                call clean_stop
           endif
           pos(idim, ipart)=xp(i,idim) * boxlen_ini
           vel(idim, ipart)=vp(i,idim) * gadgetvfact
           if (idim.eq.1) ids(ipart) = idp(i)
        end if
     end do
  end do

  call gadgetwritefile(filename, myid-1, header, pos, vel, ids)
  deallocate(pos, vel, ids)

end subroutine savegadget

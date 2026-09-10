program snrt_checkpoint_smoke
  use amr_parameters, only: dp, amr_block_size, ngridmax, twotondim
  use amr_commons, only: active, ncoarse, son
  use iso_c_binding, only: c_float
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use snrt_spectral_contract, only: &
       snrt_spectral_contract_load, snrt_spectral_contract_load_from_environment, &
       snrt_spectral_contract_ok, snrt_spectral_contract_loaded, &
       snrt_spectral_contract_runtime_allowed, snrt_group_mean_energy_ev,snrt_band_enabled,snrt_group_edges_ev,snrt_band_kind
  use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment, &
       snrt_secondary_tables_loaded, snrt_secondary_loaded_manifest_sha256, &
       snrt_thermochemistry_ok
  use snrt_state, only: snrt_nslot, snrt_ndirection, snrt_ngroups, &
       snrt_intensity, snrt_neutral_fraction, snrt_hydrogen_ii, &
       snrt_helium_ii, snrt_helium_iii, snrt_state_sync_level, &
       snrt_state_get_cell, snrt_state_get_slot, &
       snrt_state_checkpoint_write, snrt_state_checkpoint_read, snrt_energy_shift, &
       snrt_checkpoint_version, snrt_checkpoint_cell_width, snrt_checkpoint_number_width, &
       snrt_state_pack_cell, snrt_state_restore_cell, snrt_state_clear_cell, validate_cell_payload
  implicit none

  character(len=1024) :: checkpoint_file, candidate_file
  integer :: ierr, ios, failures, nleaf, nnew, unit
  integer :: idir, igroup, islot
  real(dp) :: expected_value
  logical :: payload_match
  real(c_float) :: bad_photons(3), first_photon
  integer :: case_id, stream_size, intensity_pos, write_pos

  checkpoint_file = ''
  candidate_file = ''
  call get_command_argument(1, checkpoint_file)
  call get_command_argument(2, candidate_file)
  failures = 0
  if (len_trim(checkpoint_file) == 0 .or. len_trim(candidate_file) == 0) then
     write(*,'(a)') 'checkpoint smoke requires checkpoint and candidate paths'
     error stop 2
  end if

  call snrt_spectral_contract_load_from_environment(ierr)
  call expect(ierr == snrt_spectral_contract_ok .and. &
       snrt_spectral_contract_loaded .and. snrt_spectral_contract_runtime_allowed, &
       'reference contract is loaded for checkpoint write', failures)
  call snrt_secondary_tables_load_from_environment(ierr)
  call expect(ierr == snrt_thermochemistry_ok .and. snrt_secondary_tables_loaded, &
       'FS2010 contract is loaded for checkpoint identity binding', failures)

  ! Construct one tiny native leaf layout.  This exercises the real state
  ! allocation and checkpoint payload records without launching RAMSES.
  ncoarse = 0
  ngridmax = 1
  amr_block_size = 1
  allocate(active(1))
  active(1)%ngrid = 1
  allocate(active(1)%igrid(1))
  active(1)%igrid(1) = 1
  allocate(son(twotondim))
  son = 0
  call snrt_state_sync_level(1, nleaf, nnew)
  call expect(nleaf == twotondim .and. nnew == twotondim .and. &
       snrt_nslot == twotondim, &
       'native checkpoint smoke allocates the expected leaf payload', failures)
  call expect(all(snrt_energy_shift==0.0_dp),'new slots have zero energy correction',failures)
  do islot = 1, snrt_nslot
     snrt_hydrogen_ii(islot) = 0.10d0 + 0.01d0 * islot
     snrt_neutral_fraction(islot) = 0.90d0 - 0.01d0 * islot
     snrt_helium_ii(islot) = 0.01d0 * islot
     snrt_helium_iii(islot) = 0.005d0 * islot
     do igroup = 1, size(snrt_intensity,2)
        do idir = 1, size(snrt_intensity,1)
           expected_value = 0.001d0 * idir + 0.01d0 * igroup + 0.1d0 * islot
           snrt_intensity(idir,igroup,islot) = real(expected_value,c_float)
           snrt_energy_shift(idir,igroup,islot)=(-1.0_dp)**idir* &
                0.125_dp*snrt_group_mean_energy_ev(igroup)*real(snrt_intensity(idir,igroup,islot),dp)
        end do
     end do
  end do

  if(snrt_band_enabled())then
     call band_checkpoint()
     if(failures/=0)error stop 1
     write(*,*)'SNRT_BAND_CHECKPOINT_OK'
     stop
  endif
  open(newunit=unit, file=trim(checkpoint_file), status='replace', &
       access='stream', form='unformatted', action='readwrite', iostat=ios)
  call expect(ios == 0, 'checkpoint stream opens for writing', failures)
  if (ios /= 0) error stop 3
  call snrt_state_checkpoint_write(unit, ierr)
  if (ierr /= 0) write(*,'(a,i0)') 'checkpoint_write_ierr=', ierr
  call expect(ierr == 0.and.snrt_checkpoint_version==7, 'version-7 checkpoint write succeeds', failures)
  close(unit)

  ! Unique angular families must never reinterpret another resolution's
  ! directional payload; reject at header before changing populated state.
  first_photon=snrt_intensity(1,1,1)
  open(newunit=unit,status='scratch',access='stream',form='unformatted',action='readwrite')
  write(unit)snrt_checkpoint_version,snrt_ndirection+1,snrt_ngroups,snrt_nslot
  rewind(unit)
  call snrt_state_checkpoint_read(unit,ierr)
  call expect(ierr==2.and.snrt_nslot==twotondim.and.snrt_intensity(1,1,1)==first_photon, &
       'angular dimension mismatch rejects before payload commit',failures)
  close(unit)

  ! First ensure the fixed model cannot consume an evolved band spectrum.
  open(newunit=unit,status='scratch',access='stream',form='unformatted',action='readwrite')
  write(unit)8,snrt_ndirection,snrt_ngroups,snrt_nslot
  rewind(unit)
  call snrt_state_checkpoint_read(unit,ierr)
  call expect(ierr==2.and.snrt_intensity(1,1,1)==first_photon, &
       'band-to-fixed model change rejects before payload commit',failures)
  close(unit)

  ! Loading a candidate changes the declared runtime identity but must not
  ! make the existing checkpoint payload interpretable.
  call snrt_spectral_contract_load(trim(candidate_file), ierr)
  call expect(ierr == snrt_spectral_contract_ok .and. &
       snrt_spectral_contract_loaded .and. &
       .not. snrt_spectral_contract_runtime_allowed, &
       'candidate contract loads but remains inadmissible', failures)
  open(newunit=unit, file=trim(checkpoint_file), status='old', &
       access='stream', form='unformatted', action='read', iostat=ios)
  call expect(ios == 0, 'checkpoint stream opens for identity rejection', failures)
  if (ios /= 0) error stop 4
  call snrt_state_checkpoint_read(unit, ierr)
  call expect(ierr == 4 .and. snrt_nslot == twotondim, &
       'checkpoint identity mismatch is rejected before state mutation', failures)
  close(unit)

  call snrt_spectral_contract_load_from_environment(ierr)
  call expect(ierr == snrt_spectral_contract_ok .and. &
       snrt_spectral_contract_runtime_allowed, &
       'reference contract is restored before checkpoint read', failures)
  snrt_secondary_loaded_manifest_sha256 = 'deliberate-mismatch'
  open(newunit=unit, file=trim(checkpoint_file), status='old', &
       access='stream', form='unformatted', action='read', iostat=ios)
  call expect(ios == 0, 'checkpoint stream opens for secondary identity rejection', failures)
  if (ios /= 0) error stop 5
  call snrt_state_checkpoint_read(unit, ierr)
  call expect(ierr == 5 .and. snrt_nslot == twotondim, &
       'secondary-table identity mismatch is rejected before state mutation', failures)
  close(unit)
  call snrt_secondary_tables_load_from_environment(ierr)
  call expect(ierr == 0 .and. snrt_secondary_tables_loaded, &
       'FS2010 identity is restored before checkpoint round trip', failures)
  snrt_intensity = 0.0_c_float
  snrt_energy_shift = 0.0_dp
  snrt_neutral_fraction = 0.0d0
  snrt_hydrogen_ii = 0.0d0
  snrt_helium_ii = 0.0d0
  snrt_helium_iii = 0.0d0
  snrt_nslot = 0
  open(newunit=unit, file=trim(checkpoint_file), status='old', &
       access='stream', form='unformatted', action='read', iostat=ios)
  call expect(ios == 0, 'checkpoint stream reopens for round trip', failures)
  if (ios /= 0) error stop 6
  call snrt_state_checkpoint_read(unit, ierr)
  call expect(ierr == 0 .and. snrt_nslot == twotondim, &
       'version-7 checkpoint payload round-trips', failures)
  payload_match = .true.
  do islot = 1, snrt_nslot
     if (abs(real(snrt_neutral_fraction(islot),dp) - &
          (0.90d0 - 0.01d0 * islot)) > 1.0d-6) payload_match = .false.
     if (abs(real(snrt_hydrogen_ii(islot),dp) - &
          (0.10d0 + 0.01d0 * islot)) > 1.0d-6) payload_match = .false.
     if (abs(real(snrt_helium_ii(islot),dp) - &
          (0.01d0 * islot)) > 1.0d-6) payload_match = .false.
     if (abs(real(snrt_helium_iii(islot),dp) - &
          (0.005d0 * islot)) > 1.0d-6) payload_match = .false.
     do igroup = 1, size(snrt_intensity,2)
        do idir = 1, size(snrt_intensity,1)
           expected_value = 0.001d0 * idir + 0.01d0 * igroup + 0.1d0 * islot
           if (abs(real(snrt_intensity(idir,igroup,islot),dp) - expected_value) > &
                1.0d-6) payload_match = .false.
           if(snrt_energy_shift(idir,igroup,islot)/=(-1.0_dp)**idir* &
                0.125_dp*snrt_group_mean_energy_ev(igroup)*real(snrt_intensity(idir,igroup,islot),dp)) &
                payload_match=.false.
        end do
     end do
  end do
  call expect(payload_match, &
       'round-tripped intensity, signed shift and H/He fractions are preserved for every entry', failures)
  call expect(snrt_state_get_cell(1) == 1 .and. snrt_state_get_slot(1) == 1, &
       'round-tripped cell-to-slot identity is preserved', failures)
  close(unit)

  ! Corrupt the serialized intensity, leaving valid identities and chemistry.
  ! Its offset is derived from the payload sizes, not a fixed header length.
  bad_photons = [-1.0_c_float, ieee_value(0.0_c_float,ieee_quiet_nan), &
       ieee_value(0.0_c_float,ieee_positive_inf)]
  first_photon = snrt_intensity(1,1,1)
  do case_id = 1, 3
     open(newunit=unit, file=trim(checkpoint_file), status='old', &
          access='stream', form='unformatted', action='readwrite')
     inquire(unit=unit,size=stream_size)
     intensity_pos = stream_size - 4*snrt_nslot*storage_size(0.0_dp)/8 - &
          size(snrt_intensity,1)*size(snrt_intensity,2)*snrt_nslot*storage_size(first_photon)/8 + 1
     write(unit,pos=intensity_pos) bad_photons(case_id)
     rewind(unit)
     call snrt_state_checkpoint_read(unit,ierr)
     call expect(ierr == 10 .and. snrt_intensity(1,1,1) == first_photon .and. &
          snrt_nslot == twotondim, 'invalid photon read rejected without state replacement', failures)
     write(unit,pos=intensity_pos) first_photon
     close(unit)

     snrt_intensity(1,1,1) = bad_photons(case_id)
     open(newunit=unit,status='scratch',access='stream',form='unformatted')
     call snrt_state_checkpoint_write(unit,ierr)
     inquire(unit=unit,pos=write_pos)
     call expect(ierr == 10 .and. write_pos == 1, &
          'invalid photon write rejected before header publication', failures)
     close(unit)
     snrt_intensity(1,1,1) = first_photon
  end do

  call energy_shift_cases()
  if (failures == 0) then
     write(*,'(a)') 'SNRT_CHECKPOINT_OK'
  else
     write(*,'(a,i0)') 'SNRT_CHECKPOINT_FAIL count=', failures
     error stop 1
  end if

contains
  subroutine band_checkpoint()
    use dust_composition_optics, only: d03_band_sha256
    use dust_iron_optics, only: fe_band_sha256
    use snrt_spectral_contract, only: snrt_chimes_bank_sha256,snrt_chimes_molecular_sha256
    real(c_float),allocatable::saved(:,:,:)
    real(dp),allocatable::saved_shift(:,:,:)
    integer::g,version,nd,ng,ns,bytes,offset,j
    character(len=:),allocatable::original,changed
    do g=1,snrt_ngroups
       snrt_energy_shift(:,g,:)=(.5d0*(snrt_group_edges_ev(g)+snrt_group_edges_ev(g+1))- &
            snrt_group_mean_energy_ev(g))*real(snrt_intensity(:,g,:),dp)
    enddo
    saved=snrt_intensity;saved_shift=snrt_energy_shift
    open(newunit=unit,file=trim(checkpoint_file),status='replace',access='stream',form='unformatted',action='readwrite')
    call snrt_state_checkpoint_write(unit,ierr)
    call expect(ierr==0,'band native checkpoint write',failures)
    rewind(unit);read(unit)version,nd,ng,ns
    call expect(version==7+snrt_band_kind().and.nd==snrt_ndirection.and.ng==9, &
         'selected band model native version header',failures)
    rewind(unit);snrt_intensity=0;snrt_energy_shift=0
    call snrt_state_checkpoint_read(unit,ierr)
    call expect(ierr==0.and.all(snrt_intensity==saved).and.all(snrt_energy_shift==saved_shift), &
         'band N/E exact restart',failures)
    close(unit)
    if(snrt_band_kind()>=3)then
       inquire(file=trim(checkpoint_file),size=bytes)
       allocate(character(len=bytes)::original)
       open(newunit=unit,file=trim(checkpoint_file),access='stream',form='unformatted',status='old')
       read(unit)original;close(unit)
       do j=3,snrt_band_kind()
          if(snrt_band_kind()==5.and.j/=5)cycle
          if(snrt_band_kind()==6.and.j==4)cycle
          if(j==3)offset=index(original,d03_band_sha256)
          if(j==4)offset=index(original,fe_band_sha256)
          if(j==5)offset=index(original,snrt_chimes_bank_sha256)
          if(j==6)offset=index(original,snrt_chimes_molecular_sha256)
          call expect(offset>0,'compiled spectral node hash serialized',failures)
          if(offset<=0)cycle
          changed=original;changed(offset:offset)='X'
          open(newunit=unit,file=trim(checkpoint_file),access='stream',form='unformatted', &
               status='replace',action='readwrite')
          write(unit)changed;rewind(unit)
          call snrt_state_checkpoint_read(unit,ierr)
          call expect(ierr==5.and.all(snrt_intensity==saved).and.all(snrt_energy_shift==saved_shift), &
               'altered grain node identity rejects before mutation',failures)
          close(unit)
       enddo
    endif
    open(newunit=unit,status='scratch',access='stream',form='unformatted',action='readwrite')
    write(unit)7,snrt_ndirection,snrt_ngroups,snrt_nslot
    rewind(unit);call snrt_state_checkpoint_read(unit,ierr)
    call expect(ierr==2.and.all(snrt_intensity==saved).and.all(snrt_energy_shift==saved_shift), &
         'fixed-to-band model change rejects before mutation',failures)
    close(unit)
    open(newunit=unit,status='scratch',access='stream',form='unformatted',action='readwrite')
    write(unit)merge(9,8,snrt_band_kind()==1),snrt_ndirection,snrt_ngroups,snrt_nslot
    rewind(unit);call snrt_state_checkpoint_read(unit,ierr)
    call expect(ierr==2.and.all(snrt_intensity==saved).and.all(snrt_energy_shift==saved_shift), &
         'mean/node secondary model change rejects before mutation',failures)
    close(unit)
    snrt_energy_shift(1,5,1)=(13.5d0-snrt_group_mean_energy_ev(5))*real(snrt_intensity(1,5,1),dp)
    open(newunit=unit,status='scratch',access='stream',form='unformatted',action='write')
    call snrt_state_checkpoint_write(unit,ierr)
    call expect(ierr/=0,'out-of-band positive energy rejects before checkpoint write',failures)
    close(unit)
    snrt_energy_shift=saved_shift
  end subroutine

  subroutine energy_shift_cases()
    real(dp) :: payload(snrt_checkpoint_cell_width), saved(snrt_checkpoint_cell_width)
    real(dp) :: bad(3), old_shift, actual, target
    integer :: first_shift, u, v, nbytes, shift_pos, shift_bytes, marker, c
    character(len=:), allocatable :: bytes
    first_shift=5+snrt_checkpoint_number_width
    call snrt_state_pack_cell(1,saved,ierr)
    call expect(ierr==0.and.saved(first_shift)<0,'signed cell pack validates actual energy',failures)
    call snrt_state_clear_cell(1)
    call expect(all(snrt_energy_shift(:,:,1)==0).and.all(snrt_intensity(:,:,1)==0), &
         'clear retires photons and correction together',failures)
    call snrt_state_restore_cell(1,saved,ierr)
    call snrt_state_pack_cell(1,payload,ierr)
    call expect(ierr==0.and.all(payload==saved),'cell correction round trip is exact',failures)
    payload=saved
    payload(5)=0.0_dp
    payload(first_shift)=1.0_dp
    call snrt_state_restore_cell(1,payload,ierr)
    call expect(ierr/=0.and.snrt_energy_shift(1,1,1)==saved(first_shift), &
         'zero photons with nonzero energy rejected before restore',failures)
    payload=saved
    payload(first_shift)=-snrt_group_mean_energy_ev(1)*payload(5)-1.0_dp
    call validate_cell_payload(payload,ierr)
    call expect(ierr/=0,'negative actual energy rejected',failures)
    payload=saved
    payload(5)=1.0_dp+2.0_dp**(-25)
    target=0.5_dp*snrt_group_mean_energy_ev(1)*payload(5)
    payload(first_shift)=target-snrt_group_mean_energy_ev(1)*payload(5)
    call snrt_state_restore_cell(1,payload,ierr)
    actual=snrt_group_mean_energy_ev(1)*real(snrt_intensity(1,1,1),dp)+snrt_energy_shift(1,1,1)
    call expect(ierr==0.and.abs(actual-target)<1.0d-14*target, &
         'FP32 count rounding rebases nonzero correction preserving energy',failures)
    payload(first_shift)=-snrt_group_mean_energy_ev(1)*payload(5)
    call snrt_state_restore_cell(1,payload,ierr)
    actual=snrt_group_mean_energy_ev(1)*real(snrt_intensity(1,1,1),dp)+snrt_energy_shift(1,1,1)
    call expect(ierr==0.and.actual==0.0_dp,'zero actual energy survives FP32 rounding',failures)
    payload(first_shift:)=0.0_dp
    call snrt_state_restore_cell(1,payload,ierr)
    call expect(ierr==0.and.all(snrt_energy_shift(:,:,1)==0),'number-only fixtures retain zero shift',failures)
    call snrt_state_restore_cell(1,saved,ierr)

    open(newunit=u,file=trim(checkpoint_file),status='old',access='stream',form='unformatted')
    inquire(unit=u,size=nbytes)
    allocate(character(len=nbytes)::bytes)
    read(u,pos=1)bytes
    close(u)
    shift_bytes=snrt_checkpoint_number_width*snrt_nslot*storage_size(0.0_dp)/8
    shift_pos=nbytes-shift_bytes-snrt_checkpoint_number_width*snrt_nslot*storage_size(0.0_c_float)/8- &
         4*snrt_nslot*storage_size(0.0_dp)/8+1
    bad=[-huge(0.0_dp),ieee_value(0.0_dp,ieee_quiet_nan),ieee_value(0.0_dp,ieee_positive_inf)]
    old_shift=snrt_energy_shift(1,1,1)
    do c=1,size(bad)
       open(newunit=u,status='scratch',access='stream',form='unformatted')
       write(u)bytes
       write(u,pos=shift_pos)bad(c)
       rewind(u)
       call snrt_state_checkpoint_read(u,ierr)
       call expect(ierr==10.and.snrt_energy_shift(1,1,1)==old_shift, &
            'invalid checkpoint shift read leaves live state intact',failures)
       close(u)
       snrt_energy_shift(1,1,1)=bad(c)
       open(newunit=u,status='scratch',access='stream',form='unformatted')
       call snrt_state_checkpoint_write(u,ierr)
       inquire(unit=u,pos=marker)
       call expect(ierr==10.and.marker==1,'invalid checkpoint shift rejected before header',failures)
       close(u)
       snrt_energy_shift(1,1,1)=old_shift
    end do
    ! Construct the exact v6 layout: same identities, cell IDs, FP32 photons
    ! and chemistry, with neither a correction record nor inferred dust work.
    open(newunit=u,status='scratch',access='stream',form='unformatted')
    write(u)bytes(:shift_pos-1),bytes(shift_pos+shift_bytes:)
    write(u)123456
    write(u,pos=1)6
    rewind(u)
    call snrt_state_checkpoint_read(u,ierr)
    read(u,iostat=ios)marker
    call expect(ierr==0.and.ios==0.and.marker==123456.and.all(snrt_energy_shift==0), &
         'v6 restart defaults to zero shift and consumes exactly its payload',failures)
    ! A truncated v7 must never be reinterpreted as the old number-only format.
    write(u,pos=1)7
    rewind(u)
    call snrt_state_checkpoint_read(u,ierr)
    call expect(ierr/=0.and.all(snrt_energy_shift==0),'truncated v7 cannot fall back to v6',failures)
    close(u)
    ! Sequential raw checkpoints use the same record ordering as streams.
    call snrt_state_restore_cell(1,saved,ierr)
    open(newunit=v,status='scratch',form='unformatted')
    call snrt_state_checkpoint_write(v,ierr)
    call expect(ierr==0,'sequential v7 write',failures)
    snrt_energy_shift=0.0_dp
    rewind(v)
    call snrt_state_checkpoint_read(v,ierr)
    call expect(ierr==0.and.snrt_energy_shift(1,1,1)==saved(first_shift),'sequential v7 read',failures)
    close(v)
    ! Cross the initial 1024-slot capacity without evolving RAMSES.
    ngridmax=129
    deallocate(active(1)%igrid,son)
    active(1)%ngrid=ngridmax
    allocate(active(1)%igrid(ngridmax),son(twotondim*ngridmax))
    active(1)%igrid=[(c,c=1,ngridmax)]
    son=0
    call snrt_state_sync_level(1,nleaf,nnew)
    call expect(snrt_nslot>1024.and.snrt_energy_shift(1,1,1)==saved(first_shift).and. &
         all(snrt_energy_shift(:,:,twotondim+1:)==0), &
         'growth preserves correction and zero initializes new capacity',failures)
  end subroutine energy_shift_cases

  subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (condition) then
       write(*,'(a)') 'PASS: ' // trim(label)
    else
       failures = failures + 1
       write(*,'(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine expect

end program snrt_checkpoint_smoke

program snrt_checkpoint_smoke
  use amr_parameters, only: dp, amr_block_size, ngridmax, twotondim
  use amr_commons, only: active, ncoarse, son
  use iso_c_binding, only: c_float
  use snrt_spectral_contract, only: &
       snrt_spectral_contract_load, snrt_spectral_contract_load_from_environment, &
       snrt_spectral_contract_ok, snrt_spectral_contract_loaded, &
       snrt_spectral_contract_runtime_allowed
  use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment, &
       snrt_secondary_tables_loaded, snrt_secondary_loaded_manifest_sha256, &
       snrt_thermochemistry_ok
  use snrt_state, only: snrt_nslot, &
       snrt_intensity, snrt_neutral_fraction, snrt_hydrogen_ii, &
       snrt_helium_ii, snrt_helium_iii, snrt_state_sync_level, &
       snrt_state_get_cell, snrt_state_get_slot, &
       snrt_state_checkpoint_write, snrt_state_checkpoint_read
  implicit none

  character(len=1024) :: checkpoint_file, candidate_file
  integer :: ierr, ios, failures, nleaf, nnew, unit
  integer :: idir, igroup, islot
  real(dp) :: expected_value
  logical :: payload_match

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
  do islot = 1, snrt_nslot
     snrt_hydrogen_ii(islot) = 0.10d0 + 0.01d0 * islot
     snrt_neutral_fraction(islot) = 0.90d0 - 0.01d0 * islot
     snrt_helium_ii(islot) = 0.01d0 * islot
     snrt_helium_iii(islot) = 0.005d0 * islot
     do igroup = 1, size(snrt_intensity,2)
        do idir = 1, size(snrt_intensity,1)
           expected_value = 0.001d0 * idir + 0.01d0 * igroup + 0.1d0 * islot
           snrt_intensity(idir,igroup,islot) = real(expected_value,c_float)
        end do
     end do
  end do

  open(newunit=unit, file=trim(checkpoint_file), status='replace', &
       access='stream', form='unformatted', action='readwrite', iostat=ios)
  call expect(ios == 0, 'checkpoint stream opens for writing', failures)
  if (ios /= 0) error stop 3
  call snrt_state_checkpoint_write(unit, ierr)
  if (ierr /= 0) write(*,'(a,i0)') 'checkpoint_write_ierr=', ierr
  call expect(ierr == 0, 'version-6 checkpoint write succeeds', failures)
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
       'version-6 checkpoint payload round-trips', failures)
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
        end do
     end do
  end do
  call expect(payload_match, &
       'round-tripped intensity and H/He fractions are preserved for every entry', failures)
  call expect(snrt_state_get_cell(1) == 1 .and. snrt_state_get_slot(1) == 1, &
       'round-tripped cell-to-slot identity is preserved', failures)
  close(unit)

  if (failures == 0) then
     write(*,'(a)') 'SNRT_CHECKPOINT_OK'
  else
     write(*,'(a,i0)') 'SNRT_CHECKPOINT_FAIL count=', failures
     error stop 1
  end if

contains

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

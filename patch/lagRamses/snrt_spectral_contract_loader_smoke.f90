program snrt_spectral_contract_loader_smoke
  use snrt_spectral_contract, only: &
       snrt_spectral_contract_load, snrt_spectral_contract_load_from_environment, &
       snrt_spectral_contract_loaded, snrt_spectral_contract_runtime_allowed, &
       snrt_spectral_contract_ok, snrt_spectral_contract_err_missing, &
       snrt_spectral_contract_err_open, snrt_spectral_contract_err_read, &
       snrt_spectral_contract_err_version, snrt_spectral_contract_err_identity, &
       snrt_spectral_contract_err_edges, &
       snrt_spectral_contract_err_fraction_semantics
  implicit none

  character(len=32) :: mode
  character(len=1024) :: path
  integer :: ierr, failures

  mode = ''
  path = ''
  call get_command_argument(1, mode)
  call get_command_argument(2, path)
  failures = 0

  select case (trim(mode))
  case ('unset_env')
     call snrt_spectral_contract_load_from_environment(ierr)
     call expect(ierr == snrt_spectral_contract_err_missing, &
          'unset SNRT_GROUP_CONTRACT is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded .and. &
          .not. snrt_spectral_contract_runtime_allowed, &
          'unset environment leaves the contract unloaded', failures)
  case ('missing_file')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_open, &
          'missing contract file is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'missing file leaves the contract unloaded', failures)
  case ('malformed_namelist')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_read, &
          'malformed namelist is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'malformed namelist leaves the contract unloaded', failures)
  case ('wrong_version')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_version, &
          'unsupported contract version is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'unsupported version leaves the contract unloaded', failures)
  case ('bad_identity')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_identity, &
          'malformed source identity is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'malformed identity leaves the contract unloaded', failures)
  case ('bad_edges')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_edges, &
          'edge digest mismatch at load is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'edge mismatch leaves the contract unloaded', failures)
  case ('bad_fraction_semantics')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_err_fraction_semantics, &
          'unknown fraction semantics is rejected', failures)
     call expect(.not. snrt_spectral_contract_loaded, &
          'unknown fraction semantics leaves the contract unloaded', failures)
  case ('candidate')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_ok .and. &
          snrt_spectral_contract_loaded, &
          'candidate contract can be loaded for inspection', failures)
     call expect(.not. snrt_spectral_contract_runtime_allowed, &
          'candidate contract is not runtime-admissible', failures)
  case ('intrinsic')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_ok .and. &
          snrt_spectral_contract_loaded, &
          'intrinsic fraction contract loads for inspection', failures)
     call expect(.not. snrt_spectral_contract_runtime_allowed, &
          'intrinsic fraction contract is blocked at resolved-domain runtime', failures)
  case ('reference_no_opt_in')
     call snrt_spectral_contract_load(trim(path), ierr)
     call expect(ierr == snrt_spectral_contract_ok .and. &
          snrt_spectral_contract_loaded, &
          'reference contract loads without execution opt-in', failures)
     call expect(.not. snrt_spectral_contract_runtime_allowed, &
          'reference contract requires explicit execution opt-in', failures)
  case default
     write(*,'(a,a)') 'unknown loader smoke mode: ', trim(mode)
     error stop 2
  end select

  if (failures == 0) then
     write(*,'(a,a)') 'SNRT_SPECTRAL_CONTRACT_LOADER_OK mode=', trim(mode)
  else
     write(*,'(a,i0)') 'SNRT_SPECTRAL_CONTRACT_LOADER_FAIL count=', failures
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

end program snrt_spectral_contract_loader_smoke

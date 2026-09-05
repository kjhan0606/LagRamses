program snrt_thermochemistry_loader_smoke
  use snrt_thermochemistry, only: &
       snrt_secondary_tables_load_from_environment, snrt_secondary_tables_load, &
       snrt_secondary_tables_loaded, snrt_thermochemistry_ok, &
       snrt_thermochemistry_err_missing, snrt_thermochemistry_err_open, &
       snrt_thermochemistry_err_read, snrt_thermochemistry_err_identity
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
     call snrt_secondary_tables_load_from_environment(ierr)
     call expect(ierr == snrt_thermochemistry_err_missing .and. &
          .not. snrt_secondary_tables_loaded, &
          'unset secondary contract is rejected', failures)
  case ('missing_file')
     call snrt_secondary_tables_load(trim(path), ierr)
     call expect(ierr == snrt_thermochemistry_err_open .and. &
          .not. snrt_secondary_tables_loaded, &
          'missing secondary contract file is rejected', failures)
  case ('malformed')
     call snrt_secondary_tables_load(trim(path), ierr)
     call expect(ierr == snrt_thermochemistry_err_read .and. &
          .not. snrt_secondary_tables_loaded, &
          'malformed secondary contract is rejected', failures)
  case ('bad_identity')
     call snrt_secondary_tables_load(trim(path), ierr)
     call expect(ierr == snrt_thermochemistry_err_identity .and. &
          .not. snrt_secondary_tables_loaded, &
          'wrong secondary source identity is rejected', failures)
  case default
     write(*,'(a,a)') 'unknown thermochemistry loader mode: ', trim(mode)
     error stop 2
  end select

  if (failures == 0) then
     write(*,'(a,a)') 'SNRT_NATIVE_THERMOCHEMISTRY_LOADER_OK mode=', trim(mode)
  else
     write(*,'(a,i0)') 'SNRT_NATIVE_THERMOCHEMISTRY_LOADER_FAIL count=', failures
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

end program snrt_thermochemistry_loader_smoke

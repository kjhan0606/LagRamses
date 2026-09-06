program snrt_dust_contract_smoke
  use snrt_dust_contract
  implicit none

  integer :: ierr
  character(len=2048) :: valid_path, invalid_path
  character(len=32) :: error_name

  call get_command_argument(1, valid_path)
  call get_command_argument(2, invalid_path)
  if (len_trim(valid_path) == 0 .or. len_trim(invalid_path) == 0) error stop 1

  call snrt_dust_contract_load(trim(valid_path), ierr)
  if (ierr /= snrt_dust_contract_ok .or. .not. snrt_dust_contract_loaded .or. &
       snrt_dust_contract_runtime_allowed .or. &
       snrt_dust_contract_number_groups /= 3 .or. &
       snrt_dust_contract_number_temperature /= 4 .or. &
       abs(snrt_dust_contract_group_edges_ev(4) - 1000.0d0) > 1.0d-12 .or. &
       abs(snrt_dust_contract_absorption_per_h_cm2(2) - 2.0d-21) > 1.0d-32 .or. &
       abs(snrt_dust_contract_temperature_k(4) - 80.0d0) > 1.0d-12 .or. &
       len_trim(snrt_dust_contract_source_id) == 0) error stop 2
  write(*,'(a,i0,a,l1)') 'SNRT_DUST_CONTRACT_CANDIDATE_OK groups=', &
       snrt_dust_contract_number_groups, ' runtime_allowed=', &
       snrt_dust_contract_runtime_allowed

  call snrt_dust_contract_load_from_environment(ierr)
  if (ierr /= snrt_dust_contract_ok .or. .not. snrt_dust_contract_loaded) error stop 3
  write(*,'(a)') 'SNRT_DUST_CONTRACT_ENVIRONMENT_OK'

  call snrt_dust_contract_load(trim(invalid_path), ierr)
  error_name = snrt_dust_contract_error_name(ierr)
  if (ierr /= snrt_dust_contract_err_status .or. snrt_dust_contract_loaded .or. &
       snrt_dust_contract_runtime_allowed .or. snrt_dust_contract_number_groups /= 0 .or. &
       trim(error_name) /= 'status') error stop 4
  write(*,'(a,a)') 'SNRT_DUST_CONTRACT_INVALID_RESET_OK error=', trim(error_name)

  write(*,'(a)') 'SNRT_NATIVE_DUST_CONTRACT_ADMISSION_OK candidate=1 environment=1 reset=1 runtime_gate=1'
end program snrt_dust_contract_smoke

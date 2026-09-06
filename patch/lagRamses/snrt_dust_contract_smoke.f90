program snrt_dust_contract_smoke
  use snrt_dust_contract
  use snrt_dust_ir
  implicit none

  integer :: ierr, ng, nt
  type(dust_ir_table) :: ir_table
  type(dust_ir_diagnostics) :: diagnostic
  real(dust_dp) :: field(2,2,1), photons(2,1), temperature(1), energy(1), capacity(1), initial
  real(dust_dp) :: rays(3,2)
  integer :: links(6,1)
  character(len=2048) :: valid_path, invalid_path, reference_path
  character(len=16) :: expected_reference
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

  call get_command_argument(3, reference_path)
  call get_command_argument(4, expected_reference)
  if (len_trim(reference_path) > 0) then
     call snrt_dust_contract_load(trim(reference_path), ierr)
     if (ierr /= 0 .or. .not. snrt_dust_contract_reference_control .or. &
          snrt_dust_contract_number_groups /= 9 .or. snrt_dust_contract_version < 2) error stop 5
     if (snrt_dust_contract_runtime_allowed .neqv. (trim(expected_reference) == '1')) error stop 6
     if (len_trim(snrt_dust_contract_approval_id) /= 0) error stop 7
     if (snrt_dust_contract_version == 3) then
        ng=snrt_dust_contract_number_ir; nt=snrt_dust_contract_number_temperature
        if (ng /= 2) error stop 10
        call snrt_dust_ir_initialize(ir_table,snrt_dust_contract_ir_energy_ev(1:ng), &
             snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
             snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr)
        if (ierr /= dust_ok) error stop 11
        ! Isolated operator check is allowed without runtime opt-in; it never
        ! starts RAMSES. All spectral and material inputs come from the file.
        field=0; photons=0; temperature=20
        capacity=snrt_dust_contract_heat_capacity_per_h_erg_k
        energy=capacity*temperature; initial=sum(energy)*1d36
        rays(:,1)=[1d0,0d0,0d0]; rays(:,2)=[-1d0,0d0,0d0]; links=0
        call snrt_dust_ir_advance(ir_table,rays,[.5d0,.5d0],links,1d12,1d6,1d5,[1d0],[0d0], &
             field,temperature,photons,diagnostic,ierr,1d-10,128,energy,capacity)
        if (ierr /= dust_ok .or. sum(photons)<=0 .or. temperature(1)>=20) error stop 12
        if (abs((sum(energy)+sum(field)*.5d0)*1d36+diagnostic%escaped_erg-initial) &
             /initial>1d-10) error stop 13
        write(*,'(a)') 'SNRT_DUST_IR_CONTRACT_COOLING_PASS'
     end if
     call snrt_dust_contract_load(trim(valid_path), ierr)
     if (ierr /= 0 .or. snrt_dust_contract_reference_control .or. &
          snrt_dust_contract_runtime_allowed) error stop 8
     if (snrt_dust_contract_number_ir /= 0 .or. any(snrt_dust_contract_ir_energy_ev /= 0)) error stop 14
     call snrt_dust_contract_load(trim(invalid_path), ierr)
     if (snrt_dust_contract_reference_control .or. snrt_dust_contract_runtime_allowed) error stop 9
     write(*,'(a,a)') 'SNRT_DUST_REFERENCE_OPT_IN_PASS expected=', trim(expected_reference)
  end if

  write(*,'(a)') 'SNRT_NATIVE_DUST_CONTRACT_ADMISSION_OK candidate=1 environment=1 reset=1 runtime_gate=1'
end program snrt_dust_contract_smoke

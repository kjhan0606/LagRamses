program snrt_dust_contract_smoke
  use snrt_dust_contract
  use snrt_dust_ir
  implicit none

  integer :: ierr, ng, nt
  type(dust_ir_table) :: ir_table
  type(dust_ir_diagnostics) :: diagnostic
  real(dust_dp) :: field(2,2,1), photons(2,1), temperature(1), energy(1), capacity(1), initial
  real(dust_dp) :: rays(3,2)
  real(dust_dp) :: before_field(2,2,1), before_energy(1), before_photons(2,1), before_temperature(1), recovered
  integer :: links(6,1)
  character(len=2048) :: valid_path, invalid_path, reference_path
  character(len=2048) :: scattering_path, exchange_path
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
     if (snrt_dust_contract_version >= 3) then
        ng=snrt_dust_contract_number_ir; nt=snrt_dust_contract_number_temperature
        if (ng /= 2) error stop 10
        if(snrt_dust_contract_version==4)then
           call snrt_dust_ir_initialize(ir_table,snrt_dust_contract_ir_energy_ev(1:ng), &
                snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
                snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr, &
                snrt_dust_contract_internal_energy_per_h_erg(1:nt))
        else
        call snrt_dust_ir_initialize(ir_table,snrt_dust_contract_ir_energy_ev(1:ng), &
             snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
             snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr)
        endif
        if (ierr /= dust_ok) error stop 11
        ! Isolated operator check is allowed without runtime opt-in; it never
        ! starts RAMSES. All spectral and material inputs come from the file.
        field=0; photons=0; temperature=20
        capacity=snrt_dust_contract_heat_capacity_per_h_erg_k
        energy=capacity*temperature
        if(snrt_dust_contract_version==4)then
           capacity=1d0 ! intentionally unrelated to physical U(T)
           energy=snrt_dust_contract_internal_energy_per_h_erg(2)
        endif
        initial=sum(energy)*1d36
        rays(:,1)=[1d0,0d0,0d0]; rays(:,2)=[-1d0,0d0,0d0]; links=0
        call snrt_dust_ir_advance(ir_table,rays,[.5d0,.5d0],links,1d12,1d6,1d5,[1d0],[0d0], &
             field,temperature,photons,diagnostic,ierr,1d-10,128,energy,capacity)
        if (ierr /= dust_ok .or. sum(photons)<=0 .or. temperature(1)>=20) error stop 12
        if (abs((sum(energy)+sum(field)*.5d0)*1d36+diagnostic%escaped_erg-initial) &
             /initial>1d-10) error stop 13
        write(*,'(a)') 'SNRT_DUST_IR_CONTRACT_COOLING_PASS'
        if(snrt_dust_contract_version==4)then
           call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
                snrt_dust_contract_internal_energy_per_h_erg(1:nt),energy(1),recovered,ierr)
           if(ierr/=dust_ok.or.abs(recovered-temperature(1))>1d-12)error stop 15
           ! Known midpoint of U(log T), independent of the time integrator.
           call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
                snrt_dust_contract_internal_energy_per_h_erg(1:nt),4.5d-22,recovered,ierr)
           if(ierr/=dust_ok.or.abs(recovered-sqrt(200d0))>1d-12)error stop 16
           before_field=field; before_energy=energy; before_photons=photons; before_temperature=temperature
           call snrt_dust_ir_advance(ir_table,rays,[.5d0,.5d0],links,1d12,1d6,1d5,[1d0],[1d0], &
                field,temperature,photons,diagnostic,ierr,1d-10,128,energy,capacity)
           if(ierr/=dust_err_range.or.any(field/=before_field).or.any(energy/=before_energy).or. &
                any(photons/=before_photons).or.any(temperature/=before_temperature))error stop 17
           field=0; photons=0; temperature=10; energy=snrt_dust_contract_internal_energy_per_h_erg(1)
           initial=(sum(energy)+1d6*1d-22)*1d36
           call snrt_dust_ir_advance(ir_table,rays,[.5d0,.5d0],links,1d12,1d6,1d5,[1d0],[1d-22], &
                field,temperature,photons,diagnostic,ierr,1d-10,128,energy,capacity)
           if(ierr/=dust_ok.or.temperature(1)<=10.or.temperature(1)>=20)error stop 20
           if(abs((sum(energy)+sum(field)*.5d0)*1d36+diagnostic%escaped_erg-initial)/initial>1d-10)error stop 21
           ! A duplicate/non-increasing material knot must not create a ready table.
           call snrt_dust_ir_initialize(ir_table,[.001d0,.01d0],[.001d0,.01d0],[1d-21,1d-21], &
                [10d0,20d0],10d0,ierr,[1d-22,1d-22])
           if(ierr/=dust_err_table)error stop 18
           write(*,'(a)') 'SNRT_DUST_U_TABLE_PASS inverse=1 cooling=1 heating=1 closure=1 range_rollback=1 invalid=1'
        endif
     end if
     call snrt_dust_contract_load(trim(valid_path), ierr)
     if (ierr /= 0 .or. snrt_dust_contract_reference_control .or. &
          snrt_dust_contract_runtime_allowed) error stop 8
     if (snrt_dust_contract_number_ir /= 0 .or. any(snrt_dust_contract_ir_energy_ev /= 0)) error stop 14
     if (any(snrt_dust_contract_internal_energy_per_h_erg/=0).or. &
          len_trim(snrt_dust_contract_material_sha256)/=0)error stop 19
     call snrt_dust_contract_load(trim(invalid_path), ierr)
     if (snrt_dust_contract_reference_control .or. snrt_dust_contract_runtime_allowed) error stop 9
     write(*,'(a,a)') 'SNRT_DUST_REFERENCE_OPT_IN_PASS expected=', trim(expected_reference)
  end if

  call get_command_argument(5, scattering_path)
  if(len_trim(scattering_path)>0)then
     call snrt_dust_contract_load(trim(scattering_path),ierr)
     if(ierr/=0.or..not.snrt_dust_contract_scattering_enabled.or. &
          .not.any(snrt_dust_contract_scattering_per_h_cm2>0))error stop 22
     if(snrt_dust_contract_runtime_allowed.neqv.(trim(expected_reference)=='1'))error stop 23
     call snrt_dust_contract_load(trim(valid_path),ierr)
     if(ierr/=0.or.snrt_dust_contract_scattering_enabled.or. &
          any(snrt_dust_contract_scattering_per_h_cm2/=0))error stop 24
     write(*,'(a)')'SNRT_SCATTER_CONTRACT_OPT_IN_AND_RESET_PASS'
  endif
  call get_command_argument(6,exchange_path)
  if(len_trim(exchange_path)>0)then
     call snrt_dust_contract_load(trim(exchange_path),ierr)
     if(ierr/=0.or..not.snrt_dust_contract_exchange_enabled.or. &
          snrt_dust_contract_collision_area_per_h<=0.or.snrt_dust_contract_accommodation/=.5d0)error stop 25
     if(snrt_dust_contract_runtime_allowed.neqv.(trim(expected_reference)=='1'))error stop 26
     call snrt_dust_contract_load(trim(valid_path),ierr)
     if(ierr/=0.or.snrt_dust_contract_exchange_enabled.or. &
          snrt_dust_contract_collision_area_per_h/=0.or.snrt_dust_contract_accommodation/=0)error stop 27
     write(*,'(a)')'SNRT_GAS_DUST_EXCHANGE_CONTRACT_OPT_IN_AND_RESET_PASS'
  endif
  write(*,'(a)') 'SNRT_NATIVE_DUST_CONTRACT_ADMISSION_OK candidate=1 environment=1 reset=1 runtime_gate=1'
end program snrt_dust_contract_smoke

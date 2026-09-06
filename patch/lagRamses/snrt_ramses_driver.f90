! P3 AMR leaf-cell state and CUDA conservation diagnostic for S_N.
! It is enabled only by SNRT_P1_DIAGNOSTIC=1 at runtime.
module snrt_ramses_driver
  use amr_parameters, only: MAXLEVEL
  implicit none

  private
  public :: snrt_ramses_diagnose_level, snrt_ramses_advance_level

contains

  subroutine snrt_ramses_diagnose_level(ilevel)
    use amr_commons, only: levelmin, nstep_coarse, myid
    use snrt_state, only: snrt_state_sync_level, snrt_ndirection
    use snrt_cuda_interface, only: snrt_cuda_available, &
         snrt_cuda_angular_reduce_tf32
    use snrt_cuda_ledger_interface, only: snrt_cuda_weighted_sum_fp32
    use iso_c_binding, only: c_float
    implicit none

    integer, intent(in) :: ilevel
    integer, parameter :: nbin = 16
    character(len=16) :: env_value
    integer :: env_length, env_status
    integer :: i, idir, ibin, nrow, nnew, ierr
    integer, save :: last_checked(MAXLEVEL) = -1
    logical, save :: enabled_resolved = .false.
    logical, save :: enabled = .false.
    real(c_float) :: max_tensor_error, max_ledger_error, reference
    real(c_float), allocatable :: directional(:), weights(:), projection(:), &
         binned(:), scalar(:)

    if (.not. enabled_resolved) then
       env_value = ''
       call get_environment_variable('SNRT_P1_DIAGNOSTIC', env_value, &
            length=env_length, status=env_status)
       enabled = env_status == 0 .and. env_length == 1 .and. &
            env_value(1:1) == '1'
       enabled_resolved = .true.
       if (enabled .and. myid == 1) then
          write(*,'(A)') ' SNRT P3 diagnostic enabled: persistent AMR leaf state'
       endif
    endif
    if (.not. enabled) return
    if (ilevel < levelmin .or. ilevel > MAXLEVEL) return
    if (last_checked(ilevel) == nstep_coarse) return

    call snrt_state_sync_level(ilevel, nrow, nnew)
    last_checked(ilevel) = nstep_coarse
    if (nrow == 0) return
    if (snrt_cuda_available() <= 0) then
       if (myid == 1) write(*,'(A)') ' SNRT P3 disabled: no CUDA device is visible'
       enabled = .false.
       return
    endif

    allocate(directional(nrow * snrt_ndirection), weights(snrt_ndirection), &
         projection(snrt_ndirection * nbin), binned(nrow * nbin), scalar(nrow))
    directional = 1.0_c_float
    weights = 1.0_c_float / real(snrt_ndirection, c_float)
    projection = 0.0_c_float
    do idir = 1, snrt_ndirection
       ibin = 1 + mod(idir - 1, nbin)
       projection((idir - 1) * nbin + ibin) = weights(idir)
    enddo

    call snrt_cuda_angular_reduce_tf32(directional, projection, binned, nrow, &
         snrt_ndirection, nbin, ierr)
    if (ierr /= 0) then
       if (myid == 1) write(*,'(A,I0)') ' SNRT P3 TF32 reduction failed, code=', ierr
       deallocate(directional, weights, projection, binned, scalar)
       enabled = .false.
       return
    endif
    call snrt_cuda_weighted_sum_fp32(directional, weights, scalar, nrow, &
         snrt_ndirection, ierr)
    if (ierr /= 0) then
       if (myid == 1) write(*,'(A,I0)') ' SNRT P3 ledger reduction failed, code=', ierr
       deallocate(directional, weights, projection, binned, scalar)
       enabled = .false.
       return
    endif

    reference = sum(weights)
    max_tensor_error = 0.0_c_float
    max_ledger_error = 0.0_c_float
    do i = 1, nrow
       max_tensor_error = max(max_tensor_error, &
            abs(sum(binned((i - 1) * nbin + 1:i * nbin)) - reference))
       max_ledger_error = max(max_ledger_error, abs(scalar(i) - reference))
    enddo
    if (myid == 1) write(*,'(A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
         ' SNRT P3 leaf diagnostic level=', ilevel, ' rows=', nrow, &
         ' new_slots=', nnew, ' tensor_abs=', max_tensor_error, &
         ' ledger_abs=', max_ledger_error

    deallocate(directional, weights, projection, binned, scalar)
  end subroutine snrt_ramses_diagnose_level

  subroutine snrt_ramses_advance_level(ilevel)
    use amr_parameters, only: sink, sink_AGN
    use amr_commons, only: levelmin, nstep_coarse, myid, dtnew, boxlen, &
         icoarse_min, icoarse_max, ncpu, nrestart
    use hydro_commons, only: uold
    use pm_commons, only: nsink, xsink, idsink, agn_pending_erg
    use snrt_state, only: snrt_ndirection, snrt_ngroups, snrt_intensity, &
         snrt_nslot, &
         snrt_neutral_fraction, snrt_hydrogen_ii, snrt_helium_ii, &
         snrt_helium_iii, snrt_state_get_slot
    use snrt_spectral_contract, only: &
         snrt_group_mean_energy_ev, snrt_group_energy_fraction, &
         snrt_group_cross_section_cm2, snrt_group_cross_section_hei_cm2, &
         snrt_group_cross_section_heii_cm2, &
         snrt_group_photoelectron_excess_energy_ev, &
         snrt_group_photoelectron_excess_hei_ev, &
         snrt_group_photoelectron_excess_heii_ev, &
         snrt_group_energy_fraction_sum, &
         snrt_group_unrepresented_energy_fraction, &
         snrt_spectral_contract_load_from_environment, &
         snrt_spectral_contract_error_name, snrt_spectral_contract_loaded, &
         snrt_spectral_contract_error_message, &
         snrt_spectral_contract_runtime_allowed, snrt_spectral_contract_status, &
         snrt_spectral_contract_source_id
    use snrt_amr_topology, only: snrt_amr_build_same_level_neighbors
    use snrt_transport_step, only: snrt_transport_absorb_multigroup_prepared_dust_trial
    use snrt_rt_transaction, only: snrt_rt_iteration_config, &
         snrt_transaction_contract_version, &
         snrt_rt_transaction_snapshot, snrt_transaction_load_config, &
         snrt_transaction_begin, snrt_transaction_restore, &
         snrt_transaction_commit_level, snrt_transaction_check_convergence, &
         snrt_transaction_reduce_decision, snrt_transaction_reduce_sum, &
         snrt_transaction_failure_requested, &
         snrt_transaction_failure_name, snrt_transaction_error_message, &
         snrt_failure_none, snrt_failure_partition, snrt_failure_chemistry, &
         snrt_failure_receiver, snrt_failure_transport, snrt_failure_convergence
    use snrt_dust_transaction, only: snrt_dust_validate_ledgers, &
         snrt_dust_transaction_ok
    use snrt_angular_quadrature, only: snrt_angular_init
    use snrt_agn_locator, only: snrt_agn_find_local_leaf
    use snrt_agn_source, only: snrt_c_cgs, snrt_agn_photon_budget_energy, &
         snrt_agn_deposit_transaction, snrt_agn_source_commit
    use snrt_agn_efficiency, only: snrt_agn_rt_requested, snrt_agn_reference_active
    use snrt_nlte_coupling, only: snrt_nlte_primordial_optical_depth_groups
    use snrt_thermochemistry, only: snrt_thermochemistry_result, &
         snrt_secondary_tables_load_from_environment, snrt_secondary_tables_loaded, &
         snrt_thermochemistry_ok, snrt_thermochemistry_error_name, &
         snrt_thermochemistry_error_message, snrt_thermochemistry_advance_cell, &
         snrt_nhelium_per_hydrogen, &
         snrt_mean_molecular_weight, snrt_inventory_tolerance
    use snrt_cuda_interface, only: snrt_cuda_available
    use amr_parameters, only: dp, ndim, spin_bh, mad_jet, X_floor
    use hydro_parameters, only: gamma
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_c_binding, only: c_float
    use omp_lib, only: omp_get_wtime
    implicit none

    integer, intent(in) :: ilevel
    character(len=1024) :: env_value
    integer :: env_length, env_status, read_status
    integer :: i, isink, igroup, ierr, nleaf, n_interface_face
    integer :: icell, islot, ilevel_found
    integer :: energy_index
    integer :: spectral_status, thermochemistry_status
    integer :: transaction_status, transaction_iteration, ledger_status
    integer :: local_transaction_failure, local_transaction_converged
    integer :: global_transaction_failure, global_transaction_converged
    integer :: convergence_status
    integer :: chemistry_failures
    integer :: n_locator_calls, n_active_sources
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), angular_weight(snrt_ndirection)
    real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(dp) :: dt_s, dx_code, cell_volume_code, cdt_over_dx
    real(dp) :: rho_code, hydrogen_ionized_fraction, helium_ionized_fraction
    real(dp) :: helium_double_ionized_fraction
    real(dp) :: n_hydrogen_cm3, n_helium_cm3
    real(dp) :: neutral_hydrogen_code, neutral_helium_i_code, neutral_helium_ii_code
    real(dp) :: temperature_k, internal_energy, kinetic_energy, molecular_weight
    real(dp) :: tau_dp(snrt_ngroups), tau_hi_dp(snrt_ngroups)
    real(dp) :: tau_hei_dp(snrt_ngroups), tau_heii_dp(snrt_ngroups)
    real(dp) :: unassigned_absorption_total, ledger_relative_error
    real(dp) :: global_unassigned_absorption
    real(dp) :: excess_energy_ev(3,snrt_ngroups)
    real(dp) :: deposited_density
    real(dp) :: wall_start
    real(dp) :: transaction_residual, global_transaction_residual
    real(dp) :: wall_sub
    real(dp) :: t_setup, t_topology, t_nlte, t_source, t_transport, t_coupling
    real(dp) :: t_source_overhead, t_locator, t_budget, t_deposit
    real(c_float), allocatable :: optical_depth(:,:), optical_depth_species(:,:,:), &
         optical_depth_dust(:,:)
    real(c_float), allocatable :: available_species_transport(:,:)
    real(c_float), allocatable :: optical_depth_hydrogen(:,:), &
         optical_depth_helium_i(:,:), optical_depth_helium_ii(:,:)
    real(c_float), allocatable :: absorbed_group(:,:), raw_group(:,:), &
         absorbed_hhe_group_species(:,:,:), absorbed_dust_group(:,:), returned_group(:,:)
    real(c_float), allocatable :: incoming_intensity(:,:,:), trial_intensity(:,:,:)
    real(c_float), allocatable :: coarse_flux_trial(:,:,:)
    real(c_float), allocatable :: iteration_tau(:,:), iteration_species_tau(:,:,:)
    real(c_float), allocatable :: target_tau(:,:), target_species_tau(:,:,:)
    real(dp), allocatable :: start_hydrogen_ii(:), start_helium_ii(:), &
         start_helium_iii(:), start_neutral_hydrogen(:), level_thermal(:)
    real(dp), allocatable :: current_hydrogen_ii(:), current_helium_ii(:), &
         current_helium_iii(:), relaxed_hydrogen_ii(:), relaxed_helium_ii(:), &
         relaxed_helium_iii(:), trial_hydrogen_ii(:), trial_helium_ii(:), &
         trial_helium_iii(:), trial_neutral_hydrogen(:), trial_thermal(:)
    real(dp), allocatable :: rho_level(:), temperature_level(:)
    real(dp), allocatable :: trial_heating_rate(:), trial_unassigned(:)
    real(dp), allocatable :: trial_absorbed_species(:,:,:)
    real(dp), allocatable :: current_fraction(:,:), target_fraction(:,:)
    real(dp), allocatable :: emitted_groups(:), luminosity_groups(:)
    logical :: enabled, source_ok, accounting_identity_ok
    logical :: transaction_converged, transaction_active
    logical :: transaction_diagnostic_mode
    logical :: hydro_state_invalid
    type(snrt_rt_iteration_config) :: transaction_config
    type(snrt_rt_transaction_snapshot) :: transaction
    type(snrt_thermochemistry_result) :: chemistry_result
    logical, save :: enabled_resolved = .false.
    logical, save :: enabled_latched = .false.
    logical, save :: spectral_contract_resolved = .false.
    logical, save :: spectral_contract_ok = .false.
    logical, save :: thermochemistry_resolved = .false.
    logical, save :: thermochemistry_contract_ok = .false.
    logical, save :: transaction_diagnostic_mode_latched = .false.
    logical, save :: transaction_config_reported = .false.
    logical, save :: dust_scaffold_reported = .false.
    real(dp), save :: reduced_c = 0.01d0
    integer, save :: level_filter = -1
    character(len=13), parameter :: dust_tau_mode = 'ZERO_SCAFFOLD'

    enabled = .false.
    if (.not. enabled_resolved) then
       enabled_latched = snrt_agn_rt_requested()
       if (enabled_latched .and. sink .and. sink_AGN .and. .not.snrt_agn_reference_active()) then
          if(myid==1)write(*,*)'AGN source ownership conflict: legacy feedback plus live SNRT is not approved'
          call clean_stop
       end if
       env_value = ''
       call get_environment_variable('SNRT_REDUCED_C', env_value, &
            length=env_length, status=env_status)
       if (env_status == 0 .and. env_length > 0) then
          read(env_value(1:env_length),*,iostat=read_status) reduced_c
          if (read_status /= 0 .or. reduced_c <= 0.0d0 .or. reduced_c > 1.0d0) &
               reduced_c = 0.01d0
       end if
       env_value = ''
       call get_environment_variable('SNRT_RT_LEVEL', env_value, &
            length=env_length, status=env_status)
       if (env_status == 0 .and. env_length > 0) then
          read(env_value(1:env_length),*,iostat=read_status) level_filter
          if (read_status /= 0 .or. level_filter < levelmin .or. &
               level_filter > MAXLEVEL) level_filter = -1
       end if
       env_value = ''
       call get_environment_variable('SNRT_RT_TX_DIAGNOSTIC_MODE', env_value, &
            length=env_length, status=env_status)
       transaction_diagnostic_mode_latched = env_status == 0 .and. env_length == 1 .and. &
            env_value(1:1) == '1'
       enabled_resolved = .true.
       if (enabled_latched .and. myid == 1) write(*,'(A,F8.4)') &
            ' SNRT S_N RT enabled; reduced speed factor=', reduced_c
       if (enabled_latched .and. myid == 1 .and. level_filter > 0) &
            write(*,'(A,I0)') ' SNRT RT level filter=', level_filter
    end if
    ! Runtime control is latched once per process.  Re-reading the
    ! environment mid-step could strand an already committed accounting
    ! marker while leaving the transport state in a different mode.
    enabled = enabled_latched
    transaction_diagnostic_mode = transaction_diagnostic_mode_latched
    if (.not. enabled) return
    if (myid == 1 .and. .not. dust_scaffold_reported) then
       write(*,'(A,A)') ' SNRT dust optical-depth mode=', trim(dust_tau_mode)
       write(*,'(A)') '   dust ledger is trial-only; no dust thermal/momentum/abundance commit'
       dust_scaffold_reported = .true.
    end if
    ! Do not guess unknown historical photon ownership by either replaying
    ! the ledger or silently rebasing away real accretion. Startup rejects
    ! these modes with sink enabled; defend restored sink arrays here too.
    if (nsink>0 .and. (ncpu>1 .or. nrestart>0)) then
       if(myid==1)write(*,*)'Live SNRT AGN requires serial fresh start until pending energy and photon state support restart/migration'
       call clean_stop
    end if
    if (.not. spectral_contract_resolved) then
       call snrt_spectral_contract_load_from_environment(spectral_status)
       spectral_contract_resolved = .true.
       spectral_contract_ok = spectral_status == 0 .and. &
            snrt_spectral_contract_loaded .and. &
            snrt_spectral_contract_runtime_allowed
       if(snrt_agn_reference_active())spectral_contract_ok=spectral_contract_ok.and. &
            trim(snrt_spectral_contract_status)=='reference_control'
       if (myid == 1) then
          if (spectral_contract_ok) then
             write(*,'(A,I0,A,A)') ' SNRT spectral contract loaded: groups=', &
                  snrt_ngroups, ' status=', trim(snrt_spectral_contract_status)
             write(*,'(A,A)') '   source: ', trim(snrt_spectral_contract_source_id)
             write(*,'(A,F12.8,A,F12.8)') '   represented energy fraction=', &
                  snrt_group_energy_fraction_sum, ' unrepresented=', &
                  snrt_group_unrepresented_energy_fraction
          else if (spectral_status == 0 .and. snrt_spectral_contract_loaded) then
             write(*,'(A,A)') ' SNRT RT disabled: spectral contract status is not runtime-admissible: ', &
                  trim(snrt_spectral_contract_status)
          else
             write(*,'(A,I0,A,A)') ' SNRT RT disabled: spectral contract error=', &
                  spectral_status, ' (', trim(snrt_spectral_contract_error_name(spectral_status))//')'
             if (len_trim(snrt_spectral_contract_error_message) > 0) &
                  write(*,'(A,A)') '   detail: ', trim(snrt_spectral_contract_error_message)
          end if
       end if
    end if
    if (.not. spectral_contract_ok) return
    if (.not. thermochemistry_resolved) then
       call snrt_secondary_tables_load_from_environment(thermochemistry_status)
       thermochemistry_resolved = .true.
       thermochemistry_contract_ok = thermochemistry_status == snrt_thermochemistry_ok .and. &
            snrt_secondary_tables_loaded
       if (myid == 1) then
          if (thermochemistry_contract_ok) then
             write(*,'(A)') ' SNRT native H/He thermochemistry contract loaded'
          else
             write(*,'(A,I0,A,A)') ' SNRT RT disabled: thermochemistry contract error=', &
                  thermochemistry_status, ' (', &
                  trim(snrt_thermochemistry_error_name(thermochemistry_status))//')'
             if (len_trim(snrt_thermochemistry_error_message) > 0) &
                  write(*,'(A,A)') '   detail: ', trim(snrt_thermochemistry_error_message)
          end if
       end if
    end if
    if (.not. thermochemistry_contract_ok) return
    call snrt_transaction_load_config(spectral_contract_ok .and. &
         .not. transaction_diagnostic_mode, &
         transaction_config, transaction_status, env_value)
    if (transaction_status /= 0) then
       if (myid == 1) write(*,'(A,I0,A,A)') &
            ' SNRT RT disabled: transaction configuration error=', transaction_status, ' (', &
            trim(snrt_transaction_error_message(transaction_status))//')'
       if (myid == 1 .and. len_trim(env_value) > 0) &
            write(*,'(A,A)') '   detail: ', trim(env_value)
       return
    end if
    if (myid == 1 .and. .not. transaction_config_reported) then
       write(*,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,F6.3)') &
            ' SNRT RT transaction contract=', snrt_transaction_contract_version, &
            ' max_iter=', transaction_config%max_iterations, &
            ' frac_tol=', transaction_config%fraction_absolute_tolerance, &
            ' tau_tol=', transaction_config%tau_relative_tolerance, &
            ' tau_floor=', transaction_config%tau_floor, &
            ' relaxation=', transaction_config%relaxation
       if (transaction_config%failure_stage /= 0) &
            write(*,'(A,I0,A,I0)') '   test failure injection stage=', &
            transaction_config%failure_stage, ' leaf=', transaction_config%failure_leaf
       if (transaction_diagnostic_mode) write(*,'(A)') &
            '   NONPRODUCTION diagnostic mode: failure injection is permitted'
       transaction_config_reported = .true.
    end if
    if (ilevel < levelmin .or. dtnew(ilevel) <= 0.0d0) return
    if (level_filter > 0 .and. ilevel /= level_filter) return
    local_transaction_failure = snrt_failure_none
    if (snrt_cuda_available() <= 0) local_transaction_failure = snrt_failure_transport
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0 .or. &
         global_transaction_failure /= snrt_failure_none) then
       if (myid == 1 .and. global_transaction_failure == snrt_failure_transport) &
            write(*,'(A)') ' SNRT RT preflight failed closed: no CUDA device is visible on every rank'
       call clean_stop
       return
    end if
    local_transaction_failure = snrt_failure_none
    if (.not. allocated(uold)) then
       local_transaction_failure = snrt_failure_transport
    else if (size(uold,1) < 1 .or. size(uold,2) < ndim+2) then
       local_transaction_failure = snrt_failure_transport
    end if
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0 .or. &
         global_transaction_failure /= snrt_failure_none) then
       if (myid == 1) write(*,'(A)') &
            ' SNRT RT preflight failed closed: hydro thermal field is unavailable on one or more ranks'
       call clean_stop
       return
    end if

    wall_start = omp_get_wtime()
    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
    dt_s = dtnew(ilevel) * scale_t
    dx_code = boxlen / dble(icoarse_max - icoarse_min + 1) * 0.5d0**ilevel
    local_transaction_failure = snrt_failure_none
    if (scale_l <= 0.0d0 .or. scale_d <= 0.0d0 .or. scale_v <= 0.0d0 .or. &
         scale_nH <= 0.0d0 .or. scale_T2 <= 0.0d0 .or. dt_s <= 0.0d0 .or. &
         dx_code <= 0.0d0) local_transaction_failure = snrt_failure_transport
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0 .or. &
         global_transaction_failure /= snrt_failure_none) then
       if (myid == 1) write(*,'(A)') &
            ' SNRT RT preflight failed closed: invalid unit conversion or timestep on one or more ranks'
       call clean_stop
       return
    end if
    cell_volume_code = dx_code**ndim
    cdt_over_dx = snrt_c_cgs * reduced_c * dt_s / (dx_code * scale_l)
    energy_index = ndim + 2
    t_setup = omp_get_wtime() - wall_start

    wall_start = omp_get_wtime()
    call snrt_amr_build_same_level_neighbors(ilevel, leaf_cell, leaf_slot, &
         neighbor, nleaf, n_interface_face)
    t_topology = omp_get_wtime() - wall_start
    allocate(optical_depth(nleaf,snrt_ngroups), &
         optical_depth_species(nleaf,snrt_ngroups,3), &
         optical_depth_dust(nleaf,snrt_ngroups), &
         available_species_transport(nleaf,3), &
         optical_depth_hydrogen(nleaf,snrt_ngroups), &
         optical_depth_helium_i(nleaf,snrt_ngroups), &
         optical_depth_helium_ii(nleaf,snrt_ngroups), &
         absorbed_group(nleaf,snrt_ngroups), raw_group(nleaf,snrt_ngroups), &
         absorbed_hhe_group_species(nleaf,snrt_ngroups,3), &
         absorbed_dust_group(nleaf,snrt_ngroups), returned_group(nleaf,snrt_ngroups), &
         incoming_intensity(snrt_ndirection,snrt_ngroups,nleaf), &
         trial_intensity(snrt_ndirection,snrt_ngroups,nleaf), &
         coarse_flux_trial(snrt_ndirection,snrt_ngroups,size(snrt_intensity,3)), &
         iteration_tau(nleaf,snrt_ngroups), &
         iteration_species_tau(nleaf,snrt_ngroups,3), &
         target_tau(nleaf,snrt_ngroups), &
         target_species_tau(nleaf,snrt_ngroups,3), &
         start_hydrogen_ii(nleaf), start_helium_ii(nleaf), &
         start_helium_iii(nleaf), start_neutral_hydrogen(nleaf), &
         level_thermal(nleaf), current_hydrogen_ii(nleaf), &
         current_helium_ii(nleaf), current_helium_iii(nleaf), &
         relaxed_hydrogen_ii(nleaf), relaxed_helium_ii(nleaf), &
         relaxed_helium_iii(nleaf), trial_hydrogen_ii(nleaf), &
         trial_helium_ii(nleaf), trial_helium_iii(nleaf), &
         trial_neutral_hydrogen(nleaf), trial_thermal(nleaf), &
         rho_level(nleaf), temperature_level(nleaf), &
         trial_heating_rate(nleaf), trial_unassigned(nleaf), &
         trial_absorbed_species(nleaf,3,snrt_ngroups), &
         current_fraction(nleaf,3), target_fraction(nleaf,3))
    optical_depth_dust = 0.0_c_float
    raw_group = 0.0_c_float
    absorbed_hhe_group_species = 0.0_c_float
    absorbed_dust_group = 0.0_c_float
    returned_group = 0.0_c_float
    call snrt_angular_init(direction_dp, angular_weight)
    excess_energy_ev = 0.0d0
    excess_energy_ev(1,:) = snrt_group_photoelectron_excess_energy_ev
    excess_energy_ev(2,:) = snrt_group_photoelectron_excess_hei_ev
    excess_energy_ev(3,:) = snrt_group_photoelectron_excess_heii_ev

    wall_start = omp_get_wtime()
    hydro_state_invalid = .false.
    do i = 1, nleaf
       icell = leaf_cell(i)
       islot = leaf_slot(i)
       rho_code = 0.0d0
       level_thermal(i) = 0.0d0
       temperature_level(i) = 0.0d0
       if (icell >= 1 .and. icell <= size(uold,1) .and. size(uold,2) >= ndim+2) then
          if (ieee_is_finite(uold(icell,1))) then
             rho_code = max(0.0d0,uold(icell,1))
          else
             hydro_state_invalid = .true.
          end if
          if (ieee_is_finite(uold(icell,ndim+2))) then
             level_thermal(i) = uold(icell,ndim+2)
          else
             hydro_state_invalid = .true.
          end if
       else
          hydro_state_invalid = .true.
       end if
       rho_level(i) = rho_code
       hydrogen_ionized_fraction = 0.0d0
       helium_ionized_fraction = 0.0d0
       helium_double_ionized_fraction = 0.0d0
       if (allocated(snrt_hydrogen_ii)) then
          hydrogen_ionized_fraction = max(0.0d0, min(1.0d0, &
               real(snrt_hydrogen_ii(islot),dp)))
       else if (allocated(snrt_neutral_fraction)) then
          hydrogen_ionized_fraction = 1.0d0 - max(0.0d0, min(1.0d0, &
               real(snrt_neutral_fraction(islot),dp)))
       end if
       if (allocated(snrt_helium_ii)) helium_ionized_fraction = max(0.0d0, &
            min(1.0d0, real(snrt_helium_ii(islot),dp)))
       if (allocated(snrt_helium_iii)) helium_double_ionized_fraction = &
            max(0.0d0, min(1.0d0, real(snrt_helium_iii(islot),dp)))
       if (helium_ionized_fraction + helium_double_ionized_fraction > 1.0d0) &
            helium_double_ionized_fraction = max(0.0d0, &
            1.0d0-helium_ionized_fraction)
       start_hydrogen_ii(i) = hydrogen_ionized_fraction
       start_helium_ii(i) = helium_ionized_fraction
       start_helium_iii(i) = helium_double_ionized_fraction
       start_neutral_hydrogen(i) = max(0.0d0,1.0d0-hydrogen_ionized_fraction)
       if (icell >= 1 .and. icell <= size(uold,1) .and. &
            size(uold,2) >= ndim+2) then
          if (rho_code > 0.0d0 .and. level_thermal(i) > 0.0d0 .and. &
               all(ieee_is_finite(uold(icell,2:ndim+1)))) then
             kinetic_energy = 0.5d0 * sum(uold(icell,2:ndim+1)**2) / rho_code
             internal_energy = level_thermal(i) - kinetic_energy
             molecular_weight = snrt_mean_molecular_weight(hydrogen_ionized_fraction, &
                  helium_ionized_fraction, helium_double_ionized_fraction)
             if (ieee_is_finite(internal_energy) .and. internal_energy > 0.0d0) &
                  temperature_level(i) = max(1.0d0,(gamma-1.0d0)*internal_energy/rho_code * &
                  scale_T2*molecular_weight)
          end if
       end if
       neutral_hydrogen_code = rho_code * (1.0d0-hydrogen_ionized_fraction)
       neutral_helium_i_code = rho_code * snrt_nhelium_per_hydrogen * &
            (1.0d0-helium_ionized_fraction-helium_double_ionized_fraction)
       neutral_helium_ii_code = rho_code * snrt_nhelium_per_hydrogen * &
            helium_ionized_fraction
       available_species_transport(i,1) = real(max(0.0d0, neutral_hydrogen_code), c_float)
       available_species_transport(i,2) = real(max(0.0d0, neutral_helium_i_code), c_float)
       available_species_transport(i,3) = real(max(0.0d0, neutral_helium_ii_code), c_float)
       call snrt_nlte_primordial_optical_depth_groups(neutral_hydrogen_code, &
            neutral_helium_i_code, neutral_helium_ii_code, scale_nH, dt_s, &
            snrt_group_cross_section_cm2, snrt_group_cross_section_hei_cm2, &
            snrt_group_cross_section_heii_cm2, tau_dp, tau_hi_dp, tau_hei_dp, &
            tau_heii_dp, ierr, reduced_c)
       if (ierr /= 0) then
          optical_depth(i,:) = 0.0_c_float
          optical_depth_hydrogen(i,:) = 0.0_c_float
          optical_depth_helium_i(i,:) = 0.0_c_float
          optical_depth_helium_ii(i,:) = 0.0_c_float
       else
          optical_depth(i,:) = real(max(tau_dp,0.0d0),c_float)
          optical_depth_hydrogen(i,:) = real(max(tau_hi_dp,0.0d0),c_float)
          optical_depth_helium_i(i,:) = real(max(tau_hei_dp,0.0d0),c_float)
          optical_depth_helium_ii(i,:) = real(max(tau_heii_dp,0.0d0),c_float)
       end if
       optical_depth_species(i,:,1) = optical_depth_hydrogen(i,:)
       optical_depth_species(i,:,2) = optical_depth_helium_i(i,:)
       optical_depth_species(i,:,3) = optical_depth_helium_ii(i,:)
    end do
    t_nlte = omp_get_wtime() - wall_start

    ! A non-finite hydro receiver value must never be replaced by the local
    ! zero initializer and written back to uold on an untouched leaf.  Make
    ! this a collective pre-source failure so no rank can deposit photons or
    ! enter the RT transaction with an invalid thermal baseline.
    local_transaction_failure = snrt_failure_none
    if (hydro_state_invalid) local_transaction_failure = snrt_failure_transport
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0 .or. &
         global_transaction_failure /= snrt_failure_none) then
       if (myid == 1) write(*,'(A)') &
            ' SNRT RT disabled: non-finite hydro state on one or more ranks'
       deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, optical_depth_species, &
            optical_depth_dust, available_species_transport, optical_depth_hydrogen, &
            optical_depth_helium_i, optical_depth_helium_ii, absorbed_group, raw_group, &
            absorbed_hhe_group_species, absorbed_dust_group, returned_group, &
            incoming_intensity, trial_intensity, coarse_flux_trial, iteration_tau, &
            iteration_species_tau, target_tau, target_species_tau, start_hydrogen_ii, &
            start_helium_ii, start_helium_iii, start_neutral_hydrogen, level_thermal, &
            current_hydrogen_ii, current_helium_ii, current_helium_iii, &
            relaxed_hydrogen_ii, relaxed_helium_ii, relaxed_helium_iii, &
            trial_hydrogen_ii, trial_helium_ii, trial_helium_iii, &
            trial_neutral_hydrogen, trial_thermal, rho_level, temperature_level, &
            trial_heating_rate, trial_unassigned, trial_absorbed_species, &
            current_fraction, target_fraction)
       call clean_stop
       return
    end if

    wall_start = omp_get_wtime()
    t_locator = 0.0d0
    t_budget = 0.0d0
    t_deposit = 0.0d0
    n_locator_calls = 0
    n_active_sources = 0
    ! Accepted event energy follows sink creation/merger arrays. Coarse rate
    ! estimates and present-day spin efficiency do not fund these photons.
    accounting_identity_ok = .true.
    if (nsink > 0) then
       do i = 1, nsink
          if (idsink(i) <= 0) accounting_identity_ok = .false.
          do isink = i + 1, nsink
             if (idsink(i) == idsink(isink)) accounting_identity_ok = .false.
          end do
       end do
    end if
    if (.not. accounting_identity_ok) then
       if (myid == 1) write(*,'(A)') &
            ' SNRT AGN source skipped: idsink identity map is invalid'
    end if

    if (accounting_identity_ok .and. nsink > 0) then
       if (allocated(agn_pending_erg) .and. allocated(xsink)) then
          if (size(agn_pending_erg)>=nsink .and. size(xsink,1)>=nsink) then
             allocate(emitted_groups(snrt_ngroups), luminosity_groups(snrt_ngroups))
             do isink = 1, nsink
                if (.not.ieee_is_finite(agn_pending_erg(isink)) .or. agn_pending_erg(isink)<0d0) then
                   if(myid==1)write(*,*) 'Invalid accepted AGN radiative energy for sink ',idsink(isink)
                   call clean_stop
                endif
                if(agn_pending_erg(isink)==0d0)cycle

                wall_sub = omp_get_wtime()
                ! The locator returns a leaf owned by this MPI rank or zero.
                ! The source loop is intentionally serial: local intensity and
                ! the sink-carried pending energy are shared mutable state, not
                ! OpenMP-threadprivate data.  A valid leaf therefore has one
                ! MPI owner for this source transaction.
                call snrt_agn_find_local_leaf(xsink(isink,1:ndim), icell, ilevel_found)
                t_locator = t_locator + omp_get_wtime() - wall_sub
                n_locator_calls = n_locator_calls + 1
                if (icell == 0 .or. ilevel_found /= ilevel) cycle
                islot = snrt_state_get_slot(icell)
                if (islot <= 0) cycle
                n_active_sources = n_active_sources + 1
                source_ok = .true.
                do igroup = 1, snrt_ngroups
                   wall_sub = omp_get_wtime()
                   call snrt_agn_photon_budget_energy(agn_pending_erg(isink), dt_s, &
                        snrt_group_energy_fraction(igroup), &
                        snrt_group_mean_energy_ev(igroup), luminosity_groups(igroup), &
                        emitted_groups(igroup),ierr)
                   if(ierr/=0)source_ok=.false.
                   t_budget = t_budget + omp_get_wtime() - wall_sub
                   if (.not. ieee_is_finite(luminosity_groups(igroup)) .or. &
                        .not. ieee_is_finite(emitted_groups(igroup)) .or. &
                        emitted_groups(igroup) < 0.0d0) source_ok = .false.
                end do
                if (source_ok) then
                   wall_sub = omp_get_wtime()
                   call snrt_agn_deposit_transaction(snrt_intensity, islot, &
                        emitted_groups, cell_volume_code, scale_l, scale_nH, &
                        angular_weight, deposited_density, ierr)
                   t_deposit = t_deposit + omp_get_wtime() - wall_sub
                   if (ierr /= 0) source_ok = .false.
                end if
                ! Failed source deposition retains fuel; successful source
                ! commit remains consumed even if subsequent transport retries.
                call snrt_agn_source_commit(agn_pending_erg(isink),source_ok)
             end do
             deallocate(emitted_groups, luminosity_groups)
          end if
       end if
    end if
    t_source = omp_get_wtime() - wall_start
    t_source_overhead = t_source - t_locator - t_budget - t_deposit

    ! The source phase above has already committed its own photon/accounting
    ! transaction.  The RT/chemistry transaction starts here, so a failed
    ! coupled level never erases a valid source deposit.
    do i = 1, nleaf
       do igroup = 1, snrt_ngroups
          incoming_intensity(:,igroup,i) = snrt_intensity(:,igroup,leaf_slot(i))
       end do
    end do
    call snrt_transaction_begin(transaction, snrt_intensity, leaf_slot, &
         snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
         snrt_neutral_fraction, level_thermal, transaction_status)
    transaction_active = transaction%active
    local_transaction_failure = snrt_failure_none
    if (transaction_status /= 0) local_transaction_failure = snrt_failure_transport
    call snrt_transaction_reduce_decision(local_transaction_failure, 0, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0) global_transaction_failure = snrt_failure_transport
    if (global_transaction_failure /= snrt_failure_none) then
       if (transaction_active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, transaction_status)
       if (myid == 1) write(*,'(A,A,A,I0)') &
            ' SNRT RT transaction could not start: class=', &
            trim(snrt_transaction_failure_name(global_transaction_failure)), &
            ' level=', ilevel
       deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, optical_depth_species, &
            optical_depth_dust, available_species_transport, optical_depth_hydrogen, &
            optical_depth_helium_i, optical_depth_helium_ii, absorbed_group, raw_group, &
            absorbed_hhe_group_species, absorbed_dust_group, returned_group, &
            incoming_intensity, trial_intensity, coarse_flux_trial, iteration_tau, &
            iteration_species_tau, target_tau, target_species_tau, start_hydrogen_ii, &
            start_helium_ii, start_helium_iii, start_neutral_hydrogen, level_thermal, &
            current_hydrogen_ii, current_helium_ii, current_helium_iii, &
            relaxed_hydrogen_ii, relaxed_helium_ii, relaxed_helium_iii, &
            trial_hydrogen_ii, trial_helium_ii, trial_helium_iii, &
            trial_neutral_hydrogen, trial_thermal, rho_level, temperature_level, &
            trial_heating_rate, trial_unassigned, trial_absorbed_species, &
            current_fraction, target_fraction)
       call clean_stop
       return
    end if

    current_hydrogen_ii = start_hydrogen_ii
    current_helium_ii = start_helium_ii
    current_helium_iii = start_helium_iii
    current_fraction(:,1) = current_hydrogen_ii
    current_fraction(:,2) = current_helium_ii
    current_fraction(:,3) = current_helium_iii
    iteration_species_tau = real(optical_depth_species,c_float)
    ! DUST-7 validates the total against the FP32 component sum.  Rebuild the
    ! total in the same precision/order at this boundary instead of relying
    ! on a separately rounded FP64 total from the NLTE helper.
    iteration_tau = sum(iteration_species_tau, dim=3) + optical_depth_dust
    chemistry_failures = 0
    unassigned_absorption_total = 0.0d0
    t_transport = 0.0d0
    t_coupling = 0.0d0
    transaction_converged = .false.
    global_transaction_converged = 0
    global_transaction_failure = snrt_failure_none
    wall_start = omp_get_wtime()

    do transaction_iteration = 1, transaction_config%max_iterations
       local_transaction_failure = snrt_failure_none
       transaction_residual = huge(1.0d0)
       transaction_converged = .false.
       trial_hydrogen_ii = start_hydrogen_ii
       trial_helium_ii = start_helium_ii
       trial_helium_iii = start_helium_iii
       trial_neutral_hydrogen = start_neutral_hydrogen
       trial_thermal = level_thermal
       trial_heating_rate = 0.0d0
       trial_unassigned = 0.0d0
       trial_absorbed_species = 0.0d0
       wall_sub = omp_get_wtime()
       call snrt_transport_absorb_multigroup_prepared_dust_trial(leaf_slot, neighbor, &
            cdt_over_dx, iteration_tau, iteration_species_tau, optical_depth_dust, &
            available_species_transport, incoming_intensity, trial_intensity, &
            coarse_flux_trial, raw_group, absorbed_hhe_group_species, &
            absorbed_dust_group, returned_group, absorbed_group, ierr, leaf_cell, ilevel)
       t_transport = t_transport + omp_get_wtime() - wall_sub
       if (ierr /= 0) then
          local_transaction_failure = snrt_failure_transport
       end if
       if (local_transaction_failure == snrt_failure_none) then
          ! CUDA is an external FP32 boundary.  Reject a corrupt trial before
          ! it can enter the chemistry or transaction commit path.
          if (any(.not. ieee_is_finite(trial_intensity)) .or. &
               any(trial_intensity < 0.0_c_float) .or. &
               any(.not. ieee_is_finite(absorbed_group)) .or. &
               any(absorbed_group < 0.0_c_float) .or. &
               any(.not. ieee_is_finite(coarse_flux_trial))) &
               local_transaction_failure = snrt_failure_transport
       end if
       if (local_transaction_failure == snrt_failure_none) then
          call snrt_dust_validate_ledgers(raw_group, absorbed_hhe_group_species, &
               absorbed_dust_group, returned_group, absorbed_group, &
               ledger_relative_error, ledger_status)
          if (ledger_status /= snrt_dust_transaction_ok) &
               local_transaction_failure = snrt_failure_receiver
       end if

       if (local_transaction_failure == snrt_failure_none) then
          unassigned_absorption_total = 0.0d0
          do i = 1, nleaf
             icell = leaf_cell(i)
             rho_code = rho_level(i)
             if (rho_code <= 0.0d0) cycle
             ! With no absorbed photons this bundle has no local chemistry
             ! source to advance.  Do not reject an otherwise untouched cell
             ! merely because its hydro internal-energy reconstruction cannot
             ! provide a positive temperature.
             if (sum(abs(real(absorbed_group(i,:),dp))) <= 0.0d0) cycle

             available_species_code = real(available_species_transport(i,:),dp)
             chemistry_cell_ok = .true.
             trial_unassigned(i) = 0.0d0
             do igroup = 1, snrt_ngroups
                opacity_species(1) = real(iteration_species_tau(i,igroup,1),dp)
                opacity_species(2) = real(iteration_species_tau(i,igroup,2),dp)
                opacity_species(3) = real(iteration_species_tau(i,igroup,3),dp)
                call snrt_partition_absorption(real(absorbed_group(i,igroup),dp), &
                     opacity_species, available_species_code, &
                     trial_absorbed_species(i,:,igroup), ierr, &
                     unassigned_absorption_code, inventory_scale_code=&
                     sum(abs(real(available_species_transport(i,:),dp))))
                trial_unassigned(i) = trial_unassigned(i) + &
                     max(0.0d0, unassigned_absorption_code)
                unassigned_absorption_total = unassigned_absorption_total + &
                     max(0.0d0, unassigned_absorption_code)
                if (ierr /= snrt_thermochemistry_ok) then
                   local_transaction_failure = snrt_failure_partition
                   chemistry_cell_ok = .false.
                   exit
                end if
                if (max(real(absorbed_group(i,igroup),dp), &
                     sum(abs(real(available_species_transport(i,:),dp)))) > 0.0d0) then
                   if (max(0.0d0,unassigned_absorption_code) > &
                        snrt_inventory_tolerance(max(real(absorbed_group(i,igroup),dp), &
                        sum(abs(real(available_species_transport(i,:),dp)))))) then
                      local_transaction_failure = snrt_failure_unassigned
                      chemistry_cell_ok = .false.
                      exit
                   end if
                end if
             end do
             if (.not. chemistry_cell_ok) then
                chemistry_failures = chemistry_failures + 1
                cycle
             end if

             trial_unassigned(i) = 0.0d0
             if (sum(abs(real(absorbed_hhe_group_species(i,:,:),dp))) <= 0.0d0 .and. &
                  sum(abs(real(absorbed_group(i,:),dp))) > 0.0d0) then
                ! A future nonzero dust source must supply a dust receiver;
                ! fail closed instead of pretending dust absorption heats gas.
                local_transaction_failure = snrt_failure_receiver
                exit
             end if
             do igroup = 1, snrt_ngroups
                trial_absorbed_species(i,:,igroup) = &
                     real(absorbed_hhe_group_species(i,igroup,:),dp)
             end do
             if (snrt_transaction_failure_requested(transaction_config, &
                  snrt_failure_partition, i)) then
                local_transaction_failure = snrt_failure_partition
                exit
             end if

             if (temperature_level(i) <= 0.0d0 .or. &
                  .not. ieee_is_finite(temperature_level(i)) .or. &
                  .not. ieee_is_finite(level_thermal(i))) then
                local_transaction_failure = snrt_failure_chemistry
                chemistry_failures = chemistry_failures + 1
                cycle
             end if
       n_hydrogen_cm3 = rho_code * scale_nH
       n_helium_cm3 = n_hydrogen_cm3 * snrt_nhelium_per_hydrogen
       call snrt_thermochemistry_advance_cell(n_hydrogen_cm3, n_helium_cm3, &
            scale_nH, temperature_level(i), dt_s, start_hydrogen_ii(i), &
            start_helium_ii(i), start_helium_iii(i), &
            trial_absorbed_species(i,:,:), excess_energy_ev, chemistry_result)
             if (chemistry_result%ierr /= snrt_thermochemistry_ok .or. &
                  .not. ieee_is_finite(chemistry_result%x_hydrogen_ii) .or. &
                  .not. ieee_is_finite(chemistry_result%x_helium_ii) .or. &
                  .not. ieee_is_finite(chemistry_result%x_helium_iii)) then
                local_transaction_failure = snrt_failure_chemistry
                chemistry_failures = chemistry_failures + 1
                cycle
             end if
             if (snrt_transaction_failure_requested(transaction_config, &
                  snrt_failure_chemistry, i)) then
                local_transaction_failure = snrt_failure_chemistry
                exit
             end if

             trial_hydrogen_ii(i) = chemistry_result%x_hydrogen_ii
             trial_helium_ii(i) = chemistry_result%x_helium_ii
             trial_helium_iii(i) = chemistry_result%x_helium_iii
             trial_neutral_hydrogen(i) = max(0.0d0, &
                  1.0d0-trial_hydrogen_ii(i))
             trial_heating_rate(i) = chemistry_result%heating_rate_erg_cm3_s
             if (.not. ieee_is_finite(trial_heating_rate(i)) .or. &
                  trial_heating_rate(i) < 0.0d0) then
                local_transaction_failure = snrt_failure_receiver
                chemistry_failures = chemistry_failures + 1
                cycle
             end if
             trial_thermal(i) = level_thermal(i) + trial_heating_rate(i)*dt_s / &
                  (scale_d*scale_v**2)
             if (.not. ieee_is_finite(trial_thermal(i)) .or. &
                  trial_thermal(i) <= 0.0d0) then
                local_transaction_failure = snrt_failure_receiver
                chemistry_failures = chemistry_failures + 1
                cycle
             end if
             if (snrt_transaction_failure_requested(transaction_config, &
                  snrt_failure_receiver, i)) then
                local_transaction_failure = snrt_failure_receiver
                exit
             end if
          end do
       end if

       ! Validate the arrays after partition and thermochemistry have filled
       ! them, immediately before they can feed the fixed-point predictor or
       ! a successful commit.  The transport-only finite checks above remain
       ! before chemistry; this block covers the actual candidate receiver.
       if (local_transaction_failure == snrt_failure_none) then
          if (any(.not. ieee_is_finite(trial_absorbed_species)) .or. &
               any(trial_absorbed_species < 0.0d0) .or. &
               any(.not. ieee_is_finite(trial_heating_rate)) .or. &
               any(trial_heating_rate < 0.0d0) .or. &
               any(.not. ieee_is_finite(trial_thermal))) then
             local_transaction_failure = snrt_failure_receiver
          else
             do i = 1, nleaf
                if (trial_hydrogen_ii(i) < 0.0d0 .or. trial_hydrogen_ii(i) > 1.0d0 .or. &
                     trial_helium_ii(i) < 0.0d0 .or. trial_helium_ii(i) > 1.0d0 .or. &
                     trial_helium_iii(i) < 0.0d0 .or. &
                     trial_helium_ii(i) + trial_helium_iii(i) > 1.0d0 + &
                     snrt_inventory_tolerance(1.0d0) .or. &
                     trial_neutral_hydrogen(i) < 0.0d0 .or. &
                     trial_neutral_hydrogen(i) > 1.0d0) then
                   local_transaction_failure = snrt_failure_chemistry
                   exit
                end if
             end do
          end if
       end if

       if (local_transaction_failure == snrt_failure_none) then
          relaxed_hydrogen_ii = (1.0d0-transaction_config%relaxation) * &
               current_hydrogen_ii + transaction_config%relaxation * trial_hydrogen_ii
          relaxed_helium_ii = (1.0d0-transaction_config%relaxation) * &
               current_helium_ii + transaction_config%relaxation * trial_helium_ii
          relaxed_helium_iii = (1.0d0-transaction_config%relaxation) * &
               current_helium_iii + transaction_config%relaxation * trial_helium_iii
          do i = 1, nleaf
             relaxed_hydrogen_ii(i) = min(max(relaxed_hydrogen_ii(i),0.0d0),1.0d0)
             relaxed_helium_ii(i) = min(max(relaxed_helium_ii(i),0.0d0),1.0d0)
             relaxed_helium_iii(i) = min(max(relaxed_helium_iii(i),0.0d0), &
                  max(0.0d0,1.0d0-relaxed_helium_ii(i)))
             current_fraction(i,1) = current_hydrogen_ii(i)
             current_fraction(i,2) = current_helium_ii(i)
             current_fraction(i,3) = current_helium_iii(i)
             target_fraction(i,1) = relaxed_hydrogen_ii(i)
             target_fraction(i,2) = relaxed_helium_ii(i)
             target_fraction(i,3) = relaxed_helium_iii(i)
             neutral_hydrogen_code = rho_level(i) * 0.5d0 * &
                  (start_neutral_hydrogen(i) + &
                  max(0.0d0,1.0d0-relaxed_hydrogen_ii(i)))
             neutral_helium_i_code = rho_level(i) * snrt_nhelium_per_hydrogen * 0.5d0 * &
                  (max(0.0d0,1.0d0-start_helium_ii(i)-start_helium_iii(i)) + &
                  max(0.0d0,1.0d0-relaxed_helium_ii(i)-relaxed_helium_iii(i)))
             neutral_helium_ii_code = rho_level(i) * snrt_nhelium_per_hydrogen * 0.5d0 * &
                  (start_helium_ii(i) + relaxed_helium_ii(i))
             call snrt_nlte_primordial_optical_depth_groups(neutral_hydrogen_code, &
                  neutral_helium_i_code, neutral_helium_ii_code, scale_nH, dt_s, &
                  snrt_group_cross_section_cm2, snrt_group_cross_section_hei_cm2, &
                  snrt_group_cross_section_heii_cm2, tau_dp, tau_hi_dp, tau_hei_dp, &
                  tau_heii_dp, ierr, reduced_c)
             if (ierr /= 0) then
                local_transaction_failure = snrt_failure_convergence
                exit
             end if
             target_tau(i,:) = real(max(tau_dp,0.0d0),c_float)
             target_species_tau(i,:,1) = real(max(tau_hi_dp,0.0d0),c_float)
             target_species_tau(i,:,2) = real(max(tau_hei_dp,0.0d0),c_float)
             target_species_tau(i,:,3) = real(max(tau_heii_dp,0.0d0),c_float)
             target_tau(i,:) = sum(target_species_tau(i,:,:),dim=2) + &
                  optical_depth_dust(i,:)
          end do
       end if
       if (local_transaction_failure == snrt_failure_none) then
          call snrt_transaction_check_convergence(current_fraction, target_fraction, &
               iteration_tau, target_tau, transaction_config, transaction_residual, &
               transaction_converged, convergence_status)
          if (convergence_status /= 0) then
             local_transaction_failure = snrt_failure_convergence
             transaction_converged = .false.
          end if
       end if
       local_transaction_converged = 0
       if (local_transaction_failure == snrt_failure_none .and. &
            transaction_converged) local_transaction_converged = 1
       call snrt_transaction_reduce_decision(local_transaction_failure, &
            local_transaction_converged, transaction_residual, &
            global_transaction_failure, global_transaction_converged, &
            global_transaction_residual, convergence_status)
       if (convergence_status /= 0) then
          global_transaction_failure = snrt_failure_transport
          global_transaction_converged = 0
       end if
       if (global_transaction_failure /= snrt_failure_none) exit
       if (global_transaction_converged == 1) then
          transaction_converged = .true.
          exit
       end if
       current_hydrogen_ii = relaxed_hydrogen_ii
       current_helium_ii = relaxed_helium_ii
       current_helium_iii = relaxed_helium_iii
       current_fraction(:,1) = current_hydrogen_ii
       current_fraction(:,2) = current_helium_ii
       current_fraction(:,3) = current_helium_iii
       iteration_species_tau = target_species_tau
       iteration_tau = sum(iteration_species_tau, dim=3) + optical_depth_dust
    end do

    if (global_transaction_failure /= snrt_failure_none .or. &
         .not. transaction_converged .or. global_transaction_converged == 0) then
       if (transaction_active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, transaction_status)
       if (myid == 1) then
          if (global_transaction_failure /= snrt_failure_none) then
             write(*,'(A,A,A,I0,A,I0,A,ES12.4)') &
                  ' SNRT RT transaction rollback: class=', &
                  trim(snrt_transaction_failure_name(global_transaction_failure)), &
                  ' level=', ilevel, ' iteration=', transaction_iteration, &
                  ' residual=', global_transaction_residual
          else
             write(*,'(A,I0,A,I0,A,ES12.4)') &
                  ' SNRT RT transaction non-converged: level=', ilevel, &
                  ' max_iter=', transaction_config%max_iterations, &
                  ' residual=', global_transaction_residual
          end if
       end if
       deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, optical_depth_species, &
            optical_depth_dust, available_species_transport, optical_depth_hydrogen, &
            optical_depth_helium_i, optical_depth_helium_ii, absorbed_group, raw_group, &
            absorbed_hhe_group_species, absorbed_dust_group, returned_group, &
            incoming_intensity, trial_intensity, coarse_flux_trial, iteration_tau, &
            iteration_species_tau, target_tau, target_species_tau, start_hydrogen_ii, &
            start_helium_ii, start_helium_iii, start_neutral_hydrogen, level_thermal, &
            current_hydrogen_ii, current_helium_ii, current_helium_iii, &
            relaxed_hydrogen_ii, relaxed_helium_ii, relaxed_helium_iii, &
            trial_hydrogen_ii, trial_helium_ii, trial_helium_iii, &
            trial_neutral_hydrogen, trial_thermal, rho_level, temperature_level, &
            trial_heating_rate, trial_unassigned, trial_absorbed_species, &
            current_fraction, target_fraction)
       call clean_stop
       return
    end if

    call snrt_transaction_commit_level(transaction, snrt_intensity, leaf_slot, &
         snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, snrt_neutral_fraction, &
         trial_intensity, coarse_flux_trial, trial_hydrogen_ii, trial_helium_ii, &
         trial_helium_iii, trial_neutral_hydrogen, level_thermal, trial_thermal, &
         transaction_status)
    local_transaction_failure = snrt_failure_none
    if (transaction_status /= 0) local_transaction_failure = snrt_failure_receiver
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0) global_transaction_failure = snrt_failure_receiver
    if (global_transaction_failure /= snrt_failure_none) then
       if (transaction%active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, ierr)
       if (myid == 1) write(*,'(A,A,A,I0)') &
            ' SNRT RT transaction commit failed: class=', &
            trim(snrt_transaction_failure_name(global_transaction_failure)), &
            ' level=', ilevel
       deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, optical_depth_species, &
            optical_depth_dust, available_species_transport, optical_depth_hydrogen, &
            optical_depth_helium_i, optical_depth_helium_ii, absorbed_group, raw_group, &
            absorbed_hhe_group_species, absorbed_dust_group, returned_group, &
            incoming_intensity, trial_intensity, coarse_flux_trial, iteration_tau, &
            iteration_species_tau, target_tau, target_species_tau, start_hydrogen_ii, &
            start_helium_ii, start_helium_iii, start_neutral_hydrogen, level_thermal, &
            current_hydrogen_ii, current_helium_ii, current_helium_iii, &
            relaxed_hydrogen_ii, relaxed_helium_ii, relaxed_helium_iii, &
            trial_hydrogen_ii, trial_helium_ii, trial_helium_iii, &
            trial_neutral_hydrogen, trial_thermal, rho_level, temperature_level, &
            trial_heating_rate, trial_unassigned, trial_absorbed_species, &
            current_fraction, target_fraction)
       call clean_stop
       return
    end if
    call snrt_transaction_reduce_sum(unassigned_absorption_total, &
         global_unassigned_absorption, convergence_status)
    if (convergence_status /= 0) global_unassigned_absorption = &
         unassigned_absorption_total
    do i = 1, nleaf
       icell = leaf_cell(i)
       if (icell >= 1 .and. icell <= size(uold,1) .and. &
            size(uold,2) >= energy_index) uold(icell,energy_index) = level_thermal(i)
    end do
    t_coupling = omp_get_wtime() - wall_start

    if (myid == 1) then
       write(*,'(A,I0)') ' SNRT source internals level=', ilevel
       write(*,'(A,F10.3)') '   locator  : ', t_locator
       write(*,'(A,F10.3)') '   photon   : ', t_budget
       write(*,'(A,F10.3)') '   deposit  : ', t_deposit
       write(*,'(A,F10.3)') '   overhead : ', t_source_overhead
       write(*,'(A,I0)') '   locator calls: ', n_locator_calls
       write(*,'(A,I0)') '   active sources: ', n_active_sources
       write(*,'(A,I0)') '   chemistry failures: ', chemistry_failures
       write(*,'(A,ES12.4)') '   unassigned absorption code (global): ', &
            global_unassigned_absorption
       write(*,'(A,ES12.4)') '   dust ledger relative error (last trial): ', &
            ledger_relative_error
       write(*,'(A,I0,A,I0,6(A,F10.3,1X))') &
         ' SNRT stage timings level=', ilevel, ' leaves=', nleaf, &
         ' setup=', t_setup, ' topology=', t_topology, ' nlte=', t_nlte, &
         ' source=', t_source, ' transport=', t_transport, &
         ' coupling=', t_coupling
    endif

    deallocate(leaf_cell, leaf_slot, neighbor, optical_depth, optical_depth_species, &
         optical_depth_dust, available_species_transport, &
         optical_depth_hydrogen, optical_depth_helium_i, optical_depth_helium_ii, &
         absorbed_group, raw_group, absorbed_hhe_group_species, absorbed_dust_group, &
         returned_group)
  end subroutine snrt_ramses_advance_level

end module snrt_ramses_driver

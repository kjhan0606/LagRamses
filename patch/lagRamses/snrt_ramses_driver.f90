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

  subroutine snrt_ramses_advance_level(ilevel,step_start_proper)
    use amr_parameters, only: sink, sink_AGN
    use amr_commons, only: levelmin, nstep_coarse, myid, dtnew, boxlen, &
         icoarse_min, icoarse_max, ncpu, nrestart, texp, aexp, active
    use hydro_commons, only: uold,magnetic_energy
    use dust_mass_physics, only: dust_atomic_cooling_enabled,dust_chimes_enabled,dust_gas_elements
#ifdef SNRT_CHIMES
    use snrt_chimes_runtime, only: chimes_live_capacity,chimes_live_stage,chimes_cell_state,chimes_live_band_stage, &
         chimes_live_cold_stage
    use snrt_chimes, only: chimes_ns,chimes_group_binding,chimes_boltzmann,chimes_round_subnormal_survivors
#endif
    use snrt_atomic_cooling, only: atomic_mh,atomic_temperature,atomic_heat_capacity,atomic_advance
    use pm_commons, only: nsink, xsink, idsink, agn_pending_erg, nindsink, msink, vsink, jsink, &
         dMBH_coarse, dMEd_coarse, dMsmbh, Esave, spinmag, agn_checkpoint_restored, &
         headp, numbp, nextp, ptypep, PTYPE_STAR, xp, mp0, tpp, zp
    use snrt_stellar_source, only: stellar_sed_enabled, stellar_photon_interval,stellar_sed_has_energy
    use snrt_runtime_backend, only: snrt_backend_initialize
    use snrt_state, only: snrt_ndirection, snrt_nmu, snrt_nphi, snrt_ngroups, snrt_intensity, &
         snrt_nslot, snrt_energy_shift, &
         snrt_neutral_fraction, snrt_hydrogen_ii, snrt_helium_ii, &
         snrt_helium_iii, snrt_state_get_slot
    use snrt_spectral_contract, only: &
         snrt_nedges, snrt_group_edges_ev, snrt_group_edges_sha256, snrt_band_enabled,snrt_band_model, &
         snrt_d03_band_enabled,snrt_fe_band_enabled,snrt_grain_band_bins, &
         snrt_node_secondaries_enabled,snrt_chimes_band_enabled,snrt_chimes_cold_enabled,snrt_chimes_transition_enabled, &
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
#ifdef DUST_LIVE
    use snrt_dust_live, only: snrt_dust_live_stage, snrt_dust_live_commit, dust_live_coarse_trial
    use dust_composition_material, only: dust_material_composition_enabled,dust_composition_curve,dust_composition_area
    use dust_mass_physics, only: dust_optics_enabled,dust_sublimation_rt_enabled,dust_iron_enabled,dust_fe_max_primary_ev, &
         dust_sublimation_enabled
    use dust_iron_compare, only: iron_compare_curve,iron_compare_weights,iron_compare_temperature,fe_six_opacity_basis
    use dust_mass_physics, only: dust_pah_enabled,dust_pah_nstate,dust_pah_hc,dust_pah_molecule_g,dust_pah_charged, &
         dust_pah_inventory,dust_pah_state_mass,dust_pah_hydrogenated,dust_pah_h2_enabled
    use dust_pah_live_model, only: pah_live_prepare,pah_primary_sigma,pah_primary_alpha,pah_max_primary_ev
    use dust_composition_optics, only: d03_ng,d03_nir,d03_cell_weights,d03_opacity_basis
    use snrt_dust_ir, only: dust_ir_diagnostics
    use snrt_dust_contract, only: snrt_dust_contract_version, &
         snrt_dust_contract_scattering_enabled, snrt_dust_contract_scattering_per_h_cm2, &
         snrt_dust_contract_exchange_enabled, snrt_dust_contract_collision_area_per_h, &
         snrt_dust_contract_accommodation, snrt_dust_contract_ir_background_k
    use snrt_runtime_backend, only: snrt_runtime_isotropic_scatter
    use dust_phase_state, only: dust_phase_read,dust_phase_kinetic
    use snrt_moving_scatter, only: snrt_moving_scatter_cell
    use hydro_parameters, only: dust_relative_motion,ndust_phase,idust_momentum,nvar
#endif
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
    use snrt_dust_contract, only: snrt_dust_contract_loaded, &
         snrt_dust_contract_runtime_allowed, snrt_dust_contract_number_groups, &
         snrt_dust_contract_number_temperature, snrt_dust_contract_group_edges_ev, &
         snrt_dust_contract_group_edges_sha256, &
         snrt_dust_contract_absorption_per_h_cm2, &
         snrt_dust_contract_absorption_mean_energy_ev, &
         snrt_dust_contract_temperature_k, &
         snrt_dust_contract_mass_per_h_g, &
         snrt_dust_contract_heat_capacity_per_h_erg_k, snrt_dust_contract_internal_energy_per_h_erg
    use snrt_dust_ir, only: snrt_dust_material_temperature
    use snrt_dust_receiver, only: snrt_dust_prepare_cell_optical_depth, &
         snrt_dust_receiver_stage
    use snrt_angular_quadrature, only: snrt_angular_init
    use snrt_agn_locator, only: snrt_agn_find_local_leaf
    use snrt_agn_source, only: snrt_c_cgs, snrt_ev_to_erg, snrt_agn_photon_budget_energy, &
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
    use hydro_parameters, only: gamma, idust, idust_energy, inener, idust_species, idust_bins,ichem,ichimes,idust_iron,idust_pah
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_c_binding, only: c_float
    use omp_lib, only: omp_get_wtime
    implicit none

    integer, intent(in) :: ilevel
    real(dp), intent(in), optional :: step_start_proper
    real(dp) :: stellar_photons(snrt_ngroups), stellar_time_scale
    real(dp),allocatable :: stellar_energy(:),stellar_injection_mean(:)
    integer :: igrid, ipart, ip, ig
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
    logical, save :: test_seed_checked = .false.
    integer, allocatable :: leaf_cell(:), leaf_slot(:), neighbor(:,:)
    real(dp) :: direction_dp(snrt_ndirection,3), angular_weight(snrt_ndirection)
    real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(dp) :: dt_s, dx_code, cell_volume_code, cdt_over_dx
    real(dp) :: rho_code, hydrogen_ionized_fraction, helium_ionized_fraction
    real(dp) :: helium_double_ionized_fraction
    real(dp) :: n_hydrogen_cm3, n_helium_cm3
    real(dp) :: neutral_hydrogen_code, neutral_helium_i_code, neutral_helium_ii_code
    real(dp) :: temperature_k, internal_energy, kinetic_energy, molecular_weight
#ifdef DUST_LIVE
    real(dp) :: dust_mass_code, dust_energy_code, dust_energy_scale
#endif
    real(dp) :: tau_dp(snrt_ngroups), tau_hi_dp(snrt_ngroups)
    real(dp) :: tau_hei_dp(snrt_ngroups), tau_heii_dp(snrt_ngroups)
    real(dp) :: unassigned_absorption_total, ledger_relative_error
    real(dp) :: global_unassigned_absorption
    real(dp) :: excess_energy_ev(3,snrt_ngroups),cell_excess_ev(3,snrt_ngroups)
    real(dp),parameter::band_threshold(3)=[13.60d0,24.59d0,54.42d0]
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
    real(dp),allocatable::incoming_energy_shift(:,:,:),trial_energy_shift(:,:,:),coarse_energy_shift(:,:,:)
    real(dp),allocatable::hhe_energy(:,:,:),absorbed_dust_energy_ev(:,:),dust_energy_moment(:,:,:)
    real(dp),allocatable::secondary_xi(:),band_deposition(:,:)
    real(c_float), allocatable :: iteration_tau(:,:), iteration_species_tau(:,:,:)
    real(c_float), allocatable :: target_tau(:,:), target_species_tau(:,:,:)
    real(dp), allocatable :: start_hydrogen_ii(:), start_helium_ii(:), &
         start_helium_iii(:), start_neutral_hydrogen(:), level_thermal(:)
    real(dp), allocatable :: current_hydrogen_ii(:), current_helium_ii(:), &
         current_helium_iii(:), relaxed_hydrogen_ii(:), relaxed_helium_ii(:), &
         relaxed_helium_iii(:), trial_hydrogen_ii(:), trial_helium_ii(:), &
         trial_helium_iii(:), trial_neutral_hydrogen(:), trial_thermal(:)
    real(dp), allocatable :: rho_level(:), temperature_level(:)
    real(dp), allocatable :: h_number_code(:),he_number_code(:),atomic_gas_x(:,:)
    logical :: atomic_cooling_on,band_on,paired_transport
    real(dp),allocatable::species_columns(:,:),target_columns(:,:),grain_columns(:,:)
    logical :: chimes_on
#ifdef SNRT_CHIMES
    real(dp),allocatable::chemical_trial(:,:)
    real(dp)::chemical_absorbed_ev,chemical_events,chemical_dissociation_ev
#endif
    real(dp), allocatable :: trial_heating_rate(:), trial_unassigned(:)
    real(dp), allocatable :: trial_absorbed_species(:,:,:)
    real(dp), allocatable :: current_fraction(:,:), target_fraction(:,:)
    real(dp), allocatable :: emitted_groups(:), luminosity_groups(:)
#ifdef DUST_LIVE
    real(dp), allocatable :: dust_relative_abundance(:), dust_heat_capacity(:)
    real(dp), allocatable :: dust_old_energy(:), dust_old_temperature(:)
    real(dp), allocatable :: dust_trial_energy(:), dust_trial_temperature(:)
    real(dp), allocatable :: dust_absorbed_photons(:,:), dust_absorbed_energy(:)
    real(dp),allocatable::pah_number(:,:),pah_primary_energy(:,:),dust_receiver_abundance(:)
    real(dp)::pah_solid_hc(2),pah_n
    real(dp), allocatable :: dust_n_hydrogen_cm3(:), dust_path_cm(:)
    real(dp), allocatable :: dust_tau_dp(:,:)
    real(dp), allocatable :: dust_ir_trial(:,:,:)
    real(dp), allocatable :: dust_cell_u(:,:),dust_cell_area(:)
    real(dp),allocatable::sublimation_bins(:,:),sublimation_next(:,:)
    real(dp),allocatable::dust_weights(:,:),dust_primary_sigma(:,:),dust_scatter_sigma(:,:)
    real(dp)::d03_pa(d03_ng,6),d03_ps(d03_ng,6),d03_pg(d03_ng,6)
    real(dp)::d03_ia(d03_nir,6),d03_isc(d03_nir,6),d03_ig(d03_nir,6),solid_fe
    integer::dust_nb
    real(dp),allocatable::phase_rows(:,:),phase_mass(:,:),phase_pold(:,:,:),phase_pnext(:,:,:)
    real(dp),allocatable::phase_ir_work(:,:)
    real(dp),allocatable::primary_heat(:,:),primary_pah_heat(:,:),primary_pah_captures(:,:)
    type(dust_ir_diagnostics) :: dust_ir_result
    type(dust_live_coarse_trial) :: dust_ir_coarse
#endif
    logical :: enabled, source_ok, accounting_identity_ok
    logical, allocatable :: source_transaction_ok(:)
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
#ifdef DUST_LIVE
    logical, save :: dust_contract_checked = .false.
    logical, save :: dust_contract_ok = .false.
#endif
    real(dp), save :: reduced_c = 0.01d0
    integer, save :: level_filter = -1
    character(len=20), parameter :: dust_tau_mode = &
#ifdef DUST_LIVE
         'LIVE_ABSORPTION'
#else
         'ZERO_SCAFFOLD'
#endif

    band_on=snrt_band_enabled();paired_transport=band_on
#ifdef DUST_LIVE
    paired_transport=paired_transport.or.dust_relative_motion
    if(snrt_d03_band_enabled())then
       if(.not.dust_optics_enabled().or.(dust_iron_enabled().neqv.snrt_fe_band_enabled()).or.dust_pah_enabled().or. &
            dust_sublimation_enabled().or.dust_relative_motion.or. &
            (dust_chimes_enabled().and..not.snrt_chimes_cold_enabled()))then
          if(myid==1)write(*,*)'Grain band RT requires matching static D03/Fe mode or explicit cold CHIMES; no PAH/drift'
          call clean_stop;return
       endif
    endif
#endif
    enabled = .false.
    if (.not. enabled_resolved) then
       enabled_latched = snrt_agn_rt_requested()
       if(enabled_latched.and.band_on.and.myid==1)write(*,'(A,A)')' SNRT spectral absorption model=',snrt_band_model
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
#ifdef DUST_LIVE
    if (.not. dust_contract_checked) then
       dust_contract_ok = snrt_dust_contract_loaded .and. &
            snrt_dust_contract_runtime_allowed .and. &
            snrt_dust_contract_number_groups == snrt_ngroups .and. &
            snrt_dust_contract_number_temperature >= 2 .and. &
            trim(snrt_dust_contract_group_edges_sha256) == trim(snrt_group_edges_sha256)
       if (dust_contract_ok) then
          dust_contract_ok = maxval(abs(snrt_dust_contract_group_edges_ev(1:snrt_nedges) - &
               snrt_group_edges_ev)) <= 1.0d-12 * max(1.0d0, &
               maxval(abs(snrt_group_edges_ev)))
       end if
       if(dust_sublimation_rt_enabled())dust_contract_ok=dust_contract_ok.and. &
            snrt_dust_contract_exchange_enabled.and.snrt_dust_contract_version==4
       if(snrt_d03_band_enabled())dust_contract_ok=dust_contract_ok.and. &
            snrt_dust_contract_version==4.and.snrt_dust_contract_scattering_enabled
       dust_contract_checked = .true.
       if (myid == 1 .and. .not. dust_contract_ok) then
          write(*,'(A)') ' SNRT DUST_LIVE disabled: dust contract does not match the canonical nine-group runtime'
       endif
    endif
    if (.not. dust_contract_ok) then
       call clean_stop
       return
    endif
#endif
    if (myid == 1 .and. .not. dust_scaffold_reported) then
       write(*,'(A,A)') ' SNRT dust optical-depth mode=', trim(dust_tau_mode)
#ifdef DUST_LIVE
       write(*,'(A)') '   dust mass/thermal-energy fields are staged with the coupled RT transaction'
#else
       write(*,'(A)') '   dust ledger is trial-only; no dust thermal/momentum/abundance commit'
#endif
       dust_scaffold_reported = .true.
    end if
    ! Do not guess unknown historical photon ownership by either replaying
    ! the ledger or silently rebasing away real accretion. Startup rejects
    ! these modes with sink enabled; defend restored sink arrays here too.
    if (nsink>0 .and. nrestart>0.and..not.agn_checkpoint_restored) then
       if(myid==1)write(*,*)'Live SNRT AGN requires a restored restart energy ledger'
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
             write(*,'(A,I0,A,I0,A,I0)') ' SNRT angular quadrature: mu=',snrt_nmu, &
                  ' phi=',snrt_nphi,' directions=',snrt_ndirection
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
    if (level_filter > 0 .and. ilevel /= level_filter.and..not.dust_chimes_enabled()) return
    local_transaction_failure = snrt_failure_none
    call snrt_backend_initialize(ierr)
    if (ierr /= 0) local_transaction_failure = snrt_failure_transport
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0 .or. &
         global_transaction_failure /= snrt_failure_none) then
       if (myid == 1 .and. global_transaction_failure == snrt_failure_transport) &
            write(*,'(A)') ' SNRT RT preflight failed: invalid backend controls or forced CUDA unavailable'
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

    ! Explicitly opt-in harness seed.  The normal runtime never sets this
    ! variable; it exists only to exercise the production driver with a
    ! fresh-start source when no restart payload is available.  Keep the
    ! seed in the driver so the tested source->RT call sequence is unchanged.
    if (.not. test_seed_checked) then
       test_seed_checked = .true.
       env_value = ''
       call get_environment_variable('SNRT_DRIVER_TEST_SEED_SOURCE', env_value, &
            length=env_length, status=env_status)
       if (env_status == 0 .and. trim(env_value) == '1') then
          if (.not. transaction_diagnostic_mode) then
             if (myid == 1) write(*,'(A)') &
                  ' SNRT_DRIVER_TEST_SEED_SOURCE rejected: diagnostic mode required'
          else if (nsink /= 0) then
             if (myid == 1) write(*,'(A)') &
                  ' SNRT_DRIVER_TEST_SEED_SOURCE rejected: sinks already exist'
          else if (.not. allocated(xsink) .or. .not. allocated(idsink) .or. &
               .not. allocated(agn_pending_erg) .or. .not. allocated(msink) .or. &
               .not. allocated(vsink) .or. .not. allocated(jsink) .or. &
               .not. allocated(dMBH_coarse) .or. .not. allocated(dMEd_coarse) .or. &
               .not. allocated(dMsmbh) .or. .not. allocated(Esave) .or. &
               .not. allocated(spinmag)) then
             if (myid == 1) write(*,'(A)') &
                  ' SNRT_DRIVER_TEST_SEED_SOURCE rejected: sink arrays unavailable'
          else if (size(xsink,1) < 1 .or. size(xsink,2) < ndim .or. &
               size(idsink) < 1 .or. size(agn_pending_erg) < 1 .or. &
               size(msink) < 1 .or. size(vsink,1) < 1 .or. size(vsink,2) < ndim .or. &
               size(jsink,1) < 1 .or. size(jsink,2) < ndim .or. &
               size(dMBH_coarse) < 1 .or. size(dMEd_coarse) < 1 .or. &
               size(dMsmbh) < 1 .or. size(Esave) < 1 .or. size(spinmag) < 1) then
             if (myid == 1) write(*,'(A)') &
                  ' SNRT_DRIVER_TEST_SEED_SOURCE rejected: sink arrays too small'
          else
             nsink = 1
             nindsink = 1
             idsink(1) = 1
             xsink(1,1:ndim) = 0.5d0 * boxlen
             msink(1) = 1.0d0
             vsink(1,1:ndim) = 0.0d0
             jsink(1,1:ndim) = 0.0d0
             dMBH_coarse(1) = 0.0d0
             dMEd_coarse(1) = 0.0d0
             dMsmbh(1) = 0.0d0
             Esave(1) = 0.0d0
             spinmag(1) = 0.0d0
             agn_pending_erg(1) = 1.0d-6
             if (myid == 1) write(*,'(A)') &
                  ' SNRT_DRIVER_TEST_SEED_SOURCE applied: NONPRODUCTION'
          end if
       end if
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
    if(band_on)allocate(species_columns(nleaf,3),target_columns(nleaf,3))
    if(snrt_d03_band_enabled())allocate(grain_columns(nleaf,snrt_grain_band_bins()))
    atomic_cooling_on=dust_atomic_cooling_enabled()
    chimes_on=dust_chimes_enabled()
    if((band_on.and.chimes_on.and..not.snrt_chimes_band_enabled()).or. &
         (snrt_chimes_band_enabled().and..not.chimes_on))then
       if(myid==1)write(*,*)'Band H/He closure has no matched CHIMES coefficients'
       call clean_stop;return
    endif
#ifdef SNRT_CHIMES
    if(chimes_on)allocate(chemical_trial(chimes_ns,nleaf))
#endif
    allocate(h_number_code(nleaf),he_number_code(nleaf))
    if(atomic_cooling_on.or.chimes_on)then
       allocate(atomic_gas_x(11,nleaf));atomic_gas_x=0
    endif
    hydro_state_invalid = .false.
#ifdef SNRT_CHIMES
    if(chimes_on.and..not.snrt_chimes_band_enabled())then
       ierr=chimes_group_binding(snrt_ngroups,snrt_group_edges_ev,snrt_group_mean_energy_ev)
       if(ierr/=0)hydro_state_invalid=.true.
    endif
#endif
#ifdef DUST_LIVE
    allocate(dust_relative_abundance(nleaf), dust_heat_capacity(nleaf), &
         dust_old_energy(nleaf), dust_old_temperature(nleaf), &
         dust_trial_energy(nleaf), dust_trial_temperature(nleaf), &
         dust_absorbed_photons(snrt_ngroups,nleaf), dust_absorbed_energy(nleaf), &
         dust_n_hydrogen_cm3(nleaf), dust_path_cm(nleaf), &
         dust_tau_dp(nleaf,snrt_ngroups))
    if(dust_material_composition_enabled())then
       allocate(dust_cell_u(snrt_dust_contract_number_temperature,nleaf),dust_cell_area(nleaf))
       dust_cell_u=0;dust_cell_area=0
    endif
    if(dust_pah_enabled())then
       call pah_live_prepare(ierr)
       if(ierr/=0)hydro_state_invalid=.true.
       allocate(pah_number(dust_pah_nstate(),nleaf),pah_primary_energy(snrt_ngroups,nleaf))
       do i=1,nleaf
          block
            integer::k
            do k=1,dust_pah_nstate()
               pah_number(k,i)=uold(leaf_cell(i),idust_pah+k-1)*scale_d/dust_pah_state_mass(k)
            enddo
          end block
       enddo
       pah_primary_energy=0
    endif
    if(dust_sublimation_rt_enabled())then
       allocate(sublimation_bins(4,nleaf),sublimation_next(4,nleaf))
       do i=1,nleaf
          sublimation_bins(:,i)=uold(leaf_cell(i),idust_bins:idust_bins+3)*scale_d
       enddo
    endif
    if(dust_optics_enabled())then
       dust_nb=merge(6,4,dust_iron_enabled())
       allocate(dust_weights(dust_nb,nleaf),dust_primary_sigma(snrt_ngroups,nleaf),dust_scatter_sigma(snrt_ngroups,nleaf))
       dust_weights=.25d0;dust_primary_sigma=0;dust_scatter_sigma=0
       if(dust_iron_enabled())then
          call fe_six_opacity_basis(snrt_dust_contract_mass_per_h_g,d03_pa,d03_ps,d03_pg,d03_ia,d03_isc,d03_ig,ierr)
       else
          call d03_opacity_basis(snrt_dust_contract_mass_per_h_g,d03_pa(:,1:4),d03_ps(:,1:4),d03_pg(:,1:4), &
               d03_ia(:,1:4),d03_isc(:,1:4),d03_ig(:,1:4),ierr)
       endif
       if(ierr/=0.or.snrt_ngroups/=d03_ng)hydro_state_invalid=.true.
    endif
#endif
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
#ifdef DUST_LIVE
    dust_energy_scale = scale_d * scale_v**2
    if (.not. ieee_is_finite(dust_energy_scale) .or. dust_energy_scale <= 0.0d0) &
         hydro_state_invalid = .true.
#endif
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
       h_number_code(i)=rho_code
       he_number_code(i)=rho_code*snrt_nhelium_per_hydrogen
       if((atomic_cooling_on.or.chimes_on).and.rho_code>0)then
          call dust_gas_elements(uold(icell,ichem:ichem+10), &
               uold(icell,idust_species:idust_species+1),atomic_gas_x(:,i),ierr)
#ifdef DUST_LIVE
          if(dust_iron_enabled())call dust_gas_elements(uold(icell,ichem:ichem+10), &
               uold(icell,idust_species:idust_species+1),atomic_gas_x(:,i),ierr, &
               sum(uold(icell,idust_iron:idust_iron+1)))
#endif
          if(ierr/=0)hydro_state_invalid=.true.
#ifdef DUST_LIVE
          if(dust_pah_enabled())then
             pah_solid_hc=dust_pah_inventory(uold(icell,idust_pah:idust_pah+dust_pah_nstate()-1))
             solid_fe=0
             if(dust_iron_enabled())solid_fe=sum(uold(icell,idust_iron:idust_iron+1))
             call dust_gas_elements(uold(icell,ichem:ichem+10),uold(icell,idust_species:idust_species+1), &
                  atomic_gas_x(:,i),ierr,solid_fe,pah_solid_hc)
             if(ierr/=0)hydro_state_invalid=.true.
          endif
#endif
          atomic_gas_x(:,i)=atomic_gas_x(:,i)/rho_code
          h_number_code(i)=rho_code*scale_d*atomic_gas_x(1,i)/(atomic_mh*scale_nH)
          he_number_code(i)=rho_code*scale_d*atomic_gas_x(2,i)/(4*atomic_mh*scale_nH)
       endif
#ifdef DUST_LIVE
       dust_n_hydrogen_cm3(i) = h_number_code(i) * scale_nH
       dust_path_cm(i) = dx_code * scale_l
       dust_relative_abundance(i) = 0.0d0
       dust_heat_capacity(i) = 1.0d0
       dust_old_energy(i) = 0.0d0
       dust_old_temperature(i) = max(1.0d0, snrt_dust_contract_temperature_k(1))
       dust_mass_code = -1.0d0
       dust_energy_code = -1.0d0
       if (icell >= 1 .and. icell <= size(uold,1) .and. &
            size(uold,2) >= idust_energy) then
          dust_mass_code = uold(icell,idust)
          dust_energy_code = uold(icell,idust_energy)
          solid_fe=0
          if(dust_iron_enabled())solid_fe=sum(uold(icell,idust_iron:idust_iron+1))
          if(dust_optics_enabled())then
             if(dust_iron_enabled())then
                call iron_compare_weights([uold(icell,idust_bins:idust_bins+3),uold(icell,idust_iron:idust_iron+1)], &
                     dust_weights(:,i),dust_cell_area(i),snrt_dust_contract_mass_per_h_g,ierr)
             else
                call d03_cell_weights(uold(icell,idust_bins:idust_bins+3),dust_weights(:,i),ierr)
             endif
             if(ierr/=0)hydro_state_invalid=.true.
             dust_primary_sigma(:,i)=matmul(d03_pa(:,1:dust_nb),dust_weights(:,i))
             dust_scatter_sigma(:,i)=matmul(d03_ps(:,1:dust_nb)-d03_pg(:,1:dust_nb),dust_weights(:,i))
          endif
          if(dust_iron_enabled())then
             call iron_compare_curve(snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature), &
                  uold(icell,idust_species:idust_species+1),solid_fe, &
                  snrt_dust_contract_mass_per_h_g,dust_cell_u(:,i),ierr)
             if(ierr/=0)hydro_state_invalid=.true.
          else if(dust_material_composition_enabled())then
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature), &
                  uold(icell,idust_species:idust_species+1),snrt_dust_contract_mass_per_h_g,dust_cell_u(:,i),ierr)
             if(ierr/=0)hydro_state_invalid=.true.
             call dust_composition_area(uold(icell,idust_bins:idust_bins+3), &
                  snrt_dust_contract_mass_per_h_g,dust_cell_area(i),ierr)
             if(ierr/=0)hydro_state_invalid=.true.
          endif
          if (.not. ieee_is_finite(dust_mass_code) .or. &
               .not. ieee_is_finite(dust_energy_code) .or. &
               dust_mass_code < 0.0d0 .or. dust_energy_code < 0.0d0) then
             hydro_state_invalid = .true.
          else if (dust_mass_code > 0.0d0) then
             if (dust_n_hydrogen_cm3(i) <= 0.0d0 .or. &
                  dust_energy_code <= 0.0d0) then
                hydro_state_invalid = .true.
             else
                dust_relative_abundance(i) = dust_mass_code * scale_d / &
                     (dust_n_hydrogen_cm3(i) * snrt_dust_contract_mass_per_h_g)
                dust_heat_capacity(i) = dust_n_hydrogen_cm3(i) * &
                     dust_relative_abundance(i) * &
                     snrt_dust_contract_heat_capacity_per_h_erg_k
                dust_old_energy(i) = dust_energy_code * dust_energy_scale
                if(snrt_dust_contract_version==4)then
                   ! Legacy argument retained for ABI; not a physical capacity
                   ! in v4. The IR solver uses density*U(T) instead.
                   dust_heat_capacity(i)=1d0
                   if(dust_iron_enabled().or.dust_pah_enabled())then
                      call iron_compare_temperature(uold(icell,idust_species:idust_species+1),solid_fe, &
                           dust_energy_code*scale_v**2,dust_old_temperature(i),ierr)
                   else if(dust_material_composition_enabled())then
                      call snrt_dust_material_temperature( &
                           snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature), &
                           dust_cell_u(:,i),dust_old_energy(i)/(dust_n_hydrogen_cm3(i)*dust_relative_abundance(i)), &
                           dust_old_temperature(i),ierr)
                   else
                   call snrt_dust_material_temperature( &
                        snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature), &
                        snrt_dust_contract_internal_energy_per_h_erg(1:snrt_dust_contract_number_temperature), &
                        dust_old_energy(i)/(dust_n_hydrogen_cm3(i)*dust_relative_abundance(i)), &
                        dust_old_temperature(i),ierr)
                   endif
                   if(ierr/=0)hydro_state_invalid=.true.
                else
                   dust_old_temperature(i) = dust_old_energy(i) / dust_heat_capacity(i)
                endif
                if (.not. ieee_is_finite(dust_relative_abundance(i)) .or. &
                     .not. ieee_is_finite(dust_heat_capacity(i)) .or. &
                     .not. ieee_is_finite(dust_old_energy(i)) .or. &
                     .not. ieee_is_finite(dust_old_temperature(i)) .or. &
                     dust_relative_abundance(i) <= 0.0d0 .or. &
                     dust_heat_capacity(i) <= 0.0d0 .or. &
                     dust_old_energy(i) <= 0.0d0 .or. &
                     dust_old_temperature(i) <= 0.0d0) hydro_state_invalid = .true.
             end if
          else if (dust_energy_code > 0.0d0) then
             ! A thermal dust energy without a dust mass is not a valid state.
             hydro_state_invalid = .true.
          end if
       else
          hydro_state_invalid = .true.
       end if
#endif
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
#ifdef DUST_LIVE
             if(dust_relative_motion)kinetic_energy=dust_phase_kinetic(uold(icell,:))
#endif
             internal_energy = level_thermal(i) - kinetic_energy-magnetic_energy(uold(icell,:))
#if NENER>0
             internal_energy=internal_energy-sum(uold(icell,inener:inener+NENER-1))
#endif
             molecular_weight = snrt_mean_molecular_weight(hydrogen_ionized_fraction, &
                  helium_ionized_fraction, helium_double_ionized_fraction)
             if (ieee_is_finite(internal_energy) .and. internal_energy > 0.0d0) &
                  temperature_level(i) = max(1.0d0,(gamma-1.0d0)*internal_energy/rho_code * &
                  scale_T2*molecular_weight)
             if(atomic_cooling_on.and.internal_energy>0)temperature_level(i)= &
                  atomic_temperature(rho_code*scale_d,atomic_gas_x(:,i), &
                  [hydrogen_ionized_fraction,helium_ionized_fraction,helium_double_ionized_fraction], &
                  gamma,internal_energy*scale_d*scale_v**2)
#ifdef SNRT_CHIMES
             if(chimes_on.and.internal_energy>0)temperature_level(i)=internal_energy*scale_d*scale_v**2/ &
                  chimes_live_capacity(uold(icell,ichimes:ichimes+chimes_ns-1),scale_d)
#endif
          end if
       end if
       neutral_hydrogen_code = h_number_code(i) * (1.0d0-hydrogen_ionized_fraction)
       neutral_helium_i_code = he_number_code(i) * &
            (1.0d0-helium_ionized_fraction-helium_double_ionized_fraction)
       neutral_helium_ii_code = he_number_code(i) * &
            helium_ionized_fraction
       if(band_on)species_columns(i,:)=[neutral_hydrogen_code,neutral_helium_i_code,neutral_helium_ii_code]* &
            (scale_nH*snrt_c_cgs*reduced_c*dt_s)
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
    if(chimes_on)then
       ! The full native receiver owns ALL atomic and molecular absorption.
       ! Transport still handles dust absorption/scattering; never apply a
       ! second H/He photon sink or the legacy secondary/recombination step.
       optical_depth=0;optical_depth_species=0
       optical_depth_hydrogen=0;optical_depth_helium_i=0;optical_depth_helium_ii=0
       if(band_on)species_columns=0
    endif
#ifdef DUST_LIVE
    call snrt_dust_prepare_cell_optical_depth(dust_n_hydrogen_cm3, dust_path_cm, &
         dust_relative_abundance, snrt_dust_contract_absorption_per_h_cm2(1:snrt_ngroups), &
         dust_tau_dp, ierr)
    if(dust_optics_enabled())then
       do i=1,nleaf
          dust_tau_dp(i,:)=dust_n_hydrogen_cm3(i)*dust_path_cm(i)*dust_relative_abundance(i)*dust_primary_sigma(:,i)
       enddo
    endif
    dust_receiver_abundance=dust_relative_abundance
    if(dust_pah_enabled().and.allocated(pah_primary_sigma))then
       do i=1,nleaf
          pah_n=sum(pah_number(:,i))
          ! Sum physical opacities before the ONE primary photon debit.
          ! A PAH-bearing cell need not contain bulk grains. Bulk opacity,
          ! collision area and IR mass normalization remain bulk-only.
          dust_tau_dp(i,:)=dust_tau_dp(i,:)+dust_path_cm(i)*pah_primary_alpha(pah_number(:,i))
          ! This deferred receiver checks only absence/presence, never uses
          ! this value as a material mass, opacity or heat capacity.
          if(pah_n>0) dust_receiver_abundance(i)=max(dust_receiver_abundance(i),1d0)
       enddo
    endif
    if (ierr /= 0 .or. any(.not. ieee_is_finite(dust_tau_dp)) .or. &
         any(dust_tau_dp < 0.0d0)) then
       hydro_state_invalid = .true.
       dust_tau_dp = 0.0d0
    end if
    optical_depth_dust = real(dust_tau_dp, c_float)
    if(allocated(grain_columns))then
       do i=1,nleaf
          grain_columns(i,1:4)=uold(leaf_cell(i),idust_bins:idust_bins+3)*scale_d*dust_path_cm(i)
          if(snrt_fe_band_enabled()) &
               grain_columns(i,5:6)=uold(leaf_cell(i),idust_iron:idust_iron+1)*scale_d*dust_path_cm(i)
       enddo
       ! The native node operator owns grain extinction. Legacy group tau
       ! remains exactly zero and must not supply a second opacity sink.
       optical_depth_dust=0
    endif
#endif
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
#ifdef DUST_LIVE
       deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
            dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
            dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
            dust_path_cm, dust_tau_dp)
#endif
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

    ! Validate the accepted-event ledger before entering the transaction
    ! window.  Invalid energy is terminal input corruption, not a coupled
    ! RT failure; checking it here ensures that no source has been deposited
    ! and no rollback/deallocation branch can be bypassed.
    if (accounting_identity_ok .and. nsink > 0 .and. allocated(agn_pending_erg)) then
       if (size(agn_pending_erg) >= nsink) then
          do isink = 1, nsink
             if (.not. ieee_is_finite(agn_pending_erg(isink)) .or. &
                  agn_pending_erg(isink) < 0.0d0) then
                if (myid == 1) write(*,'(A,I0)') &
                     'Invalid accepted AGN radiative energy for sink ', idsink(isink)
                call clean_stop
                return
             end if
          end do
       end if
    end if

    ! Begin the enclosing transaction before any AGN source photon is added.
    ! The pending-energy marker is deliberately not consumed in the source
    ! loop: it is cleared only after RT/chemistry/dust has passed its global
    ! commit.  Before terminal failure, this preserves the accepted event in
    ! memory without duplicating or losing it; durable restart retry is not
    ! provided by this routine.
    call snrt_transaction_begin(transaction, snrt_intensity, leaf_slot, &
         snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
         snrt_neutral_fraction, level_thermal, transaction_status, &
         persistent_energy_shift=snrt_energy_shift,group_mean_energy_ev=snrt_group_mean_energy_ev)
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
            snrt_neutral_fraction, level_thermal, transaction_status,persistent_energy_shift=snrt_energy_shift)
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
#ifdef DUST_LIVE
       deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
            dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
            dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
            dust_path_cm, dust_tau_dp)
#endif
       call clean_stop
       return
    end if
    allocate(source_transaction_ok(max(1,nsink)))
    source_transaction_ok = .false.

    ! Proper-time endpoints belong to this recursive level step, not to the
    ! feedback indtab marker (which was already advanced). Birth mass/epoch
    ! are carried by the native stellar particle and HDF5 particle payload.
    if(stellar_sed_enabled)then
       if(stellar_sed_has_energy)allocate(stellar_energy(snrt_ngroups),stellar_injection_mean(snrt_ngroups))
       local_transaction_failure=snrt_failure_none
       if(.not.present(step_start_proper))then
          local_transaction_failure=snrt_failure_transport
       else
          stellar_time_scale=scale_t/aexp**2/31557600d6
          do ig=1,active(ilevel)%ngrid
             igrid=active(ilevel)%igrid(ig)
             ipart=headp(igrid)
             do ip=1,numbp(igrid)
                if(ptypep(ipart)==PTYPE_STAR)then
                   call stellar_photon_interval((step_start_proper-tpp(ipart))*stellar_time_scale, &
                        (texp-tpp(ipart))*stellar_time_scale,zp(ipart), &
                        mp0(ipart)*scale_d*scale_l**3/1.98847d33,stellar_photons,ierr,energy_ev=stellar_energy)
                   if(ierr/=0)then
                      write(*,*)'SNRT stellar interval rejected: rank, code, Z, initial mass=',myid,ierr,zp(ipart),mp0(ipart)
                      local_transaction_failure=snrt_failure_transport
                   else if(any(stellar_photons>0d0))then
                      call snrt_agn_find_local_leaf(xp(ipart,1:ndim),icell,ilevel_found,igrid,ilevel)
                      islot=0
                      if(icell>0.and.ilevel_found==ilevel)islot=snrt_state_get_slot(icell)
                      if(islot<=0)then
                         local_transaction_failure=snrt_failure_transport
                      else
                         if(paired_transport)then
                         if(stellar_sed_has_energy)then
                            stellar_injection_mean=snrt_group_mean_energy_ev
                            where(stellar_photons>0d0)stellar_injection_mean=stellar_energy/stellar_photons
                         endif
                         call snrt_agn_deposit_transaction(snrt_intensity,islot,stellar_photons, &
                              cell_volume_code,scale_l,scale_nH,angular_weight,deposited_density,ierr, &
                              persistent_energy_shift=snrt_energy_shift,group_mean_energy_ev=snrt_group_mean_energy_ev, &
                              quantize_source=band_on,source_mean_energy_ev=stellar_injection_mean)
                         else
                         call snrt_agn_deposit_transaction(snrt_intensity,islot,stellar_photons, &
                              cell_volume_code,scale_l,scale_nH,angular_weight,deposited_density,ierr)
                         endif
                         if(ierr/=0)then
                            write(*,*)'SNRT stellar deposition rejected: rank, code, min/max photons=', &
                                 myid,ierr,minval(stellar_photons),maxval(stellar_photons)
                            local_transaction_failure=snrt_failure_transport
                         endif
                         if(ierr==0)n_active_sources=n_active_sources+1
                      endif
                   endif
                endif
                ipart=nextp(ipart)
             enddo
          enddo
       endif
       call snrt_transaction_reduce_decision(local_transaction_failure,0,0d0, &
            global_transaction_failure,global_transaction_converged,global_transaction_residual,convergence_status)
       if(global_transaction_failure/=snrt_failure_none.or.convergence_status/=0)then
          call snrt_transaction_restore(transaction,snrt_intensity,leaf_slot,snrt_hydrogen_ii, &
               snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction,level_thermal,transaction_status, &
               persistent_energy_shift=snrt_energy_shift)
          if(myid==1)write(*,*)'SNRT stellar source rejected: age/Z coverage, birth mass or local ownership'
          call clean_stop
          return
       endif
    endif

    local_transaction_failure=snrt_failure_none
    if (accounting_identity_ok .and. nsink > 0) then
       if (allocated(agn_pending_erg) .and. allocated(xsink)) then
          if (size(agn_pending_erg)>=nsink .and. size(xsink,1)>=nsink) then
             allocate(emitted_groups(snrt_ngroups), luminosity_groups(snrt_ngroups))
             do isink = 1, nsink
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
                   if(paired_transport)then
                   call snrt_agn_deposit_transaction(snrt_intensity,islot,emitted_groups, &
                        cell_volume_code,scale_l,scale_nH,angular_weight,deposited_density,ierr, &
                        persistent_energy_shift=snrt_energy_shift,group_mean_energy_ev=snrt_group_mean_energy_ev, &
                        quantize_source=band_on)
                   else
                   call snrt_agn_deposit_transaction(snrt_intensity, islot, &
                        emitted_groups, cell_volume_code, scale_l, scale_nH, &
                        angular_weight, deposited_density, ierr)
                   endif
                   t_deposit = t_deposit + omp_get_wtime() - wall_sub
                   if (ierr /= 0) source_ok = .false.
                end if
                ! Defer fuel consumption until the enclosing RT/chemistry/dust
                ! transaction has passed its global commit.  A failed source
                ! deposition, or a later coupled rollback, leaves this event
                ! pending in memory until a successful coupled commit.
                if (source_ok) source_transaction_ok(isink) = .true.
                if (.not.source_ok) local_transaction_failure=snrt_failure_transport
             end do
             deallocate(emitted_groups, luminosity_groups)
          end if
       end if
    end if
    ! All ranks carry the same accepted sink receipts. Exactly one spatial
    ! owner injects each event; publish its success before the global coupled
    ! solve, then clear the replicated receipt everywhere only after commit.
    block
#ifndef WITHOUTMPI
      use mpi_mod, only: MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD
#endif
      integer :: owned(max(1,nsink)),owners(max(1,nsink)),info
      owned=merge(1,0,source_transaction_ok)
      owners=owned
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(owned,owners,max(1,nsink),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
      if(info/=0)local_transaction_failure=snrt_failure_transport
#endif
      if(any(owners>1))local_transaction_failure=snrt_failure_transport
      source_transaction_ok=owners==1
    end block
    call snrt_transaction_reduce_decision(local_transaction_failure,0,0d0, &
         global_transaction_failure,global_transaction_converged,global_transaction_residual,convergence_status)
    if(global_transaction_failure/=snrt_failure_none.or.convergence_status/=0)then
       call snrt_transaction_restore(transaction,snrt_intensity,leaf_slot,snrt_hydrogen_ii, &
            snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction,level_thermal,transaction_status, &
            persistent_energy_shift=snrt_energy_shift)
       if(myid==1)write(*,*)'SNRT AGN source failed: invalid budget/deposition or duplicate MPI owner'
       call clean_stop
       return
    endif
    t_source = omp_get_wtime() - wall_start
    t_source_overhead = t_source - t_locator - t_budget - t_deposit

    ! The enclosing RT/chemistry/dust transaction now consumes the source
    ! photons staged above.  Its snapshot predates source injection, so every
    ! coupled rollback restores both the pre-source RT state and the pending
    ! AGN event marker.
    allocate(incoming_energy_shift(snrt_ndirection,snrt_ngroups,nleaf), &
         trial_energy_shift(snrt_ndirection,snrt_ngroups,nleaf), &
         coarse_energy_shift(snrt_ndirection,snrt_ngroups,size(snrt_intensity,3)))
    trial_energy_shift=0;coarse_energy_shift=0
    if(paired_transport)allocate(hhe_energy(nleaf,snrt_ngroups,3), &
         absorbed_dust_energy_ev(nleaf,snrt_ngroups),dust_energy_moment(nleaf,snrt_ngroups,3))
    if(snrt_node_secondaries_enabled())allocate(secondary_xi(nleaf),band_deposition(nleaf,8))
#ifdef DUST_LIVE
    if(dust_relative_motion)then
       allocate(phase_rows(nvar,nleaf),phase_mass(ndust_phase,nleaf), &
            phase_pold(3,ndust_phase,nleaf),phase_pnext(3,ndust_phase,nleaf),primary_heat(snrt_ngroups,nleaf), &
            phase_ir_work(ndust_phase,nleaf))
       phase_ir_work=0
       if(dust_pah_enabled())allocate(primary_pah_heat(snrt_ngroups,nleaf),primary_pah_captures(snrt_ngroups,nleaf))
       do i=1,nleaf
          block
            real(dp)::gas_mass,gas_momentum(3),kinetic
            phase_rows(:,i)=uold(leaf_cell(i),:)
            call dust_phase_read(phase_rows(:,i),phase_mass(:,i),phase_pold(:,:,i), &
                 gas_mass,gas_momentum,kinetic,ierr)
            if(ierr/=0)then
               call clean_stop;return
            endif
          end block
       enddo
       phase_mass=phase_mass*scale_d;phase_pold=phase_pold*(scale_d*scale_v)
    endif
#endif
    do i = 1, nleaf
       do igroup = 1, snrt_ngroups
          incoming_intensity(:,igroup,i) = snrt_intensity(:,igroup,leaf_slot(i))
          incoming_energy_shift(:,igroup,i)=snrt_energy_shift(:,igroup,leaf_slot(i))
       end do
    end do

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
#ifdef DUST_LIVE
       dust_trial_energy = dust_old_energy
       dust_trial_temperature = dust_old_temperature
#endif
       wall_sub = omp_get_wtime()
       if(allocated(secondary_xi))secondary_xi=start_hydrogen_ii
       if(paired_transport)then
          call snrt_transport_absorb_multigroup_prepared_dust_trial(leaf_slot, neighbor, &
               cdt_over_dx, iteration_tau, iteration_species_tau, optical_depth_dust, &
               available_species_transport, incoming_intensity, trial_intensity, &
               coarse_flux_trial, raw_group, absorbed_hhe_group_species, &
               absorbed_dust_group, returned_group, absorbed_group, ierr, leaf_cell, ilevel, &
               incoming_shift=incoming_energy_shift,trial_shift=trial_energy_shift,coarse_shift=coarse_energy_shift, &
               absorbed_hhe_energy=hhe_energy,absorbed_dust_energy=absorbed_dust_energy_ev, &
               dust_energy_moment=dust_energy_moment,species_columns=species_columns, &
               secondary_xi=secondary_xi,band_deposition=band_deposition,grain_columns=grain_columns)
       else
       call snrt_transport_absorb_multigroup_prepared_dust_trial(leaf_slot, neighbor, &
            cdt_over_dx, iteration_tau, iteration_species_tau, optical_depth_dust, &
            available_species_transport, incoming_intensity, trial_intensity, &
            coarse_flux_trial, raw_group, absorbed_hhe_group_species, &
            absorbed_dust_group, returned_group, absorbed_group, ierr, leaf_cell, ilevel)
       endif
       if(ierr/=0)write(*,*)'SNRT transport failure: rank=',myid,' level=',ilevel,' code=',ierr
#ifdef DUST_LIVE
       ! Lie split after transport/absorption, rebuilt from the same incoming
       ! state on every nonlinear trial. Scatter only owned leaves; each group
       ! conserves photons/energy locally and adds no absorption/heating ledger.
       if(ierr==0.and.dust_relative_motion)then
          phase_pnext=phase_pold;primary_heat=0
          if(allocated(primary_pah_heat))then
             primary_pah_heat=0;primary_pah_captures=0
          endif
          do i=1,nleaf
             block
               real(dp)::sigma(snrt_ngroups,ndust_phase),fractions(snrt_ngroups,ndust_phase)
               real(dp)::impulse(3,snrt_ngroups),kick(3),velocity(3),heat(snrt_ngroups),work(ndust_phase)
               real(dp)::tau(ndust_phase,snrt_ngroups),number(snrt_ndirection,snrt_ngroups)
               real(dp)::radiation_energy(snrt_ndirection,snrt_ngroups),total_sigma,ke_change
               integer::b,g,k
               sigma=0;fractions=0;tau=0
               do b=1,dust_nb
                  sigma(:,b)=d03_pa(:,b)*phase_mass(b,i)/snrt_dust_contract_mass_per_h_g
                  tau(b,:)=(d03_ps(:,b)-d03_pg(:,b))*phase_mass(b,i)/ &
                       snrt_dust_contract_mass_per_h_g*snrt_c_cgs*reduced_c*dt_s
               enddo
               if(dust_pah_enabled())sigma(:,ndust_phase)=sum(pah_number(:,i))*pah_primary_sigma
               do g=1,snrt_ngroups
                  total_sigma=sum(sigma(g,:))
                  if(total_sigma>0)fractions(g,:)=sigma(g,:)/total_sigma
                  if(total_sigma==0.and.absorbed_dust_energy_ev(i,g)>0)ierr=1
               enddo
               if(ierr/=0)exit
               ke_change=0
               do b=1,ndust_phase
                  if(phase_mass(b,i)==0)cycle
                  do g=1,snrt_ngroups
                     impulse(:,g)=dust_energy_moment(i,g,:)*fractions(g,b)*scale_nH*snrt_ev_to_erg/snrt_c_cgs
                  enddo
                  kick=sum(impulse,dim=2)
                  velocity=(phase_pold(:,b,i)+.5d0*kick)/phase_mass(b,i)
                  if(sqrt(sum(velocity**2))>.01d0*snrt_c_cgs)then
                     ierr=1;exit
                  endif
                  heat=absorbed_dust_energy_ev(i,:)*fractions(:,b)*scale_nH*snrt_ev_to_erg- &
                       matmul(velocity,impulse)
                  if(any(.not.ieee_is_finite(heat)).or.any(heat<0))then
                     ierr=1;exit
                  endif
                  phase_pnext(:,b,i)=phase_pold(:,b,i)+kick
                  ke_change=ke_change+dot_product(velocity,kick)
                  primary_heat(:,i)=primary_heat(:,i)+heat
                  if(dust_pah_enabled().and.b==ndust_phase)then
                     primary_pah_heat(:,i)=heat
                     primary_pah_captures(:,i)=real(absorbed_dust_group(i,:),dp)*scale_nH*fractions(:,b)
                  endif
               enddo
               if(ierr/=0)exit
               if(snrt_dust_contract_scattering_enabled)then
                  number=real(trial_intensity(:,:,i),dp)*scale_nH
                  do g=1,snrt_ngroups
                     radiation_energy(:,g)=(real(trial_intensity(:,g,i),dp)*snrt_group_mean_energy_ev(g)+ &
                          trial_energy_shift(:,g,i))*scale_nH*snrt_ev_to_erg
                  enddo
                  work=0
                  call snrt_moving_scatter_cell(number,radiation_energy,phase_pnext(:,:,i), &
                       phase_mass(:,i),tau,transpose(direction_dp),angular_weight,snrt_c_cgs,work,ierr)
                  if(ierr/=0)exit
                  trial_intensity(:,:,i)=real(number/scale_nH,c_float)
                  do g=1,snrt_ngroups
                     trial_energy_shift(:,g,i)=radiation_energy(:,g)/(scale_nH*snrt_ev_to_erg)- &
                          snrt_group_mean_energy_ev(g)*real(trial_intensity(:,g,i),dp)
                  enddo
                  ke_change=ke_change+sum(work)
               endif
               phase_rows(:,i)=uold(leaf_cell(i),:)
               phase_rows(2:4,i)=phase_rows(2:4,i)+sum(phase_pnext(:,:,i)-phase_pold(:,:,i),dim=2)/(scale_d*scale_v)
               do b=1,ndust_phase
                  k=idust_momentum+3*(b-1);phase_rows(k:k+2,i)=phase_pnext(:,b,i)/(scale_d*scale_v)
               enddo
               phase_rows(5,i)=phase_rows(5,i)+ke_change/(scale_d*scale_v**2)
               trial_thermal(i)=trial_thermal(i)+ke_change/(scale_d*scale_v**2)
             end block
          enddo
       endif
       if (ierr==0.and.snrt_dust_contract_scattering_enabled.and..not.dust_relative_motion.and. &
            .not.snrt_d03_band_enabled()) &
            call snrt_runtime_isotropic_scatter(trial_intensity,angular_weight, &
            dust_n_hydrogen_cm3*dust_relative_abundance, &
            snrt_dust_contract_scattering_per_h_cm2(1:snrt_ngroups), &
            snrt_c_cgs*reduced_c*dt_s,ierr,dust_scatter_sigma)
#endif
       t_transport = t_transport + omp_get_wtime() - wall_sub
#ifdef DUST_LIVE
       if(dust_pah_enabled().and.ierr==0)then
          do i=1,nleaf
             if(sum(pah_number(:,i))<=0)cycle
             do igroup=1,snrt_ngroups
                if(snrt_dust_contract_absorption_mean_energy_ev(igroup)<=pah_max_primary_ev)cycle
                ! A zero unsupported PAH coefficient is a rejection mask,
                ! NOT permission to propagate hard light through PAH-only
                ! cells as if those grains were physically transparent.
                if(any(trial_intensity(:,igroup,i)>0).or.absorbed_dust_group(i,igroup)>0)ierr=1
             enddo
          enddo
       endif
#endif
       if (ierr /= 0) then
          local_transaction_failure = snrt_failure_transport
          write(*,*)'SNRT transport/scatter failure: rank=',myid,' code=',ierr
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
#ifdef DUST_LIVE
       ! Transport's absorption ledger was validated above. In the cold
       ! CHIMES mode it is zero: gas and dust now compete before the material
       ! stage, always starting chemistry from the unchanged incoming cell.
#ifdef SNRT_CHIMES
       if(snrt_chimes_cold_enabled().and.local_transaction_failure==snrt_failure_none)then
          chemical_absorbed_ev=0
          chemical_events=0;chemical_dissociation_ev=0
!$omp parallel do default(shared) private(i,icell,ierr,igroup) reduction(+:chemical_absorbed_ev) &
!$omp reduction(max:local_transaction_failure) reduction(+:chemical_events,chemical_dissociation_ev)
          do i=1,nleaf
             block
               real(dp)::rn(snrt_ndirection,9),re(snrt_ndirection,9),nn(snrt_ndirection,9),ne(snrt_ndirection,9)
               real(dp)::ledger(11),gn(9),ge(9),nh_code,he_code,events(2)
               icell=leaf_cell(i)
               rn=real(trial_intensity(:,:,i),dp)*scale_nH
               do igroup=1,9
                  re(:,igroup)=(snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)+ &
                       trial_energy_shift(:,igroup,i))*scale_nH
               enddo
               call chimes_live_cold_stage(icell,scale_d,scale_v,dt_s,dx_code*scale_l,dust_old_temperature(i), &
                    reduced_c,snrt_ndirection,rn,re,chemical_trial(:,i),trial_thermal(i),nn,ne,ledger,gn,ge,ierr,events)
               if(ierr==0)ierr=chimes_round_subnormal_survivors(scale_nH,sum(re),nn,ne)
               if(ierr/=0)then
!$omp critical(chimes_cold_failure)
                  write(*,*)'CHIMES cold cell rejection: rank/cell/status=',myid,icell,ierr
!$omp end critical(chimes_cold_failure)
                  local_transaction_failure=snrt_failure_chemistry;cycle
               endif
               trial_intensity(:,:,i)=real(nn/scale_nH,c_float)
               do igroup=1,9
                  trial_energy_shift(:,igroup,i)=ne(:,igroup)/scale_nH- &
                       snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)
               enddo
               absorbed_dust_group(i,:)=real(gn/scale_nH,c_float)
               absorbed_dust_energy_ev(i,:)=ge/scale_nH
               nh_code=h_number_code(i)*scale_nH*atomic_mh/scale_d
               he_code=he_number_code(i)*scale_nH*atomic_mh/scale_d
               trial_hydrogen_ii(i)=chemical_trial(3,i)/nh_code
               trial_neutral_hydrogen(i)=1-trial_hydrogen_ii(i)
               trial_helium_ii(i)=0;trial_helium_iii(i)=0
               if(he_code>0)then
                  trial_helium_ii(i)=chemical_trial(6,i)/he_code
                  trial_helium_iii(i)=chemical_trial(7,i)/he_code
               endif
               chemical_absorbed_ev=chemical_absorbed_ev+ledger(7)*cell_volume_code*scale_l**3
               chemical_events=chemical_events+events(1)
               chemical_dissociation_ev=chemical_dissociation_ev+events(2)*cell_volume_code*scale_l**3
             end block
          enddo
!$omp end parallel do
       endif
#endif
       if (local_transaction_failure == snrt_failure_none) then
          do i = 1, nleaf
             do igroup = 1, snrt_ngroups
                dust_absorbed_photons(igroup,i) = &
                     real(absorbed_dust_group(i,igroup),dp) * scale_nH
             end do
          end do
          if(snrt_d03_band_enabled())then
          call snrt_dust_receiver_stage(dust_absorbed_photons, &
               snrt_dust_contract_absorption_mean_energy_ev(1:snrt_ngroups),dt_s, &
               dust_receiver_abundance,dust_heat_capacity,dust_old_energy,dust_old_temperature, &
               dust_trial_energy,dust_trial_temperature,dust_absorbed_energy,ierr, &
               defer_temperature=.true., &
               deposited_spectrum_erg_cm3=transpose(absorbed_dust_energy_ev)*scale_nH*snrt_ev_to_erg)
          else if(.not.dust_relative_motion)then
          call snrt_dust_receiver_stage(dust_absorbed_photons, &
               snrt_dust_contract_absorption_mean_energy_ev(1:snrt_ngroups), dt_s, &
               dust_receiver_abundance, dust_heat_capacity, dust_old_energy, &
               dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
               dust_absorbed_energy, ierr, defer_temperature=snrt_dust_contract_version==4)
          else
          call snrt_dust_receiver_stage(dust_absorbed_photons, &
               snrt_dust_contract_absorption_mean_energy_ev(1:snrt_ngroups),dt_s, &
               dust_receiver_abundance,dust_heat_capacity,dust_old_energy,dust_old_temperature, &
               dust_trial_energy,dust_trial_temperature,dust_absorbed_energy,ierr, &
               defer_temperature=.true.,deposited_spectrum_erg_cm3=primary_heat)
          endif
          if (ierr /= 0) local_transaction_failure = snrt_failure_receiver
       end if
#endif

       if (local_transaction_failure == snrt_failure_none) then
          unassigned_absorption_total = 0.0d0
          do i = 1, nleaf
             icell = leaf_cell(i)
             rho_code = rho_level(i)
             if (rho_code <= 0.0d0) cycle
             if(chimes_on)cycle
             ! With no absorbed photons this bundle has no local chemistry
             ! source to advance.  Do not reject an otherwise untouched cell
             ! merely because its hydro internal-energy reconstruction cannot
             ! provide a positive temperature.  The receiver failure hook is
             ! checked here as well so the initialized, source-free RAMSES
             ! smoke can exercise the full rollback path without inventing a
             ! physical photon source.
             if (sum(abs(real(absorbed_group(i,:),dp))) <= 0.0d0.and..not.atomic_cooling_on) then
                if (snrt_transaction_failure_requested(transaction_config, &
                     snrt_failure_receiver, i)) then
                   local_transaction_failure = snrt_failure_receiver
                   exit
                end if
                cycle
             end if

             trial_unassigned(i) = 0.0d0
             if (sum(abs(real(absorbed_hhe_group_species(i,:,:),dp))) <= 0.0d0 .and. &
                  sum(abs(real(absorbed_group(i,:),dp))) > 0.0d0) then
#ifndef DUST_LIVE
                ! A future nonzero dust source must supply a dust receiver;
                ! fail closed instead of pretending dust absorption heats gas.
                local_transaction_failure = snrt_failure_receiver
                exit
#else
                ! Dust-only absorption has already been staged by the live
                ! dust receiver.  It must not be reinterpreted as gas heating.
                if (snrt_transaction_failure_requested(transaction_config, &
                     snrt_failure_receiver, i)) then
                   local_transaction_failure = snrt_failure_receiver
                   exit
                end if
                if(.not.atomic_cooling_on)cycle
#endif
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
       n_hydrogen_cm3 = h_number_code(i) * scale_nH
       n_helium_cm3 = n_hydrogen_cm3 * snrt_nhelium_per_hydrogen
       if(atomic_cooling_on)n_helium_cm3=he_number_code(i)*scale_nH
       cell_excess_ev=excess_energy_ev
       if(band_on)then
          block
            integer::s,g
            real(dp)::count,mean
            cell_excess_ev=0
            do g=1,snrt_ngroups
               do s=1,3
                  count=trial_absorbed_species(i,s,g)
                  if(count==0)then
                     if(hhe_energy(i,g,s)/=0)local_transaction_failure=snrt_failure_chemistry
                     cycle
                  endif
                  mean=hhe_energy(i,g,s)/count
                  if(.not.ieee_is_finite(mean).or.mean<band_threshold(s)*(1-8*epsilon(0.0_c_float))) &
                       local_transaction_failure=snrt_failure_chemistry
                  cell_excess_ev(s,g)=max(0.0d0,mean-band_threshold(s))
               enddo
            enddo
          end block
          if(local_transaction_failure/=snrt_failure_none)exit
       endif
       if(allocated(band_deposition))then
          call snrt_thermochemistry_advance_cell(n_hydrogen_cm3,n_helium_cm3, &
               scale_nH,temperature_level(i),dt_s,start_hydrogen_ii(i),start_helium_ii(i),start_helium_iii(i), &
               trial_absorbed_species(i,:,:),cell_excess_ev,chemistry_result, &
               defer_recombination=atomic_cooling_on,band_deposition=band_deposition(i,:))
          ! Independent absorption-energy closure. Excitation represents
          ! escaping line energy here, NOT a transported radiation carrier.
          if(abs(chemistry_result%absorbed_photon_energy_ev_cm3-sum(hhe_energy(i,:,:))*scale_nH)> &
               2d-12*max(sum(hhe_energy(i,:,:))*scale_nH,tiny(1.0_dp))) &
               chemistry_result%ierr=snrt_failure_chemistry
       else
       call snrt_thermochemistry_advance_cell(n_hydrogen_cm3, n_helium_cm3, &
            scale_nH, temperature_level(i), dt_s, start_hydrogen_ii(i), &
            start_helium_ii(i), start_helium_iii(i), &
            trial_absorbed_species(i,:,:), cell_excess_ev, chemistry_result, &
            defer_recombination=atomic_cooling_on)
       endif
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
             if(atomic_cooling_on)then
                block
                  real(dp)::nonthermal,e0,e1,x0(3),x1(3),net_loss
                  nonthermal=.5d0*sum(uold(icell,2:ndim+1)**2)/rho_code+magnetic_energy(uold(icell,:))
#if NENER>0
                  nonthermal=nonthermal+sum(uold(icell,inener:inener+NENER-1))
#endif
                  e0=(trial_thermal(i)-nonthermal)*scale_d*scale_v**2
                  x0=[trial_hydrogen_ii(i),trial_helium_ii(i),trial_helium_iii(i)]
                  call atomic_advance(rho_code*scale_d,atomic_gas_x(:,i),gamma,aexp,dt_s, &
                       e0,x0,e1,x1,net_loss,ierr)
                  if(ierr/=0)then
                     local_transaction_failure=snrt_failure_chemistry
                     chemistry_failures=chemistry_failures+1
                     cycle
                  endif
                  trial_thermal(i)=nonthermal+e1/(scale_d*scale_v**2)
                  trial_hydrogen_ii(i)=x1(1);trial_helium_ii(i)=x1(2);trial_helium_iii(i)=x1(3)
                  trial_neutral_hydrogen(i)=1-x1(1)
                end block
             endif
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
             neutral_hydrogen_code = h_number_code(i) * 0.5d0 * &
                  (start_neutral_hydrogen(i) + &
                  max(0.0d0,1.0d0-relaxed_hydrogen_ii(i)))
             neutral_helium_i_code = he_number_code(i) * 0.5d0 * &
                  (max(0.0d0,1.0d0-start_helium_ii(i)-start_helium_iii(i)) + &
                  max(0.0d0,1.0d0-relaxed_helium_ii(i)-relaxed_helium_iii(i)))
             neutral_helium_ii_code = he_number_code(i) * 0.5d0 * &
                  (start_helium_ii(i) + relaxed_helium_ii(i))
             if(band_on)target_columns(i,:)=[neutral_hydrogen_code,neutral_helium_i_code,neutral_helium_ii_code]* &
                  (scale_nH*snrt_c_cgs*reduced_c*dt_s)
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
             if(chimes_on)then
                target_species_tau(i,:,:)=0
                if(band_on)target_columns(i,:)=0
             endif
             target_tau(i,:) = sum(target_species_tau(i,:,:),dim=2) + &
                  optical_depth_dust(i,:)
          end do
       end if
       if (local_transaction_failure == snrt_failure_none) then
          if(snrt_chimes_cold_enabled())then
             ! This declared Lie split has no opacity feedback into transport:
             ! transport is scattering-only, and CVODE already solves gas and
             ! grain competition simultaneously. Repeating the identical solve
             ! until relaxed diagnostic fractions agree changes no physics.
             transaction_converged=.true.;transaction_residual=0;convergence_status=0
          else
          call snrt_transaction_check_convergence(current_fraction, target_fraction, &
               iteration_tau, target_tau, transaction_config, transaction_residual, &
               transaction_converged, convergence_status)
          endif
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
       if(band_on)species_columns=target_columns
       iteration_tau = sum(iteration_species_tau, dim=3) + optical_depth_dust
    end do

    if (global_transaction_failure /= snrt_failure_none .or. &
         .not. transaction_converged .or. global_transaction_converged == 0) then
       if (transaction_active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, transaction_status,persistent_energy_shift=snrt_energy_shift)
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
          if (transaction_diagnostic_mode .and. &
               global_transaction_failure /= snrt_failure_none) &
               write(*,'(A,A,A,I0)') ' SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=', &
               trim(snrt_transaction_failure_name(global_transaction_failure)), &
               ' level=', ilevel
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
       deallocate(source_transaction_ok)
#ifdef DUST_LIVE
       deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
            dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
            dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
            dust_path_cm, dust_tau_dp)
#endif
       call clean_stop
       return
    end if

#ifdef DUST_LIVE
    ! Validate the dust state and its code-unit conversion before committing
    ! any persistent RT/chemistry state.  The dust fields are intentionally
    ! outside the existing transaction snapshot; they remain untouched until
    ! this collective pre-commit check has passed everywhere.
    local_transaction_failure = snrt_failure_none
#ifdef SNRT_CHIMES
    if(chimes_on.and..not.snrt_chimes_cold_enabled())then
       chemical_absorbed_ev=0
       ! Tables are read-only; cell abundances, rates, CVODE workspaces and
       ! grain temperature are private to each native call.
!$omp parallel do default(shared) private(i,icell,ierr,igroup) reduction(+:chemical_absorbed_ev) &
!$omp reduction(max:local_transaction_failure)
       do i=1,nleaf
          block
            real(dp)::photons(9),next_photons(9),area,nh_code,he_code,fraction,actual(snrt_ndirection),delta_heat
            real(dp)::ray_number(snrt_ndirection,9),ray_energy(snrt_ndirection,9)
            real(dp)::next_number(snrt_ndirection,9),next_energy(snrt_ndirection,9),photo_ledger(9)
            icell=leaf_cell(i)
            do igroup=1,9
               ! Source/transport stores direction-INTEGRATED counts, not I.
               ! Applying quadrature weights again loses photons here.
               photons(igroup)=sum(real(trial_intensity(:,igroup,i),dp))*scale_nH
            enddo
            area=dust_cell_area(i)*dust_relative_abundance(i)
            if(dust_iron_enabled().or.dust_pah_charged())then
               call dust_composition_area(uold(icell,idust_bins:idust_bins+3), &
                    snrt_dust_contract_mass_per_h_g,area,ierr)
               if(ierr/=0)then
                  local_transaction_failure=snrt_failure_receiver
                  cycle
               endif
               area=area*sum(uold(icell,idust_species:idust_species+1))*scale_d/ &
                    (dust_n_hydrogen_cm3(i)*snrt_dust_contract_mass_per_h_g)
            endif
            if(snrt_chimes_band_enabled())then
            ray_number=real(trial_intensity(:,:,i),dp)*scale_nH
            do igroup=1,9
               ray_energy(:,igroup)=(snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)+ &
                    trial_energy_shift(:,igroup,i))*scale_nH
            enddo
            call chimes_live_band_stage(icell,scale_d,scale_v,dt_s,dx_code*scale_l,dust_old_temperature(i), &
                 reduced_c,snrt_ndirection,ray_number,ray_energy,chemical_trial(:,i),trial_thermal(i), &
                 next_number,next_energy,photo_ledger,ierr)
            else if(dust_relative_motion)then
            call chimes_live_stage(icell,scale_d,scale_v,dt_s,dx_code*scale_l,dust_old_temperature(i), &
                 area,reduced_c,photons,chemical_trial(:,i),trial_thermal(i),next_photons,ierr,staged_row=phase_rows(:,i))
            else
            call chimes_live_stage(icell,scale_d,scale_v,dt_s,dx_code*scale_l,dust_old_temperature(i), &
                 area,reduced_c,photons,chemical_trial(:,i),trial_thermal(i),next_photons,ierr)
            endif
            if(ierr/=0)then
               local_transaction_failure=snrt_failure_chemistry
               cycle
            endif
            nh_code=h_number_code(i)*scale_nH*atomic_mh/scale_d
            he_code=he_number_code(i)*scale_nH*atomic_mh/scale_d
            trial_hydrogen_ii(i)=chemical_trial(3,i)/nh_code
            ! Preserve the legacy checkpoint complement field. In this mode
            ! true HI and molecular hydrogen are in the 157 passive species;
            ! this compatibility value must not be used as a HI opacity.
            trial_neutral_hydrogen(i)=1d0-trial_hydrogen_ii(i)
            trial_helium_ii(i)=0;trial_helium_iii(i)=0
            if(he_code>0)then
               trial_helium_ii(i)=chemical_trial(6,i)/he_code
               trial_helium_iii(i)=chemical_trial(7,i)/he_code
            endif
            if(snrt_chimes_band_enabled())then
               ! Transported N is FP32; store the actual outgoing FP64 E,
               ! rebased to the rounded N, not a uniform group attenuation.
               trial_intensity(:,:,i)=real(next_number/scale_nH,c_float)
               do igroup=1,9
                  trial_energy_shift(:,igroup,i)=next_energy(:,igroup)/scale_nH- &
                       snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)
               enddo
               chemical_absorbed_ev=chemical_absorbed_ev+photo_ledger(7)*cell_volume_code*scale_l**3
            else
            do igroup=1,9
               if(dust_relative_motion.and.photons(igroup)>0)then
                  fraction=next_photons(igroup)/photons(igroup)
                  actual=snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)+ &
                       trial_energy_shift(:,igroup,i)
                  ! Frozen group-grey reaction tables retain their reference
                  ! yields. Split off the actual-minus-reference energy of
                  ! accepted captures; this is not a second photon debit.
                  delta_heat=(1-fraction)*sum(trial_energy_shift(:,igroup,i))*scale_nH*snrt_ev_to_erg
                  trial_thermal(i)=trial_thermal(i)+delta_heat/(scale_d*scale_v**2)
                  chemical_absorbed_ev=chemical_absorbed_ev+delta_heat/snrt_ev_to_erg*cell_volume_code*scale_l**3
                  trial_intensity(:,igroup,i)=real(real(trial_intensity(:,igroup,i),dp)*fraction,c_float)
                  trial_energy_shift(:,igroup,i)=actual*fraction- &
                       snrt_group_mean_energy_ev(igroup)*real(trial_intensity(:,igroup,i),dp)
                  cycle
               endif
               if(photons(igroup)>0)trial_intensity(:,igroup,i)=real( &
                    real(trial_intensity(:,igroup,i),dp)*next_photons(igroup)/photons(igroup),c_float)
            enddo
            chemical_absorbed_ev=chemical_absorbed_ev+ &
                 sum((photons-next_photons)*snrt_group_mean_energy_ev)*cell_volume_code*scale_l**3
            endif
            if(dust_relative_motion)then
               delta_heat=trial_thermal(i)-dust_phase_kinetic(phase_rows(:,i))
#if NENER>0
               delta_heat=delta_heat-sum(phase_rows(inener:inener+NENER-1,i))
#endif
               if(.not.ieee_is_finite(delta_heat).or.delta_heat<=0) &
                    local_transaction_failure=snrt_failure_chemistry
            endif
          end block
       enddo
!$omp end parallel do
       call snrt_transaction_reduce_decision(local_transaction_failure,1,0d0,global_transaction_failure, &
            global_transaction_converged,global_transaction_residual,convergence_status)
       if(global_transaction_failure/=snrt_failure_none.or.convergence_status/=0)then
          if(transaction_active)call snrt_transaction_restore(transaction,snrt_intensity,leaf_slot, &
               snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction,level_thermal,transaction_status, &
               persistent_energy_shift=snrt_energy_shift)
          if(myid==1)write(*,*)'ERROR: native CHIMES trial rejected; RT/chemistry not committed'
          call clean_stop
          return
       endif
       call snrt_transaction_reduce_sum(chemical_absorbed_ev,global_unassigned_absorption,convergence_status)
       if(myid==1)write(*,'(A,ES18.10)')' SNRT_CHIMES_PRIMARY_ABSORBED_EV=',global_unassigned_absorption
    endif
#endif
    if(snrt_chimes_cold_enabled())then
       call snrt_transaction_reduce_sum(chemical_absorbed_ev,global_unassigned_absorption,convergence_status)
       if(myid==1)write(*,'(A,ES18.10)')' SNRT_CHIMES_PRIMARY_ABSORBED_EV=',global_unassigned_absorption
       call snrt_transaction_reduce_sum(sum(absorbed_dust_energy_ev)*scale_nH*cell_volume_code*scale_l**3, &
            global_unassigned_absorption,convergence_status)
       if(myid==1)write(*,'(A,ES18.10)')' SNRT_CHIMES_GRAIN_ABSORBED_EV=',global_unassigned_absorption
       if(snrt_chimes_transition_enabled())then
          call snrt_transaction_reduce_sum(chemical_events,global_unassigned_absorption,convergence_status)
          if(myid==1)write(*,'(A,ES18.10)')' SNRT_CHIMES_ATOMIZATION_EVENTS=',global_unassigned_absorption
          call snrt_transaction_reduce_sum(chemical_dissociation_ev,global_unassigned_absorption,convergence_status)
          if(myid==1)write(*,'(A,ES18.10)')' SNRT_CHIMES_DISSOCIATION_COST_EV=',global_unassigned_absorption
       endif
    endif
    if(dust_iron_enabled().and..not.snrt_fe_band_enabled())then
       do i=1,nleaf
          if(sum(uold(leaf_cell(i),idust_iron:idust_iron+1))<=0)cycle
          if(any(dust_absorbed_photons(:,i)>0.and. &
               snrt_dust_contract_absorption_mean_energy_ev(1:snrt_ngroups)>dust_fe_max_primary_ev)) &
               local_transaction_failure=snrt_failure_receiver
       enddo
       call snrt_transaction_reduce_decision(local_transaction_failure,1,0d0,global_transaction_failure, &
            global_transaction_converged,global_transaction_residual,convergence_status)
       if(global_transaction_failure/=snrt_failure_none.or.convergence_status/=0)then
          if(transaction_active)call snrt_transaction_restore(transaction,snrt_intensity,leaf_slot, &
               snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction,level_thermal,transaction_status, &
               persistent_energy_shift=snrt_energy_shift)
          if(myid==1)write(*,*)'ERROR: Fe comparison primary absorption exceeds 4 eV; no trial committed'
          call clean_stop;return
       endif
    endif
    if(snrt_dust_contract_version>=3)then
             ! Start from the pre-primary material energy and inject exactly
             ! the accepted primary absorption. Receiver-stage energy already
             ! includes that absorption and must NOT be fed as old energy.
             ! Primary quadrature integrates over 4*pi. IR stores energy per
             ! normalized direction and therefore requires weights summing 1.
             ! Stage validates faces and reduces errors collectively before
             ! halo exchange. Never skip this call on a rank-local condition.
       if(dust_pah_enabled())then
          do i=1,snrt_ngroups
             pah_primary_energy(i,:)=dust_absorbed_photons(i,:)* &
                  snrt_dust_contract_absorption_mean_energy_ev(i)*snrt_ev_to_erg
          enddo
          if(dust_relative_motion)pah_primary_energy=primary_heat
       endif
       if(snrt_dust_contract_exchange_enabled)then
          block
            real(dp) :: gas_energy(nleaf),gas_capacity(nleaf),transfer(nleaf)
            real(dp),allocatable :: pah_electrons(:),pah_gas_h(:),pah_gas_h2(:)
            real(dp)::electron_cv
            electron_cv=0
#ifdef SNRT_CHIMES
            if(dust_pah_charged())then
               pah_electrons=chemical_trial(1,:)*scale_d/atomic_mh
               electron_cv=1.5d0*chimes_boltzmann()
               if(dust_pah_hydrogenated())pah_gas_h=chemical_trial(2,:)*scale_d/atomic_mh
               ! Pinned CHIMES157: H2 is index138, one native weight per
               ! molecule, not twice that weight as an H-mass fraction.
               if(dust_pah_h2_enabled())pah_gas_h2=chemical_trial(138,:)*scale_d/atomic_mh
            endif
#endif
            do i=1,nleaf
               icell=leaf_cell(i)
               kinetic_energy=.5d0*sum(uold(icell,2:ndim+1)**2)/rho_level(i)
               if(dust_relative_motion)kinetic_energy=dust_phase_kinetic(phase_rows(:,i))
               gas_energy(i)=(trial_thermal(i)-kinetic_energy-magnetic_energy(uold(icell,:)))*dust_energy_scale
#if NENER>0
               gas_energy(i)=gas_energy(i)-sum(uold(icell,inener:inener+NENER-1))*dust_energy_scale
#endif
               molecular_weight=snrt_mean_molecular_weight(trial_hydrogen_ii(i),trial_helium_ii(i),trial_helium_iii(i))
               gas_capacity(i)=rho_level(i)*dust_energy_scale/((gamma-1)*scale_T2*molecular_weight)
               if(atomic_cooling_on)gas_capacity(i)=atomic_heat_capacity(rho_level(i)*scale_d, &
                    atomic_gas_x(:,i),[trial_hydrogen_ii(i),trial_helium_ii(i),trial_helium_iii(i)],gamma)
#ifdef SNRT_CHIMES
               if(chimes_on)gas_capacity(i)=chimes_live_capacity(chemical_trial(:,i),scale_d)
#endif
            enddo
            ! Solve exchange and IR emission together, inside the IR implicit
            ! material solve, not as a post-radiation dust temperature kick.
            call snrt_dust_live_stage(ilevel,leaf_cell,leaf_slot,neighbor, &
                 transpose(direction_dp),angular_weight/sum(angular_weight),dx_code*scale_l,dt_s,snrt_c_cgs*reduced_c, &
                 dust_n_hydrogen_cm3*dust_relative_abundance,dust_absorbed_energy,dust_old_energy, &
                 dust_heat_capacity,dust_ir_trial,dust_trial_energy,dust_trial_temperature,dust_ir_result,ierr, &
                 dust_ir_coarse,gas_energy,gas_capacity,dust_n_hydrogen_cm3,transfer,dust_cell_u,dust_cell_area,dust_weights, &
                 sublimation_bins,sublimation_next,pah_number,pah_primary_energy, &
                 primary_pah_heat=primary_pah_heat,primary_pah_captures=primary_pah_captures, &
                 phase_density=phase_mass,phase_momentum=phase_pnext,phase_work=phase_ir_work, &
                 gas_electrons=pah_electrons,electron_capacity=electron_cv,gas_atomic_h=pah_gas_h, &
                 gas_molecular_h2=pah_gas_h2)
            if(ierr==0)then
#ifdef SNRT_CHIMES
               if(dust_pah_charged())chemical_trial(1,:)=pah_electrons*atomic_mh/scale_d
               if(dust_pah_hydrogenated())chemical_trial(2,:)=pah_gas_h*atomic_mh/scale_d
               if(dust_pah_h2_enabled())chemical_trial(138,:)=pah_gas_h2*atomic_mh/scale_d
#endif
               trial_thermal=trial_thermal-transfer/dust_energy_scale
               if(dust_relative_motion)then
                  trial_thermal=trial_thermal+sum(phase_ir_work,dim=1)/dust_energy_scale
                  do i=1,nleaf
                     block
                       integer::b,k
                       phase_rows(2:4,i)=uold(leaf_cell(i),2:4)+ &
                            sum(phase_pnext(:,:,i)-phase_pold(:,:,i),dim=2)/(scale_d*scale_v)
                       do b=1,ndust_phase
                          k=idust_momentum+3*(b-1)
                          phase_rows(k:k+2,i)=phase_pnext(:,b,i)/(scale_d*scale_v)
                       enddo
                     end block
                  enddo
               endif
               if(any(.not.ieee_is_finite(trial_thermal)).or.any(trial_thermal<=0))ierr=1
            endif
          end block
       else
             call snrt_dust_live_stage(ilevel,leaf_cell,leaf_slot,neighbor, &
                  transpose(direction_dp),angular_weight/sum(angular_weight), &
                  dx_code*scale_l,dt_s,snrt_c_cgs*reduced_c, &
                  dust_n_hydrogen_cm3*dust_relative_abundance,dust_absorbed_energy,dust_old_energy, &
                  dust_heat_capacity,dust_ir_trial,dust_trial_energy,dust_trial_temperature,dust_ir_result,ierr, &
                  dust_ir_coarse,cell_material_u=dust_cell_u,cell_collision_area=dust_cell_area,cell_weights=dust_weights)
       endif
       if(ierr/=0)then
          local_transaction_failure=snrt_failure_receiver
          if(myid==1)write(*,'(A,I0)')' SNRT live IR staging failed: error=',ierr
       end if
    end if
    if (any(.not. ieee_is_finite(dust_trial_energy)) .or. &
         any(.not. ieee_is_finite(dust_trial_temperature)) .or. &
         any(dust_trial_energy < 0.0d0) .or. &
         any(dust_trial_temperature <= 0.0d0)) then
       local_transaction_failure = snrt_failure_receiver
    else
       do i = 1, nleaf
          dust_energy_code = dust_trial_energy(i) / dust_energy_scale
          if (.not. ieee_is_finite(dust_energy_code) .or. &
               dust_energy_code < 0.0d0) then
             local_transaction_failure = snrt_failure_receiver
             exit
          end if
       end do
    end if
#ifdef SNRT_CHIMES
    if(allocated(sublimation_next).and.local_transaction_failure==snrt_failure_none)then
       block
         real(dp)::next_chemical(chimes_ns),elements(11),grains(2)
         do i=1,nleaf
            grains=[sum(sublimation_next(1:2,i)),sum(sublimation_next(3:4,i))]/scale_d
            call chimes_cell_state(leaf_cell(i),grains,next_chemical,elements,ierr,chemical_trial(:,i))
            if(ierr/=0)then
               local_transaction_failure=snrt_failure_receiver
               exit
            endif
            chemical_trial(:,i)=next_chemical
         enddo
       end block
    endif
#endif
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0) global_transaction_failure = snrt_failure_receiver
    if (global_transaction_failure /= snrt_failure_none) then
       if (transaction_active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, transaction_status,persistent_energy_shift=snrt_energy_shift)
       if (myid == 1) write(*,'(A,I0)') &
            ' SNRT DUST_LIVE pre-commit validation failed at level=', ilevel
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
       deallocate(source_transaction_ok)
       deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
            dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
            dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
            dust_path_cm, dust_tau_dp)
       call clean_stop
       return
    end if
#endif

    if(paired_transport)then
       block
         integer::slot,g,leaf_of_slot(size(snrt_intensity,3))
         real(c_float)::base(snrt_ndirection),combined(snrt_ndirection)
         leaf_of_slot=0
         do i=1,nleaf
            leaf_of_slot(leaf_slot(i))=i
         enddo
         do slot=1,size(snrt_intensity,3)
            if(all(coarse_flux_trial(:,:,slot)==0.0_c_float))cycle
            do g=1,snrt_ngroups
               base=snrt_intensity(:,g,slot)
               if(leaf_of_slot(slot)>0)base=trial_intensity(:,g,leaf_of_slot(slot))
               combined=base+coarse_flux_trial(:,g,slot)
               ! The topology already rebased accumulation into the flux
               ! register; this is the separate final addition to cell N.
               coarse_energy_shift(:,g,slot)=coarse_energy_shift(:,g,slot)+snrt_group_mean_energy_ev(g)* &
                    (real(base,dp)+real(coarse_flux_trial(:,g,slot),dp)-real(combined,dp))
            enddo
         enddo
       end block
    endif
    ! Validate every proposed leaf AND coarse-slot energy before any rank
    ! publishes and discards its snapshot. A local successful commit cannot
    ! be rolled back after another rank rejects; this collective comes first.
    call snrt_transaction_commit_level(transaction,snrt_intensity,leaf_slot, &
         snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction, &
         trial_intensity,coarse_flux_trial,trial_hydrogen_ii,trial_helium_ii, &
         trial_helium_iii,trial_neutral_hydrogen,level_thermal,trial_thermal,transaction_status, &
         persistent_energy_shift=snrt_energy_shift,trial_energy_shift=trial_energy_shift, &
         coarse_energy_shift_trial=coarse_energy_shift,validate_only=.true.)
    local_transaction_failure=snrt_failure_none
    if(transaction_status/=0)local_transaction_failure=snrt_failure_receiver
    call snrt_transaction_reduce_decision(local_transaction_failure,1,0d0,global_transaction_failure, &
         global_transaction_converged,global_transaction_residual,convergence_status)
    if(global_transaction_failure/=snrt_failure_none.or.convergence_status/=0)then
       if(transaction%active)call snrt_transaction_restore(transaction,snrt_intensity,leaf_slot, &
            snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction,level_thermal,transaction_status, &
            persistent_energy_shift=snrt_energy_shift)
       if(myid==1)write(*,*)'ERROR: collective SNRT commit validation rejected; no primary/material state published'
       call clean_stop;return
    endif
    call snrt_transaction_commit_level(transaction, snrt_intensity, leaf_slot, &
         snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, snrt_neutral_fraction, &
         trial_intensity, coarse_flux_trial, trial_hydrogen_ii, trial_helium_ii, &
         trial_helium_iii, trial_neutral_hydrogen, level_thermal, trial_thermal, &
         transaction_status,persistent_energy_shift=snrt_energy_shift,trial_energy_shift=trial_energy_shift, &
         coarse_energy_shift_trial=coarse_energy_shift)
    local_transaction_failure = snrt_failure_none
    if (transaction_status /= 0) local_transaction_failure = snrt_failure_receiver
    call snrt_transaction_reduce_decision(local_transaction_failure, 1, 0.0d0, &
         global_transaction_failure, global_transaction_converged, &
         global_transaction_residual, convergence_status)
    if (convergence_status /= 0) global_transaction_failure = snrt_failure_receiver
    if (global_transaction_failure /= snrt_failure_none) then
       if (transaction%active) call snrt_transaction_restore(transaction, snrt_intensity, &
            leaf_slot, snrt_hydrogen_ii, snrt_helium_ii, snrt_helium_iii, &
            snrt_neutral_fraction, level_thermal, ierr,persistent_energy_shift=snrt_energy_shift)
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
       deallocate(source_transaction_ok)
#ifdef DUST_LIVE
       deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
            dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
            dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
            dust_path_cm, dust_tau_dp)
#endif
       call clean_stop
       return
    end if
    call snrt_transaction_reduce_sum(unassigned_absorption_total, &
         global_unassigned_absorption, convergence_status)
    if (convergence_status /= 0) global_unassigned_absorption = &
         unassigned_absorption_total
    if (myid == 1) then
       write(*,'(A,I0,A,I0,A,ES12.4)') &
            ' SNRT_RT_TRANSACTION_COMMIT_PASS level=', ilevel, &
            ' iteration=', transaction_iteration, &
            ' residual=', global_transaction_residual
       write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.4)') &
            ' SNRT_RT_CLOSURE_PASS level=', ilevel, ' leaves=', nleaf, &
            ' photon_nonnegative=', 1, ' species_simplex=', 1, &
            ' thermal_finite=', 1, ' unassigned_code=', global_unassigned_absorption
    end if
    do i = 1, nleaf
       icell = leaf_cell(i)
       if (icell >= 1 .and. icell <= size(uold,1) .and. &
            size(uold,2) >= energy_index) uold(icell,energy_index) = level_thermal(i)
#ifdef SNRT_CHIMES
       if(chimes_on)uold(icell,ichimes:ichimes+chimes_ns-1)=chemical_trial(:,i)
#endif
#ifdef DUST_LIVE
       if(dust_relative_motion)then
          uold(icell,2:4)=phase_rows(2:4,i)
          uold(icell,idust_momentum:idust_momentum+3*ndust_phase-1)= &
               phase_rows(idust_momentum:idust_momentum+3*ndust_phase-1,i)
       endif
       if (icell >= 1 .and. icell <= size(uold,1) .and. &
            size(uold,2) >= idust_energy) then
          dust_energy_code = dust_trial_energy(i) / dust_energy_scale
          uold(icell,idust_energy) = dust_energy_code
          if(dust_pah_enabled())then
             block
               integer::k
               do k=1,dust_pah_nstate()
                  uold(icell,idust_pah+k-1)=pah_number(k,i)*dust_pah_state_mass(k)/scale_d
               enddo
             end block
          endif
          if(allocated(sublimation_next))then
             uold(icell,idust_bins:idust_bins+3)=sublimation_next(:,i)/scale_d
             uold(icell,idust_species:idust_species+1)= &
                  [sum(sublimation_next(1:2,i)),sum(sublimation_next(3:4,i))]/scale_d
             uold(icell,idust)=sum(sublimation_next(:,i))/scale_d
          endif
       end if
#endif
    end do
#ifdef DUST_LIVE
    if(dust_relative_motion)then
       do i=2,4
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
       do i=idust_momentum,idust_momentum+3*ndust_phase-1
          call make_virtual_fine_dp(uold(1,i),ilevel)
       enddo
       call make_virtual_fine_dp(uold(1,energy_index),ilevel)
    endif
    if(snrt_dust_contract_version>=3)then
       call snrt_dust_live_commit(leaf_slot,dust_ir_trial,dust_ir_coarse)
       if(dust_pah_enabled())then
          do i=idust_pah,idust_pah+dust_pah_nstate()-1
             call make_virtual_fine_dp(uold(1,i),ilevel)
          enddo
       endif
       if(allocated(sublimation_next))then
          do i=idust,idust_bins+3
             call make_virtual_fine_dp(uold(1,i),ilevel)
          enddo
          call make_virtual_fine_dp(uold(1,energy_index),ilevel)
#ifdef SNRT_CHIMES
          do i=ichimes,ichimes+chimes_ns-1
             call make_virtual_fine_dp(uold(1,i),ilevel)
          enddo
#endif
       endif
       if(myid==1)write(*,'(A,ES12.4,A,ES12.4)') &
            ' SNRT_DUST_IR_COMMIT_PASS balance=',dust_ir_result%balance_relative, &
            ' escaped_erg=',dust_ir_result%escaped_erg
    end if
#endif
    ! Fuel consumption is the final accounting action of the coupled source
    ! -> RT/chemistry -> dust commit.  No rollback path can observe a cleared
    ! marker because every failure branch returns before this point.
    if (allocated(agn_pending_erg) .and. allocated(source_transaction_ok)) then
       if (size(agn_pending_erg) >= nsink .and. size(source_transaction_ok) >= nsink) then
          do isink = 1, nsink
             if (source_transaction_ok(isink)) &
                  call snrt_agn_source_commit(agn_pending_erg(isink), .true.)
          end do
       end if
    end if
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
    deallocate(source_transaction_ok)
#ifdef DUST_LIVE
    deallocate(dust_relative_abundance, dust_heat_capacity, dust_old_energy, &
         dust_old_temperature, dust_trial_energy, dust_trial_temperature, &
         dust_absorbed_photons, dust_absorbed_energy, dust_n_hydrogen_cm3, &
         dust_path_cm, dust_tau_dp)
#endif
  end subroutine snrt_ramses_advance_level

end module snrt_ramses_driver

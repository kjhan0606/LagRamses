! Patch change:
! - added parameters for Kimm feedback and star formation
subroutine read_hydro_params(nml_ok)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_parameters, only: grafic_nreaders
  use amr_commons
  use hydro_commons
  use cosmic_ray_physics
  use dust_mass_physics
#ifdef DUST_DYNAMICS
  use snrt_dust_contract, only: snrt_dust_contract_loaded,snrt_dust_contract_load_from_environment, &
       snrt_dust_contract_version,snrt_dust_contract_exchange_enabled
#endif
#ifdef SNRT
  use snrt_agn_efficiency, only: snrt_agn_rt_requested
#endif
  use eunha_cooling_mod, only: eunha_load_multi_z
#ifdef PHASE0_STELLAR_ENRICHMENT
  use stellar_enrichment_config, only: read_enrichment_namelist, &
       stellar_feedback_mode, use_channel_resolved_feedback, default_imf_id, active_element, &
       population_model_id, yield_source_basis_name, configured_imf_mass_min, &
       configured_imf_mass_max, configured_binary_fraction, stellar_fate_policy, &
       stellar_fate_map_sha256, stellar_fate_approval_id, production_fate_policy_supported, &
       high_mass_model, high_mass_max_remnant_adjust_fraction, user_source_model_requested, high_mass_history_file
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  logical::cr_ok,dust_ok
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  integer::i,idim,ivar,nboundary_true=0,dum
#ifdef PHASE0_STELLAR_ENRICHMENT
  integer::stellar_nml_iostat
#endif
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::scale,ek_bound
  character(LEN=1)::a1

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,multiple,grafic_nreaders,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
#if NENER>0
       & ,prad_region &
#endif
#if NVAR>NDIM+2+NENER
       & ,var_region &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
       & ,mhd_enabled,mhd_omp,mhd_gpu_faces,mhd_seed,mhd_initial_condition,riemann2d,slope_mag_type,interpol_mag_type &
#if NENER>0
       & ,gamma_rad &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,void_refine,void_refine_min_level &
       & ,void_web_refine,void_web_env_level &
       & ,void_web_base_level,void_web_wall_level &
       & ,void_web_hydro_max_level &
       & ,void_web_scope_ivar,void_web_scope_cut &
       & ,void_web_lambda_on,void_web_lambda_off &
       & ,void_web_update_interval &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,err_jump_u,ekin_flux_refine,d_keflux_max,floor_keflux &
       & ,void_web_jump_compression_gate,void_web_jump_pressure_min &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type,sink_refine,d_jeans_thre &
       & ,q_refine_holdback,m_refine_effective,ref_fall_rate & !(ONS)
       & ,dr_proper
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
#if NENER>0
       & ,prad_bound &
#endif
#if NVAR>NDIM+2+NENER
       & ,var_bound &
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow
  namelist/physics_params/cooling,haardt_madau,metal,isothermal &
       & ,dust_relative_motion,dust_drag_collision_cross_section_cm2 &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star,jeans_ncells &
       & ,rbubble,f_ek,ndebris,f_w,mass_gmc,kappa_IR &
       & ,J21,a_spec,z_ave,z_reion,ind_rsink,delayed_cooling,T2max &
       & ,self_shielding,smbh,agn,yieldtablefilename &
       & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff,t_diss,t_sne &
       & ,T2thres_SF, sn2_real_delay, fstar_min, M_SNII       &
       & ,star_imf, n_dc, n_gmc, nsn2mass, star_maker, SN_dT2_min,sf_lamjt&
       & ,A_SN, expN_SN, E_SNII, SF_kick_kms, write_stellar_densities              &
       & ,n_sink,sink_AGN,rAGN,eAGN_K,eAGN_T,TAGN,X_floor,r_gal,boost_acc,Mseed    &
       & ,sigmav_max,star_ratio_floor,mloadAGN,T2maxAGN,f_bondi,random_jet         &
       & ,drag,boost_drag,selfgrav,spin_bh,bhspinmerge,vrel_merge,rmerge,omega_b   &
       & ,smbh_capture_ledger,smbh_capture_ledger_file                            &
       & ,cdm_zoom_plan_manifest_sha256,cdm_zoom_capture_event_sha256             &
       & ,cdm_zoom_host_orbit_initial_conditions_sha256                           &
       & ,cdm_zoom_initial_conditions_sha256,cdm_zoom_sink_initial_conditions_sha256 &
       & ,model_zoom_manifest_sha256,model_zoom_case_id                           &
       & ,model_zoom_capture_event_sha256,model_zoom_initial_conditions_sha256    &
       & ,model_zoom_baryon_configuration_sha256                                  &
       & ,model_zoom_sink_initial_conditions_sha256                               &
       & ,agn_coarse_dump,agn_coarse_dump_file                                    &
       & ,bondi,mad_jet,eps_sn1,eps_sn2,tol                             &
       & ,sf_virial,sf_trelax,sf_model,sf_birth_properties &
       & ,cr_enabled,cr_transport,cr_sn_fraction,cr_snia_fraction,cr_sf_support &
       & ,dust_mass_enabled,dust_mass_model,dust_cooling,dust_growth,dust_sputtering,dust_condensation &
       & ,dust_grain_radius_cm,dust_grain_density,dust_sticking,dust_growth_max_temperature &
       & ,dust_metal_atom_mass,dust_injection_temperature &
       & ,dust_size_radius_cm,dust_size_density,dust_small_injection_fraction,dust_coagulation,dust_shattering &
       & ,dust_sn_shocks &
       & ,dust_material_model,dust_optics_model,dust_sublimation &
       & ,dust_iron_model,dust_fe_condensation,dust_fe_kinetics,dust_fe_sticking &
       & ,dust_pah_model,dust_pah_condensation &
       & ,cooling_method,grackle_table
#ifdef grackle
  namelist/grackle_params/grackle_comoving_coordinates,grackle_with_radiative_cooling,grackle_primordial_chemistry &
       & ,grackle_metal_cooling,grackle_UVbackground,grackle_h2_on_dust,grackle_cmb_temperature_floor &
       & ,grackle_data_file
#endif

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
#ifndef SOLVERmhd
  if(mhd_omp.or.mhd_gpu_faces)then
     write(*,*)'ERROR: MHD execution flags require SOLVER=mhd'
     call clean_stop
  endif
  if(mhd_enabled.or.any(mhd_seed/=0d0).or.riemann=='hlld'.or.mhd_initial_condition/='uniform')then
     write(*,*)'ERROR: magnetic fields require a SOLVER=mhd executable'
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  if(hydro)then
     if(mhd_gpu_faces)then
#ifndef HYDRO_CUDA
        write(*,*)'ERROR: mhd_gpu_faces requires USE_CUDA=1'
        call clean_stop
#endif
        if(nener/=0)then
           write(*,*)'ERROR: CUDA MHD face solver currently requires NENER=0'
           call clean_stop
        endif
     endif
     select case(trim(mhd_initial_condition))
     case('uniform')
     case('alfven_x','brio_wu_x')
        if(cosmo.or.any(initfile/=' ').or.nener/=0)then
           write(*,*)'ERROR: analytic MHD ICs require noncosmo, empty initfile and NENER=0'
           call clean_stop
        endif
        if(mhd_initial_condition=='brio_wu_x'.and.abs(gamma-2d0)>1d-12)then
           write(*,*)'ERROR: Brio-Wu reference requires gamma=2'
           call clean_stop
        endif
     case default
        write(*,*)'ERROR: unknown mhd_initial_condition'
        call clean_stop
     end select
     if(any(.not.ieee_is_finite(mhd_seed)))then
        write(*,*)'ERROR: mhd_seed must contain finite code-unit components'
        call clean_stop
     endif
     if(.not.mhd_enabled)then
        write(*,*)'ERROR: SOLVER=mhd requires explicit mhd_enabled=.true.'
        call clean_stop
     endif
     if(scheme/='muscl'.or.riemann/='hlld'.or.riemann2d/='hlld')then
        write(*,*)'ERROR: gas MHD requires muscl / hlld / hlld'
        call clean_stop
     endif
     if(gpu_hydro.or.use_sgs.or.dust_relative_motion)then
        write(*,*)'ERROR: gas MHD currently requires CPU hydro, no SGS or separate dust dynamics'
        call clean_stop
     endif
#ifndef HDF5
     write(*,*)'ERROR: gas MHD requires HDF5 checkpoint support'
     call clean_stop
#endif
     if(outformat/='hdf5'.or.(nrestart>0.and.informat/='hdf5'))then
        write(*,*)'ERROR: gas MHD requires HDF5 output/restart (face fields included)'
        call clean_stop
     endif
     ischeme=0;iriemann=3;iriemann2d=5
  endif
#endif
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue
#ifdef SOLVERmhd
  if(hydro.and.(use_sgs.or.dust_relative_motion))then
     write(*,*)'ERROR: gas MHD does not yet admit SGS or separate dust dynamics'
     call clean_stop
  endif
  if(hydro.and.nener>0.and.(sink.or.sink_AGN.or.agn))then
     write(*,*)'ERROR: MHD sink/AGN profile requires NENER=0 (no accreted CR closure)'
     call clean_stop
  endif
#endif
#ifdef PHASE0_STELLAR_ENRICHMENT
  rewind(1)
  call read_enrichment_namelist(1,stellar_nml_iostat)
  if(stellar_nml_iostat==0.and.stellar_fate_policy=='user_selected_model_v1')then
     if(.not.user_source_model_requested().or..not.pic)stellar_nml_iostat=1013
  endif
  if(stellar_nml_iostat/=0)then
     if(myid==1)then
        write(*,*) 'ERROR: invalid &STELLAR_ENRICHMENT_PARAMS namelist'
        if(stellar_nml_iostat==1001) &
             write(*,*) '  feedback_mode must be channel_resolved or legacy'
        if(stellar_nml_iostat==1002) &
             write(*,*) '  imf_id must be in the supported range 0:4'
        if(stellar_nml_iostat==1012) &
             write(*,*) '  invalid high-mass preset/adjustment limit, or override requested in legacy mode'
        if(stellar_nml_iostat==1013) &
             write(*,*) '  user source requires PIC, wind+SNII; SNIa requires matching binary SSP and AGB WD source; PISN off'
        if(stellar_nml_iostat==1003) &
             write(*,*) '  population_model must be single_star_ssp or binary_ssp'
        if(stellar_nml_iostat==1004) &
             write(*,*) '  every channel mass window must be finite and increasing'
        if(stellar_nml_iostat==1005) &
             write(*,*) '  &STELLAR_ENRICHMENT_PARAMS is required in this build'
        if(stellar_nml_iostat==1006) &
             write(*,*) '  yield_source_basis is required and must be recognized'
        if(stellar_nml_iostat==1007) &
             write(*,*) '  IMF mass bounds must be finite, positive, and increasing'
        if(stellar_nml_iostat==1008) &
             write(*,*) '  binary_fraction is invalid for the selected population model'
        if(stellar_nml_iostat==1009) &
             write(*,*) '  an enabled channel mass window lies outside the IMF support'
     end if
     nml_ok=.false.
  end if
  if(myid==1 .and. stellar_nml_iostat==0)then
     write(*,'(A,A)') ' Stellar feedback mode: ',trim(stellar_feedback_mode)
     write(*,'(A,I0)') ' Stellar IMF id: ',default_imf_id
     write(*,'(A,I0)') ' Stellar population model id: ',population_model_id
     write(*,'(A,A)') ' Stellar yield source basis: ',trim(yield_source_basis_name())
     write(*,'(A,2(1X,ES12.4))') ' Stellar IMF mass support:', &
          configured_imf_mass_min,configured_imf_mass_max
     write(*,'(A,1X,ES12.4)') ' Stellar binary fraction:', &
          configured_binary_fraction
     write(*,'(A,A)') ' Stellar terminal-fate policy: ',trim(stellar_fate_policy)
     write(*,'(A,A)') ' Stellar fate-map SHA256: ',trim(stellar_fate_map_sha256)
     write(*,'(A,A)') ' Stellar fate approval id: ',trim(stellar_fate_approval_id)
     write(*,'(A,A)') ' High-mass endpoint preset: ',trim(high_mass_model)
     write(*,'(A,ES16.7)') ' Maximum remnant correction / initial mass: ',high_mass_max_remnant_adjust_fraction
     if(user_source_model_requested())then
        write(*,*) 'Stellar source: USER_SELECTED_MODEL, not an automatic physical approval'
        write(*,*) 'High-mass history: ',trim(high_mass_history_file)
     else if(.not.production_fate_policy_supported())then
        write(*,'(A)') ' Stellar terminal-fate policy is review-only; production admission is blocked'
     endif
  end if
#endif
#ifdef grackle
  rewind(1)
  read(1,NML=grackle_params)
#endif
  ! Validate and initialize Eunha cooling
  if(cooling_method=='exact' .or. cooling_method=='compare') then
     if(len_trim(grackle_table)==0) then
        if(myid==1) write(*,*) 'ERROR: grackle_table required for cooling_method=',trim(cooling_method)
        call clean_stop
     endif
     if(myid==1) write(*,'(A,A,A,A)') ' Cooling method: ', trim(cooling_method), &
          '  table: ', trim(grackle_table)
     call eunha_load_multi_z(trim(grackle_table))
  else
     if(myid==1 .and. cooling) write(*,'(A)') ' Cooling method: original (RAMSES solve_cooling)'
  endif
#ifdef ATON
  if(aton)call read_radiation_params(1)
#endif

  select case(scheme)
  case('muscl','plmde')
  case('weno3')
#if NDIM != 3
     if(myid==1)write(*,*)'scheme=weno3 is currently available only with NDIM=3'
     nml_ok=.false.
#endif
     if(slope_type.ne.2)then
        if(myid==1)write(*,*)'scheme=weno3 requires slope_type=2 for the time predictor'
        nml_ok=.false.
     end if
  case('weno5','weno5ppm','ppm')
#if NDIM != 3
     if(myid==1)write(*,*)'scheme=',trim(scheme),' is currently available only with NDIM=3'
     nml_ok=.false.
#endif
     if(slope_type.ne.2)then
        if(myid==1)write(*,*)'scheme=',trim(scheme),' requires slope_type=2 for the predictor'
        nml_ok=.false.
     end if
     if(levelmin.ne.nlevelmax)then
        if(myid==1)write(*,*)'scheme=',trim(scheme),' prototype currently requires a uniform grid'
        nml_ok=.false.
     end if
     if(poisson)then
        if(myid==1)write(*,*)'scheme=',trim(scheme),' prototype currently requires poisson=.false.'
        nml_ok=.false.
     end if
     if(nboundary.ne.0)then
        if(myid==1)write(*,*)'scheme=',trim(scheme),' prototype currently requires periodic boundaries'
        nml_ok=.false.
     end if
     ! A five-point reconstruction needs two fine-cell layers across MPI
     ! ownership boundaries.  Expand the virtual mesh accordingly.
     nexpand_bound=max(nexpand_bound,2)
  case default
     if(myid==1)write(*,*)'unknown hydro scheme: ',trim(scheme)
     nml_ok=.false.
  end select

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(t_star>0)then
     star=.true.
     pic=.true.
  else if(eps_star>0)then
     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
  if(metal.and.nvar<(ndim+3))then
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+3'
     if(myid==1)write(*,*)'Modify hydro_parameters.f90 and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
  if(nvar<(ndim+2+nener))then
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  dust_ok=dust_mass_parameters_ok()
#ifdef SNRT
  if(dust_sublimation_rt_enabled().and..not.snrt_agn_rt_requested())dust_ok=.false.
#endif
  if(dust_chimes_enabled())then
#ifndef SNRT_CHIMES
     dust_ok=.false.
#endif
     ! CHIMES thermal evolution uses translational 3/2 n_tot k_B T.
     if(abs(gamma-5d0/3)>1d-12)dust_ok=.false.
  endif
  if(dust_mass_enabled)then
#if !defined(SNRT) || !defined(DUST_LIVE) || !defined(HDF5) || !defined(PHASE0_STELLAR_ENRICHMENT)
     dust_ok=.false.
#else
     if(.not.use_channel_resolved_feedback())dust_ok=.false.
     if(dust_composition_enabled().and..not.all(active_element))dust_ok=.false.
#endif
     if(.not.hydro.or..not.metal.or.cosmo.or.nboundary>0) dust_ok=.false.
     ! Total-metal reservoir closure; no element-resolved depleted cooling or sink removal yet.
     if(sink.or.sink_AGN.or.agn.or.neq_chem.or.delayed_cooling) dust_ok=.false.
     if(cooling.neqv.(trim(dust_cooling)/='none'))dust_ok=.false.
     if(trim(dust_cooling)/='none')then
        ! Explicit collisional scalar-Z comparison; no duplicate UV heating.
        if(trim(cooling_method)/='original'.or.haardt_madau.or.J21/=0.or.self_shielding) dust_ok=.false.
#ifdef grackle
        dust_ok=.false.
#endif
     endif
     if(dust_composition_enabled().and.gpu_hydro)dust_ok=.false.
     if(trim(outformat)/='hdf5'.or.(nrestart>0.and.trim(informat)/='hdf5'))dust_ok=.false.
  endif
  if(.not.dust_ok)then
     if(myid==1)write(*,*)'ERROR: dust mass requires valid bulk parameters, SNRT/DUST_LIVE/HDF5 channel feedback,'
     if(myid==1)write(*,*)'noncosmo periodic metal hydro; no sinks/neq/delayed cooling; cooling needs an explicit dust closure'
     nml_ok=.false.
  else if(dust_mass_enabled.and.myid==1)then
     write(*,*)'DUST_MASS model=',trim(dust_mass_model),'; condensation/growth/sputtering, total-metal budget'
     if(dust_composition_enabled())write(*,*)'DUST_COMPOSITION C/MgFeSiO4: source-segment C/O'
     if(dust_two_size_enabled())write(*,*)'DUST_SIZE: four masses; resolved-density coagulation/shattering'
     if(dust_material_composition_enabled())write(*,*)'DUST_MATERIAL: local DL01 composition, common T, geometric size area'
     if(dust_sublimation_enabled())write(*,*)'DUST_SUBLIMATION: ',trim(dust_sublimation),'; vacuum, fixed-radius BE'
     if(dust_sublimation_rt_enabled())write(*,*) &
          'DUST_SUBLIMATION_RT: adaptive IR/material/phase; lagged opacities/Cv; native CPU/OpenMP'
     if(dust_silicate_sublimation_enabled())write(*,*) &
          'DUST_OLIVINE: Xu crystalline surface rates; RH95 ideal-mixture atomic phase reference; congruent vapor'
     if(dust_optics_enabled())then
        write(*,*)'DUST_OPTICS: D03 four populations; primary/IR delta-isotropic Qsca*(1-g); common T'
     else
        write(*,*)'DUST_OPTICS: fixed reference mixture'
     endif
     if(dust_sn_shocks)write(*,*)'DUST_SHOCK: energy-equivalent unresolved SN comparison; fresh ejecta protected; not resolution-calibrated'
     if(trim(dust_cooling)=='depleted_scalar') &
          write(*,*)'DUST_COOLING depleted_scalar: solar-mixture curve at (total metal - dust)/rho, NOT element cooling'
     if(trim(dust_cooling)=='wss09_cie') &
          write(*,*)'DUST_COOLING WSS09_CIE: gas-phase H/He + nine metals; CIE comparison, NOT radiation-dependent NEQ'
     if(trim(dust_cooling)=='snrt_hhe_cie_metals') &
          write(*,*)'DUST_COOLING SNRT_HHE_CIE_METALS: actual H/He NEQ atomic cooling; WSS09 metals remain CIE'
  endif
  call cr_validate(nener,hydro,gpu_hydro,gamma_rad(1),cr_ok)
  if(cr_enabled)then
#ifndef HDF5
     cr_ok=.false.
#endif
     ! No cosmological super-comoving CR source or inflow CR state yet.
     ! simple_boundary is still the "BOUNDARY_PARAMS present" flag here;
     ! nboundary=0 is normalized to periodic below, even with that block.
     if(cosmo.or.nboundary>0)cr_ok=.false.
#ifdef SOLVERmhd
     if(trim(riemann)/='hlld')cr_ok=.false.
#else
     if(trim(riemann)/='hllc'.and.trim(riemann)/='hll'.and.trim(riemann)/='llf')cr_ok=.false.
#endif
     if(trim(outformat)/='hdf5'.or.(nrestart>0.and.trim(informat)/='hdf5'))cr_ok=.false.
     if(sink.or.sink_AGN.or.agn)cr_ok=.false. ! sink/accretion CR partition not yet qualified
     if(delayed_cooling)cr_ok=.false. ! no duplicate delayed-SN reservoir
#ifdef PHASE0_STELLAR_ENRICHMENT
     if(.not.use_channel_resolved_feedback())cr_ok=.false.
#else
     cr_ok=.false.
#endif
     if(cr_sf_support.and.(.not.sf_virial.or.(sf_model/=1.and.sf_model/=2.and.sf_model/=4)))cr_ok=.false.
  endif
  if(.not.cr_ok)then
     if(myid==1)write(*,*)'CR requires NENER=1, gamma_rad=4/3, CPU hydro, advective transport, HDF5, no sink/delayed cooling'
     if(myid==1)write(*,*)'CR SF support requires sf_virial and sf_model=1,2,4; source fractions in [0,1]'
     if(myid==1)write(*,*)'CR comparison requires noncosmo periodic hydro with hllc/hll/llf'
     nml_ok=.false.
  else if(cr_enabled.and.myid==1)then
     write(*,*)'CR trapped-fluid pressure/work/advection enabled; SNII/Ia energy fractions=',cr_sn_fraction,cr_snia_fraction
     write(*,*)'CR SF effective pressure support=',cr_sf_support,'; diffusion, streaming and collisional losses absent'
  endif

  !--------------------------------------------------
  ! Check ind_rsink
  !--------------------------------------------------
  if(ind_rsink<=0.0d0)then
     if(myid==1)write(*,*)'Error in the namelist'
     if(myid==1)write(*,*)'Check ind_rsink'
     nml_ok=.false.
  end if

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
#ifdef SOLVERmhd
     boundary_var(i,:)=0d0
#endif
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
#ifdef SOLVERmhd
     boundary_var(i,6:8)=mhd_seed
     boundary_var(i,nvar+1:nvar+3)=mhd_seed
     boundary_var(i,5)=boundary_var(i,5)+.5d0*sum(mhd_seed**2)
#if NENER>0
     do ivar=1,nener
        boundary_var(i,nhydro+ivar)=prad_bound(i,ivar)/(gamma_rad(ivar)-1d0)
        boundary_var(i,5)=boundary_var(i,5)+boundary_var(i,nhydro+ivar)
     enddo
#endif
     do ivar=nhydro+nener+1,nvar
        boundary_var(i,ivar)=boundary_var(i,1)*var_bound(i,ivar-nhydro-nener)
     enddo
#endif
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do


  !------------------------------------------------------
  !READING THE YIELD TABLE TO SET nvar AUTOMATICALLY
  !------------------------------------------------------
  !Read the yield table                                                                                                              
!  open(33,file=trim(yieldtablefilename),status='old', form='formatted')
!  read(33,'(a8,I)')a1,dum
!  read(33,'(a8,I)')a1,dum
!  read(33,'(a11,I)')a1,nelt

!  nvar = ndim+2+nener+nelt
!  if(metal) nvar=nvar+1
!  if(delayed_cooling) nvar=nvar+1
!  if(sf_virial) nvar=nvar+1
!  if(aton) nvar=nvar+1
!  if(myid==1) write(*,*) 'nvar = ',nvar

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
  inener=nhydro+1
  imetal=nener+nhydro+1
  idelay=imetal
  if(metal)idelay=imetal+1
  ivirial=idelay
  if(delayed_cooling)ivirial=idelay+1
  ixion=ivirial
  if(sf_virial)ixion=ivirial+1
  ! Chemical fields follow every enabled passive variable.  In particular,
  ! they must not overlap the delayed-cooling reservoir at idelay.
  ichem=ixion
  if(aton)ichem=ixion+1
  isgs=ichem
  if(use_sgs)isgs=ichem+1
#ifdef DUST_LIVE
  ! Reserve the dust fields after the full channel-resolved element window.
  ! This prevents collision with ichem:ichem+10 even when that optional
  ! feedback path is disabled in a particular runtime namelist.
  idust=ichem+11
  idust_energy=idust+1
  idust_iron=-1
  idust_pah=-1
  ichimes=-1
  if(dust_chimes_enabled())then
     ichimes=idust+11 ! Reserve the complete existing dust/source field window.
     if(nvar<ichimes+156)then
        if(myid==1)write(*,*)'ERROR: CHIMES live chemistry requires NVAR >= ',ichimes+156
        nml_ok=.false.
     endif
  endif
  idust_species=-1
  if(dust_iron_enabled())then
     idust_iron=ichimes+157
     if(.not.dust_chimes_enabled())idust_iron=idust+11
     if(.not.dust_chimes_enabled().and.dust_relative_motion)nml_ok=.false.
     if(nvar<idust_iron+1)then
        if(myid==1)write(*,*)'ERROR: Fe comparison requires NVAR >= ',idust_iron+1
        nml_ok=.false.
     endif
     if(myid==1)then
        write(*,*)'Fe ELECTRIC-ONLY comparison: grain T<=300 K, primary<=4 eV'
        if(.not.dust_chimes_enabled())write(*,*)'Fe static comparison: CHIMES off; all grain mass reactions disabled'
        if(dust_fe_kinetics)then
           write(*,*)'Fe kinetics: geometric seed growth + Choban26/Nozawa06 thermal sputtering; sticking=',dust_fe_sticking
        else
           write(*,*)'Fe kinetics OFF: fixed mass after injection'
        endif
     endif
  endif
  idust_bins=-1
  if(dust_pah_enabled())then
     if(dust_pah_charged().and.dust_relative_motion)then
        if(myid==1)write(*,*)'ERROR: charged PAH comparisons require coadvected grains'
        call clean_stop
     endif
     idust_pah=ichimes+157+merge(2,0,dust_iron_enabled())
     if(ichimes<1.or.nvar<idust_pah+dust_pah_nstate()-1.or.cosmo)then
        if(myid==1)write(*,*)'ERROR: PAH comparison requires noncosmo CHIMES and NVAR >= ',idust_pah+dust_pah_nstate()-1
        nml_ok=.false.
     endif
     if(myid==1)then
        if(dust_pah_hydrogenated())then
           write(*,*)'PAH H/charge M13-DL01 comparison: 3584 states, H=0--13, <=13.6 eV; no carbon destruction'
           if(dust_pah_h2_enabled())write(*,*) &
                'PAH H2 vacancy-refilling: cation H0--10, M13 bound rate; not a bound on the total H2 effect'
        else if(dust_pah_charged())then
           write(*,*)'PAH FIXED-H charge comparison: 256 states, <=13.6 eV, gas 10--10000 K; NO H loss/destruction'
        else
           write(*,*)'PAH neutral C24H12 comparison: absolute IR; Rayleigh long-wave tail; primary<=4 eV'
        endif
     endif
  endif
  idust_shock=-1;idust_fresh=-1
  if(dust_composition_enabled())then
     idust_species=idust_energy+1
     if(nvar<idust_species+1)then
        if(myid==1)write(*,*)'ERROR: dust composition requires NVAR >= ',idust_species+1
        nml_ok=.false.
     endif
  endif
  if(dust_two_size_enabled())then
     idust_bins=idust_species+2
     if(nvar<idust_bins+3)then
        if(myid==1)write(*,*)'ERROR: two-size dust requires NVAR >= ',idust_bins+3
        nml_ok=.false.
     endif
     if(dust_sn_shocks)then
        idust_shock=idust_bins+4;idust_fresh=idust_shock+1
        if(nvar<idust_fresh+1)then
           if(myid==1)write(*,*)'ERROR: dust SN source transaction requires NVAR >= ',idust_fresh+1
           nml_ok=.false.
        else if(any(var_region(:,idust_shock-imetal+1:idust_fresh-imetal+2)/=0d0))then
           if(myid==1)write(*,*)'ERROR: transient dust SN/fresh IC fields must be zero'
           nml_ok=.false.
        endif
     endif
  endif
  if(.not.hydro)then
     if(myid==1)write(*,*) 'ERROR: DUST_LIVE requires hydro=.true.'
     nml_ok=.false.
  endif
  if(nvar<idust_energy)then
     if(myid==1)then
        write(*,*) 'ERROR: DUST_LIVE requires NVAR >= ',idust_energy
        write(*,*) '  Current NVAR=',nvar,'. Recompile with the live profile or a larger NVAR'
     endif
     call clean_stop
  endif
#endif
  idust_momentum=-1;ndust_phase=0
  if(dust_relative_motion)then
#if defined(DUST_DYNAMICS) && defined(DUST_LIVE) && defined(SNRT_CHIMES) && !defined(SOLVERmhd)
     ndust_phase=4+merge(2,0,dust_iron_enabled())+merge(1,0,dust_pah_enabled())
     idust_momentum=ichimes+157
     if(idust_iron>0)idust_momentum=max(idust_momentum,idust_iron+2)
     if(idust_pah>0)idust_momentum=max(idust_momentum,idust_pah+dust_pah_nstate())
     if(.not.dust_chimes_enabled().or..not.dust_two_size_enabled().or..not.hydro.or.cosmo) nml_ok=.false.
     if(.not.dust_material_composition_enabled().or..not.dust_optics_enabled().or.dust_sublimation_rt_enabled())then
        if(myid==1)write(*,*)'ERROR: relative dust requires DL01/D03; coupled RT sublimation is not admitted'
        nml_ok=.false.
     endif
     if(dust_iron_enabled().and.dust_sublimation_enabled())then
        if(myid==1)write(*,*)'ERROR: relative Fe sublimation needs a C/S/Fe common-enthalpy phase solver'
        nml_ok=.false.
     endif
     dum=0
     if(.not.snrt_dust_contract_loaded)call snrt_dust_contract_load_from_environment(dum)
     if(dum/=0.or.snrt_dust_contract_version/=4.or..not.snrt_dust_contract_exchange_enabled)then
        if(myid==1)write(*,*)'ERROR: relative dust requires version-4 live dust contract with gas exchange'
        nml_ok=.false.
     endif
     if(idust_momentum+3*ndust_phase-1>nvar)then
        if(myid==1)write(*,*)'ERROR: relative dust requires NVAR >= ',idust_momentum+3*ndust_phase-1
        nml_ok=.false.
     endif
     if(ndim/=3.or.nboundary>0.or.pressure_fix.or.interpol_var/=0.or.use_sgs.or.isothermal) nml_ok=.false.
     if(gpu_hydro)then
        if(myid==1)write(*,*)'ERROR: relative dust requires CPU/OMP hydro; RT CUDA may remain enabled'
        nml_ok=.false.
     endif
     if(.not.ieee_is_finite(dust_drag_collision_cross_section_cm2).or.dust_drag_collision_cross_section_cm2<=0)then
        if(myid==1)write(*,*)'ERROR: relative dust needs explicit gas collision cross section for Epstein validity'
        nml_ok=.false.
     endif
     if(sink.or.sink_AGN.or.agn.or.delayed_cooling.or.T2_star>0) nml_ok=.false.
     if(myid==1)write(*,*)'Relative dust: first-order conserved Rusanov, absolute phase momenta, common-T solids'
     if(myid==1)write(*,*)'  idust_momentum, ndust_phase = ',idust_momentum,ndust_phase
#else
     if(myid==1)write(*,*)'ERROR: relative dust requires DUST_DYNAMICS/DUST_LIVE/CHIMES hydro build'
     nml_ok=.false.
#endif
  endif
  if(myid==1) then
     write(*,*) 'Hydro var indices:'
#if NENER>0
     write(*,*) '   inener  = ',inener
#endif
     if(metal)           write(*,*) '   imetal  = ',imetal
     if(delayed_cooling) write(*,*) '   idelay  = ',idelay
     if(sf_virial)       write(*,*) '   ivirial = ',ivirial
     if(aton)            write(*,*) '   ixion   = ',ixion
     write(*,*) '   ichem   = ',ichem
     if(use_sgs)         write(*,*) '   isgs    = ',isgs
#ifdef DUST_LIVE
     write(*,*) '   idust   = ',idust
     write(*,*) '   idust_energy = ',idust_energy
     write(*,*) '   DUST_LIVE fields are primitive dust mass/energy per gas mass in IC input'
#endif
  endif
  ! Last variable is isgs (or ichem if use_sgs=.false.)
  ! Runtime check: make sure NVAR is large enough
  if(use_sgs .and. isgs > nvar) then
     if(myid==1) then
        write(*,*) 'ERROR: use_sgs=T requires NVAR >= ', isgs
        write(*,*) '  Current NVAR=', nvar, '. Recompile with NVAR=', isgs
     end if
     call clean_stop
  end if

#ifdef PHASE0_STELLAR_ENRICHMENT
  if(use_channel_resolved_feedback() .and. metal .and. ichem+10 > nvar) then
     if(myid==1) then
        write(*,*) 'ERROR: channel_resolved feedback requires eleven element fields'
        write(*,*) '  Last element index=',ichem+10,' Current NVAR=',nvar
     end if
     call clean_stop
  end if
#endif

end subroutine read_hydro_params

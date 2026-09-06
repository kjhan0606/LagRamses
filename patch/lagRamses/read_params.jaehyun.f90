subroutine read_params
  use amr_commons
  use snrt_agn_efficiency, only: snrt_agn_rt_requested, snrt_agn_model, snrt_agn_model_reference, &
       snrt_agn_reference_config_ok, snrt_agn_admit_reference
#ifdef SNRT
  use snrt_spectral_contract, only: snrt_spectral_contract_load_from_environment, &
       snrt_spectral_contract_status, snrt_spectral_contract_runtime_allowed, &
       snrt_spectral_contract_error_name, snrt_spectral_contract_error_message, &
       snrt_spectral_contract_source_id
  use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment, &
       snrt_secondary_tables_loaded, snrt_secondary_loaded_source_id, &
       snrt_secondary_loaded_upstream_commit, snrt_thermochemistry_ok, &
       snrt_thermochemistry_error_name, snrt_thermochemistry_error_message
  use snrt_dust_contract, only: snrt_dust_contract_load_from_environment, &
       snrt_dust_contract_loaded, snrt_dust_contract_runtime_allowed, &
       snrt_dust_contract_error_name, snrt_dust_contract_error_message
#endif
  use pm_parameters
  use pm_commons, only: npartmax_auto
  use poisson_parameters
  use hydro_parameters
  use pbh_commons
  use fdm_commons, only: fdm_ghost2, fdm_ghost2_rev
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax,sink_nml_iostat
  integer::agn_model_local,agn_model_min,agn_model_max,agn_contract_error
  logical::agn_snrt_built
#ifdef SNRT
  integer :: snrt_requested_local, snrt_requested_min, snrt_requested_max
  integer :: snrt_thermochemistry_error, snrt_dust_contract_error
  integer :: snrt_dust_contract_env_length
  character(len=1024) :: snrt_dust_contract_env
#endif
  integer::nlevelmax_sink=0
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  character(LEN=15)::accretion_scheme_saved
  character(LEN=15)::agn_acc_method='mass',agn_inj_method='volume'
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok,check_energies=.true.,bondi_use_vrel=.true.
  logical::sink_descent=.false.
  real(dp)::mass_smbh_seed=0d0,mass_merger_vel_check=1d100
  real(dp)::eddington_cap=1d0,AGN_fbk_frac_ener=1d0,AGN_fbk_frac_mom=0d0
  real(dp)::boost_threshold_density=0.1d0,epsilon_kin=1d0
  real(dp)::AGN_fbk_mode_switch_threshold=0.01d0,kin_mass_loading=100d0
  real(dp)::mass_star_AGN=0d0,max_mass_nsc=1d15
  real(dp)::gamma_grad_descent=0d0,fudge_graddescent=1d0
  integer,parameter::tag=1134
  integer::dummy_io,info2
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
!jhshin1
  namelist/run_params/clumpfind,cosmo,pic,sink,sinkprops,lightcone,poisson,hydro,rt,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,remap_thresh,ordering &
       & ,bisec_tol,static,geom,overload,cost_weighting,aton,varcpu_chunk_nfile &
  & ,memory_balance,memory_balance_fast_particles,particle_tree_fast_relink &
  & ,particle_tree_fast_relink_maxlevel &
  & ,mem_weight_grid,mem_weight_part,mem_weight_sink &
       & ,lb_grid_headroom &
  & ,work_weight_grid,work_weight_part,work_weight_sidm_pair &
       & ,time_balance_alpha,lb_timing_interval,timer_report_interval,lb_timing_ema_alpha &
       & ,ksec_level_balance_alpha,ksec_level_min_fraction,ksec_level_bins &
       & ,lb_remap_min_interval,lb_remap_horizon,lb_remap_safety &
       & ,aexp_step_limit &
  & ,jobcontrolfile &
       & ,gpu_hydro,gpu_poisson,gpu_fft,gpu_sink,gpu_scalar,gpu_particle,pm_gpu_min_part,gpu_auto_tune,n_cuda_streams &
       & ,use_fftw &
       & ,mg_merged_rb &
       & ,dump_pk &
       & ,exchange_method &
       & ,use_neutrino &
       & ,sidm &
       & ,de_perturb &
       & ,use_mond &
       & ,use_fR &
       & ,use_nDGP &
       & ,use_symmetron &
       & ,use_dilaton &
       & ,use_galileon &
       & ,scalar_solver_strict &
       & ,use_coupled_de &
       & ,use_quintessence &
       & ,use_kessence &
       & ,use_chaplygin &
       & ,use_rvm &
       & ,use_horndeski &
       & ,use_ede &
       & ,use_sgs &
       & ,use_adm &
       & ,use_fdm &
       & ,use_pbh
  ! Non-standard model namelists (read only when enabled)
  namelist/cpl_params/w0,wa,cs2_de,de_table
  namelist/neutrino_params/omega_nu,neutrino_table
  namelist/fR_params/fR0,fR_n,n_iter_fR,fR_eps
  namelist/nDGP_params/omega_rc,nDGP_branch,n_iter_nDGP,nDGP_eps
  namelist/symmetron_params/a_ssb,beta_symmetron,L_symmetron, &
       & n_iter_symmetron,symmetron_eps
  namelist/dilaton_params/beta_dilaton,L_dilaton,a0_dilaton, &
       & n_iter_dilaton,dilaton_eps
  namelist/galileon_params/galileon_tracker,c2_galileon,c3_galileon, &
       & n_iter_galileon,galileon_eps
  namelist/coupled_de_params/beta_cde,cde_friction,cde_vary_mass
  namelist/quint_params/quint_pot,quint_ic_mode,quint_alpha,quint_lambda,quint_phi_ini
  namelist/kessence_params/kes_x0
  namelist/chaplygin_params/chaplygin_As,chaplygin_alpha
  namelist/rvm_params/rvm_nu
  namelist/horndeski_params/hs_mu0,hs_mass
  namelist/ede_params/omega_ede,z_ede,w_ede
  namelist/sgs_params/sgs_C_prod,sgs_C_diss,sgs_C_smag,sgs_floor,sgs_cap,sgs_e_init,sgs_hydro
  namelist/sidm_params/sidm,sidm_cross_section,sidm_npart_min, &
       & sidm_type,sidm_v0,sidm_power, &
       & sidm_courant,sidm_estimator_diagnostics, &
       & sidm_angular,sidm_epsilon, &
       & sidm_inelastic,sidm_delta,sidm_frac_excited, &
       & sidm_nstates,sidm_energy,sidm_frac_init,sidm_mchi, &
       & sidm_a_type,sidm_a_transition,sidm_sigma_ratio,sidm_a_width, &
       & sidm_fdiss, &
       & sidm_baryon,sidm_baryon_sigma,sidm_baryon_power, &
       & sidm_vrel_max
namelist/adm_params/adm_alpha,adm_mp,adm_me_ratio,adm_xi, &
       & adm_cross_section,adm_T_init,adm_T_floor,adm_adiabatic, &
       & adm_hpm,adm_hpm_gamma,adm_hpm_courant,adm_mol,adm_fH2
  namelist/fdm_params/m_axion,fdm_courant,fdm_nrefine_dB,fdm_hybrid,fdm_split_order,fdm_kinetic, &
       & fdm_outer_ledger, &
       & fdm_dual_soliton_ic,fdm_dual_soliton_rho0,fdm_dual_soliton_rc_box, &
       & fdm_dual_soliton_center_box,fdm_dual_soliton_velocity,fdm_dual_soliton_phase, &
       & fdm_dual_soliton_profile_c, &
       & fdm_cost_mode,fdm_use_hjm,fdm_first_wave_level,fdm_hjm_C1,fdm_hjm_C2,fdm_refine_rho_min, &
       & fdm_nla,fdm_match_aout,fdm_hjm_qp,fdm_qp_c1max,fdm_cn_tol,fdm_refine_matched, &
       & fdm_ghost2,fdm_ghost2_rev
  namelist/pbh_params/pbh_table_file,pbh_fraction,pbh_boost, &
       & pbh_energy_sink,pbh_bkg_warn,pbh_check_provenance, &
       & pbh_mf_model,pbh_spin_model,pbh_hawking_model, &
       & pbh_epsdep_model,pbh_fheat_model,pbh_cr_ivar
  namelist/mond_params/a0_mond,mond_mu_type,mond_type, &
       & n_iter_mond,mond_eps,g_ext_mond
  namelist/cosmo_params/omega_b,omega_m,omega_l,h0
  namelist/output_params/noutput,foutput,fbackup,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump &
       & ,informat,outformat,match_aout
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,ngridmax_auto,npartmax_auto &
       & ,nexpand,boxlen,nsinkmax,nlevel_collapse
  namelist/poisson_params/epsilon,maxiter_fine,restart_phi_warm_start &
       & ,abort_on_mg_nonconvergence,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax
  ! The active lagRamses VPATH selects sink_particle.kjhan.f90 and the
  ! cuRamses pm_parameters module, so the upstream read_sink_params routine
  ! is not linked.  Keep the standard SINK_PARAMS spelling here and connect
  ! the fields used by the active sink implementation below.
  namelist/sink_params/n_sink,rho_sink,d_sink,accretion_scheme,merging_timescale &
       & ,ir_cloud_massive,sink_soft,mass_sink_direct_force,ir_cloud,nsinkmax &
       & ,create_sinks,check_energies,mass_sink_seed,mass_smbh_seed,c_acc &
       & ,nlevelmax_sink,eddington_limit,eddington_cap,acc_sink_boost &
       & ,mass_merger_vel_check,clump_core,verbose_AGN,T2_AGN,T2_min &
       & ,cone_opening,mass_halo_AGN,mass_clump_AGN,mass_star_AGN &
       & ,AGN_fbk_frac_ener,AGN_fbk_frac_mom,T2_max,v_max &
       & ,boost_threshold_density,epsilon_kin,AGN_fbk_mode_switch_threshold &
       & ,kin_mass_loading,bondi_use_vrel,smbh,agn,max_mass_nsc &
       & ,agn_acc_method,agn_inj_method,sink_descent,gamma_grad_descent &
       & ,fudge_graddescent
  namelist/lightcone_params/zmax_cone  &
       & ,elongated_axis_cone1,observer_cone1,minboxr_cone1,maxboxr_cone1 &
       & ,elongated_axis_cone2,observer_cone2,minboxr_cone2,maxboxr_cone2
  namelist/spherical_region_params/spherical_region &
       & ,scenter1,scenter2,scenter3,scenter4,scenter5 &
       & ,sradius1,sradius2,sradius3,sradius4,sradius5
!yonghwi (changed all from previous namelist/movie_params/)
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov &
       & ,proj_axis,movie_vars,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera &
       & ,perspective_camera,smooth_frame,shader_frame &
       & ,tstart_theta_camera,tstart_phi_camera &
       & ,tend_theta_camera,tend_phi_camera,dist_camera,ddist_camera
!yonghwi

!jhshin2
  ! MPI initialization
#ifndef WITHOUTMPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
  !--------------------------------------------------
  ! Advertise cuRAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'  _/_/_/   _/    _/   _/_/_/         _/_/     _/    _/    _/_/_/    _/_/_/_/    _/_/_/ '
  write(*,*)' _/        _/    _/    _/    _/     _/  _/    _/_/_/_/   _/    _/   _/         _/    _/'
  write(*,*)'_/         _/    _/    _/    _/    _/    _/   _/ _/ _/   _/         _/         _/      '
  write(*,*)'_/         _/    _/    _/_/_/     _/_/_/_/    _/    _/     _/_/     _/_/_/       _/_/  '
  write(*,*)'_/         _/    _/    _/    _/   _/    _/    _/    _/         _/   _/               _/'
  write(*,*)' _/        _/    _/    _/    _/   _/    _/    _/    _/   _/    _/   _/         _/    _/'
  write(*,*)'  _/_/_/    _/_/_/     _/    _/   _/    _/    _/    _/    _/_/_/    _/_/_/_/    _/_/_/ '
  write(*,*)'                              Version 3.0                                             '
  write(*,*)'             written by Romain Teyssier (University of Zurich)                        '
  write(*,*)'             GPU acceleration by Juhan Kim (KIAS)                                       '
  write(*,*)'                     (c) CEA 1999-2007, UZH 2008-2014                                 '
  write(*,*)'        GPU & optimization by Juhan Kim (KIAS) 2026                                   '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<8)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif
  
  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token                                                                                                                                                                                
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  if(aexp_step_limit<=0.0d0)then
     if(myid==1)write(*,*)'ERROR: aexp_step_limit must be positive'
     call clean_stop
  endif
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=lightcone_params,END=84)
84 continue
  rewind(1)
  read(1,NML=spherical_region_params,END=83)
83 continue
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue
  if(myid==1) write(*,'(A,I0)') ' Fine MG max iterations = ',maxiter_fine
  if(myid==1 .and. restart_phi_warm_start) write(*,'(A)') &
       ' Restart phi policy = valid-marker warm start (explicit opt-in)'
  if(myid==1 .and. .not.restart_phi_warm_start) write(*,'(A)') &
       ' Restart phi policy = predictor (default)'
  if(myid==1) write(*,'(A,L1)') ' Abort on fine MG nonconvergence = ', &
       abort_on_mg_nonconvergence
  rewind(1)
  read(1,NML=cosmo_params,END=80)
80 continue
  rewind(1)
  read(1,NML=cpl_params,END=79)
79 continue
  rewind(1)
  read(1,NML=neutrino_params,END=78)
78 continue
  rewind(1)
  read(1,NML=sidm_params,END=77)
77 continue
  ! aDM parameters
  if(use_adm) then
     rewind(1)
     read(1,NML=adm_params,END=76)
76   continue
  end if
  ! FDM parameters
  if(use_fdm) then
     rewind(1)
     read(1,NML=fdm_params,END=176)
176  continue
  end if
  ! Evaporating-PBH parameters
  if(use_pbh) then
     rewind(1)
     read(1,NML=pbh_params,END=175)
175  continue
  end if
  rewind(1)
  read(1,NML=mond_params,END=75)
75 continue
  ! f(R) parameters
  if(use_fR) then
     rewind(1)
     read(1,NML=fR_params,END=74)
74   continue
  end if
  ! nDGP parameters
  if(use_nDGP) then
     rewind(1)
     read(1,NML=nDGP_params,END=73)
73   continue
  end if
  ! Symmetron parameters
  if(use_symmetron) then
     rewind(1)
     read(1,NML=symmetron_params,END=72)
72   continue
  end if
  ! Dilaton parameters
  if(use_dilaton) then
     rewind(1)
     read(1,NML=dilaton_params,END=71)
71   continue
  end if
  ! Galileon parameters
  if(use_galileon) then
     rewind(1)
     read(1,NML=galileon_params,END=70)
70   continue
  end if
  ! Coupled DE parameters
  if(use_coupled_de) then
     rewind(1)
     read(1,NML=coupled_de_params,END=69)
69   continue
  end if
  ! Early DE parameters
  if(use_ede) then
     rewind(1)
     read(1,NML=ede_params,END=68)
68   continue
  end if
  ! Quintessence parameters
  if(use_quintessence) then
     rewind(1)
     read(1,NML=quint_params,END=561)
561  continue
  end if
  ! k-essence parameters
  if(use_kessence) then
     rewind(1)
     read(1,NML=kessence_params,END=562)
562  continue
  end if
  ! Chaplygin gas parameters
  if(use_chaplygin) then
     rewind(1)
     read(1,NML=chaplygin_params,END=564)
564  continue
  end if
  ! Running vacuum parameters
  if(use_rvm) then
     rewind(1)
     read(1,NML=rvm_params,END=565)
565  continue
  end if
  ! Horndeski parameters
  if(use_horndeski) then
     rewind(1)
     read(1,NML=horndeski_params,END=563)
563  continue
  end if
  ! SGS turbulence parameters
  if(use_sgs) then
     rewind(1)
     read(1,NML=sgs_params,END=67)
67   continue
  end if

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (myid==1 .and. narg == 2) then
     CALL getarg(2,cmdarg)
     read(cmdarg,*) nrestart
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! GPU acceleration: disable if not compiled with USE_CUDA
  !-------------------------------------------------
#ifndef HYDRO_CUDA
  ! Most gpu_* switches default on, so a build without USE_CUDA turns them
  ! off silently; warning here would fire on every CPU run.
  gpu_hydro = .false.
  gpu_poisson = .false.
  gpu_fft = .false.
  gpu_sink = .false.
  gpu_scalar = .false.
  gpu_particle = .false.
#else
  if(myid==1 .and. (gpu_hydro .or. gpu_poisson .or. gpu_fft .or. gpu_sink &
       & .or. gpu_scalar .or. gpu_particle)) then
     write(*,'(A,L1,A,L1,A,L1,A,L1,A,L1,A,L1,A,I0)') &
          ' GPU acceleration: hydro=',gpu_hydro, &
          ' poisson=',gpu_poisson,' fft=',gpu_fft,' sink=',gpu_sink, &
          ' scalar=',gpu_scalar,' particle=',gpu_particle, &
          ' streams=',n_cuda_streams
  end if
#endif

  !-------------------------------------------------
  ! FFTW3 CPU Poisson solver: disable if not compiled with USE_FFTW
  !-------------------------------------------------
#ifndef USE_FFTW
  if(use_fftw) then
     if(myid==1) write(*,*) 'WARNING: use_fftw ignored (not compiled with USE_FFTW)'
     use_fftw = .false.
  end if
#else
  if(myid==1 .and. use_fftw) then
     write(*,'(A)') ' FFTW3 CPU direct Poisson solver enabled (use_fftw=T)'
  end if
#endif

  !-------------------------------------------------
  ! Multigrid merged red/black smoother diagnostic
  !-------------------------------------------------
  if(myid==1) then
     if(mg_merged_rb) then
        write(*,'(A)') ' Multigrid GPU smoother in merged red/black mode (mg_merged_rb=T)'
     else
        write(*,'(A)') ' Multigrid GPU smoother in strict red/black mode (mg_merged_rb=F, bitwise NCPU independent)'
     end if
  end if

  !-------------------------------------------------
  ! Auto-compute mem_weight_grid from nvar if sentinel (0)
  ! Per cell (×twotondim per grid):
  !   Hydro:   2*nvar*8 (uold+unew)
  !   Topo:    5*4 (son,flag1,flag2,cpu_map,cpu_map2) + sizeof(qdp) (hilbert_key)
  !   PFix:    2*8 (enew,divu)  [pressure_fix, default for cosmo]
  !   Poisson: 7*8 (rho,rho_star,phi,phi_old,f*3)
  !   FDM:     2*8 (psi_re,psi_im)  [use_fdm only]
  ! Per grid:
  !   AMR:  3*8+3*4+6*4 (xg,father/next/prev,nbor) = 48
  !   Part: 3*4 (headp,tailp,numbp) = 12
  !-------------------------------------------------
  if(fdm_cost_mode < 0 .or. fdm_cost_mode > 1) then
     if(myid==1) write(*,'(A,I0)') ' ERROR: fdm_cost_mode must be 0 or 1; got ',fdm_cost_mode
     call clean_stop
  end if
  ! Backward-compatible alias for FDM input files written before the generic
  ! runtime load-balance controls existed.  Mode 0 leaves RUN_PARAMS intact;
  ! mode 1 selects the timed work model for the initial mode only.
  if(use_fdm .and. fdm_cost_mode==1) then
     memory_balance=.false.
     cost_weighting=.true.
     if(time_balance_alpha<=0d0) time_balance_alpha=0.3d0
     if(myid==1) write(*,'(A)') &
          ' FDM fdm_cost_mode=1 (deprecated alias): selecting timed work balance'
  end if
  ! Compute this even when the initial mode is work-based: job control may
  ! switch to memory mode later without rereading the namelist.
  if(mem_weight_grid <= 0) then
     mem_weight_grid = twotondim * (2*nvar*8 + 20 + storage_size(0.0_qdp)/8 + 16 + 56) + 48 + 12
     if(use_fdm) mem_weight_grid = mem_weight_grid + twotondim * 2 * 8
     if(myid==1) write(*,'(A,I6,A,I3,A)') &
          ' Memory balance: mem_weight_grid=',mem_weight_grid,' (nvar=',nvar,')'
  end if
  if(myid==1) write(*,'(A,F6.3,A)') &
       ' Load-balance grid headroom=',lb_grid_headroom,' x ngridmax'
  if(myid==1 .and. memory_balance)then
     if(memory_balance_fast_particles)then
        write(*,'(A)') ' Memory balance particle placement: fast grid-total mode'
     else
        write(*,'(A)') ' Memory balance particle placement: exact linked-list mode'
     end if
  end if
  particle_tree_fast_relink_maxlevel=max(0,min(MAXLEVEL-1, &
       particle_tree_fast_relink_maxlevel))
  if(myid==1 .and. pic)then
     if(particle_tree_fast_relink)then
        write(*,'(A,I0)') ' Particle-tree rebuild: one-pass relink through level ', &
             particle_tree_fast_relink_maxlevel
     else
        write(*,'(A)') ' Particle-tree rebuild: exact remove/add path'
     end if
  end if

  !-------------------------------------------------
  ! Work-balance model validation
  !-------------------------------------------------
  work_weight_grid=max(0,work_weight_grid)
  work_weight_part=max(0,work_weight_part)
  work_weight_sidm_pair=max(0,work_weight_sidm_pair)
  ksec_level_balance_alpha=max(0d0,min(1d0,ksec_level_balance_alpha))
  ksec_level_min_fraction=max(0d0,min(1d0,ksec_level_min_fraction))
  ksec_level_bins=max(64,ksec_level_bins)
  lb_timing_interval=max(0,lb_timing_interval)
  timer_report_interval=max(0,timer_report_interval)
  lb_timing_ema_alpha=max(0d0,min(1d0,lb_timing_ema_alpha))
  lb_remap_min_interval=max(0,lb_remap_min_interval)
  lb_remap_horizon=max(1,lb_remap_horizon)
  lb_remap_safety=max(0d0,lb_remap_safety)
  if((.not.memory_balance).and.myid==1)then
     write(*,'(A,3(I0,1X))') ' Work balance weights grid/part/SIDM-pair: ', &
          work_weight_grid,work_weight_part,merge(work_weight_sidm_pair,0,sidm)
     write(*,'(A,I0,A,F5.2,A,I0,A,I0,A,F5.2)') &
          ' Work timing: every ',lb_timing_interval,' steps, EMA=', &
          lb_timing_ema_alpha,', remap min/horizon=',lb_remap_min_interval, &
          '/',lb_remap_horizon,', safety=',lb_remap_safety
  end if
  if(myid==1 .and. timer_report_interval>0) write(*,'(A,I0,A)') &
       ' Phase timer report: every ',timer_report_interval,' coarse steps'

  !-------------------------------------------------
  ! Exchange method auto-tune
  !-------------------------------------------------
  if(ordering=='ksection' .and. myid==1) then
     write(*,'(A,A)') ' Exchange method: ', trim(exchange_method)
     write(*,'(A,F5.2,A,F6.3,A,I0)') &
          ' K-section level balance: alpha=',ksec_level_balance_alpha, &
          ', min fraction=',ksec_level_min_fraction, &
          ', bins=',ksec_level_bins
  end if

  !-------------------------------------------------
  ! DE perturbation (CPL)
  ! Three modes:
  !   1) de_table provided → table-based linear response (any cs2_de)
  !   2) no de_table, cs2_de>0 → kappa2/alpha quasi-static method
  !   3) no de_table, cs2_de<=0 → unsupported, disable
  !-------------------------------------------------
  if(de_perturb) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'WARNING: de_perturb=T but not cosmo run, disabling'
        de_perturb = .false.
     else if(len_trim(de_table) > 0) then
        ! Table-based linear response (works for any cs2_de)
        if(myid==1) then
           write(*,'(A,A)') ' DE perturbation (table): ', trim(de_table)
           write(*,'(A,ES10.3,A,F6.3,A,F6.3)') &
                '   cs2_de=', cs2_de, ' w0=', w0, ' wa=', wa
        end if
     else if(use_quintessence .or. use_kessence) then
        ! Scalar-field DE provides model w(a), cs2(a) internally
        if(myid==1) write(*,'(A)') &
             ' DE perturbation (kappa2/alpha): cs2(a), w(a) from scalar-field DE model'
     else if(cs2_de <= 0.0d0) then
        if(myid==1) write(*,*) 'WARNING: de_perturb=T, no de_table, cs2_de<=0 -> disabling'
        de_perturb = .false.
     else
        ! Fallback: kappa2/alpha quasi-static method
        if(myid==1) then
           write(*,'(A,ES10.3,A,F6.3,A,F6.3)') &
                ' DE perturbation (kappa2/alpha): cs2_de=', cs2_de, &
                ' w0=', w0, ' wa=', wa
        end if
     end if
  end if

  !-------------------------------------------------
  ! Chaplygin gas / running-vacuum smooth dark energy
  !-------------------------------------------------
  if(use_chaplygin .or. use_rvm) then
     if(use_chaplygin .and. use_rvm) then
        if(myid==1) write(*,*) 'ERROR: use_chaplygin and use_rvm are mutually exclusive'
        call clean_stop
     end if
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: use_chaplygin/use_rvm require cosmo=.true.'
        call clean_stop
     end if
     if(use_quintessence .or. use_kessence .or. use_coupled_de .or. use_ede &
          & .or. use_galileon) then
        if(myid==1) write(*,*) &
             & 'ERROR: Chaplygin/running vacuum redefine the DE background; disable ', &
             & 'quintessence, k-essence, coupled DE, EDE and the Galileon'
        call clean_stop
     end if
     if(w0 /= -1.0d0 .or. wa /= 0.0d0) then
        if(myid==1) write(*,*) &
             & 'ERROR: Chaplygin/running vacuum set their own w(a); leave w0=-1, wa=0'
        call clean_stop
     end if
     if(use_fR .or. use_nDGP .or. use_symmetron .or. use_dilaton .or. use_horndeski) then
        if(myid==1) write(*,*) &
             & 'ERROR: screened/Horndeski gravity assumes a LCDM background; ', &
             & 'incompatible with Chaplygin/running vacuum'
        call clean_stop
     end if
  end if
  if(use_chaplygin) then
     if(chaplygin_As <= 0.0d0 .or. chaplygin_As >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: chaplygin_As must lie in (0,1), got', chaplygin_As
        call clean_stop
     end if
     if(chaplygin_alpha < 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: chaplygin_alpha must be >= 0, got', chaplygin_alpha
        call clean_stop
     end if
     if(myid==1) write(*,'(A,F7.4,A,ES10.3)') &
          & ' Generalized Chaplygin gas DE: A_s=', chaplygin_As, ' alpha=', chaplygin_alpha
  end if
  if(use_rvm) then
     if(abs(rvm_nu) >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: rvm_nu must satisfy |nu|<1, got', rvm_nu
        call clean_stop
     end if
     if(myid==1) write(*,'(A,ES12.4)') ' Running vacuum Lambda(H^2)=c0+nu*H^2: nu=', rvm_nu
  end if

  !-------------------------------------------------
  ! Neutrino linear response
  !-------------------------------------------------
  if(use_neutrino) then
     if(omega_nu <= 0.0d0) then
        if(myid==1) write(*,*) 'WARNING: use_neutrino=T but omega_nu<=0, disabling'
        use_neutrino = .false.
     else if(len_trim(neutrino_table) == 0) then
        if(myid==1) write(*,*) 'WARNING: use_neutrino=T but neutrino_table not set, disabling'
        use_neutrino = .false.
     else
        if(myid==1) then
           write(*,'(A,F7.4,A,F7.4)') &
                ' Neutrino linear response: omega_nu=', omega_nu, &
                ' omega_cb=', omega_m - omega_nu
           write(*,'(A,A)') '   table: ', trim(neutrino_table)
        end if
     end if
  end if

  !-------------------------------------------------
  ! SIDM (Self-Interacting Dark Matter) scattering
  !-------------------------------------------------
  if(sidm) then
     if(.not. pic) then
        if(myid==1) write(*,*) 'ERROR: sidm=T requires pic=T'
        call clean_stop
     end if
     if(sidm_cross_section <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm=T but sidm_cross_section<=0'
        call clean_stop
     end if
     if(sidm_inelastic .and. sidm_delta <= 0.0d0 &
          .and. sidm_energy(1) <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_inelastic=T but no energy splitting set'
        call clean_stop
     end if
     if(sidm_inelastic .and. sidm_mchi <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_inelastic=T but sidm_mchi<=0 (DM mass in GeV)'
        call clean_stop
     end if
     ! Auto-populate multi-state arrays from 2-state shortcut
     if(sidm_inelastic .and. sidm_nstates==2 &
          .and. sidm_energy(1)==0.0d0 .and. sidm_delta>0.0d0) then
        sidm_energy(0) = 0.0d0
        sidm_energy(1) = sidm_delta
     end if
     if(sidm_inelastic .and. sidm_nstates==2 &
          .and. sidm_frac_init(0)==0.0d0 .and. sidm_frac_init(1)==0.0d0 &
          .and. sidm_frac_excited>0.0d0) then
        sidm_frac_init(0) = 1.0d0 - sidm_frac_excited
        sidm_frac_init(1) = sidm_frac_excited
     end if
     if(myid==1) then
        write(*,'(A,ES10.3,A)') ' SIDM enabled: sigma/m=', sidm_cross_section, ' cm^2/g'
        write(*,'(A,A)')        '   cross-section type: ', trim(sidm_type)
        if(trim(sidm_type) /= 'constant') then
           write(*,'(A,F8.1,A)') '   v0=', sidm_v0, ' km/s'
           if(trim(sidm_type) == 'power_law') &
                write(*,'(A,F6.2)') '   power=', sidm_power
        end if
        write(*,'(A,I4)')        '   npart_min=', sidm_npart_min
        write(*,'(A,F5.2)')      '   courant=', sidm_courant
        write(*,'(A,ES10.3,A)')  '   vrel_max=', sidm_vrel_max, ' cm/s'
        write(*,'(A,A)')         '   angular: ', trim(sidm_angular)
        if(trim(sidm_angular) == 'rutherford') &
             write(*,'(A,ES10.3)') '   epsilon=', sidm_epsilon
        if(sidm_inelastic) then
           write(*,'(A,I3)')       '   iSIDM: nstates=', sidm_nstates
           write(*,'(A,10ES10.3)') '   energies [keV]=', &
                sidm_energy(0:sidm_nstates-1)
           write(*,'(A,10F7.3)')   '   frac_init=', &
                sidm_frac_init(0:sidm_nstates-1)
        end if
        if(trim(sidm_a_type) /= 'none') then
           write(*,'(A,A)')         '   Phase transition type: ', trim(sidm_a_type)
           write(*,'(A,F8.4)')      '   a_transition=', sidm_a_transition
           write(*,'(A,F8.2)')      '   sigma_ratio=', sidm_sigma_ratio
           if(trim(sidm_a_type) == 'sigmoid') &
                write(*,'(A,F8.4)') '   a_width=', sidm_a_width
        end if
        if(sidm_fdiss > 0.0d0) then
           write(*,'(A,F6.3)')     '   dSIDM: fdiss=', sidm_fdiss
        end if
        if(sidm_baryon) then
           write(*,'(A,ES10.3,A,F5.1)') &
                '   IDM: sigma_DM-b/m=', sidm_baryon_sigma, &
                ' cm^2/g  power=', sidm_baryon_power
        end if
     end if
     if(sidm_fdiss < 0.0d0 .or. sidm_fdiss >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_fdiss must be in [0,1)'
        call clean_stop
     end if
     if(sidm_baryon .and. .not.hydro) then
        if(myid==1) write(*,*) 'ERROR: sidm_baryon requires hydro=T'
        call clean_stop
     end if
     if(sidm_baryon .and. sidm_baryon_sigma <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: sidm_baryon_sigma must be > 0'
        call clean_stop
     end if
  end if

  !-------------------------------------------------
  ! Atomic Dark Matter (aDM)
  !-------------------------------------------------
  if(use_adm) then
     if(.not. sidm) then
        if(myid==1) write(*,*) 'ERROR: use_adm=T requires sidm=T'
        call clean_stop
     end if
     if(adm_alpha <= 0.0d0 .or. adm_alpha >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_alpha must be in (0,1)'
        call clean_stop
     end if
     if(adm_mp <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_mp must be > 0'
        call clean_stop
     end if
     if(adm_me_ratio <= 0.0d0 .or. adm_me_ratio >= 1.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_me_ratio must be in (0,1)'
        call clean_stop
     end if
     if(adm_xi <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_xi must be > 0'
        call clean_stop
     end if
     if(adm_T_init < adm_T_floor) then
        if(myid==1) write(*,*) 'ERROR: adm_T_init must be >= adm_T_floor'
        call clean_stop
     end if
     if(adm_T_floor <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: adm_T_floor must be positive'
        call clean_stop
     end if
     if(adm_hpm) then
        if(.not.poisson) then
           if(myid==1) write(*,*) 'ERROR: adm_hpm requires poisson=T'
           call clean_stop
        end if
        if(adm_hpm_gamma <= 1.0d0) then
           if(myid==1) write(*,*) 'ERROR: adm_hpm_gamma must exceed 1'
           call clean_stop
        end if
        if(adm_hpm_courant <= 0.0d0 .or. adm_hpm_courant > 1.0d0) then
           if(myid==1) write(*,*) 'ERROR: adm_hpm_courant must be in (0,1]'
           call clean_stop
        end if
        ! The first HPM implementation applies the mesh kick through f(),
        ! which is shared by the collisionless particle population.  Do not
        ! silently accelerate stellar or sink particles as ADM.
        if(star .or. sink) then
           if(myid==1) write(*,*) 'ERROR: adm_hpm currently requires star=F and sink=F'
           call clean_stop
        end if
     end if
     if(myid==1) then
        write(*,'(A)') ' Atomic Dark Matter (aDM) enabled:'
        write(*,'(A,ES10.3)') '   alpha_D  =', adm_alpha
        write(*,'(A,F8.3,A)') '   m_p''     =', adm_mp, ' GeV'
        write(*,'(A,ES10.3)') '   m_e''/m_p''=', adm_me_ratio
        write(*,'(A,F6.3)')   '   xi       =', adm_xi
        write(*,'(A,ES10.3,A)') '   sigma/m  =', adm_cross_section, ' cm^2/g'
        write(*,'(A,ES10.3,A)') '   T_init   =', adm_T_init, ' K (new runs)'
        write(*,'(A,ES10.3,A,L1)') '   T_floor  =', adm_T_floor, &
             & ' K; adiabatic=',adm_adiabatic
        if(adm_hpm) write(*,'(A,F6.3,A,F6.3)') &
             & '   HPM pressure ON: gamma=',adm_hpm_gamma, &
             & '; CFL=',adm_hpm_courant
        if(adm_mol) then
           write(*,'(A,ES10.3)') '   dark H2 cooling ON, fH2 =', adm_fH2
        else
           write(*,'(A)')        '   dark H2 cooling OFF'
        end if
     end if
  end if

  !-------------------------------------------------
  ! Fuzzy Dark Matter (FDM)
  !-------------------------------------------------
  if(use_fdm) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_fdm=T requires poisson=T'
        call clean_stop
     end if
     if(m_axion <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: m_axion must be > 0'
        call clean_stop
     end if
     if(sidm) then
        if(myid==1) write(*,*) 'ERROR: FDM and SIDM are mutually exclusive'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)')         ' Fuzzy Dark Matter (FDM) enabled:'
        write(*,'(A,ES10.3,A)') '   m_axion  =', m_axion, ' eV'
        write(*,'(A,F5.2)')    '   courant  =', fdm_courant
        write(*,'(A,I3)')      '   nrefine_dB=', fdm_nrefine_dB
        write(*,'(A,L1)')      '   hybrid   =', fdm_hybrid
        write(*,'(A,L1)')      '   outer wave provenance=', fdm_outer_ledger
        write(*,'(A,L1)')      '   dual soliton IC=', fdm_dual_soliton_ic
        if(fdm_dual_soliton_ic)then
           write(*,'(A,2(1X,ES10.3))') '     rho0=',fdm_dual_soliton_rho0
           write(*,'(A,2(1X,ES10.3))') '     rc/box=',fdm_dual_soliton_rc_box
        end if
        write(*,'(A,I2,A)')    '   split_order=', fdm_split_order, &
             & merge(' (Strang DKD)   ', merge(' (Yoshida 4th)  ', ' (UNKNOWN->DKD) ', fdm_split_order==4), fdm_split_order==2)
        write(*,'(A,I2,A)')    '   kinetic   =', fdm_kinetic, &
             & merge(' (explicit subcyc)', ' (Crank-Nicolson) ', fdm_kinetic==0)
        write(*,'(A,I2)')      '   cost_mode (deprecated alias)=',fdm_cost_mode
        if(memory_balance) then
           write(*,'(A)')      '   effective LB=memory'
        else if(time_balance_alpha>0d0) then
           write(*,'(A,F5.2,A)') '   effective LB=timed work (alpha=',time_balance_alpha,')'
        else
           write(*,'(A)')      '   effective LB=work'
        end if
        write(*,'(A,L1)')      '   use_hjm   =', fdm_use_hjm
        write(*,'(A,L1)')      '   refine matched=', fdm_refine_matched
        write(*,'(A,L1)')      '   match aout=', fdm_match_aout
        write(*,'(A,ES10.3)')  '   CN tolerance=', fdm_cn_tol
        if(fdm_use_hjm) then
           if(fdm_kinetic /= 1) then
              fdm_kinetic = 1
              write(*,'(A)')      '   HJM: forcing fdm_kinetic=1 (Crank-Nicolson) for wave levels'
           end if
           write(*,'(A)')         '   HJM hybrid mode:'
           write(*,'(A,I3)')      '     first_wave_level=', fdm_first_wave_level
           write(*,'(A,ES10.3)')  '     C1 thresh  =', fdm_hjm_C1
           write(*,'(A,ES10.3)')  '     C2 thresh  =', fdm_hjm_C2
           write(*,'(A,ES10.3)')  '     rho_min(dB)=', fdm_refine_rho_min
           write(*,'(A,L1)')      '     qp on fluid=', fdm_hjm_qp
           if(fdm_hjm_qp) write(*,'(A,ES10.3)')  '     qp C1 gate =', fdm_qp_c1max
        end if
     end if
  end if

  !-------------------------------------------------
  ! MOND (Modified Newtonian Dynamics)
  !-------------------------------------------------
  if(use_mond) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_mond=T requires poisson=T'
        call clean_stop
     end if
     if(a0_mond <= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: use_mond=T but a0_mond<=0'
        call clean_stop
     end if
     if(mond_mu_type < 1 .or. mond_mu_type > 2) then
        if(myid==1) write(*,*) 'ERROR: mond_mu_type must be 1 (simple) or 2 (standard)'
        call clean_stop
     end if
     if(mond_type < 0 .or. mond_type > 2) then
        if(myid==1) write(*,*) 'ERROR: mond_type must be 0, 1 (QUMOND), or 2 (AQUAL)'
        call clean_stop
     end if
     if(mond_type == 2) then
        if(n_iter_mond < 1) then
           if(myid==1) write(*,*) 'ERROR: n_iter_mond must be >= 1'
           call clean_stop
        end if
        if(mond_eps <= 0d0) then
           if(myid==1) write(*,*) 'ERROR: mond_eps must be > 0'
           call clean_stop
        end if
     end if
     ! Warn if DM particles likely present
     if(pic) then
        if(myid==1) then
           write(*,'(A)') ' WARNING: use_mond=T with pic=T (DM particles likely present)'
           write(*,'(A)') '   MOND replaces DM — running both gives excess gravity.'
           write(*,'(A)') '   Use DM-free ICs or set pic=.false. for pure MOND.'
        end if
     end if
     if(myid==1) then
        write(*,'(A)') ' MOND (QUMOND) enabled:'
        write(*,'(A,ES12.4,A)') '   a0 = ', a0_mond, ' cm/s^2'
        if(mond_mu_type==1) then
           write(*,'(A)') '   mu-function: simple [mu=x/(1+x)]'
        else
           write(*,'(A)') '   mu-function: standard [mu=x/sqrt(1+x^2)]'
        end if
        if(mond_type==0) then
           write(*,'(A)') '   mode: algebraic QUMOND (Phase 0)'
        else if(mond_type==1) then
           write(*,'(A)') '   mode: full QUMOND with phantom density (Phase 1)'
        else
           write(*,'(A,I3,A,ES10.3)') &
                '   mode: AQUAL iterative (Phase 2), max_iter=', n_iter_mond, &
                ' eps=', mond_eps
        end if
        if(g_ext_mond(1)**2+g_ext_mond(2)**2+g_ext_mond(3)**2 > 0d0) then
           write(*,'(A,3ES12.4,A)') '   g_ext = (', g_ext_mond, ') cm/s^2'
        end if
     end if
  end if

  !-------------------------------------------------
  ! f(R) Hu-Sawicki gravity
  !-------------------------------------------------
  if(use_fR) then
     ! Mutual exclusion checks
     if(use_nDGP) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both f(R) and nDGP simultaneously'
        call clean_stop
     end if
     if(use_mond) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both f(R) and MOND simultaneously'
        call clean_stop
     end if
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: f(R) gravity requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_fR=T requires poisson=T'
        call clean_stop
     end if
     if(fR0 >= 0d0) then
        if(myid==1) write(*,*) 'ERROR: fR0 must be negative'
        call clean_stop
     end if
     if(fR_n < 1) then
        if(myid==1) write(*,*) 'ERROR: fR_n must be >= 1'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' f(R) Hu-Sawicki gravity enabled'
        write(*,'(A,ES10.3,A,I2)') '   fR0=', fR0, '  n=', fR_n
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_fR, '  eps=', fR_eps
     end if
  end if

  !-------------------------------------------------
  ! nDGP gravity
  !-------------------------------------------------
  if(use_nDGP) then
     if(use_mond) then
        if(myid==1) write(*,*) 'ERROR: Cannot use both nDGP and MOND simultaneously'
        call clean_stop
     end if
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: nDGP gravity requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_nDGP=T requires poisson=T'
        call clean_stop
     end if
     if(omega_rc <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: omega_rc must be > 0'
        call clean_stop
     end if
     if(abs(nDGP_branch) /= 1) then
        if(myid==1) write(*,*) 'ERROR: nDGP_branch must be 1 or -1'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' nDGP gravity enabled'
        write(*,'(A,ES10.3,A,I2)') '   omega_rc=', omega_rc, '  branch=', nDGP_branch
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_nDGP, '  eps=', nDGP_eps
     end if
  end if

  !-------------------------------------------------
  ! Scalar field mutual exclusion check
  ! fR, nDGP, MOND, symmetron, dilaton, galileon: at most 1
  !-------------------------------------------------
  i = 0
  if(use_fR)        i = i + 1
  if(use_nDGP)      i = i + 1
  if(use_mond)      i = i + 1
  if(use_symmetron)  i = i + 1
  if(use_dilaton)    i = i + 1
  if(use_galileon)   i = i + 1
  if(i > 1) then
     if(myid==1) write(*,*) 'ERROR: Only one scalar field gravity model allowed at a time'
     if(myid==1) write(*,*) '  (fR, nDGP, MOND, symmetron, dilaton, galileon)'
     call clean_stop
  end if

  !-------------------------------------------------
  ! MG solver backgrounds assume LCDM; DE-boost paths break
  ! the QUMOND phantom cancellation in the FFT solver
  !-------------------------------------------------
  if(use_fR .or. use_nDGP .or. use_symmetron .or. use_dilaton .or. use_galileon) then
     if(use_quintessence .or. use_kessence) then
        if(myid==1) write(*,*) 'ERROR: scalar-field DE cannot be combined with the MG solvers (LCDM background assumed)'
        call clean_stop
     end if
     if(w0 /= -1.0d0 .or. wa /= 0.0d0 .or. use_ede) then
        if(myid==1) write(*,*) 'WARNING: MG solver backgrounds assume LCDM; CPL/EDE combination is inconsistent'
     end if
  end if
  if(use_mond .and. (de_perturb .or. use_horndeski .or. use_coupled_de)) then
     if(myid==1) write(*,*) 'ERROR: use_mond cannot be combined with de_perturb/use_horndeski/use_coupled_de'
     if(myid==1) write(*,*) '  (DE source boosts break the QUMOND phantom-density cancellation)'
     call clean_stop
  end if
  if(use_galileon .and. .not. galileon_tracker .and. myid==1) then
     write(*,*) 'WARNING: cubic Galileon LEGACY template (galileon_tracker=F):'
     write(*,*) '  simplified non-tracker coefficients, NOT Barreira+13; experimental'
  end if
  if(use_galileon .and. galileon_tracker .and. myid==1) then
     write(*,'(A)') ' Cubic Galileon: Barreira+13 tracker (parameter-free)'
     write(*,'(A,F8.4,A,F8.4)') '   xi=sqrt(6(1-Om))=', sqrt(6d0*(1d0-omega_m)), &
          & '  Geff/G(a=1)=', 1d0+sqrt(6d0*(1d0-omega_m))**3 &
          & /(18d0*(-(sqrt(6d0*(1d0-omega_m))/3d0) &
          & *(2d0*(-1.5d0*omega_m/(2d0-omega_m))-1d0+(1d0-omega_m))))
  end if
  if(use_galileon .and. .not. galileon_tracker .and. c2_galileon == 0d0) then
     if(myid==1) write(*,*) 'ERROR: c2_galileon must be nonzero'
     call clean_stop
  end if

  !-------------------------------------------------
  ! Symmetron gravity
  !-------------------------------------------------
  if(use_symmetron) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Symmetron requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_symmetron=T requires poisson=T'
        call clean_stop
     end if
     if(a_ssb <= 0d0 .or. a_ssb >= 1d0) then
        if(myid==1) write(*,*) 'ERROR: a_ssb must be in (0,1)'
        call clean_stop
     end if
     if(L_symmetron <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: L_symmetron must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Symmetron gravity enabled'
        write(*,'(A,F6.3,A,F6.3,A,F8.3,A)') &
             '   a_ssb=', a_ssb, '  beta=', beta_symmetron, '  L=', L_symmetron, ' Mpc/h'
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_symmetron, '  eps=', symmetron_eps
     end if
  end if

  !-------------------------------------------------
  ! Dilaton gravity
  !-------------------------------------------------
  if(use_dilaton) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Dilaton requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_dilaton=T requires poisson=T'
        call clean_stop
     end if
     if(L_dilaton <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: L_dilaton must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Dilaton gravity enabled (Brax+12 environmentally-damped)'
        write(*,'(A,F6.3,A,F8.3,A)') &
             '   beta0=', beta_dilaton, '  range L=', L_dilaton, ' Mpc/h'
        write(*,'(A,ES10.3,A,F7.4,A)') '   A2=', &
             & 1d0/(3d0*(L_dilaton/2997.92458d0)**2), '  s=3*Om=', 3d0*omega_m, &
             & '  (a0_dilaton is ignored)'
        write(*,'(A,I3,A,ES10.3)') '   max_iter=', n_iter_dilaton, '  eps=', dilaton_eps
     end if
  end if

  !-------------------------------------------------
  ! Galileon (cubic) gravity
  !-------------------------------------------------
  if(use_galileon) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: Galileon requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_galileon=T requires poisson=T'
        call clean_stop
     end if
     if(abs(c3_galileon) < 1d-30) then
        if(myid==1) write(*,*) 'ERROR: c3_galileon must be non-zero'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Cubic Galileon gravity enabled'
        write(*,'(A,ES10.3,A,ES10.3)') '   c2=', c2_galileon, '  c3=', c3_galileon
        write(*,'(A,I5,A,ES10.3)') '   max_iter=', n_iter_galileon, '  eps=', galileon_eps
     end if
  end if

  !-------------------------------------------------
  ! Coupled Dark Energy
  !-------------------------------------------------
  if(use_coupled_de) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_coupled_de=T requires poisson=T'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A,F8.4)') ' Coupled Dark Energy enabled: beta_cde=', beta_cde
        write(*,'(A,F10.6)') '   G_eff/G = ', 1d0 + 2d0*beta_cde**2
        if(use_quintessence) then
           write(*,'(A,L2,A,L2)') '   coupled quintessence: friction=', &
                & cde_friction, '  vary_mass=', cde_vary_mass
        else
           write(*,'(A)') '   (no use_quintessence: constant G_eff boost only)'
        end if
     end if
  end if

  !-------------------------------------------------
  ! Quintessence / k-essence scalar-field dark energy
  !-------------------------------------------------
  if(use_quintessence .and. use_kessence) then
     if(myid==1) write(*,*) 'ERROR: use_quintessence and use_kessence are mutually exclusive'
     call clean_stop
  end if
  if(use_quintessence .or. use_kessence) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: quintessence/k-essence requires cosmo=.true.'
        call clean_stop
     end if
     if(w0 /= -1.0d0 .or. wa /= 0.0d0) then
        if(myid==1) write(*,*) 'ERROR: scalar-field DE replaces CPL; keep w0=-1, wa=0'
        call clean_stop
     end if
     if(use_ede) then
        if(myid==1) write(*,*) 'ERROR: use_ede cannot be combined with scalar-field DE'
        call clean_stop
     end if
  end if
  if(use_quintessence) then
     if(quint_pot /= 1 .and. quint_pot /= 2) then
        if(myid==1) write(*,*) 'ERROR: quint_pot must be 1 (Ratra-Peebles) or 2 (exponential)'
        call clean_stop
     end if
     if(quint_ic_mode < 0 .or. quint_ic_mode > 1) then
        if(myid==1) write(*,*) 'ERROR: quint_ic_mode must be 0 (frozen) or 1 (RP tracker)'
        call clean_stop
     end if
     if(quint_ic_mode == 1 .and. quint_pot /= 1) then
        if(myid==1) write(*,*) 'ERROR: quint_ic_mode=1 is defined only for Ratra-Peebles phiCDM'
        call clean_stop
     end if
     if(quint_ic_mode == 1 .and. use_coupled_de) then
        if(myid==1) write(*,*) 'ERROR: RP tracker IC is not defined for coupled quintessence'
        call clean_stop
     end if
     if(quint_alpha <= 0d0 .or. quint_lambda <= 0d0 .or. &
          & (quint_ic_mode == 0 .and. quint_phi_ini <= 0d0)) then
        if(myid==1) write(*,*) 'ERROR: active quintessence potential/initial parameters must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        if(quint_pot == 1) then
           write(*,'(A,F8.4)') ' Quintessence enabled: V=A*phi^-alpha, alpha=', quint_alpha
        else
           write(*,'(A,F8.4)') ' Quintessence enabled: V=A*exp(-lambda*phi), lambda=', quint_lambda
        end if
        if(quint_ic_mode == 1) then
           write(*,'(A)') '   initial condition: matter-era Ratra-Peebles tracker (phiCDM)'
        else
           write(*,'(A,ES10.3)') '   frozen phi_ini [Mpl] =', quint_phi_ini
        end if
     end if
  end if
  if(use_kessence) then
     if(kes_x0 <= 0.5d0) then
        if(myid==1) write(*,*) 'ERROR: kes_x0 must be > 0.5'
        call clean_stop
     end if
     if(myid==1) write(*,'(A,F12.8)') &
          & ' k-essence enabled: P(X)=-X+X^2, X(a=1)/M^4=', kes_x0
  end if

  !-------------------------------------------------
  ! Horndeski quasi-static mu(a,k) gravity
  !-------------------------------------------------
  if(use_horndeski) then
     if(.not. poisson) then
        if(myid==1) write(*,*) 'ERROR: use_horndeski=T requires poisson=T'
        call clean_stop
     end if
     if(hs_mass < 0d0) then
        if(myid==1) write(*,*) 'ERROR: hs_mass must be >= 0'
        call clean_stop
     end if
     if(use_fR .or. use_nDGP .or. use_symmetron .or. use_dilaton .or. use_galileon) then
        if(myid==1) write(*,*) 'ERROR: use_horndeski cannot be combined with another MG solver'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A,F8.4,A,ES10.3,A)') ' Horndeski mu(a,k) gravity enabled: mu0=', &
             & hs_mu0, '  Compton mass=', hs_mass, ' h/Mpc'
        if(hs_mass > 0d0) write(*,'(A)') &
             & '   NOTE: k-dependence exact only in FFT Poisson paths; MG/CG use k->inf limit'
     end if
  end if

  !-------------------------------------------------
  ! Early Dark Energy (Doran-Robbers)
  !-------------------------------------------------
  if(use_ede) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: EDE requires cosmo=.true.'
        call clean_stop
     end if
     if(omega_ede < 0d0 .or. omega_ede >= 1d0) then
        if(myid==1) write(*,*) 'ERROR: omega_ede must be in [0,1)'
        call clean_stop
     end if
     if(z_ede <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: z_ede must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' Early Dark Energy (Poulin+19 fluid form) enabled'
        write(*,'(A,F8.4,A,F10.1,A,F6.3)') &
             '   omega_ede=', omega_ede, '  z_ede=', z_ede, '  w_ede=', w_ede
     end if
  end if

  !-------------------------------------------------
  ! Evaporating primordial black hole dark matter
  !-------------------------------------------------
  if(use_pbh) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'ERROR: use_pbh requires cosmo=.true.'
        call clean_stop
     end if
     if(.not. pic) then
        if(myid==1) write(*,*) 'ERROR: use_pbh requires pic=.true.'
        call clean_stop
     end if
     if(pbh_fraction < 0d0 .or. pbh_fraction > 1d0) then
        if(myid==1) write(*,*) 'ERROR: pbh_fraction must be in [0,1]'
        call clean_stop
     end if
     if(pbh_boost <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: pbh_boost must be > 0'
        call clean_stop
     end if
     if(trim(pbh_energy_sink)/='local_heat' .and. &
          & trim(pbh_energy_sink)/='uniform_heat' .and. &
          & trim(pbh_energy_sink)/='removed') then
        if(myid==1) write(*,*) &
             & 'ERROR: pbh_energy_sink must be local_heat, uniform_heat, or removed'
        call clean_stop
     end if
     if(pbh_cr_ivar/=0 .and. (pbh_cr_ivar<=ndim+2 .or. pbh_cr_ivar>nvar)) then
        if(myid==1) write(*,*) 'ERROR: pbh_cr_ivar must be 0 or a passive slot in (', &
             & ndim+2, ',', nvar, ']'
        call clean_stop
     end if
     call pbh_read_table(myid)
     call pbh_validate_models(myid)
     if(myid==1) then
        write(*,'(A)') ' Evaporating-PBH dark matter enabled'
        write(*,'(A,ES10.3,A,ES10.3,A)') '   pbh_fraction=', pbh_fraction, &
             & '  pbh_boost=', pbh_boost, '  sink='//trim(pbh_energy_sink)
     end if
  end if

  !-------------------------------------------------
  ! SGS (Sub-Grid Scale) Turbulence model
  !-------------------------------------------------
  if(use_sgs) then
     if(.not. hydro) then
        if(myid==1) write(*,*) 'ERROR: use_sgs=T requires hydro=T'
        call clean_stop
     end if
     if(sgs_C_diss <= 0d0) then
        if(myid==1) write(*,*) 'ERROR: sgs_C_diss must be > 0'
        call clean_stop
     end if
     if(myid==1) then
        write(*,'(A)') ' SGS turbulence model enabled'
        write(*,'(A,F6.3,A,F6.3,A,F6.3)') &
             '   C_prod=', sgs_C_prod, '  C_diss=', sgs_C_diss, '  C_smag=', sgs_C_smag
        write(*,'(A,ES10.3,A,F6.2,A,ES10.3,A,L1)') '   floor=', sgs_floor, &
             '  cap=', sgs_cap, '  e_init=', sgs_e_init, '  hydro=', sgs_hydro
        write(*,'(A,I3,A,I3)') '   isgs=', isgs, '  nvar=', nvar
     end if
  end if

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     allocate(tmovout(1:imovout))
     allocate(amovout(1:imovout))
     tmovout=1d100
     amovout=1d100
     if(tendmov>0)then
        do i=1,imovout
           tmovout(i)=tendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(aendmov>0)then
        do i=1,imovout
           amovout(i)=aendmov*dble(i)/dble(imovout)
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  ! [RESIZABLE] ngridmax/ngridtot set the initial capacity.  Runtime growth is
  ! enabled by default and can be disabled with ngridmax_auto=.false.
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  ! Phase 2 chunk 2: fixed block size for now; make this a namelist key later.
  amr_block_size=64
  if(mod(ngridmax,amr_block_size)/=0) &
       ngridmax=((ngridmax/amr_block_size)+1)*amr_block_size
  if(npartmax==0)then
     if(nparttot==0)then
        ! Keep zero as the initial capacity; the bundle grows on first demand.
        npartmax=0
     else
        npartmax=nparttot/int(ncpu,kind=8)
     endif
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  if(clumpfind.and.(.not.pic))then
     pic=.true.
  endif
  !if(pic.and.(.not.poisson))then
  !   poisson=.true.
  !endif

  call read_hydro_params(nml_ok)
#ifdef RT
  call rt_read_hydro_params(nml_ok)
#endif
  if(sink)then
     ! Preserve the active default when a SINK_PARAMS group omits the
     ! optional accretion_scheme key; a sentinel distinguishes that case.
     accretion_scheme_saved=accretion_scheme
     accretion_scheme='__unspecified__'
     rewind(1)
     read(1,NML=sink_params,IOSTAT=sink_nml_iostat)
     if(sink_nml_iostat>0)then
        if(myid==1)write(*,*)'ERROR: invalid &SINK_PARAMS namelist, iostat=', &
             & sink_nml_iostat
        nml_ok=.false.
        accretion_scheme=accretion_scheme_saved
     else if(sink_nml_iostat<0)then
        ! SINK_PARAMS is optional for legacy lagRamses inputs.
        accretion_scheme=accretion_scheme_saved
     else
        if(trim(accretion_scheme)=='__unspecified__')then
           accretion_scheme=accretion_scheme_saved
        else
           flux_accretion=.false.
           threshold_accretion=.false.
           bondi_accretion=.false.
           select case(trim(accretion_scheme))
           case('none')
              bondi=.false.
           case('flux')
              flux_accretion=.true.
              bondi=.false.
           case('threshold')
              threshold_accretion=.true.
              bondi=.false.
           case('bondi')
              bondi_accretion=.true.
              bondi=.true.
           case default
              if(myid==1)write(*,*)'ERROR: unsupported sink accretion_scheme=', &
                   & trim(accretion_scheme)
              nml_ok=.false.
           end select
        endif
        ! kjhan_make_sink uses Mseed rather than the upstream mass_sink_seed.
        if(mass_sink_seed>0d0)Mseed=mass_sink_seed
        if(myid==1)write(*,'(A,L1,A,A,A,L1)') &
             & ' Sink parameters: create_sinks=',create_sinks, &
             & ' accretion_scheme=',trim(accretion_scheme),' bondi=',bondi
     endif
     rewind(1)
  endif
  if (clumpfind)call read_clumpfind_params
  if (movie)call set_movie_vars


  close(1)

  !-------------------------------------------------
  ! FPR (Fixed Proper Resolution, Gnedin 2016)
  ! Must be after read_hydro_params (reads dr_proper)
  !-------------------------------------------------
  if(dr_proper > 0.0d0) then
     if(.not. cosmo) then
        if(myid==1) write(*,*) 'WARNING: dr_proper>0 but not cosmo run, disabling FPR'
        dr_proper = 0.0d0
     else
        if(myid==1) then
           write(*,'(A)') ' FPR (Fixed Proper Resolution, Gnedin 2016) enabled'
           write(*,'(A,F8.3,A)') '   dr_proper = ', dr_proper, ' kpc'
           if(q_refine_holdback) then
              write(*,'(A)') '   mode: FPR + binary holdback'
           else
              write(*,'(A)') '   mode: FPR only (no holdback)'
           end if
        end if
     end if
  end if

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif
  


  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
  
  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     m_basic_refine(i) = m_refine  (i) !(ONS)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do
  ! Initialize m_refine_eff from m_refine (FPR adjusts at runtime)
  m_refine_eff = m_refine

  ! Validate the opt-in mesh-level floor for the void target.  Refining a
  ! level creates cells at ilevel+1, hence the geometric region must be
  ! defined from levelmin through void_refine_min_level-1.
  if(void_refine)then
     if(void_refine_min_level<=levelmin)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'void_refine_min_level must exceed levelmin'
        nml_ok=.false.
     else if(void_refine_min_level>nlevelmax)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'void_refine_min_level must not exceed levelmax'
        nml_ok=.false.
     else
        do i=levelmin,void_refine_min_level-1
           if(r_refine(i)<=0.0d0)then
              if(myid==1)write(*,*)'Error in the namelist:'
              if(myid==1)write(*,*)'void refinement requires r_refine>0 at level ',i
              nml_ok=.false.
           end if
           if(a_refine(i)<=0.0d0 .or. b_refine(i)<=0.0d0 .or. exp_refine(i)<=0.0d0)then
              if(myid==1)write(*,*)'Error in the namelist:'
              if(myid==1)write(*,*)'void refinement shape parameters must be positive at level ',i
              nml_ok=.false.
           end if
        end do
        if(myid==1 .and. nml_ok)then
           write(*,'(A,I3)')' Void refinement floor enabled at level ',void_refine_min_level
        end if
     end if
  end if

  ! The V-web path is deliberately restricted to cosmological void runs.
  ! Ordinary simulations do not allocate its state or execute its analysis.
  if(void_web_refine)then
     if(void_web_env_level<0)void_web_env_level=levelmin
     if(.not.hydro)then
        if(myid==1)write(*,*)'void_web_refine requires hydro=.true.'
        nml_ok=.false.
     end if
     if(.not.cosmo)then
        if(myid==1)write(*,*)'void_web_refine is only available for cosmo=.true.'
        nml_ok=.false.
     end if
     if(simple_boundary)then
        if(myid==1)write(*,*)'void_web_refine currently requires periodic boundaries'
        nml_ok=.false.
     end if
#if NDIM != 3
     if(myid==1)write(*,*)'void_web_refine requires NDIM=3'
     nml_ok=.false.
#endif
     if(void_web_env_level/=levelmin)then
        if(myid==1)write(*,*)'void_web_env_level must equal levelmin in this implementation'
        nml_ok=.false.
     end if
     if(void_web_base_level<=void_web_env_level .or. &
          & void_web_base_level>nlevelmax)then
        if(myid==1)write(*,*)'void_web_base_level must be in (env_level, levelmax]'
        nml_ok=.false.
     end if
     if(void_web_wall_level<void_web_base_level .or. &
          & void_web_wall_level>nlevelmax)then
        if(myid==1)write(*,*)'void_web_wall_level must be in [base_level, levelmax]'
        nml_ok=.false.
     end if
     if(void_web_hydro_max_level<0) &
          & void_web_hydro_max_level=min(nlevelmax,void_web_wall_level+1)
     if(void_web_hydro_max_level<void_web_wall_level .or. &
          & void_web_hydro_max_level>nlevelmax)then
        if(myid==1)write(*,*)'void_web_hydro_max_level must be in [wall_level, levelmax]'
        nml_ok=.false.
     end if
     if(void_web_scope_ivar>nvar .or. void_web_scope_ivar<-1)then
        if(myid==1)write(*,*)'void_web_scope_ivar must be -1, 0, or a valid hydro variable'
        nml_ok=.false.
     else if(void_web_scope_ivar>0 .and. void_web_scope_ivar<=ndim+2)then
        if(myid==1)write(*,*)'void_web_scope_ivar must select a passive scalar, not a hydro variable'
        nml_ok=.false.
     end if
     if((err_grad_d>=0.0d0 .or. err_grad_p>=0.0d0 .or. &
          & err_grad_u>=0.0d0 .or. err_jump_u>=0.0d0 .or. &
          & ekin_flux_refine>=0.0d0) .and. &
          & void_web_scope_ivar<=ndim+2)then
        if(myid==1)write(*,*)'Void-web hydro triggers require a passive-scalar scope'
        nml_ok=.false.
     end if
     if(void_web_scope_cut<0.0d0 .or. void_web_scope_cut>1.0d0)then
        if(myid==1)write(*,*)'void_web_scope_cut must lie between zero and one'
        nml_ok=.false.
     end if
     if(void_web_lambda_off<0.0d0 .or. &
          & void_web_lambda_on<=void_web_lambda_off)then
        if(myid==1)write(*,*)'Require 0 <= void_web_lambda_off < void_web_lambda_on'
        nml_ok=.false.
     end if
     if(void_web_update_interval<1)then
        if(myid==1)write(*,*)'void_web_update_interval must be positive'
        nml_ok=.false.
     end if
     if(void_web_jump_pressure_min>=0.0d0 .and. err_jump_u<0.0d0)then
        if(myid==1)write(*,*)'void_web_jump_pressure_min requires err_jump_u>=0'
        nml_ok=.false.
     end if
     if(void_web_jump_compression_gate .and. err_jump_u<0.0d0)then
        if(myid==1)write(*,*)'void_web_jump_compression_gate requires err_jump_u>=0'
        nml_ok=.false.
     end if
     if(myid==1 .and. nml_ok)then
        write(*,'(A)')' Void-only V-web refinement enabled'
        write(*,'(A,I3,A,I3,A,I3)')'   environment/base/wall levels = ', &
             & void_web_env_level,' / ',void_web_base_level,' / ',void_web_wall_level
        write(*,'(A,I3)')'   scoped hydro-trigger maximum level = ', &
             & void_web_hydro_max_level
        write(*,'(A,2F8.3)')'   lambda on/off = ',void_web_lambda_on,void_web_lambda_off
        if(void_web_scope_ivar==-1)then
           write(*,'(A)')'   WARNING: global scope selected; the base floor covers the full box'
        else if(void_web_scope_ivar==0)then
           write(*,'(A)')'   scope = IC Lagrangian refinement map (cpu_map2)'
        else
           write(*,'(A,I3,A,F8.3)')'   scope passive variable = ', &
                & void_web_scope_ivar,' cut = ',void_web_scope_cut
        end if
        if(err_jump_u>=0.0d0)then
           write(*,'(A,L1,A,F8.3)')'   jump compression gate = ', &
                & void_web_jump_compression_gate,' pressure minimum = ', &
                & void_web_jump_pressure_min
        end if
     end if
  end if

#ifdef SNRT
  ! Check before any sink feedback or RT state mutation, including before
  ! sinks are first created. Do not allow rank-local environment selection
  ! to split collective RT execution or the source ownership decision.
  snrt_requested_local=merge(1,0,snrt_agn_rt_requested())
  snrt_requested_min=snrt_requested_local
  snrt_requested_max=snrt_requested_local
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(snrt_requested_local,snrt_requested_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(snrt_requested_local,snrt_requested_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
  if (snrt_requested_min/=snrt_requested_max) then
     if(myid==1)write(*,*)'SNRT_RT_ENABLE must agree across MPI ranks'
     nml_ok=.false.
  end if
  if (snrt_requested_max==1 .and. sink .and. sink_AGN .and. &
       snrt_agn_model()/=snrt_agn_model_reference) then
     if(myid==1)write(*,*)'AGN source ownership conflict: legacy feedback plus live SNRT is not approved'
     nml_ok=.false.
  end if
  if (snrt_requested_max==1 .and. sink .and. (ncpu>1 .or. nrestart>0)) then
     if(myid==1)write(*,*)'Live SNRT AGN requires serial fresh start until pending energy and photon state support restart/migration'
     nml_ok=.false.
  end if
  if (snrt_requested_max==1) then
     ! Admit the runtime contracts during namelist initialization.  The
     ! driver also keeps a defensive per-process loader, but waiting until
     ! the first AMR level would allow a missing or non-admissible contract
     ! to look like a successful no-op run.
     call snrt_spectral_contract_load_from_environment(agn_contract_error)
     if (agn_contract_error/=0 .or. .not.snrt_spectral_contract_runtime_allowed) then
        if(myid==1)write(*,'(A,I0,A,A)') &
             'SNRT startup rejected spectral contract: error=',agn_contract_error, &
             ' (',trim(snrt_spectral_contract_error_name(agn_contract_error))//')'
        if(myid==1 .and. len_trim(snrt_spectral_contract_error_message)>0) &
             write(*,'(A,A)') '  detail: ',trim(snrt_spectral_contract_error_message)
        nml_ok=.false.
     else if(myid==1)then
        write(*,'(A,A,A,A)') 'SNRT startup spectral contract admitted: status=', &
             trim(snrt_spectral_contract_status),' source=',trim(snrt_spectral_contract_source_id)
     end if
     call snrt_secondary_tables_load_from_environment(snrt_thermochemistry_error)
     if (snrt_thermochemistry_error/=snrt_thermochemistry_ok .or. &
          .not.snrt_secondary_tables_loaded) then
        if(myid==1)write(*,'(A,I0,A,A)') &
             'SNRT startup rejected thermochemistry contract: error=', &
             snrt_thermochemistry_error,' (', &
             trim(snrt_thermochemistry_error_name(snrt_thermochemistry_error))//')'
        if(myid==1 .and. len_trim(snrt_thermochemistry_error_message)>0) &
             write(*,'(A,A)') '  detail: ',trim(snrt_thermochemistry_error_message)
        nml_ok=.false.
     else if(myid==1)then
        write(*,'(A,A,A,A)') 'SNRT startup thermochemistry contract admitted: source=', &
             trim(snrt_secondary_loaded_source_id),' upstream=', &
             trim(snrt_secondary_loaded_upstream_commit)
     end if
     ! A supplied dust contract must never be silently ignored.  Candidate
     ! contracts remain available for inspection in the legacy profile.  The
     ! DUST_LIVE profile requires a version-2, physically approved contract;
     ! otherwise the executable must not enter the live driver.
     snrt_dust_contract_env = ''
     call get_environment_variable('SNRT_DUST_CONTRACT', snrt_dust_contract_env, &
          length=snrt_dust_contract_env_length)
     if (snrt_dust_contract_env_length > 0) then
        call snrt_dust_contract_load_from_environment(snrt_dust_contract_error)
        if (snrt_dust_contract_error/=0 .or. .not.snrt_dust_contract_loaded) then
           if(myid==1)write(*,'(A,I0,A,A)') &
                'SNRT startup rejected dust contract: error=', &
                snrt_dust_contract_error,' (', &
                trim(snrt_dust_contract_error_name(snrt_dust_contract_error))//')'
           if(myid==1 .and. len_trim(snrt_dust_contract_error_message)>0) &
                write(*,'(A,A)') '  detail: ',trim(snrt_dust_contract_error_message)
           nml_ok=.false.
        else if (snrt_dust_contract_runtime_allowed) then
#ifdef DUST_LIVE
           if(myid==1)write(*,'(A)') &
                'SNRT startup admitted approved DUST_LIVE contract'
#else
           if(myid==1)write(*,'(A)') &
                'SNRT startup rejected approved dust contract: live dust receiver is not wired'
           nml_ok=.false.
#endif
        else
#ifdef DUST_LIVE
           if(myid==1)write(*,'(A)') &
                'SNRT startup rejected dust contract: DUST_LIVE requires runtime-approved version-2 opacity/thermal data'
           nml_ok=.false.
#else
           if(myid==1)write(*,'(A)') &
                'SNRT dust contract loaded for inspection only; live dust remains ZERO_SCAFFOLD'
#endif
        end if
#ifdef DUST_LIVE
     else
        if(myid==1)write(*,'(A)') &
             'SNRT startup rejected: DUST_LIVE requires SNRT_DUST_CONTRACT'
        nml_ok=.false.
#endif
     end if
  end if
#endif

  ! The comparison profile is opt-in, rank-uniform, serial and fresh only.
  agn_model_local=snrt_agn_model()
  agn_model_min=agn_model_local; agn_model_max=agn_model_local
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(agn_model_local,agn_model_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(agn_model_local,agn_model_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
  if(agn_model_min<0.or.agn_model_min/=agn_model_max)then
     if(myid==1)write(*,*) 'Unknown or rank-inconsistent SNRT_AGN_MODEL'
     nml_ok=.false.
  endif
  agn_snrt_built=.false.
#ifdef SNRT
  agn_snrt_built=.true.
#endif
  if(agn_model_max==snrt_agn_model_reference)then
     if(.not.snrt_agn_reference_config_ok(agn_snrt_built,snrt_agn_rt_requested(),hydro,sink,bondi, &
          sink_AGN,mad_jet,ndim,nener,ncpu,nrestart,X_floor))then
        if(myid==1)write(*,*) 'partition_reference_v1 requires serial fresh 3D NENER=0 non-MAD Bondi + RT + AGN'
        nml_ok=.false.
     endif
#ifdef SNRT
     call snrt_spectral_contract_load_from_environment(agn_contract_error)
     if(agn_contract_error/=0.or..not.snrt_spectral_contract_runtime_allowed.or. &
          trim(snrt_spectral_contract_status)/='reference_control')then
        if(myid==1)write(*,*) 'partition_reference_v1 requires explicitly enabled reference_control SED'
        nml_ok=.false.
     endif
#endif
  endif
  call snrt_agn_admit_reference(nml_ok.and.agn_model_local==snrt_agn_model_reference)
  if(nml_ok.and.agn_model_local==snrt_agn_model_reference.and.myid==1)write(*,*) &
       'SNRT_AGN_MODEL=partition_reference_v1 comparison only; mechanical shares high=0.15 low=1; MAD excluded'

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params

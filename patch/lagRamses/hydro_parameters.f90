module hydro_parameters
  use amr_parameters

  integer,parameter::neul=ndim+2
#ifdef SOLVERmhd
  integer,parameter::nhydro=8
#else
  integer,parameter::nhydro=ndim+2
#endif

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=ndim+2+nener
#else
  integer,parameter::nvar=NVAR
#endif
  integer,parameter::nvar_all=nvar &
#ifdef SOLVERmhd
       +3 &
#endif
       +0
  logical::mhd_enabled=.false.
  ! MHD face-batch execution. CUDA HLLD is distinct from hydro-only gpu_hydro.
  logical::mhd_omp=.false.,mhd_gpu_faces=.false.
  ! Code B absorbs sqrt(4*pi): magnetic energy density is |B|^2/2.
  real(dp)::mhd_seed(3)=0d0
  ! Analytic gas-only ICs exercise the same production CT solver.
  character(len=24)::mhd_initial_condition='uniform'
  integer::slope_mag_type=1,interpol_mag_type=1
  character(len=10)::riemann2d='hlld'
#ifdef SOLVERmhd
  integer::ischeme=0,iriemann=3,iriemann2d=5
  real(dp)::eta_mag=0d0
  real(dp)::err_grad_A=-1d0,err_grad_B=-1d0,err_grad_C=-1d0,err_grad_B2=-1d0
  real(dp)::floor_A=1d-10,floor_B=1d-10,floor_C=1d-10,floor_B2=1d-10
  logical::allow_switch_solver=.true.,allow_switch_solver2D=.true.
  real(dp)::switch_solv_B=1d20,switch_solv_dens=1d20,switch_solv_min_dens=1d-20
#endif
  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar_all)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-NDIM-2-NENER)::var_bound=0.0
#endif
  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
  ! Void numerical-dissipation control (velocity-jump / KE-flux refinement).
  ! Unlike err_grad_u, these are NOT bulk-velocity-suppressed, so they flag
  ! shear/shocks and supersonic flow in fast, low-density void gas where
  ! Eulerian advection (truncation) dissipation is largest. Off by default.
  real(dp)::err_jump_u=-1.0     ! Velocity-jump (Mach-of-jump) refinement: |dv|/c_s threshold
  real(dp)::ekin_flux_refine=-1.0 ! KE-flux refinement: |v|^2/c_s^2 (kinetic/thermal) threshold
  real(dp)::d_keflux_max=-1.0   ! Density gate: jump/KE-flux act only where rho<d_keflux_max (<=0 disables gate)
  real(dp)::floor_keflux=1.d-30 ! Floor for KE-flux denominator
  logical::void_web_jump_compression_gate=.false. ! Void-only: require convergent normal flow
  real(dp)::void_web_jump_pressure_min=-1.0 ! Void-only: minimum symmetric pressure jump (<0 disables)
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:NVAR-NDIM-2)::err_grad_var=-1.0
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-NDIM-2-NENER)::var_region=0.0
#endif
  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::slope_theta=1.5d0
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
#ifdef SOLVERmhd
  character(LEN=10)::riemann='hlld'
#else
  character(LEN=10)::riemann='llf'
#endif

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6
  integer::ivirial=6
  integer::inener=6
  integer::iHydrogen=-1
  integer::iHelium=-1
  ! DUST_LIVE reserves two fields after the complete existing passive map.
  ! They are assigned by read_hydro_params only in the live profile.
  integer::idust=-1
  integer::idust_energy=-1
  integer::idust_species=-1 ! first of C, MgFeSiO4, only for explicit composition model
  integer::idust_bins=-1 ! C-small/large, silicate-small/large; two-size model only
  integer::idust_shock=-1,idust_fresh=-1 ! transient SN energy and fresh C/sil densities; cleared before RT/SF/output
  integer::ichimes=-1 ! 157 transported m_H*n_species densities; NOT additional baryon mass
  integer::idust_iron=-1 ! Fe after CHIMES or reserved dust window (static only); subset of rho/Fe/idust
  integer::idust_pah=-1 ! 128 neutral or 256 neutral/cation mass states; separate from bulk idust/energy
  ! Internal admission only: main exposes the namelist after radiation wiring.
  logical::dust_relative_motion=.false.
  integer::idust_momentum=-1,ndust_phase=0 ! three ABSOLUTE momenta per phase
  ! Explicit effective gas collision cross section for the Epstein-domain
  ! check. No implicit physical default; main owns final namelist exposure.
  real(dp)::dust_drag_collision_cross_section_cm2=0d0

contains
  pure real(dp) function magnetic_energy(row) result(e)
    real(dp),intent(in)::row(:)
    e=0d0
#ifdef SOLVERmhd
    e=.125d0*sum((row(6:8)+row(nvar+1:nvar+3))**2)
#endif
  end function magnetic_energy
end module hydro_parameters

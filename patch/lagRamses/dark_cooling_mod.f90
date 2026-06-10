!################################################################
!################################################################
! Dark Cooling Module for Atomic Dark Matter (aDM)
!
! Implements dark-sector cooling processes following
! Rosenberg & Fan (2017, PDU 15, 136) and Shen+ (2023).
!
! Dark sector: dark proton (p'), dark electron (e'), dark photon
! Dark U(1) Coulomb-like force with coupling alpha_D
!
! Five cooling processes:
!   1. Dark bremsstrahlung
!   2. Dark Compton cooling (off dark CMB)
!   3. Dark recombination cooling
!   4. Dark collisional ionization cooling
!   5. Dark collisional excitation cooling
!################################################################
module dark_cooling_mod
  use amr_parameters, only: dp
  implicit none
  private
  public :: dark_net_cooling, dark_cool_implicit

  ! Physical constants (CGS)
  real(dp),parameter::pi_dc      = 3.141592653589793d0
  real(dp),parameter::kB_cgs     = 1.380649d-16     ! Boltzmann [erg/K]
  real(dp),parameter::hbar_cgs   = 1.054571817d-27  ! hbar [erg s]
  real(dp),parameter::c_cgs      = 2.997924580d+10  ! speed of light [cm/s]
  real(dp),parameter::eV_to_erg  = 1.602176634d-12  ! 1 eV [erg]
  real(dp),parameter::GeV_to_g   = 1.78266192d-24   ! 1 GeV [g]
  real(dp),parameter::GeV_to_erg = 1.602176634d-3   ! 1 GeV [erg]
  real(dp),parameter::m_e_cgs    = 9.1093837d-28    ! electron mass [g]
  real(dp),parameter::alpha_em   = 7.2973525693d-3  ! SM fine-structure constant
  real(dp),parameter::sigma_T_sm = 6.6524587d-25    ! SM Thomson cross-section [cm^2]
  real(dp),parameter::a_rad      = 7.5657d-15       ! radiation constant [erg/cm^3/K^4]
  real(dp),parameter::T_CMB0     = 2.7255d0         ! CMB temperature today [K]

contains

!================================================================
! Dark sector derived quantities
!================================================================
subroutine dark_params(alpha_D, mp_GeV, me_ratio, xi, &
     me_D_g, E_ion_erg, T_ion_K, sigma_T_D)
  implicit none
  real(dp),intent(in)::alpha_D, mp_GeV, me_ratio, xi
  real(dp),intent(out)::me_D_g, E_ion_erg, T_ion_K, sigma_T_D

  real(dp)::me_D_GeV, mu_D_GeV

  me_D_GeV = mp_GeV * me_ratio
  me_D_g   = me_D_GeV * GeV_to_g

  ! Reduced mass for hydrogen-like atom
  mu_D_GeV = mp_GeV * me_D_GeV / (mp_GeV + me_D_GeV)

  ! Ionization energy: E_ion = 0.5 * alpha_D^2 * mu_D * c^2
  E_ion_erg = 0.5d0 * alpha_D**2 * mu_D_GeV * GeV_to_erg
  T_ion_K   = E_ion_erg / kB_cgs

  ! Dark Thomson cross-section: sigma_T' = (8pi/3)*(alpha_D/alpha_em)^2 * (m_e/m_e')^2 * sigma_T
  sigma_T_D = sigma_T_sm * (alpha_D/alpha_em)**2 * (m_e_cgs/me_D_g)**2

end subroutine dark_params

!================================================================
! Saha equilibrium ionization fraction x_e(T_D, n_D)
! n_D = total dark baryon number density [cm^{-3}]
! Returns x_e in [0,1]
!================================================================
function saha_ionization(T_D, n_D, me_D_g, E_ion_erg) result(x_e)
  implicit none
  real(dp),intent(in)::T_D, n_D, me_D_g, E_ion_erg
  real(dp)::x_e
  real(dp)::lambda_th, n_Q, ratio, disc

  ! Thermal de Broglie wavelength: lambda = hbar*sqrt(2*pi/(m_e'*kB*T))
  if(T_D < 1.0d0) then
     x_e = 0.0d0
     return
  end if
  lambda_th = hbar_cgs * sqrt(2.0d0*pi_dc / (me_D_g * kB_cgs * T_D))
  ! Quantum concentration
  n_Q = 1.0d0 / lambda_th**3

  ! Saha: x_e^2/(1-x_e) = n_Q/n_D * exp(-E_ion/(kB*T))
  ratio = (n_Q / n_D) * exp(-E_ion_erg / (kB_cgs * T_D))

  ! Solve x_e^2 + ratio*x_e - ratio = 0
  disc = ratio**2 + 4.0d0*ratio
  x_e = 0.5d0*(-ratio + sqrt(disc))
  x_e = max(0.0d0, min(1.0d0, x_e))

end function saha_ionization

!================================================================
! Dark net cooling rate Lambda_dark [erg/s/cm^3]
!
! T_D    : dark temperature [K]
! n_D    : dark baryon number density [cm^{-3}]
! aexp   : scale factor (for dark CMB temperature)
!
! Uses module-level parameters from amr_parameters:
!   adm_alpha, adm_mp, adm_me_ratio, adm_xi
!================================================================
function dark_net_cooling(T_D, n_D, aexp) result(Lambda)
  use amr_parameters, only: adm_alpha, adm_mp, adm_me_ratio, adm_xi, adm_mol
  implicit none
  real(dp),intent(in)::T_D, n_D, aexp
  real(dp)::Lambda

  real(dp)::me_D_g, E_ion_erg, T_ion_K, sigma_T_D
  real(dp)::x_e, n_e, n_i, n_n
  real(dp)::T_D_CMB, mp_D_g
  real(dp)::Lambda_brem, Lambda_compton, Lambda_recomb
  real(dp)::Lambda_ci, Lambda_cx, Lambda_mol
  real(dp)::g_ff, kT

  call dark_params(adm_alpha, adm_mp, adm_me_ratio, adm_xi, &
       me_D_g, E_ion_erg, T_ion_K, sigma_T_D)

  mp_D_g = adm_mp * GeV_to_g
  kT = kB_cgs * T_D

  ! Ionization fraction (Saha equilibrium)
  x_e = saha_ionization(T_D, n_D, me_D_g, E_ion_erg)

  n_e = x_e * n_D         ! free dark electron density
  n_i = x_e * n_D         ! free dark proton density
  n_n = (1.0d0 - x_e)*n_D ! neutral dark atom density

  ! Dark CMB temperature
  T_D_CMB = adm_xi * T_CMB0 / aexp  ! T_dark_CMB = xi * T_CMB / a

  !----- 1. Dark Bremsstrahlung -----
  ! Lambda_ff = (2^5 pi / 3) * (alpha_D^3 / m_e'^2) * sqrt(2*pi*m_e'/(3*kT))
  !             * g_ff * n_e * n_i * kT   (Rosenberg & Fan 2017, eq 3)
  ! Gaunt factor g_ff ~ 1 (for simplicity)
  g_ff = 1.0d0
  if(T_D > 1.0d0 .and. x_e > 1.0d-10) then
     Lambda_brem = (32.0d0*pi_dc/3.0d0) &
          * (adm_alpha**3 * hbar_cgs**2 / me_D_g**2) &
          * sqrt(2.0d0*pi_dc*me_D_g / (3.0d0*kT)) &
          * g_ff * n_e * n_i * kT / hbar_cgs
     ! Convert: factor has units. Use exact formula:
     ! Λ_ff = 1.426e-27 * (α'/α)^3 * (m_e/m_e')^{1/2} * T^{1/2} * g_ff * n_e * n_i * Z^2
     ! Scaling from SM bremsstrahlung (Rybicki & Lightman):
     Lambda_brem = 1.426d-27 * (adm_alpha/alpha_em)**3 &
          * sqrt(m_e_cgs/me_D_g) * sqrt(T_D) * g_ff * n_e * n_i
  else
     Lambda_brem = 0.0d0
  end if

  !----- 2. Dark Compton Cooling -----
  ! Lambda_C = (4*sigma_T'*a_rad'*T_DCMB^4)/(m_e'*c) * n_e * kB*(T_D - T_DCMB)
  ! where a_rad' = a_rad (same Stefan-Boltzmann for dark photon)
  if(x_e > 1.0d-10 .and. T_D > T_D_CMB) then
     Lambda_compton = 4.0d0 * sigma_T_D * a_rad * T_D_CMB**4 &
          / (me_D_g * c_cgs) * n_e * kB_cgs * (T_D - T_D_CMB)
  else
     Lambda_compton = 0.0d0
  end if

  !----- 3. Dark Recombination Cooling -----
  ! Lambda_rec ~ alpha_rec * n_e * n_i * E_ion
  ! alpha_rec ~ alpha_D^5 * (m_e'/m_e)^{3/2} * alpha_B(T)
  ! Use case-B-like scaling: alpha_B ~ 2.06e-11 * (alpha'/alpha)^3 * (m_e'/m_e)^{3/2}
  !                          * T_ion^{1/2} * (T/T_ion)^{-0.7}
  if(T_D > 1.0d0 .and. x_e > 1.0d-10) then
     Lambda_recomb = 2.06d-11 * (adm_alpha/alpha_em)**3 &
          * (me_D_g/m_e_cgs)**1.5d0 &
          * sqrt(T_ion_K) * (T_D/T_ion_K)**(-0.7d0) &
          * n_e * n_i * kB_cgs * T_D
     ! Ensure reasonable: recombination cools by releasing E_ion but
     ! net cooling = recomb_rate * kB*T (kinetic energy removed)
  else
     Lambda_recomb = 0.0d0
  end if

  !----- 4. Dark Collisional Ionization -----
  ! Lambda_ci ~ n_e * n_n * beta_ci * E_ion
  ! beta_ci ~ sigma_0 * sqrt(kT/m_e') * exp(-E_ion/kT)
  ! Using SM scaling: beta_ci ~ 5.85e-11 * (alpha'/alpha)^2 * sqrt(T)
  !                   * (E_ion/kT)^{-1} * exp(-E_ion/kT) * (1 + (E_ion/kT)^{-1/2})^{-1}
  if(T_D > 1.0d0 .and. x_e > 1.0d-10 .and. (1.0d0-x_e) > 1.0d-10) then
     Lambda_ci = 5.85d-11 * (adm_alpha/alpha_em)**2 * sqrt(T_D) &
          * exp(-E_ion_erg/kT) / (1.0d0 + sqrt(E_ion_erg/kT)) &
          * n_e * n_n * E_ion_erg
  else
     Lambda_ci = 0.0d0
  end if

  !----- 5. Dark Collisional Excitation -----
  ! Dominant transition: 1s -> 2p, E_21 = (3/4)*E_ion
  ! Lambda_cx ~ n_e * n_n * q_12 * E_21
  ! q_12 ~ (alpha_D/alpha_em)^2 * 2.41e-6 / sqrt(T) * (E_21/kT)^{-1/3} * exp(-E_21/kT)
  if(T_D > 1.0d0 .and. x_e > 1.0d-10 .and. (1.0d0-x_e) > 1.0d-10) then
     Lambda_cx = 2.41d-6 * (adm_alpha/alpha_em)**2 / sqrt(T_D) &
          * (0.75d0*E_ion_erg/kT)**(-1.0d0/3.0d0) &
          * exp(-0.75d0*E_ion_erg/kT) &
          * n_e * n_n * 0.75d0*E_ion_erg
  else
     Lambda_cx = 0.0d0
  end if

  !----- 6. Dark Molecular (H2') Line Cooling -----
  ! Re-scaled SM H2 rovibrational cooling (Ryan, Gurian, Shandera & Jeong 2021).
  ! Dominant below the atomic floor; needs neutral dark hydrogen.
  Lambda_mol = 0.0d0
  if(adm_mol .and. T_D > 1.0d0) then
     Lambda_mol = dark_mol_cooling(T_D, n_D, x_e)
  end if

  Lambda = Lambda_brem + Lambda_compton + Lambda_recomb &
       + Lambda_ci + Lambda_cx + Lambda_mol

end function dark_net_cooling

!================================================================
! Dark molecular (H2') line cooling rate [erg/s/cm^3]
!
! Re-scaled Standard-Model H2 rovibrational cooling following
! Ryan, Gurian, Shandera & Jeong (2021, "Molecular Chemistry for
! Dark Matter"). SM base rates: Galli & Palla (1998) low-density
! limit + Hollenbach & McKee (1979)/Glover & Abel (2008) LTE fits.
!
! Re-scaling (their Table 1, Eqs 3.15, 3.21, 3.23, 3.31-3.34):
!   r_alpha = alpha_D/alpha_SM, r_m = m_e'/m_e, r_M = m_p'/m_p
!   low-density (rotational):  amp = r_a r_m r_M^-2 ; T~_r = (r_M/(r_a^2 r_m^2)) T
!   LTE rotational:            amp = r_a^9 r_m^8 r_M^-6 ; T~_r
!   LTE vibrational:           amp = r_a^9 r_m^5 r_M^-3 ; T~_v = (r_M^1/2/(r_m^3/2 r_a^2)) T
!   density interpolation:     C = C_LTE / (1 + C_LTE/(C_n0 * n_X))
!
! PHASE-1 (equilibrium) approximations:
!   * dark-H2 fraction fixed by namelist adm_fH2 (n_H2'/n_neutral)
!   * single rotational re-scaling applied to the GP98 total-rovib
!     low-density fit (T0 -> infinity fallback, paper Sec 3.6)
! Phase-2 will replace adm_fH2 with a tracked non-equilibrium x_H2'.
!================================================================
function dark_mol_cooling(T_D, n_D, x_e) result(Lambda_mol)
  use amr_parameters, only: adm_alpha, adm_mp, adm_me_ratio, adm_fH2
  implicit none
  real(dp),intent(in)::T_D, n_D, x_e
  real(dp)::Lambda_mol

  ! SM reference micro-physics
  real(dp),parameter::mp_SM    = 0.9382720813d0    ! proton mass [GeV]
  real(dp),parameter::me_SM    = 0.5109989461d-3   ! electron mass [GeV]
  real(dp),parameter::alpha_SM = 7.2973525693d-3   ! SM fine-structure constant

  real(dp)::r_alf, r_me, r_mp
  real(dp)::n_neutral, n_H2, n_HI
  real(dp)::Tr, Tv, T3r, T3v
  real(dp)::amp_LD, amp_rotLTE, amp_vibLTE
  real(dp)::C_n0, C_rotLTE, C_vibLTE, C_LTE, C_per_mol

  Lambda_mol = 0.0d0
  if(adm_fH2 <= 0.0d0) return

  ! Neutral dark hydrogen: molecules + atomic collider
  n_neutral = (1.0d0 - x_e) * n_D
  if(n_neutral <= 0.0d0) return
  n_H2 = 0.5d0 * adm_fH2 * n_neutral          ! 2 atoms per molecule
  n_HI = (1.0d0 - adm_fH2) * n_neutral         ! atomic collider density
  if(n_H2 <= 0.0d0 .or. n_HI <= 0.0d0) return

  ! Dark/SM parameter ratios
  r_alf = adm_alpha / alpha_SM
  r_me  = (adm_mp * adm_me_ratio) / me_SM       ! dark electron / SM electron
  r_mp  = adm_mp / mp_SM                         ! dark proton / SM proton

  ! Re-scaled SM evaluation temperatures
  Tr = (r_mp / (r_alf**2 * r_me**2)) * T_D                ! rotational
  Tv = (sqrt(r_mp) / (r_me**1.5d0 * r_alf**2)) * T_D      ! vibrational

  ! Amplitude factors
  amp_LD     = r_alf * r_me / r_mp**2           ! Eq 3.21 (rotational)
  amp_rotLTE = r_alf**9 * r_me**8 / r_mp**6     ! Eq 3.23
  amp_vibLTE = r_alf**9 * r_me**5 / r_mp**3     ! A_vib*E_vib re-scaling

  ! --- Low-density limit: GP98 total-rovib fit [erg cm^3/s, per H2 per H] ---
  C_n0 = amp_LD * sm_h2_lowden(Tr)

  ! --- LTE limit: HM79/GP98 rot + vib fits [erg/s, per H2] ---
  T3r = Tr * 1.0d-3
  T3v = Tv * 1.0d-3
  C_rotLTE = amp_rotLTE * ( 9.5d-22 * T3r**3.76d0 &
                 / (1.0d0 + 0.12d0 * T3r**2.1d0) &
                 * exp(-(0.13d0/max(T3r,1.0d-30))**3) &
                 + 3.0d-24 * exp(-0.51d0/max(T3r,1.0d-30)) )
  C_vibLTE = amp_vibLTE * ( 6.7d-19 * exp(-5.86d0/max(T3v,1.0d-30)) &
                 + 1.6d-18 * exp(-11.7d0/max(T3v,1.0d-30)) )
  C_LTE = C_rotLTE + C_vibLTE

  ! Density interpolation (Eq 3.15): C_n0*n_HI is the low-density rate/molecule
  if(C_n0 <= 0.0d0 .or. C_LTE <= 0.0d0) return
  C_per_mol = C_LTE / (1.0d0 + C_LTE / (C_n0 * n_HI))

  Lambda_mol = n_H2 * C_per_mol

end function dark_mol_cooling

!================================================================
! Standard-Model H2 cooling, low-density limit (Galli & Palla 1998)
! Returns cooling-rate coefficient [erg cm^3 s^-1] (per H2 per H).
! Fit valid 10 K <= T <= 1e4 K; clamp outside that range.
!================================================================
function sm_h2_lowden(T_in) result(L)
  implicit none
  real(dp),intent(in)::T_in
  real(dp)::L
  real(dp)::Tc, x, logL

  Tc = min(max(T_in, 10.0d0), 1.0d4)
  x  = log10(Tc)
  logL = -103.0d0 + 97.59d0*x - 48.05d0*x**2 &
       + 10.80d0*x**3 - 0.9032d0*x**4
  L = 10.0d0**logL

end function sm_h2_lowden

!================================================================
! Backward Euler implicit update for dark internal energy
!
! edp_old : dark internal energy per unit mass [code units, erg/g]
! rho_D   : dark matter density [g/cm^3]
! n_D     : dark baryon number density [cm^{-3}]
! dt_phys : physical timestep [s]
! aexp    : scale factor
!
! Returns edp_new [same units as edp_old]
!================================================================
function dark_cool_implicit(edp_old, rho_D, n_D, dt_phys, aexp) result(edp_new)
  use amr_parameters, only: adm_mp
  implicit none
  real(dp),intent(in)::edp_old, rho_D, n_D, dt_phys, aexp
  real(dp)::edp_new

  real(dp)::T_D, Lambda, dE
  real(dp)::T_floor, edp_floor
  integer::iter
  integer,parameter::max_iter = 10
  real(dp),parameter::tol = 1.0d-4

  ! Temperature floor: 1 K
  T_floor = 1.0d0
  ! edp = (3/2) * n_D * kB * T_D / rho_D  (energy per unit mass)
  edp_floor = 1.5d0 * n_D * kB_cgs * T_floor / rho_D

  edp_new = edp_old
  if(edp_new < edp_floor) edp_new = edp_floor

  ! Simple backward Euler with subcycling if needed
  do iter = 1, max_iter
     ! Current dark temperature from internal energy
     ! edp = (3/2) * n_D * kB * T / rho   =>  T = (2/3) * edp * rho / (n_D * kB)
     T_D = (2.0d0/3.0d0) * edp_new * rho_D / (n_D * kB_cgs)
     if(T_D < T_floor) exit

     ! Cooling rate [erg/s/cm^3]
     Lambda = dark_net_cooling(T_D, n_D, aexp)
     if(Lambda <= 0.0d0) exit  ! no cooling (could be heating from Compton if T < T_DCMB)

     ! Energy loss per unit mass: dE/dt = -Lambda / rho
     dE = Lambda * dt_phys / rho_D

     edp_new = edp_new - dE
     if(edp_new < edp_floor) then
        edp_new = edp_floor
        exit
     end if

     ! Check convergence
     if(abs(dE) < tol * edp_new) exit
  end do

  if(edp_new < edp_floor) edp_new = edp_floor

end function dark_cool_implicit

end module dark_cooling_mod

module scalar_de_commons
  !--------------------------------------------------------------
  ! Scalar-field dark energy backgrounds beyond CPL:
  !
  !  - Quintessence (use_quintessence): Klein-Gordon field phi
  !    integrated from a=1e-6, potential selected by quint_pot:
  !      1 = Ratra-Peebles  V = A * phi^(-quint_alpha)
  !      2 = Exponential    V = A * exp(-quint_lambda*phi)
  !    Amplitude A fixed by shooting so that Omega_phi(a=1)=omega_l.
  !    Rest-frame sound speed cs2 = 1.
  !
  !  - Purely kinetic k-essence (use_kessence): P(X) = -X + X^2
  !    (X in units of M^4). Shift-symmetry integral
  !      P_X * sqrt(2X) * a^3 = const
  !    gives X(a) algebraically from X(a=1)=kes_x0 (>1/2).
  !    w -> 1/3 early (dark-radiation-like), w -> -1 late.
  !    cs2(a) = P_X/(P_X + 2X P_XX) = (2X-1)/(6X-1).
  !
  !  - Coupled quintessence (use_coupled_de + use_quintessence):
  !    DM mass m(phi) = m0*exp(-beta_cde*phi/Mpl) (Amendola 2000), i.e.
  !      rho_dm = rho_dm0 * a^-3 * exp(-beta_cde*(phi-phi0))
  !    KG source +beta_cde*rho_dm (pushes phi up the RP potential),
  !    fifth force G_eff=G(1+2 beta^2) (applied elsewhere), particle
  !    velocity term +beta*phidot*v (momentum conservation as m drops).
  !
  !  - Horndeski quasi-static parametrized gravity (use_horndeski):
  !      mu(a)   = 1 + hs_mu0 * Omega_de(a)/Omega_de(a=1)
  !      mu(a,k) = 1 + (mu(a)-1) * k^2/(k^2 + (a*hs_mass)^2)
  !    (hs_mass in h/Mpc; hs_mass=0 -> scale-independent)
  !
  ! All backgrounds are tabulated on a uniform ln(a) grid from
  ! a=1e-6 to a~1.55 and accessed via interpolators:
  !   sde_fde_of_a(a)    = rho_de(a)/rho_de(1)   (CPL f_de analogue)
  !   sde_w_of_a(a)      = w_de(a)
  !   sde_cs2_of_a(a)    = rest-frame cs^2(a)
  !   sde_dmcorr_of_a(a) = rho_dm(a)*a^3/rho_dm0 (coupled DE, else 1)
  !   sde_phip_of_a(a)   = dphi/dlna [Mpl]       (coupled DE friction)
  !
  ! Initialization is lazy and self-healing: omega_m/omega_l may be
  ! overridden by IC headers (grafic/gadget) after read_params, so
  ! every accessor re-tabulates if the cosmology it was built with
  ! has changed. First call happens in serial init code.
  !--------------------------------------------------------------
  use amr_parameters, only: dp
  implicit none

  integer, parameter :: sde_i0 = 3969        ! index of a=1 (N=0)
  integer, parameter :: sde_na = 4096        ! total table size
  real(dp), parameter :: sde_lna_min = -13.815510557964274d0  ! ln(1e-6)

  logical  :: sde_loaded = .false.
  real(dp) :: sde_dlna
  real(dp) :: sde_tab_fde(sde_na)     ! rho_de(a)/rho_de(1)
  real(dp) :: sde_tab_w(sde_na)       ! w_de(a)
  real(dp) :: sde_tab_cs2(sde_na)     ! rest-frame cs^2
  real(dp) :: sde_tab_dmcorr(sde_na)  ! rho_dm(a)*a^3/rho_dm0
  real(dp) :: sde_tab_phip(sde_na)    ! dphi/dlna [Mpl]

  ! Cosmology/model snapshot used to build the tables (self-healing)
  real(dp) :: sde_om_used = -1d0, sde_ol_used = -1d0, sde_ok_used = -1d0
  real(dp) :: sde_par1_used = -1d99, sde_par2_used = -1d99
  real(dp) :: sde_par3_used = -1d99, sde_beta_used = -1d99
  integer  :: sde_pot_used = -1

contains

  !--------------------------------------------------------------
  logical function sde_active()
    use amr_parameters, only: use_quintessence, use_kessence
    sde_active = use_quintessence .or. use_kessence
  end function sde_active

  !--------------------------------------------------------------
  ! Lazy init; rebuild if cosmology or model parameters changed
  ! (IC headers may override omega_m/omega_l after read_params).
  !--------------------------------------------------------------
  subroutine sde_check_init()
    use amr_parameters, only: omega_m, omega_l, omega_k, &
         use_quintessence, use_coupled_de, quint_pot, quint_alpha, &
         quint_lambda, quint_phi_ini, kes_x0, beta_cde
    real(dp) :: p1, p2, p3, bet
    integer  :: ip
    if(use_quintessence) then
       ip = quint_pot; p1 = quint_alpha; p2 = quint_lambda; p3 = quint_phi_ini
    else
       ip = 0; p1 = kes_x0; p2 = 0d0; p3 = 0d0
    end if
    if(use_coupled_de .and. use_quintessence) then
       bet = beta_cde
    else
       bet = 0d0
    end if
    if(sde_loaded) then
       if(omega_m == sde_om_used .and. omega_l == sde_ol_used .and. &
            omega_k == sde_ok_used .and. ip == sde_pot_used .and. &
            p1 == sde_par1_used .and. p2 == sde_par2_used .and. &
            p3 == sde_par3_used .and. bet == sde_beta_used) return
    end if
    call sde_init()
    sde_om_used = omega_m; sde_ol_used = omega_l; sde_ok_used = omega_k
    sde_pot_used = ip; sde_par1_used = p1; sde_par2_used = p2
    sde_par3_used = p3; sde_beta_used = bet
  end subroutine sde_check_init

  !--------------------------------------------------------------
  subroutine sde_init()
    use amr_parameters, only: use_quintessence, use_kessence
    use amr_commons, only: myid
    if(use_quintessence) then
       call quint_tabulate()
    else if(use_kessence) then
       call kessence_tabulate()
    else
       return
    end if
    sde_loaded = .true.
    if(myid==1) call sde_report()
  end subroutine sde_init

  !--------------------------------------------------------------
  ! Quintessence background: RK4 in N=ln(a), state (u,v)=(phi,dphi/dN)
  ! in Mpl units. Densities in units of rho_crit0 = 3 H0^2 Mpl^2:
  !   rho_m = Cm * exp(-3N) * exp(-beta*u)       (coupled matter)
  !   rho_k = omega_k * exp(-2N)
  !   E^2   = (rho_m + rho_k + V) / (1 - v^2/6)
  !   dE2/dN = -3 rho_m - 2 rho_k - E^2 v^2
  !   v'    = -(3 + dE2/(2E2)) v - 3 (V_u - beta*rho_m)/E^2
  ! Amplitude ln(A) fixed by bisection so rho_phi(N=0) = omega_l.
  ! For coupled DE, outer fixed-point on Cm = omega_m*exp(+beta*phi0).
  !--------------------------------------------------------------
  subroutine quint_tabulate()
    use amr_parameters, only: omega_m, omega_l, use_coupled_de, &
         use_quintessence, beta_cde
    use amr_commons, only: myid
    real(dp) :: lnA_lo, lnA_hi, lnA, rphi_end, u0, cm, bet, u0_prev
    real(dp) :: f_lo, f_hi, fmid
    integer  :: it, ic
    logical  :: ok

    bet = 0d0
    if(use_coupled_de .and. use_quintessence) bet = beta_cde
    u0_prev = 0d0
    cm = omega_m

    do ic = 1, 10   ! outer fixed-point for the coupling (1 pass if beta=0)
       cm = omega_m * exp(bet*u0_prev)
       lnA_lo = -60d0
       lnA_hi =  60d0
       call quint_integrate(lnA_lo, cm, bet, .false., rphi_end, u0, ok)
       f_lo = rphi_end - omega_l
       call quint_integrate(lnA_hi, cm, bet, .false., rphi_end, u0, ok)
       f_hi = rphi_end - omega_l
       if(f_lo*f_hi > 0d0) then
          if(myid==1) write(*,*) 'ERROR: quintessence shooting not bracketed', f_lo, f_hi
          call clean_stop
       end if
       do it = 1, 200
          lnA = 0.5d0*(lnA_lo + lnA_hi)
          call quint_integrate(lnA, cm, bet, .false., rphi_end, u0, ok)
          fmid = rphi_end - omega_l
          if(fmid > 0d0) then
             lnA_hi = lnA
          else
             lnA_lo = lnA
          end if
          if(abs(fmid) < 1d-13 .or. lnA_hi-lnA_lo < 1d-14) exit
       end do
       if(bet == 0d0) exit
       if(abs(u0 - u0_prev) < 1d-10) then
          u0_prev = u0
          exit
       end if
       u0_prev = u0
    end do

    ! Final pass: fill the module tables
    call quint_integrate(0.5d0*(lnA_lo+lnA_hi), cm, bet, .true., rphi_end, u0, ok)
    if(.not. ok) then
       if(myid==1) write(*,*) 'ERROR: quintessence background integration failed'
       call clean_stop
    end if
  end subroutine quint_tabulate

  !--------------------------------------------------------------
  subroutine quint_integrate(lnA, cm, bet, fill, rphi_end, u_end, ok)
    use amr_parameters, only: quint_pot, quint_alpha, quint_lambda, &
         quint_phi_ini
    real(dp), intent(in)  :: lnA, cm, bet
    logical,  intent(in)  :: fill
    real(dp), intent(out) :: rphi_end, u_end
    logical,  intent(out) :: ok
    real(dp) :: u, v, xn, h, e2, rphi, rphi0, u0
    real(dp) :: k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v
    integer  :: i

    sde_dlna = -sde_lna_min/dble(sde_i0-1)
    h = sde_dlna
    u = max(quint_phi_ini, 1d-8)
    v = 0d0
    ok = .true.
    rphi0 = 1d0

    do i = 1, sde_na
       xn = sde_lna_min + dble(i-1)*sde_dlna
       call quint_state(xn, u, v, lnA, cm, bet, e2, rphi, ok)
       if(.not. ok) then
          rphi_end = huge(1d0)*1d-10  ! poison: drives bisection downward
          u_end = u
          return
       end if
       if(fill) then
          sde_tab_w(i)      = (rphi - 2d0*quint_pot_v(u,lnA)) / rphi
          sde_tab_cs2(i)    = 1d0
          sde_tab_dmcorr(i) = exp(-bet*u)  ! /exp(-bet*u0) applied below
          sde_tab_phip(i)   = v
          sde_tab_fde(i)    = rphi         ! /rphi(a=1) applied below
       end if
       if(i == sde_i0) then
          rphi_end = rphi
          u_end = u
          u0 = u
       end if
       if(i == sde_na) exit
       ! RK4 step from xn to xn+h
       call quint_rhs(xn,        u,           v,           lnA, cm, bet, k1u, k1v, ok); if(.not.ok) exit
       call quint_rhs(xn+0.5d0*h, u+0.5d0*h*k1u, v+0.5d0*h*k1v, lnA, cm, bet, k2u, k2v, ok); if(.not.ok) exit
       call quint_rhs(xn+0.5d0*h, u+0.5d0*h*k2u, v+0.5d0*h*k2v, lnA, cm, bet, k3u, k3v, ok); if(.not.ok) exit
       call quint_rhs(xn+h,      u+h*k3u,     v+h*k3v,     lnA, cm, bet, k4u, k4v, ok); if(.not.ok) exit
       u = u + h*(k1u + 2d0*k2u + 2d0*k3u + k4u)/6d0
       v = v + h*(k1v + 2d0*k2v + 2d0*k3v + k4v)/6d0
    end do

    if(.not. ok) then
       rphi_end = huge(1d0)*1d-10
       u_end = u
       return
    end if

    if(fill) then
       rphi0 = sde_tab_fde(sde_i0)
       do i = 1, sde_na
          sde_tab_fde(i)    = sde_tab_fde(i) / rphi0
          sde_tab_dmcorr(i) = sde_tab_dmcorr(i) / exp(-bet*u0)
       end do
    end if
  end subroutine quint_integrate

  !--------------------------------------------------------------
  ! Potential V(phi)/rho_crit0 and its derivative, overflow-safe
  !--------------------------------------------------------------
  function quint_pot_v(u, lnA) result(vh)
    use amr_parameters, only: quint_pot, quint_alpha, quint_lambda
    real(dp), intent(in) :: u, lnA
    real(dp) :: vh, lnv
    if(quint_pot == 2) then
       lnv = lnA - quint_lambda*u
    else
       lnv = lnA - quint_alpha*log(max(u,1d-30))
    end if
    vh = exp(min(lnv, 200d0))
  end function quint_pot_v

  function quint_pot_dv(u, lnA) result(vp)
    use amr_parameters, only: quint_pot, quint_alpha, quint_lambda
    real(dp), intent(in) :: u, lnA
    real(dp) :: vp
    if(quint_pot == 2) then
       vp = -quint_lambda * quint_pot_v(u, lnA)
    else
       vp = -quint_alpha/max(u,1d-30) * quint_pot_v(u, lnA)
    end if
  end function quint_pot_dv

  !--------------------------------------------------------------
  subroutine quint_state(xn, u, v, lnA, cm, bet, e2, rphi, ok)
    use amr_parameters, only: omega_k
    real(dp), intent(in)  :: xn, u, v, lnA, cm, bet
    real(dp), intent(out) :: e2, rphi
    logical,  intent(out) :: ok
    real(dp) :: rhom, rhok, vh, denom
    ok = .true.
    rhom = cm * exp(-3d0*xn - bet*u)
    rhok = omega_k * exp(-2d0*xn)
    vh   = quint_pot_v(u, lnA)
    denom = 1d0 - v*v/6d0
    if(denom < 1d-8) then
       ok = .false.; e2 = 1d0; rphi = 0d0
       return
    end if
    e2 = (rhom + rhok + vh)/denom
    if(e2 <= 0d0) then
       ok = .false.; e2 = 1d0; rphi = 0d0
       return
    end if
    rphi = e2*v*v/6d0 + vh
  end subroutine quint_state

  !--------------------------------------------------------------
  subroutine quint_rhs(xn, u, v, lnA, cm, bet, du, dv, ok)
    use amr_parameters, only: omega_k
    real(dp), intent(in)  :: xn, u, v, lnA, cm, bet
    real(dp), intent(out) :: du, dv
    logical,  intent(out) :: ok
    real(dp) :: e2, rphi, rhom, rhok, de2
    call quint_state(xn, u, v, lnA, cm, bet, e2, rphi, ok)
    if(.not. ok) then
       du = 0d0; dv = 0d0
       return
    end if
    rhom = cm * exp(-3d0*xn - bet*u)
    rhok = omega_k * exp(-2d0*xn)
    de2  = -3d0*rhom - 2d0*rhok - e2*v*v
    du = v
    dv = -(3d0 + de2/(2d0*e2))*v - 3d0*(quint_pot_dv(u,lnA) - bet*rhom)/e2
  end subroutine quint_rhs

  !--------------------------------------------------------------
  ! Purely kinetic k-essence P(X) = -X + X^2 (X in M^4 units):
  !   (2X-1)*sqrt(X) = C0 * a^-3,  C0 = (2*x0-1)*sqrt(x0)
  !   rho = 3X^2 - X,  w = (X^2-X)/(3X^2-X),  cs2 = (2X-1)/(6X-1)
  ! Amplitude cancels in f_de = rho(X(a))/rho(x0). No shooting.
  !--------------------------------------------------------------
  subroutine kessence_tabulate()
    use amr_parameters, only: kes_x0
    use amr_commons, only: myid
    real(dp) :: c0, x0, x, targ, g, gp, dxs, rho0, xn
    integer  :: i, it

    x0 = kes_x0
    if(x0 <= 0.5d0) then
       if(myid==1) write(*,*) 'ERROR: kes_x0 must be > 0.5, got', x0
       call clean_stop
    end if
    c0 = (2d0*x0 - 1d0)*sqrt(x0)
    rho0 = 3d0*x0*x0 - x0
    sde_dlna = -sde_lna_min/dble(sde_i0-1)

    x = x0
    ! March from late times (largest i) to early: X grows monotonically
    do i = sde_na, 1, -1
       xn = sde_lna_min + dble(i-1)*sde_dlna
       targ = c0 * exp(-3d0*xn)
       ! Newton on g(X) = (2X-1)*sqrt(X) - targ, monotone for X>1/2
       x = max(x, 0.5d0 + 0.25d0*(targ/2d0)**(2d0/3d0))
       if(targ > 10d0) x = max(x, (targ/2d0)**(2d0/3d0))
       do it = 1, 100
          g  = (2d0*x - 1d0)*sqrt(x) - targ
          gp = 2d0*sqrt(x) + (2d0*x - 1d0)/(2d0*sqrt(x))
          dxs = g/gp
          if(x - dxs <= 0.5d0) dxs = 0.5d0*(x - 0.5d0)  ! damped step
          x = x - dxs
          if(abs(dxs) < 1d-14*max(x,1d0)) exit
       end do
       sde_tab_fde(i)    = (3d0*x*x - x)/rho0
       sde_tab_w(i)      = (x*x - x)/(3d0*x*x - x)
       sde_tab_cs2(i)    = (2d0*x - 1d0)/(6d0*x - 1d0)
       sde_tab_dmcorr(i) = 1d0
       sde_tab_phip(i)   = 0d0
    end do
  end subroutine kessence_tabulate

  !--------------------------------------------------------------
  ! Table interpolators (linear in ln a, clamped at both ends)
  !--------------------------------------------------------------
  function sde_interp(tab, a) result(val)
    real(dp), intent(in) :: tab(sde_na), a
    real(dp) :: val, xn, t
    integer  :: i
    xn = log(max(a, 1d-30))
    t  = (xn - sde_lna_min)/sde_dlna + 1d0
    if(t <= 1d0) then
       val = tab(1)
    else if(t >= dble(sde_na)) then
       val = tab(sde_na)
    else
       i = int(t)
       val = tab(i) + (t - dble(i))*(tab(i+1) - tab(i))
    end if
  end function sde_interp

  function sde_fde_of_a(a) result(val)
    real(dp), intent(in) :: a
    real(dp) :: val
    call sde_check_init()
    val = sde_interp(sde_tab_fde, a)
  end function sde_fde_of_a

  function sde_w_of_a(a) result(val)
    real(dp), intent(in) :: a
    real(dp) :: val
    call sde_check_init()
    val = sde_interp(sde_tab_w, a)
  end function sde_w_of_a

  function sde_cs2_of_a(a) result(val)
    real(dp), intent(in) :: a
    real(dp) :: val
    call sde_check_init()
    val = sde_interp(sde_tab_cs2, a)
  end function sde_cs2_of_a

  function sde_dmcorr_of_a(a) result(val)
    real(dp), intent(in) :: a
    real(dp) :: val
    call sde_check_init()
    val = sde_interp(sde_tab_dmcorr, a)
  end function sde_dmcorr_of_a

  function sde_phip_of_a(a) result(val)
    real(dp), intent(in) :: a
    real(dp) :: val
    call sde_check_init()
    val = sde_interp(sde_tab_phip, a)
  end function sde_phip_of_a

  !--------------------------------------------------------------
  ! Horndeski quasi-static mu(a): 1 + hs_mu0*Omega_de(a)/Omega_de(1)
  ! Works with LCDM, CPL, or scalar-field DE backgrounds.
  ! The k-dependent form mu(a,k)=1+(mu(a)-1)*k^2/(k^2+(a*hs_mass)^2)
  ! is applied in the FFT Poisson paths; MG/CG paths use this
  ! scale-independent limit (exact when hs_mass=0).
  !--------------------------------------------------------------
  function horndeski_mu_of_a(a) result(mu)
    use amr_parameters, only: use_horndeski, hs_mu0, omega_m, omega_l, &
         omega_k, w0, wa
    real(dp), intent(in) :: a
    real(dp) :: mu, fde, e2
    mu = 1d0
    if(.not. use_horndeski) return
    if(sde_active()) then
       fde = sde_fde_of_a(a)
    else if(wa == 0d0 .and. w0 == -1d0) then
       fde = 1d0
    else if(wa == 0d0) then
       fde = a**(-3d0*(1d0+w0))
    else
       fde = a**(-3d0*(1d0+w0+wa)) * exp(-3d0*wa*(1d0-a))
    end if
    e2 = omega_m/a**3 + omega_k/a**2 + omega_l*fde
    mu = 1d0 + hs_mu0 * fde/e2
  end function horndeski_mu_of_a

  !--------------------------------------------------------------
  subroutine sde_report()
    use amr_parameters, only: use_quintessence, use_kessence, &
         use_coupled_de, beta_cde, omega_m, omega_l
    use amr_commons, only: aexp_ini
    real(dp) :: w0_eff, wz1, fde_i, ode_i, ai
    ai = max(aexp_ini, 1d-6)
    w0_eff = sde_interp(sde_tab_w, 1d0)
    wz1    = sde_interp(sde_tab_w, 0.5d0)
    fde_i  = sde_interp(sde_tab_fde, ai)
    ode_i  = omega_l*fde_i*ai**3 / (omega_m + omega_l*fde_i*ai**3)
    if(use_quintessence) then
       write(*,'(" Quintessence background tabulated:")')
    else if(use_kessence) then
       write(*,'(" k-essence background tabulated:")')
    end if
    write(*,'("   w(z=0) =",F9.5,"   w(z=1) =",F9.5)') w0_eff, wz1
    write(*,'("   rho_de(a_ini)/rho_m(a_ini) =",1pe11.3," at a_ini =",1pe11.3)') &
         & ode_i/max(1d0-ode_i,1d-30), ai
    if(use_coupled_de .and. use_quintessence) then
       write(*,'("   coupled DE: beta =",F8.4,"  dm mass corr(a_ini) =",F10.6)') &
            & beta_cde, sde_interp(sde_tab_dmcorr, ai)
    end if
  end subroutine sde_report

end module scalar_de_commons

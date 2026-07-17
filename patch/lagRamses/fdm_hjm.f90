!################################################################
!################################################################
! Hamilton-Jacobi-Madelung (HJM) fluid solver for FDM
!
! Evolves (rho, S) on coarse AMR levels where lambda_dB >> dx.
! psi_re stores rho, psi_im stores S on fluid levels.
! Wave levels (>= fdm_first_wave_level) use existing Schrodinger solver.
!
! Reference: Kunkel+2024 (arXiv:2411.17288), Chan+2025 (arXiv:2504.10387)
!################################################################
!################################################################

!################################################################
! Main HJM step: SSP-RK3 time integration
! Called from fdm_step when fdm_use_hjm and ilevel < fdm_first_wave_level
!################################################################
subroutine fdm_hjm_step(ilevel, dt_loc)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dt_loc

  integer::igrid,ind,iskip,icell
  real(dp)::dx,scale,dx_loc
  integer::nx_loc
  integer::idim,icL,icR
  real(dp)::sqrho_c,sqrho_L,sqrho_R,d2,sum_d2,c1_cell,qp

  dx = 0.5d0**ilevel
  nx_loc = icoarse_max - icoarse_min + 1
  scale = boxlen / dble(nx_loc)
  dx_loc = dx * scale

  ! Fluid levels hold (rho, S) PERSISTENTLY: psi_re=rho, psi_im=S with S an
  ! unwrapped real action. No psi<->rhoS conversion per step — atan2 would
  ! wrap S to [-pi*hbar,pi*hbar] and break |dS| > pi*hbar flows (large boxes).
  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  ! ============================================================
  ! SSP-RK3 (Shu 2007):
  !   u1 = un + dt*L(un)
  !   u2 = 3/4*un + 1/4*u1 + 1/4*dt*L(u1)
  !   u^{n+1} = 1/3*un + 2/3*u2 + 2/3*dt*L(u2)
  !
  ! We use rho(1:ngrids*8) = psi_re as working array.
  ! rho_n/S_n are saved in f(:,1)/f(:,2) scratch, then the
  ! RK stages update psi_re/psi_im in-place via fdm_hjm_rk.
  ! ============================================================

  call fdm_hjm_rk(ilevel, dx_loc, dt_loc)

  ! --- Apply source kick to S: dS/dt += -Phi + QP ---
  ! Same operator-split stage as the gravitational source. When fdm_hjm_qp is
  ! on, we also add the Madelung quantum pressure QP = +(hbar^2/2) lap(sqrt rho)/sqrt rho.
  ! Sync rho ghosts first so the QP Laplacian reads post-RK neighbour amplitudes.
  if(fdm_hjm_qp) call make_virtual_fine_dp(psi_re(1), ilevel)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           psi_im(icell) = psi_im(icell) - phi(icell) * dt_loc
           if(fdm_hjm_qp) then
              sqrho_c = sqrt(max(psi_re(icell), 1.0d-8))
              sum_d2  = 0.0d0
              c1_cell = 0.0d0
              do idim=1,ndim
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icL)
                 call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icR)
                 ! son-guard Neumann mirror: missing or refined neighbour -> centre value
                 if(icL > 0 .and. son(icL) == 0) then
                    sqrho_L = sqrt(max(psi_re(icL), 1.0d-8))
                 else
                    sqrho_L = sqrho_c
                 end if
                 if(icR > 0 .and. son(icR) == 0) then
                    sqrho_R = sqrt(max(psi_re(icR), 1.0d-8))
                 else
                    sqrho_R = sqrho_c
                 end if
                 d2      = sqrho_R - 2.0d0*sqrho_c + sqrho_L
                 sum_d2  = sum_d2 + d2
                 c1_cell = max(c1_cell, abs(d2)/sqrho_c)
              end do
              ! validity gate: caustic cells (C1>max) belong to the wave solver
              if(c1_cell <= fdm_qp_c1max) then
                 qp = 0.5d0 * hbar_code**2 * sum_d2 / (dx_loc**2 * sqrho_c)
                 psi_im(icell) = psi_im(icell) + qp * dt_loc
              end if
           end if
        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_hjm_step

!################################################################
!################################################################
! SSP-RK3 integrator for HJM equations.
! Stores (rho_n, S_n) in allocated scratch arrays, applies 3 RK
! stages with ghost syncs between them.
! Loops directly over (ind, igrid) to avoid cell_list overhead.
!################################################################
subroutine fdm_hjm_rk(ilevel, dx_loc, dt_loc)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dx_loc,dt_loc

  real(dp),allocatable,dimension(:)::rho_n,S_n,drho,dS_arr
  integer::igrid,ind,iskip,icell
  integer::ntot
  real(dp)::rho1,S1,inv_a2,dx_inv,dx2_inv,hbar2_over_2

  inv_a2  = 1.0d0   ! supercomoving time absorbs all a factors
  dx_inv  = 1.0d0 / dx_loc
  dx2_inv = 1.0d0 / (dx_loc**2)
  hbar2_over_2 = 0.5d0 * hbar_code**2

  ! Allocate scratch arrays (same size as psi_re/psi_im)
  ntot = ncoarse + twotondim*ngridmax
  allocate(rho_n(1:ntot), S_n(1:ntot))
  allocate(drho(1:ntot), dS_arr(1:ntot))

  ! Save u^n
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rho_n(icell) = psi_re(icell)
           S_n(icell)   = psi_im(icell)
        end if
        igrid = next(igrid)
     end do
  end do

  ! --- Stage 1: u1 = un + dt*L(un) ---
  call fdm_hjm_rhs_grid(ilevel, dx_inv, dx2_inv, inv_a2, hbar2_over_2, drho, dS_arr)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           psi_re(icell) = rho_n(icell) + dt_loc * drho(icell)
           psi_im(icell) = S_n(icell)   + dt_loc * dS_arr(icell)
           if(psi_re(icell) < 1.0d-10) psi_re(icell) = 1.0d-10
        end if
        igrid = next(igrid)
     end do
  end do
  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  ! --- Stage 2: u2 = 3/4*un + 1/4*u1 + 1/4*dt*L(u1) ---
  call fdm_hjm_rhs_grid(ilevel, dx_inv, dx2_inv, inv_a2, hbar2_over_2, drho, dS_arr)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rho1 = psi_re(icell)
           S1   = psi_im(icell)
           psi_re(icell) = 0.75d0*rho_n(icell) + 0.25d0*rho1 &
                         + 0.25d0*dt_loc*drho(icell)
           psi_im(icell) = 0.75d0*S_n(icell)   + 0.25d0*S1 &
                         + 0.25d0*dt_loc*dS_arr(icell)
           if(psi_re(icell) < 1.0d-10) psi_re(icell) = 1.0d-10
        end if
        igrid = next(igrid)
     end do
  end do
  call make_virtual_fine_dp(psi_re(1), ilevel)
  call make_virtual_fine_dp(psi_im(1), ilevel)

  ! --- Stage 3: u^{n+1} = 1/3*un + 2/3*u2 + 2/3*dt*L(u2) ---
  call fdm_hjm_rhs_grid(ilevel, dx_inv, dx2_inv, inv_a2, hbar2_over_2, drho, dS_arr)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rho1 = psi_re(icell)
           S1   = psi_im(icell)
           psi_re(icell) = (1.0d0/3.0d0)*rho_n(icell) + (2.0d0/3.0d0)*rho1 &
                         + (2.0d0/3.0d0)*dt_loc*drho(icell)
           psi_im(icell) = (1.0d0/3.0d0)*S_n(icell)   + (2.0d0/3.0d0)*S1 &
                         + (2.0d0/3.0d0)*dt_loc*dS_arr(icell)
           if(psi_re(icell) < 1.0d-10) psi_re(icell) = 1.0d-10
        end if
        igrid = next(igrid)
     end do
  end do

  deallocate(rho_n, S_n, drho, dS_arr)

end subroutine fdm_hjm_rk

!################################################################
!################################################################
! RHS of HJM equations computed directly on AMR grid structure.
! Loops over (ind, igrid) and uses fdm_neighbor_cell.
!
! Continuity: drho/dt = -(1/a^2) * div(rho * grad(S))
! Hamilton-Jacobi: dS/dt = -(1/(2a^2))*|grad(S)|^2
!                        + (hbar^2/(2a^2)) * Q
!   Q = 0.5*lap(ln rho) + 0.25*|grad(ln rho)|^2
!################################################################
subroutine fdm_hjm_rhs_grid(ilevel, dx_inv, dx2_inv, inv_a2, hbar2_over_2, drho, dS_arr)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dx_inv,dx2_inv,inv_a2,hbar2_over_2
  real(dp),dimension(*),intent(out)::drho,dS_arr

  integer::igrid,ind,iskip,icell,idim
  integer::icL,icR
  real(dp)::rho_c,S_c
  real(dp)::rho_L,rho_R,S_L,S_R
  real(dp)::vel_L,vel_R
  real(dp)::flux_L,flux_R
  real(dp)::grad2S
  real(dp)::dS_bk,dS_fw

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           rho_c = psi_re(icell)
           S_c   = psi_im(icell)

           grad2S  = 0.0d0
           drho(icell) = 0.0d0

           do idim=1,ndim
              call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icL)
              call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icR)

              if(icL > 0 .and. son(icL) == 0) then
                 rho_L = psi_re(icL); S_L = psi_im(icL)
              else
                 rho_L = rho_c; S_L = S_c
              end if
              if(icR > 0 .and. son(icR) == 0) then
                 rho_R = psi_re(icR); S_R = psi_im(icR)
              else
                 rho_R = rho_c; S_R = S_c
              end if

              ! S is an unwrapped real action — plain differences are exact
              dS_bk = S_c - S_L
              dS_fw = S_R - S_c

              ! --- Continuity: upwind flux ---
              vel_L = dS_bk * dx_inv
              vel_R = dS_fw * dx_inv

              if(vel_R >= 0.0d0) then
                 flux_R = rho_c * vel_R
              else
                 flux_R = rho_R * vel_R
              end if
              if(vel_L >= 0.0d0) then
                 flux_L = rho_L * vel_L
              else
                 flux_L = rho_c * vel_L
              end if

              drho(icell) = drho(icell) - inv_a2 * (flux_R - flux_L) * dx_inv

              ! --- Hamilton-Jacobi: Sethian-Osher ---
              grad2S = grad2S + max(max(dS_bk*dx_inv, 0.0d0)**2, &
                                    min(dS_fw*dx_inv, 0.0d0)**2)

              ! --- Quantum pressure: log-density form ---
              ! QP omitted at fluid levels: by Madelung criterion C1<threshold,
              ! QP is O(C1/dx^2) — a small correction handled by wave solver
              ! at fine levels. Including it here destabilizes oct boundaries.
           end do

           dS_arr(icell) = inv_a2 * (-0.5d0 * grad2S)

        end if
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_hjm_rhs_grid

!################################################################
!################################################################
! Convert psi (Re, Im) -> (rho, S) for fluid levels
! After this: psi_re(icell) = rho = |psi|^2
!             psi_im(icell) = S   = hbar * atan2(Im, Re)
!################################################################
subroutine fdm_psi_to_rhoS(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell
  real(dp)::re_val,im_val

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        re_val = psi_re(icell)
        im_val = psi_im(icell)
        psi_re(icell) = re_val**2 + im_val**2          ! rho = |psi|^2
        psi_im(icell) = hbar_code * atan2(im_val, re_val)  ! S = hbar*theta
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_psi_to_rhoS

!################################################################
!################################################################
! Convert (rho, S) -> psi (Re, Im) for fluid levels
! After this: psi_re(icell) = sqrt(rho) * cos(S/hbar)
!             psi_im(icell) = sqrt(rho) * sin(S/hbar)
!################################################################
subroutine fdm_rhoS_to_psi(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell
  real(dp)::rho_val,S_val,amp,theta

  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        rho_val = psi_re(icell)
        S_val   = psi_im(icell)
        amp   = sqrt(max(rho_val, 0.0d0))
        theta = S_val / hbar_code
        psi_re(icell) = amp * cos(theta)
        psi_im(icell) = amp * sin(theta)
        igrid = next(igrid)
     end do
  end do

end subroutine fdm_rhoS_to_psi

!################################################################
!################################################################
! Fluid -> wave refinement criterion (detect Madelung breakdown).
!
! The fluid (rho,S) form is exact only in the single-stream regime.
! It fails at caustics / multi-streaming, where the wave solver is
! mandatory. We flag a fluid cell where either
!
!   (1) C_Q = dx^2 |lap sqrt(rho)| / sqrt(rho)  >  fdm_hjm_C1
!       Normalised amplitude curvature. Interference fringes raise it
!       ahead of the first node. Dimensionless, boost- and box-size
!       invariant (no m or dx in the threshold), and identical to the
!       C1 CFL-exclusion quantity. Primary detector.
!
!   (2) min_axis( d2S )  <  -dx^2 / (fdm_nla * dt)
!       Time-to-caustic look-ahead. v = grad S, so the second
!       difference of S is dx^2 dv/dx; a converging axis reaches a
!       caustic in ~ -1/(dv/dx). We refine fdm_nla steps early.
!       Sign-gated (converging axes only). Predictive backstop.
!
! A mass gate rho > 0.1*rho_mean drops dynamically irrelevant void
! interference. NOTE the old phase-Laplacian C2 = |d2S|/hbar was
! removed: its threshold hbar/(m dx^2) tests WAVE RESOLVABILITY, not
! fluid validity, so it fired on ordinary single-stream bulk
! convergence and flooded ~5% of the box.
!################################################################
subroutine fdm_madelung_refine_flag(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,idim,icL,icR
  integer::nflag_loc,nflag_cq_loc,nflag_cs_loc,nflag_tot,info
  integer::nflag_cq_tot,nflag_cs_tot
  integer::nx_loc
  real(dp)::sqrho_c,sqrho_L,sqrho_R,S_c,S_L,S_R
  real(dp)::lap_sq,d2S,d2S_min,CQ
  real(dp)::dx,scale,dx_loc,dt,dm_share,rho_gate,cs_thresh,rho_c
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(.not.use_fdm) return
  if(.not.fdm_use_hjm) return
  if(hbar_code <= 0.0d0) return
  if(ilevel >= fdm_first_wave_level) return

  ! Geometry and time step
  nx_loc = icoarse_max - icoarse_min + 1
  scale  = boxlen / dble(nx_loc)
  dx     = 0.5d0**ilevel
  dx_loc = dx * scale
  dt     = dtnew(ilevel)

  ! Mean fluid density: |psi|^2 averages the DM share when gas is present
  dm_share = 1.0d0
  if(hydro .and. omega_m > 0.0d0) dm_share = 1.0d0 - omega_b/omega_m
  ! Nonlinearity gate: C_Q must not fire on the linear cutoff-scale
  ! amplitude structure of the IC (delta ~ 0.1 at high z), only where the
  ! flow is entering the nonlinear / pre-caustic regime. rho > 2*rho_mean
  ! is delta > 1. This makes fluid->wave refinement density-driven, so it
  ! unlocks progressively with structure growth exactly as CDM does.
  rho_gate = 2.0d0 * dm_share

  ! Time-to-caustic threshold (v = grad S, look-ahead fdm_nla steps)
  if(dt > 0.0d0) then
     cs_thresh = -dx_loc*dx_loc / (fdm_nla * dt)
  else
     cs_thresh = -huge(1.0d0)
  end if

  nflag_loc = 0; nflag_cq_loc = 0; nflag_cs_loc = 0

  ! Fluid levels store (rho, S) directly: psi_re=rho, psi_im=S (unwrapped)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        rho_c = max(psi_re(icell),0.0d0)
        sqrho_c = sqrt(rho_c)

        if(sqrho_c > 1.0d-15 .and. rho_c > rho_gate) then
           S_c = psi_im(icell)
           lap_sq  = 0.0d0
           d2S_min = 0.0d0

           do idim=1,ndim
              call fdm_neighbor_cell(igrid, ilevel, ind, idim, 1, icL)
              call fdm_neighbor_cell(igrid, ilevel, ind, idim, 2, icR)

              if(icL > 0) then
                 sqrho_L = sqrt(max(psi_re(icL),0.0d0)); S_L = psi_im(icL)
              else
                 sqrho_L = sqrho_c; S_L = S_c
              end if
              if(icR > 0) then
                 sqrho_R = sqrt(max(psi_re(icR),0.0d0)); S_R = psi_im(icR)
              else
                 sqrho_R = sqrho_c; S_R = S_c
              end if

              ! isotropic 7-point Laplacian of sqrt(rho)
              lap_sq = lap_sq + (sqrho_L - sqrho_c) + (sqrho_R - sqrho_c)
              ! per-axis second difference of the action (converging = negative)
              d2S = S_L + S_R - 2.0d0*S_c
              if(d2S < d2S_min) d2S_min = d2S
           end do

           CQ = abs(lap_sq) / sqrho_c

           if(CQ >= fdm_hjm_C1 .or. d2S_min < cs_thresh) then
              flag1(icell) = 1
              nflag_loc = nflag_loc + 1
              if(CQ >= fdm_hjm_C1)      nflag_cq_loc = nflag_cq_loc + 1
              if(d2S_min < cs_thresh)   nflag_cs_loc = nflag_cs_loc + 1
           end if
        end if
        igrid = next(igrid)
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag_loc,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nflag_cq_loc,nflag_cq_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nflag_cs_loc,nflag_cs_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  nflag_tot = nflag_loc
  nflag_cq_tot = nflag_cq_loc
  nflag_cs_tot = nflag_cs_loc
#endif
  if(myid==1 .and. nstep_coarse_old < 80) then
     write(*,'(" MAD_FLAG lv=",I2," total=",I10," CQ=",I10," caustic=",I10)') &
          ilevel, nflag_tot, nflag_cq_tot, nflag_cs_tot
  end if

end subroutine fdm_madelung_refine_flag

!################################################################
!################################################################
! Auto-detect fdm_first_wave_level if set to 0
! Default: levelmin + 2 (2 fluid levels above base)
! Called from fdm_compute_hbar (initialization)
!################################################################
subroutine fdm_hjm_init()
  use amr_commons
  use fdm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(.not.fdm_use_hjm) return

  if(fdm_first_wave_level <= 0) then
     fdm_first_wave_level = levelmin + 2
  end if

  ! Clamp to valid range
  if(fdm_first_wave_level <= levelmin) fdm_first_wave_level = levelmin + 1
  if(fdm_first_wave_level > nlevelmax) fdm_first_wave_level = nlevelmax + 1

  if(myid == 1) then
     write(*,'(A)')       ' ============================================'
     write(*,'(A)')       ' HJM Hybrid FDM enabled'
     write(*,'(A,I3,A,I3)') '   Fluid levels: ', levelmin, ' to ', &
          fdm_first_wave_level - 1
     write(*,'(A,I3,A,I3)') '   Wave  levels: ', fdm_first_wave_level, &
          ' to ', nlevelmax
     write(*,'(A,ES9.2)')   '   C1 threshold: ', fdm_hjm_C1
     write(*,'(A,ES9.2)')   '   C2 threshold: ', fdm_hjm_C2
     write(*,'(A)')       ' ============================================'
  end if

end subroutine fdm_hjm_init

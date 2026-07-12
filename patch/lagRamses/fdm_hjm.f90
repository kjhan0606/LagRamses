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

  ! --- Apply potential kick: S -> S - Phi*dt ---
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        if(son(icell) == 0) then
           psi_im(icell) = psi_im(icell) - phi(icell) * dt_loc
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
! Madelung refinement criterion for fluid levels
! C1: quantum pressure curvature of sqrt(rho)
! C2: phase curvature of S
! If either exceeds threshold, flag cell for refinement
!################################################################
subroutine fdm_madelung_refine_flag(ilevel)
  use amr_commons
  use poisson_commons
  use fdm_commons
  implicit none
  integer,intent(in)::ilevel

  integer::igrid,ind,iskip,icell,idim,icL,icR
  integer::nflag_loc,nflag_C1_loc,nflag_C2_loc,nflag_tot,info
  integer::nflag_C1_tot,nflag_C2_tot
  real(dp)::sqrho_c,sqrho_L,sqrho_R,S_c,S_L,S_R
  real(dp)::C1_dim,C2_dim,C1_max,C2_max
  real(dp)::re_c,im_c,re_L,im_L,re_R,im_R
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if(.not.use_fdm) return
  if(.not.fdm_use_hjm) return
  if(hbar_code <= 0.0d0) return
  if(ilevel >= fdm_first_wave_level) return

  nflag_loc = 0
  nflag_C1_loc = 0
  nflag_C2_loc = 0

  ! Fluid levels store (rho, S) directly: psi_re=rho, psi_im=S (unwrapped)
  do ind=1,twotondim
     iskip = ncoarse + (ind-1)*ngridmax
     igrid = headl(myid, ilevel)
     do while(igrid > 0)
        icell = igrid + iskip
        sqrho_c = sqrt(max(psi_re(icell),0.0d0))

        if(sqrho_c > 1.0d-15) then
           C1_max = 0.0d0
           C2_max = 0.0d0
           S_c = psi_im(icell)

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

              C1_dim = abs(sqrho_R - 2.0d0*sqrho_c + sqrho_L) / sqrho_c
              C2_dim = abs(S_R - 2.0d0*S_c + S_L) / hbar_code

              C1_max = max(C1_max, C1_dim)
              C2_max = max(C2_max, C2_dim)
           end do

           if(C1_max >= fdm_hjm_C1 .or. C2_max >= fdm_hjm_C2) then
              flag1(icell) = 1
              nflag_loc = nflag_loc + 1
              if(C1_max >= fdm_hjm_C1) nflag_C1_loc = nflag_C1_loc + 1
              if(C2_max >= fdm_hjm_C2) nflag_C2_loc = nflag_C2_loc + 1
           end if
        end if
        igrid = next(igrid)
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nflag_loc,nflag_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nflag_C1_loc,nflag_C1_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nflag_C2_loc,nflag_C2_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  nflag_tot = nflag_loc
  nflag_C1_tot = nflag_C1_loc
  nflag_C2_tot = nflag_C2_loc
#endif
  if(myid==1 .and. nstep_coarse_old < 80) then
     write(*,'(" MAD_FLAG lv=",I2," total=",I10," C1=",I10," C2=",I10)') &
          ilevel, nflag_tot, nflag_C1_tot, nflag_C2_tot
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

!################################################################
!################################################################
! Dark Cooling Fine: Apply dark-sector cooling to DM particles
!
! For Atomic Dark Matter (aDM): each DM particle carries a dark
! internal energy edp. Per-particle cooling via backward Euler.
!
! Pattern follows cooling_fine.kjhan.f90 (grid traversal)
!################################################################
subroutine dark_cooling_fine(ilevel)
  use amr_commons
  use pm_commons
  use amr_parameters, only: adm_adiabatic
  use dark_cooling_mod, only: dark_cool_implicit
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel

  integer::ncache,igrid,ngrid
  real(dp)::aexp_next

  if(.not. use_adm) return
  if(numbtot(1,ilevel)==0) return
  if(verbose) write(*,111) ilevel

  aexp_next = aexp
  if(cosmo .and. adm_adiabatic) call adm_step_aexp(dtnew(ilevel),aexp_next)

  ncache = active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid) schedule(dynamic)
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     call sub_dark_cooling_fine(ilevel, igrid, ngrid, aexp_next)
  end do

111 format('   Entering dark_cooling_fine for level ',I2)
end subroutine dark_cooling_fine
!################################################################
!################################################################
subroutine sub_dark_cooling_fine(ilevel, igrid_start, ngrid, aexp_next)
  use amr_commons
  use pm_commons
  use amr_parameters, only: adm_adiabatic
  use dark_cooling_mod, only: dark_cool_implicit, dark_adiabatic_expand
#include "amr_index.h"
  implicit none

  integer,intent(in)::ilevel, igrid_start, ngrid
  real(dp),intent(in)::aexp_next

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_phys,dt_phys
  real(dp)::rho_D,n_D,mp_D_g
  real(dp),dimension(1:twotondim)::dm_mass_cell
  integer::nx_loc,i,ipart,jpart,ind,icell
  integer,dimension(1:nvector)::ind_grid
  real(dp)::edp_new,edp_phys
  integer,external::cell_index_from_part

  ! Physical constants
  real(dp),parameter::GeV_to_g = 1.78266192d-24

  ! Unit conversions
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Cell size at this level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max-icoarse_min+1)
  scale = boxlen/dble(nx_loc)
  dx_loc = dx*scale

  ! Physical leaf-cell volume [cm^3].  units() already returns the
  ! proper length conversion for the current cosmological epoch.
  vol_phys = (dx_loc*scale_l)**3

  ! Physical timestep [s]
  dt_phys = dtnew(ilevel)*scale_t

  ! Dark proton mass [g]
  mp_D_g = adm_mp * GeV_to_g

  ! Fill grid index array
  do i = 1, ngrid
     ind_grid(i) = active(ilevel)%igrid(igrid_start + i - 1)
  end do

  ! Traverse particles per grid.  All macro-particles in one leaf cell
  ! use the same local density, obtained from the summed DM mass in
  ! that leaf.
  do i = 1, ngrid
     dm_mass_cell(:) = 0.0d0

     ipart = headp(ind_grid(i))
     do jpart = 1, numbp(ind_grid(i))
        if(idp(ipart) > 0 .and. ptypep(ipart) /= PTYPE_STAR .and. &
             ptypep(ipart) /= PTYPE_SINK) then
           ind = cell_index_from_part(ipart, ind_grid(i), ilevel)
           icell = ICELL_OF(ind_grid(i), ind)
           if(son(icell) == 0) then
              dm_mass_cell(ind) = dm_mass_cell(ind) + &
                   mp(ipart) * scale_d * scale_l**3
           end if
        end if
        ipart = nextp(ipart)
     end do

     ipart = headp(ind_grid(i))
     do jpart = 1, numbp(ind_grid(i))
        ! Only DM-like particles (idp>0, ptypep is DM or excited iSIDM — not star/sink)
        if(idp(ipart) > 0 .and. ptypep(ipart) /= PTYPE_STAR .and. ptypep(ipart) /= PTYPE_SINK) then
           ind = cell_index_from_part(ipart, ind_grid(i), ilevel)
           icell = ICELL_OF(ind_grid(i), ind)
           if(son(icell) == 0) then
              rho_D = dm_mass_cell(ind) / vol_phys
              n_D = rho_D / mp_D_g
           else
              n_D = 0.0d0
           end if

           if(n_D > 0.0d0) then
              edp_phys = edp(ipart) * scale_v**2
              if(cosmo .and. adm_adiabatic) then
                 edp_phys = dark_adiabatic_expand(edp_phys,aexp,aexp_next)
              end if
              ! Apply dark cooling (edp in code energy/mass units)
              edp_new = dark_cool_implicit( &
                   edp_phys, &                 ! physical [erg/g]
                   rho_D, n_D, dt_phys, aexp)
              edp(ipart) = edp_new / scale_v**2  ! back to code units
           end if
        end if

        ipart = nextp(ipart)
     end do
  end do

end subroutine sub_dark_cooling_fine

!################################################################
! The ADM operator executes before update_time. Interpolate the same
! Friedmann table at the end of this level step so the a^-2 factor is exact
! for the code's cosmological time coordinate.
!################################################################
subroutine adm_step_aexp(dt_step,a_next)
  use amr_commons
  implicit none
  real(dp),intent(in)::dt_step
  real(dp),intent(out)::a_next
  real(dp)::t_next
  integer::i

  if(.not.cosmo) then
     a_next=aexp
     return
  end if
  t_next=t+dt_step
  i=1
  do while(tau_frw(i)>t_next .and. i<n_frw)
     i=i+1
  end do
  a_next=aexp_frw(i)*(t_next-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1)) + &
       & aexp_frw(i-1)*(t_next-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
end subroutine adm_step_aexp

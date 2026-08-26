!################################################################
! Atomic-dark-matter hydro-particle-mesh pressure force.
!
! rho_star is an already-reserved AMR scalar scratch field.  The HPM branch
! is rejected for star/sink runs, so it cannot overwrite stellar-density
! state.  Reusing it avoids adding an O(ncell) pressure allocation.
!################################################################
subroutine adm_hpm_force_fine(ilevel)
  use amr_commons
  use amr_parameters, only: adm_hpm,adm_hpm_gamma
  use pm_commons
  use poisson_commons, only: rho_star,f
  use adm_hpm_mod, only: adm_hpm_pressure
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  integer::igrid,ipart,jpart,ind,icell,i,icpu,ibound
  integer::ncache,nx_loc
  integer,external::cell_index_from_part
  real(dp)::dx,scale,dx_loc,vol_loc

  if(.not.adm_hpm) return
  if(numbtot(1,ilevel)==0) return

  dx=0.5d0**ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Start a fresh pressure field only on this level. Parent-level pressure
  ! remains available to the AMR interpolation used at coarse-fine faces.
  do ind=1,twotondim
     do i=1,active(ilevel)%ngrid
        rho_star(ICELL_OF(active(ilevel)%igrid(i),ind))=0.0d0
     end do
     do icpu=1,ncpu
        do i=1,reception(icpu,ilevel)%ngrid
           rho_star(ICELL_OF(reception(icpu,ilevel)%igrid(i),ind))=0.0d0
        end do
     end do
  end do

  ! NGP deposition is intentionally matched to the leaf-cell density used by
  ! dark_cooling_fine: the accumulator is rho_D e_D before conversion to P_D.
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ipart=headp(igrid)
     do jpart=1,numbp(igrid)
        if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
             & ptypep(ipart)/=PTYPE_SINK) then
           ind=cell_index_from_part(ipart,igrid,ilevel)
           icell=ICELL_OF(igrid,ind)
           if(son(icell)==0) rho_star(icell)=rho_star(icell)+mp(ipart)*edp(ipart)
        end if
        ipart=nextp(ipart)
     end do
  end do
  do ind=1,twotondim
     do i=1,active(ilevel)%ngrid
        icell=ICELL_OF(active(ilevel)%igrid(i),ind)
        rho_star(icell)=adm_hpm_pressure(1.0d0,rho_star(icell)/vol_loc, &
             & adm_hpm_gamma)
     end do
  end do

  ! Match the ordinary particle-density exchange: reverse accumulation first,
  ! then populate virtual grids for face gradients on MPI boundaries.
  call make_virtual_reverse_dp(rho_star(1),ilevel)
  call make_virtual_fine_dp(rho_star(1),ilevel)
  do ibound=1,nboundary
     do ind=1,twotondim
        do i=1,boundary(ibound,ilevel)%ngrid
           rho_star(ICELL_OF(boundary(ibound,ilevel)%igrid(i),ind))=0.0d0
        end do
     end do
  end do

  call adm_hpm_add_pressure_force(ilevel,dx_loc,vol_loc)
end subroutine adm_hpm_force_fine

!################################################################
subroutine adm_hpm_add_pressure_force(ilevel,dx_loc,vol_loc)
  use amr_commons
  use pm_commons
  use poisson_commons, only: rho_star,f
  use morton_hash
  use adm_hpm_mod, only: adm_hpm_acceleration
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dx_loc,vol_loc
  integer::igrid,ngrid,ncache,i,ind,idim,ipart,jpart,icell
  integer::ig_left,ig_right,ih_left,ih_right
  integer,external::cell_index_from_part
  real(dp)::rho_D,p_cen,p_left,p_right
  real(dp),dimension(1:twotondim)::dm_mass_cell
  integer,dimension(1:3,1:2,1:8)::ggg,hhh
  integer,dimension(1:nvector)::ind_grid_w,ind_cell_w
  integer,dimension(1:nvector,0:twondim)::igridn_w

  ncache=active(ilevel)%ngrid
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  do igrid=1,ncache,nvector
     ngrid=min(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid_w(i)=active(ilevel)%igrid(igrid+i-1)
        igridn_w(i,0)=ind_grid_w(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           igridn_w(i,2*idim-1)=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim-1)
           igridn_w(i,2*idim  )=morton_nbor_grid(ind_grid_w(i),ilevel,2*idim  )
        end do
     end do

     do i=1,ngrid
        dm_mass_cell(:)=0.0d0
        ipart=headp(ind_grid_w(i))
        do jpart=1,numbp(ind_grid_w(i))
           if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
                & ptypep(ipart)/=PTYPE_SINK) then
              ind=cell_index_from_part(ipart,ind_grid_w(i),ilevel)
              icell=ICELL_OF(ind_grid_w(i),ind)
              if(son(icell)==0) dm_mass_cell(ind)=dm_mass_cell(ind)+mp(ipart)
           end if
           ipart=nextp(ipart)
        end do

        do ind=1,twotondim
           ind_cell_w(i)=ICELL_OF(ind_grid_w(i),ind)
           rho_D=dm_mass_cell(ind)/vol_loc
           if(rho_D<=0.0d0 .or. son(ind_cell_w(i))/=0) cycle
           p_cen=rho_star(ind_cell_w(i))
           do idim=1,ndim
              ig_left=ggg(idim,1,ind); ig_right=ggg(idim,2,ind)
              ih_left=ICELL_OF(1,hhh(idim,1,ind))-1
              ih_right=ICELL_OF(1,hhh(idim,2,ind))-1
              if(igridn_w(i,ig_left)>0) then
                 p_left=rho_star(igridn_w(i,ig_left)+ih_left)
              else
                 call adm_hpm_sample_axis(ind_grid_w(i),ind,ilevel,idim,-1,p_cen,p_left)
              end if
              if(igridn_w(i,ig_right)>0) then
                 p_right=rho_star(igridn_w(i,ig_right)+ih_right)
              else
                 call adm_hpm_sample_axis(ind_grid_w(i),ind,ilevel,idim,1,p_cen,p_right)
              end if
              f(ind_cell_w(i),idim)=f(ind_cell_w(i),idim) + &
                   & adm_hpm_acceleration(p_left,p_right,rho_D,dx_loc)
           end do
        end do
     end do
  end do

  do idim=1,ndim
     call make_virtual_fine_dp(f(1,idim),ilevel)
  end do
end subroutine adm_hpm_add_pressure_force

!################################################################
! Sample pressure at a same-level neighbour, falling back to a parent-level
! CIC value at a coarse-fine interface.  This follows the force solver's
! scalar-field boundary convention and avoids a zero-pressure AMR edge.
!################################################################
subroutine adm_hpm_lookup_cell(ilevel,ix_in,iy_in,iz_in,value,found)
  use amr_commons
  use poisson_commons, only: rho_star
  use morton_hash
  use morton_keys, only: mkey_t,morton_encode
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  integer(8),intent(in)::ix_in,iy_in,iz_in
  real(dp),intent(out)::value
  logical,intent(out)::found
  integer(8)::ix,iy,iz,ncx,ncy,ncz,gx,gy,gz
  integer::igrid,ind,icell
  type(mkey_t)::key

  found=.false.; value=0.0d0
  if(.not.allocated(mort_table)) return
  if(ilevel<1 .or. ilevel>size(mort_table)) return
  ncx=int(nx,8)*2_8**ilevel; ncy=int(ny,8)*2_8**ilevel; ncz=int(nz,8)*2_8**ilevel
  ix=modulo(ix_in,ncx); iy=modulo(iy_in,ncy); iz=modulo(iz_in,ncz)
  gx=ix/2_8; gy=iy/2_8; gz=iz/2_8
  key=morton_encode(gx,gy,gz)
  igrid=morton_hash_lookup(mort_table(ilevel),key)
  if(igrid<=0) return
  ind=1+int(modulo(ix,2_8))+2*int(modulo(iy,2_8))+4*int(modulo(iz,2_8))
  icell=ICELL_OF(igrid,ind)
  value=rho_star(icell)
  found=.true.
end subroutine adm_hpm_lookup_cell

subroutine adm_hpm_sample_offset(igrid,ind,ilevel,ox,oy,oz,fallback,value)
  use amr_commons
  use morton_keys, only: mkey_t,grid_to_morton,morton_decode
  implicit none
  integer,intent(in)::igrid,ind,ilevel,ox,oy,oz
  real(dp),intent(in)::fallback
  real(dp),intent(out)::value
  integer(8)::gx,gy,gz,cx,cy,cz,tx,ty,tz,x0,y0,z0
  integer::bx,by,bz
  real(dp)::wx,wy,wz,w,val
  logical::found
  type(mkey_t)::key

  key=grid_to_morton(igrid,ilevel)
  call morton_decode(key,gx,gy,gz)
  cx=2_8*gx+int(mod(ind-1,2),8)
  cy=2_8*gy+int(mod((ind-1)/2,2),8)
  cz=2_8*gz+int((ind-1)/4,8)
  tx=cx+int(ox,8); ty=cy+int(oy,8); tz=cz+int(oz,8)
  call adm_hpm_lookup_cell(ilevel,tx,ty,tz,value,found)
  if(found) return
  if(ilevel<=1) then
     value=fallback
     return
  end if

  if(modulo(tx,2_8)==0_8) then; x0=tx/2_8-1_8; wx=0.75d0
  else; x0=(tx-1_8)/2_8; wx=0.25d0; end if
  if(modulo(ty,2_8)==0_8) then; y0=ty/2_8-1_8; wy=0.75d0
  else; y0=(ty-1_8)/2_8; wy=0.25d0; end if
  if(modulo(tz,2_8)==0_8) then; z0=tz/2_8-1_8; wz=0.75d0
  else; z0=(tz-1_8)/2_8; wz=0.25d0; end if
  value=0.0d0
  do bz=0,1
     do by=0,1
        do bx=0,1
           call adm_hpm_lookup_cell(ilevel-1,x0+bx,y0+by,z0+bz,val,found)
           if(.not.found) then
              value=fallback
              return
           end if
           w=merge(wx,1.0d0-wx,bx==1)*merge(wy,1.0d0-wy,by==1) &
                & *merge(wz,1.0d0-wz,bz==1)
           value=value+w*val
        end do
     end do
  end do
end subroutine adm_hpm_sample_offset

subroutine adm_hpm_sample_axis(igrid,ind,ilevel,idim,side,fallback,value)
  use amr_parameters, only: dp
  implicit none
  integer,intent(in)::igrid,ind,ilevel,idim,side
  real(dp),intent(in)::fallback
  real(dp),intent(out)::value
  integer::ox,oy,oz

  ox=0; oy=0; oz=0
  select case(idim)
  case(1); ox=side
  case(2); oy=side
  case(3); oz=side
  end select
  call adm_hpm_sample_offset(igrid,ind,ilevel,ox,oy,oz,fallback,value)
end subroutine adm_hpm_sample_axis

!################################################################
subroutine adm_hpm_timestep(ilevel)
  use amr_commons
  use amr_parameters, only: adm_hpm,adm_hpm_gamma,adm_hpm_courant
  use pm_commons
  use adm_hpm_mod, only: adm_hpm_sound_speed
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel
  integer::igrid,i,ipart,jpart,nx_loc,info
  real(dp)::dx,scale,dx_loc,cs,dt_loc,dt_all

  if(.not.adm_hpm) return
  if(numbtot(1,ilevel)==0) return
  dx=0.5d0**ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dt_loc=dtnew(ilevel)
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ipart=headp(igrid)
     do jpart=1,numbp(igrid)
        if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
             & ptypep(ipart)/=PTYPE_SINK) then
           cs=adm_hpm_sound_speed(edp(ipart),adm_hpm_gamma)
           if(cs>0.0d0) dt_loc=min(dt_loc,adm_hpm_courant*dx_loc/cs)
        end if
        ipart=nextp(ipart)
     end do
  end do
  dt_all=dt_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
  dtnew(ilevel)=min(dtnew(ilevel),dt_all)
end subroutine adm_hpm_timestep

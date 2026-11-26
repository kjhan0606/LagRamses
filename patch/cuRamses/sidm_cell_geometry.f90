! Production SIDM particle-to-child classification in RAMSES mesh coordinates.
integer function cell_index_from_part(ipart,igrid,ilevel)
  use pm_commons, only: xp
  use amr_commons, only: xg,boxlen,icoarse_min,icoarse_max, &
       & jcoarse_min,kcoarse_min,dp
  implicit none
  integer,intent(in)::ipart,igrid,ilevel
  integer::ix,iy,iz,nx_loc
  real(dp)::scale,skip_loc(3),x_mesh(3)

  ! move1 maps physical particle coordinates to the normalized mesh as
  ! xp/scale + skip_loc before comparing them with xg.
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  skip_loc=(/dble(icoarse_min),dble(jcoarse_min),dble(kcoarse_min)/)
  x_mesh=0.0d0
  x_mesh(1)=xp(ipart,1)/scale+skip_loc(1)
#if NDIM>1
  x_mesh(2)=xp(ipart,2)/scale+skip_loc(2)
#endif
#if NDIM==3
  x_mesh(3)=xp(ipart,3)/scale+skip_loc(3)
#endif

  ix=0
  iy=0
  iz=0
  if(x_mesh(1)>=xg(igrid,1))ix=1
#if NDIM>1
  if(x_mesh(2)>=xg(igrid,2))iy=1
#endif
#if NDIM==3
  if(x_mesh(3)>=xg(igrid,3))iz=1
#endif
  cell_index_from_part=1+ix+2*iy+4*iz
end function cell_index_from_part

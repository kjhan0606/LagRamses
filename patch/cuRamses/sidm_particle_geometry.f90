module sidm_particle_geometry
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: sidm_particle_subcell_index

contains

  pure integer function sidm_particle_subcell_index(xp_part,xg_center, &
       scale,skip_loc,ndim) result(subcell)
    real(real64),intent(in)::xp_part(:),xg_center(:),scale,skip_loc(:)
    integer,intent(in)::ndim
    integer::idim
    real(real64)::xcoord

    subcell = 1
    if(ndim<1 .or. ndim>size(xp_part) .or. &
         ndim>size(xg_center) .or. ndim>size(skip_loc)) then
       subcell = 0
       return
    end if

    do idim=1,ndim
       xcoord=xp_part(idim)/scale+skip_loc(idim)
       if(xcoord>=xg_center(idim)) subcell=subcell+2**(idim-1)
    end do
  end function sidm_particle_subcell_index

end module sidm_particle_geometry

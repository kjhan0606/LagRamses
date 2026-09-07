module snrt_agn_locator
  ! Local AMR ownership lookup for an xsink-style position.
  ! The returned cell is a leaf owned by this MPI rank, or zero when the
  ! source belongs to another rank.  The caller is responsible for using a
  ! per-update accretion increment rather than a sink/BH mass.
  use amr_parameters
  use amr_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none

contains

  subroutine snrt_agn_find_local_leaf(xsource, icell, ilevel_found, grid_hint, level_hint)
    real(dp), intent(in) :: xsource(1:ndim)
    integer, intent(out) :: icell, ilevel_found
    integer, intent(in), optional :: grid_hint, level_hint
    integer :: ilevel, i, igrid, idim, ind, child
    integer :: nx_loc
    integer :: first_level,last_level,ngrid_search
    logical :: hinted
    integer, dimension(1:ndim) :: skip_loc
    real(dp) :: dx, scale, q
    logical :: in_grid

    icell = 0
    ilevel_found = 0
    if(any(.not.ieee_is_finite(xsource)))return
    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc = 0
    if (ndim >= 1) skip_loc(1) = icoarse_min
    if (ndim >= 2) skip_loc(2) = jcoarse_min
    if (ndim >= 3) skip_loc(3) = kcoarse_min
    scale = boxlen / dble(nx_loc)
    hinted=present(grid_hint).and.present(level_hint)
    first_level=nlevelmax;last_level=levelmin
    if(hinted)then
       if(level_hint<levelmin.or.level_hint>nlevelmax)return
       if(grid_hint<1.or.grid_hint>size(xg,1))return
       first_level=level_hint;last_level=level_hint
    endif

    ! A source can be known on every rank after sink synchronization, while
    ! its hydrodynamic leaf is local to exactly one rank.
    do ilevel = first_level, last_level, -1
       dx = 0.5d0**ilevel
       ngrid_search=active(ilevel)%ngrid
       if(hinted)ngrid_search=1
       do i = 1, ngrid_search
          if(hinted)then
             ! Native stars already belong to an owned particle grid. Do
             ! not scan the entire mesh separately for every stellar source.
             igrid=grid_hint
          else
             igrid = active(ilevel)%igrid(i)
          endif
          ind = 0
          in_grid = .true.
          do idim = 1, ndim
             q = (xsource(idim) / scale + dble(skip_loc(idim)) - &
                  (xg(igrid,idim) - dx)) / dx
             ! INT truncates negative fractions toward zero, falsely making
             ! a source just outside a grid belong to that grid (and rank).
             if (q < 0d0 .or. q >= 2d0) then
                in_grid = .false.
                exit
             end if
             child = floor(q)
             ind = ind + child * 2**(idim-1)
          end do
          if (.not. in_grid) cycle

          icell = ICELL_OF(igrid,ind+1)
          if (son(icell) == 0) then
             ilevel_found = ilevel
             return
          end if
       end do
    end do
    icell = 0
  end subroutine snrt_agn_find_local_leaf

end module snrt_agn_locator

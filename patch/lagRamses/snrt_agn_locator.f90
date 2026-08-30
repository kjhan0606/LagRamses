module snrt_agn_locator
  ! Local AMR ownership lookup for an xsink-style position.
  ! The returned cell is a leaf owned by this MPI rank, or zero when the
  ! source belongs to another rank.  The caller is responsible for using a
  ! per-update accretion increment rather than a sink/BH mass.
  use amr_parameters
  use amr_commons
  implicit none

contains

  subroutine snrt_agn_find_local_leaf(xsource, icell, ilevel_found)
    real(dp), intent(in) :: xsource(1:ndim)
    integer, intent(out) :: icell, ilevel_found
    integer :: ilevel, i, igrid, idim, ind, child
    integer :: nx_loc
    integer, dimension(1:ndim) :: skip_loc
    real(dp) :: dx, scale, q
    logical :: in_grid

    icell = 0
    ilevel_found = 0
    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc = 0
    if (ndim >= 1) skip_loc(1) = icoarse_min
    if (ndim >= 2) skip_loc(2) = jcoarse_min
    if (ndim >= 3) skip_loc(3) = kcoarse_min
    scale = boxlen / dble(nx_loc)

    ! A source can be known on every rank after sink synchronization, while
    ! its hydrodynamic leaf is local to exactly one rank.
    do ilevel = nlevelmax, levelmin, -1
       dx = 0.5d0**ilevel
       do i = 1, active(ilevel)%ngrid
          igrid = active(ilevel)%igrid(i)
          ind = 0
          in_grid = .true.
          do idim = 1, ndim
             q = (xsource(idim) / scale + dble(skip_loc(idim)) - &
                  (xg(igrid,idim) - dx)) / dx
             child = int(q)
             if (child < 0 .or. child > 1) then
                in_grid = .false.
                exit
             end if
             ind = ind + child * 2**(idim-1)
          end do
          if (.not. in_grid) cycle

          icell = ncoarse + ind * ngridmax + igrid
          if (son(icell) == 0) then
             ilevel_found = ilevel
             return
          end if
       end do
    end do
    icell = 0
  end subroutine snrt_agn_find_local_leaf

end module snrt_agn_locator

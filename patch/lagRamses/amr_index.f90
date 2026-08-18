module amr_index
  use amr_commons, only: ncoarse, ngridmax, twotondim
  implicit none
contains
  pure elemental integer function icell_of(igrid, ichild)
    integer, intent(in) :: igrid, ichild
#ifdef AMR_INDEX_CHECK
    if (ichild < 1 .or. ichild > twotondim) error stop 'amr_index: invalid child'
    if (igrid < 1 .or. igrid > ngridmax) error stop 'amr_index: invalid grid'
#endif
    icell_of = ncoarse + (ichild-1)*ngridmax + igrid
  end function
  pure elemental integer function igrid_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) error stop 'amr_index: invalid cell'
#endif
    igrid_of = mod(icell - ncoarse - 1, ngridmax) + 1
  end function
  pure elemental integer function ichild_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) error stop 'amr_index: invalid cell'
#endif
    ichild_of = (icell - ncoarse - 1)/ngridmax + 1
  end function
  pure elemental logical function is_coarse_cell(icell)
    integer, intent(in) :: icell
    is_coarse_cell = icell <= ncoarse
  end function
end module

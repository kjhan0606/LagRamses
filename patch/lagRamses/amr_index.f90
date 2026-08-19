module amr_index
  use amr_commons, only: ncoarse, ngridmax, twotondim, amr_block_size
  implicit none
contains
  pure elemental integer function icell_of(igrid, ichild)
    integer, intent(in) :: igrid, ichild
#ifdef AMR_INDEX_CHECK
    if (ichild < 1 .or. ichild > twotondim) error stop 'amr_index: invalid child'
    if (igrid < 1 .or. igrid > ngridmax) error stop 'amr_index: invalid grid'
#endif
    icell_of = ncoarse + ((igrid-1)/amr_block_size)*(twotondim*amr_block_size) &
         & + (ichild-1)*amr_block_size + mod(igrid-1,amr_block_size) + 1
  end function
  pure elemental integer function igrid_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) error stop 'amr_index: invalid cell'
#endif
    igrid_of = ((icell-ncoarse-1)/(twotondim*amr_block_size))*amr_block_size &
         & + mod(mod(icell-ncoarse-1,twotondim*amr_block_size),amr_block_size) + 1
  end function
  pure elemental integer function ichild_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) error stop 'amr_index: invalid cell'
#endif
    ichild_of = (mod(icell-ncoarse-1,(twotondim*amr_block_size))/amr_block_size) + 1
  end function
  pure elemental logical function is_coarse_cell(icell)
    integer, intent(in) :: icell
    is_coarse_cell = icell <= ncoarse
  end function
end module

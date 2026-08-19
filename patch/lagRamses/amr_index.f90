module amr_index
  use amr_commons, only: ncoarse, ngridmax, twotondim, amr_block_size
  implicit none
contains
  pure elemental integer function icell_legacy(icell)
    integer, intent(in) :: icell
    integer :: g, c
    if (icell <= ncoarse) then
       icell_legacy = icell
    else
       c = (mod(icell-ncoarse-1, twotondim*amr_block_size)/amr_block_size) + 1
       g = (((icell-ncoarse-1)/(twotondim*amr_block_size))*amr_block_size) &
            & + mod(mod(icell-ncoarse-1, twotondim*amr_block_size), amr_block_size) + 1
       icell_legacy = ncoarse + (c-1)*ngridmax + g
    end if
  end function
  pure elemental integer function icell_of(igrid, ichild)
    integer, intent(in) :: igrid, ichild
#ifdef AMR_INDEX_CHECK
    if (ichild < 1 .or. ichild > twotondim) error stop 'amr_index: invalid child'
    if (igrid < 0 .or. igrid > ngridmax) error stop 'amr_index: invalid grid'
    if (igrid == 0) then
       icell_of = ncoarse
       return
    end if
#endif
    icell_of = ncoarse + ((igrid-1)/amr_block_size)*(twotondim*amr_block_size) &
         & + (ichild-1)*amr_block_size + mod(igrid-1,amr_block_size) + 1
  end function
  pure elemental integer function igrid_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) then
       igrid_of = icell
       return
    end if
#endif
    igrid_of = ((icell-ncoarse-1)/(twotondim*amr_block_size))*amr_block_size &
         & + mod(mod(icell-ncoarse-1,twotondim*amr_block_size),amr_block_size) + 1
  end function
  pure elemental integer function ichild_of(icell)
    integer, intent(in) :: icell
#ifdef AMR_INDEX_CHECK
    if (icell <= ncoarse) then
       ichild_of = icell
       return
    end if
#endif
    ichild_of = (mod(icell-ncoarse-1,(twotondim*amr_block_size))/amr_block_size) + 1
  end function
  pure elemental logical function is_coarse_cell(icell)
    integer, intent(in) :: icell
    is_coarse_cell = icell <= ncoarse
  end function
end module

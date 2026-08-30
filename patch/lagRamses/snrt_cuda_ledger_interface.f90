module snrt_cuda_ledger_interface
  use iso_c_binding, only: c_float, c_int, c_loc, c_ptr
  implicit none

  private
  public :: snrt_cuda_weighted_sum_fp32

  interface
     subroutine snrt_cuda_weighted_sum_fp32_c(directional, weights, scalar, &
          nrow, ndirection, ierr) bind(C, name='snrt_cuda_weighted_sum_fp32_c')
       import :: c_ptr, c_int
       type(c_ptr), value :: directional, weights, scalar
       integer(c_int), value :: nrow, ndirection
       integer(c_int) :: ierr
     end subroutine snrt_cuda_weighted_sum_fp32_c
  end interface

contains

  subroutine snrt_cuda_weighted_sum_fp32(directional, weights, scalar, nrow, &
       ndirection, ierr)
    real(c_float), contiguous, target, intent(in) :: directional(:), weights(:)
    real(c_float), contiguous, target, intent(out) :: scalar(:)
    integer, intent(in) :: nrow, ndirection
    integer, intent(out) :: ierr
    integer(c_int) :: c_ierr

    if (nrow <= 0 .or. ndirection <= 0 .or. size(directional) /= nrow * ndirection &
         .or. size(weights) /= ndirection .or. size(scalar) /= nrow) then
       ierr = 1
       return
    endif
    call snrt_cuda_weighted_sum_fp32_c(c_loc(directional(1)), c_loc(weights(1)), &
         c_loc(scalar(1)), int(nrow, c_int), int(ndirection, c_int), c_ierr)
    ierr = int(c_ierr)
  end subroutine snrt_cuda_weighted_sum_fp32

end module snrt_cuda_ledger_interface

module snrt_cuda_sparse_transport_interface
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  implicit none

  interface
     function snrt_cuda_upwind_sparse(state, direction, neighbor, ncell, &
          ndirection, cdt_over_dx) bind(C, name='snrt_cuda_upwind_sparse_c') &
          result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*)
       integer(c_int), intent(in) :: neighbor(*)
       integer(c_int), value :: ncell, ndirection
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_upwind_sparse
  end interface

end module snrt_cuda_sparse_transport_interface

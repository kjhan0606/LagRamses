module snrt_cuda_transport_interface
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  implicit none

  interface
     function snrt_cuda_upwind_periodic(state, direction, nx, ny, nz, &
          ndirection, cdt_over_dx) bind(C, name='snrt_cuda_upwind_periodic_c') &
          result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*)
       integer(c_int), value :: nx, ny, nz, ndirection
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_upwind_periodic
  end interface

end module snrt_cuda_transport_interface

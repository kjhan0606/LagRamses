module snrt_cuda_absorption_interface
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  implicit none

  interface
     function snrt_cuda_absorb(state, optical_depth, absorbed, ncell, &
          ndirection) bind(C, name='snrt_cuda_absorb_c') result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: optical_depth(*)
       real(c_float), intent(out) :: absorbed(*)
       integer(c_int), value :: ncell, ndirection
       integer(c_int) :: ierr
     end function snrt_cuda_absorb
  end interface

end module snrt_cuda_absorption_interface

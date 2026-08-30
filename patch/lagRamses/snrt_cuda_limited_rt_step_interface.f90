module snrt_cuda_limited_rt_step_interface
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  implicit none

  interface
     function snrt_cuda_transport_absorb_limited(state, direction, neighbor, &
          optical_depth, neutral_hydrogen, absorbed, ncell, ndirection, &
          cdt_over_dx) bind(C, name='snrt_cuda_transport_absorb_limited_c') &
          result(ierr)
       import :: c_int, c_float
       real(c_float), intent(inout) :: state(*)
       real(c_float), intent(in) :: direction(*), optical_depth(*), neutral_hydrogen(*)
       integer(c_int), intent(in) :: neighbor(*)
       real(c_float), intent(out) :: absorbed(*)
       integer(c_int), value :: ncell, ndirection
       real(c_float), value :: cdt_over_dx
       integer(c_int) :: ierr
     end function snrt_cuda_transport_absorb_limited
  end interface

end module snrt_cuda_limited_rt_step_interface

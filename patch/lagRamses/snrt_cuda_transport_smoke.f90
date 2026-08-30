program snrt_cuda_transport_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use snrt_cuda_transport_interface, only: snrt_cuda_upwind_periodic
  implicit none

  integer, parameter :: nx = 11, ny = 7, nz = 5, ndirection = 80
  integer, parameter :: ncell = nx * ny * nz
  integer(c_int) :: ierr
  integer :: icell, idir
  real(c_float) :: state(ncell, ndirection), direction(3, ndirection)
  real(c_double) :: sum_before, sum_after

  direction = 0.0_c_float
  do idir = 1, ndirection
     direction(1,idir) = 0.30_c_float
     direction(2,idir) = -0.20_c_float
     direction(3,idir) = 0.40_c_float
  end do

  state = 1.0_c_float
  ierr = snrt_cuda_upwind_periodic(state, direction, int(nx,c_int), int(ny,c_int), &
       int(nz,c_int), int(ndirection,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 1
  if (maxval(abs(real(state,c_double) - 1.0d0)) > 1.0d-6) error stop 2

  do idir = 1, ndirection
     do icell = 1, ncell
        state(icell,idir) = real(1 + modulo(17*icell + 13*idir, 29), c_float)
     end do
  end do
  sum_before = sum(real(state,c_double))
  ierr = snrt_cuda_upwind_periodic(state, direction, int(nx,c_int), int(ny,c_int), &
       int(nz,c_int), int(ndirection,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 3
  sum_after = sum(real(state,c_double))
  if (abs(sum_after-sum_before) / sum_before > 1.0d-6) error stop 4
  if (minval(state) < 0.0_c_float) error stop 5

  write(*,'(a,es14.6,a,es14.6)') 'SNRT_CUDA_TRANSPORT_OK relative_sum_error=', &
       abs(sum_after-sum_before) / sum_before, ' minimum=', real(minval(state),c_double)
end program snrt_cuda_transport_smoke

program snrt_cuda_absorption_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use snrt_cuda_absorption_interface, only: snrt_cuda_absorb
  implicit none

  integer, parameter :: ncell = 319, ndirection = 80
  integer :: icell, idir
  integer(c_int) :: ierr
  real(c_float) :: state(ncell,ndirection), tau(ncell), absorbed(ncell)
  real(c_double) :: before, after, removed, budget_error

  do icell = 1, ncell
     tau(icell) = real(modulo(7*icell,11),c_float) * 0.1_c_float
     do idir = 1, ndirection
        state(icell,idir) = real(1 + modulo(19*icell + 7*idir,31),c_float)
     end do
  end do
  before = sum(real(state,c_double))
  ierr = snrt_cuda_absorb(state, tau, absorbed, int(ncell,c_int), int(ndirection,c_int))
  if (ierr /= 0_c_int) error stop 1
  after = sum(real(state,c_double))
  removed = sum(real(absorbed,c_double))
  budget_error = abs(before-after-removed) / before
  if (budget_error > 1.0d-6) error stop 2
  if (minval(state) < 0.0_c_float .or. minval(absorbed) < 0.0_c_float) error stop 3

  tau = 0.0_c_float
  before = sum(real(state,c_double))
  ierr = snrt_cuda_absorb(state, tau, absorbed, int(ncell,c_int), int(ndirection,c_int))
  if (ierr /= 0_c_int) error stop 4
  if (abs(sum(real(state,c_double))-before) / before > 1.0d-7) error stop 5
  if (maxval(abs(absorbed)) > 1.0e-6_c_float) error stop 6

  write(*,'(a,es14.6)') 'SNRT_CUDA_ABSORPTION_OK relative_budget_error=', &
       budget_error
end program snrt_cuda_absorption_smoke

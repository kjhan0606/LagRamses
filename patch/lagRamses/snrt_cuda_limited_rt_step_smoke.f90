program snrt_cuda_limited_rt_step_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use snrt_cuda_limited_rt_step_interface, only: snrt_cuda_transport_absorb_limited
  implicit none

  integer, parameter :: nx = 9, ny = 7, nz = 5, ndirection = 80
  integer, parameter :: ncell = nx * ny * nz
  integer :: ix, iy, iz, icell, idir
  integer(c_int) :: ierr, neighbor(6,ncell)
  real(c_float) :: state(ncell,ndirection), direction(3,ndirection)
  real(c_float) :: tau(ncell), neutral_hydrogen(ncell), absorbed(ncell)
  real(c_double) :: before, after, removed

  do iz = 0, nz-1
     do iy = 0, ny-1
        do ix = 0, nx-1
           icell = 1 + ix + nx * (iy + ny * iz)
           neighbor(1,icell) = 1 + modulo(ix-1,nx) + nx * (iy + ny * iz)
           neighbor(2,icell) = 1 + modulo(ix+1,nx) + nx * (iy + ny * iz)
           neighbor(3,icell) = 1 + ix + nx * (modulo(iy-1,ny) + ny * iz)
           neighbor(4,icell) = 1 + ix + nx * (modulo(iy+1,ny) + ny * iz)
           neighbor(5,icell) = 1 + ix + nx * (iy + ny * modulo(iz-1,nz))
           neighbor(6,icell) = 1 + ix + nx * (iy + ny * modulo(iz+1,nz))
        end do
     end do
  end do
  do idir = 1, ndirection
     direction(:,idir) = (/0.30_c_float, -0.20_c_float, 0.40_c_float/)
  end do
  do icell = 1, ncell
     tau(icell) = 10.0_c_float
     neutral_hydrogen(icell) = real(1 + modulo(3*icell,5),c_float)
     do idir = 1, ndirection
        state(icell,idir) = real(1 + modulo(11*icell + 17*idir,23),c_float)
     end do
  end do
  before = sum(real(state,c_double))
  ierr = snrt_cuda_transport_absorb_limited(state, direction, neighbor, tau, &
       neutral_hydrogen, absorbed, int(ncell,c_int), int(ndirection,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 1
  after = sum(real(state,c_double))
  removed = sum(real(absorbed,c_double))
  if (abs(before-after-removed) / before > 1.0d-6) error stop 2
  if (maxval(real(absorbed,c_double) - real(neutral_hydrogen,c_double)) > 1.0d-5) error stop 3
  if (minval(state) < 0.0_c_float .or. minval(absorbed) < 0.0_c_float) error stop 4

  write(*,'(a,es14.6)') 'SNRT_CUDA_LIMITED_RT_STEP_OK relative_budget_error=', &
       abs(before-after-removed) / before
end program snrt_cuda_limited_rt_step_smoke

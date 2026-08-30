program snrt_source_transport_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection
  use snrt_angular_quadrature, only: snrt_angular_init
  use snrt_agn_source, only: snrt_agn_deposit_isotropic
  use snrt_cuda_sparse_transport_interface, only: snrt_cuda_upwind_sparse
  implicit none

  integer, parameter :: nx = 9, ny = 7, nz = 5, ncell = nx * ny * nz
  integer :: ix, iy, iz, icell, ierr
  integer(c_int) :: cuda_ierr, neighbor(6,ncell)
  real(dp) :: direction_dp(snrt_ndirection,3), weight(snrt_ndirection)
  real(dp) :: deposited_density, sum_before, sum_after
  real(c_float) :: state(snrt_ndirection,1,ncell)
  real(c_float) :: direction_c(3,snrt_ndirection), packed(ncell,snrt_ndirection)

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

  call snrt_angular_init(direction_dp, weight)
  direction_c = real(transpose(direction_dp),c_float)
  state = 0.0_c_float
  icell = 1 + nx/2 + nx * (ny/2 + ny * (nz/2))
  call snrt_agn_deposit_isotropic(state, icell, 1, 1.0d63, 125.0d0, 1.0d21, &
       1.0d-3, weight, deposited_density, ierr)
  if (ierr /= 0) error stop 1
  if (abs(sum(real(state,dp)) - deposited_density) > 1.0d-6) error stop 2

  do icell = 1, ncell
     packed(icell,:) = state(:,1,icell)
  end do
  sum_before = sum(real(packed,dp))
  cuda_ierr = snrt_cuda_upwind_sparse(packed, direction_c, neighbor, &
       int(ncell,c_int), int(snrt_ndirection,c_int), 0.05_c_float)
  if (cuda_ierr /= 0_c_int) error stop 3
  sum_after = sum(real(packed,dp))
  if (abs(sum_after-sum_before) / sum_before > 1.0d-6) error stop 4
  if (minval(packed) < 0.0_c_float) error stop 5

  write(*,'(a,es14.6,a,es14.6)') 'SNRT_SOURCE_TRANSPORT_OK density=', &
       deposited_density, ' relative_sum_error=', abs(sum_after-sum_before) / sum_before
end program snrt_source_transport_smoke

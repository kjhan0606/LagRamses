module snrt_angular_quadrature
  ! 3-D product S_N quadrature: 8 Gauss-Legendre polar ordinates times
  ! 10 uniform azimuthal ordinates.  It supplies 80 directions and exactly
  ! integrates a constant over 4*pi.
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection
  implicit none

  integer, parameter :: snrt_nmu = 8
  integer, parameter :: snrt_nphi = 10

contains

  subroutine snrt_angular_init(direction, weight)
    real(dp), intent(out) :: direction(snrt_ndirection, 3)
    real(dp), intent(out) :: weight(snrt_ndirection)
    real(dp), parameter :: mu(snrt_nmu) = (/ &
         -0.9602898564975363d0, -0.7966664774136267d0, &
         -0.5255324099163290d0, -0.1834346424956498d0, &
          0.1834346424956498d0,  0.5255324099163290d0, &
          0.7966664774136267d0,  0.9602898564975363d0 /)
    real(dp), parameter :: wmu(snrt_nmu) = (/ &
         0.1012285362903763d0, 0.2223810344533745d0, &
         0.3137066458778873d0, 0.3626837833783620d0, &
         0.3626837833783620d0, 0.3137066458778873d0, &
         0.2223810344533745d0, 0.1012285362903763d0 /)
    real(dp) :: phi, sin_theta, pi
    integer :: imu, iphi, idir

    if (snrt_ndirection /= snrt_nmu * snrt_nphi) error stop &
         'snrt_ndirection does not match the S_N quadrature'

    pi = acos(-1.0d0)
    idir = 0
    do imu = 1, snrt_nmu
       sin_theta = sqrt(max(0.0d0, 1.0d0 - mu(imu)**2))
       do iphi = 1, snrt_nphi
          idir = idir + 1
          phi = 2.0d0 * pi * (dble(iphi) - 0.5d0) / dble(snrt_nphi)
          direction(idir,1) = sin_theta * cos(phi)
          direction(idir,2) = sin_theta * sin(phi)
          direction(idir,3) = mu(imu)
          weight(idir) = wmu(imu) * 2.0d0 * pi / dble(snrt_nphi)
       end do
    end do
  end subroutine snrt_angular_init

end module snrt_angular_quadrature

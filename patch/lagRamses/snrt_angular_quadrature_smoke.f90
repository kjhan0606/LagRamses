program snrt_angular_quadrature_smoke
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection
  use snrt_angular_quadrature, only: snrt_angular_init
  implicit none

  real(dp) :: direction(snrt_ndirection, 3), weight(snrt_ndirection)
  real(dp) :: pi

  call snrt_angular_init(direction, weight)
  pi = acos(-1.0d0)
  if (abs(sum(weight) - 4.0d0*pi) > 1.0d-13) error stop 1
  if (maxval(abs(matmul(weight, direction))) > 1.0d-13) error stop 2
  if (maxval(abs(sum(direction**2, dim=2) - 1.0d0)) > 1.0d-13) error stop 3

  write(*,'(a,i0,a,es14.6)') 'SNRT_QUADRATURE_OK ndirection=', &
       snrt_ndirection, ' weight_sum=', sum(weight)
end program snrt_angular_quadrature_smoke

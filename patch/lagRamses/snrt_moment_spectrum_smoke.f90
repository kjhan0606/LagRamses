program snrt_moment_spectrum_smoke
  use snrt_moment_transport, only: mn_dp,mn_ok
  use snrt_moment_spectrum, only: mn_band_split_angular
  implicit none

  real(mn_dp),parameter :: lo=11.2_mn_dp,hi=13.6_mn_dp
  real(mn_dp) :: n(3),e(3),a(3),b(3),tiny_positive
  integer :: status,k

  ! Reproduce one representable step outside each band endpoint. The energy
  ! and number are subnormal, where epsilon*scale rounds the old tolerance to
  ! zero despite the discrepancy being only one ULP.
  tiny_positive=nearest(0._mn_dp,1._mn_dp)
  n=[64._mn_dp*tiny_positive,64._mn_dp*tiny_positive,2._mn_dp]
  e=[nearest(lo*n(1),-1._mn_dp),nearest(hi*n(2),1._mn_dp),24._mn_dp]
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status/=mn_ok)error stop 'one-ULP subnormal endpoint was rejected'
  if(any(a<0._mn_dp).or.any(b<0._mn_dp))error stop 'negative endpoint packet'
  if(a(1)/=0._mn_dp.or.b(1)/=n(1))error stop 'lower endpoint was not conserved'
  if(a(2)/=n(2).or.b(2)/=0._mn_dp)error stop 'upper endpoint was not conserved'
  if(any(a+b/=n))error stop 'split did not conserve photon number'
  if(8*epsilon(1._mn_dp)*max(e(1),hi*n(1))/(hi-lo)/=0._mn_dp) &
       error stop 'regression fixture does not underflow the former tolerance'

  ! More than the declared ULP envelope remains a hard spectral failure.
  e(1)=lo*n(1)
  do k=1,16
     e(1)=nearest(e(1),-1._mn_dp)
  enddo
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status==mn_ok)error stop 'multi-ULP subnormal band violation was accepted'

  ! The arithmetic allowance must not admit physically out-of-band energy.
  e(3)=lo*n(3)-1.e-8_mn_dp
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status==mn_ok)error stop 'finite lower-band violation was accepted'
  e(3)=hi*n(3)+1.e-8_mn_dp
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status==mn_ok)error stop 'finite upper-band violation was accepted'

  write(*,'(a)')'SNRT_MOMENT_SPECTRUM_ENDPOINTS_OK'
end program snrt_moment_spectrum_smoke

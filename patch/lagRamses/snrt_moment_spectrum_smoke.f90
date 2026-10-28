program snrt_moment_spectrum_smoke
  use snrt_moment_transport, only: mn_dp,mn_ok
  use snrt_moment_spectrum, only: mn_band_split_angular,mn_band_err_outside_edges,mn_band_err_zero_n_energy
  implicit none

  real(mn_dp),parameter :: lo=11.2_mn_dp,hi=13.6_mn_dp
  real(mn_dp) :: n(3),e(3),a(3),b(3),tiny_positive
  integer :: status,k

  ! Reproduce one representable step outside each band endpoint. The energy
  ! and number are subnormal, where epsilon*scale rounds the old tolerance to
  ! zero despite the discrepancy being only one ULP.
  tiny_positive=nearest(0._mn_dp,1._mn_dp)
  n=[4096._mn_dp*tiny_positive,4096._mn_dp*tiny_positive,2._mn_dp]
  e=[nearest(lo*n(1),-1._mn_dp),nearest(hi*n(2),1._mn_dp),24._mn_dp]
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status/=mn_ok)error stop 'one-ULP subnormal endpoint was rejected'
  if(any(a<0._mn_dp).or.any(b<0._mn_dp))error stop 'negative endpoint packet'
  if(a(1)/=0._mn_dp.or.b(1)/=n(1))error stop 'lower endpoint was not conserved'
  if(a(2)/=n(2).or.b(2)/=0._mn_dp)error stop 'upper endpoint was not conserved'
  if(any(a+b/=n))error stop 'split did not conserve photon number'
  if(8*epsilon(1._mn_dp)*max(e(1),hi*n(1))/(hi-lo)/=0._mn_dp) &
       error stop 'regression fixture does not underflow the former tolerance'

  ! Independent N/E transport can leave a meaningless subnormal ray: its
  ! apparent mean is below the 2 keV group edge. Keep N at the lower endpoint;
  ! the caller records the tiny E projection residual, never gas heating.
  n=[3*tiny_positive,0._mn_dp,2._mn_dp]
  e=[4000*tiny_positive,0._mn_dp,8000._mn_dp]
  call mn_band_split_angular(n,e,2000._mn_dp,10000._mn_dp,a,b,status)
  if(status/=mn_ok)error stop 'subnormal unresolved ray was rejected'
  if(a(1)/=0._mn_dp.or.b(1)/=n(1))error stop 'subnormal ray did not project to lower endpoint'
  if(10000._mn_dp*a(1)+2000._mn_dp*b(1)<=e(1)) &
       error stop 'subnormal projection did not retain its positive energy residual'
  n=[4096*tiny_positive,4096*tiny_positive,2._mn_dp]
  e=[lo*n(1),hi*n(2),24._mn_dp]

  ! The measured nine-ULP endpoint envelope is admitted, but a 1,024-ULP
  ! violation at a resolved packet remains a hard spectral failure.
  e(3)=lo*n(3)
  do k=1,9
     e(3)=nearest(e(3),-1._mn_dp)
  enddo
  call mn_band_split_angular(n(3:3),e(3:3),lo,hi,a(3:3),b(3:3),status)
  if(status/=mn_ok)error stop 'measured endpoint roundoff was rejected'
  if(a(3)/=0._mn_dp.or.b(3)/=n(3))error stop 'near-edge lower packet was not conserved'

  e(3)=lo*n(3)
  do k=1,1024
     e(3)=nearest(e(3),-1._mn_dp)
  enddo
  call mn_band_split_angular(n(3:3),e(3:3),lo,hi,a(3:3),b(3:3),status)
  if(status/=mn_band_err_outside_edges)error stop 'multi-ULP resolved band violation was not classified'

  e(1)=tiny_positive
  n(1)=0._mn_dp
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status/=mn_band_err_zero_n_energy)error stop 'zero-number energy was not classified'
  n(1)=4096._mn_dp*tiny_positive

  ! The arithmetic allowance must not admit physically out-of-band energy.
  e(3)=lo*n(3)-1.e-8_mn_dp
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status/=mn_band_err_outside_edges)error stop 'finite lower-band violation was not classified'
  e(3)=hi*n(3)+1.e-8_mn_dp
  call mn_band_split_angular(n,e,lo,hi,a,b,status)
  if(status/=mn_band_err_outside_edges)error stop 'finite upper-band violation was not classified'

  write(*,'(a)')'SNRT_MOMENT_SPECTRUM_ENDPOINTS_OK'
end program snrt_moment_spectrum_smoke

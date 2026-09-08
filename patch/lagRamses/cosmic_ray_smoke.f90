program cosmic_ray_smoke
  use cosmic_ray_physics
  implicit none
  logical::ok
  integer::ierr
  real(cr_dp)::energy,thermal,total,c2,cs,alpha_without,alpha_with
  call cr_validate(0,.true.,.false.,cr_gamma,ok)
  if(.not.ok)stop 1
  cr_enabled=.true.;cr_sn_fraction=.1d0;cr_snia_fraction=.2d0
  call cr_validate(1,.true.,.false.,cr_gamma,ok)
  if(.not.ok)stop 2
  call cr_source_partition(1d51,2d50,energy,ierr)
  if(ierr/=0.or.abs(energy-1.4d50)>1d-14*1d51)stop 3
  total=1.2d51;thermal=total-energy
  if(abs(thermal+energy-total)>1d-14*total.or.thermal<0)stop 4
  call cr_validate(0,.true.,.false.,cr_gamma,ok)
  if(ok)stop 5
  call cr_validate(1,.true.,.true.,cr_gamma,ok)
  if(ok)stop 6
  call cr_validate(1,.true.,.false.,5d0/3,ok)
  if(ok)stop 7
  call cr_source_partition(-1d0,0d0,energy,ierr)
  if(ierr==0)stop 8
  cs=1d0;c2=cr_support_speed2(9d0,2d0)
  if(abs(c2-2d0)>1d-14)stop 9
  alpha_without=5*cs;alpha_with=5*(cs+c2)
  if(alpha_with<=alpha_without)stop 10
  ! Compression Ecr~rho^gamma has the pressure derivative used by support.
  if(abs((cr_pressure(9d0*(1d0+1d-5)**cr_gamma)-cr_pressure(9d0))/2d-5-c2)>2d-5)stop 11
  cr_sn_fraction=1.01d0
  call cr_validate(1,.true.,.false.,cr_gamma,ok)
  if(ok)stop 12
  write(*,*)'CR_PRESSURE_SOURCE_PARTITION_SUPPORT_VALIDATION_PASS'
end program

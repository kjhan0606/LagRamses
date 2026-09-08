program dust_mass_smoke
  use dust_mass_physics
  implicit none
  integer::ierr
  real(real64)::d,e,q,a,b,half,full
  if(.not.dust_mass_parameters_ok())stop 1
  dust_mass_enabled=.true.
  call dust_condense([1d0,2d0,3d0],[.7d0,1d0,1d0],[.2d0,.5d0,1d0],d,ierr)
  if(ierr/=0.or.abs(d-.25d0)>1d-14)stop 2
  call dust_mass_step(1d0,.1d0,2d0,1d0,0d0,d,ierr)
  if(ierr/=0.or.abs(d-1/(1+9*exp(-2d0)))>1d-14)stop 3
  call dust_mass_step(1d0,.1d0,2d0,0d0,1d0,d,ierr)
  if(ierr/=0.or.abs(d-.1d0*exp(-2d0))>1d-14)stop 4
  call dust_mass_step(1d0,.1d0,2d0,1d0,1d0,d,ierr)
  if(ierr/=0.or.abs(d-.1d0/1.2d0)>1d-14)stop 5
  call dust_mass_step(1d0,.1d0,1d0,.7d0,.3d0,half,ierr)
  call dust_mass_step(1d0,half,1d0,.7d0,.3d0,d,ierr)
  call dust_mass_step(1d0,.1d0,2d0,.7d0,.3d0,full,ierr)
  if(ierr/=0.or.abs(d-full)>1d-14)stop 6
  call dust_mass_step(1d0,.1d0,1d4,1d0,0d0,d,ierr)
  if(ierr/=0.or.d/=1d0)stop 7
  call dust_mass_step(1d0,0d0,1d4,1d0,0d0,d,ierr)
  if(ierr/=0.or.d/=0d0)stop 8
  call dust_mass_step(1d0,1.1d0,1d0,1d0,0d0,d,ierr)
  if(ierr==0)stop 9
  call dust_mass_exchange(.1d0,.2d0,2d0,10d0,e,q,ierr)
  if(ierr/=0.or.e/=4.or.q/=2.or.(10-q+e)/=12)stop 10
  call dust_mass_exchange(.1d0,0d0,2d0,10d0,e,q,ierr)
  if(ierr/=0.or.e/=0.or.q/=-2)stop 11
  call dust_mass_exchange(.1d0,.2d0,2d0,1d0,e,q,ierr)
  if(ierr==0)stop 12
  call dust_mass_rates(1d-24,1d-26,100d0,a,b,ierr)
  if(ierr/=0.or.a<=0.or.b<0)stop 13
  call dust_mass_rates(1d-24,1d-26,2d6,a,b,ierr)
  if(ierr/=0.or.a/=0.or.abs(b-3*3.2d-18*(1d-24/dust_mp)*.5d0/dust_grain_radius_cm)>1d-25)stop 14
  write(*,*)'DUST_MASS_CONDENSATION_GROWTH_SPUTTERING_ENERGY_PASS'
end program

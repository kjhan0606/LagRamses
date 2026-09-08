! One-moment, fixed-characteristic-size dust mass closure. Metal fields
! denote total (gas+dust) metals; dust is a subset, not additional gas mass.
module dust_mass_physics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: dust_mass_dp=real64
  logical :: dust_mass_enabled=.false.,dust_growth=.true.,dust_sputtering=.true.
  real(real64) :: dust_condensation(3)=[0d0,.2d0,.15d0] ! wind, AGB, SNII; no Ia
  real(real64) :: dust_grain_radius_cm=1d-5,dust_grain_density=3d0
  real(real64) :: dust_sticking=.3d0,dust_growth_max_temperature=300d0
  real(real64) :: dust_metal_atom_mass=24d0,dust_injection_temperature=20d0
  real(real64),parameter :: dust_mp=1.67262192369d-24,dust_kb=1.380649d-16
contains
  logical function dust_mass_parameters_ok() result(ok)
    real(real64)::v(9)
    v=[dust_condensation,dust_grain_radius_cm,dust_grain_density,dust_sticking, &
       dust_growth_max_temperature,dust_metal_atom_mass,dust_injection_temperature]
    ok=all(ieee_is_finite(v)).and.all(v>=0)
    ok=ok.and.all(dust_condensation<=1).and.dust_sticking<=1
    ok=ok.and.dust_grain_radius_cm>0.and.dust_grain_density>0.and.dust_metal_atom_mass>0
    ok=ok.and.dust_growth_max_temperature>0.and.dust_injection_temperature>0
  end function

  function dust_mass_identity() result(v)
    real(real64)::v(13)
    v=[1d0,merge(1d0,0d0,dust_mass_enabled),merge(1d0,0d0,dust_growth), &
       merge(1d0,0d0,dust_sputtering),dust_condensation,dust_grain_radius_cm, &
       dust_grain_density,dust_sticking,dust_growth_max_temperature, &
       dust_metal_atom_mass,dust_injection_temperature]
  end function

  subroutine dust_condense(returned,hydrogen,helium,dust,ierr)
    real(real64),intent(in)::returned(3),hydrogen(3),helium(3)
    real(real64),intent(out)::dust
    integer,intent(out)::ierr
    real(real64)::metals(3)
    ierr=1;dust=0
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite(returned)).or..not.all(ieee_is_finite(hydrogen)).or. &
       .not.all(ieee_is_finite(helium)))return
    if(any(returned<0).or.any(hydrogen<0).or.any(helium<0))return
    metals=returned-hydrogen-helium
    if(any(metals < -1d-10*max(returned,tiny(1d0))))return
    if(dust_mass_enabled)dust=sum(dust_condensation*max(metals,0d0))
    if(.not.ieee_is_finite(dust))return
    ierr=0
  end subroutine

  subroutine dust_mass_rates(rho,metal_density,temperature,growth_rate,loss_rate,ierr)
    real(real64),intent(in)::rho,metal_density,temperature ! physical cgs
    real(real64),intent(out)::growth_rate,loss_rate
    integer,intent(out)::ierr
    real(real64)::speed,x
    growth_rate=0;loss_rate=0;ierr=1
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([rho,metal_density,temperature])))return
    if(min(rho,metal_density,temperature)<0.or.metal_density>rho)return
    if(dust_growth.and.temperature>0.and.temperature<=dust_growth_max_temperature)then
       speed=sqrt(8*dust_kb*temperature/(acos(-1d0)*dust_metal_atom_mass*dust_mp))
       ! Geometric collisions: pi*a^2*rho_Z*v*S / (4*pi*a^3*rho_s/3).
       ! The separate (1-D/Z) factor debits gas-phase metal availability.
       growth_rate=.75d0*dust_sticking*metal_density*speed/(dust_grain_density*dust_grain_radius_cm)
    endif
    if(dust_sputtering.and.temperature>0)then
       ! Tsai & Mathews erosion, McKinnon+2018 eq50; bulk mass lifetime a/(3|adot|).
       if(temperature<2d6)then
          x=(temperature/2d6)**2.5d0;x=x/(1+x)
       else
          x=1/(1+(2d6/temperature)**2.5d0)
       endif
       loss_rate=3*3.2d-18*(rho/dust_mp)*x/dust_grain_radius_cm
    endif
    if(.not.all(ieee_is_finite([growth_rate,loss_rate])))return
    ierr=0
  end subroutine

  subroutine dust_mass_step(metal,dust,dt,a,b,new_dust,ierr)
    ! Exact frozen-rate dD/dt = (A-B)D - A D^2/Z. Preserves D=0 and 0<=D<=Z.
    real(real64),intent(in)::metal,dust,dt,a,b
    real(real64),intent(out)::new_dust
    integer,intent(out)::ierr
    real(real64)::r,x,phi,decay,denominator
    ierr=1;new_dust=dust
    if(.not.all(ieee_is_finite([metal,dust,dt,a,b])))return
    if(min(metal,dust,dt,a,b)<0.or.dust>metal)return
    if(dust==0.or.dt==0.or.(a==0.and.b==0))then
       ierr=0;return
    endif
    r=a-b;x=r*dt
    if(.not.ieee_is_finite(x).or..not.ieee_is_finite(a*dt))return
    if(abs(x)<1d-5)then
       phi=dt*(1+x/2+x*x/6+x*x*x/24)
       new_dust=dust*(1+r*phi)/(1+a*(dust/metal)*phi)
    else if(r>0)then
       decay=exp(-x)
       denominator=decay+(a/r)*(dust/metal)*(1-decay)
       new_dust=dust/denominator
    else
       decay=exp(x);phi=(decay-1)/r
       new_dust=dust*decay/(1+a*(dust/metal)*phi)
    endif
    if(.not.ieee_is_finite(new_dust))return
    if(new_dust<0.or.new_dust>metal*(1+32*epsilon(1d0)))return
    new_dust=min(new_dust,metal);ierr=0
  end subroutine

  subroutine dust_mass_exchange(old_mass,new_mass,old_dust_energy,gas_thermal,new_energy,transfer,ierr)
    ! Fixed dust temperature during mass exchange; no latent heat model.
    ! transfer is charged to gas thermal energy, not CR/kinetic energy.
    real(real64),intent(in)::old_mass,new_mass,old_dust_energy,gas_thermal
    real(real64),intent(out)::new_energy,transfer
    integer,intent(out)::ierr
    new_energy=old_dust_energy;transfer=0;ierr=1
    if(.not.all(ieee_is_finite([old_mass,new_mass,old_dust_energy,gas_thermal])))return
    if(min(old_mass,new_mass,old_dust_energy,gas_thermal)<0)return
    if(old_mass==0)then
       if(new_mass/=0.or.old_dust_energy/=0)return
    else
       new_energy=old_dust_energy*(new_mass/old_mass)
       transfer=new_energy-old_dust_energy
    endif
    if(.not.ieee_is_finite(new_energy).or.transfer>gas_thermal)return
    ierr=0
  end subroutine
end module

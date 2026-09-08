! Relativistic CR fluid: pressure/work/advection use RAMSES NENER=1.
! CR energy is part of TOTAL energy, not thermal energy or a metal scalar.
! This closure is the trapped/advective limit, not diffusion/streaming.
module cosmic_ray_physics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: cr_dp=real64
  real(cr_dp),parameter :: cr_gamma=4d0/3d0
  logical :: cr_enabled=.false.,cr_sf_support=.false.
  real(cr_dp) :: cr_sn_fraction=0d0,cr_snia_fraction=0d0
  character(len=24) :: cr_transport='advective'
contains
  pure real(cr_dp) function cr_pressure(energy) result(p)
    real(cr_dp),intent(in)::energy
    p=(cr_gamma-1)*energy
  end function

  pure real(cr_dp) function cr_support_speed2(energy,density) result(c2)
    real(cr_dp),intent(in)::energy,density
    ! Adiabatic restoring term for CR trapped with isothermal gas.
    c2=cr_gamma*cr_pressure(energy)/density
  end function

  subroutine cr_validate(nener,hydro,gpu_hydro,gamma_nonthermal,ok)
    integer,intent(in)::nener
    logical,intent(in)::hydro,gpu_hydro
    real(cr_dp),intent(in)::gamma_nonthermal
    logical,intent(out)::ok
    ok=ieee_is_finite(cr_sn_fraction).and.ieee_is_finite(cr_snia_fraction)
    ok=ok.and.cr_sn_fraction>=0.and.cr_sn_fraction<=1.and.cr_snia_fraction>=0.and.cr_snia_fraction<=1
    if(.not.cr_enabled)then
       ok=ok.and.cr_sn_fraction==0.and.cr_snia_fraction==0.and..not.cr_sf_support
       return
    endif
    ok=ok.and.hydro.and.nener==1.and..not.gpu_hydro.and.trim(cr_transport)=='advective'
    ok=ok.and.ieee_is_finite(gamma_nonthermal).and.abs(gamma_nonthermal-cr_gamma)<1d-9
  end subroutine

  subroutine cr_source_partition(sn_energy,snia_coupled_energy,energy,ierr)
    real(cr_dp),intent(in)::sn_energy,snia_coupled_energy
    real(cr_dp),intent(out)::energy
    integer,intent(out)::ierr
    energy=0;ierr=1
    if(.not.ieee_is_finite(sn_energy).or..not.ieee_is_finite(snia_coupled_energy))return
    if(sn_energy<0.or.snia_coupled_energy<0)return
    if(.not.ieee_is_finite(cr_sn_fraction).or..not.ieee_is_finite(cr_snia_fraction))return
    if(min(cr_sn_fraction,cr_snia_fraction)<0.or.max(cr_sn_fraction,cr_snia_fraction)>1)return
    if(cr_enabled)energy=cr_sn_fraction*sn_energy+cr_snia_fraction*snia_coupled_energy
    if(.not.ieee_is_finite(energy))return
    ! Set a component of the existing total-energy increment; no extra energy.
    ierr=0
  end subroutine

  function cr_identity() result(values)
    real(cr_dp)::values(6)
    values=[1d0,merge(1d0,0d0,cr_enabled),cr_gamma,cr_sn_fraction,cr_snia_fraction,merge(1d0,0d0,cr_sf_support)]
  end function
end module

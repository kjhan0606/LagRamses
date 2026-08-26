program test_adm_thermal_evolution
  use amr_parameters
  use dark_cooling_mod, only: dark_adiabatic_expand, dark_net_cooling, &
       & dark_cool_implicit
  implicit none

  real(dp),parameter::kB_cgs=1.380649d-16
  real(dp),parameter::GeV_to_g=1.78266192d-24
  real(dp),parameter::alpha_em=7.2973525693d-3
  real(dp)::edp_old,edp_new,temperature_old,temperature_new
  real(dp)::rho_D,n_D,lambda

  ! A decoupled non-relativistic gas has T_D proportional to a^-2.
  adm_mp=40.0d0
  adm_T_floor=1.0d-6
  temperature_old=100.0d0
  edp_old=1.5d0*kB_cgs*temperature_old/(adm_mp*GeV_to_g)
  edp_new=dark_adiabatic_expand(edp_old,0.01d0,0.02d0)
  temperature_new=(2.0d0/3.0d0)*edp_new*adm_mp*GeV_to_g/kB_cgs
  if(abs(temperature_new/temperature_old-0.25d0)>1.0d-14) then
     write(*,*) 'FAIL: adiabatic a^-2 evolution',temperature_new
     stop 1
  end if

  ! In the ionised low-density limit a gas below the dark-radiation bath
  ! receives energy. This catches a regression to the old cooling-only sign.
  adm_alpha=alpha_em
  adm_mp=0.9382720813d0
  adm_me_ratio=0.5109989461d-3/adm_mp
  adm_xi=1.0d4
  n_D=1.0d-10
  rho_D=n_D*adm_mp*GeV_to_g
  temperature_old=1.0d4
  edp_old=1.5d0*kB_cgs*temperature_old/(adm_mp*GeV_to_g)
  lambda=dark_net_cooling(temperature_old,n_D,1.0d0)
  if(lambda>=0.0d0) then
     write(*,*) 'FAIL: Compton heating sign',lambda
     stop 2
  end if
  edp_new=dark_cool_implicit(edp_old,rho_D,n_D,1.0d10,1.0d0)
  if(edp_new<=edp_old) then
     write(*,*) 'FAIL: implicit Compton heating',edp_old,edp_new
     stop 3
  end if

  ! Empty AMR leaves are a valid no-op, not a division-by-zero path.
  edp_new=dark_cool_implicit(edp_old,0.0d0,0.0d0,1.0d10,1.0d0)
  if(edp_new/=edp_old) then
     write(*,*) 'FAIL: empty-cell thermal guard',edp_old,edp_new
     stop 4
  end if

  write(*,*) 'ADM thermal-evolution unit test passed'
end program test_adm_thermal_evolution

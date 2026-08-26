!################################################################
! Atomic-dark-matter hydro-particle-mesh (HPM) closure utilities.
!
! This is deliberately an approximate pressure closure, not a dark-fluid
! Riemann solver.  Its conservative variables remain the collisionless
! macro-particles; P_D=(gamma_D-1) rho_D e_D is evaluated on the AMR mesh.
!################################################################
module adm_hpm_mod
  use amr_parameters, only: dp
  implicit none
  private
  public :: adm_hpm_pressure, adm_hpm_acceleration, adm_hpm_sound_speed

contains

pure function adm_hpm_pressure(rho_D, edp_D, gamma_D) result(pressure_D)
  real(dp),intent(in)::rho_D,edp_D,gamma_D
  real(dp)::pressure_D

  pressure_D=max(0.0d0,(gamma_D-1.0d0)*rho_D*max(edp_D,0.0d0))
end function adm_hpm_pressure

pure function adm_hpm_acceleration(p_left,p_right,rho_D,dx) result(accel_D)
  real(dp),intent(in)::p_left,p_right,rho_D,dx
  real(dp)::accel_D

  if(rho_D>0.0d0 .and. dx>0.0d0) then
     accel_D=-(p_right-p_left)/(2.0d0*dx*rho_D)
  else
     accel_D=0.0d0
  end if
end function adm_hpm_acceleration

pure function adm_hpm_sound_speed(edp_D,gamma_D) result(cs_D)
  real(dp),intent(in)::edp_D,gamma_D
  real(dp)::cs_D

  cs_D=sqrt(max(0.0d0,gamma_D*(gamma_D-1.0d0)*edp_D))
end function adm_hpm_sound_speed

end module adm_hpm_mod

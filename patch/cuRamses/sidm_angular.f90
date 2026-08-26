!################################################################
! SIDM angular-scattering helpers
!################################################################
!
! These external routines are shared by the production collision operator
! and the standalone angular-kernel regression test.
subroutine sidm_sample_rutherford_cosine(uniform, eps2, cos_theta)
  implicit none

  real(kind(1.0d0)),intent(in)::uniform,eps2
  real(kind(1.0d0)),intent(out)::cos_theta
  real(kind(1.0d0))::cdf_lower,cdf_upper

  cdf_lower = 1.0d0/(2.0d0+eps2)
  cdf_upper = 1.0d0/eps2
  cos_theta = 1.0d0+eps2 &
       - 1.0d0/(cdf_lower+uniform*(cdf_upper-cdf_lower))
  cos_theta = max(-1.0d0,min(1.0d0,cos_theta))
end subroutine sidm_sample_rutherford_cosine
!################################################################
subroutine sidm_rotate_scattered_direction(v_rel_vec,v_rel_mag, &
     & cos_theta,sin_theta,phi,direction)
  implicit none

  real(kind(1.0d0)),intent(in)::v_rel_vec(3),v_rel_mag
  real(kind(1.0d0)),intent(in)::cos_theta,sin_theta,phi
  real(kind(1.0d0)),intent(out)::direction(3)
  real(kind(1.0d0))::v_hat(3),e1(3),e2(3),e1_mag

  v_hat = v_rel_vec/v_rel_mag
  if(abs(v_hat(3))<0.9d0) then
     e1(1) =  v_hat(2)
     e1(2) = -v_hat(1)
     e1(3) =  0.0d0
  else
     e1(1) =  0.0d0
     e1(2) =  v_hat(3)
     e1(3) = -v_hat(2)
  end if
  e1_mag = sqrt(sum(e1**2))
  e1 = e1/e1_mag
  e2(1) = v_hat(2)*e1(3)-v_hat(3)*e1(2)
  e2(2) = v_hat(3)*e1(1)-v_hat(1)*e1(3)
  e2(3) = v_hat(1)*e1(2)-v_hat(2)*e1(1)
  direction = sin_theta*cos(phi)*e1+sin_theta*sin(phi)*e2 &
       +cos_theta*v_hat
end subroutine sidm_rotate_scattered_direction

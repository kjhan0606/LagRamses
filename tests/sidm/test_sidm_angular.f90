program test_sidm_angular
  implicit none

  integer,parameter::ncase=3,n_sample=50000
  integer::case_index,index,seed_size
  integer,allocatable::seed(:)
  real(kind(1.0d0))::eps2(ncase),uniform,cos_theta,sin_theta,phi
  real(kind(1.0d0))::v_rel(3),direction(3),v_norm,direction_norm
  real(kind(1.0d0))::cos_check
  external sidm_sample_rutherford_cosine,sidm_rotate_scattered_direction

  eps2 = (/2.0d-3,2.0d-2,2.0d-1/)
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  do index=1,seed_size
     seed(index) = 97*index+23
  end do
  call random_seed(put=seed)

  v_rel = (/2.0d0,-3.0d0,5.0d0/)
  v_norm = sqrt(sum(v_rel**2))
  do case_index=1,ncase
     do index=1,n_sample
        call random_number(uniform)
        call sidm_sample_rutherford_cosine(uniform,eps2(case_index),cos_theta)
        sin_theta = sqrt(max(0.0d0,1.0d0-cos_theta**2))
        phi = 2.0d0*acos(-1.0d0)*dble(index)/dble(n_sample)
        call sidm_rotate_scattered_direction(v_rel,v_norm,cos_theta, &
             & sin_theta,phi,direction)
        direction_norm = sqrt(sum(direction**2))
        cos_check = sum(direction*v_rel)/v_norm
        if(abs(direction_norm-1.0d0)>2.0d-13 .or. &
             & abs(cos_check-cos_theta)>2.0d-13) error stop 1
        write(*,'(I2,1X,ES24.16)') case_index,cos_theta
     end do
  end do
end program test_sidm_angular

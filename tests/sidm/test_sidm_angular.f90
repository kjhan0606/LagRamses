program test_sidm_angular
  implicit none

  integer,parameter::ncase=5,n_sample=50000
  integer::case_index,index,seed_size
  integer,allocatable::seed(:)
  real(kind(1.0d0))::eps2(ncase),uniform,cos_theta,sin_theta,phi
  real(kind(1.0d0))::v_rel(3),direction(3),v_norm,direction_norm
  real(kind(1.0d0))::cos_check
  external sidm_sample_rutherford_cosine,sidm_rotate_scattered_direction

  eps2 = (/2.0d-3,2.0d-2,2.0d-1,2.0d0,2.0d1/)
  call verify_rotation_geometry()
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

contains

  subroutine verify_rotation_geometry()
    implicit none
    real(kind(1.0d0))::axes(3,10),cosines(5),azimuths(4),z_values(4)
    real(kind(1.0d0))::result(3),norm,cosine,theta_sine
    integer::axis,angle,azimuth

    axes=0d0
    axes(1,1)=1d0; axes(1,2)=-1d0
    axes(2,3)=1d0; axes(2,4)=-1d0
    axes(3,5)=1d0; axes(3,6)=-1d0
    z_values=(/0.899999d0,0.900001d0,-0.899999d0,-0.900001d0/)
    do axis=7,10
       axes(3,axis)=z_values(axis-6)
       axes(1,axis)=sqrt(1d0-axes(3,axis)**2)
    end do
    cosines=(/-1d0,-0.7d0,0d0,0.6d0,1d0/)
    azimuths=(/0d0,0.4d0,acos(-1d0),2d0*acos(-1d0)-0.2d0/)
    do axis=1,10
       norm=sqrt(sum(axes(:,axis)**2))
       do angle=1,5
          cosine=cosines(angle)
          theta_sine=sqrt(max(0d0,1d0-cosine**2))
          do azimuth=1,4
             call sidm_rotate_scattered_direction(axes(:,axis),norm, &
                  & cosine,theta_sine,azimuths(azimuth),result)
             if(abs(sum(result**2)-1d0)>2d-13) error stop 2
             if(abs(sum(result*axes(:,axis))/norm-cosine)>2d-13) error stop 3
          end do
       end do
    end do
  end subroutine verify_rotation_geometry
end program test_sidm_angular

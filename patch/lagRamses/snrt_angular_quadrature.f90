module snrt_angular_quadrature
  ! 3-D product S_N: Gauss-Legendre mu times uniform half-shifted phi.
  ! Default 8x10 is preserved bit-for-bit; refined families are 16x20/24x30.
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection, snrt_nmu, snrt_nphi
  implicit none

contains

  subroutine snrt_angular_init(direction, weight)
    real(dp), intent(out) :: direction(snrt_ndirection, 3)
    real(dp), intent(out) :: weight(snrt_ndirection)
    integer :: ierr
    call snrt_angular_product(snrt_nmu,snrt_nphi,direction,weight,ierr)
    if(ierr/=0)error stop 'invalid SNRT angular quadrature'
  end subroutine snrt_angular_init

  subroutine snrt_angular_product(nmu,nphi,direction,weight,ierr)
    ! Three fixed families only; changing the order changes checkpoint width.
    ! NIST DLMF 3.5: Legendre roots and 2/[(1-mu^2) Pn'(mu)^2].
    integer,intent(in)::nmu,nphi
    real(dp),intent(out)::direction(:,:),weight(:)
    integer,intent(out)::ierr
    real(dp)::mu(24),wmu(24),phi,sin_theta,pi,z,previous,p0,p1,p2,derivative
    integer::imu,iphi,idir,iteration,k
    ierr=1;direction=0;weight=0
    if(.not.((nmu==8.and.nphi==10).or.(nmu==16.and.nphi==20).or. &
         (nmu==24.and.nphi==30)))return
    if(size(weight)/=nmu*nphi.or.any(shape(direction)/=[nmu*nphi,3]))return
    pi = acos(-1.0d0)
    if(nmu==8)then
       call snrt_angular_legacy(direction,weight)
       ierr=0;return
    endif
       do imu=1,nmu/2
          z=cos(pi*(real(imu,dp)-.25d0)/(real(nmu,dp)+.5d0))
          do iteration=1,64
             p0=1;p1=z
             do k=2,nmu
                p2=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p2
             enddo
             derivative=nmu*(z*p1-p0)/(z*z-1)
             previous=z;z=previous-p1/derivative
             if(abs(z-previous)<=2*epsilon(z))exit
          enddo
          if(iteration>64)return
          ! Re-evaluate derivative at the accepted root for consistent weights.
          p0=1;p1=z
          do k=2,nmu
             p2=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p2
          enddo
          derivative=nmu*(z*p1-p0)/(z*z-1)
          mu(imu)=-z;mu(nmu+1-imu)=z
          wmu(imu)=2/((1-z*z)*derivative**2);wmu(nmu+1-imu)=wmu(imu)
       enddo
    idir = 0
    do imu = 1, nmu
       sin_theta = sqrt(max(0.0d0, 1.0d0 - mu(imu)**2))
       do iphi = 1, nphi
          idir = idir + 1
          phi = 2.0d0 * pi * (dble(iphi) - 0.5d0) / dble(nphi)
          direction(idir,1) = sin_theta * cos(phi)
          direction(idir,2) = sin_theta * sin(phi)
          direction(idir,3) = mu(imu)
          weight(idir) = wmu(imu) * 2.0d0 * pi / dble(nphi)
       end do
    end do
    ierr=0
  end subroutine snrt_angular_product

  ! Preserve the original fixed-bound loop as well as its literal nodes.
  ! A variable-bound equivalent loop changes compiler rounding at the last bit.
  subroutine snrt_angular_legacy(direction, weight)
    real(dp), intent(out) :: direction(80, 3)
    real(dp), intent(out) :: weight(80)
    real(dp), parameter :: mu(8) = (/ &
         -0.9602898564975363d0, -0.7966664774136267d0, &
         -0.5255324099163290d0, -0.1834346424956498d0, &
          0.1834346424956498d0,  0.5255324099163290d0, &
          0.7966664774136267d0,  0.9602898564975363d0 /)
    real(dp), parameter :: wmu(8) = (/ &
         0.1012285362903763d0, 0.2223810344533745d0, &
         0.3137066458778873d0, 0.3626837833783620d0, &
         0.3626837833783620d0, 0.3137066458778873d0, &
         0.2223810344533745d0, 0.1012285362903763d0 /)
    real(dp) :: phi, sin_theta, pi
    integer :: imu, iphi, idir

    if (80 /= 8 * 10) error stop &
         '80 does not match the S_N quadrature'

    pi = acos(-1.0d0)
    idir = 0
    do imu = 1, 8
       sin_theta = sqrt(max(0.0d0, 1.0d0 - mu(imu)**2))
       do iphi = 1, 10
          idir = idir + 1
          phi = 2.0d0 * pi * (dble(iphi) - 0.5d0) / dble(10)
          direction(idir,1) = sin_theta * cos(phi)
          direction(idir,2) = sin_theta * sin(phi)
          direction(idir,3) = mu(imu)
          weight(idir) = wmu(imu) * 2.0d0 * pi / dble(10)
       end do
    end do
  end subroutine snrt_angular_legacy

end module snrt_angular_quadrature

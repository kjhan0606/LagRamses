program snrt_angular_quadrature_smoke
  use amr_parameters, only: dp
  use snrt_state, only: snrt_ndirection
  use snrt_angular_quadrature, only: snrt_angular_init, snrt_angular_product
  implicit none

  real(dp) :: direction(snrt_ndirection, 3), weight(snrt_ndirection)
  real(dp) :: pi,axis(3),kappa,truth,flux_truth,error(3),moment_error(3),flux(3),second
  real(dp),allocatable::rays(:,:),w(:),lobe(:)
  integer::level,n,imu,iphi,opposite,d,j,ierr,nmu,nphi

  call snrt_angular_init(direction, weight)
  pi = acos(-1.0d0)
  if (abs(sum(weight) - 4.0d0*pi) > 1.0d-13) error stop 1
  if (maxval(abs(matmul(weight, direction))) > 1.0d-13) error stop 2
  if (maxval(abs(sum(direction**2, dim=2) - 1.0d0)) > 1.0d-13) error stop 3

  write(*,'(a,i0,a,es14.6)') 'SNRT_QUADRATURE_OK ndirection=', &
       snrt_ndirection, ' weight_sum=', sum(weight)
  axis=[1d0,2d0,3d0]/sqrt(14d0);kappa=20
  truth=2*pi*(1-exp(-2*kappa))/kappa
  flux_truth=truth*((1+exp(-2*kappa))/(1-exp(-2*kappa))-1/kappa)
  do level=0,2
     nmu=8*(level+1);nphi=10*(level+1);n=nmu*nphi
     allocate(rays(n,3),w(n),lobe(n))
     call snrt_angular_product(nmu,nphi,rays,w,ierr)
     if(ierr/=0.or.any(w<=0))error stop 4
     if(abs(sum(w)/(4*pi)-1)>2d-14)error stop 5
     if(maxval(abs(matmul(w,rays)))>2d-13)error stop 6
     if(maxval(abs(sum(rays*rays,dim=2)-1))>2d-14)error stop 7
     if(n==snrt_ndirection)then
        if(any(rays/=direction).or.any(w/=weight))error stop 8
     endif
     do d=1,3
        do j=1,3
           second=0;if(d==j)second=4*pi/3
           if(abs(sum(w*rays(:,d)*rays(:,j))-second)>2d-13)error stop 9
        enddo
     enddo
     do imu=1,nmu
        do iphi=1,nphi
           d=(imu-1)*nphi+iphi
           opposite=(nmu-imu)*nphi+modulo(iphi-1+nphi/2,nphi)+1
           if(maxval(abs(rays(d,:)+rays(opposite,:)))>2d-14)error stop 10
           if(w(d)/=w(opposite))error stop 11
        enddo
     enddo
     ! Oriented smooth narrow radiation lobe: analytic zeroth/first moments.
     ! Angular quadrature only, not a claim of spatial ray-effect convergence.
     lobe=exp(kappa*(matmul(rays,axis)-1))
     error(level+1)=abs(sum(w*lobe)/truth-1)
     flux=matmul(w*lobe,rays)
     moment_error(level+1)=sqrt(sum((flux-flux_truth*axis)**2))/flux_truth
     write(*,'(A,I4,A,2ES18.9)')'SNRT_ANGULAR_LOBE directions=',n, &
          ' number/flux_relative_error=',error(level+1),moment_error(level+1)
     deallocate(rays,w,lobe)
  enddo
  if(error(2)>=error(1).or.error(3)>=error(2).or.error(3)>1d-8)error stop 12
  if(moment_error(2)>=moment_error(1).or.moment_error(3)>=moment_error(2))error stop 13
  call snrt_angular_product(10,8,direction,weight,ierr)
  if(ierr==0.or.any(direction/=0).or.any(weight/=0))error stop 14
  write(*,'(A)')'SNRT_ANGULAR_REFINEMENT_SYMMETRY_MOMENTS_PASS'
end program snrt_angular_quadrature_smoke

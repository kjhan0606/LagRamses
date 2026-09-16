! Compact angular M_N transport. No RAMSES state, MPI or accelerator globals.
! The entropy is classical Boltzmann entropy on a positive spherical quadrature.
! u(1) is angular-integrated density; u(2:) are real spherical-harmonic moments.
! This is NOT a P_N polynomial reconstruction or a low-rank S_N array.
module snrt_moment_transport
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: mn_dp=real64, mn_ok=0, mn_bad_input=1, mn_no_closure=2, mn_bad_step=3
  type, public :: mn_basis
     integer :: order=0, nm=0, nq=0
     integer :: spatial_limiter=1 ! 1=minmod, 2=MC, 3=MC-beta(1.5), 4=MLP, 5=CWENO3
     real(mn_dp), allocatable :: direction(:,:), weight(:), harmonic(:,:)
  end type
  public :: mn_initialize, mn_project, mn_reconstruct, mn_flux, mn_advance, mn_material_exchange
  public :: mn_stress_tensor
  public :: mn_configure_spatial, mn_limited_slope, mn_vertex_gradient
  public :: mn_cweno_faces
contains
  subroutine mn_configure_spatial(b,name,ierr)
    type(mn_basis),intent(inout) :: b
    character(*),intent(in) :: name
    integer,intent(out) :: ierr
    ierr=mn_bad_input
    select case(trim(name))
    case('minmod');b%spatial_limiter=1
    case('mc');b%spatial_limiter=2
    case('mc15');b%spatial_limiter=3
    case('mlp');b%spatial_limiter=4
    case('cweno3');b%spatial_limiter=5
    case default;return
    end select
    ierr=mn_ok
  end subroutine

  pure subroutine mn_cweno_faces(stencil,face)
    ! Cell-average preserving 3-D CWENO: one quadratic central polynomial and
    ! eight octant linear polynomials. Coordinates are in cell widths. The
    ! quadratic basis x_a**2-1/12 and x_a*x_b has zero volume mean. Positive
    ! nonlinear weights sum to one, so reconstruction never changes the mean.
    ! A Bernstein convex-hull bound scales the WHOLE polynomial about that
    ! mean into [0,2*mean] throughout the cube; no pointwise clipping or floor.
    ! This supplies a conservative sufficient Euler bound at the existing
    ! c*dt/dx<=1/12. Closure success still requires a separate stage check.
    ! For constant face normal and linear kinetic advection the full quadratic
    ! face integral is analytic; cross terms enter weights/limiting but their
    ! full-face integrals vanish. A face-centre sample would NOT be exact.
    real(mn_dp),intent(in) :: stencil(-1:1,-1:1,-1:1)
    real(mn_dp),intent(out) :: face(6) ! -x,+x,-y,+y,-z,+z, face averages
    real(mn_dp) :: s(-1:1,-1:1,-1:1),p(9,0:8),beta(0:8),w(0:8),c(9)
    real(mn_dp),parameter :: bx(0:2)=[-.5d0,0d0,.5d0],bq(0:2)=[1d0/6,-1d0/3,1d0/6]
    real(mn_dp) :: scale,v,lo,hi,value,theta
    integer :: sx,sy,sz,k,ix,iy,iz,a
    face=0
    scale=maxval(stencil)
    if(scale<=0.or.stencil(0,0,0)<=0)return
    s=stencil/scale;v=s(0,0,0);p=0
    p(1:3,0)=.5d0*[s(1,0,0)-s(-1,0,0),s(0,1,0)-s(0,-1,0),s(0,0,1)-s(0,0,-1)]
    ! Central polynomial is (P_opt-sum(d_k P_k))/d_0 with d_0=1/2,
    ! d_1..8=1/16. Its linear coefficients equal those of P_opt.
    p(4:6,0)=[s(1,0,0)-2*v+s(-1,0,0),s(0,1,0)-2*v+s(0,-1,0),s(0,0,1)-2*v+s(0,0,-1)]
    p(7,0)=.5d0*(s(1,1,0)-s(1,-1,0)-s(-1,1,0)+s(-1,-1,0))
    p(8,0)=.5d0*(s(1,0,1)-s(1,0,-1)-s(-1,0,1)+s(-1,0,-1))
    p(9,0)=.5d0*(s(0,1,1)-s(0,1,-1)-s(0,-1,1)+s(0,-1,-1))
    k=0
    do sz=-1,1,2
       do sy=-1,1,2
          do sx=-1,1,2
             k=k+1;p(1:3,k)=[sx*(s(sx,0,0)-v),sy*(s(0,sy,0)-v),sz*(s(0,0,sz)-v)]
          enddo
       enddo
    enddo
    do k=0,8
       beta(k)=sum(p(1:3,k)**2)+(13d0/3)*sum(p(4:6,k)**2)+(7d0/6)*sum(p(7:9,k)**2)
    enddo
    w=1d0/16;w(0)=.5d0
    w=w/(1d-12+beta)**2;w=w/sum(w);c=matmul(p,w)
    lo=v;hi=v
    do iz=0,2
       do iy=0,2
          do ix=0,2
             value=v+c(1)*bx(ix)+c(2)*bx(iy)+c(3)*bx(iz)+ &
                  c(4)*bq(ix)+c(5)*bq(iy)+c(6)*bq(iz)+ &
                  c(7)*bx(ix)*bx(iy)+c(8)*bx(ix)*bx(iz)+c(9)*bx(iy)*bx(iz)
             lo=min(lo,value);hi=max(hi,value)
          enddo
       enddo
    enddo
    theta=1
    if(lo<0)theta=min(theta,v/(v-lo))
    if(hi>2*v)theta=min(theta,v/(hi-v))
    if(theta<1)theta=theta*(1-16*epsilon(1d0))
    c=theta*c
    do a=1,3
       face(2*a-1)=scale*(v-.5d0*c(a)+c(a+3)/6)
       face(2*a)=scale*(v+.5d0*c(a)+c(a+3)/6)
    enddo
  end subroutine

  elemental real(mn_dp) function mn_limited_slope(a,b,limiter) result(c)
    real(mn_dp),intent(in) :: a,b
    integer,intent(in) :: limiter
    c=minmod(a,b)
    if(limiter==2)c=minmod(.5d0*a+.5d0*b,2*c)
    ! beta=1.5 leaves every reconstructed face >= one quarter of its positive
    ! cell average; full MC can make a face exactly zero. No intensity floor.
    if(limiter==3)c=minmod(.5d0*a+.5d0*b,1.5d0*c)
  end function

  pure subroutine mn_vertex_gradient(stencil,gradient)
    ! Cartesian vertex-neighbour multidimensional limiting. Each vertex uses
    ! the eight cell averages sharing it, not three independent 1-D bounds.
    ! One theta scales the full central gradient. The additional corner bound
    ! keeps the reconstructed linear intensity >= I_cell/4 throughout the
    ! cell, without adding photons. Caller supplies finite nonnegative values.
    real(mn_dp),intent(in) :: stencil(-1:1,-1:1,-1:1)
    real(mn_dp),intent(out) :: gradient(3)
    real(mn_dp) :: value,theta,lo,hi,change,excursion
    integer :: sx,sy,sz,ix,iy,iz
    value=stencil(0,0,0)
    if(value==0)then
       gradient=0;return
    endif
    gradient=.5d0*[stencil(1,0,0)-stencil(-1,0,0), &
         stencil(0,1,0)-stencil(0,-1,0),stencil(0,0,1)-stencil(0,0,-1)]
    if(all(gradient==0))return
    theta=1
    do sz=-1,1,2
       do sy=-1,1,2
          do sx=-1,1,2
             lo=value;hi=value
             do iz=0,1
                do iy=0,1
                   do ix=0,1
                      lo=min(lo,stencil(ix*sx,iy*sy,iz*sz))
                      hi=max(hi,stencil(ix*sx,iy*sy,iz*sz))
                   enddo
                enddo
             enddo
             change=.5d0*(sx*gradient(1)+sy*gradient(2)+sz*gradient(3))
             if(change>0)theta=min(theta,(hi-value)/change)
             if(change<0)theta=min(theta,(lo-value)/change)
          enddo
       enddo
    enddo
    excursion=.5d0*sum(abs(gradient))
    if(excursion>0)theta=min(theta,.75d0*value/excursion)
    gradient=max(0d0,theta)*gradient
  end subroutine

  subroutine mn_stress_tensor(b,u,tensor,ierr)
    ! Local orthonormal frame: R00=E, R0i=F_i/c, Rij=P_ij.
    ! Higher harmonics do not add independent entries to this 4x4 tensor.
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(in) :: u(:)
    real(mn_dp), intent(out) :: tensor(4,4)
    integer, intent(out) :: ierr
    real(mn_dp) :: angular(b%nq)
    real(mn_dp) :: vector(4)
    integer :: q,j
    tensor=0
    call mn_reconstruct(b,u,angular,ierr)
    if(ierr/=mn_ok)return
    if(b%nm==4)then
       do q=1,b%nq
          vector=[1d0,b%direction(:,q)]
          do j=1,4
             tensor(:,j)=tensor(:,j)+b%weight(q)*angular(q)*vector*vector(j)
          enddo
       enddo
       return
    endif
    tensor(1,1)=u(1)
    tensor(1,2:4)=[-u(3),-u(4),u(2)]/sqrt(3d0)
    tensor(2:4,1)=tensor(1,2:4)
    tensor(4,4)=(u(1)+2*u(5)/sqrt(5d0))/3
    tensor(2,2)=(u(1)-tensor(4,4))/2+u(8)/sqrt(15d0)
    tensor(3,3)=(u(1)-tensor(4,4))/2-u(8)/sqrt(15d0)
    tensor(2,3)=u(9)/sqrt(15d0);tensor(3,2)=tensor(2,3)
    tensor(2,4)=-u(6)/sqrt(15d0);tensor(4,2)=tensor(2,4)
    tensor(3,4)=-u(7)/sqrt(15d0);tensor(4,3)=tensor(3,4)
  end subroutine
  subroutine mn_initialize(b,order,ierr,angular_mu,angular_phi)
    type(mn_basis), intent(out) :: b
    integer, intent(in) :: order
    integer, intent(out) :: ierr
    integer, optional, intent(in) :: angular_mu,angular_phi
    ! More azimuthal nodes than 2*order: the legacy 10-phi rule is singular
    ! for one of the l=5 harmonics. This rule is independent of stored width.
    integer :: nmu,nphi
    real(mn_dp), allocatable :: mu(:),wm(:)
    real(mn_dp) :: z,previous,p0,p1,p2,derivative,phi,pi,s,leg(0:5,0:5),norm
    integer :: i,j,k,it,q,l,m,h
    ierr=mn_bad_input
    if(order<1.or.order>5)return
    nmu=16;nphi=24
    if(present(angular_mu))nmu=angular_mu
    if(present(angular_phi))nphi=angular_phi
    if(nmu<order+1.or.nphi<=2*order.or.mod(nmu,2)/=0)return
    if(nmu>64.or.nphi>128)return
    allocate(mu(nmu),wm(nmu))
    pi=acos(-1d0)
    do i=1,nmu/2
       z=cos(pi*(i-.25d0)/(nmu+.5d0))
       do it=1,64
          p0=1;p1=z
          do k=2,nmu
             p2=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p2
          enddo
          derivative=nmu*(z*p1-p0)/(z*z-1)
          previous=z;z=z-p1/derivative
          if(abs(z-previous)<4*epsilon(z))exit
       enddo
       if(it>64)return
       p0=1;p1=z
       do k=2,nmu
          p2=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=p2
       enddo
       derivative=nmu*(z*p1-p0)/(z*z-1)
       mu(i)=-z;mu(nmu+1-i)=z
       wm(i)=1/((1-z*z)*derivative**2);wm(nmu+1-i)=wm(i)
    enddo
    b%order=order;b%nm=(order+1)**2;b%nq=nmu*nphi
    allocate(b%direction(3,b%nq),b%weight(b%nq),b%harmonic(b%nm,b%nq))
    q=0
    do i=1,nmu
       z=mu(i);s=sqrt(1-z*z);leg=0;leg(0,0)=1
       do m=1,order
          leg(m,m)=-(2*m-1)*s*leg(m-1,m-1)
       enddo
       do m=0,order-1
          leg(m+1,m)=(2*m+1)*z*leg(m,m)
          do l=m+2,order
             leg(l,m)=((2*l-1)*z*leg(l-1,m)-(l+m-1)*leg(l-2,m))/(l-m)
          enddo
       enddo
       do j=1,nphi
          q=q+1;phi=2*pi*(j-.5d0)/nphi
          b%direction(:,q)=[s*cos(phi),s*sin(phi),z];b%weight(q)=wm(i)/nphi
          h=0
          do l=0,order
             h=h+1;b%harmonic(h,q)=sqrt(real(2*l+1,mn_dp))*leg(l,0)
             do m=1,l
                norm=2*(2*l+1)
                do k=l-m+1,l+m
                   norm=norm/k
                enddo
                norm=sqrt(norm)*leg(l,m)
                h=h+1;b%harmonic(h,q)=norm*cos(m*phi)
                h=h+1;b%harmonic(h,q)=norm*sin(m*phi)
             enddo
          enddo
       enddo
    enddo
    b%weight=b%weight/sum(b%weight)
    ierr=mn_ok
  end subroutine

  subroutine mn_project(b,intensity,u,ierr)
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(in) :: intensity(:)
    real(mn_dp), intent(out) :: u(:)
    integer, intent(out) :: ierr
    ierr=mn_bad_input
    if(b%nm==0.or.size(intensity)/=b%nq.or.size(u)/=b%nm)return
    if(any(.not.ieee_is_finite(intensity)))return
    if(any(intensity<0))return
    u=matmul(b%harmonic,b%weight*intensity)
    if(any(.not.ieee_is_finite(u)))return
    ierr=mn_ok
  end subroutine

  subroutine mn_reconstruct(b,u,intensity,ierr,residual,diagnostic,dual_hint,stream_cache,max_iterations)
    ! Normalized dual Newton solve, avoiding overflow with log-sum-exp.
    ! No clipping of intensity, isotropic moment mixing or hidden SN fallback.
    ! Vacuum is exact; non-realizable/boundary states fail without mutation.
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(in) :: u(:)
    real(mn_dp), intent(out) :: intensity(:)
    integer, intent(out) :: ierr
    real(mn_dp), optional, intent(out) :: residual
    integer, optional, intent(out) :: diagnostic(2) ! iterations, failure reason (1 input/2 solve/3 search/4 limit)
    real(mn_dp), optional, intent(inout) :: dual_hint(:) ! numerical initial guess only; updated on success
    ! Optional bounded-memory replay coefficients: alpha(1:nm-1), logit max,
    ! normalization. Replay uses the SAME input moments as its target. This
    ! is a stage-local numerical cache, not checkpoint/persistent physics.
    real(mn_dp), optional, intent(inout) :: stream_cache(:)
    integer,optional,intent(in) :: max_iterations
    real(mn_dp) :: alpha(b%nm-1),trial(b%nm-1),target(b%nm-1),mean(b%nm-1),nextmean(b%nm-1)
    real(mn_dp) :: gradient(b%nm-1),hessian(b%nm-1,b%nm-1),step(b%nm-1)
    real(mn_dp) :: centered(b%nm-1,b%nq)
    real(mn_dp) :: eta(b%nq),trial_eta(b%nq),eta_step(b%nq)
    real(mn_dp) :: probability(b%nq),nextprob(b%nq),objective,nextobjective,scale,err
    real(mn_dp) :: radius,decrement,predicted,ratio,merit,nextmerit,merit_slope,change,top,partition
    integer :: it,q,ls,status,iteration_limit
    ierr=mn_bad_input;intensity=0
    ! A realizable M5 tail from live AMR/defrag needs 687 Newton steps
    ! at the unchanged 2e-12 residual. Budget exhaustion is not infeasibility.
    iteration_limit=1024
    if(present(max_iterations))iteration_limit=max_iterations
    if(iteration_limit<1)return
    if(present(diagnostic))diagnostic=[0,1]
    if(present(residual))residual=huge(1d0)
    if(b%nm==0.or.size(u)/=b%nm.or.size(intensity)/=b%nq)return
    if(present(stream_cache))then
       if(size(stream_cache)/=b%nm+1)return
    endif
    if(any(.not.ieee_is_finite(u)))return
    if(u(1)<0)return
    if(u(1)==0)then
       if(any(u/=0))return
       ierr=mn_ok
       if(present(diagnostic))diagnostic=0
       if(present(residual))residual=0
       if(present(stream_cache))stream_cache=0
       return
    endif
    target=u(2:)/u(1)
    if(any(.not.ieee_is_finite(target)))return
    if(any(abs(target)>maxval(abs(b%harmonic(2:,:)),dim=2)+1d-12))return
    alpha=0;ierr=mn_no_closure;radius=1d0
    if(present(dual_hint))then
       if(size(dual_hint)==b%nm-1)then
          if(all(ieee_is_finite(dual_hint)))alpha=dual_hint
       endif
    endif
    if(present(diagnostic))diagnostic=[0,4]
    ! Work with target-centered logits and update them incrementally. Re-forming
    ! alpha.Y and subtracting alpha.target at every trial loses significant
    ! digits when large dual coefficients cancel near the realizable boundary.
    ! In exact arithmetic eta_q=alpha.(Y_q-target): identical dual problem.
    do q=1,b%nq
       eta(q)=dot_product(alpha,b%harmonic(2:,q)-target)
    enddo
    ! Front tails can approach the realizable boundary (large dual variables)
    ! even with tiny physical density. Keep the moment tolerance unchanged;
    ! allow the safeguarded solve to finish rather than mixing in isotropy.
    do it=1,iteration_limit
       if(present(diagnostic))diagnostic(1)=it
       call dual(eta,probability,mean,objective)
       gradient=mean-target;err=maxval(abs(gradient))
       if(err<=2d-12.and.present(stream_cache))then
          ! Incrementally updated logits need not equal a replayed dot product
          ! to the last bit near the cone boundary. Validate the replay itself
          ! before exporting; if needed continue Newton at the unchanged tol.
          do q=1,b%nq
             eta(q)=dot_product(alpha,b%harmonic(2:,q)-target)
          enddo
          call dual(eta,probability,mean,objective)
          gradient=mean-target;err=maxval(abs(gradient))
       endif
       if(present(residual))residual=err
       if(err<=2d-12)then
          intensity=u(1)*(probability/sum(probability))/b%weight
          if(present(stream_cache))then
             top=maxval(eta);partition=sum(b%weight*exp(eta-top))
             intensity=u(1)*(exp(eta-top)/partition)
          endif
          if(any(.not.ieee_is_finite(intensity)))return
          if(present(stream_cache))stream_cache=[alpha,top,partition]
          if(present(dual_hint))then
             if(size(dual_hint)==b%nm-1)dual_hint=alpha
          endif
          if(present(diagnostic))diagnostic(2)=0
          ierr=mn_ok;return
       endif
       ! Weighted covariance as a matrix contraction, enabling optimized
       ! BLAS builds without changing the moments or adding a library API.
       do q=1,b%nq
          centered(:,q)=(b%harmonic(2:,q)-mean)*sqrt(probability(q))
       enddo
       hessian=matmul(centered,transpose(centered))
       call whitened_newton(hessian,gradient,step,status)
       if(status/=0)then
          if(present(diagnostic))diagnostic(2)=2
          return
       endif
       ! Radial trust-region step in the covariance metric, not an arbitrary
       ! bound on individual dual coefficients. The physical target is fixed.
       decrement=dot_product(gradient,step)
       if(.not.ieee_is_finite(decrement).or.decrement<=0)then
          if(present(diagnostic))diagnostic(2)=2
          return
       endif
       scale=min(1d0,radius/sqrt(decrement))
       merit=dot_product(gradient,gradient)
       merit_slope=dot_product(gradient,matmul(hessian,step))
       do q=1,b%nq
          eta_step(q)=-dot_product(step,b%harmonic(2:,q)-target)
       enddo
       do ls=1,48
          trial=alpha-scale*step
          trial_eta=eta+scale*eta_step
          call dual(trial_eta,nextprob,nextmean,nextobjective)
          if(maxval(abs(nextmean-target))<=2d-12)exit
          change=local_change(probability,eta,scale*eta_step,objective,nextobjective)
          if(change<0.and.change<=-1d-4*scale*decrement)exit
          if(err<1d-7)then
             ! Supplement stable objective descent with residual descent in
             ! its roundoff-limited regime. Do not require gradient norm to
             ! decrease on EVERY Newton step of a nonquadratic objective.
             nextmerit=sum((nextmean-target)**2)
             if(nextmerit<merit.and.nextmerit<=merit- &
                  1d-4*scale*max(0d0,merit_slope))exit
          endif
          scale=scale*.5d0
       enddo
       if(ls>48)then
          if(present(diagnostic))diagnostic(2)=3
          return
       endif
       predicted=scale*decrement-.5d0*scale**2*dot_product(step,matmul(hessian,step))
       ! Objective differences are cancellation dominated near the solution;
       ! retain the radius there and rely on the residual acceptance above.
       if(err>=1d-7.and.predicted>0)then
          change=local_change(probability,eta,scale*eta_step,objective,nextobjective)
          ratio=-change/predicted
          if(ratio<.25d0)radius=max(1d-8,.25d0*radius)
          if(ratio>.75d0.and.scale*sqrt(decrement)>.8d0*radius)radius=min(16d0,2d0*radius)
       endif
       alpha=trial;eta=trial_eta
    enddo
  contains
    real(mn_dp) function local_change(p,logits,increment,value,nextvalue) result(result)
      ! Exact identity: delta Phi = log(sum p*exp(delta alpha.(Y-target))).
      ! expm1/log1p corrections and compensated summation preserve small
      ! changes. Log-space evaluation for large increments revives nodes
      ! whose current probability has underflowed; never multiply p in place.
      real(mn_dp), intent(in) :: p(:),logits(:),increment(:),value,nextvalue
      real(mn_dp) :: s,compensation,term,y,t,v,postlog
      integer :: node
      s=0;compensation=0
      do node=1,b%nq
         if(increment(node)>50d0)then
            postlog=log(b%weight(node))+logits(node)-value+increment(node)
            if(postlog>50d0)then
               result=huge(1d0);return ! necessarily an uphill trial
            endif
            term=exp(postlog)-p(node)
         else if(increment(node)<-50d0)then
            ! expm1 rounds to -1 here. Avoid log(exp(x)) corrections on
            ! subnormal exp(x), where relative rounding is no longer small.
            term=-p(node)
         else
            v=exp(increment(node))
            if(v==1)then
               term=p(node)*increment(node)
            else if(v==0)then
               term=-p(node)
            else
               term=p(node)*(v-1)*(increment(node)/log(v))
            endif
         endif
         y=term-compensation;t=s+y;compensation=(t-s)-y;s=t
      enddo
      if(s<=-.5d0)then
         result=nextvalue-value ! large decrease: no cancellation risk
      else
         v=1+s
         if(v==1)then
            result=s
         else
            result=log(v)*(s/(v-1))
         endif
      endif
    end function
    subroutine dual(logits,p,avg,value)
      real(mn_dp), intent(in) :: logits(:)
      real(mn_dp), intent(out) :: p(:),avg(:),value
      real(mn_dp) :: top,partition
      top=maxval(logits)
      p=b%weight*exp(logits-top);partition=sum(p);p=p/partition
      avg=matmul(b%harmonic(2:,:),p)
      value=top+log(partition)
    end subroutine
  end subroutine

  subroutine whitened_newton(matrix,rhs,x,ierr)
    ! H=L L^T changes the local dual coordinates to a unit covariance basis.
    ! Triangular solves avoid explicitly forming an ill-conditioned inverse.
    ! A diagonal shift may stabilize the STEP near rank deficiency. It is
    ! never added to the dual objective or the requested moment constraints.
    real(mn_dp), intent(in) :: matrix(:,:),rhs(:)
    real(mn_dp), intent(out) :: x(:)
    integer, intent(out) :: ierr
    real(mn_dp) :: l(size(rhs),size(rhs)),v(size(rhs)),pivot,damping,matrix_scale
    integer :: i,j,n,attempt
    logical :: positive
    n=size(rhs);l=0;x=0;ierr=mn_no_closure
    matrix_scale=max(1d0,maxval(abs(matrix)));damping=0
    do attempt=1,10
       l=0;positive=.true.
       do i=1,n
          do j=1,i
             pivot=matrix(i,j)-dot_product(l(i,1:j-1),l(j,1:j-1))
             if(i==j)then
                pivot=pivot+damping
                if(.not.ieee_is_finite(pivot).or.pivot<=epsilon(1d0)*matrix_scale)then
                   positive=.false.;exit
                endif
                l(i,j)=sqrt(pivot)
             else
                l(i,j)=pivot/l(j,j)
             endif
          enddo
          if(.not.positive)exit
       enddo
       if(positive)exit
       damping=max(1d-15*matrix_scale,10*damping)
    enddo
    if(.not.positive)return
    do i=1,n
       v(i)=(rhs(i)-dot_product(l(i,1:i-1),v(1:i-1)))/l(i,i)
    enddo
    do i=n,1,-1
       x(i)=(v(i)-dot_product(l(i+1:n,i),x(i+1:n)))/l(i,i)
    enddo
    if(any(.not.ieee_is_finite(x)))return
    ierr=mn_ok
  end subroutine

  subroutine mn_flux(b,left,right,normal,flux,ierr)
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(in) :: left(:),right(:),normal(3)
    real(mn_dp), intent(out) :: flux(:)
    integer, intent(out) :: ierr
    real(mn_dp) :: il(b%nq),ir(b%nq),mu(b%nq)
    ierr=mn_bad_input
    if(size(flux)/=b%nm.or.any(.not.ieee_is_finite(normal)))return
    if(abs(sum(normal**2)-1d0)>1d-12)return
    call mn_reconstruct(b,left,il,ierr)
    if(ierr/=mn_ok)return
    call mn_reconstruct(b,right,ir,ierr)
    if(ierr/=mn_ok)return
    mu=matmul(normal,b%direction)
    flux=matmul(b%harmonic,b%weight*(max(mu,0d0)*il+min(mu,0d0)*ir))
  end subroutine

  subroutine mn_advance(b,u,volume,left,right,normal,area,cdt,boundary_integral,projection_integral,ierr, &
       centres,face_centres,face_integral,stream_angles)
    ! One face, ONE flux, opposite updates to its two cells, including unequal
    ! AMR volumes. right=0 denotes vacuum, not an MPI ghost. MPI callers must
    ! supply exchanged states and use a shared face ownership convention.
    ! SSPRK2 supplies a local time predictor/corrector without storing history.
    ! Optional Cartesian centres enable minmod MUSCL reconstruction; slopes
    ! are angular work arrays, never persistent radiation state. A failed
    ! realizability/closure/CFL check leaves u AND the receipt unchanged.
    ! Numerical projection is explicit: sum(V*(u_new-u_old)) + boundary
    ! equals projection_integral (up to roundoff), NOT zero for every moment.
    ! The caller must retain this receipt, not deposit it as physical exchange.
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(inout) :: u(:,:),boundary_integral(:),projection_integral(:)
    real(mn_dp), intent(in) :: volume(:),normal(:,:),area(:),cdt
    integer, intent(in) :: left(:),right(:)
    integer, intent(out) :: ierr
    real(mn_dp), optional, intent(in) :: centres(:,:),face_centres(:,:)
    ! Time/area-integrated moments per oriented face for AMR reflux registers.
    real(mn_dp), optional, intent(inout) :: face_integral(:,:)
    ! Stream one angular ordinate across the cell/face list using compact
    ! closure coefficients. No nq-by-ncell scratch. The old workspace path
    ! stays available for comparison; no change to its default arithmetic.
    logical, optional, intent(in) :: stream_angles
    real(mn_dp) :: stage(size(u,1),size(u,2)),next(size(u,1),size(u,2)),first(b%nm),second(b%nm)
    real(mn_dp) :: outgoing(size(u,2)),check(b%nq)
    real(mn_dp) :: face_first(b%nm,size(left)),face_second(b%nm,size(left))
    real(mn_dp) :: base(size(u,1),size(u,2)),project_first(b%nm),project_second(b%nm)
    integer :: f,i,nc,nf
    ierr=mn_bad_input;nc=size(u,2);nf=size(left)
    ! MLP needs complete fixed-grid vertex connectivity, provided by the
    ! benchmark transport. Do not silently use a 1-D fallback on this face list.
    if(b%spatial_limiter<1.or.b%spatial_limiter>3)return
    if(b%nm==0.or.size(u,1)/=b%nm.or.size(boundary_integral)/=b%nm.or.size(projection_integral)/=b%nm)return
    if(size(volume)/=nc.or.size(right)/=nf.or.size(area)/=nf.or.any(shape(normal)/=[3,nf]))return
    if(present(face_integral))then
       if(any(shape(face_integral)/=[b%nm,nf]))return
    endif
    if(any(left<1).or.any(left>nc).or.any(right<0).or.any(right>nc).or.any(left==right))return
    if(any(.not.ieee_is_finite(u)).or.any(.not.ieee_is_finite(volume)))return
    if(any(.not.ieee_is_finite(normal)).or.any(.not.ieee_is_finite(area)))return
    if(.not.ieee_is_finite(cdt))return
    if(any(volume<=0).or.any(area<=0).or.cdt<0)return
    if(present(centres).neqv.present(face_centres))return
    if(present(centres).and.present(face_centres))then
       if(any(shape(centres)/=[3,nc]).or.any(shape(face_centres)/=[3,nf]))return
       if(any(.not.ieee_is_finite(centres)).or.any(.not.ieee_is_finite(face_centres)))return
    endif
    outgoing=0
    do f=1,nf
       if(abs(sum(normal(:,f)**2)-1d0)>1d-12)return
       outgoing(left(f))=outgoing(left(f))+area(f)
       if(right(f)>0)outgoing(right(f))=outgoing(right(f))+area(f)
    enddo
    ierr=mn_bad_step
    ! Conservative sufficient bound, including all incident face areas.
    if(any(cdt*outgoing/volume>.5d0+8*epsilon(1d0)))return
    if(cdt==0)then
       boundary_integral=0;projection_integral=0;ierr=mn_ok
       if(present(face_integral))face_integral=0
       return
    endif
    call dispatch_euler(u,stage,first,project_first,face_first,.true.,ierr)
    if(ierr/=mn_ok)return
    call dispatch_euler(stage,next,second,project_second,face_second,.false.,ierr)
    if(ierr/=mn_ok)return
    next=.5d0*base+.5d0*next
    do i=1,nc
       call mn_reconstruct(b,next(:,i),check,ierr)
       if(ierr/=mn_ok)return
    enddo
    u=next;boundary_integral=.5d0*(first+second);ierr=mn_ok
    projection_integral=project_first+.5d0*project_second
    if(present(face_integral))face_integral=.5d0*(face_first+face_second)
  contains
    subroutine dispatch_euler(input,output,escape,projection,face_receipt,save_base,status)
      real(mn_dp),intent(in) :: input(:,:)
      real(mn_dp),intent(out) :: output(:,:),escape(:),projection(:),face_receipt(:,:)
      logical,intent(in) :: save_base
      integer,intent(out) :: status
      if(present(stream_angles))then
         if(stream_angles)then
            call euler_stream(input,output,escape,projection,face_receipt,save_base,status)
            return
         endif
      endif
      call euler(input,output,escape,projection,face_receipt,save_base,status)
    end subroutine

    subroutine euler_stream(input,output,escape,projection,face_receipt,save_base,status)
      real(mn_dp),intent(in) :: input(:,:)
      real(mn_dp),intent(out) :: output(:,:),escape(:),projection(:),face_receipt(:,:)
      logical,intent(in) :: save_base
      integer,intent(out) :: status
      real(mn_dp) :: cache(b%nm+1,nc),angular(b%nq),projected(b%nm),theta(nc)
      real(mn_dp) :: intensity(nc),slope(3,nc),lo(3,nc),hi(3,nc),transported(nc)
      real(mn_dp) :: d,distance,change,il,ir,mu,angular_flux,flux(b%nm)
      logical :: seenlo(3,nc),seenhi(3,nc),muscl
      integer :: cell,face,a,l,r,q,side,pass
      muscl=present(centres).and.present(face_centres)
      projection=0;output=0;face_receipt=0;escape=0;theta=1
      do cell=1,nc
         call mn_reconstruct(b,input(:,cell),angular,status,stream_cache=cache(:,cell))
         if(status/=mn_ok)return
         call mn_project(b,angular,projected,status)
         if(status/=mn_ok)return
         projection=projection+volume(cell)*(projected-input(:,cell))
         if(save_base)base(:,cell)=projected
      enddo
      ! Pass one computes the SAME single theta per cell over every ordinate
      ! and incident face. Pass two replays intensities/slopes and transports.
      ! A separate theta for each ordinate would change the old method.
      do pass=1,2
         if(pass==1.and..not.muscl)cycle
         do q=1,b%nq
            do cell=1,nc
               intensity(cell)=0
               if(input(1,cell)>0)intensity(cell)=input(1,cell)* &
                    (exp(dot_product(cache(1:b%nm-1,cell),b%harmonic(2:,q)- &
                    input(2:,cell)/input(1,cell))-cache(b%nm,cell))/cache(b%nm+1,cell))
            enddo
            slope=0;lo=0;hi=0;seenlo=.false.;seenhi=.false.
            if(muscl)then
               do face=1,nf
                  a=maxloc(abs(normal(:,face)),dim=1)
                  if(abs(abs(normal(a,face))-1d0)>1d-12)then
                     status=mn_bad_input;return
                  endif
                  l=left(face);r=right(face)
                  if(r==0)cycle
                  distance=centres(a,r)-centres(a,l)
                  if(abs(distance)<tiny(1d0))then
                     status=mn_bad_input;return
                  endif
                  do side=1,2
                     cell=l
                     if(side==2)cell=r
                     d=(intensity(r)-intensity(l))/distance
                     if((distance>0).eqv.(side==1))then
                        if(seenhi(a,cell))d=minmod(hi(a,cell),d)
                        hi(a,cell)=d;seenhi(a,cell)=.true.
                     else
                        if(seenlo(a,cell))d=minmod(lo(a,cell),d)
                        lo(a,cell)=d;seenlo(a,cell)=.true.
                     endif
                  enddo
               enddo
               do cell=1,nc
                  do a=1,3
                     if(seenlo(a,cell).and.seenhi(a,cell)) &
                          slope(a,cell)=mn_limited_slope(lo(a,cell),hi(a,cell),b%spatial_limiter)
                  enddo
               enddo
               if(pass==1)then
                  do face=1,nf
                     do side=1,2
                        cell=left(face)
                        if(side==2)cell=right(face)
                        if(cell==0)cycle
                        change=dot_product(slope(:,cell),face_centres(:,face)-centres(:,cell))
                        if(change<0)theta(cell)=min(theta(cell),.5d0*intensity(cell)/(-change))
                     enddo
                  enddo
               else
                  do cell=1,nc
                     slope(:,cell)=theta(cell)*slope(:,cell)
                  enddo
               endif
            endif
            if(pass==1)cycle
            transported=intensity
            do face=1,nf
               l=left(face);r=right(face);il=intensity(l);ir=0
               if(r>0)ir=intensity(r)
               if(muscl)then
                  il=il+dot_product(slope(:,l),face_centres(:,face)-centres(:,l))
                  if(r>0)ir=ir+dot_product(slope(:,r),face_centres(:,face)-centres(:,r))
               endif
               mu=dot_product(normal(:,face),b%direction(:,q))
               angular_flux=cdt*area(face)*(max(mu,0d0)*il+min(mu,0d0)*ir)
               flux=b%harmonic(:,q)*(b%weight(q)*angular_flux)
               face_receipt(:,face)=face_receipt(:,face)+flux
               transported(l)=transported(l)-angular_flux/volume(l)
               if(r>0)then
                  transported(r)=transported(r)+angular_flux/volume(r)
               else
                  escape=escape+flux
               endif
            enddo
            if(any(.not.ieee_is_finite(transported)).or.any(transported<0))then
               status=mn_bad_step;return
            endif
            do cell=1,nc
               output(:,cell)=output(:,cell)+b%harmonic(:,q)*(b%weight(q)*transported(cell))
            enddo
         enddo
      enddo
      status=mn_ok
    end subroutine

    subroutine euler(input,output,escape,projection,face_receipt,save_base,status)
      real(mn_dp), intent(in) :: input(:,:)
      real(mn_dp), intent(out) :: output(:,:),escape(:),projection(:),face_receipt(:,:)
      logical, intent(in) :: save_base
      integer, intent(out) :: status
      real(mn_dp) :: intensity(b%nq,nc),slope(b%nq,3,nc),lo(b%nq,3,nc),hi(b%nq,3,nc)
      real(mn_dp) :: il(b%nq),ir(b%nq),mu(b%nq),flux(b%nm),d(b%nq),offset(3),theta(nc),change(b%nq)
      real(mn_dp) :: distance
      real(mn_dp) :: transported(b%nq,nc),angular_flux(b%nq),projected(b%nm)
      logical :: seenlo(3,nc),seenhi(3,nc)
      integer :: cell,face,a,l,r,q,side
      projection=0
      do cell=1,nc
         call mn_reconstruct(b,input(:,cell),intensity(:,cell),status)
         if(status/=mn_ok)return
         call mn_project(b,intensity(:,cell),projected,status)
         if(status/=mn_ok)return
         projection=projection+volume(cell)*(projected-input(:,cell))
         if(save_base)base(:,cell)=projected
      enddo
      transported=intensity
      slope=0;lo=0;hi=0;seenlo=.false.;seenhi=.false.
      if(present(centres).and.present(face_centres))then
         ! Cartesian minmod. Coarse/fine multiplicities take the most
         ! restrictive one-sided slope; unsupported oblique normals reject.
         do face=1,nf
            a=maxloc(abs(normal(:,face)),dim=1)
            if(abs(abs(normal(a,face))-1d0)>1d-12)then
               status=mn_bad_input;return
            endif
            l=left(face);r=right(face)
            if(r==0)cycle
            distance=centres(a,r)-centres(a,l)
            if(abs(distance)<tiny(1d0))then
               status=mn_bad_input;return
            endif
            d=(intensity(:,r)-intensity(:,l))/distance
            do side=1,2
               cell=l
               if(side==2)cell=r
               if((distance>0).eqv.(side==1))then
                  if(seenhi(a,cell))d=minmod(hi(:,a,cell),d)
                  hi(:,a,cell)=d;seenhi(a,cell)=.true.
               else
                  if(seenlo(a,cell))d=minmod(lo(:,a,cell),d)
                  lo(:,a,cell)=d;seenlo(a,cell)=.true.
               endif
               d=(intensity(:,r)-intensity(:,l))/distance
            enddo
         enddo
         do cell=1,nc
            do a=1,3
               if(seenlo(a,cell).and.seenhi(a,cell)) &
                    slope(:,a,cell)=mn_limited_slope(lo(:,a,cell),hi(:,a,cell),b%spatial_limiter)
            enddo
         enddo
         theta=1
         do face=1,nf
            do side=1,2
               cell=left(face)
               if(side==2)cell=right(face)
               if(cell==0)cycle
               offset=face_centres(:,face)-centres(:,cell)
               change=matmul(slope(:,:,cell),offset)
               do q=1,b%nq
                  if(change(q)<0)theta(cell)=min(theta(cell),.5d0*intensity(q,cell)/(-change(q)))
               enddo
            enddo
         enddo
         do cell=1,nc
            slope(:,:,cell)=theta(cell)*slope(:,:,cell)
         enddo
      endif
      escape=0
      do face=1,nf
         l=left(face);r=right(face);il=intensity(:,l);ir=0
         if(r>0)ir=intensity(:,r)
         if(present(centres).and.present(face_centres))then
            il=il+matmul(slope(:,:,l),face_centres(:,face)-centres(:,l))
            if(r>0)ir=ir+matmul(slope(:,:,r),face_centres(:,face)-centres(:,r))
         endif
         mu=matmul(normal(:,face),b%direction)
         angular_flux=cdt*area(face)*(max(mu,0d0)*il+min(mu,0d0)*ir)
         flux=matmul(b%harmonic,b%weight*angular_flux)
         face_receipt(:,face)=flux
         transported(:,l)=transported(:,l)-angular_flux/volume(l)
         if(r>0)then
            transported(:,r)=transported(:,r)+angular_flux/volume(r)
         else
            escape=escape+flux
         endif
      enddo
      do cell=1,nc
         call mn_project(b,transported(:,cell),output(:,cell),status)
         if(status/=mn_ok)return
      enddo
      status=mn_ok
    end subroutine
  end subroutine

  elemental real(mn_dp) function minmod(a,b) result(c)
    real(mn_dp), intent(in) :: a,b
    c=0
    if((a>0.and.b>0).or.(a<0.and.b<0))c=sign(min(abs(a),abs(b)),a)
  end function

  subroutine mn_material_exchange(b,u,tau_abs,tau_scatter,emitted,receipt,ierr)
    ! Isotropic absorption/emission and isotropic elastic scattering.
    ! emission is the already time-integrated injection, not a rate. This
    ! operator is source-split; no claim of a stiff coupled material solver.
    ! receipt = OLD radiation + emitted monopole - NEW radiation. Components
    ! 2:4 encode material momentum times c with the known l=1 normalization.
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(inout) :: u(:),receipt(:)
    real(mn_dp), intent(in) :: tau_abs,tau_scatter,emitted
    integer, intent(out) :: ierr
    real(mn_dp) :: candidate(b%nm),check(b%nq),balance(b%nm)
    ierr=mn_bad_input
    if(size(u)/=b%nm.or.size(receipt)/=b%nm)return
    if(.not.all(ieee_is_finite([tau_abs,tau_scatter,emitted])))return
    if(min(tau_abs,tau_scatter,emitted)<0)return
    call mn_reconstruct(b,u,check,ierr)
    if(ierr/=mn_ok)return
    candidate=u*exp(-tau_abs);candidate(2:)=candidate(2:)*exp(-tau_scatter)
    candidate(1)=candidate(1)+emitted
    call mn_reconstruct(b,candidate,check,ierr)
    if(ierr/=mn_ok)return
    balance=u-candidate;balance(1)=balance(1)+emitted
    u=candidate;receipt=balance
  end subroutine
end module

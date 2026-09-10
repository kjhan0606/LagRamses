! Standalone, local moving elastic-scattering receiver; no RAMSES state.
! N(d,g) is direction-INTEGRATED photon number density; E(d,g)
! is independent direction-integrated LAB energy density.
! rho(b), p(:,b) are absolute phase mass and momentum densities.
! tau(b,g) is an already integrated, frozen, nonnegative scattering depth.
! A caller may form tau with its interaction clock (including reduced c),
! but recoil, Doppler factors and work ALWAYS use the required physicalc.
! Unit contract: E, rho, p and physicalc MUST use one consistent unit system
! (e.g. erg/cm^3, g/cm^3, g/cm^2/s, cm/s). E in eV/cm^3 requires matching
! mass and momentum units, or conversion to cgs before calling. No energy
! conversion is inferred. physicalc is NEVER an interaction/reduced speed;
! arbitrary code units prevent this routine from numerically identifying a
! caller's reduced-c substitution. snrt_scatter_c is physical c in cgs.
!
! Physics: isotropic coherent events in each phase's comoving frame,
! grey within each group, no group-boundary transfer, thermal recoil,
! absorption, grain heat or gas reservoir. Nonrelativistic (|v|/c <= .01).
! The retained Doppler/aberration factors are NOT a full relativistic model:
! nu0 = gamma*nu*q, q=1-beta.n; event probability a ~ w/q^2. See
! Tominaga et al. 2015, sec II.2, https://arxiv.org/abs/1507.05141 .
! The common gamma cancels in event probability and coherent energy ratios.
!
! Discretization (derived here): the continuous frozen-beta collision
! generator for both X=N and X=q*E is dX/dtau=-q*X+a*sum(q*X).
! Backward Euler gives removal f=tau*q/(1+tau*q), and redistributes with
! A ~ a/(1+tau*q). Thus Xnew=(1-f)*X+A*sum(f*X). This is positive,
! conserves sum(X), and preserves the moving isotropic equilibrium
! N~w/q^3, E~w/q^4 at ANY tau. Using exponential removal followed by a
! single redistribution with a instead would give the wrong finite-step
! equilibrium X~a/(1-exp(-tau*q)). BE is deliberately first order.
!
! Solve beta_mid=(pold+pnew)/(2*rho*c) over ALL groups of one phase.
! sum(q*(Enew-Eold))=0 implies delta Erad=c*beta_mid.delta prad;
! pnew=pold-delta prad then closes Newtonian kinetic+rad energy exactly
! at the discrete level (up to the nonlinear solve and roundoff).
! Phases split in input order (first order). No signed work is thermalized.
! ALL outputs publish only after every phase succeeds. Caller must not
! alias N,E,p or work with each other or with any input.
module snrt_moving_scatter
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(real64),parameter,public :: snrt_scatter_c=2.99792458d10
  real(real64),parameter,public :: snrt_scatter_beta_limit=1d-2
  integer,parameter,public :: moving_scatter_ok=0,moving_scatter_config=1
  integer,parameter,public :: moving_scatter_state=2,moving_scatter_range=3
  integer,parameter,public :: moving_scatter_convergence=4
  integer,parameter,public :: moving_scatter_allocation=5
  public :: snrt_moving_scatter_cell,snrt_moving_scatter_energy_cell
contains
  subroutine snrt_moving_scatter_energy_cell(energy,momenta,rho,tau,direction,weight,physicalc,work,ierr,max_iterations)
    ! IR/grey energy-only API. The internal auxiliary scalar merely reuses
    ! the positive collision map; it is NOT a physical photon population.
    ! No photon state is supplied, published or claimed conserved by this API.
    real(real64),intent(inout) :: energy(:,:),momenta(:,:),work(:)
    real(real64),intent(in) :: rho(:),tau(:,:),direction(:,:),weight(:),physicalc
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: max_iterations
    real(real64),allocatable :: auxiliary(:,:)
    integer :: stat
    allocate(auxiliary(size(energy,1),size(energy,2)),stat=stat)
    ierr=moving_scatter_allocation
    if(stat/=0)return
    auxiliary=energy
    call snrt_moving_scatter_cell(auxiliary,energy,momenta,rho,tau,direction,weight,physicalc,work,ierr,max_iterations)
  end subroutine snrt_moving_scatter_energy_cell

  subroutine snrt_moving_scatter_cell(number,energy,momenta,rho,tau,direction,weight,physicalc,work,ierr,max_iterations)
    real(real64),intent(inout) :: number(:,:),energy(:,:),momenta(:,:),work(:)
    real(real64),intent(in) :: rho(:),tau(:,:),direction(:,:),weight(:),physicalc
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: max_iterations
    real(real64),allocatable :: nn(:,:),ee(:,:),pp(:,:),ww(:),w(:),en(:,:),out_n(:,:),out_e(:,:)
    real(real64) :: beta0(3),beta(3),res(3),shift(3),jac(3,3),probe(3),rp(3),candidate(3),rc(3)
    real(real64) :: impulse(3),pnext(3),actual_beta(3),delta_e,esc,inertia,ratio,h,alpha,rnorm,oldnorm
    real(real64) :: residual_tol,energy_tol,phase_work,scale,bound,check_p(3)
    integer :: nd,ng,nb,b,d,g,k,it,ls,niter,stat
    logical :: solved,accepted,linear_ok
    ierr=moving_scatter_config
    nd=size(number,1);ng=size(number,2);nb=size(rho)
    if(min(nd,ng,nb)<1)return
    if(any(shape(energy)/=[nd,ng]).or.any(shape(momenta)/=[3,nb]))return
    if(any(shape(tau)/=[nb,ng]).or.any(shape(direction)/=[3,nd]))return
    if(size(weight)/=nd.or.size(work)/=nb)return
    if(.not.ieee_is_finite(physicalc))return
    if(physicalc<=sqrt(tiny(1d0)).or.physicalc>=sqrt(huge(1d0)/32))return
    niter=64
    if(present(max_iterations))niter=max_iterations
    if(niter<1.or.niter>256)return
    if(any(.not.ieee_is_finite(direction)).or.any(.not.ieee_is_finite(weight)))return
    if(any(weight<=0).or.any(abs(direction)>1d0+1d-12))return
    allocate(w(nd),stat=stat)
    ierr=moving_scatter_allocation
    if(stat/=0)return
    ierr=moving_scatter_config
    w=weight/maxval(weight);w=w/sum(w)
    do d=1,nd
       if(abs(sum(direction(:,d)**2)-1d0)>1d-12)return
    enddo
    ! A stationary isotropic field must not acquire quadrature momentum.
    if(maxval(abs(matmul(direction,w)))>1d-12)return
    ierr=moving_scatter_state
    if(any(.not.ieee_is_finite(number)).or.any(.not.ieee_is_finite(energy)))return
    if(any(.not.ieee_is_finite(momenta)).or.any(.not.ieee_is_finite(rho)))return
    if(any(.not.ieee_is_finite(tau)))return
    if(any(number<0).or.any(energy<0).or.any(rho<0).or.any(tau<0))return
    if(any((number==0).neqv.(energy==0)))return ! no energy without photons or zero-energy photons
    bound=huge(1d0)/(32d0*real(nd,real64)*real(ng,real64)*real(nb,real64))
    if(maxval(number)>bound.or.maxval(energy)>bound)return
    do b=1,nb
       if(rho(b)==0)then
          ierr=moving_scatter_state
          if(any(momenta(:,b)/=0).or.any(tau(b,:)/=0))return
       else
          ierr=moving_scatter_range
          if(rho(b)>huge(1d0)/(32*max(physicalc,physicalc**2,1d0)))return
          scale=rho(b)*physicalc
          if(scale<=tiny(1d0))return
          if(any(abs(momenta(:,b))>snrt_scatter_beta_limit*scale))return
          if(norm2(momenta(:,b)/scale)>snrt_scatter_beta_limit)return
       endif
    enddo
    allocate(nn(nd,ng),ee(nd,ng),pp(3,nb),ww(nb),en(nd,ng),out_n(nd,ng),out_e(nd,ng),stat=stat)
    ierr=moving_scatter_allocation
    if(stat/=0)return
    nn=number;ee=energy;pp=momenta;ww=0
    do b=1,nb
       if(rho(b)==0.or.all(tau(b,:)==0).or.all(ee==0))cycle
       esc=maxval(ee);inertia=rho(b)*physicalc**2
       ! Guard arithmetic before the nonlinear solve, including FP-trap builds.
       ! Ratios beyond inverse roundoff cannot resolve grain backreaction.
       ierr=moving_scatter_range
       if(inertia<=tiny(1d0).or.inertia<esc*epsilon(1d0))return
       if(physicalc<1)then
          if(esc>(huge(1d0)*physicalc)/(8d0*real(nd,real64)*real(ng,real64)))return
       endif
       ratio=esc/inertia;en=ee/esc
       beta0=pp(:,b)/(rho(b)*physicalc);beta=beta0
       residual_tol=32*epsilon(1d0)*snrt_scatter_beta_limit
       solved=.false.
       do it=1,niter
          call residual(beta,res)
          if(any(.not.ieee_is_finite(res)))return
          rnorm=norm2(res)
          if(rnorm<=residual_tol)then
             solved=.true.;exit
          endif
          ! A 3x3 finite-difference Newton solve; line search keeps q positive
          ! and the midpoint inside the stated nonrelativistic domain.
          do k=1,3
             h=1d-6
             if(beta(k)>0)h=-h
             probe=beta;probe(k)=probe(k)+h
             call residual(probe,rp)
             jac(:,k)=(rp-res)/h
          enddo
          rp=-res
          call solve_three(jac,rp,shift,linear_ok)
          if(.not.linear_ok)exit
          alpha=1;accepted=.false.;oldnorm=rnorm
          do ls=1,32
             candidate=beta+alpha*shift
             if(norm2(candidate)<=snrt_scatter_beta_limit)then
                call residual(candidate,rc)
                if(norm2(rc)<oldnorm.or.norm2(rc)<=residual_tol)then
                   beta=candidate;accepted=.true.;exit
                endif
             endif
             alpha=.5d0*alpha
          enddo
          if(.not.accepted)exit
       enddo
       ierr=moving_scatter_convergence
       if(.not.solved)return
       call collision(beta,nn,en,tau(b,:),direction,w,out_n,out_e,check_p,delta_e)
       impulse=-(esc/physicalc)*check_p
       pnext=pp(:,b)+impulse
       actual_beta=.5d0*(pp(:,b)/(rho(b)*physicalc)+pnext/(rho(b)*physicalc))
       ierr=moving_scatter_range
       if(norm2(pnext/(rho(b)*physicalc))>snrt_scatter_beta_limit)return
       if(any(.not.ieee_is_finite(out_n)).or.any(.not.ieee_is_finite(out_e)))return
       if(any(out_n<0).or.any(out_e<0))return
       if(any((out_n==0).neqv.(out_e==0)))return
       ierr=moving_scatter_convergence
       if(norm2(actual_beta-beta)>2*residual_tol)return
       phase_work=dot_product(actual_beta,impulse)*physicalc
       energy_tol=1024*epsilon(1d0)*sum(en)
       ! Test both stable event differences and the actual stored candidate.
       ! No repair, clipping, energy rescaling, or invisible heat reservoir.
       if(abs(delta_e+phase_work/esc)>energy_tol)return
       if(abs(sum(out_e-en)+phase_work/esc)>energy_tol)return
       do g=1,ng
          if(abs(sum(out_n(:,g))-sum(nn(:,g)))>1024*epsilon(1d0)*sum(nn(:,g)))return
       enddo
       if(maxval(out_e)>huge(1d0)/max(esc,1d0))return
       nn=out_n;ee=out_e*esc;pp(:,b)=pnext;ww(b)=phase_work
    enddo
    number=nn;energy=ee;momenta=pp;work=ww;ierr=moving_scatter_ok
  contains
    subroutine residual(at,value)
      real(real64),intent(in) :: at(3)
      real(real64),intent(out) :: value(3)
      real(real64) :: dp_rad(3),de
      call collision(at,nn,en,tau(b,:),direction,w,out_n,out_e,dp_rad,de)
      value=(at-beta0)+(.5d0*ratio)*dp_rad
    end subroutine
  end subroutine snrt_moving_scatter_cell

  subroutine collision(beta,number,energy,tau,direction,w,next_n,next_e,prad,de)
    real(real64),intent(in) :: beta(3),number(:,:),energy(:,:),tau(:),direction(:,:),w(:)
    real(real64),intent(out) :: next_n(:,:),next_e(:,:),prad(3),de
    real(real64) :: q(size(w)),a(size(w)),prob(size(w)),f(size(w)),keep(size(w))
    real(real64) :: removed(size(w)),emitted(size(w)),change(size(w)),denom(size(w)),captured,comoving
    integer :: g
    q=1-matmul(beta,direction);a=w/q**2;a=a/sum(a)
    prad=0;de=0
    do g=1,size(tau)
       if(tau(g)==0)then
          next_n(:,g)=number(:,g);next_e(:,g)=energy(:,g);cycle
       endif
       if(tau(g)>1)then
          denom=q+1/tau(g);keep=(1/tau(g))/denom;f=q/denom
          prob=a/denom
       else
          denom=1+tau(g)*q;keep=1/denom;f=tau(g)*q/denom
          prob=a/denom
       endif
       prob=prob/sum(prob)
       captured=sum(f*number(:,g));removed=f*energy(:,g)
       comoving=sum(q*removed);emitted=prob*comoving/q
       next_n(:,g)=keep*number(:,g)+prob*captured
       next_e(:,g)=keep*energy(:,g)+emitted
       change=emitted-removed
       prad=prad+matmul(direction,change);de=de+sum(change)
    enddo
  end subroutine collision

  subroutine solve_three(matrix,rhs,x,ok)
    real(real64),intent(in) :: matrix(3,3),rhs(3)
    real(real64),intent(out) :: x(3)
    logical,intent(out) :: ok
    real(real64) :: a(3,4),row(4),factor,scale
    integer :: k,j,pivot
    ok=.false.;x=0;a(:,1:3)=matrix;a(:,4)=rhs
    if(any(.not.ieee_is_finite(a)))return
    scale=maxval(abs(matrix))
    do k=1,3
       pivot=k-1+maxloc(abs(a(k:3,k)),dim=1)
       if(abs(a(pivot,k))<=64*epsilon(1d0)*scale)return
       row=a(k,:);a(k,:)=a(pivot,:);a(pivot,:)=row
       do j=k+1,3
          factor=a(j,k)/a(k,k);a(j,k:4)=a(j,k:4)-factor*a(k,k:4)
       enddo
    enddo
    do k=3,1,-1
       x(k)=(a(k,4)-sum(a(k,k+1:3)*x(k+1:3)))/a(k,k)
    enddo
    ok=all(ieee_is_finite(x))
  end subroutine solve_three
end module snrt_moving_scatter

! Focused native test, independent of RAMSES and its Makefile.
! Compile in a separate build directory, for example:
! gfortran -std=f2008 -O0 -g -Wall -Wextra -Wno-compare-reals \
!   -fcheck=all -ffpe-trap=invalid,zero,overflow \
!   /absolute/source/snrt_moving_scatter.f90 \
!   /absolute/source/snrt_moving_scatter_smoke.f90 -o moving_scatter_smoke
program snrt_moving_scatter_smoke
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use snrt_moving_scatter
  implicit none
  integer,parameter :: nd=6,ng=3,nb=3
  real(real64) :: n(nd,ng),e(nd,ng),p(3,nb),rho(nb),tau(nb,ng),dirs(3,nd),w(nd),work(nb)
  real(real64) :: n0(nd,ng),e0(nd,ng),p0(3,nb),work0(nb),t0(nb,ng)
  real(real64) :: na(nd,ng),ea(nd,ng),pa(3,nb),wa(nb),q(nd),beta(3),c,expected,ceff,dt,opacity
  real(real64) :: rotated(3,nd),pr(3,nb),rel,max_number,max_energy,max_momentum,ev,err(3),exact(nd)
  integer :: ierr,b,g,k,j,count,steps,s
  c=snrt_scatter_c;count=0;max_number=0;max_energy=0;max_momentum=0
  dirs=0
  do k=1,3
     dirs(k,2*k-1)=1;dirs(k,2*k)=-1
  enddo
  w=1d0/nd

  call reset()
  tau=0
  call save()
  call advance()
  call require(ierr==0,'zero optical depth accepted')
  call require(all(n==n0).and.all(e==e0).and.all(p==p0).and.all(work==0),'zero-depth identity')

  ! Comoving-isotropic angular number/energy are DIFFERENT Doppler powers.
  ! Test transparent, moderate and very thick groups, all three phases.
  call reset()
  beta=[.003d0,-.002d0,.001d0];q=1-matmul(beta,dirs)
  do g=1,ng
     n(:,g)=1d9*g*w/q**3;e(:,g)=1d-3*g*w/q**4
  enddo
  do b=1,nb
     p(:,b)=rho(b)*c*beta;tau(b,:)=[1d-8,1d0,1d6]
  enddo
  call save();call advance()
  call require(ierr==0,'moving isotropic equilibrium accepted')
  call require(maxval(abs(n-n0))/maxval(n0)<2d-13,'moving number equilibrium')
  call require(maxval(abs(e-e0))/maxval(e0)<2d-13,'moving energy equilibrium')
  call require(maxval(abs(p-p0))*c/sum(e0)<2d-12,'moving equilibrium no force')
  call closure('moving equilibrium')

  ! Initially stationary finite-mass scatterer: beam produces recoil AND
  ! loses precisely the finite quadratic Newtonian kinetic energy.
  call reset();call beam()
  tau=0;tau(1,:)=[.2d0,1d0,4d0]
  call save();call advance()
  call require(ierr==0,'stationary beam finite recoil accepted')
  call require(p(1,1)>0.and.work(1)>0.and.sum(e)<sum(e0),'beam work debited from radiation')
  expected=sum(p(:,1)**2)/(2*rho(1))
  call require(abs(work(1)-expected)<1d-11*expected,'finite quadratic recoil work')
  call require(abs(sum(e0-e)-expected)<1d-8*expected,'beam independent energy debit')
  call require(abs(sum(e(:,1))/sum(n(:,1))-sum(e0(:,1))/sum(n0(:,1)))>0, &
       'photon energies evolve independently of photon number')
  call closure('stationary beam')

  ! Massive-grain limit gives an independent BE beam impulse reference.
  ! Effective c changes exposure tau only; equal exposure must give the
  ! same recoil. A mistaken 1/c_eff recoil would be wrong by a factor 100.
  call reset();call beam();rho=1d-10
  opacity=1d-10;dt=1d0;ceff=c
  tau=0;tau(1,:)=opacity*ceff*dt
  call save();call advance()
  call require(ierr==0,'physical-c beam reference')
  expected=sum(e0)*t0(1,1)/(1+t0(1,1))/c
  call require(abs(p(1,1)/expected-1)<1d-10,'recoil uses physical c')
  na=n;ea=e;pa=p;wa=work
  n=n0;e=e0;p=p0;work=work0;ceff=.01d0*c;dt=100d0
  tau(1,:)=opacity*ceff*dt
  call advance()
  call require(ierr==0,'equal reduced-c exposure accepted')
  call require(maxval(abs(p-pa))*c/sum(e0)<2d-13,'c_eff absent from recoil')
  call require(maxval(abs(e-ea))/sum(e0)<2d-13,'equal exposure energy')
  call closure('physical c')

  ! First-order collision-time convergence to the independent static
  ! exponential solution. Finite-step BE is not mislabeled as exponential.
  do j=1,3
     call reset();call beam();rho=1d-5;steps=16*2**(j-1)
     tau=0;tau(1,:)=2d0/steps
     call save()
     do s=1,steps
        call advance()
        call require(ierr==0,'time refinement accepted')
     enddo
     exact=sum(e0(:,1))*w+(e0(:,1)-sum(e0(:,1))*w)*exp(-2d0)
     err(j)=sum(abs(e(:,1)-exact))/sum(e0(:,1))
  enddo
  call require(err(1)/err(2)>1.9d0.and.err(1)/err(2)<2.1d0,'first-order time refinement 16/32')
  call require(err(2)/err(3)>1.9d0.and.err(2)/err(3)<2.1d0,'first-order time refinement 32/64')
  write(*,'(A,3ES14.5)')'BE time refinement L1 errors (16,32,64): ',err

  ! Lab-isotropic radiation is NOT equilibrium for a moving grain.
  ! Thin-limit force is -(4/3)*tau*Erad*beta/c on an isotropic quadrature.
  call reset();p(1,1)=rho(1)*c*1d-4;tau=0;tau(1,:)=1d-5
  call save();call advance()
  call require(ierr==0,'moving grain in lab-isotropic radiation')
  call require(p(1,1)<p0(1,1).and.work(1)<0.and.sum(e)>sum(e0),'drag gives energy back to photons')
  expected=-(4d0/3)*1d-5*sum(e0)*1d-4/c
  call require(abs((p(1,1)-p0(1,1))/expected-1)<5d-5,'mixed-frame thin drag coefficient')
  call closure('radiative drag')

  ! Each material sees the updated radiation but keeps its own absolute p.
  call reset()
  n(1,:)=2*n(1,:);e(1,:)=3*e(1,:);e(3,:)=2*e(3,:)
  do b=1,nb
     p(:,b)=rho(b)*c*[real(b,real64)*3d-4,-2d-4,1d-4]
     tau(b,:)=[.1d0,.7d0,2d0]*b
  enddo
  call save();call advance()
  call require(ierr==0,'multiple phases accepted');call closure('multiple phases')
  na=n;ea=e;pa=p;wa=work
  n=n0;e=e0;p=p0;work=work0
  do b=1,nb
     call snrt_moving_scatter_cell(n,e,p(:,b:b),rho(b:b),tau(b:b,:),dirs,w,c,work(b:b),ierr)
     call require(ierr==0,'individual split phase accepted')
  enddo
  call require(all(n==na).and.all(e==ea).and.all(p==pa).and.all(work==wa),'ordered phase split reference')

  ! IR API publishes energy/momenta only; its internal scalar is not an IR
  ! photon ledger. Energy and recoil cannot depend on the primary N array.
  e=e0;p=p0;work=work0
  call snrt_moving_scatter_energy_cell(e,p,rho,tau,dirs,w,c,work,ierr)
  call require(ierr==0,'energy-only IR entry accepted')
  call require(all(e==ea).and.all(p==pa).and.all(work==wa),'energy-only/paired recoil parity')

  ! eV energy units require the SAME energy-unit conversion of rho and p.
  ! This is a unit covariance test, not permission to mix eV with cgs mass.
  ev=1.602176634d-12;n=n0;e=e0/ev;p=p0/ev;work=work0
  call snrt_moving_scatter_cell(n,e,p,rho/ev,tau,dirs,w,c,work,ierr)
  call require(ierr==0,'consistent eV units accepted')
  call require(maxval(abs(e*ev-ea))/sum(e0)<3d-12,'eV energy covariance')
  call require(maxval(abs(p*ev-pa))*c/sum(e0)<3d-12,'eV momentum covariance')

  ! c=1 physical code units: rho_code=rho*c_cgs^2, p_code=p*c_cgs.
  n=n0;e=e0;p=p0*c;work=work0
  call snrt_moving_scatter_cell(n,e,p,rho*c*c,tau,dirs,w,1d0,work,ierr)
  call require(ierr==0,'physical c=1 units accepted')
  call require(maxval(abs(e-ea))/sum(e0)<3d-12,'c=1 energy covariance')
  call require(maxval(abs(p/c-pa))*c/sum(e0)<3d-12,'c=1 momentum covariance')

  ! Rotate the full problem: no distinguished force axis in the solver.
  rotated(1,:)=dirs(2,:);rotated(2,:)=dirs(3,:);rotated(3,:)=dirs(1,:)
  pr(1,:)=p0(2,:);pr(2,:)=p0(3,:);pr(3,:)=p0(1,:)
  n=n0;e=e0;work=work0
  call snrt_moving_scatter_cell(n,e,pr,rho,tau,rotated,w,c,work,ierr)
  call require(ierr==0,'rotated state accepted')
  call require(maxval(abs(pr(1,:)-pa(2,:)))*c/sum(e0)<1d-10,'rotation x')
  call require(maxval(abs(pr(2,:)-pa(3,:)))*c/sum(e0)<1d-10,'rotation y')
  call require(maxval(abs(pr(3,:)-pa(1,:)))*c/sum(e0)<1d-10,'rotation z')

  ! Static isotropic limit, with extreme finite tau, has no energy/force.
  call reset();tau=huge(1d0)
  call save();call advance()
  call require(ierr==0,'extreme finite optical depth accepted')
  call require(maxval(abs(p))*c/sum(e0)<1d-12,'static isotropic no recoil')
  call closure('thick isotropic')

  ! Empty radiation and an absent (zero density, zero opacity) phase.
  call reset();n=0;e=0;rho(2)=0;tau(2,:)=0
  call save();call advance()
  call require(ierr==0.and.all(work==0),'empty radiation/absent phase')
  call require(all(n==n0).and.all(e==e0).and.all(p==p0),'empty identity')

  ! Deterministic directional/group/phase mixtures spanning thin to thick.
  do j=1,24
     call reset()
     do g=1,ng
        do k=1,nd
           n(k,g)=1d8*(1+mod(k*g+j,11))
           e(k,g)=1d-5*(1+mod(3*k+g*j,17))
        enddo
     enddo
     do b=1,nb
        p(:,b)=rho(b)*c*[sin(real(j+b,real64)),cos(real(j+2*b,real64)),.3d0]*1d-3
        tau(b,:)=10d0**[real(mod(j+b,7)-4,real64),-2d0,2d0]
     enddo
     call save();call advance()
     call require(ierr==0,'deterministic mixture accepted');call closure('mixture sweep')
  enddo

  call reset();call beam();tau=1
  call save()
  call snrt_moving_scatter_cell(n,e,p,rho,tau,dirs,w,c,work,ierr,max_iterations=1)
  call require(ierr==moving_scatter_convergence,'nonlinear exhaustion rejected');call unchanged('solve rollback')

  ! First phase demonstrably succeeds; a much lighter last phase cannot
  ! receive the remaining beam within |v|/c <= .01. No earlier phase leaks.
  call reset();call beam();tau=0;tau(1,:)=.1d0;tau(3,:)=1d0;rho(3)=1d-24
  call save()
  call snrt_moving_scatter_cell(n,e,p(:,1:1),rho(1:1),tau(1:1,:),dirs,w,c,work(1:1),ierr)
  call require(ierr==0.and.p(1,1)>0,'first phase of late rejection succeeds')
  n=n0;e=e0;p=p0;work=work0
  call advance()
  call require(ierr/=0,'late nonrelativistic rejection');call unchanged('whole-cell rollback')

  call reset();tau(2,2)=-1;call save();call advance()
  call require(ierr/=0,'negative depth rejected');call unchanged('negative-depth rollback')
  call reset();tau(2,2)=ieee_value(0d0,ieee_quiet_nan);call save();call advance()
  call require(ierr/=0,'NaN depth rejected');call unchanged('NaN rollback')
  call reset();rho(3)=0;call save();call advance()
  call require(ierr/=0,'opacity without phase mass rejected');call unchanged('absent-mass rollback')
  call reset();p(1,2)=rho(2)*c*.02d0;call save();call advance()
  call require(ierr==moving_scatter_range,'super-admissible velocity rejected');call unchanged('velocity rollback')
  call reset();n(1,1)=0;call save();call advance()
  call require(ierr/=0,'energy without photons rejected');call unchanged('number-energy support rollback')
  call reset();e(1,1)=-1;call save();call advance()
  call require(ierr/=0,'negative energy rejected');call unchanged('negative-energy rollback')
  call reset();call save()
  call snrt_moving_scatter_cell(n,e,p,rho,tau,dirs,w(:5),c,work,ierr)
  call require(ierr==moving_scatter_config,'shape rejected');call unchanged('shape rollback')
  call reset();call save();rotated=dirs;rotated(:,1)=2*rotated(:,1)
  call snrt_moving_scatter_cell(n,e,p,rho,tau,rotated,w,c,work,ierr)
  call require(ierr==moving_scatter_config,'nonunit direction rejected');call unchanged('quadrature rollback')
  call reset();call save()
  call snrt_moving_scatter_energy_cell(e,p,rho,tau,dirs,w,0d0,work,ierr)
  call require(ierr==moving_scatter_config,'invalid physical c rejected');call unchanged('IR rollback')
  call reset();call save();tau(2,1)=-1
  call snrt_moving_scatter_energy_cell(e,p,rho,tau,dirs,w,c,work,ierr)
  call require(ierr/=0,'IR negative depth rejected');call unchanged('IR state rollback')

  write(*,'(A,I0)')'MOVING_SCATTER_SMOKE_PASS assertions=',count
  write(*,'(A,3ES14.5)')'max relative number / energy / momentum residuals: ',max_number,max_energy,max_momentum
contains
  subroutine reset()
    integer :: g
    do g=1,ng
       n(:,g)=1d8*g*w;e(:,g)=1d-3*g*w
    enddo
    rho=[1d-20,2d-20,3d-20];p=0;tau=.3d0;work=-777
  end subroutine
  subroutine beam()
    integer :: g
    n=0;e=0
    do g=1,ng
       n(1,g)=1d8*g;e(1,g)=1d-3*g
    enddo
  end subroutine
  subroutine save()
    n0=n;e0=e;p0=p;work0=work;t0=tau
  end subroutine
  subroutine advance()
    call snrt_moving_scatter_cell(n,e,p,rho,tau,dirs,w,c,work,ierr)
  end subroutine
  subroutine require(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if(.not.condition)then
       write(*,'(A,A,A,I0)')'FAIL: ',label,' ierr=',ierr
       error stop 1
    endif
    count=count+1
  end subroutine
  subroutine unchanged(label)
    character(len=*),intent(in) :: label
    call require(all(n==n0).and.all(e==e0).and.all(p==p0).and.all(work==work0),label)
  end subroutine
  subroutine closure(label)
    character(len=*),intent(in) :: label
    real(real64) :: dk,total,pscale,dp_rad(3),dp_grain(3),part
    integer :: b,g
    dk=0;total=sum(e0)
    do b=1,nb
       if(rho(b)==0)cycle
       ! Independent kinetic difference; stable for moving phases.
       part=dot_product(.5d0*(p(:,b)+p0(:,b))/rho(b),p(:,b)-p0(:,b))
       dk=dk+part
       call require(abs(part-work(b))<3d-12*max(total,abs(part)),label//' phase work')
    enddo
    rel=abs(sum(e-e0)+dk)/total
    max_energy=max(max_energy,rel)
    call require(rel<3d-12,label//' radiation+kinetic energy')
    dp_rad=matmul(dirs,sum(e-e0,dim=2))/c;dp_grain=sum(p-p0,dim=2)
    pscale=total/c;rel=maxval(abs(dp_rad+dp_grain))/pscale
    max_momentum=max(max_momentum,rel)
    call require(rel<3d-12,label//' momentum')
    do g=1,ng
       rel=abs(sum(n(:,g))-sum(n0(:,g)))/sum(n0(:,g))
       max_number=max(max_number,rel)
       call require(rel<3d-12,label//' group photon number')
    enddo
    call require(all(n>=0).and.all(e>=0),label//' positivity')
  end subroutine
end program snrt_moving_scatter_smoke

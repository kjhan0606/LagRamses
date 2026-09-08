! H/He non-equilibrium thermal receiver for the dust comparison.
! Rosdahl+2013 Appendix E, case B; existing native recombination fits.
! WSS09 metal terms remain an explicitly separate CIE approximation.
module snrt_atomic_cooling
  use amr_parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_thermochemistry, only: snrt_alpha_hydrogen_case_b,snrt_alpha_helium_ii_case_b, &
       snrt_alpha_helium_ii_radiative_case_b,snrt_alpha_helium_iii_case_b
  use dust_element_cooling, only: wss09_metals_rate
  implicit none
  private
  real(dp),parameter,public::atomic_mh=1.66d-24,atomic_kb=1.380649d-16
  real(dp),parameter::ev=1.602176634d-12,tolerance=1d-3,fraction_floor=1d-5
  real(dp),parameter::atom(11)=[1d0,4d0,12d0,14d0,16d0,20d0,24d0,28d0,32d0,40d0,56d0]
  public::atomic_temperature,atomic_heat_capacity,atomic_rates,atomic_advance,atomic_identity
contains
  function atomic_identity() result(v)
    real(dp)::v(18)
    ! Version changes when fits/algorithm/normalization change. The WSS09
    ! numerical tables are independently bound by their existing identity.
    v=[1d0,atomic_mh,atomic_kb,ev,tolerance,fraction_floor,atom,1d0]
  end function

  real(dp) function atomic_heat_capacity(rho,gas_x,x,gamma) result(cv)
    real(dp),intent(in)::rho,gas_x(11),x(3),gamma
    real(dp)::ne
    ne=gas_x(1)*x(1)+gas_x(2)/4*(x(2)+2*x(3))
    cv=rho/atomic_mh*(sum(gas_x/atom)+ne)*atomic_kb/(gamma-1)
  end function

  real(dp) function atomic_temperature(rho,gas_x,x,gamma,energy) result(t)
    real(dp),intent(in)::rho,gas_x(11),x(3),gamma,energy
    t=energy/atomic_heat_capacity(rho,gas_x,x,gamma)
  end function

  subroutine atomic_rates(t,beta,alpha,recomb,excitation,brems,dielectronic,ierr)
    real(dp),intent(in)::t
    real(dp),intent(out)::beta(3),alpha(3),recomb(3),excitation(2),brems,dielectronic
    integer,intent(out)::ierr
    real(dp)::f,lambda
    ierr=1;beta=0;alpha=0;recomb=0;excitation=0;brems=0;dielectronic=0
    if(.not.ieee_is_finite(t).or.t<1d0.or.t>1d9)return
    f=1+sqrt(t/1d5)
    beta=[5.85d-11,2.38d-11,5.68d-12]*sqrt(t)/f*exp(-[157809.1d0,285335.4d0,631515d0]/t)
    alpha=[snrt_alpha_hydrogen_case_b(t),snrt_alpha_helium_ii_case_b(t),snrt_alpha_helium_iii_case_b(t)]
    lambda=315614d0/t
    recomb(1)=3.435d-30*t*lambda**1.97/(1+(lambda/2.25d0)**.376d0)**3.72
    recomb(2)=atomic_kb*t*snrt_alpha_helium_ii_radiative_case_b(t)
    lambda=1263030d0/t
    recomb(3)=27.48d-30*t*lambda**1.97/(1+(lambda/2.25d0)**.376d0)**3.72
    excitation=[7.5d-19*exp(-118348d0/t),5.54d-17*t**(-.397d0)*exp(-473638d0/t)]/f
    brems=1.42d-27*sqrt(t) ! Rosdahl+2013 Appendix E (no extra Gaunt factor).
    dielectronic=1.24d-13*t**(-1.5d0)*exp(-470000d0/t)*(1+.3d0*exp(-94000d0/t))
    if(.not.all(ieee_is_finite([beta,alpha,recomb,excitation,brems,dielectronic])))return
    ierr=0
  end subroutine

  subroutine atomic_trial(rho,gas_x,gamma,aexp,h,e,x,next_e,next_x,ierr)
    real(dp),intent(in)::rho,gas_x(11),gamma,aexp,h,e,x(3)
    real(dp),intent(out)::next_e,next_x(3)
    integer,intent(out)::ierr
    real(dp)::t,nh,nhe,ne,q,beta(3),alpha(3),rec(3),exc(2),bre,die,neutral(3),ion(3)
    real(dp)::a,b,y,rate,metal,tbg
    next_e=e;next_x=x;ierr=1
    t=atomic_temperature(rho,gas_x,x,gamma,e)
    call atomic_rates(t,beta,alpha,rec,exc,bre,die,ierr)
    if(ierr/=0)return
    nh=rho/atomic_mh*gas_x(1);nhe=rho/atomic_mh*gas_x(2)/4
    ne=nh*x(1)+nhe*(x(2)+2*x(3));q=ne*h
    ! Positive backward-Euler tridiagonal H/He system, coefficients and ne
    ! frozen for this substep. Step doubling controls the resulting error.
    next_x(1)=(x(1)+q*beta(1))/(1+q*(beta(1)+alpha(1)))
    a=1+q*beta(2);b=1+q*alpha(3)
    y=(x(2)+(1-x(2)-x(3))*q*beta(2)/a+x(3)*q*alpha(3)/b)/ &
         (1+q*alpha(2)/a+q*beta(3)/b)
    next_x(2)=y;next_x(3)=(x(3)+q*beta(3)*y)/b
    ierr=1
    if(any(.not.ieee_is_finite(next_x)).or.any(next_x<0).or.next_x(1)>1)return
    if(next_x(2)+next_x(3)>1+32*epsilon(1d0))return
    neutral=[nh*(1-next_x(1)),nhe*max(0d0,1-next_x(2)-next_x(3)),nhe*next_x(2)]
    ion=[nh*next_x(1),nhe*next_x(2),nhe*next_x(3)]
    ! Retain this substep's frozen ne in both reaction counts and cooling.
    ! Ionization cooling uses the SAME beta and thresholds as the chemistry.
    ! This charges the kinetic-to-ionization-potential transfer exactly once.
    rate=ne*(sum(beta*[13.60d0,24.59d0,54.42d0]*ev*neutral)+sum(rec*ion)+ &
         exc(1)*neutral(1)+exc(2)*neutral(3)+bre*(ion(1)+ion(2)+4*ion(3))+die*ion(2))
    tbg=2.727d0/aexp
    rate=rate+ne*1.017d-37*tbg**4*(t-tbg)
    call wss09_metals_rate(t,gas_x,metal,ierr)
    if(ierr/=0)return
    rate=rate+nh*nh*metal
    next_e=e-h*rate;ierr=1
    if(.not.ieee_is_finite(next_e).or.next_e<=0)return
    if(abs(next_e-e)>.1d0*e)return
    t=atomic_temperature(rho,gas_x,next_x,gamma,next_e)
    call wss09_metals_rate(t,gas_x,metal,ierr)
    if(t<1d0.or.t>1d9)ierr=1
  end subroutine

  subroutine atomic_advance(rho,gas_x,gamma,aexp,dt,energy,x,next_energy,next_x,net_loss,ierr)
    real(dp),intent(in)::rho,gas_x(11),gamma,aexp,dt,energy,x(3)
    real(dp),intent(out)::next_energy,next_x(3),net_loss
    integer,intent(out)::ierr
    real(dp)::left,h,e,e1,e2,ef,y(3),y1(3),y2(3),yf(3),err,denom(5),s0(5),sf(5),s2(5)
    integer::k,status
    ierr=1;next_energy=energy;next_x=x;net_loss=0
    if(.not.all(ieee_is_finite([rho,gas_x,gamma,aexp,dt,energy,x])))return
    if(rho<=0.or.gamma<=1.or.aexp<=0.or.dt<0.or.energy<=0.or.any(gas_x<0))return
    if(gas_x(1)<=0.or.sum(gas_x)>1+128*epsilon(1d0))return
    if(any(x<0).or.x(1)>1.or.x(2)+x(3)>1)return
    e=energy;y=x;left=dt;h=dt
    do k=1,20000
       if(left<=0)exit
       h=min(h,left)
       if(h<=0.or.(left-h==left))return
       call atomic_trial(rho,gas_x,gamma,aexp,h,e,y,ef,yf,status)
       if(status==0)call atomic_trial(rho,gas_x,gamma,aexp,h/2,e,y,e1,y1,status)
       if(status==0)call atomic_trial(rho,gas_x,gamma,aexp,h/2,e1,y1,e2,y2,status)
       if(status/=0)then
          h=h/2;cycle
       endif
       s0=[1-y(1),y(1),1-y(2)-y(3),y(2),y(3)]
       sf=[1-yf(1),yf(1),1-yf(2)-yf(3),yf(2),yf(3)]
       s2=[1-y2(1),y2(1),1-y2(2)-y2(3),y2(2),y2(3)]
       denom=max(abs(s0),abs(s2),fraction_floor)
       err=max(abs(e2-ef)/max(e,e2),maxval(abs(s2-sf)/denom))
       if(err>tolerance.or.maxval(abs(s2-s0)/denom)>.1d0)then
          h=h/2;cycle
       endif
       e=e2;y=y2;left=left-h
       if(err<tolerance/4)h=h*2
    enddo
    if(left>0)return
    next_energy=e;next_x=y;net_loss=energy-e;ierr=0
  end subroutine
end module

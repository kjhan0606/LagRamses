! Wiersma, Schaye & Smith (2009), author-supplied z_collis.txt.
! Opt-in, low-density CIE comparison, NOT local-radiation/NEQ metal cooling.
module dust_element_cooling
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter,public :: cie_n=352
  ! Match the active RAMSES units() convention exactly (mH is its mass unit).
  real(real64),parameter :: cie_mh=1.6600000d-24,cie_kb=1.3806200d-16
  real(real64),parameter :: atom(11)=[1d0,4d0,12d0,14d0,16d0,20d0,24d0,28d0,32d0,40d0,56d0]
  real(real64),parameter :: solar(9)=[2.46d-4,8.51d-5,4.90d-4,1d-4,3.47d-5,3.47d-5,1.86d-5,2.29d-6,2.82d-5]
  real(real64),parameter :: he(7)=[.0786528d0,.0830474d0,.0875605d0,.0921968d0,.0969616d0,.10186d0,.106898d0]
  real(real64),save :: cie_data(25,cie_n)
  include 'dust_wss09_cie_data.inc'
  public :: wss09_curve,wss09_step,wss09_rate,wss09_identity,wss09_metals_rate
contains
  subroutine wss09_metals_rate(temperature,gas_x,rate,ierr)
    ! Equation (3), metal terms only. No H/He CIE cooling/electron heat
    ! capacity enters this receiver. Metal ions still assume CIE; NOT NEQ
    ! metal cooling or an invented equation-(4) solar-electron correction.
    real(real64),intent(in)::temperature,gas_x(11)
    real(real64),intent(out)::rate
    integer,intent(out)::ierr
    real(real64)::w,weights(9)
    integer::j
    ierr=1;rate=0
    if(.not.all(ieee_is_finite(gas_x)).or.any(gas_x<0))return
    if(gas_x(1)<=0.or.sum(gas_x)>1+128*epsilon(1d0))return
    if(.not.ieee_is_finite(temperature).or.temperature<=0)return
    if(all(gas_x(3:11)==0))then
       ierr=0;return
    endif
    if(temperature<cie_data(1,1).or.temperature>cie_data(1,cie_n))return
    j=interval(cie_data(1,:),temperature)
    w=(temperature-cie_data(1,j))/(cie_data(1,j+1)-cie_data(1,j))
    weights=gas_x(3:11)/(atom(3:11)*gas_x(1)*solar)
    rate=sum(weights*((1-w)*cie_data(17:25,j)+w*cie_data(17:25,j+1)))
    if(.not.ieee_is_finite(rate))return
    ierr=0
  end subroutine

  subroutine wss09_identity(first,second)
    ! Split below the HDF5 per-attribute size limit; bind all numerical data
    ! and conversion coefficients, not merely a path or user-provided label.
    real(real64),intent(out)::first(15*cie_n+30),second(10*cie_n)
    first=[1d0,cie_mh,cie_kb,atom,solar,he,reshape(cie_data(1:15,:),[15*cie_n])]
    second=reshape(cie_data(16:25,:),[10*cie_n])
  end subroutine

  subroutine wss09_curve(gas_x,tau,rate,ierr)
    ! x = gas-phase elemental mass / total mixture mass, order H..Fe.
    ! Tabulated H/He net rates and electrons are interpolated in nHe/nH.
    ! Metals use WSS09 equation (3); the ASCII table does not give the solar
    ! electron reference needed by equation (4). Do not invent that ratio.
    real(real64),intent(in)::gas_x(11)
    real(real64),intent(out)::tau(cie_n),rate(cie_n)
    integer,intent(out)::ierr
    real(real64)::y,w,ne(cie_n),particles(cie_n),weight(9)
    integer::j,k
    ierr=1;tau=0;rate=0
    if(.not.all(ieee_is_finite(gas_x)).or.any(gas_x<0))return
    if(gas_x(1)<=0.or.sum(gas_x)>1+128*epsilon(1d0))return
    y=gas_x(2)/(4*gas_x(1))
    if(y<he(1).or.y>he(7))return ! No silent He extrapolation.
    j=1
    do while(j<6)
       if(y<he(j+1))exit
       j=j+1
    enddo
    w=(y-he(j))/(he(j+1)-he(j))
    rate=(1-w)*cie_data(2*j,:)+w*cie_data(2*j+2,:)
    ne=(1-w)*cie_data(2*j+1,:)+w*cie_data(2*j+3,:)
    weight=gas_x(3:11)/(atom(3:11)*gas_x(1)*solar)
    do k=1,9
       rate=rate+weight(k)*cie_data(16+k,:)
    enddo
    ! CIE H/He electrons + one particle per gas-phase metal nucleus.
    ! Trace-metal approximation: no invented metal ion fractions/electrons.
    particles=gas_x(1)*(1+y+ne)+sum(gas_x(3:11)/atom(3:11))
    tau=cie_data(1,:)*particles ! T/mu = T * particles per mixture mass unit.
    if(.not.all(ieee_is_finite(rate)).or..not.all(ieee_is_finite(tau)))return
    if(any(tau(2:)<=tau(:cie_n-1)))return
    ierr=0
  end subroutine

  subroutine wss09_rate(t2,gas_x,rate,temperature,ierr)
    real(real64),intent(in)::t2,gas_x(11)
    real(real64),intent(out)::rate,temperature
    integer,intent(out)::ierr
    real(real64)::tau(cie_n),rates(cie_n),w
    integer::j
    rate=0;temperature=0
    call wss09_curve(gas_x,tau,rates,ierr)
    if(ierr/=0)return
    ierr=2
    if(.not.ieee_is_finite(t2).or.t2<tau(1).or.t2>tau(cie_n))return
    j=interval(tau,t2);w=(t2-tau(j))/(tau(j+1)-tau(j))
    rate=(1-w)*rates(j)+w*rates(j+1)
    temperature=(1-w)*cie_data(1,j)+w*cie_data(1,j+1)
    ierr=0
  end subroutine

  subroutine wss09_step(rho,gas_x,t2,dt,gamma,delta,ierr)
    ! Exact integration of a piecewise-linear NET rate in thermal energy
    ! (T/mu), frozen density and elemental abundances. Signed rates and
    ! equilibrium zeros are retained; no scalar-Z curve is added a second time.
    real(real64),intent(in)::rho,gas_x(11),t2,dt,gamma
    real(real64),intent(out)::delta
    integer,intent(out)::ierr
    real(real64)::tau(cie_n),rates(cie_n),x,left,k,r,s,boundary,rb,edge_time,z,fac,w
    integer::j,n
    delta=0;ierr=1
    if(.not.all(ieee_is_finite([rho,t2,dt,gamma])))return
    if(rho<=0.or.dt<0.or.gamma<=1)return
    call wss09_curve(gas_x,tau,rates,ierr)
    if(ierr/=0)return
    ierr=2
    if(t2<tau(1).or.t2>tau(cie_n))return
    k=(gamma-1)*gas_x(1)**2*rho/(cie_mh*cie_kb)
    x=t2;left=dt
    do n=1,cie_n+1
       if(left<=0)exit
       j=interval(tau,x)
       w=(x-tau(j))/(tau(j+1)-tau(j))
       r=(1-w)*rates(j)+w*rates(j+1)
       if(r==0)then
          left=0
          exit
       endif
       if(r>0.and.x==tau(j))then
          if(j==1)return ! Would leave the published domain: fail, no clamp.
          j=j-1
       else if(r<0.and.x==tau(cie_n))then
          return
       endif
       s=(rates(j+1)-rates(j))/(tau(j+1)-tau(j))
       boundary=tau(j);rb=rates(j)
       if(r<0)then
          boundary=tau(j+1);rb=rates(j+1)
       endif
       edge_time=huge(1d0)
       if(rb*r>0)then
          if(abs(s*(boundary-x))<1d-8*abs(r))then
             edge_time=abs(boundary-x)/(k*abs(r))
          else
             edge_time=-log(rb/r)/(k*s)
          endif
       endif
       if(left>=edge_time)then
          x=boundary;left=left-edge_time
       else
          z=-k*s*left
          if(abs(z)<1d-5)then
             fac=1+z/2+z*z/6+z*z*z/24
          else
             fac=(exp(z)-1)/z
          endif
          x=x-k*r*left*fac
          left=0
       endif
    enddo
    if(.not.ieee_is_finite(x).or.left>0.or.x<tau(1).or.x>tau(cie_n))return
    delta=x-t2;ierr=0
  end subroutine

  integer function interval(tau,x) result(j)
    real(real64),intent(in)::tau(cie_n),x
    integer::lo,hi,m
    lo=1;hi=cie_n
    do while(hi-lo>1)
       m=(hi+lo)/2
       if(x>=tau(m))then
          lo=m
       else
          hi=m
       endif
    enddo
    j=lo
  end function
end module

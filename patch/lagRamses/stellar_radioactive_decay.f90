! Named approximate abundance comparison; no RAMSES globals or dispatch.
module stellar_radioactive_decay
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter,public :: radioactive_nparent=2
  character(len=*),parameter,public :: radioactive_model='lc18_al26_fe60_transparent_v1'
  ! NUBASE2020 ground states, pinned by the companion converter. LC18's
  ! single Al26 label is the long-lived inventory; no invented isomer split.
  real(real64),parameter,public :: radioactive_year_s=31557600d0
  real(real64),parameter,public :: radioactive_half_life_s(2)=[717d3,2.62d6]*radioactive_year_s
  real(real64),parameter :: lambda(2)=log(2d0)/radioactive_half_life_s
  public :: radioactive_decay,radioactive_release,radioactive_gas_decay,radioactive_element_shift
contains
  subroutine kernel(dt,survive,lost,mean_survive,mean_lost)
    real(real64),intent(in)::dt
    real(real64),intent(out)::survive(2),lost(2),mean_survive(2),mean_lost(2)
    real(real64)::x
    integer::j
    do j=1,2
       x=lambda(j)*dt
       survive(j)=exp(-x)
       if(x<1d-4)then
          lost(j)=x*(1-x*(.5d0-x*(1d0/6-x*(1d0/24-x/120))))
          mean_lost(j)=x*(.5d0-x*(1d0/6-x*(1d0/24-x*(1d0/120-x/720))))
          mean_survive(j)=1-mean_lost(j)
       else
          lost(j)=1-survive(j)
          mean_survive(j)=lost(j)/x
          mean_lost(j)=1-mean_survive(j)
       endif
    enddo
  end subroutine

  subroutine radioactive_decay(parent,dt_s,surviving,decayed,ierr)
    ! Old GAS only: call once per local physical timestep, not per source.
    real(real64),intent(in)::parent(2),dt_s
    real(real64),intent(out)::surviving(2),decayed(2)
    integer,intent(out)::ierr
    real(real64)::s(2),d(2),ms(2),md(2)
    surviving=parent;decayed=0;ierr=1
    if(.not.all(ieee_is_finite([parent,dt_s])))return
    if(min(minval(parent),dt_s)<0)return
    call kernel(dt_s,s,d,ms,md)
    surviving=parent*s;decayed=parent*d;ierr=0
  end subroutine

  subroutine radioactive_release(emitted,duration_s,delay_s,surviving,decayed,ierr)
    ! Constant release over duration, then free decay over delay to step end.
    ! duration=0 is an instantaneous terminal release at its actual lifetime.
    ! emitted is the fresh MASS in the clipped interval, not a rate or SSP age.
    real(real64),intent(in)::emitted(2),duration_s,delay_s
    real(real64),intent(out)::surviving(2),decayed(2)
    integer,intent(out)::ierr
    real(real64)::s(2),d(2),ms(2),md(2),sd(2),dd(2),unused(2),unused2(2)
    surviving=0;decayed=0;ierr=1
    if(.not.all(ieee_is_finite([emitted,duration_s,delay_s])))return
    if(min(minval(emitted),duration_s,delay_s)<0)return
    call kernel(duration_s,s,d,ms,md)
    call kernel(delay_s,sd,dd,unused,unused2)
    surviving=emitted*(sd*ms)
    decayed=emitted*(dd+sd*md)
    ierr=0
  end subroutine

  function radioactive_element_shift(decayed) result(delta)
    ! A*m_u baryonic convention. Al -> tracked Mg; Fe -> untracked Ni
    ! (Co60 delay collapsed). Total metal/rho, momentum and energy unchanged.
    real(real64),intent(in)::decayed(2)
    real(real64)::delta(11)
    delta=0;delta(7)=decayed(1);delta(11)=-decayed(2)
  end function

  subroutine radioactive_gas_decay(metal,elements,parent,dt_s,ierr)
    ! Atomic abundance transaction. Parent carriers are SUBSETS, never added
    ! to rho/metals. MeV gamma/lepton energy escapes: no thermal/CR/RT source.
    ! Already-aged NEW ejecta must be added after this old-gas operation.
    real(real64),intent(in)::metal,dt_s
    real(real64),intent(inout)::elements(11),parent(2)
    integer,intent(out)::ierr
    real(real64)::surviving(2),decayed(2),trial(11),tol,untracked
    ierr=1
    if(.not.all(ieee_is_finite([metal,elements,parent,dt_s])))return
    if(min(metal,minval(elements),minval(parent),dt_s)<0)return
    untracked=metal-sum(elements(3:11))
    tol=64*epsilon(1d0)*max(metal,sum(elements(3:11)),tiny(1d0))
    if(parent(1)>untracked+tol.or.parent(2)>elements(11)+tol)return
    call radioactive_decay(parent,dt_s,surviving,decayed,ierr)
    if(ierr/=0)return
    ierr=1;trial=elements+radioactive_element_shift(decayed)
    if(.not.all(ieee_is_finite(trial)).or.any(trial<0))return
    if(sum(trial(3:11))>metal+tol)return
    elements=trial;parent=surviving;ierr=0
  end subroutine
end module

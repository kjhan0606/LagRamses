! Bulk metallic Fe enthalpy and common-temperature C/olivine/Fe mixtures.
! Native material building block; no live Fe selector/carriers implied.
module dust_iron_material
  use dust_composition_material, only: dust_composition_curve,dust_composition_cold,dl01_t,dl01_hot_max
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'dust_iron_janaf_data.inc'
  real(real64),parameter :: iron_theta=415d0,iron_gamma=6d-4,iron_join=298.15d0
  real(real64),parameter :: iron_r=8.31446261815324d7,iron_max_t=3000d0
  integer,parameter,public :: iron_identity_n=10+2*iron_n
  public :: iron_enthalpy_bounds,iron_mixture_enthalpy,iron_mixture_state,iron_material_identity
contains
  pure real(real64) function iron_cold_raw(t) result(u)
    ! Hensley & Draine 2017, arXiv:1611.08607 eq(1), erg/g, no zero-point U.
    ! Debye phonons plus electronic heat. Used only below 298.15 K.
    ! Integral_0^z x^3/(exp(x)-1) dx = pi^4/15 - exponential tail series.
    ! Here z>=415/298.15: fast convergence, no small-z cancellation.
    real(real64),intent(in)::t
    real(real64)::z,integral,term,r
    integer::n
    u=0
    if(t<=0)return
    integral=acos(-1d0)**4/15
    if(t>iron_theta/100)then
       z=iron_theta/t
       do n=1,128
          r=real(n,real64)
          term=exp(-r*z)*(z**3/r+3*z*z/r**2+6*z/r**3+6/r**4)
          integral=integral-term
          if(term<epsilon(1d0)*integral)exit
       enddo
    endif
    u=iron_r/iron_molar_mass*(9*t*(t/iron_theta)**3*integral+iron_gamma*t*t/2)
  end function

  pure real(real64) function iron_cold_scale() result(f)
    ! Explicit normalization of the approximate HD17 low-T shape to JANAF's
    ! measured zero-to-298.15K enthalpy. Not a fit to its 100/200K entries.
    f=(-iron_h_zero)*1d10/iron_molar_mass/iron_cold_raw(iron_join)
  end function

  subroutine iron_material_identity(v)
    ! Includes source table (also unused 100/200/3100K entries), conventions,
    ! splice and numerical approximation. No external mutable table at runtime.
    real(real64),intent(out)::v(iron_identity_n)
    ! Version 2 admits the analytic DL01 cold mixture down to U(0)=0.
    v=[2d0,iron_molar_mass,iron_h_zero,iron_theta,iron_gamma,iron_join, &
         iron_r,iron_max_t,iron_cold_scale(),1d0,iron_t,iron_h]
  end subroutine

  subroutine iron_enthalpy_bounds(t,lo,hi,ierr)
    ! U ~= H-H(0) per gram of Fe: condensed pV neglected, equilibrium bulk
    ! phase transitions retained. At a latent plateau [lo,hi] is multivalued.
    ! Piecewise-linear H away from plateaus, NOT interpolated divergent Cp.
    ! Curie anomaly is resolved only at JANAF's tabulated T resolution.
    real(real64),intent(in)::t
    real(real64),intent(out)::lo,hi
    integer,intent(out)::ierr
    integer::i,k
    real(real64)::w
    lo=0;hi=0;ierr=1
    if(.not.ieee_is_finite(t))return
    if(t<0.or.t>iron_max_t)return
    if(t<=iron_join)then
       lo=iron_cold_scale()*iron_cold_raw(t);hi=lo;ierr=0;return
    endif
    do k=1,3
       if(t/=iron_transition_t(k))cycle
       i=iron_transition_index(k)
       lo=(iron_h(i)-iron_h_zero)*1d10/iron_molar_mass
       hi=(iron_h(i+1)-iron_h_zero)*1d10/iron_molar_mass
       ierr=0;return
    enddo
    do i=4,iron_n-1
       if(t<iron_t(i).or.t>iron_t(i+1).or.iron_t(i)==iron_t(i+1))cycle
       w=(t-iron_t(i))/(iron_t(i+1)-iron_t(i))
       lo=((1-w)*iron_h(i)+w*iron_h(i+1)-iron_h_zero)*1d10/iron_molar_mass
       hi=lo;ierr=0;return
    enddo
  end subroutine

  subroutine iron_mixture_enthalpy(t,mass,lo,hi,ierr)
    ! mass = [graphite, MgFeSiO4, metallic Fe] in g (or densities consistently).
    ! Output erg (or erg/cm3). Do not normalize away the actual Fe reservoir.
    real(real64),intent(in)::t,mass(3)
    real(real64),intent(out)::lo,hi
    integer,intent(out)::ierr
    real(real64)::c(1),s(1),a,b,base,cold(2)
    integer::status
    lo=0;hi=0;ierr=1
    if(.not.ieee_is_finite(t).or.any(.not.ieee_is_finite(mass)))return
    if(any(mass<0).or.t<0.or.t>iron_max_t)return
    base=0
    if(any(mass(1:2)>0))then
       if(t>dl01_hot_max)return
       if(t<dl01_t(1))then
          call dust_composition_cold(t,cold,status)
          if(status/=0)return
          base=dot_product(mass(1:2),cold)
       else
       if(mass(1)>0)then
          call dust_composition_curve([t],[1d0,0d0],1d0,c,status)
          if(status/=0)return
          base=mass(1)*c(1)
       endif
       if(mass(2)>0)then
          call dust_composition_curve([t],[0d0,1d0],1d0,s,status)
          if(status/=0)return
          base=base+mass(2)*s(1)
       endif
       endif
    endif
    call iron_enthalpy_bounds(t,a,b,status)
    if(status/=0)return
    a=base+mass(3)*a;b=base+mass(3)*b
    if(.not.all(ieee_is_finite([a,b])))return
    lo=a;hi=b;ierr=0
  end subroutine

  subroutine iron_mixture_state(mass,energy,t,phase,ierr)
    ! Enthalpy inverse with phase mass fractions [alpha,gamma,delta,liquid].
    ! Keep energy itself as state: at melting it fixes liquid fraction, not T.
    ! No hysteresis, supercooling, nanoparticle melting depression or vapor.
    ! Domain failure publishes zero outputs, never a clipped hot/cold state.
    real(real64),intent(in)::mass(3),energy
    real(real64),intent(out)::t,phase(4)
    integer,intent(out)::ierr
    real(real64)::left,right,mid,lo,hi,f,answer,p(4)
    integer::k,j,status
    t=0;phase=0;ierr=1
    if(any(.not.ieee_is_finite(mass)).or..not.ieee_is_finite(energy))return
    if(any(mass<0).or.energy<0.or..not.any(mass>0))return
    left=0
    right=iron_max_t
    call iron_mixture_enthalpy(left,mass,lo,hi,status)
    if(status/=0.or.energy<lo)return
    if(energy==lo)then
       t=left;if(mass(3)>0)phase(1)=1
       ierr=0;return
    endif
    call iron_mixture_enthalpy(right,mass,lo,hi,status)
    if(status/=0.or.energy>hi)return
    if(energy==hi)then
       t=right;if(mass(3)>0)phase(4)=1
       ierr=0;return
    endif
    p=0
    do k=1,3
       call iron_mixture_enthalpy(iron_transition_t(k),mass,lo,hi,status)
       if(status/=0)return
       if(mass(3)>0.and.energy>=lo.and.energy<=hi.and.hi>lo)then
          f=(energy-lo)/(hi-lo);p(k)=1-f;p(k+1)=f
          t=iron_transition_t(k);phase=p;ierr=0;return
       endif
       if(energy<lo)then
          right=iron_transition_t(k);exit
       endif
       left=iron_transition_t(k)
    enddo
    do j=1,100
       mid=left+(right-left)/2
       if(mid==left.or.mid==right)exit
       call iron_mixture_enthalpy(mid,mass,lo,hi,status)
       if(status/=0)return
       if(energy<lo)then
          right=mid
       else
          left=mid
       endif
    enddo
    answer=left+(right-left)/2
    if(mass(3)>0)then
       k=1
       do j=1,3
          if(answer>=iron_transition_t(j))k=j+1
       enddo
       p(k)=1
    endif
    t=answer;phase=p;ierr=0
  end subroutine
end module

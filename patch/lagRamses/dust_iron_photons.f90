! Explicit UV catalytic comparison. No hydro/global state, file IO, charge
! passives, equilibrium gas replacement, or keV/Auger yield extrapolation.
module dust_iron_photons
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_mass_physics, only: dust_fe_trace_charge,dust_fe_trace_energy,dust_fe_relax_ratio, &
       dust_fe_max_temperature
  use dust_iron_optics, only: fe_radius_cm,fe_density,fe_primary_ev,fe_electric_basis,fe_nir
  implicit none
  private
  public :: fe_uv_step,fe_uv_photo,fe_oml,fe_charge_energy,fe_photon_data_identity,fe_photon_data_identity_n
  public :: fe_gas_ip_ev
  integer,parameter :: fe_photon_data_identity_n=95
  real(real64),parameter :: ev=1.602176634d-12,kb=1.380649d-16,mp=1.67262192369d-24
  real(real64),parameter :: pi=3.14159265358979323846d0,e2=14.3996454784255d0
  integer,parameter :: zlo(2)=[-110,-8604],zhi(2)=[63,632],ns=157
  ! Pinned CHIMES157, Fortran indices. Other ions/molecules are unchanged.
  integer,parameter :: ie=1,ih=2,ihp=3,ihe=5,ihep=6,ihepp=7,ic=8,icp=9
  include 'dust_fe_photon_data.inc'
  real(real64),parameter :: kjmol_ev=1000d0/(6.02214076d23*1.602176634d-19)
  ! H/C use the SAME ATcT ground-state reference as CHIMES atomization.
  ! Helium is unchanged and retains the existing CHIMES secondary convention.
  ! FS fractions are energy fractions: events divide by these binding costs,
  ! not by the rounded energy weights internal to the FS interpolation.
  real(real64),parameter :: fe_gas_ip_ev(4)=[ &
       (fe_photon_atct_formation(2)-fe_photon_atct_formation(1))*kjmol_ev,24.59d0,54.42d0, &
       (fe_photon_atct_formation(4)-fe_photon_atct_formation(3))*kjmol_ev]
  real(real64),parameter :: hd_work=4.5d0,hd_ip_positive=.23d0,hd_ip_negative=.34d0
  real(real64),parameter :: hd_bulk=.0021d0,hd_bulk_denom=.0068d0,hd_bulk_power=5d0,hd_escape_length=1d-7
  real(real64),parameter :: hd_barrier_coeff=.3d0,hd_barrier_radius_power=-.45d0,hd_barrier_charge_power=-.26d0
  real(real64),parameter :: hd_det_sigma=1.2d-17,hd_det_width=3d0,hd_det_denom=3d0,hd_det_power=2d0
  real(real64),parameter :: oml_e_stick=.5d0,oml_e_mass=9.1093837139d-28,oml_c_mass=12d0,oml_return=2d0
  real(real64),parameter :: step_fraction=.05d0,tail_bound=1d-12,gas_tmax=1d4
  integer,parameter :: step_cap=10000
  real(real64),parameter :: geom_switch=.01d0,energy_tolerance=2d-10,nuclear_tolerance=2d-12
  real(real64),parameter :: partition_tolerance=1d-10,relax_factor=2d0,variance_floor=1d0
  real(real64),parameter :: uv_max_ev=13.6d0,light_speed=2.99792458d10,gas_cv_factor=1.5d0
  real(real64),parameter :: quad_x(8)=[-.9602898564975363d0,-.7966664774136267d0,-.5255324099163290d0, &
       -.1834346424956498d0,.1834346424956498d0,.5255324099163290d0,.7966664774136267d0,.9602898564975363d0]
  real(real64),parameter :: quad_w(8)=[.1012285362903763d0,.2223810344533745d0,.3137066458778873d0, &
       .3626837833783620d0,.3626837833783620d0,.3137066458778873d0,.2223810344533745d0,.1012285362903763d0]
  abstract interface
     subroutine fe_secondary_callback(energy,state,f,ierr)
       import real64
       real(real64),intent(in)::energy,state(157)
       real(real64),intent(out)::f(5)
       integer,intent(out)::ierr
     end subroutine
  end interface
contains
  function fe_photon_data_identity() result(v)
    real(real64)::v(fe_photon_data_identity_n)
    v=[2d0,e2,hd_work,hd_escape_length,fe_primary_ev,reshape(fe_photon_attenuation_cm,[18]), &
         step_fraction,tail_bound,gas_tmax,real(step_cap,real64),real(zlo,real64),real(zhi,real64), &
         hd_work,hd_ip_positive,hd_ip_negative,hd_bulk,hd_bulk_denom,hd_bulk_power,hd_escape_length, &
         hd_barrier_coeff,hd_barrier_radius_power,hd_barrier_charge_power, &
         hd_det_sigma,hd_det_width,hd_det_denom,hd_det_power, &
         oml_e_stick,oml_e_mass,oml_c_mass,oml_return,kb,mp,ev,pi, &
         geom_switch,energy_tolerance,nuclear_tolerance,partition_tolerance,relax_factor,variance_floor, &
         quad_x,quad_w,fe_photon_atct_formation,kjmol_ev,fe_gas_ip_ev,uv_max_ev,light_speed,gas_cv_factor]
  end function

  real(real64) function fe_charge_energy(z,a) result(f)
    integer,intent(in)::z
    real(real64),intent(in)::a
    real(real64)::d,q
    d=e2/(a*1d8);q=real(z,real64)
    f=hd_work*q+d*(.5d0*q*q-(.5d0-merge(hd_ip_positive,hd_ip_negative,z>=0))*q)
  end function

  real(real64) function ionization(z,a) result(ip)
    integer,intent(in)::z
    real(real64),intent(in)::a
    ip=hd_work+(real(z,real64)+merge(hd_ip_positive,hd_ip_negative,z>=0))*e2/(a*1d8)
  end function

  subroutine fe_uv_photo(z,a,energy,attenuation,yield,pdt,lo,hi,ip,ierr)
    ! HD17 eq32--46; WD01 negative-charge energy interval, not the Z>=0
    ! interval. Correct dimensions of HD17 eq35 using vanHoof04 eq1.
    integer,intent(in)::z
    real(real64),intent(in)::a,energy,attenuation
    real(real64),intent(out)::yield,pdt,lo,hi,ip
    integer,intent(out)::ierr
    real(real64)::d,v,barrier,threshold,theta,y0,y1,y2,alpha,beta,x
    yield=0;pdt=0;lo=0;hi=0;ip=0;ierr=1
    if(.not.all(ieee_is_finite([a,energy,attenuation])))return
    if(a<=0.or.attenuation<=0.or.energy<0.or.energy>uv_max_ev)return
    ip=ionization(z,a);d=e2/(a*1d8);barrier=0
    if(z< -1)then
       v=real(-z-1,real64)
       barrier=d*v/(1+1/sqrt(v))*(1-hd_barrier_coeff*(a/hd_escape_length)**hd_barrier_radius_power* &
            v**hd_barrier_charge_power)
    endif
    threshold=ip+barrier;ierr=0
    if(energy<=threshold)return
    theta=energy-threshold
    if(z>=0)theta=theta+(z+1)*d
    y0=hd_bulk*(theta/hd_work)**hd_bulk_power/(1+hd_bulk_denom*(theta/hd_work)**hd_bulk_power)
    beta=a/attenuation;alpha=beta+a/hd_escape_length
    y1=(beta/alpha)**2*geom(alpha)/geom(beta)
    if(z>=0)then
       lo=-(z+1)*d;hi=energy-threshold
       y2=hi**2*(hi-3*lo)/(hi-lo)**3
    else
       lo=barrier;hi=barrier+energy-threshold;y2=1
       x=(energy-threshold)/hd_det_width
       pdt=-hd_det_sigma*z*x/(1+x*x/hd_det_denom)**hd_det_power
    endif
    yield=y2*min(y0*y1,1d0)
    if(.not.all(ieee_is_finite([yield,pdt,lo,hi,ip])).or.yield<0.or.yield>1.or.pdt<0)ierr=1
  end subroutine

  real(real64) function geom(x) result(g)
    real(real64),intent(in)::x
    ! x^2-2x+2-2exp(-x); stable thin-attenuation limit.
    if(x<geom_switch)then
       g=x**3*(1d0/3-x/12+x*x/60-x**3/360+x**4/2520-x**5/20160)
    else
       g=x*x-2*x+2-2*exp(-x)
    endif
  end function

  subroutine fe_oml(z,a,t,mass,sign,sticking,rate,kinetic,ierr)
    ! Classical OML comparison, NOT the DS87 image-potential fit.
    ! rate is cm3/s; kinetic is mean incoming energy at infinity in eV.
    integer,intent(in)::z,sign
    real(real64),intent(in)::a,t,mass,sticking
    real(real64),intent(out)::rate,kinetic
    integer,intent(out)::ierr
    real(real64)::u,kt,j
    rate=0;kinetic=0;ierr=1
    if(.not.all(ieee_is_finite([a,t,mass,sticking])))return
    if(min(a,t,mass)<=0.or.sticking<0.or.sticking>1.or.abs(sign)/=1)return
    kt=kb*t/ev;u=real(z*sign,real64)*e2/(a*1d8)/kt
    if(u>=0)then
       j=exp(-u);kinetic=kt*(u+2)
    else
       j=1-u;kinetic=kt*(2-u)/(1-u)
    endif
    rate=sticking*pi*a*a*sqrt(8*kb*t/(pi*mass))*j
    if(.not.all(ieee_is_finite([rate,kinetic])))return
    ierr=0
  end subroutine

  subroutine distribution(b,t,ne,nhp,ncp,flux,cross,p,up,down,rh,rc,ke,ki,y,pd,lo,hi,ip,ierr)
    integer,intent(in)::b
    real(real64),intent(in)::t,ne,nhp,ncp,flux(9),cross(9)
    real(real64),intent(out)::p(zlo(b):zhi(b)),up(zlo(b):zhi(b)),down(zlo(b):zhi(b))
    real(real64),intent(out)::rh(zlo(b):zhi(b)),rc(zlo(b):zhi(b)),ke(zlo(b):zhi(b)),ki(zlo(b):zhi(b))
    real(real64),intent(out)::y(4,zlo(b):zhi(b)),pd(4,zlo(b):zhi(b)),lo(4,zlo(b):zhi(b))
    real(real64),intent(out)::hi(4,zlo(b):zhi(b)),ip(zlo(b):zhi(b))
    integer,intent(out)::ierr
    integer::z,g,mode,status
    real(real64)::a,k,s,area,total
    a=fe_radius_cm(b);area=pi*a*a;s=oml_e_stick*(1-exp(-a/hd_escape_length));p=0;ierr=1
    do z=zlo(b),zhi(b)
       call fe_oml(z,a,t,oml_e_mass,-1,s,k,ke(z),status)
       if(status/=0)return
       down(z)=ne*k
       call fe_oml(z,a,t,mp,1,1d0,k,ki(z),status)
       if(status/=0)return
       rh(z)=nhp*k;rc(z)=ncp*k/sqrt(oml_c_mass)
       up(z)=rh(z)+rc(z)
       do g=1,4
          call fe_uv_photo(z,a,fe_primary_ev(g),fe_photon_attenuation_cm(g,b), &
               y(g,z),pd(g,z),lo(g,z),hi(g,z),ip(z),status)
          if(status/=0)return
          up(z)=up(z)+flux(g)*(cross(g)*y(g,z)+pd(g,z))
       enddo
    enddo
    down(zlo(b))=0 ! HD17 autoionization-limited sticking convention.
    mode=zlo(b)
    do z=zlo(b),zhi(b)
       mode=z
       if(up(z)<=down(z).and.down(z)>0)exit
    enddo
    if(mode==zhi(b).and.up(mode)>down(mode))return
    p(mode)=1
    do z=mode-1,zlo(b),-1
       if(p(z+1)==0)exit
       if(up(z)<=0)return
       p(z)=p(z+1)*down(z+1)/up(z)
    enddo
    do z=mode+1,zhi(b)
       if(p(z-1)==0)exit
       if(down(z)<=0)return
       p(z)=p(z-1)*up(z-1)/down(z)
    enddo
    total=sum(p)
    if(.not.ieee_is_finite(total).or.total<=0.or.any(p<0))return
    p=p/total
    ! Explicit negligible stationary-tail bound, not an untested reflecting
    ! ion boundary at the photon-only maximum charge.
    if(p(zhi(b))*up(zhi(b))>tail_bound*sum(p*(up+down)))return
    ierr=0
  end subroutine

  subroutine electron_partition(lo,hi,z,state,secondary,parts,ierr)
    ! Integrate the CONDITIONAL escaping parabolic distribution. Its total
    ! energy is analytic; quadrature averages normalized FS fractions with
    ! energy weights, preserving that exact first moment in the partition.
    real(real64),intent(in)::lo,hi,state(ns)
    integer,intent(in)::z
    procedure(fe_secondary_callback)::secondary
    real(real64),intent(out)::parts(5)
    integer,intent(out)::ierr
    real(real64)::left,e,weight,f(5),mean,total
    integer::j
    parts=0;ierr=0;left=max(0d0,lo)
    if(hi<=left)return
    do j=1,size(quad_x)
       e=left+(hi-left)*(1+quad_x(j))/2
       call secondary(e,state,f,ierr)
       if(ierr/=0)return
       if(any(.not.ieee_is_finite(f)).or.any(f<0).or.abs(sum(f)-1)>partition_tolerance)then
          ierr=1;return
       endif
       weight=quad_w(j)*e*(e-lo)*(hi-e)
       parts=parts+weight*f
    enddo
    if(z>=0)then
       mean=hi*(hi-2*lo)/(2*(hi-3*lo))
    else
       mean=(lo+hi)/2
    endif
    total=sum(parts)
    if(total<=0)then
       ierr=1;return
    endif
    parts=parts*(mean/total)
  end subroutine

  subroutine uv_rates(masses,number,abund,t,td,nh,chat,secondary,dn,ds,gas,solid,exc,qbound,fbound,relax,ierr)
    real(real64),intent(in)::masses(2),number(9),abund(ns),t,td,nh,chat
    procedure(fe_secondary_callback)::secondary
    real(real64),intent(out)::dn(9),ds(ns),gas,solid,exc,qbound,fbound,relax
    integer,intent(out)::ierr
    real(real64),allocatable::p(:),up(:),down(:),rh(:),rc(:),ke(:),ki(:),y(:,:),pd(:,:),lo(:,:),hi(:,:),ip(:)
    real(real64)::ng,area,flux(9),cap,pe,det,cycles,h,c,parts(5),f(5),rate,total,variance,mean,ktout,ion
    real(real64)::cross(9,2),pa(9,2),ps(9,2),pg(9,2),ia(fe_nir,2),isc(fe_nir,2),ig(fe_nir,2)
    integer::b,z,g,n,status
    dn=0;ds=0;gas=0;solid=0;exc=0;qbound=0;fbound=0;relax=0;ierr=1
    call fe_electric_basis(1d0,pa,ps,pg,ia,isc,ig,status)
    if(status/=0)return
    do b=1,2
       cross(:,b)=pa(:,b)*(4*pi/3*fe_density(1)*fe_radius_cm(b)**3)
    enddo
    flux=chat*number;ktout=oml_return*kb*td/ev
    do b=1,2
       if(masses(b)==0)cycle
       area=pi*fe_radius_cm(b)**2;ng=masses(b)/(4*pi/3*fe_density(1)*fe_radius_cm(b)**3)
       n=zhi(b)-zlo(b)+1
       allocate(p(zlo(b):zhi(b)),up(zlo(b):zhi(b)),down(zlo(b):zhi(b)),rh(zlo(b):zhi(b)), &
            rc(zlo(b):zhi(b)),ke(zlo(b):zhi(b)),ki(zlo(b):zhi(b)),y(4,zlo(b):zhi(b)), &
            pd(4,zlo(b):zhi(b)),lo(4,zlo(b):zhi(b)),hi(4,zlo(b):zhi(b)),ip(zlo(b):zhi(b)))
       call distribution(b,t,abund(ie)*nh,abund(ihp)*nh,abund(icp)*nh,flux,cross(:,b), &
            p,up,down,rh,rc,ke,ki,y,pd,lo,hi,ip,status)
       if(status/=0)return
       mean=0
       do z=zlo(b),zhi(b)
          mean=mean+p(z)*z
          qbound=qbound+ng*p(z)*abs(real(z,real64))
          fbound=fbound+ng*p(z)*abs(fe_charge_energy(z,fe_radius_cm(b)))
       enddo
       variance=0;total=0
       do z=zlo(b),zhi(b)
          variance=variance+p(z)*(z-mean)**2
          total=total+p(z)*(up(z)+down(z))
       enddo
       if(total<=0)return
       ! Charge diffusion relaxation ESTIMATE, not a certified spectral gap.
       ! Narrow populations use at least one elementary charging time.
       relax=max(relax,relax_factor*max(variance,variance_floor)/total)
       do g=1,4
          cap=ng*flux(g)*cross(g,b)
          dn(g)=dn(g)-cap;solid=solid+cap*fe_primary_ev(g)
       enddo
       do z=zlo(b),zhi(b)-1
          if(p(z)==0)cycle
          rate=ng*p(z);h=rate*rh(z);c=rate*rc(z);ion=h+c
          ds(ihp)=ds(ihp)-h/nh;ds(ih)=ds(ih)+h/nh
          ds(icp)=ds(icp)-c/nh;ds(ic)=ds(ic)+c/nh
          ds(ie)=ds(ie)-ion/nh
          gas=gas+ion*(ktout-ki(z)-ke(z+1))
          solid=solid+ion*(ki(z)+ke(z+1)-ktout)+h*fe_gas_ip_ev(1)+c*fe_gas_ip_ev(4)
          do g=1,4
             pe=rate*flux(g)*cross(g,b)*y(g,z)
             det=rate*flux(g)*pd(g,z)
             if(pe+det==0)cycle
             call electron_partition(lo(g,z),hi(g,z),z,abund,secondary,parts,status)
             if(status/=0)return
             parts=pe*parts
             if(det>0)then
                f=0
                call secondary(fe_primary_ev(g)-ip(z),abund,f,status)
                if(status/=0.or.any(.not.ieee_is_finite(f)).or.any(f<0).or. &
                     abs(sum(f)-1)>partition_tolerance)return
                parts=parts+det*(fe_primary_ev(g)-ip(z))*f
             endif
             cycles=pe+det
             dn(g)=dn(g)-det
             ! Complete up/down cycle: charge-potential terms cancel exactly.
             gas=gas+parts(1)-cycles*ke(z+1)
             solid=solid+det*fe_primary_ev(g)-sum(parts)+cycles*ke(z+1)
             exc=exc+parts(5)
             f(1:3)=parts(2:4)/fe_gas_ip_ev(1:3)/nh
             ds(ih)=ds(ih)-f(1);ds(ihp)=ds(ihp)+f(1)
             ds(ihe)=ds(ihe)-f(2);ds(ihep)=ds(ihep)+f(2)-f(3);ds(ihepp)=ds(ihepp)+f(3)
             ds(ie)=ds(ie)+sum(f(1:3))
          enddo
       enddo
       deallocate(p,up,down,rh,rc,ke,ki,y,pd,lo,hi,ip)
    enddo
    if(.not.all(ieee_is_finite([dn,ds,gas,solid,exc,qbound,fbound,relax])))return
    ierr=0
  end subroutine

  subroutine fe_uv_step(dt,chat,nh,td,bins,number,state,gas_energy,solid_energy,secondary,ledger,ierr)
    ! All outputs are untouched on failure, including the ledger.
    ! Input state is abundances/H. E are erg/cm3; number is photons/cm3.
    real(real64),intent(in)::dt,chat,nh,td,bins(6)
    real(real64),intent(inout)::number(9),state(ns),gas_energy,solid_energy,ledger(8)
    procedure(fe_secondary_callback)::secondary
    integer,intent(out)::ierr
    real(real64)::n(9),s(ns),dn(9),ds(ns),ug,ud,gd,dd,xd,qb,fb,relax,t
    real(real64)::elapsed,h,force,initial,charge0,charge1,chem0,chem1,balance,excit,captured,trial(8),scale
    integer::iteration,j,status
    ierr=1
    if(.not.all(ieee_is_finite([dt,chat,nh,td,bins,number,state,gas_energy,solid_energy])))return
    if(dt<0.or.chat<=0.or.chat>light_speed.or.nh<=0.or.any(bins<0).or.any(number<0).or.any(state<0))return
    if(gas_energy<=0.or.solid_energy<0)return
    if(td<=0.or.td>dust_fe_max_temperature)return
    if(sum(bins(5:6))==0.or.dt==0)then
       ledger=0;ierr=0;return
    endif
    if(any(number(5:9)>0))then
       ierr=2;return
    endif
    n=number;s=state;ug=gas_energy/ev;ud=solid_energy/ev;elapsed=0;excit=0;trial=0
    initial=ug+ud;charge0=charge(s);chem0=chemical(s)*nh
    do iteration=1,step_cap
       if(elapsed>=dt)exit
       t=ug*ev/(gas_cv_factor*kb*nh*sum(s))
       ! Freeze material T only for this source split. The existing implicit
       ! IR/material solve consumes the energy receipt and enforces T<=300 K.
       if(t<=0.or.t>gas_tmax)then
          ierr=3;return
       endif
       call uv_rates(bins(5:6),n,s,t,td,nh,chat,secondary,dn,ds,gd,dd,xd,qb,fb,relax,status)
       if(status/=0)then
          ierr=4;return
       endif
       if(qb>dust_fe_trace_charge*s(ie)*nh.or.fb>dust_fe_trace_energy*(ug+ud))then
          ierr=5;return
       endif
       ! Omitted-charge relaxation must beat both dt and actual forcing.
       force=dt
       do j=1,9
          if(dn(j)<0)force=min(force,n(j)/(-dn(j)))
       enddo
       do j=1,ns
          if(ds(j)<0)force=min(force,s(j)/(-ds(j)))
       enddo
       if(gd/=0)force=min(force,ug/abs(gd))
       if(dd<0)force=min(force,max(ud,tiny(1d0))/(-dd))
       if(relax>dust_fe_relax_ratio*force)then
          ierr=6;return
       endif
       h=min(dt-elapsed,step_fraction*force)
       if(h<=0.or.elapsed+h==elapsed)then
          ierr=7;return
       endif
       n=n+h*dn;s=s+h*ds;ug=ug+h*gd;ud=ud+h*dd;excit=excit+h*xd
       trial(1)=trial(1)-h*dot_product(dn,fe_primary_ev)
       trial(2)=trial(2)+h*gd;trial(3)=trial(3)+h*dd
       trial(5)=max(trial(5),qb/(s(ie)*nh));trial(6)=max(trial(6),fb/(ug+ud))
       trial(7)=max(trial(7),relax/force)
       if(.not.all(ieee_is_finite([n,s,ug,ud,excit])).or.any(n<0).or.any(s<0).or.ug<=0.or.ud<0)then
          ierr=7;return
       endif
       elapsed=elapsed+h
    enddo
    if(elapsed<dt)then
       ierr=8;return
    endif
    ! The last finite increment can cross a validity bound even when its
    ! starting state passed. Validate the actual candidate, not just stages.
    t=ug*ev/(gas_cv_factor*kb*nh*sum(s))
    if(.not.ieee_is_finite(t).or.t<=0.or.t>gas_tmax)then
       ierr=3;return
    endif
    call uv_rates(bins(5:6),n,s,t,td,nh,chat,secondary,dn,ds,gd,dd,xd,qb,fb,relax,status)
    if(status/=0)then
       ierr=4;return
    endif
    if(qb>dust_fe_trace_charge*s(ie)*nh.or.fb>dust_fe_trace_energy*(ug+ud))then
       ierr=5;return
    endif
    force=dt
    do j=1,9
       if(dn(j)<0)force=min(force,n(j)/(-dn(j)))
    enddo
    do j=1,ns
       if(ds(j)<0)force=min(force,s(j)/(-ds(j)))
    enddo
    if(gd/=0)force=min(force,ug/abs(gd))
    if(dd<0)force=min(force,max(ud,tiny(1d0))/(-dd))
    if(relax>dust_fe_relax_ratio*force)then
       ierr=6;return
    endif
    trial(5)=max(trial(5),qb/(s(ie)*nh));trial(6)=max(trial(6),fb/(ug+ud))
    trial(7)=max(trial(7),relax/force)
    captured=dot_product(number-n,fe_primary_ev);chem1=chemical(s)*nh;charge1=charge(s)
    balance=(ug+ud-initial)+(chem1-chem0)+excit-captured
    scale=max(initial,abs(captured),abs(chem1-chem0),tiny(1d0))
    if(abs(balance)>energy_tolerance*scale.or.abs(charge1-charge0)>nuclear_tolerance*max(sum(state),tiny(1d0)))then
       ierr=9;return
    endif
    if(abs(s(ih)+s(ihp)-state(ih)-state(ihp))>nuclear_tolerance.or. &
         abs(sum(s(ihe:ihepp))-sum(state(ihe:ihepp)))>nuclear_tolerance.or. &
         abs(s(ic)+s(icp)-state(ic)-state(icp))>nuclear_tolerance*max(state(ic)+state(icp),tiny(1d0)))then
       ierr=9;return
    endif
    trial(4)=excit;trial(8)=real(iteration-1,real64)
    number=n;state=s;gas_energy=ug*ev;solid_energy=ud*ev;ledger=trial;ierr=0
  end subroutine

  real(real64) function charge(s) result(q)
    real(real64),intent(in)::s(ns)
    q=s(ihp)+s(ihep)+2*s(ihepp)+s(icp)-s(ie)
  end function
  real(real64) function chemical(s) result(e)
    real(real64),intent(in)::s(ns)
    e=fe_gas_ip_ev(1)*s(ihp)+fe_gas_ip_ev(2)*s(ihep)+(fe_gas_ip_ev(2)+fe_gas_ip_ev(3))*s(ihepp)+ &
         fe_gas_ip_ev(4)*s(icp)
  end function
end module

! Local conservative drag primitive for the unified dust-extension bundle.
! Inputs are separate gas/dust component densities, NOT the barycentric
! RAMSES mixture density. No live selector is exposed by this module alone.
module dust_drag_physics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public::dust_drag_step,dust_epstein_stopping_time
  public::dust_mixture_split,dust_mixture_drag,dust_mixture_radiation_kick
contains
  subroutine dust_mixture_split(rho,momentum,solid,relative,gas_rho,gas_p,solid_p,drift_energy,ierr)
    ! Exact barycentric representation, not the terminal-velocity limit.
    ! rho and momentum are TOTAL gas+solid mass/momentum densities.
    ! relative(:,b) = rho_b*(v_b-v_bary), a SIGNED momentum carrier, not
    ! a mass fraction. It must never pass through a positive scalar limiter.
    ! Total kinetic energy = |momentum|^2/(2*rho) + drift_energy.
    real(real64),intent(in)::rho,momentum(3),solid(:),relative(:,:)
    real(real64),intent(inout)::gas_rho,gas_p(3),solid_p(:,:),drift_energy
    integer,intent(out)::ierr
    real(real64)::rg,pg(3),pd(3,size(solid)),drift,offset(3)
    integer::b,n
    ierr=1;n=size(solid)
    if(any(shape(relative)/=[3,n]).or.any(shape(solid_p)/=[3,n]))return
    if(.not.all(ieee_is_finite([rho,momentum])).or.rho<=0)return
    if(any(.not.ieee_is_finite(solid)).or.any(solid<0))return
    if(any(.not.ieee_is_finite(relative)))return
    rg=rho-sum(solid)
    if(.not.ieee_is_finite(rg).or.rg<=0)return
    offset=sum(relative,dim=2);drift=.5d0*dot_product(offset,offset)/rg
    do b=1,n
       if(solid(b)==0)then
          if(any(relative(:,b)/=0))return
       else
          drift=drift+.5d0*dot_product(relative(:,b),relative(:,b))/solid(b)
       endif
       pd(:,b)=(solid(b)/rho)*momentum+relative(:,b)
    enddo
    pg=(rg/rho)*momentum-offset
    if(.not.all(ieee_is_finite([pg,drift])).or.any(.not.ieee_is_finite(pd)))return
    gas_rho=rg;gas_p=pg;solid_p=pd;drift_energy=drift;ierr=0
  end subroutine

  subroutine dust_mixture_drag(rho,solid,relative,stopping_time,dt,next_relative,gas_heat,ierr)
    ! Internal drag leaves total mixture momentum unchanged. Solve in the
    ! barycentric rest frame so an arbitrarily large common bulk velocity
    ! cannot contaminate drift heating by subtracting two bulk energies.
    ! Frozen stopping times are single-grain times; gas backreaction is in
    ! dust_drag_step. No claim of spatial transport or live RAMSES fields.
    real(real64),intent(in)::rho,solid(:),relative(:,:),stopping_time(:),dt
    real(real64),intent(inout)::next_relative(:,:),gas_heat
    integer,intent(out)::ierr
    real(real64)::rg,pg(3),pd(3,size(solid)),ng(3),nd(3,size(solid)),old_k,new_k,heat
    ierr=1
    if(any(shape(next_relative)/=[3,size(solid)]))return
    rg=0;pg=0;pd=0;old_k=0
    call dust_mixture_split(rho,[0d0,0d0,0d0],solid,relative,rg,pg,pd,old_k,ierr)
    if(ierr/=0)return
    call dust_drag_step(rg,solid,pg,pd,stopping_time,dt,ng,nd,heat,ierr)
    if(ierr/=0)return
    new_k=0
    call dust_mixture_split(rho,[0d0,0d0,0d0],solid,nd,rg,pg,pd,new_k,ierr)
    if(ierr/=0)return
    ierr=1
    if(abs(old_k-new_k-heat)>2d-12*max(old_k,new_k,heat,tiny(1d0)))return
    if(any(abs(pg-ng)>2d-12*max(maxval(abs(ng)),maxval(abs(nd)),tiny(1d0))))return
    next_relative=nd;gas_heat=heat;ierr=0
  end subroutine

  subroutine dust_mixture_radiation_kick(rho,momentum,solid,relative,impulse,radiation_energy, &
       next_momentum,next_relative,solid_heat,ierr)
    ! Receive ALREADY DEBITED radiation momentum and energy per component.
    ! cgs impulse = photon momentum density removed (physical h*nu/c, not
    ! h*nu/c_hat). radiation_energy must include mechanical work: a pure
    ! elastic-scattering force with zero energy debit is not energy closed.
    ! Photoelectron/charge/destruction channels must first be separated by
    ! the physical radiation receiver. No unaccounted thermal energy repair.
    ! Masses fixed during the kick. Empty grains cannot receive an impulse.
    real(real64),intent(in)::rho,momentum(3),solid(:),relative(:,:),impulse(:,:),radiation_energy(:)
    real(real64),intent(inout)::next_momentum(3),next_relative(:,:),solid_heat(:)
    integer,intent(out)::ierr
    real(real64)::rg,pg(3),pd(3,size(solid)),old_k,jtotal(3),work_p(3),work_j(3,size(solid))
    real(real64)::heat(size(solid)),velocity(3),work
    integer::b,n
    ierr=1;n=size(solid)
    if(any(shape(impulse)/=[3,n]).or.any(shape(next_relative)/=[3,n]))return
    if(size(radiation_energy)/=n.or.size(solid_heat)/=n)return
    if(any(.not.ieee_is_finite(impulse)).or.any(.not.ieee_is_finite(radiation_energy)))return
    rg=0;pg=0;pd=0;old_k=0
    call dust_mixture_split(rho,momentum,solid,relative,rg,pg,pd,old_k,ierr)
    if(ierr/=0)return
    ierr=1;jtotal=sum(impulse,dim=2);work_p=momentum+jtotal;work_j=relative
    heat=0
    do b=1,n
       if(solid(b)==0)then
          if(any(impulse(:,b)/=0).or.radiation_energy(b)/=0)return
          cycle
       endif
       velocity=momentum/rho+relative(:,b)/solid(b)
       work=dot_product(velocity,impulse(:,b))+.5d0*dot_product(impulse(:,b),impulse(:,b))/solid(b)
       heat(b)=radiation_energy(b)-work
       if(.not.ieee_is_finite(heat(b)).or.heat(b)<0)return
       work_j(:,b)=relative(:,b)+impulse(:,b)-(solid(b)/rho)*jtotal
    enddo
    if(.not.all(ieee_is_finite(work_p)).or.any(.not.ieee_is_finite(work_j)))return
    next_momentum=work_p;next_relative=work_j;solid_heat=heat;ierr=0
  end subroutine

  subroutine dust_epstein_stopping_time(gas_density,temperature,particle_mass,mean_free_path, &
       radius,solid_density,relative_velocity,stopping_time,ierr)
    ! Neutral free-molecular Epstein drag, with Kwok's supersonic correction
    ! (Laibe & Price 2012b, eq. 17). c_iso^2=k_B*T/m_particle.
    ! Return the SINGLE-GRAIN time rho_d/K, NOT the mixture relative-decay
    ! time rho_g*rho_d/[K*(rho_g+rho_d)]. dust_drag_step owns back-reaction.
    ! cgs inputs. Caller must supply the gas collision mean free path; reject
    ! the Stokes domain a>9*lambda/4 instead of extrapolating this closure.
    ! Does not include Coulomb, Lorentz, sputtering or radiation forces.
    real(real64),intent(in)::gas_density,temperature,particle_mass,mean_free_path
    real(real64),intent(in)::radius(:),solid_density(:),relative_velocity(:,:)
    real(real64),intent(out)::stopping_time(:)
    integer,intent(out)::ierr
    real(real64),parameter::kb=1.380649d-16,pi=acos(-1d0)
    real(real64)::ciso,scale,speed_factor,logt,work(size(radius))
    integer::b,n
    ierr=1;stopping_time=0;n=size(radius)
    if(size(solid_density)/=n.or.size(stopping_time)/=n)return
    if(any(shape(relative_velocity)/=[3,n]))return
    if(.not.all(ieee_is_finite([gas_density,temperature,particle_mass,mean_free_path])))return
    if(min(gas_density,temperature,particle_mass,mean_free_path)<=0)return
    if(.not.all(ieee_is_finite(radius)).or..not.all(ieee_is_finite(solid_density)))return
    if(.not.all(ieee_is_finite(relative_velocity)))return
    if(any(radius<=0).or.any(solid_density<=0))return
    if(any(radius/2.25d0>mean_free_path))return
    logt=.5d0*(log(kb)+log(temperature)-log(particle_mass))
    if(logt<log(tiny(1d0)).or.logt>log(huge(1d0)))return
    ciso=exp(logt)
    do b=1,n
       scale=max(ciso,maxval(abs(relative_velocity(:,b))))
       speed_factor=sqrt((ciso/scale)**2+(9*pi/128)*sum((relative_velocity(:,b)/scale)**2))
       logt=.5d0*log(pi/8)+log(solid_density(b))+log(radius(b))-log(gas_density) &
            -log(scale)-log(speed_factor)
       if(logt<log(tiny(1d0)).or.logt>log(huge(1d0)))return
       work(b)=exp(logt)
    enddo
    stopping_time=work;ierr=0
  end subroutine

  subroutine dust_drag_step(gas_density,dust_density,gas_momentum,dust_momentum, &
       stopping_time,dt,next_gas,next_dust,gas_heat,ierr)
    ! Backward Euler of dp_d/dt=-rho_d*(v_d-v_g)/t_s, with the opposite
    ! gas force. t_s is supplied by the physical drag closure and held fixed
    ! during this step. Momenta in g cm^-2 s^-1, densities in g cm^-3,
    ! times in s. Relative kinetic energy dissipated into gas heat [erg cm^-3].
    ! This does not implement spatial transport, radiation force or charging.
    real(real64),intent(in)::gas_density,dust_density(:),gas_momentum(3),dust_momentum(:,:)
    real(real64),intent(in)::stopping_time(:),dt
    real(real64),intent(out)::next_gas(3),next_dust(:,:),gas_heat
    integer,intent(out)::ierr
    real(real64)::weight(size(dust_density)),remain(size(dust_density))
    real(real64)::vg(3),vd(3,size(dust_density)),new_vd(3),denom,x
    real(real64)::work_g(3),work_d(3,size(dust_density)),delta_g(3),heat
    integer::n,b
    next_gas=gas_momentum;gas_heat=0;ierr=1;n=size(dust_density)
    if(any(shape(dust_momentum)/=[3,n]).or.any(shape(next_dust)/=[3,n]))return
    next_dust=dust_momentum
    if(size(stopping_time)/=n)return
    if(.not.all(ieee_is_finite([gas_density,gas_momentum,dt])))return
    if(.not.all(ieee_is_finite(dust_density)).or..not.all(ieee_is_finite(dust_momentum)))return
    if(.not.all(ieee_is_finite(stopping_time)))return
    if(gas_density<=0.or.dt<0.or.any(dust_density<0).or.any(stopping_time<=0))return
    do b=1,n
       if(dust_density(b)==0.and.any(dust_momentum(:,b)/=0))return
    enddo
    if(dt==0.or.sum(dust_density)==0)then
       ierr=0;return
    endif
    vg=gas_momentum/gas_density;vd=0;weight=0;remain=1
    do b=1,n
       if(dust_density(b)==0)cycle
       vd(:,b)=dust_momentum(:,b)/dust_density(b)
       ! Evaluate dt/(dt+t_s) and t_s/(dt+t_s) without overflowing dt/t_s.
       if(dt>=stopping_time(b))then
          x=stopping_time(b)/dt;weight(b)=1/(1+x);remain(b)=x/(1+x)
       else
          x=dt/stopping_time(b);weight(b)=x/(1+x);remain(b)=1/(1+x)
       endif
    enddo
    denom=gas_density+sum(weight*dust_density)
    ! Solve in the original gas rest frame to avoid subtracting large common
    ! bulk velocities in the kinetic-energy accounting.
    delta_g=0
    do b=1,n
       delta_g=delta_g+weight(b)*dust_density(b)*(vd(:,b)-vg)
    enddo
    delta_g=delta_g/denom
    work_d=dust_momentum;work_g=gas_momentum;heat=.5d0*gas_density*sum(delta_g**2)
    do b=1,n
       if(dust_density(b)==0)cycle
       new_vd=remain(b)*(vd(:,b)-vg)+weight(b)*delta_g
       ! BE's kinetic-energy loss is a sum of nonnegative terms:
       ! dt*K*|v_d'-v_g'|^2 + 1/2 sum rho*|v'-v|^2.
       ! Rewrite the first term using weight*remain to avoid dt/t_s overflow.
       heat=heat+dust_density(b)*weight(b)*remain(b)*sum((vd(:,b)-vg-delta_g)**2) &
            +.5d0*dust_density(b)*sum((new_vd-(vd(:,b)-vg))**2)
       work_d(:,b)=dust_density(b)*(vg+new_vd)
       work_g=work_g-(work_d(:,b)-dust_momentum(:,b))
    enddo
    if(.not.all(ieee_is_finite(work_g)).or..not.all(ieee_is_finite(work_d)))return
    if(.not.ieee_is_finite(heat).or.heat<0)return
    next_gas=work_g;next_dust=work_d;gas_heat=heat;ierr=0
  end subroutine
end module

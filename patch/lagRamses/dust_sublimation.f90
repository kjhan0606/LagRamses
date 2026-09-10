! Opt-in graphite / crystalline olivine vacuum sublimation primitives.
! Used by the old pre-RT split operator and the coupled IR material solver.
! Shared-temperature, fixed-radius two-size closure: not resolved shrinking
! particles, PAH dissociation or vapor-pressure balance.
module dust_sublimation_physics
  use dust_composition_material
  use snrt_dust_ir, only: snrt_dust_material_temperature
  implicit none
contains
  subroutine dust_olivine_surface_flux(temperature,flux,derivative,ierr)
    ! Xu et al., arXiv:2509.11036v1, Table 2: MgFeSiO4 crystalline
    ! surface MASS flux [g cm^-2 s^-1], not a bulk latent heat. The physical
    ! Boltzmann sign is negative (the preprint equation has a double minus).
    ! Equal exposed areas of (100), (010), (001) are an explicit closure.
    ! This kinetic helper alone does not enable a live sublimation mode.
    real(real64),intent(in)::temperature
    real(real64),intent(out)::flux,derivative
    integer,intent(out)::ierr
    real(real64),parameter::log_prefactor(3)=[20.06d0,16.29d0,22.23d0]
    real(real64),parameter::activation_k(3)=[84780d0,74420d0,90590d0]
    real(real64)::surface(3)
    flux=0;derivative=0;ierr=1
    if(.not.ieee_is_finite(temperature).or.temperature<=0)return
    surface=exp(log_prefactor-activation_k/temperature)
    flux=sum(surface)/3
    ! Divide successively, avoiding overflow of T**2 for rejected extremes.
    derivative=sum(surface*(activation_k/temperature)/temperature)/3
    if(.not.all(ieee_is_finite([flux,derivative])))return
    ierr=0
  end subroutine

  subroutine dust_sublimation_step(nodes,bath,bins,energy,dt,next,new_energy,gas_heat,phase_energy,ierr)
    ! Masses and energy may be extensive or densities, but MUST use cgs:
    ! g (or g/cm3), erg (or erg/cm3), seconds, kelvin. No persistent writes.
    ! Conserved energy is Egas+Edust-L*M_C,solid. The equivalent nonnegative
    ! phase reservoir L*(M_C,total-M_C,solid) is derived from existing fields.
    real(real64),intent(in)::nodes(:),bath,bins(4),energy,dt
    real(real64),intent(out)::next(4),new_energy,gas_heat,phase_energy
    integer,intent(out)::ierr
    real(real64)::uc(size(nodes)),us(size(nodes)),mix(size(nodes)),initial_t,lo,hi,t
    real(real64)::trial(4),ed,heat,latent,need,lower_need,tol,removed,removed_s,latent_s
    logical::silicate
    integer::status,it,n
    next=bins;new_energy=energy;gas_heat=0;phase_energy=0;ierr=1;n=size(nodes)
    if(n<2.or..not.dust_material_domain(nodes))return
    if(any(nodes(2:n)<=nodes(1:n-1)))return
    if(.not.all(ieee_is_finite([bath,bins,energy,dt])))return
    if(any(bins<0).or.min(energy,dt)<0.or.bath<nodes(1).or.bath>nodes(n))return
    if(any(.not.ieee_is_finite(dust_size_radius_cm)).or.any(dust_size_radius_cm<=0))return
    if(any(.not.ieee_is_finite(dust_size_density)).or.any(dust_size_density<=0))return
    silicate=dust_silicate_sublimation_enabled();latent_s=0
    if(silicate)then
       call dust_olivine_phase_reference(latent_s,status)
       if(status/=0)return
    endif
    if(sum(bins)==0)then
       if(energy/=0)return
       ierr=0;return
    endif
    if(dt==0.or.(sum(bins(1:2))==0.and..not.silicate))then
       ierr=0;return
    endif
    call dust_composition_curve(nodes,[1d0,0d0],1d0,uc,status)
    if(status/=0)return
    call dust_composition_curve(nodes,[0d0,1d0],1d0,us,status)
    if(status/=0)return
    ! This inequality also ensures the residual increases monotonically
    ! when sublimation shifts mass from solid energy into binding+gas heat.
    if(any(uc>=dust_carbon_latent))return
    if(silicate)then
       if(any(us>=latent_s))return
    endif
    mix=sum(bins(1:2))*uc+sum(bins(3:4))*us
    call snrt_dust_material_temperature(nodes,mix,energy,initial_t,status)
    if(status/=0)return
    ! The validated inverse evaluates exp(log(T)); snap its roundoff at a
    ! table endpoint before passing T to the strict shared rate evaluator.
    if(initial_t<nodes(1)*(1-128*epsilon(1d0)).or.initial_t>nodes(n)*(1+128*epsilon(1d0)))return
    initial_t=max(nodes(1),min(nodes(n),initial_t))
    if(initial_t<bath)return
    lo=bath;hi=initial_t;tol=128*epsilon(1d0)*max(energy,tiny(1d0))
    call evaluate(lo,lower_need)
    if(status/=0)return
    if(lower_need>energy+tol)return ! Do not invent energy from a bath floor.
    call evaluate(hi,need)
    if(status/=0)return
    if(removed==0.and.removed_s==0)then
       ierr=0;return
    endif
    if(need<energy-tol)return
    do it=1,96
       t=lo+(hi-lo)/2
       call evaluate(t,need)
       if(status/=0)return
       if(abs(need-energy)<=tol)exit
       if(need>energy)then
          hi=t
       else
          lo=t
       endif
    enddo
    if(abs(need-energy)>tol)return
    if(.not.all(ieee_is_finite([trial,ed,heat,latent])))return
    if(minval(trial)<0.or.min(ed,heat,latent)<0)return
    next=trial;new_energy=ed;gas_heat=heat;phase_energy=latent;ierr=0
  contains
    subroutine evaluate(temp,total)
      real(real64),intent(in)::temp
      real(real64),intent(out)::total
      call dust_sublimation_at_temperature(nodes,uc,us,bins,temp,dt,latent_s,trial,ed,heat,latent,status)
      total=huge(1d0)
      if(status/=0)return
      removed=sum(bins(1:2)-trial(1:2));removed_s=sum(bins(3:4)-trial(3:4))
      total=ed+heat+latent
    end subroutine
  end subroutine
  subroutine dust_sublimation_at_temperature(nodes,uc,us,bins,temp,dt,latent_s,trial,ed,heat,latent,ierr)
    ! Exact fixed-temperature BE mass/phase evaluation shared by the isolated
    ! and radiation-coupled solves. Caller validates monotone material tables.
    real(real64),intent(in)::nodes(:),uc(:),us(:),bins(4),temp,dt,latent_s
    real(real64),intent(out)::trial(4),ed,heat,latent
    integer,intent(out)::ierr
    real(real64)::w,cc,ss,logx,x,flux,derivative,removed,removed_s,total
    integer::k,b,n,status
    logical::silicate
    ierr=1;n=size(nodes);trial=bins;ed=0;heat=0;latent=0
    if(n<2.or.size(uc)/=n.or.size(us)/=n)return
    if(.not.all(ieee_is_finite([bins,temp,dt,latent_s])))return
    if(any(bins<0).or.dt<=0.or.latent_s<0.or.temp<nodes(1).or.temp>nodes(n))return
    silicate=latent_s>0
      k=1
      do while(k<n-1)
         if(temp<=nodes(k+1))exit
         k=k+1
      enddo
      w=log(temp/nodes(k))/log(nodes(k+1)/nodes(k))
      cc=(1-w)*uc(k)+w*uc(k+1);ss=(1-w)*us(k)+w*us(k+1)
      trial=bins;removed=0;removed_s=0
      do b=1,2
         ! da/dt=-nu*(m/rho_s)^(1/3)*exp(-B/kT); fixed-bin mass rate=3|da/dt|/a.
         logx=log(dt)+log(3*dust_carbon_nu)+(log(dust_carbon_atom/dust_size_density(1)))/3 &
              -log(dust_size_radius_cm(b))-dust_carbon_binding_k/temp
         if(logx>0)then
            x=exp(-logx);trial(b)=bins(b)*x/(1+x)
         else
            x=exp(logx);trial(b)=bins(b)/(1+x)
         endif
         ! Use the represented mass change for exact shared-carrier accounting.
         removed=removed+(bins(b)-trial(b))
      enddo
      if(silicate)then
         call dust_olivine_surface_flux(temp,flux,derivative,status)
         if(status/=0)then
            return
         endif
         if(flux>0)then
            do b=3,4
               logx=log(dt)+log(3d0)+log(flux)-log(dust_size_density(2))-log(dust_size_radius_cm(b-2))
               if(logx>0)then
                  x=exp(-logx);trial(b)=bins(b)*x/(1+x)
               else
                  x=exp(logx);trial(b)=bins(b)/(1+x)
               endif
               removed_s=removed_s+(bins(b)-trial(b))
            enddo
         endif
      endif
      ed=sum(trial(1:2))*cc+sum(trial(3:4))*ss
      heat=removed*2*dust_kb*temp/dust_carbon_atom
      latent=removed*dust_carbon_latent
      if(silicate)then
         ! Congruent atomic Mg+Fe+Si+4O vapor, immediately thermalized.
         heat=heat+removed_s*2*7*dust_kb*temp/(dust_olivine_molar_mass*1.66053906660d-24)
         latent=latent+removed_s*latent_s
      endif

    if(.not.all(ieee_is_finite([trial,ed,heat,latent])))return
    ierr=0
  end subroutine
end module

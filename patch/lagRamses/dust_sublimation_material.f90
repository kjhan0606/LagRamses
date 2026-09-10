! Coupled radiative heating / emission / gas exchange / vacuum evaporation.
! The exact BE mass law is inside the temperature root, not a later kick.
module dust_sublimation_material
  use dust_sublimation_physics
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  private
  public::dust_radiative_sublimation_cell,dust_radiative_sublimation_evolve
  type::material_cache
     real(real64),allocatable::uc(:),us(:),bands(:,:),power(:)
     real(real64)::bath_power,latent_s
  end type
  interface
     function gas_transfer_native(eg,cv,kdt,td) bind(C,name='snrt_dust_gas_transfer_c') result(q)
       import c_double
       real(c_double),value::eg,cv,kdt,td
       real(c_double)::q
     end function
  end interface
contains
  subroutine dust_radiative_sublimation_evolve(nodes,bath,bins,old_energy,heating,dt,reference_mass, &
       basis_band,eg,cv,conductance,next,energy,temperature,rate,transfer,phase,ierr,accuracy,steps)
    ! Native adaptive BE step doubling. A 25% state-change pre-limit prevents
    ! both coarse and half-step solves from missing the same hot transient.
    ! Compare one full step with two half steps; accept the positive half-step
    ! solution, never an extrapolated state or renormalized emission spectrum.
    ! All optical/geometric coefficients remain at the IR-substep input;
    ! gas thermal speed follows the evolving gas reservoir. Rates returned to
    ! the IR receiver are time averages, including the resolved early pulse.
    real(real64),intent(in)::nodes(:),bath,bins(4),old_energy,heating,dt,reference_mass
    real(real64),intent(in)::basis_band(:,:,:),eg,cv,conductance
    real(real64),intent(out)::next(4),energy,temperature,rate(:),transfer,phase
    integer,intent(out)::ierr
    real(real64),optional,intent(in)::accuracy
    integer,optional,intent(out)::steps
    type(material_cache)::cache
    real(real64)::work(4),trial(4),u,ed,gas,q,latent,td,p(size(rate)),mean_p(size(rate))
    real(real64)::elapsed,step,limit,change,coefficient,total_q,total_phase,remaining
    real(real64)::half(4),twice(4),eh,et,qh,qt,lh,lt,th,tt,ph(size(rate)),pt(size(rate)),err,norm
    integer::it,status,accepted
    ierr=1;next=bins;energy=old_energy;temperature=bath;rate=0;transfer=0;phase=0
    if(present(steps))steps=0
    limit=1d-4
    if(present(accuracy))limit=accuracy
    if(.not.ieee_is_finite(limit).or.limit<=0.or.limit>.01d0)return
    if(.not.ieee_is_finite(dt).or.dt<=0.or..not.ieee_is_finite(eg).or.eg<0)return
    work=bins;u=old_energy;gas=eg;elapsed=0;step=dt;mean_p=0;total_q=0;total_phase=0;accepted=0
    do it=1,4096
       remaining=dt-elapsed;step=min(step,remaining)
       if(step<=0.or.elapsed+step==elapsed)exit
       coefficient=conductance
       if(conductance>0.and.eg>0)coefficient=conductance*sqrt(gas/eg)
       call dust_radiative_sublimation_cell(nodes,bath,work,u,heating,step,reference_mass,basis_band, &
            gas,cv,coefficient,trial,ed,td,p,q,latent,status,cache)
       if(status==1)return ! Invalid material/input, not a timestep repair.
       if(status/=0)then
          step=step*.25d0
          cycle
       endif
       change=abs(ed-u)/max(ed,u,tiny(1d0))
       change=max(change,maxval(abs(trial-work)/max(work,trial,tiny(1d0))))
       change=max(change,abs(q)/max(gas,gas-q,tiny(1d0)))
       if(change>.25d0)then
          step=step*max(.1d0,.2d0/change)
          cycle
       endif
       call dust_radiative_sublimation_cell(nodes,bath,work,u,heating,step/2,reference_mass,basis_band, &
            gas,cv,coefficient,half,eh,th,ph,qh,lh,status,cache)
       if(status==1)return
       if(status/=0)then
          step=step*.25d0
          cycle
       endif
       coefficient=conductance
       if(conductance>0.and.eg>0)coefficient=conductance*sqrt((gas-qh)/eg)
       call dust_radiative_sublimation_cell(nodes,bath,half,eh,heating,step/2,reference_mass,basis_band, &
            gas-qh,cv,coefficient,twice,et,tt,pt,qt,lt,status,cache)
       if(status==1)return
       if(status/=0)then
          step=step*.25d0
          cycle
       endif
       err=abs(et-ed)/max(u,et,ed,tiny(1d0))
       err=max(err,maxval(abs(twice-trial)/max(work,twice,tiny(1d0))))
       err=max(err,abs(q-qh-qt)/max(gas,gas-qh-qt,tiny(1d0)))
       norm=max(u,step*heating,abs(qh+qt),tiny(1d0))
       err=max(err,step*sum(abs(p-(ph+pt)/2))/norm)
       if(err>limit)then
          step=step*max(.1d0,.8d0*sqrt(limit/err))
          cycle
       endif
       work=twice;u=et;td=tt;total_q=total_q+qh+qt;gas=eg-total_q;total_phase=total_phase+lh+lt
       mean_p=mean_p+(step/dt)*(ph+pt)/2;accepted=accepted+1
       if(step>=remaining)then
          if(.not.all(ieee_is_finite([work,u,gas,mean_p,total_q,total_phase])))return
          next=work;energy=u;temperature=td;rate=mean_p;transfer=total_q;phase=total_phase;ierr=0
          if(present(steps))steps=accepted
          return
       endif
       elapsed=elapsed+step
       step=step*min(2d0,.9d0*sqrt(limit/max(err,tiny(1d0))))
    enddo
    ierr=5 ! Explicit failure; no partially subcycled state is committed.
  end subroutine

  subroutine dust_radiative_sublimation_cell(nodes,bath,bins,old_energy,heating,dt,reference_mass, &
       basis_band,eg,cv,conductance,next,energy,temperature,rate,transfer,phase,ierr,cache)
    ! cgs densities/energies, band power per reference H. Opacity and gas Cv
    ! are lagged coefficients for this IR substep; mass/U/phase use final T.
    ! energy is SENSIBLE dust energy. The IR trial ledger uses energy+phase,
    ! with transfer=collisional gas-to-dust energy minus effusive vapor heat.
    real(real64),intent(in)::nodes(:),bath,bins(4),old_energy,heating,dt,reference_mass
    real(real64),intent(in)::basis_band(:,:,:),eg,cv,conductance
    real(real64),intent(out)::next(4),energy,temperature,rate(:),transfer,phase
    integer,intent(out)::ierr
    ! Private cache is only reused by evolve with identical radiation inputs.
    ! Its optical coefficients are deliberately frozen at the initial bins.
    type(material_cache),optional,intent(inout)::cache
    real(real64)::uc(size(nodes)),us(size(nodes)),bands(size(rate),size(nodes)),background(size(rate)),power(size(nodes))
    real(real64)::trial(4),trial_rate(size(rate)),ed,heat,latent,coll,ls,lo,hi,t,residual,tol,scale,w
    real(real64)::emitted,bath_power,width,start,finish
    integer::n,k,b,it,status
    ierr=1;next=bins;energy=old_energy;temperature=bath;rate=0;transfer=0;phase=0
    n=size(nodes)
    if(n<2.or.any(shape(basis_band)/=[size(rate),n,4]))return
    if(.not.dust_material_domain(nodes).or.any(nodes(2:)<=nodes(:n-1)))return
    if(.not.all(ieee_is_finite([bath,bins,old_energy,heating,dt,reference_mass,eg,cv,conductance])))return
    if(any(bins<0).or.min(old_energy,heating,eg,conductance)<0.or.min(dt,reference_mass,cv)<=0)return
    if(bath<nodes(1).or.bath>nodes(n))return
    if(sum(bins)==0)then
       if(old_energy/=0.or.heating/=0)return
       ierr=0;return
    endif
    if(present(cache))then
       if(allocated(cache%uc))then
          uc=cache%uc;us=cache%us;bands=cache%bands;power=cache%power
          bath_power=cache%bath_power;ls=cache%latent_s
          goto 100
       endif
    endif
    if(any(.not.ieee_is_finite(basis_band)).or.any(basis_band<0))return
    call dust_composition_curve(nodes,[1d0,0d0],1d0,uc,status)
    if(status/=0)return
    call dust_composition_curve(nodes,[0d0,1d0],1d0,us,status)
    if(status/=0)return
    ls=0
    if(dust_silicate_sublimation_enabled())then
       call dust_olivine_phase_reference(ls,status)
       if(status/=0.or.any(us>=ls))return
    endif
    if(any(uc>=dust_carbon_latent))return
    bands=0
    do b=1,4
       bands=bands+basis_band(:,:,b)*(bins(b)/reference_mass)
    enddo
    if(any(bands(:,2:)<bands(:,:n-1)))return
    k=1
    do while(k<n-1)
       if(bath<=nodes(k+1))exit
       k=k+1
    enddo
    w=log(bath/nodes(k))/log(nodes(k+1)/nodes(k))
    background=(1-w)*bands(:,k)+w*bands(:,k+1)
    power=sum(bands,dim=1);bath_power=sum(background)
    if(any(power(2:)<=power(:n-1)))return
    if(present(cache))then
       cache%uc=uc;cache%us=us;cache%bands=bands;cache%power=power
       cache%bath_power=bath_power;cache%latent_s=ls
    endif
100 continue
    scale=max(old_energy,dt*heating,tiny(1d0));tol=2d-12*scale
    ! Root in NET emitted power, not temperature. A large simulation dt can
    ! require T-Tbath smaller than one representable temperature increment.
    ! Keeping the power increment explicit avoids catastrophic subtraction
    ! of two nearly equal Planck powers near the background floor.
    lo=0
    ierr=2
    call evaluate(lo,residual)
    if(status/=0.or.residual>tol)return
    hi=min(power(n)-bath_power,max(heating+(old_energy-ed-latent-heat+coll)/dt,0d0))
    ierr=3
    call evaluate(hi,residual)
    if(status/=0.or.residual < -tol)return
    ierr=4
    do it=1,96
       emitted=lo+(hi-lo)/2
       call evaluate(emitted,residual)
       if(status/=0)return
       if(abs(residual)<=tol)exit
       if(residual>0)then
          hi=emitted
       else
          lo=emitted
       endif
    enddo
    if(abs(residual)>tol.or.eg-coll+heat<0)return
    trial_rate=0
    do k=1,n-1
       start=max(power(k)-bath_power,0d0);finish=max(power(k+1)-bath_power,0d0)
       ! Above the solved power every later width is exactly zero. Do not
       ! traverse all 136 bands at those unused temperatures in each local
       ! substep; this changes no nonzero arithmetic or physical cutoff.
       if(emitted<=start)exit
       width=max(min(emitted,finish)-start,0d0)
       trial_rate=trial_rate+width*(bands(:,k+1)-bands(:,k))/(power(k+1)-power(k))
    enddo
    if(.not.all(ieee_is_finite(trial_rate)).or.any(trial_rate<0))return
    next=trial;energy=ed;temperature=t;rate=trial_rate;transfer=coll-heat;phase=latent;ierr=0
  contains
    subroutine evaluate(net_power,f)
      real(real64),intent(in)::net_power
      real(real64),intent(out)::f
      k=1
      do while(k<n-1)
         if(bath_power+net_power<=power(k+1))exit
         k=k+1
      enddo
      w=(bath_power+net_power-power(k))/(power(k+1)-power(k))
      t=max(bath,min(nodes(n),exp(log(nodes(k))+w*log(nodes(k+1)/nodes(k)))))
      call dust_sublimation_at_temperature(nodes,uc,us,bins,t,dt,ls,trial,ed,heat,latent,status)
      f=huge(1d0)
      if(status/=0)return
      coll=0
      if(conductance>0)coll=gas_transfer_native(eg,cv,dt*conductance,t)
      f=(ed-old_energy)+latent+heat+dt*(net_power-heating)-coll
      ! A cold grain can radiate far more gas energy than its own stored U.
      ! Scale the residual with the terms actually cancelled in this equation,
      ! not with that tiny stored U alone (nor the entire unused gas reservoir).
      tol=2d-12*max(scale,ed,latent,heat,dt*net_power,abs(coll))
      if(.not.all(ieee_is_finite([f,coll])))status=1
    end subroutine
  end subroutine
end module

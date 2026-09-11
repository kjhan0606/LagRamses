! Excitation-resolved H chemistry for the optional M13/DL01 live comparison.
! Number populations, gas H, vibrational/bond energies and IR are independent
! conserved inventories; fully dehydrogenated C24 is NOT gas-phase carbon.
module dust_pah_hydrogen
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_pah_radiation, only: pah_radiative_model,pah_radiative_coefficients, &
       pah_coronene_photoionize,pah_coronene_recombine
  use dust_stochastic_physics, only: dust_pah_modes,dust_stochastic_photon_rates,dust_stochastic_evolve
  implicit none
  private
  real(real64),parameter::ev=1.602176634d-12
  real(real64),parameter::er_sigma=.06d-16,er_kb=1.380649d-16,er_mh=1.6735575d-24
  ! Molecular ground-state dissociation energy, not an activation barrier.
  real(real64),parameter,public::pah_h2_binding=4.4781d0*ev
  public::pah_hydrogen_parameters,pah_harmonic_dissociation
  public::pah_hydrogen_cool,pah_hydrogen_attach
  public::pah_hydrogen_prepare,pah_hydrogen_absorbed_step
  public::pah_hydrogen_charged_step
  public::pah_h13_abstraction
  public::pah_hydrogen_identity
  type,public::pah_hydrogen_model
     private
     logical::ready=.false.
     real(real64),allocatable::level(:),down(:),band(:,:),diss(:,:)
     real(real64)::bond(0:13),attach(0:13)
     real(real64)::attach_h2(0:13)=0
     logical::h2_enabled=.false.
     logical::h2_abstraction=.false.
     real(real64)::spacing=0
  end type
contains
  subroutine pah_hydrogen_identity(model,values,ierr)
    type(pah_hydrogen_model),intent(in)::model
    real(real64),allocatable,intent(out)::values(:)
    integer,intent(out)::ierr
    ierr=1
    if(.not.model%ready)return
    values=[1d0,model%spacing,model%bond,model%attach,reshape(model%diss,[size(model%diss)])]
    if(model%h2_enabled)values=[2d0,values(2:),model%attach_h2,pah_h2_binding]
    if(model%h2_abstraction)values=[3d0,values(2:),er_sigma,13d0,12d0,0d0,1d0,er_kb,er_mh,8/acos(-1d0)]
    ierr=0
  end subroutine
  subroutine pah_hydrogen_charged_step(model,photon_ev,captured,old,dt,temperature, &
       next,electrons,gas_h,gas_heat,rate,ierr,gas_h2)
    ! Composition of charge, H capture, simultaneous photon/IR/H-loss, and
    ! recombination. Molecules/cm3 indexed (excitation,H+1,charge+1).
    ! Constant IP at every H is an explicit comparison closure. This is not
    ! an H-specific electronic spectrum, nor a multi-charge X-ray network.
    type(pah_hydrogen_model),intent(in)::model(2)
    real(real64),intent(in)::photon_ev(:),captured(:,:),old(:,0:,:),dt,temperature
    real(real64),intent(inout)::next(:,0:,:),electrons,gas_h,gas_heat,rate(:)
    real(real64),optional,intent(inout)::gas_h2
    integer,intent(out)::ierr
    real(real64),allocatable::flat(:,:),vib(:,:),p(:,:,:),trial(:,:),bands(:,:),rec(:,:)
    real(real64)::ne,gh,heat,stored,ip_before,bond_before,bond_after,change,kinetic,attach_heat,db
    real(real64)::rec_ne,rec_heat,rec_ip,np,f,flow,target,w,overflow,expected,scale,balance
    real(real64)::gh2,molecular_change,h_balance
    integer::n,ng,q,nh,j,k,status
    ierr=1
    if(.not.all(model%ready))return
    if(any(model%h2_enabled).neqv.all(model%h2_enabled))return
    if(any(model%h2_abstraction).neqv.all(model%h2_abstraction))return
    if(model(1)%h2_enabled.neqv.present(gas_h2))return
    gh2=0;molecular_change=0
    if(present(gas_h2))then
       if(.not.ieee_is_finite(gas_h2).or.gas_h2<0)return
       gh2=gas_h2
    endif
    n=size(model(1)%level);ng=size(rate)
    if(any(shape(old)/=[n,14,2]).or.any(shape(next)/=shape(old)))return
    if(size(model(2)%level)/=n.or.size(model(1)%band,1)/=ng.or.size(model(2)%band,1)/=ng)return
    if(any(model(2)%level/=model(1)%level))return
    if(any(shape(captured)/=[size(photon_ev),2]))return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    if(.not.all(ieee_is_finite([temperature,dt,electrons,gas_h])))return
    if(dt<=0.or.min(electrons,gas_h)<0)return
    if(sum(old)==0)then
       if(any(.not.ieee_is_finite(captured)).or.any(captured/=0))return
       next=old;rate=0;gas_heat=0;ierr=0;return
    endif
    if(temperature<10.or.temperature>1d4)return
    flat=reshape(old,[n*14,2]);allocate(vib(size(photon_ev),2));vib=0
    ne=electrons;gh=gas_h;heat=0;stored=0
    call pah_coronene_photoionize(photon_ev,captured,flat,ne,vib,heat,stored,status)
    if(status/=0)then
       write(*,*)'PAH H/charge photoionize rejected status=',status
       return
    endif
    allocate(p(n,0:13,2),trial(n,0:13),bands(ng,2));bands=0;p=reshape(flat,[n,14,2])
    ip_before=sum(old(:,:,2))*7.02d0*ev;bond_before=0
    do q=1,2
       bond_before=bond_before+dot_product(model(q)%bond,sum(old(:,:,q),dim=1))
    enddo
    ! Constant k implies sigma*v is velocity independent: captured atoms
    ! sample the Maxwellian mean 3kT/2 (unlike electron alpha~T^-1/2).
    kinetic=1.5d0*1.380649d-16*temperature
    do q=1,2
       trial=p(:,:,q);attach_heat=0;db=0
       call pah_hydrogen_attach(model(q)%level,model(q)%bond,model(q)%attach,p(:,:,q),dt,kinetic, &
            trial,gh,attach_heat,db,status)
       if(status/=0)then
          write(*,*)'PAH H/charge atomic attachment rejected q/status=',q,status
          return
       endif
       heat=heat+attach_heat;p(:,:,q)=trial
       if(model(q)%h2_enabled)then
          call pah_hydrogen_attach(model(q)%level,model(q)%bond,model(q)%attach_h2,p(:,:,q),dt,kinetic, &
               trial,gh2,attach_heat,db,status,capture_stride=2,donor_binding=pah_h2_binding)
          if(status/=0)then
             write(*,*)'PAH H/charge molecular attachment rejected q/status=',q,status
             return
          endif
          heat=heat+attach_heat;p(:,:,q)=trial
       endif
       call pah_hydrogen_absorbed_step(model(q),photon_ev,vib(:,q),p(:,:,q),dt,trial,gh,bands(:,q),db,status)
       if(status/=0)then
          write(*,*)'PAH H/charge absorbed step rejected q/status=',q,status
          ierr=status;return
       endif
       p(:,:,q)=trial
    enddo
    if(all(model%h2_abstraction))then
       ! Explicit frozen-T operator split. Both charges compete for ONE H
       ! donor, after the existing attachment/photo/IR operators. H2 returns
       ! in its ground state; net bond release is immediately thermalized.
       block
         real(real64)::abstraction_heat,breaking(2)
         real(real64),allocatable::abstracted(:,:,:)
         breaking=[model(1)%bond(12)-model(1)%bond(13),model(2)%bond(12)-model(2)%bond(13)]
         abstracted=p;abstraction_heat=0
         call pah_h13_abstraction(p,dt,temperature,breaking,abstracted,gh,gh2,abstraction_heat,status)
         if(status/=0)then
            write(*,*)'PAH H/charge H13 abstraction rejected status=',status
            return
         endif
         p=abstracted;heat=heat+abstraction_heat
       end block
    endif
    ! Use the existing exact bimolecular count ONCE for all H states. A
    ! sequential per-H recombination solve would bias their competing rates.
    allocate(rec(n,2));rec=0;np=sum(p(:,:,2));rec(1,2)=np
    rec_ne=ne;rec_heat=0;rec_ip=np*7.02d0*ev
    call pah_coronene_recombine(model(1)%level,temperature,dt,rec,rec_ne,rec_heat,rec_ip,status)
    if(status/=0)then
       write(*,*)'PAH H/charge recombination rejected status=',status
       return
    endif
    f=0
    if(np>0)f=sum(rec(:,1))/np
    overflow=0;change=0;kinetic=1.380649d-16*temperature
    do nh=0,13
       do j=1,n
          flow=f*p(j,nh,2)
          if(flow==0)cycle
          target=model(1)%level(j)+7.02d0*ev+kinetic
          if(target>model(1)%level(n))then
             overflow=overflow+flow;cycle
          endif
          call bracket(model(1)%level,target,k,w)
          p(j,nh,2)=(1-f)*p(j,nh,2)
          p(k,nh,1)=p(k,nh,1)+flow*(1-w);p(k+1,nh,1)=p(k+1,nh,1)+flow*w
          change=change+flow
       enddo
    enddo
    if(overflow>256*epsilon(1d0)*max(f*np,tiny(1d0)))return
    ne=ne-change;heat=heat-change*kinetic;bond_after=0
    do q=1,2
       bond_after=bond_after+dot_product(model(q)%bond,sum(p(:,:,q),dim=1))
    enddo
    expected=sum(sum(captured,dim=2)*photon_ev)*ev
    if(present(gas_h2))molecular_change=(gas_h2-gh2)*pah_h2_binding
    balance=sum(matmul(model(1)%level,sum(p-old,dim=3)))+bond_after-bond_before &
         +sum(p(:,:,2))*7.02d0*ev-ip_before+heat+dt*sum(bands)-expected+molecular_change
    scale=max(expected,abs(bond_before),abs(bond_after),ip_before,abs(heat), &
         sum(matmul(model(1)%level,sum(old,dim=3))),tiny(1d0))
    if(.not.all(ieee_is_finite([balance,ne,gh,heat,scale])).or.min(ne,gh)<0)return
    if(abs(balance)>4096*epsilon(1d0)*n*28*scale)return
    if(abs(sum(p)-sum(old))>4096*epsilon(1d0)*n*28*max(sum(old),tiny(1d0)))return
    if(abs(sum(p(:,:,2))-ne-sum(old(:,:,2))+electrons)> &
         4096*epsilon(1d0)*n*28*max(sum(old),electrons,tiny(1d0)))return
    if(present(gas_h2))then
       h_balance=gh-gas_h+2*(gh2-gas_h2)
       do nh=0,13
          h_balance=h_balance+nh*sum(p(:,nh,:)-old(:,nh,:))
       enddo
       if(abs(h_balance)>4096*epsilon(1d0)*n*28*max(gas_h,2*gas_h2,13*sum(old),tiny(1d0)))return
       gas_h2=gh2
    endif
    next=p;electrons=ne;gas_h=gh;gas_heat=heat;rate=sum(bands,dim=2);ierr=0
  end subroutine
  subroutine pah_hydrogen_prepare(model,normal,charge,spacing,ierr,h2_capture,h2_abstraction)
    ! Explicit M13/DL01 comparison closure: H-state-dependent generic DL01
    ! harmonic DOS, but normal-H cooling/optics at every H state, as assumed
    ! for cooling in Montillaud+2013 section3.3.1. Not an H-dependent AIB
    ! spectrum or a reproduction of that paper's molecule-specific modes.
    ! H13 uses H11 dissociation DOS, the source's super-H lower-rate bound.
    type(pah_hydrogen_model),intent(out)::model
    type(pah_radiative_model),intent(in)::normal
    integer,intent(in)::charge
    real(real64),intent(in)::spacing
    integer,intent(out)::ierr
    logical,optional,intent(in)::h2_capture
    logical,optional,intent(in)::h2_abstraction
    real(real64),allocatable::modes(:)
    integer::nh,ns,status,nc,normal_h
    ierr=1
    call pah_radiative_coefficients(normal,model%level,model%band,model%down,status,nc,normal_h)
    if(status/=0.or.nc/=24.or.normal_h/=12)return
    call pah_hydrogen_parameters(charge,model%bond,model%attach,status)
    if(status/=0)return
    if(present(h2_capture))model%h2_enabled=h2_capture
    if(present(h2_abstraction))model%h2_abstraction=h2_abstraction
    if(model%h2_abstraction.and..not.model%h2_enabled)return
    ! M13 bound RATE, but only the two-vacancy subset, not an effect bound.
    ! No molecular superhydrogenation or single-vacancy abstraction implied.
    if(model%h2_enabled.and.charge==1)model%attach_h2(0:10)=5d-13
    allocate(model%diss(size(model%level),0:13));model%diss=0
    do nh=1,13
       ns=nh
       if(nh==13)ns=11
       allocate(modes(3*(24+ns-2)))
       call dust_pah_modes(24,ns,modes,status)
       if(status/=0)return
       call pah_harmonic_dissociation(modes,spacing,model%level, &
            model%bond(nh-1)-model%bond(nh),6.8d17,model%diss(:,nh),status)
       if(status/=0)return
       deallocate(modes)
    enddo
    model%spacing=spacing;model%ready=.true.;ierr=0
  end subroutine

  subroutine pah_h13_abstraction(old,dt,temperature,breaking,next,gas_h,gas_h2,gas_heat,ierr)
    ! Boschman2015 Eq3/TableA1: C24H13^(0/+) + H -> C24H12^(0/+) + H2,
    ! one barrierless site, sigma=.06 Angstrom^2. Retain the M13 bond ladder:
    ! this is a named hybrid comparison, not the full Boschman H0--36 model.
    ! Exact bimolecular counts at frozen T, shared by both charges and all
    ! excitation states. No charge change, photon, new carrier or H2 pumping.
    ! PAH excitation is retained; D(H2)-breaking goes into the gas heat.
    real(real64),intent(in)::old(:,0:,:),dt,temperature,breaking(2)
    real(real64),intent(inout)::next(:,0:,:),gas_h,gas_h2,gas_heat
    integer,intent(out)::ierr
    real(real64)::np,nlim,a,delta,x,ratio,loss,denom,events,f,heat,scale,balance
    real(real64)::by_charge(2),transferred(size(old,1),2)
    real(real64),allocatable::trial(:,:,:)
    ierr=1
    if(size(old,1)<1.or.size(old,2)/=14.or.size(old,3)/=2)return
    if(any(shape(next)/=shape(old)))return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    if(.not.all(ieee_is_finite([dt,temperature,gas_h,gas_h2,breaking])))return
    if(dt<0.or.temperature<10.or.temperature>1d4.or.min(gas_h,gas_h2)<0)return
    if(any(breaking<0).or.any(breaking>pah_h2_binding))return
    np=sum(old(:,13,:));nlim=min(np,gas_h)
    if(.not.ieee_is_finite(np))return
    if(dt==0.or.nlim==0)then
       next=old;gas_heat=0;ierr=0;return
    endif
    a=er_sigma*sqrt(8*er_kb*temperature/(acos(-1d0)*er_mh))*dt
    delta=abs(np-gas_h);x=a*delta
    if(.not.all(ieee_is_finite([a,x])))return
    if(x<1d-4)then
       ratio=1-x/2+x*x/6-x**3/24
       loss=x*ratio
    else
       loss=1-exp(-x);ratio=loss/x
    endif
    denom=a*nlim*ratio
    if(.not.ieee_is_finite(denom))return
    events=nlim*(loss+denom)/(1+denom)
    if(events<0.or.events>nlim)return
    f=events/np
    if(.not.ieee_is_finite(f).or.f<0.or.f>1)return
    transferred=f*old(:,13,:)
    by_charge=sum(transferred,dim=1)
    heat=dot_product(by_charge,pah_h2_binding-breaking)
    trial=old
    ! Both factors are nonnegative. Subtracting the rounded transfer from
    ! subnormal H13 tails can acquire a negative sign under optimized
    ! reassociation; form survival directly, without clipping any donor.
    trial(:,13,:)=(1d0-f)*old(:,13,:)
    trial(:,12,:)=old(:,12,:)+transferred
    ! Debit the analytic event count. Population summation differs from it
    ! only at floating-point precision, checked below, never by donor clipping.
    balance=sum(by_charge)-events
    scale=max(np,gas_h,tiny(1d0))
    if(abs(balance)>512*epsilon(1d0)*scale)return
    if(any(.not.ieee_is_finite(trial)).or.any(trial<0))return
    if(.not.all(ieee_is_finite([heat,gas_h2+events])))return
    if(abs(sum(trial)-sum(old))>512*epsilon(1d0)*max(sum(old),tiny(1d0)))return
    next=trial;gas_h=gas_h-events;gas_h2=gas_h2+events;gas_heat=heat;ierr=0
  end subroutine

  subroutine pah_hydrogen_absorbed_step(model,photon_ev,captured,old,dt,next,gas_h,rate,bond_change,ierr)
    ! Actual absorbed photons/cm3, already partitioned by charge and debited
    ! by transport. Simultaneous excitation, IR and H loss in every H state;
    ! a grain may absorb further photons before dissociating. Monochromatic
    ! within-group closure and normal-H optical coefficients are explicit.
    type(pah_hydrogen_model),intent(in)::model
    real(real64),intent(in)::photon_ev(:),captured(:),old(:,0:),dt
    real(real64),intent(inout)::next(:,0:),gas_h,rate(:),bond_change
    integer,intent(out)::ierr
    real(real64),allocatable::up(:,:),over_rate(:),over_power(:),rhs(:,:),p(:,:),removed(:),spectrum(:)
    real(real64)::number,pa,pe,energy,stored,lost,expected,absorbed,emitted,overflow,target,f,scale
    integer::n,g,nh,j,k,status
    ierr=1
    if(.not.model%ready)return
    n=size(model%level);g=size(photon_ev)
    if(any(shape(old)/=[n,14]).or.any(shape(next)/=shape(old)))return
    if(size(rate)/=size(model%band,1).or.size(captured)/=g.or.g<1)return
    if(.not.all(ieee_is_finite([dt,gas_h])).or.dt<=0.or.gas_h<0)return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    if(any(.not.ieee_is_finite(captured)).or.any(captured<0))return
    if(any(.not.ieee_is_finite(photon_ev)).or.any(photon_ev<=0))return
    if(any(photon_ev>13.6d0.and.captured>0))return
    number=sum(old);expected=sum(photon_ev*captured)*ev
    if(.not.all(ieee_is_finite([number,expected])))return
    if(number==0)then
       if(expected/=0)return
       next=old;rate=0;bond_change=0;ierr=0;return
    endif
    allocate(up(n,n),over_rate(n),over_power(n),removed(n),spectrum(size(rate)))
    call dust_stochastic_photon_rates(model%level,photon_ev*ev,captured/number/dt,up,over_rate,over_power,status)
    if(status/=0)return
    rhs=old;p=old;spectrum=0;lost=0;stored=0;absorbed=0;emitted=0;overflow=0
    do nh=13,0,-1
       pa=0;pe=0;removed=0
       call dust_stochastic_evolve(model%level,up,model%down,rhs(:,nh),dt,p(:,nh),pa,pe,status, &
            model%diss(:,nh),removed)
       if(status/=0)return
       absorbed=absorbed+pa;emitted=emitted+pe
       overflow=overflow+dt*dot_product(over_power,p(:,nh))
       spectrum=spectrum+matmul(model%band,p(:,nh))
       do j=1,n
          if(removed(j)==0)cycle
          target=model%level(j)-model%bond(nh-1)+model%bond(nh)
          call bracket(model%level,target,k,f)
          rhs(k,nh-1)=rhs(k,nh-1)+removed(j)*(1-f)
          rhs(k+1,nh-1)=rhs(k+1,nh-1)+removed(j)*f
          lost=lost+removed(j);stored=stored+removed(j)*(model%bond(nh-1)-model%bond(nh))
       enddo
    enddo
    ierr=5
    if(overflow>256*epsilon(1d0)*max(expected,tiny(1d0)))return
    ierr=1
    energy=sum(matmul(model%level,p-old))
    scale=max(expected,emitted,stored,sum(matmul(model%level,old)),tiny(1d0))
    if(.not.all(ieee_is_finite([energy,stored,lost,gas_h+lost,absorbed,emitted,overflow])))return
    if(any(.not.ieee_is_finite(spectrum)).or.any(spectrum<0))return
    if(abs(absorbed+overflow-expected)>4096*epsilon(1d0)*n*14*scale)return
    if(abs(energy+stored+dt*sum(spectrum)-expected)>4096*epsilon(1d0)*n*14*scale)return
    if(abs(sum(p)-number)>4096*epsilon(1d0)*n*14*number)return
    next=p;gas_h=gas_h+lost;rate=spectrum;bond_change=stored;ierr=0
  end subroutine
  subroutine pah_hydrogen_parameters(charge,bond,attach,ierr)
    ! Montillaud+2013 eq3, tables3/4: coronene H=0..13, neutral/monocation.
    ! Bond offsets relative to C24H12 + free H; equal neutral/cation offsets
    ! are the source's approximation. Zero fragment kinetic energy is the
    ! declared threshold-loss closure, NOT a measured fragment spectrum.
    ! Neutral H attachment and H2 addition are the standard model's zero
    ! choices, not measured exclusions. No states above H13 are admitted.
    integer,intent(in)::charge
    real(real64),intent(out)::bond(0:13),attach(0:13)
    integer,intent(out)::ierr
    integer::nh
    ierr=1;bond=0;attach=0
    if(charge<0.or.charge>1)return
    do nh=12,1,-1
       bond(nh-1)=bond(nh)+merge(4.8d0,3.2d0,mod(nh,2)==0)*ev
    enddo
    bond(13)=-3.2d0*ev
    if(charge==1)then
       do nh=0,12
          attach(nh)=merge(1.4d-10,5d-11,mod(nh,2)==0)
       enddo
    endif
    ierr=0
  end subroutine

  subroutine pah_harmonic_dissociation(modes,spacing,levels,threshold,prefactor,rates,ierr)
    ! Beyer--Swinehart count of harmonic states on an explicit uniform mesh.
    ! modes/spacing/levels/threshold in erg. k=A*rho(E-E0)/rho(E).
    ! Rounding mode energies to this mesh is controlled by spacing; callers
    ! must compare refinements. No Arrhenius temperature surrogate or fitted
    ! mean grain temperature substitutes for the excitation distribution.
    real(real64),intent(in)::modes(:),spacing,levels(:),threshold,prefactor
    real(real64),intent(inout)::rates(:)
    integer,intent(out)::ierr
    real(real64),allocatable::rho(:),trial(:)
    real(real64)::top,a,b
    integer::n,j,k,m,quanta
    ierr=1;n=size(levels)
    if(n<2.or.size(rates)/=n.or.size(modes)<1)return
    if(any(.not.ieee_is_finite(modes)).or.any(modes<=0))return
    if(any(.not.ieee_is_finite(levels)).or.any(levels<0))return
    if(any(levels(2:)<=levels(:n-1)))return
    if(.not.all(ieee_is_finite([spacing,threshold,prefactor])))return
    if(spacing<=0.or.threshold<=0.or.prefactor<0)return
    top=levels(n)/spacing
    if(.not.ieee_is_finite(top).or.top>2000000d0)return
    if(minval(modes)/spacing<1.or.maxval(modes)/spacing>real(huge(1)-1,real64))return
    m=ceiling(top)+1;allocate(rho(0:m),trial(n));rho=0;rho(0)=1;trial=0
    do j=1,size(modes)
       quanta=nint(modes(j)/spacing)
       do k=quanta,m
          rho(k)=rho(k)+rho(k-quanta)
       enddo
       if(any(.not.ieee_is_finite(rho)))return
    enddo
    do j=1,n
       if(levels(j)<=threshold)cycle
       a=sample(levels(j));b=sample(levels(j)-threshold)
       if(a<=0)return ! unresolved density of states is not a zero reaction
       trial(j)=prefactor*(b/a)
    enddo
    if(any(.not.ieee_is_finite(trial)).or.any(trial<0))return
    rates=trial;ierr=0
  contains
    real(real64) function sample(e) result(density)
      real(real64),intent(in)::e
      real(real64)::f
      integer::lo
      lo=floor(e/spacing);f=e/spacing-lo
      ! Piecewise linear *density*, not extrapolation through unoccupied bins.
      density=(1-f)*rho(lo)+f*rho(lo+1)
    end function
  end subroutine

  subroutine pah_hydrogen_cool(levels,down,band,diss,bond,old,dt,next,gas_h,emitted,bond_change,ierr)
    ! Simultaneous backward-Euler H loss and IR cooling; no process ordering
    ! that lets every grain cool first and thereby suppresses dissociation.
    ! A downward H/excitation graph permits a positive triangular solve.
    ! A loss emits neutral atomic H and retains the parent's charge. Emission
    ! is integrated erg/cm3 per IR band. All outputs commit together.
    real(real64),intent(in)::levels(:),down(:,0:),band(:,:,0:),diss(:,0:),bond(0:),old(:,0:),dt
    real(real64),intent(inout)::next(:,0:),gas_h,emitted(:),bond_change
    integer,intent(out)::ierr
    real(real64),allocatable::rhs(:,:),p(:,:),light(:)
    real(real64)::x,total,flow,irflow,rate,barrier,target,f,lost,stored,energy,scale,number
    integer::n,nhmax,nh,j,k,g
    ierr=1;n=size(old,1);nhmax=ubound(old,2);g=size(emitted)
    if(n<2.or.nhmax<1.or.size(levels)/=n.or.size(bond)/=nhmax+1)return
    if(any(shape(next)/=shape(old)).or.any(shape(diss)/=shape(old)))return
    if(any(shape(down)/=shape(old)).or.any(shape(band)/=[g,n,nhmax+1]))return
    if(any(.not.ieee_is_finite(levels)).or.levels(1)/=0)return
    if(any(levels(2:)<=levels(:n-1)))return
    if(any(.not.ieee_is_finite(down)).or.any(down<0).or.any(down(1,:)/=0))return
    if(any(.not.ieee_is_finite(band)).or.any(band<0).or.any(band(:,1,:)/=0))return
    if(.not.all(ieee_is_finite([dt,gas_h])).or.dt<0.or.gas_h<0)return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    if(any(.not.ieee_is_finite(diss)).or.any(diss<0).or.any(diss(:,0)/=0))return
    if(any(.not.ieee_is_finite(bond)).or.any(bond(1:)>=bond(0:nhmax-1)))return
    do nh=0,nhmax
       do j=2,n
          total=sum(band(:,j,nh));x=down(j,nh)*(levels(j)-levels(j-1))
          if(abs(total-x)>256*epsilon(1d0)*max(total,x,tiny(1d0)))return
       enddo
       if(nh>0)then
          barrier=bond(nh-1)-bond(nh)
          if(any(diss(:,nh)>0.and.levels<barrier))return
       endif
    enddo
    number=sum(old)
    if(.not.ieee_is_finite(number))return
    if(dt==0.or.number==0)then
       next=old;emitted=0;bond_change=0;ierr=0;return
    endif
    rhs=old;allocate(p(n,0:nhmax),light(g));p=0;light=0;lost=0;stored=0
    do nh=nhmax,0,-1
       do j=n,1,-1
          rate=down(j,nh)+diss(j,nh);x=dt*rate
          if(.not.ieee_is_finite(x))return
          total=rhs(j,nh);p(j,nh)=total/(1+x)
          if(rate==0)cycle
          flow=total*(x/(1+x))
          irflow=flow*(down(j,nh)/rate)
          if(j>1)rhs(j-1,nh)=rhs(j-1,nh)+irflow
          flow=flow*(diss(j,nh)/rate)
          if(flow==0)cycle
          barrier=bond(nh-1)-bond(nh);target=levels(j)-barrier
          call bracket(levels,target,k,f)
          rhs(k,nh-1)=rhs(k,nh-1)+flow*(1-f)
          rhs(k+1,nh-1)=rhs(k+1,nh-1)+flow*f
          lost=lost+flow;stored=stored+flow*barrier
       enddo
       light=light+dt*matmul(band(:,:,nh),p(:,nh))
    enddo
    energy=0
    do nh=0,nhmax
       energy=energy+dot_product(levels,p(:,nh)-old(:,nh))
    enddo
    scale=max(sum(light),stored,sum(matmul(levels,old)),tiny(1d0))
    if(any(.not.ieee_is_finite(p)).or.any(p<0).or.any(.not.ieee_is_finite(light)))return
    if(.not.all(ieee_is_finite([lost,stored,gas_h+lost])))return
    if(abs(sum(p)-number)>4096*epsilon(1d0)*n*(nhmax+1)*number)return
    if(abs(energy+sum(light)+stored)>4096*epsilon(1d0)*n*(nhmax+1)*scale)return
    next=p;gas_h=gas_h+lost;emitted=light;bond_change=stored;ierr=0
  end subroutine

  subroutine pah_hydrogen_attach(level,bond,coefficient,old,dt,kinetic,next,gas_h,gas_heat,bond_change,ierr, &
       capture_stride,donor_binding)
    ! Finite atomic-H donor. Implicit reaction chain, solved together with
    ! the scalar final gas-H inventory (never a post-hoc clipped abundance).
    ! The caller supplies the capture-weighted kinetic energy per H. Bond
    ! release + that energy excite the daughter; gas loses kinetic only.
    ! This is a local frozen-T split, not a gas thermal/charge solve.
    real(real64),intent(in)::level(:),bond(0:),coefficient(0:),old(:,0:),dt,kinetic
    real(real64),intent(inout)::next(:,0:),gas_h,gas_heat,bond_change
    integer,intent(out)::ierr
    integer,optional,intent(in)::capture_stride
    real(real64),optional,intent(in)::donor_binding
    real(real64),allocatable::rhs(:,:),p(:,:)
    real(real64)::lo,hi,mid,residual,captured,heat,stored,overflow,number,energy,scale
    integer::n,nhmax,it,nh,status,stride
    real(real64)::binding
    ierr=1;n=size(level);nhmax=ubound(old,2)
    stride=1;binding=0
    if(present(capture_stride))stride=capture_stride
    if(present(donor_binding))binding=donor_binding
    if(stride<1.or.stride>2.or.stride>nhmax)return
    if(.not.ieee_is_finite(binding).or.binding<0)return
    if(n<2.or.size(old,1)/=n.or.any(shape(next)/=shape(old)).or.nhmax<1)return
    if(size(bond)/=nhmax+1.or.size(coefficient)/=nhmax+1)return
    if(any(.not.ieee_is_finite(level)).or.level(1)/=0)return
    if(any(level(2:)<=level(:n-1)))return
    if(any(.not.ieee_is_finite(bond)).or.any(bond(1:)>=bond(0:nhmax-1)))return
    if(any(.not.ieee_is_finite(coefficient)).or.any(coefficient<0))return
    if(any(coefficient(nhmax-stride+1:)/=0))return
    do nh=0,nhmax-stride
       if(coefficient(nh)>0.and.bond(nh)-bond(nh+stride)<binding)return
    enddo
    if(.not.all(ieee_is_finite([dt,kinetic,gas_h])).or.min(dt,kinetic,gas_h)<0)return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    number=sum(old)
    if(.not.ieee_is_finite(number))return
    if(dt==0.or.number==0.or.gas_h==0.or.all(coefficient==0))then
       next=old;gas_heat=0;bond_change=0;ierr=0;return
    endif
    allocate(rhs(n,0:nhmax),p(n,0:nhmax))
    lo=0;hi=gas_h
    do it=1,180
       mid=lo+(hi-lo)/2
       call trial(mid,status)
       if(status/=0)return
       residual=mid+captured-gas_h
       if(residual>0)then
          hi=mid
       else
          lo=mid
       endif
       if(residual<=0.and.abs(residual)<=128*epsilon(1d0)*max(gas_h,tiny(1d0)))exit
    enddo
    if(abs(residual)>256*epsilon(1d0)*max(gas_h,tiny(1d0)))return
    ! Use the non-overdrawing side of the root bracket, then debit the actual
    ! accepted events, not the independently rounded scalar root.
    mid=lo;call trial(mid,status)
    if(status/=0.or.captured>gas_h)return
    residual=mid+captured-gas_h
    if(abs(residual)>512*epsilon(1d0)*max(gas_h,tiny(1d0)))return
    ! Significant above-grid captures require a larger energy representation.
    ! A roundoff tail remains unreacted, including its donor H.
    if(overflow>256*epsilon(1d0)*max(captured,tiny(1d0)))return
    if(captured==0)then
       next=old;gas_heat=0;bond_change=0;ierr=0;return
    endif
    heat=-captured*kinetic;energy=0;stored=0
    do nh=0,nhmax
       energy=energy+dot_product(level,p(:,nh)-old(:,nh))
       stored=stored+bond(nh)*sum(p(:,nh)-old(:,nh))
    enddo
    scale=max(abs(heat),abs(stored),sum(matmul(level,old)),tiny(1d0))
    if(.not.all(ieee_is_finite([energy,stored,heat,scale])))return
    if(any(.not.ieee_is_finite(p)).or.any(p<0))return
    if(abs(sum(p)-number)>4096*epsilon(1d0)*n*(nhmax+1)*number)return
    if(abs(energy+stored+heat+captured*binding)>4096*epsilon(1d0)*n*(nhmax+1)*scale)return
    next=p;gas_h=gas_h-captured;gas_heat=heat;bond_change=stored;ierr=0
  contains
    subroutine trial(donor,status)
      real(real64),intent(in)::donor
      integer,intent(out)::status
      real(real64)::x,total,flow,target,f
      integer::nh,j,k
      status=1;rhs=old;p=0;captured=0;overflow=0
      do nh=0,nhmax
         x=(dt*coefficient(nh))*donor
         if(.not.ieee_is_finite(x))return
         do j=1,n
            total=rhs(j,nh);p(j,nh)=total/(1+x)
            flow=total*(x/(1+x))
            if(flow==0)cycle
            target=level(j)+bond(nh)-bond(nh+stride)+kinetic
            if(binding/=0)target=target-binding
            if(target>level(n))then
               overflow=overflow+flow;p(j,nh)=total;cycle
            endif
            call bracket(level,target,k,f)
            rhs(k,nh+stride)=rhs(k,nh+stride)+flow*(1-f)
            rhs(k+1,nh+stride)=rhs(k+1,nh+stride)+flow*f
            captured=captured+flow
         enddo
      enddo
      status=0
    end subroutine
  end subroutine

  subroutine bracket(level,target,k,f)
    real(real64),intent(in)::level(:),target
    integer,intent(out)::k
    real(real64),intent(out)::f
    k=1
    do while(k<size(level)-1)
       if(target<=level(k+1))exit
       k=k+1
    enddo
    f=(target-level(k))/(level(k+1)-level(k))
  end subroutine
end module

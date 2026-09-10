! Native PAH population <-> spectral SNRT IR coupling. This is not a
! common-temperature graphite replacement or a stationary-population reset.
! A live hydro caller must additionally own H/C carriers, source chemistry,
! charge/survival admission, co-advection and restart identity.
module dust_pah_radiation
  use dust_stochastic_physics, only: dust_pah_modes,dust_vibrational_curve, &
       dust_stochastic_photon_rates,dust_stochastic_evolve
  use snrt_dust_ir, only: dust_ir_table,dust_ir_diagnostics,snrt_dust_ir_initialize,snrt_dust_ir_advance
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(real64),parameter::ev=1.602176634d-12,kb=1.380649d-16,h=6.62607015d-27,c=2.99792458d10
  type,public::pah_radiative_model
     private
     logical::ready=.false.
     integer::nc=0,nh=0
     real(real64),allocatable::level(:),temperature(:),band(:,:),down(:),ir_ev(:),cabs(:)
     type(dust_ir_table)::transport
  end type
  public::pah_radiative_prepare,pah_radiative_advance,pah_absorbed_step,pah_solid_inventory
  public::pah_neutral_optics
  public::pah_charge_optics,pah_charged_absorbed_step
  public::pah_photon_advance
  public::pah_coronene_photoionize,pah_coronene_recombine
  public::pah_radiative_coefficients
contains
  subroutine pah_radiative_coefficients(model,levels,band,down,ierr,ncarbon,nhydrogen)
    ! Read-only copies for the coupled H-state receiver; preserve private
    ! prepared-model invariants (no external mutation of its cooling rates).
    type(pah_radiative_model),intent(in)::model
    real(real64),allocatable,intent(out)::levels(:),band(:,:),down(:)
    integer,intent(out)::ierr
    integer,optional,intent(out)::ncarbon,nhydrogen
    ierr=1
    if(present(ncarbon))ncarbon=-1
    if(present(nhydrogen))nhydrogen=-1
    if(.not.model%ready)return
    levels=model%level;band=model%band;down=model%down;ierr=0
    if(present(ncarbon))ncarbon=model%nc
    if(present(nhydrogen))nhydrogen=model%nh
  end subroutine
  subroutine pah_coronene_photoionize(photon_ev,captured,population,electrons,vibrational_captures, &
       gas_heat,ionization_energy,ierr)
    ! C24H12, neutral/monocation event receiver. Montillaud+2013 eq1,
    ! sections 3.2--3.3: ionizing absorptions leave vibrational energy
    ! unchanged; the excess above IP becomes photoelectron kinetic energy.
    ! captured(g,charge) is ALREADY absorbed photons/cm3, not incident flux.
    ! Charge-resolved opacity and excitation evolution belong to the caller.
    ! No stationary charge reset, high-energy fit extrapolation, H loss or
    ! second ionization is implied. This bounded PDR kernel is not itself
    ! admission of charged PAHs into the live neutral-only hydro model.
    real(real64),intent(in)::photon_ev(:),captured(:,:)
    real(real64),intent(inout)::population(:,:),electrons,vibrational_captures(:,:)
    real(real64),intent(inout)::gas_heat,ionization_energy
    integer,intent(out)::ierr
    real(real64),parameter::ip=7.02d0
    real(real64)::ions(size(photon_ev)),vib(size(photon_ev),2),trial(size(population,1),2)
    real(real64)::n0,ni,f,heat,stored,ne
    integer::g
    ierr=1
    if(size(population,1)<1.or.size(population,2)/=2.or.size(photon_ev)<1)return
    if(any(shape(captured)/=[size(photon_ev),2]))return
    if(any(shape(vibrational_captures)/=shape(captured)))return
    if(any(.not.ieee_is_finite(photon_ev)).or.any(photon_ev<=0))return
    if(any(.not.ieee_is_finite(captured)).or.any(captured<0))return
    if(any(photon_ev>13.6d0.and.sum(captured,dim=2)>0))return
    if(any(.not.ieee_is_finite(population)).or.any(population<0))return
    if(.not.all(ieee_is_finite([electrons,gas_heat,ionization_energy])).or.electrons<0)return
    n0=sum(population(:,1))
    if(.not.ieee_is_finite(n0))return
    if(n0==0.and.any(captured(:,1)>0))return
    if(sum(population(:,2))==0.and.any(captured(:,2)>0))return
    ions=0
    do g=1,size(photon_ev)
       if(photon_ev(g)>=ip.and.photon_ev(g)<=13.6d0) &
            ions(g)=captured(g,1)*.8d0*exp(-.00128d0*(photon_ev(g)-14.89d0)**4)
    enddo
    ni=sum(ions)
    ! If the frozen neutral population is exhausted, retry the enclosing
    ! transport/chemistry step at smaller dt. Never clip away photon energy.
    if(.not.ieee_is_finite(ni).or.ni>n0)return
    trial=population
    if(ni>0)then
       f=ni/n0
       trial(:,1)=(1-f)*population(:,1)
       trial(:,2)=population(:,2)+f*population(:,1)
    endif
    vib=captured;vib(:,1)=vib(:,1)-ions
    ne=electrons+ni;heat=gas_heat+sum(ions*(photon_ev-ip))*ev
    stored=ionization_energy+ni*ip*ev
    if(.not.all(ieee_is_finite([ne,heat,stored])).or.any(.not.ieee_is_finite(trial)))return
    population=trial;electrons=ne;vibrational_captures=vib
    gas_heat=heat;ionization_energy=stored;ierr=0
  end subroutine

  subroutine pah_coronene_recombine(levels,temperature,dt,population,electrons,gas_heat,ionization_energy,ierr)
    ! Montillaud+2013 table4: alpha=1e-5 sqrt(300/T) cm3/s (a MODEL rate,
    ! not a measured coronene recombination coefficient). Non-dissociative
    ! capture deposits IP + electron kinetic energy in the neutral grain.
    ! Maxwellian electrons with sigma*v proportional to E^(-1/2) give a
    ! capture-weighted mean kinetic energy kT. This is the explicit kinetic
    ! assumption underlying our T^(-1/2) rate, not the 3kT/2 bulk mean.
    ! Exact bimolecular cation/electron counts at frozen gas T; caller owns
    ! the thermal solve. Interpolation preserves molecule number and energy.
    real(real64),intent(in)::levels(:),temperature,dt
    real(real64),intent(inout)::population(:,:),electrons,gas_heat,ionization_energy
    integer,intent(out)::ierr
    real(real64),parameter::ip=7.02d0
    real(real64)::trial(size(levels),2),np,nlim,delta,a,x,loss,ratio,denom,nrec,f,n,target,w
    real(real64)::heat,stored,kinetic,unresolved,accepted,limit
    integer::j,k,nlevel
    ierr=1;nlevel=size(levels)
    if(nlevel<2.or.any(shape(population)/=[nlevel,2]))return
    if(any(.not.ieee_is_finite(levels)).or.levels(1)/=0)return
    if(any(levels(2:)<=levels(:nlevel-1)))return
    if(any(.not.ieee_is_finite(population)).or.any(population<0))return
    if(.not.all(ieee_is_finite([temperature,dt,electrons,gas_heat,ionization_energy])))return
    ! Explicit bounded PDR comparison domain, not a general hot-plasma fit.
    if(temperature<10.or.temperature>1d4.or.dt<0.or.electrons<0)return
    np=sum(population(:,2));nlim=min(np,electrons)
    if(.not.ieee_is_finite(np))return
    if(dt==0.or.nlim==0)then
       ierr=0;return
    endif
    a=1d-5*sqrt(300d0/temperature)*dt
    delta=abs(np-electrons);x=a*delta
    if(.not.all(ieee_is_finite([a,x])))return
    if(x<1d-4)then
       ratio=1-x/2+x*x/6-x**3/24
       loss=x*ratio
    else
       loss=1-exp(-x);ratio=loss/x
    endif
    denom=a*nlim*ratio
    if(.not.ieee_is_finite(denom))return
    ! Stable loss form also resolves small dt and exactly equal abundances.
    nrec=nlim*(loss+denom)/(1+denom)
    if(nrec<0.or.nrec>nlim)return
    f=nrec/np;kinetic=kb*temperature
    trial=population;unresolved=0;accepted=0
    do j=1,nlevel
       n=f*population(j,2)
       if(n==0)cycle
       target=levels(j)+ip*ev+kinetic
       if(target>levels(nlevel))then
          ! Do not invent an above-grid state or clip its energy. Retain
          ! these unreacted molecules; only roundoff-scale event tails may
          ! be unresolved, as in the vibrational absorption receiver.
          unresolved=unresolved+n*(ip*ev+kinetic)
          cycle
       endif
       accepted=accepted+n;trial(j,2)=(1-f)*population(j,2)
       k=1
       do while(k<nlevel-1)
          if(target<=levels(k+1))exit
          k=k+1
       enddo
       w=(target-levels(k))/(levels(k+1)-levels(k))
       trial(k,1)=trial(k,1)+n*(1-w);trial(k+1,1)=trial(k+1,1)+n*w
    enddo
    limit=256*epsilon(1d0)*max(nrec*(ip*ev+kinetic),tiny(1d0))
    if(unresolved>limit)return
    nrec=accepted
    heat=gas_heat-nrec*kinetic;stored=ionization_energy-nrec*ip*ev
    if(any(.not.ieee_is_finite(trial)).or..not.all(ieee_is_finite([heat,stored])))return
    population=trial;electrons=electrons-nrec;gas_heat=heat;ionization_energy=stored;ierr=0
  end subroutine

  subroutine pah_neutral_optics(path,ncarbon,photon_ev,weight_ev,cabs,ierr)
    character(len=*),intent(in)::path
    integer,intent(in)::ncarbon
    real(real64),allocatable,intent(out)::photon_ev(:),weight_ev(:),cabs(:)
    integer,intent(out)::ierr
    call pah_charge_optics(path,ncarbon,0,photon_ev,weight_ev,cabs,ierr)
  end subroutine

  subroutine pah_charge_optics(path,ncarbon,charge_state,photon_ev,weight_ev,cabs,ierr)
    ! Read the uncompressed ORIGINAL Draine PAHneu/ion_30 table. No extrapolation
    ! in grain radius or wavelength and no invented far-IR/X-ray tail.
    ! Source defines N_C=468*(a/.001 micron)^3. Log interpolation of actual
    ! cross sections between bracketing radii; retain ALL wavelength knots.
    ! These are phenomenological PAH-mixture optical data, NOT independent
    ! molecule-specific Malloci spectra. Charge must match the file header;
    ! neither table supplies arbitrary hydrogenation or multiple charges.
    character(len=*),intent(in)::path
    integer,intent(in)::ncarbon,charge_state
    real(real64),allocatable,intent(out)::photon_ev(:),weight_ev(:),cabs(:)
    integer,intent(out)::ierr
    integer,parameter::nr=30,nw=1201
    real(real64)::radius(nr),wave(nw),qabs(nw,nr),w,qext,qa,qs,asym,target,f,xlo,xhi
    integer::unit,ios,i,j,k,nrad,nwave
    character(len=256)::line
    ierr=1
    if(ncarbon<3.or.charge_state<0.or.charge_state>1)return
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,'(A)',iostat=ios)line
    if(ios/=0)goto 900
    if(charge_state==0.and.index(line,'neutral PAH/graphitic solid')/=1)goto 900
    if(charge_state==1.and.index(line,'ionized PAH/graphitic solid')/=1)goto 900
    do i=2,7
       read(unit,'(A)',iostat=ios)line
       if(ios/=0)goto 900
    enddo
    read(unit,*,iostat=ios)nrad,xlo,xhi
    if(ios/=0)goto 900
    if(nrad/=nr.or.xlo/=3.548d-4.or.xhi/=1d-2)goto 900
    read(unit,*,iostat=ios)nwave,xlo,xhi
    if(ios/=0)goto 900
    if(nwave/=nw.or.xlo/=1d3.or.xhi/=1d-3)goto 900
    do i=1,nr
       read(unit,'(A)',iostat=ios)line
       if(ios/=0)goto 900
       read(unit,*,iostat=ios)radius(i)
       if(ios/=0)goto 900
       read(unit,'(A)',iostat=ios)line
       if(ios/=0)goto 900
       do j=1,nw
          ! Original fixed-width Fortran output: negative g can directly
          ! follow Q_sca with NO whitespace. List-directed reads misparse it.
          read(unit,'(E9.3,3E10.3,E9.2)',iostat=ios)w,qext,qa,qs,asym
          if(ios/=0)goto 900
          if(.not.all(ieee_is_finite([w,qext,qa,qs,asym])))goto 900
          if(min(w,qa)<=0.or.min(qext,qs)<0.or.abs(asym)>1)goto 900
          if(i==1)wave(j)=w
          if(w/=wave(j))goto 900
          qabs(j,i)=qa
       enddo
    enddo
    if(any(.not.ieee_is_finite(radius)).or.any(radius<=0))goto 900
    if(any(radius(2:)<=radius(:nr-1)).or.any(wave(2:)>=wave(:nw-1)))goto 900
    target=.001d0*(real(ncarbon,real64)/468)**(1d0/3)
    if(target<radius(1).or.target>radius(nr))goto 900
    k=1
    do while(k<nr-1)
       if(target<=radius(k+1))exit
       k=k+1
    enddo
    f=log(target/radius(k))/log(radius(k+1)/radius(k))
    allocate(photon_ev(nw),weight_ev(nw),cabs(nw))
    photon_ev=(h*c/ev)*1d4/wave
    weight_ev(1)=(photon_ev(2)-photon_ev(1))/2
    weight_ev(nw)=(photon_ev(nw)-photon_ev(nw-1))/2
    weight_ev(2:nw-1)=(photon_ev(3:)-photon_ev(:nw-2))/2
    cabs=acos(-1d0)*exp((1-f)*log(qabs(:,k)*(radius(k)*1d-4)**2)+ &
         f*log(qabs(:,k+1)*(radius(k+1)*1d-4)**2))
    ierr=0
900 close(unit)
  end subroutine

  subroutine pah_radiative_prepare(model,ncarbon,nhydrogen,levels,ir_ev,weight_ev,cabs,ierr)
    ! DL01 eqs 2,32,41: thermal-continuous SPONTANEOUS cooling. Integrate
    ! only photon energies <= the grain excitation; no emission above U.
    ! No rescaling of the spectrum to a separately prescribed bolometric L.
    ! Full spectral/source coverage and survival are caller admission duties.
    ! cabs [cm2/grain], weight_ev [eV], levels [erg/grain]. Absolute IR field,
    ! not the equilibrium receiver's excess-above-a-bath convention.
    type(pah_radiative_model),intent(out)::model
    integer,intent(in)::ncarbon,nhydrogen
    real(real64),intent(in)::levels(:),ir_ev(:),weight_ev(:),cabs(:)
    integer,intent(out)::ierr
    real(real64),allocatable::modes(:),sorted(:)
    real(real64)::left,right,mid,u(1),cv(1),t,x,occupation,photon,first,twentieth
    integer::n,ng,nm,i,j,g,k,status
    ierr=1;n=size(levels);ng=size(ir_ev)
    if(ncarbon<3.or.nhydrogen<0.or.ncarbon>100000.or.nhydrogen>100000)return
    if(n<2.or.ng<1.or.size(weight_ev)/=ng.or.size(cabs)/=ng)return
    if(any(.not.ieee_is_finite(levels)).or.any(levels<0))return
    if(levels(1)/=0.or.any(levels(2:)<=levels(:n-1)))return
    if(any(.not.ieee_is_finite(ir_ev)).or.any(.not.ieee_is_finite(weight_ev)))return
    if(any(.not.ieee_is_finite(cabs)).or.any(cabs<0))return
    if(any(ir_ev<=0).or.any(weight_ev<=0).or.any(ir_ev(2:)<=ir_ev(:ng-1)))return
    nm=3*(ncarbon+nhydrogen-2)
    if(nm<20)return
    allocate(modes(nm),sorted(nm))
    call dust_pah_modes(ncarbon,nhydrogen,modes,status)
    if(status/=0)return
    ! Only the first and twentieth order statistics are needed (DL01 eq32).
    ! Partial selection avoids sorting a large molecular spectrum O(Nm^2).
    sorted=modes
    do i=1,20
       k=i-1+minloc(sorted(i:),dim=1);t=sorted(i);sorted(i)=sorted(k);sorted(k)=t
    enddo
    first=sorted(1);twentieth=sorted(20)
    allocate(model%level(n),model%temperature(n),model%band(ng,n),model%down(n),model%ir_ev(ng))
    model%level=levels;model%ir_ev=ir_ev;model%temperature=0;model%band=0;model%down=0
    model%cabs=cabs
    do j=2,n
       if(levels(j)<=twentieth)then
          t=first/(kb*log(2d0))
       else
          ! U(T)<=Nm*kT; doubling from this lower estimate brackets the inverse.
          left=0;right=max(levels(j)/(nm*kb),1d0)
          do i=1,128
             call dust_vibrational_curve(modes,[right],u,cv,status)
             if(status/=0)return
             if(u(1)>=levels(j))exit
             right=2*right
          enddo
          if(u(1)<levels(j))return
          do i=1,100
             mid=left+(right-left)/2
             if(mid==left.or.mid==right)exit
             call dust_vibrational_curve(modes,[mid],u,cv,status)
             if(status/=0)return
             if(u(1)>levels(j))then
                right=mid
             else
                left=mid
             endif
          enddo
          t=left+(right-left)/2
       endif
       model%temperature(j)=t
       do g=1,ng
          photon=ir_ev(g)*ev
          if(photon>levels(j))cycle
          x=photon/(kb*t)
          if(x<1d-3)then
             occupation=1/x-.5d0+x/12-x**3/720
          else
             occupation=exp(-x)/(1-exp(-x))
          endif
          model%band(g,j)=8*acos(-1d0)*photon**3/(h**3*c**2)*cabs(g)*weight_ev(g)*ev*occupation
       enddo
       model%down(j)=sum(model%band(:,j))/(levels(j)-levels(j-1))
    enddo
    if(any(.not.ieee_is_finite(model%band)).or.any(model%band<0))return
    if(any(.not.ieee_is_finite(model%down)).or.any(model%down(2:)<=0))return
    ! These two auxiliary thermal nodes initialize the EXISTING transport
    ! table only. Its equilibrium emission/bath are never called below.
    call snrt_dust_ir_initialize(model%transport,ir_ev,weight_ev,cabs,[10d0,20d0],10d0,status)
    if(status/=0)return
    model%nc=ncarbon;model%nh=nhydrogen;model%ready=.true.;ierr=0
  end subroutine

  subroutine pah_charged_absorbed_step(model,photon_ev,captured,old,electrons,gas_temperature,dt, &
       next,next_electrons,gas_heat,ionization_change,rate,ierr)
    ! Conservative first-order local receiver, not a transport/advection
    ! driver: photoionization -> charge-resolved vibrational absorption and
    ! IR cooling -> electron recombination. Captures have already been
    ! partitioned/debited using the TWO optical cross sections by the caller.
    ! Gas T is frozen during this split; the caller must couple the returned
    ! heat and particle-number change to its gas thermal solve. Stationary
    ! grains only: there is no mechanical-work or Doppler prescription here.
    ! Both prepared models use the same DL01 modes (explicit approximation)
    ! with distinct optical/emission cross sections. H loss is not included.
    type(pah_radiative_model),intent(in)::model(2)
    real(real64),intent(in)::photon_ev(:),captured(:,:),old(:,:),electrons,gas_temperature,dt
    real(real64),intent(inout)::next(:,:),next_electrons,gas_heat,ionization_change,rate(:)
    integer,intent(out)::ierr
    real(real64),allocatable::trial(:,:),cooled(:,:),vib(:,:),spectrum(:,:),fluence(:)
    real(real64)::ne,heat,stored,initial_stored,u,t,overflow,absorbed,balance,tolerance
    integer::n,ng,j,status
    ierr=1
    if(.not.all(model%ready))return
    if(any(model%nc/=24).or.any(model%nh/=12))return
    n=size(model(1)%level);ng=size(model(1)%ir_ev)
    if(size(model(2)%level)/=n.or.size(model(2)%ir_ev)/=ng)return
    if(any(model(1)%level/=model(2)%level).or.any(model(1)%ir_ev/=model(2)%ir_ev))return
    if(any(shape(old)/=[n,2]).or.any(shape(next)/=[n,2]).or.size(rate)/=ng)return
    if(any(shape(captured)/=[size(photon_ev),2]))return
    if(.not.ieee_is_finite(dt).or.dt<=0)return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    trial=old;ne=electrons;heat=0;initial_stored=sum(old(:,2))*7.02d0*ev;stored=initial_stored
    if(sum(old)==0)then
       if(any(.not.ieee_is_finite(captured)).or.any(captured/=0))return
       if(.not.ieee_is_finite(electrons).or.electrons<0)return
       next=old;next_electrons=electrons;gas_heat=0;ionization_change=0;rate=0;ierr=0
       return
    endif
    allocate(vib(size(photon_ev),2),spectrum(ng,2));vib=0;spectrum=0
    call pah_coronene_photoionize(photon_ev,captured,trial,ne,vib,heat,stored,status)
    if(status/=0)return
    cooled=trial
    do j=1,2
       fluence=vib(:,j)*photon_ev*ev;u=0;t=0
       call pah_absorbed_step(model(j),photon_ev,fluence,trial(:,j),dt,cooled(:,j),spectrum(:,j), &
            u,t,overflow,status,captured_photons=vib(:,j))
       if(status/=0)then
          ierr=status;return
       endif
    enddo
    call pah_coronene_recombine(model(1)%level,gas_temperature,dt,cooled,ne,heat,stored,status)
    if(status/=0)return
    absorbed=sum(sum(captured,dim=2)*photon_ev)*ev
    balance=dot_product(model(1)%level,sum(cooled-old,dim=2))+stored-initial_stored+heat &
         +dt*sum(spectrum)-absorbed
    tolerance=4096*epsilon(1d0)*n*max(absorbed,abs(heat),abs(stored),abs(initial_stored), &
         dot_product(model(1)%level,sum(old,dim=2)),tiny(1d0))
    if(.not.ieee_is_finite(balance).or.abs(balance)>tolerance)return
    if(abs(sum(cooled)-sum(old))>2048*epsilon(1d0)*n*max(sum(old),tiny(1d0)))return
    if(abs(sum(cooled(:,2))-ne-sum(old(:,2))+electrons)> &
         2048*epsilon(1d0)*n*max(sum(old),electrons,tiny(1d0)))return
    next=cooled;next_electrons=ne;gas_heat=heat;ionization_change=stored-initial_stored
    rate=sum(spectrum,dim=2);ierr=0
  end subroutine

  subroutine pah_photon_advance(model,direction,weight,neighbor,dx,dt,chat,primary_ev,primary_photons, &
       radiation,population,emitted_photons,diagnostics,ierr,tolerance,max_iterations, &
       ghost_energy,ghost_index,blocked_face)
    ! Excitation-independent PAH absorption of an already transported primary
    ! angular photon reservoir (photons/cm3 per normalized direction). Uses
    ! the SAME optical cross sections as emission, log-interpolated in E.
    ! Primary photons, IR photons and populations commit as one transaction.
    ! Other absorbers must be partitioned by the enclosing mixed-material
    ! caller; this PAH-only receiver must not debit their photons a second time.
    type(pah_radiative_model),intent(in)::model
    real(real64),intent(in)::direction(:,:),weight(:),dx,dt,chat,primary_ev(:)
    integer,intent(in)::neighbor(:,:)
    real(real64),intent(inout)::primary_photons(:,:,:),radiation(:,:,:),population(:,:),emitted_photons(:,:)
    type(dust_ir_diagnostics),intent(inout)::diagnostics
    integer,intent(out)::ierr
    real(real64),intent(in)::tolerance
    integer,intent(in)::max_iterations
    real(real64),optional,intent(in)::ghost_energy(:,:,:)
    integer,optional,intent(in)::ghost_index(:,:)
    logical,optional,intent(in)::blocked_face(:,:)
    real(real64),allocatable::trial(:,:,:),captured(:,:),number(:),sigma(:)
    real(real64)::f,tau,loss
    integer::ng,nc,nd,g,i,d,k,nopt
    ierr=1
    if(.not.model%ready)return
    ng=size(primary_ev);nc=size(population,2);nd=size(weight);nopt=size(model%ir_ev)
    if(ng<1.or.nc<1.or.nd<1.or.nopt<2)return
    if(any(shape(primary_photons)/=[ng,nd,nc]))return
    if(any(.not.ieee_is_finite(primary_photons)).or.any(primary_photons<0))return
    if(any(.not.ieee_is_finite(primary_ev)).or.any(primary_ev<=0))return
    if(any(.not.ieee_is_finite(weight)).or.any(weight<=0))return
    if(abs(sum(weight)-1)>1d-12)return
    if(any(.not.ieee_is_finite(population)).or.any(population<0))return
    if(.not.all(ieee_is_finite([dt,chat])).or.min(dt,chat)<=0)return
    allocate(sigma(ng),captured(ng,nc));sigma=0;captured=0
    number=sum(population,dim=1);trial=primary_photons
    if(any(.not.ieee_is_finite(number)))return
    do g=1,ng
       ! Zero populations do not authorize extrapolated future absorption:
       ! any occupied, nonzero primary group needs supported cross sections.
       if(.not.any(primary_photons(g,:,:)>0))cycle
       if(primary_ev(g)<model%ir_ev(1).or.primary_ev(g)>model%ir_ev(nopt))return
       k=1
       do while(k<nopt-1)
          if(primary_ev(g)<=model%ir_ev(k+1))exit
          k=k+1
       enddo
       f=log(primary_ev(g)/model%ir_ev(k))/log(model%ir_ev(k+1)/model%ir_ev(k))
       if(model%cabs(k)<=0.or.model%cabs(k+1)<=0)return
       sigma(g)=exp((1-f)*log(model%cabs(k))+f*log(model%cabs(k+1)))
       do i=1,nc
          tau=chat*dt*number(i)*sigma(g)
          if(.not.ieee_is_finite(tau))return
          if(tau<1d-4)then
             loss=tau*(1-tau/2+tau*tau/6-tau**3/24)
          else
             loss=1-exp(-tau)
          endif
          do d=1,nd
             trial(g,d,i)=primary_photons(g,d,i)*exp(-tau)
             captured(g,i)=captured(g,i)+weight(d)*primary_photons(g,d,i)*loss*primary_ev(g)*ev
          enddo
       enddo
    enddo
    call pah_radiative_advance(model,direction,weight,neighbor,dx,dt,chat,primary_ev,captured, &
         radiation,population,emitted_photons,diagnostics,ierr,tolerance,max_iterations, &
         ghost_energy,ghost_index,blocked_face)
    if(ierr==0)primary_photons=trial
  end subroutine

  subroutine pah_solid_inventory(model,population,hc,excitation,ierr)
    ! H/C solid mass densities from the SAME number distribution as the
    ! radiation receiver. Hydrogen is not metallicity; neither is new rho.
    type(pah_radiative_model),intent(in)::model
    real(real64),intent(in)::population(:)
    real(real64),intent(inout)::hc(2),excitation
    integer,intent(out)::ierr
    real(real64)::number,mass(2),u
    ierr=1
    if(.not.model%ready)return
    if(size(population)/=size(model%level))return
    if(any(.not.ieee_is_finite(population)).or.any(population<0))return
    number=sum(population)
    mass=number*1.66053906660d-24*[1.008d0*model%nh,12.011d0*model%nc]
    u=dot_product(model%level,population)
    if(.not.all(ieee_is_finite([mass,u])))return
    hc=mass;excitation=u;ierr=0
  end subroutine

  subroutine pah_absorbed_step(model,photon_ev,absorbed_fluence,old,dt,next,rate,energy,temperature,overflow,ierr, &
       captured_photons)
    ! Convert actual captured spectral energy to per-grain photon events.
    ! Excitation-independent Cabs is the DL01 assumption: every grain sees
    ! the same absorption rate. old/next are number/cm3 per energy state.
    ! absorbed_fluence is erg/cm3 PER SPECTRAL NODE, not a bolometric total.
    ! overflow is always reported (erg/cm3). A significant top-grid tail
    ! REJECTS, preserving all inout state; it is never thermalized or deleted.
    type(pah_radiative_model),intent(in)::model
    real(real64),intent(in)::photon_ev(:),absorbed_fluence(:),old(:),dt
    ! Optional actual captures / cm3, independently of deposited excitation
    ! energy. Moving-grain absorption pays mechanical work BEFORE this call.
    ! Changing that work must not manufacture or remove absorption events.
    ! Within each spectral node the retained closure is a mean excitation
    ! quantum = deposited energy / captured photons (not frequency-resolved).
    real(real64),optional,intent(in)::captured_photons(:)
    real(real64),intent(inout)::next(:),rate(:),energy,temperature
    real(real64),intent(out)::overflow
    integer,intent(out)::ierr
    real(real64),allocatable::up(:,:),over_rate(:),over_power(:),events(:),p(:),spectrum(:),quantum(:)
    real(real64)::number,pa,pe,expected,tolerance,u,t
    integer::n,ng,status
    ierr=1;overflow=0
    if(.not.model%ready)return
    n=size(model%level);ng=size(photon_ev)
    if(size(old)/=n.or.size(next)/=n.or.size(rate)/=size(model%ir_ev))return
    if(ng<1.or.size(absorbed_fluence)/=ng)return
    if(.not.ieee_is_finite(dt).or.dt<=0)return
    if(any(.not.ieee_is_finite(old)).or.any(old<0))return
    if(any(.not.ieee_is_finite(absorbed_fluence)).or.any(absorbed_fluence<0))return
    if(any(.not.ieee_is_finite(photon_ev)).or.any(photon_ev<=0))return
    if(present(captured_photons))then
       if(size(captured_photons)/=ng)return
       if(any(.not.ieee_is_finite(captured_photons)).or.any(captured_photons<0))return
       if(any(captured_photons==0.and.absorbed_fluence>0))return
    endif
    number=sum(old);expected=sum(absorbed_fluence)
    if(.not.all(ieee_is_finite([number,expected])))return
    if(number==0)then
       if(expected/=0)return
       if(present(captured_photons))then
          if(any(captured_photons/=0))return
       endif
       next=0;rate=0;energy=0;temperature=0;ierr=0;return
    endif
    allocate(up(n,n),over_rate(n),over_power(n),events(ng),p(n),spectrum(size(rate)),quantum(ng))
    events=(absorbed_fluence/number)/(photon_ev*ev)/dt
    quantum=photon_ev*ev
    if(present(captured_photons))then
       events=0
       where(captured_photons>0.and.absorbed_fluence>0)
          events=(captured_photons/number)/dt
          quantum=absorbed_fluence/max(captured_photons,tiny(1d0))
       endwhere
       ! Captures with zero excitation leave the population unchanged. Their
       ! photon debit and mechanical work remain in the enclosing receiver.
    endif
    if(any(.not.ieee_is_finite(events)).or.any(.not.ieee_is_finite(quantum)))return
    call dust_stochastic_photon_rates(model%level,quantum,events,up,over_rate,over_power,status)
    if(status/=0)return
    p=old;pa=0;pe=0
    call dust_stochastic_evolve(model%level,up,model%down,old,dt,p,pa,pe,status)
    if(status/=0)return
    overflow=dt*dot_product(over_power,p)
    if(.not.ieee_is_finite(overflow))return
    ! Numerical truncation must be below roundoff-scale absorbed energy;
    ! no user knob can conceal a missing high-energy/survival model here.
    ierr=5
    if(overflow>256*epsilon(1d0)*max(expected,tiny(1d0)))return
    spectrum=matmul(model%band,p)
    u=dot_product(model%level,p);t=dot_product(model%temperature,p)/number
    if(any(.not.ieee_is_finite(spectrum)).or.any(spectrum<0))return
    if(.not.all(ieee_is_finite([u,t])))return
    tolerance=2048*epsilon(1d0)*n*max(expected,pe,u,dot_product(model%level,old),tiny(1d0))
    if(abs(pa+overflow-expected)>tolerance)return
    if(abs(dt*sum(spectrum)-pe)>tolerance)return
    if(abs(dot_product(model%level,p-old)-expected+dt*sum(spectrum))>tolerance)return
    next=p;rate=spectrum;energy=u;temperature=t;ierr=0
  end subroutine

  subroutine pah_radiative_advance(model,direction,weight,neighbor,dx,dt,chat,primary_ev,primary_absorbed, &
       radiation,population,emitted_photons,diagnostics,ierr,tolerance,max_iterations, &
       ghost_energy,ghost_index,blocked_face)
    ! Connected SNRT transport/absorption/emission transaction. Primary
    ! absorbed photons must have been debited by the enclosing source trial;
    ! on failure that caller must discard that trial as well. This routine
    ! commits no populations/radiation/diagnostics on failure. Reabsorption
    ! retains its spectrum at every nonlinear iterate; it is not bolometric.
    type(pah_radiative_model),intent(in)::model
    real(real64),intent(in)::direction(:,:),weight(:),dx,dt,chat,primary_ev(:),primary_absorbed(:,:)
    integer,intent(in)::neighbor(:,:)
    real(real64),intent(inout)::radiation(:,:,:),population(:,:),emitted_photons(:,:)
    type(dust_ir_diagnostics),intent(inout)::diagnostics
    integer,intent(out)::ierr
    real(real64),intent(in)::tolerance
    integer,intent(in)::max_iterations
    real(real64),optional,intent(in)::ghost_energy(:,:,:)
    integer,optional,intent(in)::ghost_index(:,:)
    logical,optional,intent(in)::blocked_face(:,:)
    real(real64),allocatable::density(:),energy(:),temperature(:),capacity(:),primary(:)
    integer::nc,i
    ierr=1
    if(.not.model%ready)return
    nc=size(population,2)
    if(size(population,1)/=size(model%level).or.nc<1)return
    if(any(shape(primary_absorbed)/=[size(primary_ev),nc]))return
    if(.not.ieee_is_finite(dt).or.dt<=0)return
    if(any(.not.ieee_is_finite(primary_absorbed)).or.any(primary_absorbed<0))return
    allocate(density(nc),energy(nc),temperature(nc),capacity(nc),primary(nc))
    density=sum(population,dim=1);temperature=0;capacity=1
    do i=1,nc
       energy(i)=dot_product(model%level,population(:,i))
    enddo
    primary=sum(primary_absorbed,dim=1)/dt
    call snrt_dust_ir_advance(model%transport,direction,weight,neighbor,dx,dt,chat,density,primary, &
         radiation,temperature,emitted_photons,diagnostics,ierr,tolerance,max_iterations, &
         dust_energy=energy,heat_capacity=capacity,ghost_energy=ghost_energy,ghost_index=ghost_index, &
         blocked_face=blocked_face,thin_reabsorption=.true.,population=population,population_dispatch=material)
  contains
    subroutine material(ir_absorbed,step_dt,old,next,rate,next_energy,next_temperature,status)
      real(real64),intent(in)::ir_absorbed(:,:),step_dt,old(:,:)
      real(real64),intent(out)::next(:,:),rate(:,:),next_energy(:),next_temperature(:)
      integer,intent(out)::status
      real(real64)::overflow
      integer::cell
      next=old;rate=0;next_energy=0;next_temperature=0
      do cell=1,size(old,2)
         call pah_absorbed_step(model,[primary_ev,model%ir_ev], &
              [primary_absorbed(:,cell),ir_absorbed(:,cell)],old(:,cell),step_dt, &
              next(:,cell),rate(:,cell),next_energy(cell),next_temperature(cell),overflow,status)
         if(status/=0)return
      enddo
      status=0
    end subroutine
  end subroutine
end module

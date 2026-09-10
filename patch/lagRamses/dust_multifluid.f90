! Conservative, same-mesh gas + pressureless grain transport.
! This is a spatial operator, not yet a selectable RAMSES hydro scheme.
! State: total rho, total p(3), E; then [rho_b,p_b(3)] per grain;
! then positive carriers tagged by their transport owner (0=gas).
! E = gas thermal + nonthermal + kinetic energy of ALL components.
! Grain thermal/excitation carriers are outside E, as in DUST_LIVE.
! Absolute grain momenta are conservative. rho_b*(v_b-v_bary) is NOT.
module dust_multifluid
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_drag_physics, only: dust_drag_step
  implicit none
  private
  integer,parameter,public::dust_fv_ok=0,dust_fv_config=1,dust_fv_state=2,dust_fv_cfl=3
  type,public::dust_fv_layout
     private
     integer::nb=0,nc=0,nvar=0,base=0
     integer,allocatable::owner(:)
     real(real64),allocatable::gamma_nt(:)
     logical::ready=.false.
  end type
  public::dust_fv_initialize,dust_fv_decode,dust_fv_face,dust_fv_periodic_step,dust_fv_drag
  public::dust_fv_periodic_advance
  public::dust_fv_radiation_kick,dust_fv_absorption_kick
contains
  subroutine dust_fv_radiation_kick(layout,state,gamma,impulse,debited_energy,next,solid_heat,ierr)
    ! Fixed masses, absolute component momenta. Received energy is ALREADY
    ! removed from radiation. E_hydro receives only mechanical work; the
    ! remainder belongs to the grain enthalpy/PAH receiver outside E_hydro.
    ! Rejection leaves BOTH state and returned heating unchanged.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::state(:),gamma,impulse(:,:),debited_energy(:)
    real(real64),intent(inout)::next(:),solid_heat(:)
    integer,intent(out)::ierr
    real(real64)::rho(0:layout%nb),v(3,0:layout%nb),th,kin,p,cs,mechanical,th_after
    real(real64)::trial(size(state)),heat(layout%nb),work(layout%nb)
    integer::b,k
    ierr=dust_fv_config
    if(size(next)/=size(state).or.size(solid_heat)/=layout%nb.or.size(debited_energy)/=layout%nb)return
    if(any(shape(impulse)/=[3,layout%nb]))return
    if(any(.not.ieee_is_finite(impulse)).or.any(.not.ieee_is_finite(debited_energy)))return
    call dust_fv_decode(layout,state,gamma,rho,v,th,kin,p,cs,ierr)
    if(ierr/=0)return
    ierr=dust_fv_state;trial=state;heat=0;work=0
    do b=1,layout%nb
       if(rho(b)==0)then
          if(any(impulse(:,b)/=0).or.debited_energy(b)/=0)return
          cycle
       endif
       k=6+4*(b-1)
       work(b)=dot_product(v(:,b)+.5d0*impulse(:,b)/rho(b),impulse(:,b))
       heat(b)=debited_energy(b)-work(b)
       if(.not.ieee_is_finite(heat(b)).or.heat(b)<0)return
       trial(k+1:k+3)=state(k+1:k+3)+impulse(:,b)
    enddo
    mechanical=sum(work)
    trial(2:4)=state(2:4)+sum(impulse,dim=2)
    trial(5)=state(5)+mechanical
    call dust_fv_decode(layout,trial,gamma,rho,v,th_after,kin,p,cs,ierr)
    if(ierr/=0)return
    ierr=dust_fv_state
    if(abs(th_after-th)>1024*epsilon(1d0)*max(abs(state(5)),abs(trial(5)),tiny(1d0)))return
    next=trial;solid_heat=heat;ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_absorption_kick(layout,state,gamma,photons,moment,energy_erg,fraction, &
       number_scale,energy_scale,momentum_scale,next,solid_heat,ierr)
    ! Adapter from the primary SNRT accepted dust photon ledger and its
    ! signed first angular moment. These are angle-integrated bins, not
    ! intensity: no angular weights, dt, or reduced-c factor is applied here.
    ! fraction(b,g) is the PHYSICAL absorption share from each grain opacity.
    ! Caller must not assign scattering, returned photons, gas ionization or
    ! already diverted photoelectron energy to this absorption-only receiver.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::state(:),gamma,photons(:),moment(:,:),energy_erg(:),fraction(:,:)
    real(real64),intent(in)::number_scale,energy_scale,momentum_scale
    real(real64),intent(inout)::next(:),solid_heat(:)
    integer,intent(out)::ierr
    real(real64),parameter::physical_c=2.99792458d10
    real(real64)::impulse(3,layout%nb),energy(layout%nb),norm,scale,e
    integer::g,b,ng
    ierr=dust_fv_config;ng=size(photons)
    if(.not.layout%ready.or.ng<1)return
    if(size(energy_erg)/=ng.or.any(shape(moment)/=[3,ng]))return
    if(any(shape(fraction)/=[layout%nb,ng]))return
    if(.not.all(ieee_is_finite([number_scale,energy_scale,momentum_scale])))return
    if(min(number_scale,energy_scale,momentum_scale)<=0)return
    if(any(.not.ieee_is_finite(photons)).or.any(photons<0))return
    if(any(.not.ieee_is_finite(moment)))return
    if(any(.not.ieee_is_finite(energy_erg)).or.any(energy_erg<=0))return
    if(any(.not.ieee_is_finite(fraction)).or.any(fraction<0))return
    impulse=0;energy=0
    do g=1,ng
       scale=maxval(abs(moment(:,g)));norm=0
       if(scale>0)norm=scale*sqrt(sum((moment(:,g)/scale)**2))
       if(norm>(1d0+2d-6)*photons(g))return ! FP32 quadrature boundary only.
       if(photons(g)==0)cycle
       if(abs(sum(fraction(:,g))-1d0)>128*epsilon(1d0))return
       e=number_scale*energy_erg(g)
       do b=1,layout%nb
          energy(b)=energy(b)+(e/energy_scale)*photons(g)*fraction(b,g)
          impulse(:,b)=impulse(:,b)+(e/physical_c/momentum_scale)*moment(:,g)*fraction(b,g)
       enddo
    enddo
    if(any(.not.ieee_is_finite(energy)).or.any(.not.ieee_is_finite(impulse)))return
    call dust_fv_radiation_kick(layout,state,gamma,impulse,energy,next,solid_heat,ierr)
  end subroutine

  subroutine dust_fv_initialize(nb,owner,gamma_nt,layout,ierr)
    integer,intent(in)::nb,owner(:)
    real(real64),intent(in)::gamma_nt(:) ! 0=passive; >1=gas nonthermal energy
    type(dust_fv_layout),intent(inout)::layout
    integer,intent(out)::ierr
    type(dust_fv_layout)::trial
    ierr=dust_fv_config
    if(nb<1.or.size(owner)/=size(gamma_nt))return
    if(any(owner<0).or.any(owner>nb))return
    if(any(.not.ieee_is_finite(gamma_nt)))return
    if(any(gamma_nt/=0.and.gamma_nt<=1))return
    if(any(gamma_nt>0.and.owner/=0))return
    trial%nb=nb;trial%nc=size(owner);trial%base=5+4*nb;trial%nvar=trial%base+trial%nc
    trial%owner=owner;trial%gamma_nt=gamma_nt;trial%ready=.true.
    layout=trial;ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_decode(layout,u,gamma,rho,velocity,thermal,kinetic,pressure,sound,ierr)
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::u(:),gamma
    real(real64),intent(out)::rho(0:),velocity(:,0:),thermal,kinetic,pressure,sound
    integer,intent(out)::ierr
    real(real64)::p(3),nt,cs2
    integer::b,k,j
    ierr=dust_fv_config;rho=0;velocity=0;thermal=0;kinetic=0;pressure=0;sound=0
    if(.not.layout%ready)return
    if(size(u)/=layout%nvar.or.size(rho)/=layout%nb+1)return
    if(any(shape(velocity)/=[3,layout%nb+1]))return
    if(.not.ieee_is_finite(gamma).or.gamma<=1)return
    ierr=dust_fv_state
    if(any(.not.ieee_is_finite(u)).or.u(1)<=0)return
    p=u(2:4)
    do b=1,layout%nb
       k=6+4*(b-1);rho(b)=u(k)
       if(rho(b)<0)return
       if(rho(b)==0)then
          if(any(u(k+1:k+3)/=0))return
       else
          velocity(:,b)=u(k+1:k+3)/rho(b)
          kinetic=kinetic+.5d0*dot_product(u(k+1:k+3),velocity(:,b))
       endif
       p=p-u(k+1:k+3)
    enddo
    rho(0)=u(1)-sum(rho(1:))
    if(rho(0)<=0.or..not.ieee_is_finite(rho(0)))return
    velocity(:,0)=p/rho(0);kinetic=kinetic+.5d0*dot_product(p,velocity(:,0))
    nt=0;cs2=0
    do j=1,layout%nc
       k=layout%base+j;b=layout%owner(j)
       if(u(k)<0)return
       if(rho(b)==0.and.u(k)/=0)return
       if(layout%gamma_nt(j)>0)then
          nt=nt+u(k)
          pressure=pressure+(layout%gamma_nt(j)-1)*u(k)
          cs2=cs2+layout%gamma_nt(j)*(layout%gamma_nt(j)-1)*u(k)
       endif
    enddo
    thermal=u(5)-kinetic-nt
    if(thermal<0)return ! Never turn unresolved drift energy into a temperature floor.
    cs2=(cs2+gamma*(gamma-1)*thermal)/rho(0)
    pressure=pressure+(gamma-1)*thermal
    if(.not.all(ieee_is_finite([rho,reshape(velocity,[3*(layout%nb+1)]),thermal,kinetic,pressure,cs2])))return
    if(cs2<0)return
    sound=sqrt(cs2);ierr=dust_fv_ok
  end subroutine

  subroutine physical_flux(layout,u,axis,gamma,flux,speed,gas_velocity,ierr)
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::u(:),gamma
    integer,intent(in)::axis
    real(real64),intent(out)::flux(:),speed,gas_velocity
    integer,intent(out)::ierr
    real(real64)::rho(0:layout%nb),v(3,0:layout%nb),thermal,kinetic,p,cs,kg,kd,nt
    integer::b,k,j
    flux=0;speed=0;gas_velocity=0;ierr=dust_fv_config
    if(axis<1.or.axis>3.or.size(flux)/=size(u))return
    call dust_fv_decode(layout,u,gamma,rho,v,thermal,kinetic,p,cs,ierr)
    if(ierr/=0)return
    gas_velocity=v(axis,0);speed=max(abs(gas_velocity)+cs,maxval(abs(v(axis,1:))))
    flux(1)=u(axis+1)
    flux(2:4)=rho(0)*v(:,0)*v(axis,0);flux(axis+1)=flux(axis+1)+p
    kg=.5d0*rho(0)*sum(v(:,0)**2);nt=0
    do j=1,layout%nc
       if(layout%gamma_nt(j)>0)nt=nt+u(layout%base+j)
    enddo
    flux(5)=(thermal+nt+kg+p)*v(axis,0)
    do b=1,layout%nb
       k=6+4*(b-1)
       flux(k)=u(k+axis)
       flux(k+1:k+3)=u(k+1:k+3)*v(axis,b)
       flux(2:4)=flux(2:4)+flux(k+1:k+3)
       kd=.5d0*rho(b)*sum(v(:,b)**2)
       flux(5)=flux(5)+kd*v(axis,b)
    enddo
    do j=1,layout%nc
       flux(layout%base+j)=u(layout%base+j)*v(axis,layout%owner(j))
    enddo
    ierr=dust_fv_state
    if(any(.not.ieee_is_finite(flux)).or..not.ieee_is_finite(speed))return
    ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_face(layout,left,right,axis,gamma,flux,speed,gas_velocity,ierr)
    ! First-order shared Rusanov flux. One wave bound for gas AND all grains
    ! preserves consistency of mass, elemental and population face sums.
    ! It does not confer asymptotic-preserving strong-drag wave accuracy.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::left(:),right(:),gamma
    integer,intent(in)::axis
    real(real64),intent(inout)::flux(:),speed,gas_velocity
    integer,intent(out)::ierr
    real(real64)::fl(size(left)),fr(size(left)),work(size(left)),sl,sr,vl,vr,a
    ierr=dust_fv_config
    if(size(right)/=size(left).or.size(flux)/=size(left))return
    call physical_flux(layout,left,axis,gamma,fl,sl,vl,ierr)
    if(ierr/=0)return
    call physical_flux(layout,right,axis,gamma,fr,sr,vr,ierr)
    if(ierr/=0)return
    a=max(sl,sr);work=.5d0*fl+.5d0*fr-.5d0*a*(right-left)
    ierr=dust_fv_state
    if(any(.not.ieee_is_finite(work)))return
    flux=work;speed=a;gas_velocity=.5d0*vl+.5d0*vr;ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_periodic_step(layout,state,neighbors,dx,dt,gamma,next,ierr)
    ! Bounded native periodic FV driver using the same public face operator.
    ! No AMR/MPI/ghost or RAMSES state is accessed. All cells commit together.
    ! Per-cell sum(dt*a_face/dx)<=1 is a conservative multidimensional CFL.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::state(:,:),dx,dt,gamma
    integer,intent(in)::neighbors(:,:)
    real(real64),intent(inout)::next(:,:)
    integer,intent(out)::ierr
    real(real64),allocatable::trial(:,:),divv(:),cfl(:)
    real(real64)::flux(size(state,1)),a,vg,rho(0:layout%nb),v(3,0:layout%nb),th,kin,p,cs,f
    integer::n,i,j,axis,k,b
    ierr=dust_fv_config;n=size(state,2)
    if(.not.layout%ready.or.n<1.or.size(state,1)/=layout%nvar)return
    if(any(shape(next)/=shape(state)).or.any(shape(neighbors)/=[6,n]))return
    if(.not.all(ieee_is_finite([dx,dt])).or.dx<=0.or.dt<0)return
    if(any(neighbors<1).or.any(neighbors>n))return
    do i=1,n
       do axis=1,3
          j=neighbors(2*axis,i)
          if(neighbors(2*axis-1,j)/=i)return
          j=neighbors(2*axis-1,i)
          if(neighbors(2*axis,j)/=i)return
       enddo
    enddo
    trial=state;allocate(divv(n),cfl(n));divv=0;cfl=0;f=dt/dx
    do i=1,n
       do axis=1,3
          j=neighbors(2*axis,i);flux=0;a=0;vg=0
          call dust_fv_face(layout,state(:,i),state(:,j),axis,gamma,flux,a,vg,ierr)
          if(ierr/=0)return
          cfl(i)=cfl(i)+f*a;cfl(j)=cfl(j)+f*a
          if(j==i)cycle
          trial(:,i)=trial(:,i)-f*flux;trial(:,j)=trial(:,j)+f*flux
          divv(i)=divv(i)+vg/dx;divv(j)=divv(j)-vg/dx
       enddo
    enddo
    ierr=dust_fv_cfl
    if(any(.not.ieee_is_finite(cfl)).or.any(cfl>1d0))return
    do i=1,n
       do b=1,layout%nc
          if(layout%gamma_nt(b)==0)cycle
          k=layout%base+b
          ! Gas nonthermal pressure work. Total E is already conservative;
          ! its partition changes, never add this work to total E a second time.
          trial(k,i)=trial(k,i)-dt*(layout%gamma_nt(b)-1)*state(k,i)*divv(i)
       enddo
       call dust_fv_decode(layout,trial(:,i),gamma,rho,v,th,kin,p,cs,ierr)
       if(ierr/=0)return
    enddo
    next=trial;ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_drag(layout,state,gamma,stopping_time,dt,next,gas_heat,ierr)
    ! Connect spatially conserved component momenta to the existing implicit
    ! drag receiver. E and total p stay fixed; lost drift kinetic energy is
    ! therefore gas heat, not a second addition to the total-energy field.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::state(:),gamma,stopping_time(:),dt
    real(real64),intent(inout)::next(:),gas_heat
    integer,intent(out)::ierr
    real(real64)::rho(0:layout%nb),v(3,0:layout%nb),pd(3,layout%nb),pg(3),nd(3,layout%nb),ng(3)
    real(real64)::th,kin,p,cs,heat,work(size(state)),th_after,kin_after
    integer::b,k
    ierr=dust_fv_config
    if(size(next)/=size(state).or.size(stopping_time)/=layout%nb)return
    call dust_fv_decode(layout,state,gamma,rho,v,th,kin,p,cs,ierr)
    if(ierr/=0)return
    pg=state(2:4)
    do b=1,layout%nb
       k=6+4*(b-1);pd(:,b)=state(k+1:k+3);pg=pg-pd(:,b)
    enddo
    call dust_drag_step(rho(0),rho(1:),pg,pd,stopping_time,dt,ng,nd,heat,ierr)
    if(ierr/=0)return
    work=state
    do b=1,layout%nb
       k=6+4*(b-1);work(k+1:k+3)=nd(:,b)
    enddo
    call dust_fv_decode(layout,work,gamma,rho,v,th_after,kin_after,p,cs,ierr)
    if(ierr/=0)return
    ierr=dust_fv_state
    if(abs((th_after-th)-heat)>1024*epsilon(1d0)*max(abs(state(5)),heat,tiny(1d0)))return
    next=work;gas_heat=heat;ierr=dust_fv_ok
  end subroutine

  subroutine dust_fv_periodic_advance(layout,state,neighbors,dx,dt,gamma,stopping_time,next,gas_heat,ierr)
    ! First-order transport followed by implicit local drag, fixed input t_s.
    ! A late drag failure discards the ENTIRE spatial + drag step. The caller
    ! must resolve stiff spatial accuracy; BE stability alone is not an AP proof.
    type(dust_fv_layout),intent(in)::layout
    real(real64),intent(in)::state(:,:),dx,dt,gamma,stopping_time(:,:)
    integer,intent(in)::neighbors(:,:)
    real(real64),intent(inout)::next(:,:),gas_heat(:)
    integer,intent(out)::ierr
    real(real64),allocatable::transported(:,:),trial(:,:),heat(:)
    integer::i,n
    ierr=dust_fv_config;n=size(state,2)
    if(any(shape(next)/=shape(state)).or.size(gas_heat)/=n)return
    if(any(shape(stopping_time)/=[layout%nb,n]))return
    transported=state;trial=state;allocate(heat(n));heat=0
    call dust_fv_periodic_step(layout,state,neighbors,dx,dt,gamma,transported,ierr)
    if(ierr/=0)return
    do i=1,n
       call dust_fv_drag(layout,transported(:,i),gamma,stopping_time(:,i),dt,trial(:,i),heat(i),ierr)
       if(ierr/=0)return
    enddo
    next=trial;gas_heat=heat;ierr=dust_fv_ok
  end subroutine
end module

! First-order CPU reference adapter for the actual RAMSES flux/reflux path.
! No terminal-velocity AP or kinetic multi-stream interpretation is implied.
module dust_dynamics_runtime
  use dust_phase_state
  use hydro_parameters
  use dust_multifluid
  use dust_mass_physics, only: olivine_fraction,dust_pah_hc,dust_pah_nstate,dust_pah_molecule_g, &
       dust_size_radius_cm,dust_size_density,dust_mp,dust_kb
  use dust_drag_physics, only: dust_epstein_stopping_time
  use dust_iron_optics, only: fe_radius_cm,fe_density
  use dust_composition_material, only: dust_composition_curve
  use dust_iron_compare, only: iron_compare_curve,iron_compare_temperature
  use snrt_dust_contract, only: snrt_dust_contract_number_temperature,snrt_dust_contract_temperature_k
  use snrt_dust_ir, only: snrt_dust_material_temperature
  use snrt_chimes, only: chimes_nuclear_sums,chimes_ns
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  real(dp),parameter::atomic_mass(11)=[1d0,4d0,12d0,14d0,16d0,20d0,24d0,28d0,32d0,40d0,56d0]
  public::dust_dynamics_face,dust_dynamics_eos,dust_dynamics_unsplit,dust_dynamics_drag
  public::dust_dynamics_stopping_times
contains
  subroutine dust_dynamics_stopping_times(row,sd,sv,st,stopping,ierr)
    ! Explicit neutral-collision Epstein comparison. No Coulomb/Lorentz drag.
    ! PAH uses the equal-volume carbon sphere at the configured C density;
    ! this is a geometric molecular drag approximation, not optical Qabs.
    real(dp),intent(in)::row(:),sd,sv,st
    real(dp),intent(out)::stopping(:)
    integer,intent(out)::ierr
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke,thermal,ntot,mpart,temp,mfp
    real(dp)::radius(ndust_phase),density(ndust_phase),drift(3,ndust_phase)
    integer::b
    ierr=1;stopping=0
    if(min(sd,sv,st,dust_drag_collision_cross_section_cm2)<=0)return
    call dust_phase_read(row,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    ierr=1;thermal=row(ndim+2)-ke
    if(nener>0)thermal=thermal-sum(row(inener:inener+nener-1))
    ntot=sum(row(ichimes:ichimes+chimes_ns-1))*sd/dust_mp
    if(ntot<=0.or.thermal<=0)return
    mpart=rg*sd/ntot;temp=(gamma-1)*thermal*sd*sv**2/(ntot*dust_kb)
    mfp=1d0/(ntot*dust_drag_collision_cross_section_cm2)
    radius(1:4)=[dust_size_radius_cm,dust_size_radius_cm]
    density(1:4)=[dust_size_density(1),dust_size_density(1),dust_size_density(2),dust_size_density(2)]
    if(idust_iron>0)then
       radius(5:6)=fe_radius_cm;density(5:6)=fe_density(1)
    endif
    if(idust_pah>0)then
       density(ndust_phase)=dust_size_density(1)
       radius(ndust_phase)=(3*dust_pah_molecule_g/(4*acos(-1d0)*density(ndust_phase)))**(1d0/3)
    endif
    drift=0
    do b=1,ndust_phase
       if(mass(b)>0)drift(:,b)=(pd(:,b)/mass(b)-pg/rg)*sv
    enddo
    call dust_epstein_stopping_time(rg*sd,temp,mpart,mfp,radius,density,drift,stopping,ierr)
    if(ierr==0)stopping=stopping/st
  end subroutine

  subroutine phase_layout(layout,ierr)
    type(dust_fv_layout),intent(out)::layout
    integer,intent(out)::ierr
    integer::owner(nvar+ndust_phase),b
    real(dp)::gn(nvar+ndust_phase)
    owner=0;gn=0
    if(nener>0)gn(inener:inener+nener-1)=gamma_rad(1:nener)
    if(idust_pah>0)owner(idust_pah:idust_pah+dust_pah_nstate()-1)=ndust_phase
    do b=1,ndust_phase
       owner(nvar+b)=b
    enddo
    call dust_fv_initialize(ndust_phase,owner,gn,layout,ierr)
  end subroutine

  subroutine split_enthalpy(row,mass,sv,enthalpy,ierr)
    real(dp),intent(in)::row(:),mass(:),sv
    real(dp),intent(out)::enthalpy(:)
    integer,intent(out)::ierr
    real(dp)::grains(2),fe,td,total,curve(1),unit(2),w
    real(dp)::mix(snrt_dust_contract_number_temperature),pure(snrt_dust_contract_number_temperature)
    integer::nt,b,m,j
    ierr=1;enthalpy=0
    total=sum(mass(1:4));fe=0
    if(idust_iron>0)fe=sum(mass(5:6))
    total=total+fe
    if(row(idust_energy)<0.or.sv<=0)return
    if(total==0)then
       if(row(idust_energy)==0)ierr=0
       return
    endif
    grains=[sum(mass(1:2)),sum(mass(3:4))]
    nt=snrt_dust_contract_number_temperature
    if(nt<2)return
    if(idust_iron>0.or.idust_pah>0)then
       call iron_compare_temperature(grains,fe,row(idust_energy)*sv**2,td,ierr)
    else
       call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt),grains,1d0,mix,ierr)
       if(ierr==0)call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
            mix,row(idust_energy)*sv**2/total,td,ierr)
    endif
    if(ierr/=0)return
    j=1
    do while(j<nt-1)
       if(td<=snrt_dust_contract_temperature_k(j+1))exit
       j=j+1
    enddo
    w=0
    if(idust_iron<=0.and.idust_pah<=0)w=log(td/snrt_dust_contract_temperature_k(j))/ &
         log(snrt_dust_contract_temperature_k(j+1)/snrt_dust_contract_temperature_k(j))
    do b=1,4+merge(2,0,idust_iron>0)
       if(mass(b)==0)cycle
       m=(b+1)/2;unit=0
       if(m<=2)unit(m)=1
       if(idust_iron>0.or.idust_pah>0)then
          call iron_compare_curve([td],unit,merge(1d0,0d0,m==3),1d0,curve,ierr)
       else
          call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt),unit,1d0,pure,ierr)
          curve(1)=(1-w)*pure(j)+w*pure(j+1)
       endif
       if(ierr/=0)return
       enthalpy(b)=mass(b)*curve(1)/sv**2
    enddo
    ierr=1
    if(abs(sum(enthalpy)-row(idust_energy))>1d-9*max(row(idust_energy),tiny(1d0)))return
    ! Remove roundoff only, retaining precisely the stored common-T total.
    b=maxloc(enthalpy,dim=1);enthalpy(b)=enthalpy(b)+row(idust_energy)-sum(enthalpy)
    ierr=0
  end subroutine

  subroutine pack(row,sv,state,ierr)
    real(dp),intent(in)::row(:),sv
    real(dp),intent(out)::state(:)
    integer,intent(out)::ierr
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke,solid(11),ed(ndust_phase)
    integer::b,k,base
    state=0;base=5+4*ndust_phase
    call dust_phase_read(row,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    call split_enthalpy(row,mass,sv,ed,ierr)
    if(ierr/=0)return
    state(1:5)=row(1:5);state(base+1:base+nvar)=row
    state(base+1:base+5)=0
    state(base+idust_momentum:base+idust_momentum+3*ndust_phase-1)=0
    state(base+idust)=0;state(base+idust_energy)=0
    state(base+idust_species:base+idust_species+1)=0
    state(base+idust_bins:base+idust_bins+3)=0
    if(idust_iron>0)state(base+idust_iron:base+idust_iron+1)=0
    solid=olivine_fraction*sum(mass(3:4));solid(3)=solid(3)+sum(mass(1:2))
    if(idust_iron>0)solid(11)=solid(11)+sum(mass(5:6))
    if(idust_pah>0)then
       solid(1)=solid(1)+mass(ndust_phase)*dust_pah_hc(1)
       solid(3)=solid(3)+mass(ndust_phase)*dust_pah_hc(2)
    endif
    ! Gas nuclei come from the actual positive gas chemistry state. Subtracting
    ! two nearly equal total/solid reservoirs generates negative roundoff at
    ! complete depletion. The nuclear projection is linear and needs no face
    ! normalization, including when gas and dust mass fluxes cancel.
    ierr=chimes_nuclear_sums(row(ichimes:ichimes+chimes_ns-1),state(base+ichem:base+ichem+10))
    if(ierr/=0)return
    state(base+ichem:base+ichem+10)=state(base+ichem:base+ichem+10)*atomic_mass
    ierr=1
    if(any(abs(row(ichem:ichem+10)-solid-state(base+ichem:base+ichem+10))> &
         1d-8*max(abs(row(ichem:ichem+10)),tiny(1d0))))return
    if(abs(sum(state(base+ichem:base+ichem+10))-rg)>1d-8*row(1))return
    state(base+imetal)=sum(state(base+ichem+2:base+ichem+10))
    do b=1,ndust_phase
       k=6+4*(b-1);state(k)=mass(b);state(k+1:k+3)=pd(:,b)
    enddo
    state(base+nvar+1:)=ed
    ierr=0
  end subroutine

  subroutine unpack(state,row)
    ! Linear reconstruction, valid for signed face fluxes and states alike.
    ! Never normalize a flux by its (possibly zero) total mass flux.
    real(dp),intent(in)::state(:)
    real(dp),intent(out)::row(:)
    real(dp)::mass(ndust_phase),solid(11)
    integer::b,k,base,j
    base=5+4*ndust_phase;row=state(base+1:base+nvar);row(1:5)=state(1:5)
    do b=1,ndust_phase
       k=6+4*(b-1);j=idust_momentum+3*(b-1)
       mass(b)=state(k);row(j:j+2)=state(k+1:k+3)
    enddo
    row(idust_bins:idust_bins+3)=mass(1:4)
    row(idust_species:idust_species+1)=[sum(mass(1:2)),sum(mass(3:4))]
    row(idust)=sum(mass(1:4))
    solid=olivine_fraction*sum(mass(3:4));solid(3)=solid(3)+sum(mass(1:2))
    if(idust_iron>0)then
       row(idust_iron:idust_iron+1)=mass(5:6)
       row(idust)=row(idust)+sum(mass(5:6));solid(11)=solid(11)+sum(mass(5:6))
    endif
    if(idust_pah>0)then
       solid(1)=solid(1)+mass(ndust_phase)*dust_pah_hc(1)
       solid(3)=solid(3)+mass(ndust_phase)*dust_pah_hc(2)
    endif
    row(ichem:ichem+10)=row(ichem:ichem+10)+solid
    row(imetal)=sum(row(ichem+2:ichem+10))
    row(idust_energy)=sum(state(base+nvar+1:))
  end subroutine

  subroutine dust_dynamics_face(left,right,axis,sv,flux,speed,gas_velocity,ierr)
    real(dp),intent(in)::left(:),right(:),sv
    integer,intent(in)::axis
    real(dp),intent(inout)::flux(:),speed,gas_velocity
    integer,intent(out)::ierr
    type(dust_fv_layout)::layout
    real(dp)::l(5+5*ndust_phase+nvar),r(5+5*ndust_phase+nvar),f(5+5*ndust_phase+nvar)
    call phase_layout(layout,ierr)
    if(ierr/=0)return
    call pack(left,sv,l,ierr)
    if(ierr/=0)return
    call pack(right,sv,r,ierr)
    if(ierr/=0)return
    f=0
    call dust_fv_face(layout,l,r,axis,gamma,f,speed,gas_velocity,ierr)
    if(ierr/=0)return
    call unpack(f,flux)
  end subroutine

  subroutine dust_dynamics_eos(row,thermal,kinetic,velocity,speed,ierr)
    real(dp),intent(in)::row(:)
    real(dp),intent(out)::thermal,kinetic,velocity(3),speed(3)
    integer,intent(out)::ierr
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),cs2
    integer::b,j
    call dust_phase_read(row,mass,pd,rg,pg,kinetic,ierr)
    if(ierr/=0)return
    ierr=1;thermal=row(ndim+2)-kinetic;cs2=0
    do j=1,nener
       if(row(inener+j-1)<0)return
       thermal=thermal-row(inener+j-1)
       cs2=cs2+gamma_rad(j)*(gamma_rad(j)-1)*row(inener+j-1)
    enddo
    if(thermal<0)return
    cs2=(cs2+gamma*(gamma-1)*thermal)/rg
    velocity=pg/rg;speed=abs(velocity)+max(sqrt(cs2),smallc)
    do b=1,ndust_phase
       if(mass(b)>0)speed=max(speed,abs(pd(:,b)/mass(b)))
    enddo
    if(any(.not.ieee_is_finite(speed)))return
    ierr=0
  end subroutine

  subroutine dust_dynamics_unsplit(uin,flux,tmp,dx,dy,dz,dt,ngrid)
    real(dp),intent(in)::uin(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar),dx,dy,dz,dt
    real(dp),intent(out)::flux(nvector,if1:if2,jf1:jf2,kf1:kf2,nvar,ndim)
    real(dp),intent(out)::tmp(nvector,if1:if2,jf1:jf2,kf1:kf2,2,ndim)
    integer,intent(in)::ngrid
    real(dp)::sl,st,sd,sv,snh,st2,spacing(3),f(nvar),a,vg
    integer::axis,i,j,k,l,lo(3),hi(3),offset(3),status
    call units(sl,st,sd,sv,snh,st2)
    spacing=[dx,dy,dz];flux=0;tmp=0
    do axis=1,ndim
       lo=[1,1,1];hi=[2,2,2];hi(axis)=3;offset=0;offset(axis)=1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                do l=1,ngrid
                   f=0;a=0;vg=0
                   call dust_dynamics_face(uin(l,i-offset(1),j-offset(2),k-offset(3),:), &
                        uin(l,i,j,k,:),axis,sv,f,a,vg,status)
                   if(status/=0)then
                      write(*,*)'ERROR: relative dust conserved face rejected',axis,i,j,k,status
                      call clean_stop;return
                   endif
                   flux(l,i,j,k,:,axis)=f*dt/spacing(axis)
                   tmp(l,i,j,k,1,axis)=vg*dt/spacing(axis)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine

  subroutine dust_dynamics_drag(row,stopping_time,dt,trial,gas_heat,ierr)
    ! Supplied single-grain stopping times in code time units. The physical
    ! closure (including its validity domain) is the caller's responsibility.
    real(dp),intent(in)::row(:),stopping_time(:),dt
    real(dp),intent(inout)::trial(:),gas_heat
    integer,intent(out)::ierr
    type(dust_fv_layout)::layout
    real(dp)::state(5+4*ndust_phase+nener),next(5+4*ndust_phase+nener)
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke,heat,work(size(row))
    integer::b,k
    ierr=1
    if(size(trial)/=size(row))return
    call dust_phase_read(row,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    call dust_fv_initialize(ndust_phase,spread(0,1,nener),gamma_rad(1:nener),layout,ierr)
    if(ierr/=0)return
    state(1:5)=row(1:5)
    do b=1,ndust_phase
       k=6+4*(b-1);state(k)=mass(b);state(k+1:k+3)=pd(:,b)
    enddo
    state(6+4*ndust_phase:)=row(inener:inener+nener-1);next=state;heat=0
    call dust_fv_drag(layout,state,gamma,stopping_time,dt,next,heat,ierr)
    if(ierr/=0)return
    work=row
    do b=1,ndust_phase
       k=idust_momentum+3*(b-1);work(k:k+2)=next(7+4*(b-1):9+4*(b-1))
    enddo
    trial=work;gas_heat=heat ! total p/E unchanged; heat already follows from KE loss
  end subroutine
end module dust_dynamics_runtime

subroutine dust_dynamics_sync_level(ilevel,dteff)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use dust_phase_state, only: dust_phase_gravity
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  real(dp),intent(in)::dteff
  integer::i,ind,cell,status
  real(dp)::old(nvar),trial(nvar)
  if(.not.dust_relative_motion)then
     call synchro_hydro_fine(ilevel,dteff)
     return
  endif
  if(.not.poisson.or.numbtot(1,ilevel)==0)return
!$omp parallel do private(i,ind,cell,status,old,trial)
  do i=1,active(ilevel)%ngrid
     do ind=1,twotondim
        cell=ICELL_OF(active(ilevel)%igrid(i),ind)
        old=uold(cell,:);trial=old
        call dust_phase_gravity(trial,old,f(cell,:),dteff,status)
        if(status/=0)call clean_stop
        uold(cell,:)=trial
     enddo
  enddo
end subroutine dust_dynamics_sync_level

subroutine dust_dynamics_advance_level(ilevel)
  use amr_commons
  use hydro_commons
  use dust_dynamics_runtime, only: dust_dynamics_stopping_times,dust_dynamics_drag
#include "amr_index.h"
  implicit none
  integer,intent(in)::ilevel
  integer::i,ind,cell,k,n,status,bad,all_bad,info,ivar
  integer,allocatable::cells(:)
  real(dp),allocatable::stage(:,:)
  real(dp)::sl,st,sd,sv,snh,st2,heat,stopping(ndust_phase),trial(nvar)
  include 'mpif.h'
  if(.not.dust_relative_motion)return
  call units(sl,st,sd,sv,snh,st2)
  n=active(ilevel)%ngrid*twotondim;k=0;bad=0
  allocate(cells(n),stage(3*ndust_phase,n))
  do ind=1,twotondim
     do i=1,active(ilevel)%ngrid
        cell=ICELL_OF(active(ilevel)%igrid(i),ind)
        if(son(cell)/=0)cycle
        k=k+1;cells(k)=cell
        call dust_dynamics_stopping_times(uold(cell,:),sd,sv,st,stopping,status)
        if(status/=0)then
           bad=1;cycle
        endif
        trial=uold(cell,:);heat=0
        call dust_dynamics_drag(uold(cell,:),stopping,dtnew(ilevel),trial,heat,status)
        if(status/=0)then
           bad=1;cycle
        endif
        stage(:,k)=trial(idust_momentum:idust_momentum+3*ndust_phase-1)
     enddo
  enddo
  call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  if(all_bad/=0.or.info/=0)then
     if(myid==1)write(*,*)'ERROR: relative dust drag/state/Epstein collision-domain rejection'
     call MPI_ABORT(MPI_COMM_WORLD,11,info)
  endif
  do i=1,k
     uold(cells(i),idust_momentum:idust_momentum+3*ndust_phase-1)=stage(:,i)
  enddo
  call upload_fine(ilevel)
  do ivar=idust_momentum,idust_momentum+3*ndust_phase-1
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  enddo
end subroutine dust_dynamics_advance_level

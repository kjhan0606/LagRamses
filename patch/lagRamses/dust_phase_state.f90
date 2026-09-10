! RAMSES conservative phase interpretation. No AMR access or live commit.
module dust_phase_state
  use hydro_parameters
  use dust_mass_physics, only: dust_pah_nstate
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: dust_dynamics_enabled,dust_phase_read,dust_phase_kick,dust_phase_comoving
  public :: dust_phase_masses,dust_phase_kinetic,dust_phase_gas_velocity,dust_phase_gravity
  public :: dust_phase_erode
contains
  subroutine dust_phase_erode(before,after,gas_momentum,dust_momentum,ierr)
    ! Explicit one-way solid->gas operator (shocks/sublimation only).
    ! Growth/destruction mixtures MUST use the directed rate-substep API.
    real(dp),intent(in)::before(:),after(:)
    real(dp),intent(inout)::gas_momentum(3),dust_momentum(:,:)
    integer,intent(out)::ierr
    real(dp)::trial(3,size(before)),pg(3),lost(3)
    integer::b
    ierr=1
    if(size(after)/=size(before).or.any(shape(dust_momentum)/=[3,size(before)]))return
    if(any(.not.ieee_is_finite(before)).or.any(.not.ieee_is_finite(after)))return
    if(any(.not.ieee_is_finite(dust_momentum)).or.any(.not.ieee_is_finite(gas_momentum)))return
    if(any(after<0).or.any(after>before))return
    trial=dust_momentum;pg=gas_momentum
    do b=1,size(before)
       if(before(b)==0)then
          if(any(trial(:,b)/=0))return
          cycle
       endif
       lost=(before(b)-after(b))/before(b)*dust_momentum(:,b)
       trial(:,b)=dust_momentum(:,b)-lost
       if(after(b)==0)trial(:,b)=0
       pg=pg+lost
    enddo
    gas_momentum=pg;dust_momentum=trial;ierr=0
  end subroutine

  logical function dust_dynamics_enabled()
    dust_dynamics_enabled=dust_relative_motion
  end function

  subroutine dust_phase_masses(row,mass,ierr)
    real(dp),intent(in)::row(:)
    real(dp),intent(out)::mass(:)
    integer,intent(out)::ierr
    integer::b
    ierr=1;mass=0
    if(.not.dust_relative_motion.or.ndim/=3.or.size(row)/=nvar)return
    if(ndust_phase<4.or.ndust_phase>7.or.size(mass)/=ndust_phase)return
    if(idust_bins<1.or.idust_bins+3>nvar)return
    mass(1:4)=row(idust_bins:idust_bins+3);b=4
    if(idust_iron>0)then
       if(idust_iron+1>nvar.or.b+2>ndust_phase)return
       mass(b+1:b+2)=row(idust_iron:idust_iron+1);b=b+2
    endif
    if(idust_pah>0)then
       if(idust_pah+dust_pah_nstate()-1>nvar.or.b+1>ndust_phase)return
       if(any(row(idust_pah:idust_pah+dust_pah_nstate()-1)<0))return
       b=b+1;mass(b)=sum(row(idust_pah:idust_pah+dust_pah_nstate()-1))
    endif
    if(b/=ndust_phase.or.any(.not.ieee_is_finite(mass)).or.any(mass<0))return
    ierr=0
  end subroutine

  subroutine dust_phase_read(row,mass,momentum,gas_mass,gas_momentum,kinetic,ierr)
    real(dp),intent(in)::row(:)
    real(dp),intent(out)::mass(:),momentum(:,:),gas_mass,gas_momentum(3),kinetic
    integer,intent(out)::ierr
    integer::b,k
    ierr=1;momentum=0;gas_mass=0;gas_momentum=0;kinetic=0
    if(any(shape(momentum)/=[3,ndust_phase]))return
    if(idust_momentum<1.or.idust_momentum+3*ndust_phase-1>size(row))return
    call dust_phase_masses(row,mass,ierr)
    if(ierr/=0)return
    ierr=1
    if(any(.not.ieee_is_finite(row)))return
    gas_mass=row(1)-sum(mass);gas_momentum=row(2:4)
    if(gas_mass<=0)return
    do b=1,ndust_phase
       k=idust_momentum+3*(b-1);momentum(:,b)=row(k:k+2)
       gas_momentum=gas_momentum-momentum(:,b)
       if(mass(b)==0)then
          if(any(momentum(:,b)/=0))return
       else
          kinetic=kinetic+.5d0*sum(momentum(:,b)**2)/mass(b)
       endif
    enddo
    kinetic=kinetic+.5d0*sum(gas_momentum**2)/gas_mass
    if(.not.ieee_is_finite(kinetic))return
    ierr=0
  end subroutine

  real(dp) function dust_phase_kinetic(row) result(kinetic)
    real(dp),intent(in)::row(:)
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3)
    integer::status
    kinetic=.5d0*sum(row(2:ndim+1)**2)/max(row(1),smallr)
    if(.not.dust_relative_motion)return
    call dust_phase_read(row,mass,pd,rg,pg,kinetic,status)
    if(status/=0)then
       write(*,*)'ERROR: invalid conservative dust phase state'
       call clean_stop
    endif
  end function

  function dust_phase_gas_velocity(row) result(v)
    real(dp),intent(in)::row(:)
    real(dp)::v(3),mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke
    integer::status
    v=0;v(1:ndim)=row(2:ndim+1)/max(row(1),smallr)
    if(.not.dust_relative_motion)return
    call dust_phase_read(row,mass,pd,rg,pg,ke,status)
    if(status/=0)then
       call clean_stop;return
    endif
    v=pg/rg
  end function

  subroutine dust_phase_kick(row,delta_gas,delta_dust,delta_energy,trial,ierr)
    real(dp),intent(in)::row(:),delta_gas(3),delta_dust(:,:),delta_energy
    real(dp),intent(inout)::trial(:)
    integer,intent(out)::ierr
    real(dp)::work(size(row)),mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke,thermal
    integer::b,k
    ierr=1
    if(size(trial)/=size(row).or.any(shape(delta_dust)/=[3,ndust_phase]))return
    if(any(.not.ieee_is_finite(delta_gas)).or.any(.not.ieee_is_finite(delta_dust)))return
    if(.not.ieee_is_finite(delta_energy))return
    call dust_phase_read(row,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    work=row;work(2:4)=row(2:4)+delta_gas+sum(delta_dust,dim=2)
    work(ndim+2)=row(ndim+2)+delta_energy
    do b=1,ndust_phase
       k=idust_momentum+3*(b-1);work(k:k+2)=row(k:k+2)+delta_dust(:,b)
    enddo
    call dust_phase_read(work,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    ierr=1;thermal=work(ndim+2)-ke
    if(nener>0)thermal=thermal-sum(work(inener:inener+nener-1))
    if(.not.ieee_is_finite(thermal).or.thermal<0)return
    trial=work;ierr=0
  end subroutine

  subroutine dust_phase_comoving(row,ierr)
    real(dp),intent(inout)::row(:)
    integer,intent(out)::ierr
    real(dp)::mass(ndust_phase),v(3)
    integer::b,k
    ierr=0
    if(.not.dust_relative_motion)return
    call dust_phase_masses(row,mass,ierr)
    if(ierr/=0)return
    ierr=1
    if(idust_momentum<1.or.idust_momentum+3*ndust_phase-1>size(row))return
    if(row(1)<=sum(mass).or.any(.not.ieee_is_finite(row(1:5))))return
    v=row(2:4)/row(1)
    do b=1,ndust_phase
       k=idust_momentum+3*(b-1);row(k:k+2)=mass(b)*v
    enddo
    ierr=0 ! Co-moving KE equals the pre-existing barycentric KE exactly.
  end subroutine

  subroutine dust_phase_gravity(row,old,acceleration,dt,ierr)
    real(dp),intent(inout)::row(:)
    real(dp),intent(in)::old(:),acceleration(3),dt
    integer,intent(out)::ierr
    real(dp)::mass(ndust_phase),pd(3,ndust_phase),rg,pg(3),ke
    real(dp)::mass0(ndust_phase),pd0(3,ndust_phase),rg0,pg0(3),ke0,work(size(row))
    integer::b,k
    call dust_phase_read(old,mass0,pd0,rg0,pg0,ke0,ierr)
    if(ierr/=0)return
    call dust_phase_read(row,mass,pd,rg,pg,ke,ierr)
    if(ierr/=0)return
    work=row;work(2:4)=row(2:4)+old(1)*acceleration*dt
    do b=1,ndust_phase
       k=idust_momentum+3*(b-1);work(k:k+2)=pd(:,b)+mass0(b)*acceleration*dt
    enddo
    call dust_phase_read(work,mass,pd,rg,pg,ke0,ierr)
    if(ierr/=0)return
    work(ndim+2)=row(ndim+2)+ke0-ke
    row=work
  end subroutine
end module dust_phase_state

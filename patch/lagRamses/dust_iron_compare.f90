! Live bounded comparison adapter. Independent Fe carriers, common T and
! electric/eddy opacity only. This is not full physical admission of Fe dust.
module dust_iron_compare
  use dust_mass_physics, only: dust_fe_max_temperature,dust_fe_max_primary_ev,dust_fe_condensation
  use dust_composition_material, only: dust_composition_curve,dust_composition_area,dl01_t
  use dust_iron_material
  use dust_iron_optics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
contains
  subroutine iron_compare_neutral_area(weights,normalization,neutral_fraction,area,ierr)
    ! Remove H+ from ONLY the Fe part of the old H-nucleus collision area:
    ! the UV catalytic operator already includes those ion impacts.
    real(real64),intent(in)::weights(6),normalization,neutral_fraction
    real(real64),intent(inout)::area
    integer,intent(out)::ierr
    real(real64)::next,cs_area
    ierr=1
    if(.not.all(ieee_is_finite([weights,normalization,neutral_fraction,area])))return
    if(any(weights<0).or.normalization<=0.or.neutral_fraction<0.or.neutral_fraction>1)return
    call dust_composition_area(weights(1:4),normalization,cs_area,ierr)
    if(ierr/=0)return
    ierr=1
    next=cs_area*sum(weights(1:4))+neutral_fraction*normalization*.75d0* &
         sum(weights(5:6)/(fe_density(1)*fe_radius_cm))
    if(next<0)return
    area=next;ierr=0
  end subroutine

  subroutine iron_compare_source_receipt(delta,old,primary,ierr)
    ! Signed catalytic sensible-energy receipt, not an invented photon heat.
    ! Feed a positive source to the existing implicit IR solve; any net loss
    ! debits only the funded pre-step material reservoir.
    real(real64),intent(in)::delta
    real(real64),intent(inout)::old,primary
    integer,intent(out)::ierr
    real(real64)::source
    ierr=1
    if(.not.all(ieee_is_finite([delta,old,primary])))return
    if(min(old,primary)<0)return
    source=primary+delta
    if(old+source<0)return
    if(source>=0)then
       primary=source
    else
       old=old+source;primary=0
    endif
    ierr=0
  end subroutine

  subroutine iron_compare_curve(nodes,grains,fe,normalization,curve,ierr)
    real(real64),intent(in)::nodes(:),grains(2),fe,normalization
    real(real64),intent(out)::curve(size(nodes))
    integer,intent(out)::ierr
    real(real64)::mass(3),lo,hi
    integer::k
    ierr=1;curve=0
    if(.not.all(ieee_is_finite([grains,fe,normalization])))return
    if(minval([grains,fe])<0.or.normalization<=0)return
    if(fe==0.and.all(nodes>=dl01_t(1)))then
       call dust_composition_curve(nodes,grains,normalization,curve,ierr)
       return
    endif
    mass=[.5d0,.5d0,0d0]*normalization
    if(sum([grains,fe])>0)mass=[grains,fe]/sum([grains,fe])*normalization
    do k=1,size(nodes)
       call iron_mixture_enthalpy(nodes(k),mass,lo,hi,ierr)
       if(ierr/=0)return
       curve(k)=lo
    enddo
    ierr=0
  end subroutine

  subroutine iron_compare_weights(bins,weights,area,normalization,ierr)
    real(real64),intent(in)::bins(6),normalization
    real(real64),intent(out)::weights(6),area
    integer,intent(out)::ierr
    real(real64)::cs_area,total
    weights=0;area=0;ierr=1
    if(any(.not.ieee_is_finite(bins)).or.any(bins<0))return
    if(.not.ieee_is_finite(normalization).or.normalization<=0)return
    total=sum(bins)
    ! Same harmless equal-C/S reference as the old dust-free cells; Fe=0.
    weights(1:4)=.25d0
    if(total>0)weights=bins/total
    call dust_composition_area(bins(1:4),normalization,cs_area,ierr)
    if(ierr/=0)return
    area=cs_area*sum(weights(1:4))+normalization*.75d0* &
         sum(weights(5:6)/(fe_density(1)*fe_radius_cm))
  end subroutine

  subroutine iron_compare_temperature(grains,fe,energy,t,ierr)
    real(real64),intent(in)::grains(2),fe,energy
    real(real64),intent(out)::t
    integer,intent(out)::ierr
    real(real64)::phase(4)
    call iron_mixture_state([grains,fe],energy,t,phase,ierr)
    if(ierr==0.and.fe>0.and.t>dust_fe_max_temperature)ierr=1
  end subroutine
end module

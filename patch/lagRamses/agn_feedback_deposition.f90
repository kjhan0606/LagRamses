! Native AGN cell coupling. All inputs use RAMSES code units; the deferred
! result is an integrated energy, not an energy density. No source calibration
! or temperature-floor cooling is performed here.
module agn_feedback_deposition
  use amr_parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: agn_deposit_cell, agn_jet_delta, agn_eddington_ratio
  public :: agn_jet_geometry, agn_contains_donor, agn_scalar_map
  public :: agn_withdraw_cell, agn_deposit_material, agn_pack_load, agn_unpack_load
  public :: agn_accretion_receipt, agn_accrete_scalars
  public :: agn_merge_pending
  public :: agn_partition_release, agn_reference_event, agn_reference_commit
  public :: agn_heat_ready, agn_legacy_energy, agn_reference_receiver_ready
  public :: agn_reference_jet_speed_ok, agn_reference_max_jet_speed_fraction
  integer, parameter, public :: agn_deposit_invalid_source=1, agn_deposit_invalid_receiver=2
  ! The reference profile is coupled to the Newtonian RAMSES velocity update.
  ! Do not silently enter a superluminal regime when retained mass is used as
  ! the jet loading entitlement.  A violating event is rejected before any
  ! hydro mutation; the reference model is not relativistically corrected.
  real(dp), parameter :: agn_reference_max_jet_speed_fraction=0.9d0

contains

  pure subroutine agn_partition_release(released,retained,chi,xfloor,receipt,ierr)
    ! receipt = EM erg, mechanical heat erg, jet erg, jet loading entitlement
    ! in retained code mass. These named shares belong to the reference model.
    real(dp),intent(in)::released,retained,chi,xfloor
    real(dp),intent(out)::receipt(4)
    integer,intent(out)::ierr
    real(dp),parameter::high_mechanical_share=0.15d0
    receipt=0d0; ierr=1
    if(.not.all(ieee_is_finite([released,retained,chi,xfloor])))return
    if(released<0d0.or.retained<0d0.or.chi<0d0.or.xfloor<=0d0)return
    if(chi>=xfloor)then
       receipt(2)=high_mechanical_share*released
       receipt(1)=released-receipt(2)
    else
       receipt(3)=released
       receipt(4)=retained
    endif
    ierr=0
  end subroutine agn_partition_release

  pure subroutine agn_reference_event(pending,bh_mass,jetfrac,energy,load_mass,is_jet,is_replay)
    ! pending = heat erg, jet erg, loading code mass, deferred erg.
    ! One event per sink/coarse step; heat's temperature trigger follows
    ! the existing averaging pass. Ineligible/unselected channels stay pending.
    real(dp),intent(in)::pending(4),bh_mass,jetfrac
    real(dp),intent(out)::energy,load_mass
    logical,intent(out)::is_jet,is_replay
    real(dp)::mass_fraction
    is_jet=.false.; is_replay=.false.; energy=0d0; load_mass=0d0
    if(.not.all(ieee_is_finite(pending)).or. &
         .not.all(ieee_is_finite([bh_mass,jetfrac])).or.any(pending<0d0).or. &
         bh_mass<=0d0.or.jetfrac<0d0)return
    if(pending(4)>0d0)then
       is_replay=.true.; energy=pending(4)
       return
    endif
    mass_fraction=0d0
    if(bh_mass>pending(3))mass_fraction=pending(3)/(bh_mass-pending(3))
    if(pending(2)>0d0.and.bh_mass>pending(3).and.mass_fraction>=jetfrac)then
       is_jet=.true.; energy=pending(2); load_mass=pending(3)
    else
       energy=pending(1)
    endif
  end subroutine agn_reference_event

  pure subroutine agn_reference_commit(pending,is_jet,is_replay,fired,deferred_erg,ierr)
    real(dp),intent(inout)::pending(4)
    logical,intent(in)::is_jet,is_replay,fired
    real(dp),intent(in)::deferred_erg
    integer,intent(out)::ierr
    ierr=1
    if(.not.all(ieee_is_finite(pending)).or.any(pending<0d0))return
    if(.not.ieee_is_finite(deferred_erg).or.deferred_erg<0d0)return
    if(is_jet.and.is_replay)return
    if(fired)then
       if(.not.is_replay)then
          if(is_jet)then
             pending(2:3)=0d0 ! Unloaded entitlement is consumed, not future gas.
          else
             pending(1)=0d0
          endif
       endif
       pending(4)=deferred_erg
    endif
    ierr=0
  end subroutine agn_reference_commit

  pure logical function agn_heat_ready(energy,mass,scale_t2,gamma,threshold) result(ready)
    real(dp),intent(in)::energy,mass,scale_t2,gamma,threshold
    ready=.false.
    if(mass>0d0.and.energy>0d0)ready=(gamma-1d0)*(energy/mass)*scale_t2>threshold
  end function agn_heat_ready

  pure logical function agn_reference_receiver_ready(is_replay,selected,vol_gas,ind_blast) result(ready)
    logical,intent(in)::is_replay,selected
    real(dp),intent(in)::vol_gas
    integer,intent(in)::ind_blast
    ready=.false.
    if(.not.selected.or..not.ieee_is_finite(vol_gas))return
    if(is_replay)then
       ! Replay is deposited over the thermal bubble and does not require the
       ! single-cell donor owner used by a fresh jet/heat event.
       ready=vol_gas>0d0
    else
       ready=ind_blast>0
    endif
  end function agn_reference_receiver_ready

  pure logical function agn_reference_jet_speed_ok(energy,mass,f_ek,scale_v) result(ok)
    real(dp),intent(in)::energy,mass,f_ek,scale_v
    real(dp),parameter::c_cgs=2.99792458d10
    real(dp)::speed_code
    ok=.false.
    if(.not.all(ieee_is_finite([energy,mass,f_ek,scale_v])))return
    if(energy<0d0.or.mass<=0d0.or.f_ek<0d0.or.f_ek>1d0.or.scale_v<=0d0)return
    speed_code=sqrt(2d0*f_ek*energy/mass)
    ok=ieee_is_finite(speed_code).and.speed_code<= &
         agn_reference_max_jet_speed_fraction*c_cgs/scale_v
  end function agn_reference_jet_speed_ok

  pure real(dp) function agn_legacy_energy(is_jet,epsilon_r,spin,mad,retained,eK,eT,scale_v) result(energy)
    logical,intent(in)::is_jet,mad
    real(dp),intent(in)::epsilon_r,spin,retained,eK,eT,scale_v
    real(dp)::eff_mad
    if(is_jet)then
       if(mad)then
          eff_mad=(4.10507+0.328712*spin+76.0849*spin**2d0 &
               +47.9235*spin**3d0+3.86634*spin**4d0)/100d0
          energy=eff_mad*retained*(3d10/scale_v)**2d0
       else
          energy=eK*epsilon_r*retained*(3d10/scale_v)**2d0
       endif
    else
       energy=eT*epsilon_r*retained*(3d10/scale_v)**2d0
    endif
  end function agn_legacy_energy

  pure subroutine agn_merge_pending(pending,groups,merged,ierr)
    real(dp),intent(in)::pending(:)
    integer,intent(in)::groups(:)
    real(dp),intent(out)::merged(:)
    integer,intent(out)::ierr
    integer::k
    merged=0d0; ierr=1
    if(size(groups)/=size(pending))return
    if(any(groups<1).or.any(groups>size(merged)))return
    if(.not.all(ieee_is_finite(pending)).or.any(pending<0d0))return
    do k=1,size(groups)
       merged(groups(k))=merged(groups(k))+pending(k)
    enddo
    if(.not.all(ieee_is_finite(merged)))then
       merged=0d0
       return
    endif
    ierr=0
  end subroutine agn_merge_pending

  pure subroutine agn_accretion_receipt(rho,rho_initial,volume,requested,floor_fraction, &
       epsilon,mass_unit,gross,retained,radiated_erg,ierr)
    real(dp), intent(in) :: rho,rho_initial,volume,requested,floor_fraction,epsilon,mass_unit
    real(dp), intent(out) :: gross,retained,radiated_erg
    integer, intent(out) :: ierr
    real(dp), parameter :: c_cgs=2.99792458d10
    real(dp) :: floor_density
    gross=0d0; retained=0d0; radiated_erg=0d0; ierr=1
    if(.not.all(ieee_is_finite([rho,rho_initial,volume,requested,floor_fraction,epsilon,mass_unit])))return
    if(rho<=0d0.or.rho_initial<0d0.or.volume<=0d0.or.requested<0d0.or.mass_unit<=0d0)return
    if(floor_fraction<=0d0.or.floor_fraction>1d0.or.epsilon<0d0.or.epsilon>=1d0)return
    ! A zero unew reference can occur in an unrefreshed reception cell.
    ! Use the current donor for a conservative per-event floor, never the
    ! zero reference that would allow emptying it. The caller counts this.
    floor_density=rho_initial
    if(floor_density==0d0)floor_density=rho
    gross=max(min(requested,(rho-floor_fraction*floor_density)*volume),0d0)
    retained=(1d0-epsilon)*gross
    radiated_erg=epsilon*gross*mass_unit*c_cgs**2
    if(.not.all(ieee_is_finite([gross,retained,radiated_erg])))then
       gross=0d0; retained=0d0; radiated_erg=0d0
       return
    endif
    ierr=0
  end subroutine agn_accretion_receipt

  pure subroutine agn_accrete_scalars(row,fields,metal_slot,gross,volume,ierr,hydro_last)
    ! Only the declared constituent densities change here. The caller retains
    ! the accretion hydro/MHD energy convention, not cold jet loading.
    real(dp), intent(inout) :: row(:)
    integer, intent(in) :: fields(:),metal_slot
    real(dp), intent(in) :: gross,volume
    integer, intent(out) :: ierr
    integer, optional, intent(in) :: hydro_last
    real(dp) :: factor
    integer :: k,last
    ierr=1
    last=5
    if(present(hydro_last))last=hydro_last
    if(last<3.or.size(row)<last)return
    if(.not.all(ieee_is_finite([row(1),gross,volume])))return
    if(row(1)<=0d0.or.gross<0d0.or.volume<=0d0)return
    if(any(fields<=last).or.any(fields>size(row)))return
    if(metal_slot<0.or.metal_slot>size(fields))return
    do k=1,size(fields)
       if(any(fields(:k-1)==fields(k)))return
    enddo
    if(.not.all(ieee_is_finite(row(fields))))return
    factor=1d0-gross/(row(1)*volume)
    if(.not.ieee_is_finite(factor).or.factor<0d0.or.factor>1d0)return
    ! Finite but unphysical advected composition is not corrected here:
    ! status 2 requests a whole-event skip, not a whole-run fatal error.
    ierr=2
    if(any(row(fields)<0d0))return
    if(metal_slot>0)then
       if(row(fields(metal_slot))>row(1))return
    endif
    row(fields)=row(fields)*factor
    ierr=0
  end subroutine agn_accrete_scalars

  pure real(dp) function agn_eddington_ratio(bondi, eddington) result(ratio)
    real(dp), intent(in) :: bondi, eddington
    ! The parent rejects negative/nonfinite inputs. An idle zero cap has no
    ! accretion, and must not generate a 0/0 in the branch-selection loops.
    ratio=0d0
    if (eddington>0d0) ratio=bondi/eddington
  end function agn_eddington_ratio

  pure subroutine agn_deposit_cell(row, density_delta, momentum_delta, energy_delta, &
       volume, gamma, scale_t2, temperature_cap, deferred_energy, ierr)
    real(dp), intent(inout) :: row(5)
    real(dp), intent(in) :: density_delta, momentum_delta(3), energy_delta
    real(dp), intent(in) :: volume, gamma, scale_t2, temperature_cap
    real(dp), intent(out) :: deferred_energy
    integer, intent(out) :: ierr
    real(dp) :: staged(5), kinetic_old, kinetic_new, kinetic_input
    real(dp) :: internal_old, internal_trial, internal_limit, trial_energy, tol

    ierr = agn_deposit_invalid_receiver
    deferred_energy = 0d0
    if (.not. all(ieee_is_finite(row)) .or. row(1)<=0d0) return
    ierr = agn_deposit_invalid_source
    if (.not. all(ieee_is_finite(momentum_delta)) .or. &
         .not. all(ieee_is_finite([density_delta, energy_delta, volume, gamma, &
         scale_t2, temperature_cap]))) return
    if (density_delta < 0d0 .or. energy_delta < 0d0 .or. &
         volume <= 0d0 .or. gamma <= 1d0 .or. scale_t2 <= 0d0 .or. temperature_cap <= 0d0) return

    staged = row
    staged(1) = row(1) + density_delta
    staged(2:4) = row(2:4) + momentum_delta
    if (.not. all(ieee_is_finite(staged(1:4)))) return
    kinetic_old = 0.5d0 * sum((row(2:4)/sqrt(row(1)))**2)
    kinetic_new = 0.5d0 * sum((staged(2:4)/sqrt(staged(1)))**2)
    kinetic_input = 0d0
    if (density_delta > 0d0) then
       kinetic_input = 0.5d0 * sum((momentum_delta/sqrt(density_delta))**2)
    else if (any(momentum_delta /= 0d0)) then
       return
    end if
    internal_old = row(5) - kinetic_old
    trial_energy = row(5) + energy_delta
    internal_trial = trial_energy - kinetic_new
    if (.not. all(ieee_is_finite([kinetic_old, kinetic_new, kinetic_input, &
         internal_old, trial_energy, internal_trial]))) return
    tol = 64d0 * epsilon(1d0) * max(tiny(1d0), abs(row(5)), abs(trial_energy), kinetic_input)
    if (internal_old < -tol) then
       ierr = agn_deposit_invalid_receiver
       return
    end if
    if (energy_delta < kinetic_input-tol .or. &
         internal_trial < max(0d0,internal_old)-tol) return

    ! Limit only newly added internal energy: never cool an already-hot cell.
    internal_limit = max(internal_old, &
         (temperature_cap/scale_t2)/(gamma-1d0)*staged(1))
    if (.not. ieee_is_finite(internal_limit)) return
    staged(5) = kinetic_new + min(internal_trial, internal_limit)
    deferred_energy = max(0d0, trial_energy-staged(5))*volume
    if (.not. ieee_is_finite(staged(5)) .or. .not. ieee_is_finite(deferred_energy)) then
       deferred_energy = 0d0
       return
    end if
    row = staged
    ierr = 0
  end subroutine agn_deposit_cell

  pure subroutine agn_jet_geometry(offset, spin, radius, weights, axis)
    real(dp), intent(in) :: offset(3), spin(3), radius
    real(dp), intent(out) :: weights(2), axis(3)
    real(dp) :: norm, axial, radial, weight
    weights=0d0; axis=0d0
    norm=sqrt(sum(spin**2))
    if(norm<=0d0 .or. radius<=0d0)return
    axis=spin/norm
    axial=sum(offset*axis)
    radial=sqrt(max(0d0,sum(offset**2)-axial**2))
    if(radial>radius .or. abs(axial)>radius)return
    weight=exp(-0.5d0*(radial/radius)**2)
    if(axial>0d0)then
       weights=[weight,0d0]
    else if(axial<0d0)then
       weights=[0d0,weight]
    else
       weights=0.5d0*weight
    endif
  end subroutine agn_jet_geometry

  pure logical function agn_contains_donor(offset, width) result(inside)
    real(dp), intent(in) :: offset(3), width
    ! offset = cell centre - sink; half-open physical cell avoids face ties.
    inside=all(offset>-0.5d0*width .and. offset<=0.5d0*width)
  end function agn_contains_donor

  pure subroutine agn_jet_delta(loaded_mass, weights, volume_weight_sum, &
       bulk_velocity, axis, jet_speed, density_delta, momentum_delta, kinetic_delta)
    ! Each lobe integrates to half the loaded mass, even on unequal meshes.
    ! Weights and their volume sums come from the SAME geometry helper.
    real(dp), intent(in) :: loaded_mass, weights(2), volume_weight_sum(2)
    real(dp), intent(in) :: bulk_velocity(3), axis(3), jet_speed
    real(dp), intent(out) :: density_delta, momentum_delta(3), kinetic_delta
    real(dp) :: drho(2), plus(3), minus(3)
    ! Unsupported lobes are routed to fallback by the caller before ANY
    ! receiver deposition; never divide by a missing lobe here.
    density_delta=0d0; momentum_delta=0d0; kinetic_delta=0d0
    if(any(volume_weight_sum<=0d0))return
    drho=0.5d0*loaded_mass*(weights/volume_weight_sum)
    plus=bulk_velocity+jet_speed*axis
    minus=bulk_velocity-jet_speed*axis
    density_delta=sum(drho)
    momentum_delta=drho(1)*plus+drho(2)*minus
    kinetic_delta=0.5d0*(drho(1)*sum(plus**2)+drho(2)*sum(minus**2))
  end subroutine agn_jet_delta

  pure subroutine agn_scalar_map(nvars, metal_index, first_element, nelements, reserved, fields, ierr,hydro_last)
    integer, intent(in) :: nvars, metal_index, first_element, nelements, reserved(:)
    integer, intent(out) :: fields(:), ierr
    integer, optional, intent(in) :: hydro_last
    integer :: k, nmetal,last
    fields=0; ierr=agn_deposit_invalid_source
    last=5
    if(present(hydro_last))last=hydro_last
    if(last<3.or.last>nvars)return
    nmetal=merge(1,0,metal_index/=0)
    if(nelements<0 .or. size(fields)/=nmetal+nelements)return
    if(nmetal==1)fields(1)=metal_index
    do k=1,nelements
       fields(nmetal+k)=first_element+k-1
    enddo
    do k=1,size(fields)
       if(fields(k)<=last .or. fields(k)>nvars)return
       if(any(fields(k)==reserved))return
       if(any(fields(k)==fields(:k-1)))return
    enddo
    ierr=0
  end subroutine agn_scalar_map

  pure subroutine agn_withdraw_cell(row, fields, metal_slot, requested_mass, volume, &
       loaded_mass, velocity, fractions, ierr)
    real(dp), intent(inout) :: row(:)
    integer, intent(in) :: fields(:), metal_slot
    real(dp), intent(in) :: requested_mass, volume
    real(dp), intent(out) :: loaded_mass, velocity(3), fractions(:)
    integer, intent(out) :: ierr
    real(dp) :: staged(size(row)), rho, kinetic, internal, ratio, tol
    loaded_mass=0d0; velocity=0d0; fractions=0d0
    ierr=agn_deposit_invalid_source
    if(size(row)<5 .or. size(fields)/=size(fractions))return
    if(any(fields<=5) .or. any(fields>size(row)))return
    if(metal_slot<0 .or. metal_slot>size(fields))return
    if(.not.all(ieee_is_finite([requested_mass,volume])))return
    if(requested_mass<0d0 .or. volume<=0d0)return
    ierr=agn_deposit_invalid_receiver
    if(.not.all(ieee_is_finite(row(1:5))) .or. row(1)<=0d0)return
    if(.not.all(ieee_is_finite(row(fields))) .or. any(row(fields)<0d0))return
    rho=row(1)
    if(metal_slot>0)then
       if(row(fields(metal_slot))>rho)return
    endif
    velocity=row(2:4)/rho
    fractions=row(fields)/rho
    kinetic=0.5d0*sum((row(2:4)/sqrt(rho))**2)
    internal=row(5)-kinetic
    tol=64d0*epsilon(1d0)*max(tiny(1d0),abs(row(5)),kinetic)
    if(.not.all(ieee_is_finite([velocity,fractions,kinetic,internal])))return
    if(internal < -tol)return
    loaded_mass=min(requested_mass,0.25d0*rho*volume)
    ratio=(loaded_mass/volume)/rho
    staged=row
    staged(1:4)=row(1:4)*(1d0-ratio)
    staged(5)=internal+kinetic*(1d0-ratio)
    staged(fields)=row(fields)*(1d0-ratio)
    if(.not.all(ieee_is_finite(staged(1:5))) .or. &
         .not.all(ieee_is_finite(staged(fields))))return
    row=staged
    ierr=0
  end subroutine agn_withdraw_cell

  pure subroutine agn_deposit_material(row, fields, metal_slot, fractions, drho, momentum, energy, &
       volume, gamma, scale_t2, cap, deferred, ierr)
    real(dp), intent(inout) :: row(:)
    integer, intent(in) :: fields(:), metal_slot
    real(dp), intent(in) :: fractions(:), drho, momentum(3), energy, volume, gamma, scale_t2, cap
    real(dp), intent(out) :: deferred
    integer, intent(out) :: ierr
    real(dp) :: staged(size(row))
    deferred=0d0; ierr=agn_deposit_invalid_source
    if(size(row)<5 .or. size(fields)/=size(fractions))return
    if(any(fields<=5) .or. any(fields>size(row)))return
    if(metal_slot<0 .or. metal_slot>size(fields))return
    staged=row
    if(drho>0d0)then
       if(.not.all(ieee_is_finite(fractions)) .or. any(fractions<0d0))return
       if(metal_slot>0)then
          if(fractions(metal_slot)>1d0)return
       endif
       ierr=agn_deposit_invalid_receiver
       if(.not.all(ieee_is_finite(row(fields))) .or. any(row(fields)<0d0))return
       if(metal_slot>0)then
          if(row(fields(metal_slot))>row(1))return
       endif
       ierr=agn_deposit_invalid_source
       staged(fields)=row(fields)+drho*fractions
       if(.not.all(ieee_is_finite(staged(fields))))return
    endif
    call agn_deposit_cell(staged(1:5),drho,momentum,energy,volume,gamma,scale_t2,cap,deferred,ierr)
    if(ierr==0)row=staged
  end subroutine agn_deposit_material

  pure subroutine agn_pack_load(mass, velocity, fractions, packed)
    real(dp), intent(in) :: mass, velocity(3), fractions(:)
    real(dp), intent(out) :: packed(4+size(fractions))
    packed(1)=mass; packed(2:4)=velocity; packed(5:)=fractions
  end subroutine agn_pack_load

  pure subroutine agn_unpack_load(packed, mass, velocity, fractions)
    real(dp), intent(in) :: packed(:)
    real(dp), intent(out) :: mass, velocity(3), fractions(size(packed)-4)
    mass=packed(1); velocity=packed(2:4); fractions=packed(5:)
  end subroutine agn_unpack_load
end module agn_feedback_deposition

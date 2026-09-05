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
  integer, parameter, public :: agn_deposit_invalid_source=1, agn_deposit_invalid_receiver=2

contains

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

  pure subroutine agn_scalar_map(nvars, metal_index, first_element, nelements, reserved, fields, ierr)
    integer, intent(in) :: nvars, metal_index, first_element, nelements, reserved(:)
    integer, intent(out) :: fields(:), ierr
    integer :: k, nmetal
    fields=0; ierr=agn_deposit_invalid_source
    nmetal=merge(1,0,metal_index/=0)
    if(nelements<0 .or. size(fields)/=nmetal+nelements)return
    if(nmetal==1)fields(1)=metal_index
    do k=1,nelements
       fields(nmetal+k)=first_element+k-1
    enddo
    do k=1,size(fields)
       if(fields(k)<=5 .or. fields(k)>nvars)return
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

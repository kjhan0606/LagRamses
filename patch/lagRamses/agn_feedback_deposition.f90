! Native AGN cell coupling. All inputs use RAMSES code units; the deferred
! result is an integrated energy, not an energy density. No source calibration
! or temperature-floor cooling is performed here.
module agn_feedback_deposition
  use amr_parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: agn_deposit_cell, agn_jet_delta, agn_eddington_ratio
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

  pure subroutine agn_jet_delta(loaded_mass, weight, volume_weight_sum, &
       bulk_velocity, axis, jet_speed, axial_distance, density_delta, momentum_delta, kinetic_delta)
    ! weight is the unnormalised Gaussian. The caller sums weight*cell_volume
    ! over the SAME cylinder membership, so sum(density_delta*volume)=loaded_mass.
    real(dp), intent(in) :: loaded_mass, weight, volume_weight_sum
    real(dp), intent(in) :: bulk_velocity(3), axis(3), jet_speed, axial_distance
    real(dp), intent(out) :: density_delta, momentum_delta(3), kinetic_delta
    real(dp) :: velocity(3), speed_squared

    density_delta = loaded_mass * (weight/volume_weight_sum)
    velocity = bulk_velocity
    if (axial_distance < 0d0) velocity = bulk_velocity - jet_speed*axis
    if (axial_distance > 0d0) velocity = bulk_velocity + jet_speed*axis
    speed_squared = sum(velocity**2)
    ! Exactly on the midplane: equal opposed lobe contributions, not a stale
    ! velocity or arbitrary lobe. Their net momentum is zero in the BH frame,
    ! but their second moment (kinetic energy) is not zero.
    if (axial_distance == 0d0) speed_squared = sum(bulk_velocity**2) + jet_speed**2
    momentum_delta = density_delta*velocity
    kinetic_delta = 0.5d0*density_delta*speed_squared
  end subroutine agn_jet_delta
end module agn_feedback_deposition

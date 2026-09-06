! DUST-6: primary-photon dust absorption and heating contract.
! Dust is a non-depleting absorber.  This module deliberately does not
! change the existing three-species CUDA ABI or any RAMSES state array.
module snrt_dust_coupling
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_thermochemistry, only: snrt_partition_absorption, &
       snrt_thermochemistry_ok
  implicit none
  private

  integer, parameter, public :: dust_coupling_dp = real64
  integer, parameter, public :: dust_coupling_ok = 0
  integer, parameter, public :: dust_coupling_err_input = 1
  integer, parameter, public :: dust_coupling_err_shape = 2
  integer, parameter, public :: dust_coupling_err_closure = 3

  real(real64), parameter :: dust_ev_to_erg = 1.602176634d-12
  real(real64), parameter :: dust_closure_ulps = 256.0d0

  public :: snrt_dust_prepare_optical_depth
  public :: snrt_dust_total_optical_depth
  public :: snrt_dust_partition_group
  public :: snrt_dust_heating_from_absorbed

contains

  subroutine snrt_dust_prepare_optical_depth(n_hydrogen_cm3, &
       dust_relative_abundance, path_cm, sigma_per_h_cm2, tau_dust, ierr)
    real(real64), intent(in) :: n_hydrogen_cm3
    real(real64), intent(in) :: dust_relative_abundance
    real(real64), intent(in) :: path_cm
    real(real64), intent(in) :: sigma_per_h_cm2
    real(real64), intent(out) :: tau_dust
    integer, intent(out) :: ierr

    tau_dust = 0.0d0
    ierr = dust_coupling_ok
    if (.not. ieee_is_finite(n_hydrogen_cm3) .or. &
         .not. ieee_is_finite(dust_relative_abundance) .or. &
         .not. ieee_is_finite(path_cm) .or. &
         .not. ieee_is_finite(sigma_per_h_cm2) .or. &
         n_hydrogen_cm3 < 0.0d0 .or. dust_relative_abundance < 0.0d0 .or. &
         path_cm < 0.0d0 .or. sigma_per_h_cm2 < 0.0d0) then
       ierr = dust_coupling_err_input
       return
    end if
    tau_dust = n_hydrogen_cm3 * dust_relative_abundance * path_cm * &
         sigma_per_h_cm2
    if (.not. ieee_is_finite(tau_dust) .or. tau_dust < 0.0d0) then
       tau_dust = 0.0d0
       ierr = dust_coupling_err_input
    end if
  end subroutine snrt_dust_prepare_optical_depth


  subroutine snrt_dust_total_optical_depth(tau_hhe, tau_dust, tau_total, ierr)
    real(real64), intent(in) :: tau_hhe(:)
    real(real64), intent(in) :: tau_dust
    real(real64), intent(out) :: tau_total
    integer, intent(out) :: ierr

    tau_total = 0.0d0
    ierr = dust_coupling_ok
    if (size(tau_hhe) /= 3) then
       ierr = dust_coupling_err_shape
       return
    end if
    if (.not. ieee_is_finite(tau_dust) .or. tau_dust < 0.0d0 .or. &
         any(.not. ieee_is_finite(tau_hhe)) .or. any(tau_hhe < 0.0d0)) then
       ierr = dust_coupling_err_input
       return
    end if
    tau_total = sum(tau_hhe) + tau_dust
    if (.not. ieee_is_finite(tau_total) .or. tau_total < 0.0d0) then
       tau_total = 0.0d0
       ierr = dust_coupling_err_input
    end if
  end subroutine snrt_dust_total_optical_depth


  subroutine snrt_dust_partition_group(raw_removed, tau_hhe_species, &
       tau_dust, available_hhe, assigned_hhe, assigned_dust, returned, &
       unassigned, ierr)
    ! All absorption quantities are homogeneous photon amounts for one cell
    ! and one group.  The caller restores `returned` to the photon field.
    real(real64), intent(in) :: raw_removed
    real(real64), intent(in) :: tau_hhe_species(3)
    real(real64), intent(in) :: tau_dust
    real(real64), intent(inout) :: available_hhe(3)
    real(real64), intent(out) :: assigned_hhe(3)
    real(real64), intent(out) :: assigned_dust
    real(real64), intent(out) :: returned
    real(real64), intent(out) :: unassigned
    integer, intent(out) :: ierr

    real(real64) :: tau_hhe, tau_total, hhe_target, hhe_target_capped
    real(real64) :: eligible_inventory, hhe_excess, dust_direct
    real(real64) :: dust_transfer, dust_fraction, residual, scale, tolerance
    real(real64) :: work_available(3), hhe_unassigned, inventory_scale
    integer :: hhe_ierr

    assigned_hhe = 0.0d0
    assigned_dust = 0.0d0
    returned = 0.0d0
    unassigned = 0.0d0
    ierr = dust_coupling_ok
    if (.not. ieee_is_finite(raw_removed) .or. raw_removed < 0.0d0 .or. &
         .not. ieee_is_finite(tau_dust) .or. tau_dust < 0.0d0 .or. &
         any(.not. ieee_is_finite(tau_hhe_species)) .or. &
         any(tau_hhe_species < 0.0d0) .or. &
         any(.not. ieee_is_finite(available_hhe)) .or. &
         any(available_hhe < 0.0d0)) then
       ierr = dust_coupling_err_input
       return
    end if

    tau_hhe = sum(tau_hhe_species)
    tau_total = tau_hhe + tau_dust
    if (.not. ieee_is_finite(tau_total) .or. tau_total <= 0.0d0) then
       ! A positive upstream removal with no optical depth is inconsistent;
       ! do not silently turn that inconsistency into physical heating.
       if (raw_removed > 0.0d0) ierr = dust_coupling_err_input
       return
    end if
    if (raw_removed == 0.0d0) return

    work_available = available_hhe
    inventory_scale = sum(work_available)
    if (tau_dust == 0.0d0) then
       ! Preserve the existing three-species path exactly in the zero-dust
       ! limit, including its caller-visible partition arithmetic.
       hhe_target = raw_removed
    else
       hhe_target = min(raw_removed, raw_removed * tau_hhe / tau_total)
    end if
    eligible_inventory = sum(pack(work_available, tau_hhe_species > 0.0d0))
    hhe_target_capped = min(hhe_target, eligible_inventory)
    hhe_unassigned = 0.0d0
    hhe_ierr = snrt_thermochemistry_ok
    if (hhe_target_capped > 0.0d0) then
       call snrt_partition_absorption(hhe_target_capped, tau_hhe_species, &
            work_available, assigned_hhe, hhe_ierr, hhe_unassigned, &
            inventory_scale)
       if (hhe_ierr /= snrt_thermochemistry_ok) then
          assigned_hhe = 0.0d0
          ierr = dust_coupling_err_closure
          return
       end if
    end if

    if (.not. ieee_is_finite(hhe_unassigned) .or. hhe_unassigned < 0.0d0 .or. &
         hhe_unassigned > hhe_target) then
       assigned_hhe = 0.0d0
       ierr = dust_coupling_err_closure
       return
    end if
    ! Preserve any numerical remainder reported by the H/He callee as an
    ! unassigned ledger term; only a physical inventory deficit is eligible
    ! for the dust-transfer approximation below.
    hhe_excess = max(0.0d0, hhe_target - sum(assigned_hhe) - hhe_unassigned)
    dust_direct = max(0.0d0, raw_removed - hhe_target)
    if (tau_dust > 0.0d0) then
       if (tau_dust < 1.0d-6) then
          dust_fraction = tau_dust * (1.0d0 - 0.5d0*tau_dust + &
               tau_dust*tau_dust/6.0d0)
       else
          dust_fraction = 1.0d0 - exp(-tau_dust)
       end if
    else
       dust_fraction = 0.0d0
    end if
    dust_transfer = hhe_excess * dust_fraction
    assigned_dust = dust_direct + dust_transfer
    returned = hhe_excess - dust_transfer

    residual = raw_removed - sum(assigned_hhe) - assigned_dust - returned - &
         hhe_unassigned
    scale = max(raw_removed, sum(abs(assigned_hhe)) + abs(assigned_dust) + &
         abs(returned) + abs(hhe_unassigned))
    tolerance = dust_closure_ulps * epsilon(1.0d0) * scale
    if (.not. ieee_is_finite(assigned_dust) .or. &
         .not. ieee_is_finite(returned) .or. assigned_dust < 0.0d0 .or. &
         returned < 0.0d0 .or. abs(residual) > tolerance) then
       assigned_hhe = 0.0d0
       assigned_dust = 0.0d0
       returned = 0.0d0
       ierr = dust_coupling_err_closure
       return
    end if
    unassigned = hhe_unassigned + max(0.0d0, residual)
    available_hhe = work_available
  end subroutine snrt_dust_partition_group


  subroutine snrt_dust_heating_from_absorbed(absorbed_photons_cm3, &
       mean_energy_ev, dt_s, heating_erg_cm3_s, ierr)
    ! absorbed_photons_cm3(group,cell) is already a physical photon density.
    ! No RAMSES code-unit conversion belongs in this generic ledger routine.
    real(real64), intent(in) :: absorbed_photons_cm3(:,:)
    real(real64), intent(in) :: mean_energy_ev(:)
    real(real64), intent(in) :: dt_s
    real(real64), intent(out) :: heating_erg_cm3_s(:)
    integer, intent(out) :: ierr
    integer :: group, cell

    heating_erg_cm3_s = 0.0d0
    ierr = dust_coupling_ok
    if (size(absorbed_photons_cm3,1) /= size(mean_energy_ev) .or. &
         size(absorbed_photons_cm3,2) /= size(heating_erg_cm3_s)) then
       ierr = dust_coupling_err_shape
       return
    end if
    if (.not. ieee_is_finite(dt_s) .or. dt_s <= 0.0d0 .or. &
         any(.not. ieee_is_finite(absorbed_photons_cm3)) .or. &
         any(absorbed_photons_cm3 < 0.0d0) .or. &
         any(.not. ieee_is_finite(mean_energy_ev)) .or. &
         any(mean_energy_ev < 0.0d0)) then
       ierr = dust_coupling_err_input
       return
    end if
    do cell = 1, size(heating_erg_cm3_s)
       do group = 1, size(mean_energy_ev)
          heating_erg_cm3_s(cell) = heating_erg_cm3_s(cell) + &
               absorbed_photons_cm3(group,cell) * mean_energy_ev(group)
       end do
       heating_erg_cm3_s(cell) = heating_erg_cm3_s(cell) * &
            dust_ev_to_erg / dt_s
    end do
    if (any(.not. ieee_is_finite(heating_erg_cm3_s))) then
       heating_erg_cm3_s = 0.0d0
       ierr = dust_coupling_err_input
    end if
  end subroutine snrt_dust_heating_from_absorbed

end module snrt_dust_coupling

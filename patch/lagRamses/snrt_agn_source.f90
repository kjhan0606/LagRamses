module snrt_agn_source
  ! Photon budgets for the S_N AGN source path.
  !
  ! delta_accreted_mass_code is an increment over the current radiation
  ! update, expressed in code mass units.  It must never be a sink seed
  ! mass or the accumulated BH mass.  Sink bookkeeping and AMR-cell
  ! deposition are intentionally kept outside this module.
  use, intrinsic :: iso_c_binding, only: c_float
  use amr_parameters, only: dp
  implicit none

  real(dp), parameter :: snrt_c_cgs = 2.99792458d10
  real(dp), parameter :: snrt_ev_to_erg = 1.602176634d-12

contains

  subroutine snrt_agn_photon_budget(delta_accreted_mass_code, mass_unit_g, &
       delta_t_s, radiative_efficiency, ionizing_fraction, &
       mean_photon_energy_ev, luminosity_erg_s, emitted_photons)
    real(dp), intent(in) :: delta_accreted_mass_code, mass_unit_g, delta_t_s
    real(dp), intent(in) :: radiative_efficiency, ionizing_fraction
    real(dp), intent(in) :: mean_photon_energy_ev
    real(dp), intent(out) :: luminosity_erg_s, emitted_photons
    real(dp) :: radiated_energy_erg

    luminosity_erg_s = 0.0d0
    emitted_photons = 0.0d0

    if (delta_accreted_mass_code <= 0.0d0) return
    if (mass_unit_g <= 0.0d0 .or. delta_t_s <= 0.0d0) return
    if (radiative_efficiency <= 0.0d0 .or. ionizing_fraction <= 0.0d0) return
    if (mean_photon_energy_ev <= 0.0d0) return

    radiated_energy_erg = radiative_efficiency * delta_accreted_mass_code * &
         mass_unit_g * snrt_c_cgs**2
    luminosity_erg_s = radiated_energy_erg / delta_t_s
    emitted_photons = ionizing_fraction * radiated_energy_erg / &
         (mean_photon_energy_ev * snrt_ev_to_erg)
  end subroutine snrt_agn_photon_budget

  subroutine snrt_agn_isotropic_packet(emitted_photons, angular_weights, &
       directional_photons)
    ! Split an angle-integrated photon count using the S_N quadrature.
    ! The caller retains emitted_photons as the conservation ledger; this
    ! routine only provides its directional representation for transport.
    real(dp), intent(in) :: emitted_photons
    real(dp), intent(in) :: angular_weights(:)
    real(dp), intent(out) :: directional_photons(size(angular_weights))
    real(dp) :: total_weight

    directional_photons = 0.0d0
    if (emitted_photons <= 0.0d0) return

    total_weight = sum(max(angular_weights, 0.0d0))
    if (total_weight <= 0.0d0) return

    directional_photons = emitted_photons * max(angular_weights, 0.0d0) / &
         total_weight
  end subroutine snrt_agn_isotropic_packet

  subroutine snrt_agn_photons_to_density_code(emitted_photons, &
       cell_volume_code, length_unit_cm, n_h_unit_cm3, photon_density_code)
    ! Convert a cell-integrated physical photon count to the S_N state
    ! variable n_gamma / n_H,unit.  This keeps the c_float transport state
    ! in a usable dynamic range while leaving the integral photon ledger in
    ! double precision at the source boundary.
    real(dp), intent(in) :: emitted_photons, cell_volume_code
    real(dp), intent(in) :: length_unit_cm, n_h_unit_cm3
    real(dp), intent(out) :: photon_density_code

    photon_density_code = 0.0d0
    if (emitted_photons <= 0.0d0) return
    if (cell_volume_code <= 0.0d0 .or. length_unit_cm <= 0.0d0) return
    if (n_h_unit_cm3 <= 0.0d0) return

    photon_density_code = emitted_photons / (cell_volume_code * &
         length_unit_cm**3 * n_h_unit_cm3)
  end subroutine snrt_agn_photons_to_density_code

  subroutine snrt_agn_deposit_isotropic(state, slot, group, emitted_photons, &
       cell_volume_code, length_unit_cm, n_h_unit_cm3, angular_weights, &
       deposited_density_code, ierr)
    ! State layout is (direction, group, slot).  The source is supplied as
    ! an integrated physical photon count and is converted at this boundary
    ! to n_gamma / n_H,unit before its angular split.
    real(c_float), intent(inout) :: state(:,:,:)
    integer, intent(in) :: slot, group
    real(dp), intent(in) :: emitted_photons, cell_volume_code
    real(dp), intent(in) :: length_unit_cm, n_h_unit_cm3
    real(dp), intent(in) :: angular_weights(:)
    real(dp), intent(out) :: deposited_density_code
    integer, intent(out) :: ierr
    real(dp), allocatable :: directional_density(:)
    real(dp) :: state_limit

    ierr = 0
    deposited_density_code = 0.0d0
    if (slot < 1 .or. slot > size(state,3)) then
       ierr = 1
       return
    end if
    if (group < 1 .or. group > size(state,2)) then
       ierr = 2
       return
    end if
    if (size(state,1) /= size(angular_weights)) then
       ierr = 3
       return
    end if

    call snrt_agn_photons_to_density_code(emitted_photons, cell_volume_code, &
         length_unit_cm, n_h_unit_cm3, deposited_density_code)
    if (deposited_density_code <= 0.0d0) return

    allocate(directional_density(size(angular_weights)))
    call snrt_agn_isotropic_packet(deposited_density_code, angular_weights, &
         directional_density)
    state_limit = real(huge(0.0_c_float), dp) * 0.5d0
    if (maxval(abs(real(state(:,group,slot),dp)) + directional_density) > &
         state_limit) then
       ierr = 4
       return
    end if
    state(:,group,slot) = state(:,group,slot) + real(directional_density,c_float)
  end subroutine snrt_agn_deposit_isotropic

end module snrt_agn_source

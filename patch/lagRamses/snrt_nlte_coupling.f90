module snrt_nlte_coupling
  ! Unit conversions at the boundary between the S_N photon state and an
  ! external NLTE chemistry solver.  All code densities are normalized by
  ! n_H,unit; optical depth and heating are returned in physical units.
  use amr_parameters, only: dp
  use snrt_agn_source, only: snrt_c_cgs, snrt_ev_to_erg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

contains

  subroutine snrt_nlte_optical_depth(neutral_hydrogen_code, n_h_unit_cm3, &
       delta_t_s, cross_section_cm2, optical_depth, light_speed_factor)
    real(dp), intent(in) :: neutral_hydrogen_code, n_h_unit_cm3
    real(dp), intent(in) :: delta_t_s, cross_section_cm2
    real(dp), intent(out) :: optical_depth
    real(dp), intent(in), optional :: light_speed_factor
    real(dp) :: speed_factor

    optical_depth = 0.0d0
    if (neutral_hydrogen_code <= 0.0d0 .or. n_h_unit_cm3 <= 0.0d0) return
    if (delta_t_s <= 0.0d0 .or. cross_section_cm2 <= 0.0d0) return
    speed_factor = 1.0d0
    if (present(light_speed_factor)) speed_factor = max(0.0d0, &
         min(1.0d0, light_speed_factor))
    if (speed_factor <= 0.0d0) return
    optical_depth = snrt_c_cgs * delta_t_s * neutral_hydrogen_code * &
         n_h_unit_cm3 * cross_section_cm2 * speed_factor
  end subroutine snrt_nlte_optical_depth

  subroutine snrt_nlte_optical_depth_groups(neutral_hydrogen_code, &
       n_h_unit_cm3, delta_t_s, cross_section_cm2, optical_depth, ierr, &
       light_speed_factor)
    ! Construct all group optical depths from one neutral-H state.  The
    ! cross-section table is owned by the NLTE/chemistry layer and is not
    ! hard-coded into the transport operator.
    real(dp), intent(in) :: neutral_hydrogen_code, n_h_unit_cm3
    real(dp), intent(in) :: delta_t_s
    real(dp), intent(in) :: cross_section_cm2(:)
    real(dp), intent(out) :: optical_depth(:)
    integer, intent(out) :: ierr
    real(dp), intent(in), optional :: light_speed_factor
    real(dp) :: column_factor, speed_factor

    ierr = 0
    optical_depth = 0.0d0
    if (size(optical_depth) /= size(cross_section_cm2)) then
       ierr = 1
       return
    end if
    if (neutral_hydrogen_code <= 0.0d0 .or. n_h_unit_cm3 <= 0.0d0 .or. &
         delta_t_s <= 0.0d0) return

    speed_factor = 1.0d0
    if (present(light_speed_factor)) speed_factor = max(0.0d0, &
         min(1.0d0, light_speed_factor))
    if (speed_factor <= 0.0d0) return
    column_factor = snrt_c_cgs * delta_t_s * neutral_hydrogen_code * &
         n_h_unit_cm3 * speed_factor
    optical_depth = max(0.0d0, cross_section_cm2) * column_factor
  end subroutine snrt_nlte_optical_depth_groups

  subroutine snrt_nlte_primordial_optical_depth_groups( &
       neutral_hydrogen_code, neutral_helium_i_code, neutral_helium_ii_code, &
       n_h_unit_cm3, delta_t_s, cross_section_hydrogen, cross_section_helium_i, &
       cross_section_helium_ii, optical_depth_total, optical_depth_hydrogen, &
       optical_depth_helium_i, optical_depth_helium_ii, ierr, light_speed_factor)
    ! Construct total and species-resolved primordial optical depths.  The
    ! three input species densities are in the same code density convention as
    ! snrt_intensity (number density divided by n_H,unit).  Keeping the
    ! component depths lets the host partition the CUDA-reported absorption
    ! without re-evaluating an inconsistent opacity after transport.
    real(dp), intent(in) :: neutral_hydrogen_code, neutral_helium_i_code
    real(dp), intent(in) :: neutral_helium_ii_code, n_h_unit_cm3, delta_t_s
    real(dp), intent(in) :: cross_section_hydrogen(:), cross_section_helium_i(:)
    real(dp), intent(in) :: cross_section_helium_ii(:)
    real(dp), intent(out) :: optical_depth_total(:), optical_depth_hydrogen(:)
    real(dp), intent(out) :: optical_depth_helium_i(:), optical_depth_helium_ii(:)
    integer, intent(out) :: ierr
    real(dp), intent(in), optional :: light_speed_factor
    real(dp) :: column_factor, speed_factor

    ierr = 0
    optical_depth_total = 0.0d0
    optical_depth_hydrogen = 0.0d0
    optical_depth_helium_i = 0.0d0
    optical_depth_helium_ii = 0.0d0
    if (size(optical_depth_total) /= size(cross_section_hydrogen) .or. &
         size(optical_depth_total) /= size(cross_section_helium_i) .or. &
         size(optical_depth_total) /= size(cross_section_helium_ii) .or. &
         size(optical_depth_total) /= size(optical_depth_hydrogen) .or. &
         size(optical_depth_total) /= size(optical_depth_helium_i) .or. &
         size(optical_depth_total) /= size(optical_depth_helium_ii)) then
       ierr = 1
       return
    end if
    if (.not. ieee_is_finite(neutral_hydrogen_code) .or. &
         .not. ieee_is_finite(neutral_helium_i_code) .or. &
         .not. ieee_is_finite(neutral_helium_ii_code) .or. &
         .not. ieee_is_finite(n_h_unit_cm3) .or. .not. ieee_is_finite(delta_t_s) .or. &
         neutral_hydrogen_code < 0.0d0 .or. neutral_helium_i_code < 0.0d0 .or. &
         neutral_helium_ii_code < 0.0d0 .or. n_h_unit_cm3 <= 0.0d0 .or. &
         delta_t_s <= 0.0d0 .or. any(.not. ieee_is_finite(cross_section_hydrogen)) .or. &
         any(.not. ieee_is_finite(cross_section_helium_i)) .or. &
         any(.not. ieee_is_finite(cross_section_helium_ii)) .or. &
         any(cross_section_hydrogen < 0.0d0) .or. any(cross_section_helium_i < 0.0d0) .or. &
         any(cross_section_helium_ii < 0.0d0)) then
       ierr = 2
       return
    end if
    speed_factor = 1.0d0
    if (present(light_speed_factor)) speed_factor = max(0.0d0, &
         min(1.0d0, light_speed_factor))
    if (speed_factor <= 0.0d0) return
    column_factor = snrt_c_cgs * delta_t_s * n_h_unit_cm3 * speed_factor
    optical_depth_hydrogen = cross_section_hydrogen * &
         neutral_hydrogen_code * column_factor
    optical_depth_helium_i = cross_section_helium_i * &
         neutral_helium_i_code * column_factor
    optical_depth_helium_ii = cross_section_helium_ii * &
         neutral_helium_ii_code * column_factor
    optical_depth_total = optical_depth_hydrogen + optical_depth_helium_i + &
         optical_depth_helium_ii
  end subroutine snrt_nlte_primordial_optical_depth_groups

  subroutine snrt_nlte_photo_source(absorbed_photon_code, &
       total_hydrogen_code, ionized_fraction, n_h_unit_cm3, delta_t_s, &
       photon_energy_ev, ionization_increment, heating_rate_erg_cm3_s, &
       photoelectron_excess_energy_ev)
    real(dp), intent(in) :: absorbed_photon_code, total_hydrogen_code
    real(dp), intent(in) :: ionized_fraction, n_h_unit_cm3, delta_t_s
    real(dp), intent(in) :: photon_energy_ev
    real(dp), intent(out) :: ionization_increment, heating_rate_erg_cm3_s
    real(dp), intent(in), optional :: photoelectron_excess_energy_ev
    real(dp) :: available_neutral_code, used_photon_code, excess_energy_erg
    real(dp) :: excess_energy_ev

    ionization_increment = 0.0d0
    heating_rate_erg_cm3_s = 0.0d0
    if (.not. ieee_is_finite(absorbed_photon_code) .or. &
         .not. ieee_is_finite(total_hydrogen_code) .or. &
         .not. ieee_is_finite(ionized_fraction) .or. &
         .not. ieee_is_finite(n_h_unit_cm3) .or. &
         .not. ieee_is_finite(delta_t_s) .or. &
         .not. ieee_is_finite(photon_energy_ev)) return
    if (absorbed_photon_code <= 0.0d0 .or. total_hydrogen_code <= 0.0d0) return
    if (n_h_unit_cm3 <= 0.0d0 .or. delta_t_s <= 0.0d0) return
    if (photon_energy_ev <= 0.0d0) return

    available_neutral_code = total_hydrogen_code * max(0.0d0, &
         min(1.0d0, 1.0d0 - ionized_fraction))
    used_photon_code = min(absorbed_photon_code, available_neutral_code)
    if (present(photoelectron_excess_energy_ev)) then
       if (.not. ieee_is_finite(photoelectron_excess_energy_ev) .or. &
            photoelectron_excess_energy_ev < 0.0d0) return
       excess_energy_ev = photoelectron_excess_energy_ev
    else
       excess_energy_ev = max(0.0d0, photon_energy_ev - 13.6d0)
    end if
    ionization_increment = used_photon_code / total_hydrogen_code
    excess_energy_erg = excess_energy_ev * snrt_ev_to_erg
    heating_rate_erg_cm3_s = used_photon_code * n_h_unit_cm3 * &
         excess_energy_erg / delta_t_s
  end subroutine snrt_nlte_photo_source

end module snrt_nlte_coupling

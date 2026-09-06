! DUST-9: source-bound cell mapping and persistent thermal receiver.
!
! This module is deliberately caller-owned.  It does not access RAMSES
! hydro arrays, infer dust from legacy kappa_IR, or mutate a gas/radiation
! state.  A live adapter must supply a dedicated dust state and a validated
! source-bound opacity sidecar before calling these routines in production.
module snrt_dust_receiver
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter, public :: snrt_dust_receiver_dp = real64
  integer, parameter, public :: snrt_dust_receiver_ok = 0
  integer, parameter, public :: snrt_dust_receiver_err_shape = 1
  integer, parameter, public :: snrt_dust_receiver_err_input = 2
  integer, parameter, public :: snrt_dust_receiver_err_binding = 3
  integer, parameter, public :: snrt_dust_receiver_err_state = 4
  integer, parameter, public :: snrt_dust_receiver_err_closure = 5

  real(real64), parameter, public :: snrt_dust_receiver_ev_to_erg = &
       1.602176634d-12
  real(real64), parameter :: receiver_closure_ulps = 256.0d0

  public :: snrt_dust_validate_opacity_binding
  public :: snrt_dust_map_cell_abundance
  public :: snrt_dust_prepare_cell_optical_depth
  public :: snrt_dust_receiver_stage
  public :: snrt_dust_receiver_commit

contains

  subroutine snrt_dust_validate_opacity_binding(group_edges_ev, sigma_per_h_cm2, &
       mean_energy_ev, binding_status, source_id, source_sha256, edges_sha256, ierr)
    ! The upstream JSON reader owns file I/O and hash calculation.  This
    ! native boundary rechecks the numerical payload and the identity tokens
    ! before any cell optical depth is constructed.
    real(real64), intent(in) :: group_edges_ev(:), sigma_per_h_cm2(:)
    real(real64), intent(in) :: mean_energy_ev(:)
    character(len=*), intent(in) :: binding_status, source_id
    character(len=*), intent(in) :: source_sha256, edges_sha256
    integer, intent(out) :: ierr
    integer :: group, ng

    ierr = snrt_dust_receiver_ok
    ng = size(sigma_per_h_cm2)
    if (size(group_edges_ev) /= ng + 1 .or. size(mean_energy_ev) /= ng .or. ng < 1) then
       ierr = snrt_dust_receiver_err_shape
       return
    end if
    if (.not. all(ieee_is_finite(group_edges_ev)) .or. &
         .not. all(ieee_is_finite(sigma_per_h_cm2)) .or. &
         .not. all(ieee_is_finite(mean_energy_ev))) then
       ierr = snrt_dust_receiver_err_input
       return
    end if
    if (any(group_edges_ev <= 0.0d0) .or. &
         any(group_edges_ev(2:) <= group_edges_ev(:ng)) .or. &
         any(sigma_per_h_cm2 < 0.0d0) .or. any(mean_energy_ev <= 0.0d0)) then
       ierr = snrt_dust_receiver_err_binding
       return
    end if
    do group = 1, ng
       ! Closed endpoint admission is intentional.  The native transport
       ! group convention remains owned by the spectral contract; this check
       ! only prevents an out-of-band representative energy.
       if (mean_energy_ev(group) < group_edges_ev(group) .or. &
            mean_energy_ev(group) > group_edges_ev(group + 1)) then
          ierr = snrt_dust_receiver_err_binding
          return
       end if
    end do
    if (trim(binding_status) /= 'candidate_source_sed_matched' .and. &
         trim(binding_status) /= 'candidate_scattering_isotropic' .and. &
         trim(binding_status) /= 'reference_control' .and. &
         trim(binding_status) /= 'reference_scattering_control') then
       ierr = snrt_dust_receiver_err_binding
       return
    end if
    if (len_trim(source_id) == 0 .or. .not. is_sha256(source_sha256) .or. &
         .not. is_sha256(edges_sha256)) ierr = snrt_dust_receiver_err_binding
  end subroutine snrt_dust_validate_opacity_binding


  subroutine snrt_dust_map_cell_abundance(metallicity_solar, dust_to_metal, &
       dust_relative_abundance, ierr)
    ! Both inputs are explicit dimensionless cell fields.  The result is the
    ! factor multiplying the reference-mixture opacity per H nucleus.
    real(real64), intent(in) :: metallicity_solar(:), dust_to_metal(:)
    real(real64), intent(out) :: dust_relative_abundance(:)
    integer, intent(out) :: ierr

    ierr = snrt_dust_receiver_ok
    dust_relative_abundance = 0.0d0
    if (size(metallicity_solar) /= size(dust_to_metal) .or. &
         size(dust_relative_abundance) /= size(metallicity_solar)) then
       ierr = snrt_dust_receiver_err_shape
       return
    end if
    if (.not. all(ieee_is_finite(metallicity_solar)) .or. &
         .not. all(ieee_is_finite(dust_to_metal)) .or. &
         any(metallicity_solar < 0.0d0) .or. any(dust_to_metal < 0.0d0)) then
       ierr = snrt_dust_receiver_err_input
       return
    end if
    dust_relative_abundance = metallicity_solar * dust_to_metal
    if (.not. all(ieee_is_finite(dust_relative_abundance)) .or. &
         any(dust_relative_abundance < 0.0d0)) then
       dust_relative_abundance = 0.0d0
       ierr = snrt_dust_receiver_err_input
    end if
  end subroutine snrt_dust_map_cell_abundance


  subroutine snrt_dust_prepare_cell_optical_depth(n_hydrogen_cm3, path_cm, &
       dust_relative_abundance, sigma_per_h_cm2, tau_dust, ierr)
    real(real64), intent(in) :: n_hydrogen_cm3(:), path_cm(:)
    real(real64), intent(in) :: dust_relative_abundance(:), sigma_per_h_cm2(:)
    real(real64), intent(out) :: tau_dust(:,:)
    integer, intent(out) :: ierr
    integer :: cell, group, nc, ng

    ierr = snrt_dust_receiver_ok
    tau_dust = 0.0d0
    nc = size(n_hydrogen_cm3)
    ng = size(sigma_per_h_cm2)
    if (size(path_cm) /= nc .or. size(dust_relative_abundance) /= nc .or. &
         any(shape(tau_dust) /= [nc, ng]) .or. nc < 1 .or. ng < 1) then
       ierr = snrt_dust_receiver_err_shape
       return
    end if
    if (.not. all(ieee_is_finite(n_hydrogen_cm3)) .or. &
         .not. all(ieee_is_finite(path_cm)) .or. &
         .not. all(ieee_is_finite(dust_relative_abundance)) .or. &
         .not. all(ieee_is_finite(sigma_per_h_cm2)) .or. &
         any(n_hydrogen_cm3 < 0.0d0) .or. any(path_cm < 0.0d0) .or. &
         any(dust_relative_abundance < 0.0d0) .or. any(sigma_per_h_cm2 < 0.0d0)) then
       ierr = snrt_dust_receiver_err_input
       return
    end if
    do cell = 1, nc
       do group = 1, ng
          tau_dust(cell, group) = n_hydrogen_cm3(cell) * path_cm(cell) * &
               dust_relative_abundance(cell) * sigma_per_h_cm2(group)
       end do
    end do
    if (.not. all(ieee_is_finite(tau_dust)) .or. any(tau_dust < 0.0d0)) then
       tau_dust = 0.0d0
       ierr = snrt_dust_receiver_err_input
    end if
  end subroutine snrt_dust_prepare_cell_optical_depth


  subroutine snrt_dust_receiver_stage(absorbed_photons_cm3, mean_energy_ev, dt_s, &
       dust_relative_abundance, heat_capacity_erg_cm3_k, old_energy_erg_cm3, &
       old_temperature_k, staged_energy_erg_cm3, staged_temperature_k, &
       absorbed_energy_erg_cm3, ierr)
    ! This is a local, constant-capacity thermal step.  A future physical
    ! grain table can supply a temperature-dependent capacity at this same
    ! boundary.  The old arrays are never modified by this routine.
    real(real64), intent(in) :: absorbed_photons_cm3(:,:), mean_energy_ev(:)
    real(real64), intent(in) :: dt_s, dust_relative_abundance(:)
    real(real64), intent(in) :: heat_capacity_erg_cm3_k(:)
    real(real64), intent(in) :: old_energy_erg_cm3(:), old_temperature_k(:)
    real(real64), intent(out) :: staged_energy_erg_cm3(:), staged_temperature_k(:)
    real(real64), intent(out) :: absorbed_energy_erg_cm3(:)
    integer, intent(out) :: ierr
    integer :: cell, group, nc, ng
    real(real64) :: energy, residual, scale, tolerance

    ierr = snrt_dust_receiver_ok
    staged_energy_erg_cm3 = 0.0d0
    staged_temperature_k = 0.0d0
    absorbed_energy_erg_cm3 = 0.0d0
    nc = size(old_energy_erg_cm3)
    ng = size(mean_energy_ev)
    if (any(shape(absorbed_photons_cm3) /= [ng, nc]) .or. &
         size(dust_relative_abundance) /= nc .or. &
         size(heat_capacity_erg_cm3_k) /= nc .or. &
         size(old_temperature_k) /= nc .or. &
         size(staged_energy_erg_cm3) /= nc .or. &
         size(staged_temperature_k) /= nc .or. &
         size(absorbed_energy_erg_cm3) /= nc .or. nc < 1 .or. ng < 1) then
       ierr = snrt_dust_receiver_err_shape
       return
    end if
    if (.not. ieee_is_finite(dt_s) .or. dt_s <= 0.0d0 .or. &
         .not. all(ieee_is_finite(absorbed_photons_cm3)) .or. &
         .not. all(ieee_is_finite(mean_energy_ev)) .or. &
         .not. all(ieee_is_finite(dust_relative_abundance)) .or. &
         .not. all(ieee_is_finite(heat_capacity_erg_cm3_k)) .or. &
         .not. all(ieee_is_finite(old_energy_erg_cm3)) .or. &
         .not. all(ieee_is_finite(old_temperature_k)) .or. &
         any(absorbed_photons_cm3 < 0.0d0) .or. any(mean_energy_ev <= 0.0d0) .or. &
         any(dust_relative_abundance < 0.0d0) .or. &
         any(heat_capacity_erg_cm3_k <= 0.0d0) .or. &
         any(old_energy_erg_cm3 < 0.0d0) .or. any(old_temperature_k <= 0.0d0)) then
       ierr = snrt_dust_receiver_err_input
       return
    end if
    do cell = 1, nc
       energy = 0.0d0
       do group = 1, ng
          energy = energy + absorbed_photons_cm3(group, cell) * &
               mean_energy_ev(group)
       end do
       absorbed_energy_erg_cm3(cell) = energy * snrt_dust_receiver_ev_to_erg
       if (.not. ieee_is_finite(absorbed_energy_erg_cm3(cell)) .or. &
            absorbed_energy_erg_cm3(cell) < 0.0d0) then
          ierr = snrt_dust_receiver_err_input
          return
       end if
       if (dust_relative_abundance(cell) == 0.0d0 .and. &
            absorbed_energy_erg_cm3(cell) > 0.0d0) then
          ierr = snrt_dust_receiver_err_state
          return
       end if
       staged_energy_erg_cm3(cell) = old_energy_erg_cm3(cell) + &
            absorbed_energy_erg_cm3(cell)
       staged_temperature_k(cell) = old_temperature_k(cell) + &
            absorbed_energy_erg_cm3(cell) / heat_capacity_erg_cm3_k(cell)
       residual = staged_energy_erg_cm3(cell) - old_energy_erg_cm3(cell) - &
            absorbed_energy_erg_cm3(cell)
       scale = max(absorbed_energy_erg_cm3(cell), &
            abs(staged_energy_erg_cm3(cell)), 1.0d0)
       tolerance = receiver_closure_ulps * epsilon(1.0d0) * scale
       if (.not. ieee_is_finite(staged_energy_erg_cm3(cell)) .or. &
            .not. ieee_is_finite(staged_temperature_k(cell)) .or. &
            staged_energy_erg_cm3(cell) < 0.0d0 .or. &
            staged_temperature_k(cell) <= 0.0d0 .or. abs(residual) > tolerance) then
          ierr = snrt_dust_receiver_err_closure
          return
       end if
    end do
  end subroutine snrt_dust_receiver_stage


  subroutine snrt_dust_receiver_commit(state_energy_erg_cm3, state_temperature_k, &
       staged_energy_erg_cm3, staged_temperature_k, ierr)
    ! All validation precedes either assignment, so a rejected commit leaves
    ! both persistent state arrays byte-for-byte unchanged.
    real(real64), intent(inout) :: state_energy_erg_cm3(:), state_temperature_k(:)
    real(real64), intent(in) :: staged_energy_erg_cm3(:), staged_temperature_k(:)
    integer, intent(out) :: ierr

    ierr = snrt_dust_receiver_ok
    if (size(state_energy_erg_cm3) /= size(state_temperature_k) .or. &
         size(staged_energy_erg_cm3) /= size(state_energy_erg_cm3) .or. &
         size(staged_temperature_k) /= size(state_temperature_k)) then
       ierr = snrt_dust_receiver_err_shape
       return
    end if
    if (.not. all(ieee_is_finite(state_energy_erg_cm3)) .or. &
         .not. all(ieee_is_finite(state_temperature_k)) .or. &
         .not. all(ieee_is_finite(staged_energy_erg_cm3)) .or. &
         .not. all(ieee_is_finite(staged_temperature_k)) .or. &
         any(state_energy_erg_cm3 < 0.0d0) .or. &
         any(state_temperature_k <= 0.0d0) .or. &
         any(staged_energy_erg_cm3 < 0.0d0) .or. &
         any(staged_temperature_k <= 0.0d0)) then
       ierr = snrt_dust_receiver_err_state
       return
    end if
    state_energy_erg_cm3 = staged_energy_erg_cm3
    state_temperature_k = staged_temperature_k
  end subroutine snrt_dust_receiver_commit


  logical function is_sha256(value)
    character(len=*), intent(in) :: value
    integer :: i, code

    is_sha256 = .false.
    if (len_trim(value) /= 64) return
    do i = 1, 64
       code = iachar(value(i:i))
       if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
            (code >= iachar('a') .and. code <= iachar('f')) .or. &
            (code >= iachar('A') .and. code <= iachar('F')))) return
    end do
    is_sha256 = .true.
  end function is_sha256

end module snrt_dust_receiver

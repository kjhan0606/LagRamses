! Phase 0 runtime adapter for the RAMSES stellar-particle feedback loop.
!
! The legacy feedback routine interpolates a prompt table and overlaps the
! total-metal field with the first chemical field.  This adapter evaluates an
! age increment through the Phase 0 SSP driver and deposits total metal plus
! eleven independently tracked ejecta fields without that overlap.

module stellar_ramses_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_commons
  use pm_commons
  use hydro_commons
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       active_element, enable_pisn, &
       default_imf_id, population_model_id, yield_source_basis_id, &
       configured_imf_mass_min, configured_imf_mass_max, &
       configured_binary_fraction, configured_channel_mass_min, &
       configured_channel_mass_max, channel_owns_terminal_remnant, &
       production_source_model_supported
  use stellar_enrichment_contract, only: stellar_population_t, &
       stellar_source_t, delayed_cooling_source_mass, untracked_ejecta_mass, &
       generic_metal_ejecta_mass
  use stellar_yield_tables, only: stellar_yield_table_t, &
       set_yield_mass_assignment_mode, yield_mass_assignment_piecewise_constant, &
       yield_table_ok
  use stellar_yield_backend, only: load_yield_backend, backend_ok
  use stellar_yield_audit, only: audit_yield_table
  use stellar_enrichment_driver, only: compute_stellar_source_increment
  use stellar_population_ledger, only: stellar_population_ledger_t
  use stellar_native_units, only: code_time_to_age_gyr, &
       code_interval_to_age_gyr, mass_code_to_msun, mass_msun_to_code, &
       momentum_cgs_to_code, energy_erg_to_code, units_ok
  use stellar_progress_contract, only: stellar_progress_t, &
       progress_initialize, progress_begin, progress_commit, progress_abort, &
       progress_export, progress_ok
  use stellar_ramses_field_map, only: stellar_field_map_t, &
       clear_field_map, validate_field_map
#include "amr_index.h"
  implicit none
  private

  real(stellar_dp), parameter :: solar_mass_cgs = 1.98847d33
  real(stellar_dp), parameter :: source_tolerance = 1.0d-10
  integer, parameter :: phase0_mass_bins = 64

  type(stellar_yield_table_t), save :: yield_table
  type(stellar_field_map_t), save :: runtime_field_map
  logical, save :: initialized = .false.
  integer, save :: initialization_ierr = 0

  public :: phase0_initialize
  public :: phase0_feedback

contains

  subroutine phase0_initialize(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: status, table_ierr, audit_ierr, assignment_ierr, element
    logical :: exists

    if (initialized .or. initialization_ierr /= 0) then
       ierr = initialization_ierr
       return
    end if

    ierr = 0
    if (.not. production_source_model_supported()) then
       ierr = 3
       initialization_ierr = ierr
       if (myid == 1) write(*,*) &
            'Phase 0 source model is not implemented for production'
       return
    end if
    call get_environment_variable('PHASE0_YIELD_TABLE', filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) then
       ! Production startup must name the approved external table explicitly.
       ! The embedded backend is not a valid physical fallback.
       ierr = 1
       initialization_ierr = ierr
       if (myid == 1) write(*,*) &
            'Phase 0 requires PHASE0_YIELD_TABLE; embedded fallback is disabled'
       return
    else
       inquire(file=trim(filename), exist=exists)
       if (.not. exists) then
          ierr = 2
          initialization_ierr = ierr
          if (myid == 1) write(*,*) 'Phase 0 yield table not found: ', trim(filename)
          return
       end if
       call load_yield_backend(trim(filename), .false., yield_table, table_ierr)
    end if
    if (table_ierr /= backend_ok) then
       ierr = 10 + table_ierr
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 yield table load failed: ', table_ierr
       return
    end if
    ! Production source evaluation must not invent terminal outcomes between
    ! discrete source mass nodes.  The mode affects only mass assignment;
    ! metallicity, age, and any source-node fate axes remain independently
    ! governed by their resolver contracts.
    call set_yield_mass_assignment_mode(yield_table, &
         yield_mass_assignment_piecewise_constant, assignment_ierr)
    if (assignment_ierr /= yield_table_ok) then
       ierr = 11
       initialization_ierr = ierr
       if (myid == 1) write(*,*) &
            'Phase 0 source-cell mass assignment mode is invalid: ', assignment_ierr
       return
    end if
    call audit_yield_table(yield_table, 1.0d-8, audit_ierr, .true., &
         channel_owns_terminal_remnant)
    if (audit_ierr /= 0) then
       ierr = 20 + audit_ierr
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 yield table audit failed: ', audit_ierr
       return
    end if

    if (.not. metal) then
       ierr = 30
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 requires metal=.true.'
       return
    end if
    call clear_field_map(runtime_field_map)
    runtime_field_map%density_index = 1
    runtime_field_map%momentum_index = (/2, 3, 4/)
    runtime_field_map%energy_index = inener
    runtime_field_map%total_metal_index = imetal
    if (delayed_cooling) then
       runtime_field_map%delayed_cooling_index = idelay
    else
       runtime_field_map%delayed_cooling_index = 0
    end if
    runtime_field_map%element_index = (/ (ichem + element - 1, &
         element = 1, n_stellar_elements) /)
    call validate_field_map(runtime_field_map, nvar, ndim, audit_ierr)
    if (audit_ierr /= 0) then
       ierr = 31
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 RAMSES field map is invalid: ', audit_ierr
       return
    end if

    initialized = .true.
    initialization_ierr = 0
    if (myid == 1) then
       write(*,*) 'Phase 0 stellar enrichment enabled'
       write(*,*) '  table rows = ', yield_table%n_rows
       write(*,*) '  total-metal field = ', imetal
       write(*,*) '  first element field = ', ichem
       write(*,*) '  mass assignment = piecewise source-cell'
    end if
  end subroutine phase0_initialize

  subroutine phase0_feedback(ilevel, kgrid, subnump, ierr)
    integer, intent(in) :: ilevel, kgrid, subnump
    integer, intent(out) :: ierr
    integer :: jgrid, npart1, jpart, ipart, next_part, igrid
    integer :: local_ierr

    ierr = 0
    if (.not. initialized) call phase0_initialize(ierr)
    if (ierr /= 0 .or. subnump <= 0 .or. kgrid <= 0) return

    igrid = kgrid
    do jgrid = 1, subnump
       npart1 = numbp(igrid)
       if (npart1 > 0) then
          ipart = headp(igrid)
          do jpart = 1, npart1
             next_part = nextp(ipart)
             if (ptypep(ipart) == PTYPE_STAR) then
                call deposit_one_star(ilevel, igrid, ipart, local_ierr)
                if (local_ierr /= 0) then
                   ierr = local_ierr
                   return
                end if
             end if
             ipart = next_part
          end do
       end if
       igrid = next(igrid)
    end do
  end subroutine phase0_feedback

  subroutine deposit_one_star(ilevel, parent_grid, ipart, ierr)
    integer, intent(in) :: ilevel, parent_grid, ipart
    integer, intent(out) :: ierr
    type(stellar_population_t) :: population
    type(stellar_source_t) :: source
    type(stellar_population_ledger_t) :: population_ledger
    type(stellar_progress_t) :: progress
    real(stellar_dp) :: age_code, code_dt, previous_age_gyr, age_gyr, dt_gyr
    real(stellar_dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(stellar_dp) :: scale_mass, scale_momentum, scale_energy
    real(stellar_dp) :: returned_code, snii_returned_code, volume, energy_density
    real(stellar_dp) :: energy_code, source_momentum_code(3)
    real(stellar_dp) :: ejecta_code(n_stellar_elements), metal_ejecta_code
    real(stellar_dp) :: untracked_ejecta_msun, metal_ejecta_msun
    real(stellar_dp) :: ledger_remaining_code, ledger_scale
    real(stellar_dp) :: particle_mass_scale
    real(stellar_dp) :: bulk_energy, bulk_momentum(3)
    integer :: source_ierr, target_cell, idim, element
    integer :: units_ierr, progress_ierr
    logical :: located, should_deposit

    ierr = 0
    age_code = texp - tpp(ipart)
    if (.not. ieee_is_finite(age_code)) then
       ierr = 40
       return
    end if

    call progress_initialize(progress, indtab(ipart), progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 41
       return
    end if
    call progress_begin(progress, age_code, source_tolerance, should_deposit, &
         code_dt, progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 42
       return
    end if
    if (.not. should_deposit) return

    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
    scale_mass = scale_d * scale_l**3 / solar_mass_cgs
    scale_momentum = scale_d * scale_l**3 * scale_v
    scale_energy = scale_d * scale_l**3 * scale_v**2
    if (.not. ieee_is_finite(scale_mass) .or. &
         .not. ieee_is_finite(scale_momentum) .or. &
         .not. ieee_is_finite(scale_energy) .or. &
         scale_mass <= 0.0d0 .or. scale_momentum <= 0.0d0 .or. &
         scale_energy <= 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 43
       return
    end if

    call code_time_to_age_gyr(age_code, scale_t, aexp, age_gyr, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 44
       return
    end if
    call code_time_to_age_gyr(progress%committed_age_code, scale_t, aexp, &
         previous_age_gyr, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 45
       return
    end if
    call code_interval_to_age_gyr(code_dt, scale_t, aexp, dt_gyr, units_ierr)
    if (units_ierr /= units_ok .or. dt_gyr <= 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 46
       return
    end if

    population%formation_time = tpp(ipart)
    call mass_code_to_msun(mp0(ipart), scale_mass, population%initial_mass, &
         units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 47
       return
    end if
    call mass_code_to_msun(mp(ipart), scale_mass, population%current_mass, &
         units_ierr)
    if (units_ierr /= units_ok .or. population%initial_mass < 0.0d0 .or. &
         population%current_mass < 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 48
       return
    end if
    if (.not. ieee_is_finite(zp(ipart)) .or. zp(ipart) < 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 49
       return
    end if
    population%birth_metallicity = zp(ipart)
    population%birth_mass_fraction = 0.0d0
    population%imf_id = default_imf_id
    population%imf_mass_min = configured_imf_mass_min
    population%imf_mass_max = configured_imf_mass_max
    population%population_id = population_model_id
    population%binary_fraction = configured_binary_fraction
    population%yield_basis_id = yield_source_basis_id
    population%pisn_enabled = enable_pisn
    if (population%initial_mass <= 0.0d0) then
       call progress_commit(progress, progress_ierr)
       if (progress_ierr /= progress_ok) then
          ierr = 50
          return
       end if
       call progress_export(progress, indtab(ipart), progress_ierr)
       if (progress_ierr /= progress_ok) ierr = 50
       return
    end if

    call compute_stellar_source_increment(yield_table, population, &
         previous_age_gyr, age_gyr, configured_channel_mass_min, &
         configured_channel_mass_max, &
         phase0_mass_bins, source, source_ierr, population_ledger)
    if (source_ierr /= 0) then
       call progress_abort(progress, progress_ierr)
       ierr = 50 + source_ierr
       if (myid == 1) write(*,*) 'Phase 0 source evaluation failed: ', source_ierr
       return
    end if

    call mass_msun_to_code(source%returned_mass, scale_mass, returned_code, &
         units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 59
       return
    end if
    particle_mass_scale = max(abs(mp0(ipart)), tiny(1.0d0))
    if (returned_code < -source_tolerance * particle_mass_scale .or. &
         returned_code > mp(ipart) + source_tolerance * particle_mass_scale) then
       call progress_abort(progress, progress_ierr)
       ierr = 60
       if (myid == 1) write(*,*) 'Phase 0 returned mass violates star ledger'
       return
    end if
    returned_code = max(0.0d0, min(returned_code, mp(ipart)))
    call mass_msun_to_code(population_ledger%living_mass + &
         population_ledger%remnant_mass, scale_mass, ledger_remaining_code, &
         units_ierr)
    ledger_scale = max(particle_mass_scale, abs(mp(ipart)), &
         abs(ledger_remaining_code))
    if (units_ierr /= units_ok .or. &
         ledger_remaining_code < -source_tolerance * particle_mass_scale .or. &
         abs(mp(ipart)-returned_code-ledger_remaining_code) > &
         source_tolerance * ledger_scale) then
       call progress_abort(progress, progress_ierr)
       ierr = 72
       if (myid == 1) write(*,*) &
            'Phase 0 cumulative population ledger does not match particle mass'
       return
    end if
    call mass_msun_to_code(delayed_cooling_source_mass(source), scale_mass, &
         snii_returned_code, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 62
       return
    end if
    if (snii_returned_code < -source_tolerance * particle_mass_scale .or. &
         snii_returned_code > returned_code + &
         source_tolerance * particle_mass_scale) then
       call progress_abort(progress, progress_ierr)
       ierr = 63
       if (myid == 1) write(*,*) 'Phase 0 SNII return violates channel ledger'
       return
    end if
    snii_returned_code = max(0.0d0, min(snii_returned_code, returned_code))
    do element = 1, n_stellar_elements
       call mass_msun_to_code(source%ejected_mass(element), scale_mass, &
            ejecta_code(element), units_ierr)
       if (units_ierr /= units_ok) then
          call progress_abort(progress, progress_ierr)
          ierr = 64
          return
       end if
    end do
    if (any(ejecta_code < -source_tolerance * particle_mass_scale)) then
       call progress_abort(progress, progress_ierr)
       ierr = 65
       return
    end if
    ejecta_code = max(0.0d0, ejecta_code)
    untracked_ejecta_msun = untracked_ejecta_mass(source%returned_mass, &
         source%ejected_mass)
    if (untracked_ejecta_msun < -source_tolerance * &
         max(population%initial_mass, abs(source%returned_mass), &
         tiny(1.0d0))) then
       call progress_abort(progress, progress_ierr)
       ierr = 65
       if (myid == 1) write(*,*) &
            'Phase 0 tracked ejecta exceed returned mass'
       return
    end if
    untracked_ejecta_msun = max(0.0d0, untracked_ejecta_msun)
    metal_ejecta_msun = generic_metal_ejecta_mass(source%returned_mass, &
         source%ejected_mass)
    call mass_msun_to_code(metal_ejecta_msun, scale_mass, &
         metal_ejecta_code, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 65
       return
    end if

    call locate_star_cell(ilevel, parent_grid, ipart, target_cell, volume, located)
    if (.not. located) then
       call progress_abort(progress, progress_ierr)
       return
    end if

    if (.not. ieee_is_finite(volume) .or. volume <= 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 68
       return
    end if
    call momentum_cgs_to_code(source%momentum, scale_momentum, &
         source_momentum_code, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 66
       return
    end if
    call energy_erg_to_code(source%energy, scale_energy, energy_code, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 67
       return
    end if
    ! Rebuild the returned-particle momentum after the physical source unit
    ! conversion has succeeded.  All writes below are part of one transaction.
    ! source%energy is the non-bulk event energy, so include event momentum
    ! kinetic energy and its cross-term with the gas-frame bulk velocity.
    if (returned_code > 0.0d0) then
       bulk_energy = 0.5d0 * returned_code * sum(vp(ipart,1:3)**2) + &
            sum(vp(ipart,1:3) * source_momentum_code) + &
            0.5d0 * sum(source_momentum_code**2) / returned_code
    else if (maxval(abs(source_momentum_code)) > 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 68
       return
    else
       bulk_energy = 0.0d0
    end if
    bulk_momentum = returned_code * vp(ipart,1:3) + source_momentum_code
    energy_density = (energy_code + bulk_energy) / volume

    if (.not. ieee_is_finite(energy_density)) then
       call progress_abort(progress, progress_ierr)
       ierr = 68
       return
    end if
    unew(target_cell,runtime_field_map%density_index) = &
         unew(target_cell,runtime_field_map%density_index) + returned_code / volume
    do idim = 1, ndim
       unew(target_cell,runtime_field_map%momentum_index(idim)) = &
            unew(target_cell,runtime_field_map%momentum_index(idim)) + &
            bulk_momentum(idim) / volume
    end do
    unew(target_cell,runtime_field_map%energy_index) = &
         unew(target_cell,runtime_field_map%energy_index) + energy_density

    ! Delayed cooling represents unresolved core-collapse SN blast energy.
    ! Do not load this reservoir with winds, AGB, SNIa, or their combined
    ! mass return.  The legacy feedback_mode retains that historical rule.
    if (delayed_cooling) then
       unew(target_cell,runtime_field_map%delayed_cooling_index) = &
            unew(target_cell,runtime_field_map%delayed_cooling_index) + &
            snii_returned_code / volume
    end if

    ! The generic metal scalar must retain all metals even when an individual
    ! element field is disabled.  The residual contains source-declared
    ! elements outside the reduced eleven-species network.
    unew(target_cell,runtime_field_map%total_metal_index) = &
         unew(target_cell,runtime_field_map%total_metal_index) + &
         metal_ejecta_code / volume
    do element = 1, n_stellar_elements
       if (.not. active_element(element)) cycle
       unew(target_cell,runtime_field_map%element_index(element)) = &
            unew(target_cell,runtime_field_map%element_index(element)) + &
            ejecta_code(element) / volume
    end do

    mp(ipart) = mp(ipart) - returned_code
    call progress_commit(progress, progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 69
       return
    end if
    call progress_export(progress, indtab(ipart), progress_ierr)
    if (progress_ierr /= progress_ok) ierr = 70
  end subroutine deposit_one_star

  subroutine locate_star_cell(ilevel, parent_grid, ipart, target_cell, volume, located)
    integer, intent(in) :: ilevel, parent_grid, ipart
    integer, intent(out) :: target_cell
    real(stellar_dp), intent(out) :: volume
    logical, intent(out) :: located
    integer :: idim, nx_loc, grid_son, cell_number
    integer :: id(3), parent_id(3), cell_id(3)
    integer :: father_cell(1:nvector), kg(1:nvector)
    integer :: nbors_father_cells(1:nvector,1:threetondim)
    integer :: nbors_father_grids(1:nvector,1:twotondim)
    real(stellar_dp) :: dx, scale, dx_loc
    real(stellar_dp) :: x0(3), xpart(3), coordinate(3), skip_loc(3)

    located = .false.
    target_cell = 0
    volume = 0.0d0
    dx = 0.5d0**ilevel
    nx_loc = icoarse_max - icoarse_min + 1
    scale = boxlen / dble(nx_loc)
    dx_loc = dx * scale
    volume = dx_loc**ndim
    skip_loc = 0.0d0
    if (ndim > 0) skip_loc(1) = dble(icoarse_min)
    if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
    if (ndim > 2) skip_loc(3) = dble(kcoarse_min)

    x0 = 0.0d0
    xpart = 0.0d0
    coordinate = 0.0d0
    do idim = 1, ndim
       x0(idim) = xg(parent_grid,idim) - 3.0d0*dx
       xpart(idim) = xp(ipart,idim) / scale + skip_loc(idim)
       coordinate(idim) = (xpart(idim) - x0(idim)) / dx
       id(idim) = int(coordinate(idim))
       parent_id(idim) = id(idim) / 2
       cell_id(idim) = id(idim) - 2*parent_id(idim)
    end do
    if (ndim < 3) then
       id(3) = 0
       parent_id(3) = 0
       cell_id(3) = 0
    end if
    kg(1) = 1 + parent_id(1) + 3*parent_id(2) + 9*parent_id(3)
    father_cell(1) = father(parent_grid)
    call get3cubefather(father_cell, nbors_father_cells, nbors_father_grids, 1, ilevel)
    grid_son = son(nbors_father_cells(1,kg(1)))
    if (grid_son <= 0) return
    cell_number = 1 + cell_id(1) + 2*cell_id(2) + 4*cell_id(3)
    target_cell = ICELL_OF(grid_son,cell_number)
    located = target_cell > 0
  end subroutine locate_star_cell

end module stellar_ramses_runtime

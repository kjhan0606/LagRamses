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
       n_stellar_channels, active_element, enable_wind, enable_agb, &
       enable_snii, enable_snia, enable_pisn, default_imf_id, &
       population_model_id, yield_source_basis_id, configured_imf_mass_min, &
       configured_imf_mass_max, configured_binary_fraction, &
       configured_channel_mass_min, configured_channel_mass_max, &
       channel_owns_terminal_remnant, &
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
       code_interval_to_age_gyr, units_ok, solar_mass_cgs
  use stellar_progress_contract, only: stellar_progress_t, &
       progress_initialize, progress_begin, progress_commit, &
       progress_abort, progress_export, progress_ok
#include "amr_index.h"
  implicit none
  private

  real(stellar_dp), parameter :: source_tolerance = 1.0d-10
  integer, parameter :: phase0_mass_bins = 64

  type(stellar_yield_table_t), save :: yield_table
  logical, save :: initialized = .false.
  integer, save :: initialization_ierr = 0
  character(len=1024), save :: loaded_yield_table_path = ''
  integer, save :: loaded_yield_table_rows = 0

  public :: phase0_initialize
  public :: phase0_feedback
  public :: phase0_get_runtime_identity

contains

  subroutine phase0_initialize(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: status, table_ierr, audit_ierr, coverage_ierr, assignment_ierr
    logical :: exists

    if (initialized .or. initialization_ierr /= 0) then
       ierr = initialization_ierr
       return
    end if

    ierr = 0
    loaded_yield_table_path = ''
    loaded_yield_table_rows = 0
    if (.not. production_source_model_supported()) then
       ierr = 3
       initialization_ierr = ierr
       if (myid == 1) write(*,*) &
            'Phase 0 source model is not implemented for production'
       return
    end if
    call get_environment_variable('PHASE0_YIELD_TABLE', filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) then
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
    call audit_enabled_channel_coverage(yield_table, coverage_ierr)
    if (coverage_ierr /= 0) then
       ierr = 120 + coverage_ierr
       initialization_ierr = ierr
       if (myid == 1) write(*,*) &
            'Phase 0 yield table does not cover enabled channel: ', coverage_ierr
       return
    end if

    if (.not. metal) then
       ierr = 30
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 requires metal=.true.'
       return
    end if
    if (imetal < 1 .or. ichem < 1 .or. ichem+n_stellar_elements-1 > nvar) then
       ierr = 31
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 chemical field map exceeds NVAR'
       return
    end if

    initialized = .true.
    initialization_ierr = 0
    loaded_yield_table_path = trim(filename)
    loaded_yield_table_rows = yield_table%n_rows
    if (myid == 1) then
       write(*,*) 'Phase 0 stellar enrichment enabled'
       write(*,*) '  table rows = ', yield_table%n_rows
       write(*,*) '  total-metal field = ', imetal
       write(*,*) '  first element field = ', ichem
       write(*,*) '  mass assignment = piecewise source-cell'
    end if
  end subroutine phase0_initialize

  subroutine phase0_get_runtime_identity(path, n_rows, is_loaded)
    character(len=*), intent(out) :: path
    integer, intent(out) :: n_rows
    logical, intent(out) :: is_loaded

    path = loaded_yield_table_path
    n_rows = loaded_yield_table_rows
    is_loaded = initialized
  end subroutine phase0_get_runtime_identity

  subroutine audit_enabled_channel_coverage(table, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(out) :: ierr
    logical :: enabled(n_stellar_channels), found
    real(stellar_dp) :: table_min, table_max, scale
    integer :: channel, row

    enabled = (/enable_wind, enable_agb, enable_snii, enable_snia, enable_pisn/)
    ierr = 0
    do channel = 1, n_stellar_channels
       if (.not. enabled(channel)) cycle
       found = .false.
       table_min = huge(1.0_stellar_dp)
       table_max = -huge(1.0_stellar_dp)
       do row = 1, table%n_rows
          if (table%channel(row) /= channel) cycle
          found = .true.
          table_min = min(table_min, table%initial_mass(row))
          table_max = max(table_max, table%initial_mass(row))
       end do
       scale = max(1.0_stellar_dp, abs(configured_channel_mass_min(channel)), &
            abs(configured_channel_mass_max(channel)))
       if (.not. found .or. &
            table_min > configured_channel_mass_min(channel) + &
                 source_tolerance * scale .or. &
            table_max < configured_channel_mass_max(channel) - &
                 source_tolerance * scale) then
          ierr = channel
          return
       end if
    end do
  end subroutine audit_enabled_channel_coverage

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
    real(stellar_dp) :: age_code, code_dt, previous_age_gyr, age_gyr, dt_gyr
    real(stellar_dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(stellar_dp) :: scale_mass, scale_momentum, scale_energy
    real(stellar_dp) :: returned_code, snii_returned_code, volume, energy_density
    real(stellar_dp) :: source_momentum_code(3)
    real(stellar_dp) :: ejecta_code(n_stellar_elements), metal_ejecta_code
    real(stellar_dp) :: untracked_ejecta_msun, metal_ejecta_msun
    real(stellar_dp) :: ledger_remaining_code, ledger_scale
    real(stellar_dp) :: particle_mass_scale
    real(stellar_dp) :: bulk_energy, bulk_momentum(3)
    integer :: source_ierr, target_cell, idim, element
    integer :: progress_ierr, units_ierr
    logical :: located, should_deposit
    type(stellar_progress_t) :: progress

    ierr = 0
    age_code = texp - tpp(ipart)
    call progress_initialize(progress, indtab(ipart), progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 40
       return
    end if
    call progress_begin(progress, age_code, source_tolerance, should_deposit, &
         code_dt, progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 41
       return
    end if
    if (.not. should_deposit) return

    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
    scale_mass = scale_d * scale_l**3 / solar_mass_cgs
    scale_momentum = scale_d * scale_l**3 * scale_v
    scale_energy = scale_d * scale_l**3 * scale_v**2
    if (scale_mass <= 0.0d0 .or. scale_momentum <= 0.0d0 .or. &
         scale_energy <= 0.0d0) then
       call progress_abort(progress, progress_ierr)
       ierr = 42
       return
    end if

    call code_time_to_age_gyr(age_code, scale_t, aexp, age_gyr, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 43
       return
    end if
    call code_time_to_age_gyr(progress%committed_age_code, scale_t, aexp, &
         previous_age_gyr, units_ierr)
    if (units_ierr /= units_ok) then
       call progress_abort(progress, progress_ierr)
       ierr = 44
       return
    end if
    call code_interval_to_age_gyr(code_dt, scale_t, aexp, dt_gyr, units_ierr)
    if (units_ierr /= units_ok .or. dt_gyr <= 0.0d0 .or. &
         age_gyr <= previous_age_gyr) then
       call progress_abort(progress, progress_ierr)
       ierr = 45
       return
    end if

    population%formation_time = tpp(ipart)
    population%initial_mass = mp0(ipart) * scale_mass
    population%current_mass = mp(ipart) * scale_mass
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
          ierr = 46
          return
       end if
       call progress_export(progress, indtab(ipart), progress_ierr)
       if (progress_ierr /= progress_ok) ierr = 47
       return
    end if

    call compute_stellar_source_increment(yield_table, population, &
         previous_age_gyr, age_gyr, configured_channel_mass_min, &
         configured_channel_mass_max, &
         phase0_mass_bins, source, source_ierr, population_ledger)
    if (source_ierr /= 0) then
       ierr = 50 + source_ierr
       if (myid == 1) write(*,*) 'Phase 0 source evaluation failed: ', source_ierr
       call progress_abort(progress, progress_ierr)
       return
    end if

    returned_code = source%returned_mass / scale_mass
    particle_mass_scale = max(abs(mp0(ipart)), tiny(1.0d0))
    if (returned_code < -source_tolerance * particle_mass_scale .or. &
         returned_code > mp(ipart) + source_tolerance * particle_mass_scale) then
       ierr = 60
       if (myid == 1) write(*,*) 'Phase 0 returned mass violates star ledger'
       call progress_abort(progress, progress_ierr)
       return
    end if
    returned_code = max(0.0d0, min(returned_code, mp(ipart)))
    ledger_remaining_code = (population_ledger%living_mass + &
         population_ledger%remnant_mass) / scale_mass
    ledger_scale = max(particle_mass_scale, abs(mp(ipart)), &
         abs(ledger_remaining_code))
    if (.not. ieee_is_finite(ledger_remaining_code) .or. &
         ledger_remaining_code < -source_tolerance * particle_mass_scale .or. &
         abs(mp(ipart)-returned_code-ledger_remaining_code) > &
         source_tolerance * ledger_scale) then
       ierr = 72
       if (myid == 1) write(*,*) &
            'Phase 0 cumulative population ledger does not match particle mass'
       call progress_abort(progress, progress_ierr)
       return
    end if
    snii_returned_code = delayed_cooling_source_mass(source) / scale_mass
    if (snii_returned_code < -source_tolerance * particle_mass_scale .or. &
         snii_returned_code > returned_code + &
         source_tolerance * particle_mass_scale) then
       ierr = 62
       if (myid == 1) write(*,*) 'Phase 0 SNII return violates channel ledger'
       call progress_abort(progress, progress_ierr)
       return
    end if
    snii_returned_code = max(0.0d0, min(snii_returned_code, returned_code))
    ejecta_code = source%ejected_mass / scale_mass
    if (any(ejecta_code < -source_tolerance * particle_mass_scale)) then
       ierr = 61
       call progress_abort(progress, progress_ierr)
       return
    end if
    ejecta_code = max(0.0d0, ejecta_code)
    untracked_ejecta_msun = untracked_ejecta_mass(source%returned_mass, &
         source%ejected_mass)
    if (untracked_ejecta_msun < -source_tolerance * &
         max(population%initial_mass, abs(source%returned_mass), &
         tiny(1.0d0))) then
       ierr = 65
       if (myid == 1) write(*,*) 'Phase 0 tracked ejecta exceed returned mass'
       call progress_abort(progress, progress_ierr)
       return
    end if
    metal_ejecta_msun = generic_metal_ejecta_mass(source%returned_mass, &
         source%ejected_mass)
    metal_ejecta_code = metal_ejecta_msun / scale_mass
    if (.not. ieee_is_finite(metal_ejecta_code) .or. &
         metal_ejecta_code < -source_tolerance * particle_mass_scale) then
       ierr = 66
       call progress_abort(progress, progress_ierr)
       return
    end if
    metal_ejecta_code = max(0.0d0, metal_ejecta_code)

    call locate_star_cell(ilevel, parent_grid, ipart, target_cell, volume, located)
    if (.not. located) then
       call progress_abort(progress, progress_ierr)
       return
    end if

    source_momentum_code = source%momentum / scale_momentum
    bulk_momentum = returned_code * vp(ipart,1:3)
    ! source%energy is the non-bulk event energy.  Once a source carries
    ! event momentum, conserved total energy also contains its kinetic term
    ! and the bulk/event cross-term.  The current SSP path has zero source
    ! momentum; the expression keeps the future SNIa convention explicit.
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
    do idim = 1, 3
       bulk_momentum(idim) = bulk_momentum(idim) + &
            source%momentum(idim) / scale_momentum
    end do
    energy_density = (source%energy / scale_energy + bulk_energy) / volume

    unew(target_cell,1) = unew(target_cell,1) + returned_code / volume
    do idim = 1, ndim
       unew(target_cell,1+idim) = unew(target_cell,1+idim) + &
            bulk_momentum(idim) / volume
    end do
    unew(target_cell,5) = unew(target_cell,5) + energy_density

    ! Delayed cooling represents unresolved core-collapse SN blast energy.
    ! Do not load this reservoir with winds, AGB, SNIa, or their combined
    ! mass return.  The legacy feedback_mode retains that historical rule.
    if (delayed_cooling) then
       unew(target_cell,idelay) = unew(target_cell,idelay) + &
            snii_returned_code / volume
    end if

    ! Generic metallicity retains tracked metals even when their individual
    ! fields are disabled, plus the source-declared untracked residual.
    unew(target_cell,imetal) = unew(target_cell,imetal) + &
         metal_ejecta_code / volume
    do element = 1, n_stellar_elements
       if (.not. active_element(element)) cycle
       unew(target_cell,ichem+element-1) = &
            unew(target_cell,ichem+element-1) + ejecta_code(element) / volume
    end do

    mp(ipart) = mp(ipart) - returned_code
    call progress_commit(progress, progress_ierr)
    if (progress_ierr /= progress_ok) then
       ierr = 70
       return
    end if
    call progress_export(progress, indtab(ipart), progress_ierr)
    if (progress_ierr /= progress_ok) ierr = 71
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

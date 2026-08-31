! Phase 0 runtime adapter for the RAMSES stellar-particle feedback loop.
!
! The legacy feedback routine interpolates a prompt table and overlaps the
! total-metal field with the first chemical field.  This adapter evaluates an
! age increment through the Phase 0 SSP driver and deposits total metal plus
! eleven independently tracked ejecta fields without that overlap.

module stellar_ramses_runtime
  use amr_commons
  use pm_commons
  use hydro_commons
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, elem_c, active_element, enable_pisn
  use stellar_enrichment_contract, only: stellar_population_t, stellar_source_t
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_yield_backend, only: load_yield_backend, backend_ok
  use stellar_yield_audit, only: audit_yield_table
  use stellar_enrichment_driver, only: compute_stellar_source_increment
#include "amr_index.h"
  implicit none
  private

  real(stellar_dp), parameter :: seconds_per_gyr = 1.0d9 * 365.25d0 * 86400.0d0
  real(stellar_dp), parameter :: solar_mass_cgs = 1.98847d33
  real(stellar_dp), parameter :: source_tolerance = 1.0d-10
  integer, parameter :: phase0_mass_bins = 64

  type(stellar_yield_table_t), save :: yield_table
  logical, save :: initialized = .false.
  integer, save :: initialization_ierr = 0

  public :: phase0_initialize
  public :: phase0_feedback

contains

  subroutine phase0_initialize(ierr)
    integer, intent(out) :: ierr
    character(len=1024) :: filename
    integer :: status, table_ierr, audit_ierr
    logical :: exists

    if (initialized .or. initialization_ierr /= 0) then
       ierr = initialization_ierr
       return
    end if

    ierr = 0
    call get_environment_variable('PHASE0_YIELD_TABLE', filename, status=status)
    if (status /= 0 .or. len_trim(filename) == 0) then
       ! The embedded backend is an integration-test fallback. Production
       ! runs should set PHASE0_YIELD_TABLE to an age-resolved external table.
       call load_yield_backend('', .true., yield_table, table_ierr)
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
    call audit_yield_table(yield_table, 1.0d-8, audit_ierr)
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
    if (imetal < 1 .or. ichem < 1 .or. ichem+n_stellar_elements-1 > nvar) then
       ierr = 31
       initialization_ierr = ierr
       if (myid == 1) write(*,*) 'Phase 0 chemical field map exceeds NVAR'
       return
    end if

    initialized = .true.
    initialization_ierr = 0
    if (myid == 1) then
       write(*,*) 'Phase 0 stellar enrichment enabled'
       write(*,*) '  table rows = ', yield_table%n_rows
       write(*,*) '  total-metal field = ', imetal
       write(*,*) '  first element field = ', ichem
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
    real(stellar_dp) :: age_code, previous_age_code, age_gyr, dt_gyr
    real(stellar_dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
    real(stellar_dp) :: scale_mass, scale_momentum, scale_energy
    real(stellar_dp) :: returned_code, volume, energy_density
    real(stellar_dp) :: ejecta_code(n_stellar_elements), metal_ejecta_code
    real(stellar_dp) :: bulk_energy, bulk_momentum(3)
    real(stellar_dp) :: channel_mass_min(n_stellar_channels)
    real(stellar_dp) :: channel_mass_max(n_stellar_channels)
    integer :: source_ierr, target_cell, idim, element
    logical :: located

    ierr = 0
    age_code = texp - tpp(ipart)
    previous_age_code = max(0.0d0, indtab(ipart))
    if (age_code <= previous_age_code + source_tolerance) return
    if (aexp <= 0.0d0) then
       ierr = 40
       return
    end if

    call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
    scale_mass = scale_d * scale_l**3 / solar_mass_cgs
    scale_momentum = scale_d * scale_l**3 * scale_v
    scale_energy = scale_d * scale_l**3 * scale_v**2
    if (scale_mass <= 0.0d0 .or. scale_momentum <= 0.0d0 .or. &
         scale_energy <= 0.0d0) then
       ierr = 41
       return
    end if

    age_gyr = age_code * scale_t / (seconds_per_gyr * aexp**2)
    dt_gyr = (age_code - previous_age_code) * scale_t / &
         (seconds_per_gyr * aexp**2)
    if (age_gyr < 0.0d0 .or. dt_gyr <= 0.0d0) return

    population%formation_time = tpp(ipart)
    population%initial_mass = mp0(ipart) * scale_mass
    population%current_mass = mp(ipart) * scale_mass
    population%birth_metallicity = max(0.0d0, zp(ipart))
    population%birth_mass_fraction = 0.0d0
    population%imf_id = 1
    population%population_id = 0
    population%pisn_enabled = enable_pisn
    if (population%initial_mass <= 0.0d0) then
       indtab(ipart) = age_code
       return
    end if

    channel_mass_min = (/0.8d0, 1.0d0, 8.0d0, 3.0d0, 140.0d0/)
    channel_mass_max = (/120.0d0, 8.0d0, 40.0d0, 8.0d0, 260.0d0/)
    call compute_stellar_source_increment(yield_table, population, age_gyr, &
         dt_gyr, channel_mass_min, channel_mass_max, phase0_mass_bins, source, &
         source_ierr)
    if (source_ierr /= 0) then
       ierr = 50 + source_ierr
       if (myid == 1) write(*,*) 'Phase 0 source evaluation failed: ', source_ierr
       return
    end if

    returned_code = source%returned_mass / scale_mass
    if (returned_code < -source_tolerance .or. returned_code > &
         mp(ipart) * (1.0d0 + source_tolerance)) then
       ierr = 60
       if (myid == 1) write(*,*) 'Phase 0 returned mass violates star ledger'
       return
    end if
    returned_code = max(0.0d0, min(returned_code, mp(ipart)))
    ejecta_code = source%ejected_mass / scale_mass
    if (any(ejecta_code < -source_tolerance)) then
       ierr = 61
       return
    end if
    ejecta_code = max(0.0d0, ejecta_code)

    call locate_star_cell(ilevel, parent_grid, ipart, target_cell, volume, located)
    if (.not. located) return

    bulk_momentum = returned_code * vp(ipart,1:3)
    bulk_energy = 0.5d0 * returned_code * sum(vp(ipart,1:3)**2)
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

    metal_ejecta_code = 0.0d0
    do element = elem_c, n_stellar_elements
       if (active_element(element)) metal_ejecta_code = &
            metal_ejecta_code + ejecta_code(element)
    end do
    unew(target_cell,imetal) = unew(target_cell,imetal) + &
         metal_ejecta_code / volume
    do element = 1, n_stellar_elements
       if (.not. active_element(element)) cycle
       unew(target_cell,ichem+element-1) = &
            unew(target_cell,ichem+element-1) + ejecta_code(element) / volume
    end do

    mp(ipart) = mp(ipart) - returned_code
    indtab(ipart) = age_code
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

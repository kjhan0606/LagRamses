! Phase 0 RAMSES conserved-variable bridge.
!
! RAMSES-specific variable indices are supplied by the caller.  This keeps
! the stellar source modules independent of a particular lagRamses nvar
! layout while allowing the deposition operation to update uold directly.
! The SNIa bridge below additionally requires the code length, density, and
! velocity scales so physical event increments cannot be written into uold as
! if they were already code units.

module stellar_ramses_bridge
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_enrichment_contract, only: stellar_source_t, &
       delayed_cooling_source_mass, generic_metal_ejecta_mass
  use stellar_snia_physical_contract, only: snia_event_budget_t
  use stellar_snia_cell_deposition, only: snia_thermal_coupling_t, &
       snia_cell_increment_t, &
       snia_deposition_ok, build_snia_cell_increment
  use stellar_ramses_field_map, only: stellar_field_map_t, validate_field_map
  implicit none

  private
  integer, parameter, public :: ramses_bridge_ok = 0
  integer, parameter, public :: ramses_bridge_err_argument = 1
  integer, parameter, public :: ramses_bridge_err_source = 2
  integer, parameter, public :: ramses_bridge_err_index = 4
  integer, parameter, public :: ramses_bridge_err_closure = 8
  integer, parameter, public :: ramses_bridge_err_snia = 16
  integer, parameter, public :: ramses_bridge_err_result = 32
  integer, parameter, public :: ramses_bridge_err_target = 64

  public :: deposit_source_to_uold
  public :: build_stellar_source_unew_delta
  public :: deposit_snia_budget_to_uold
  public :: build_snia_budget_unew_delta
  public :: deposit_snia_budget_to_unew

contains

  subroutine deposit_source_to_uold(source, nvar, n_cells, cell_volume, &
       weights, density_var, energy_var, momentum_var, element_var, uold, &
       tolerance, ierr, total_metal_var)
    type(stellar_source_t), intent(in) :: source
    integer, intent(in) :: nvar, n_cells
    real(stellar_dp), intent(in) :: cell_volume(n_cells)
    real(stellar_dp), intent(in) :: weights(n_cells)
    integer, intent(in) :: density_var, energy_var
    integer, intent(in) :: momentum_var(3)
    integer, intent(in) :: element_var(n_stellar_elements)
    real(stellar_dp), intent(inout) :: uold(nvar,n_cells)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    integer, intent(in), optional :: total_metal_var

    real(stellar_dp) :: tol, weight_sum, ejected_sum, scale
    real(stellar_dp) :: normalized_weight, cell_mass, element_mass
    real(stellar_dp) :: generic_metal_mass
    integer :: cell, element

    ierr = ramses_bridge_ok
    tol = max(tolerance, 1.0e-12_stellar_dp)

    if (nvar <= 0 .or. n_cells <= 0 .or. density_var < 1 .or. &
         density_var > nvar .or. energy_var < 1 .or. energy_var > nvar .or. &
         any(momentum_var < 1) .or. any(momentum_var > nvar)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (any(cell_volume <= 0.0_stellar_dp) .or. &
         any(weights < 0.0_stellar_dp)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (any(element_var < 0) .or. any(element_var > nvar)) then
       ierr = ramses_bridge_err_index
       return
    end if
    if (present(total_metal_var)) then
       if (total_metal_var < 0 .or. total_metal_var > nvar) then
          ierr = ramses_bridge_err_index
          return
       end if
    end if

    weight_sum = sum(weights)
    if (weight_sum <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    ejected_sum = sum(source%ejected_mass)
    scale = max(tiny(1.0_stellar_dp), abs(source%returned_mass), &
         abs(ejected_sum))
    if (source%returned_mass < -tol * scale .or. &
         minval(source%ejected_mass) < -tol * scale) then
       ierr = ramses_bridge_err_source
       return
    end if
    if (ejected_sum > source%returned_mass + tol * scale) then
       ierr = ramses_bridge_err_closure
       return
    end if
    generic_metal_mass = generic_metal_ejecta_mass(source%returned_mass, &
         source%ejected_mass)

    ! All validation is completed before uold is modified.
    do cell = 1, n_cells
       normalized_weight = weights(cell) / weight_sum
       cell_mass = normalized_weight * source%returned_mass / cell_volume(cell)
       uold(density_var,cell) = uold(density_var,cell) + cell_mass
       if (present(total_metal_var)) then
          if (total_metal_var > 0) then
             uold(total_metal_var,cell) = uold(total_metal_var,cell) + &
                  normalized_weight * generic_metal_mass / cell_volume(cell)
          end if
       end if

       do element = 1, n_stellar_elements
          if (element_var(element) == 0) cycle
          element_mass = normalized_weight * source%ejected_mass(element) / &
               cell_volume(cell)
          uold(element_var(element),cell) = &
               uold(element_var(element),cell) + element_mass
       end do

       uold(energy_var,cell) = uold(energy_var,cell) + &
            normalized_weight * source%energy / cell_volume(cell)
       uold(momentum_var(1),cell) = uold(momentum_var(1),cell) + &
            normalized_weight * source%momentum(1) / cell_volume(cell)
       uold(momentum_var(2),cell) = uold(momentum_var(2),cell) + &
            normalized_weight * source%momentum(2) / cell_volume(cell)
       uold(momentum_var(3),cell) = uold(momentum_var(3),cell) + &
            normalized_weight * source%momentum(3) / cell_volume(cell)
    end do
  end subroutine deposit_source_to_uold

  subroutine build_stellar_source_unew_delta(source, bulk_velocity_code, &
       scale_mass, scale_momentum, scale_energy, volume_code, nvar, ndim, &
       field_map, delta, tolerance, ierr, channel_resolved)
    ! Convert one generic stellar source into a complete row-major RAMSES
    ! delta.  This routine is deliberately non-mutating so generic and SNIa
    ! deltas can be combined before the runtime lock/commit boundary.
    type(stellar_source_t), intent(in) :: source
    real(stellar_dp), intent(in) :: bulk_velocity_code(3)
    real(stellar_dp), intent(in) :: scale_mass, scale_momentum, scale_energy
    real(stellar_dp), intent(in) :: volume_code
    integer, intent(in) :: nvar, ndim
    type(stellar_field_map_t), intent(in) :: field_map
    real(stellar_dp), intent(out) :: delta(:)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    logical, intent(in), optional :: channel_resolved

    real(stellar_dp) :: tol, returned_code, snii_code, source_momentum_code(3)
    real(stellar_dp) :: generic_metal_code, bulk_energy, ejected_sum, scale
    real(stellar_dp) :: channel_mass, channel_mass_sum, channel_p(3), momentum_sum(3)
    logical :: resolve_channels
    integer :: map_ierr, element, idim, channel

    ierr = ramses_bridge_ok
    delta = 0.0_stellar_dp
    tol = max(tolerance, 1.0e-12_stellar_dp)
    if (nvar <= 0 .or. size(delta) < nvar .or. ndim /= 3 .or. &
         .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_stellar_dp .or. &
         .not. ieee_is_finite(volume_code) .or. volume_code <= 0.0_stellar_dp .or. &
         .not. all(ieee_is_finite(bulk_velocity_code)) .or. &
         .not. ieee_is_finite(scale_mass) .or. scale_mass <= 0.0_stellar_dp .or. &
         .not. ieee_is_finite(scale_momentum) .or. scale_momentum <= 0.0_stellar_dp .or. &
         .not. ieee_is_finite(scale_energy) .or. scale_energy <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    call validate_field_map(field_map, nvar, ndim, map_ierr)
    if (map_ierr /= 0) then
       ierr = ramses_bridge_err_index
       return
    end if
    if (.not. ieee_is_finite(source%returned_mass) .or. &
         .not. ieee_is_finite(source%energy) .or. &
         .not. all(ieee_is_finite(source%momentum)) .or. &
         .not. all(ieee_is_finite(source%ejected_mass)) .or. &
         .not. all(ieee_is_finite(source%net_yield)) .or. &
         .not. all(ieee_is_finite(source%channel_returned_mass)) .or. &
         .not. all(ieee_is_finite(source%channel_energy)) .or. &
         .not. all(ieee_is_finite(source%channel_momentum)) .or. &
         .not. all(ieee_is_finite(source%channel_ejected_mass)) .or. &
         .not. all(ieee_is_finite(source%channel_net_yield)) .or. &
         source%returned_mass < 0.0_stellar_dp .or. &
         source%energy < 0.0_stellar_dp .or. &
         minval(source%ejected_mass) < 0.0_stellar_dp) then
       ierr = ramses_bridge_err_source
       return
    end if
    ejected_sum = sum(source%ejected_mass)
    scale = max(1.0_stellar_dp, source%returned_mass, abs(ejected_sum))
    if (.not. ieee_is_finite(ejected_sum) .or. &
         ejected_sum > source%returned_mass + tol * scale) then
       ierr = ramses_bridge_err_closure
       return
    end if

    returned_code = source%returned_mass / scale_mass
    source_momentum_code = source%momentum / scale_momentum
    resolve_channels = .false.
    if (present(channel_resolved)) resolve_channels = channel_resolved
    if (resolve_channels) then
       ! Production SSP sources carry a complete channel mass/momentum ledger.
       ! Retain each component's kinetic energy before mixing in the cell:
       ! sum_c |m_c*v_star + p_c|^2/(2*m_c), not |sum_c p_c|^2/(2*sum_c m_c).
       ! A nonzero residual becomes gas internal energy on conservative mixing.
       channel_mass_sum = sum(source%channel_returned_mass)
       momentum_sum = sum(source%channel_momentum, dim=1)
       if (any(source%channel_returned_mass < 0.0_stellar_dp)) then
          ierr = ramses_bridge_err_source
          return
       end if
       scale = max(tiny(1.0_stellar_dp), source%returned_mass, channel_mass_sum)
       if (.not. ieee_is_finite(channel_mass_sum) .or. &
            .not. all(ieee_is_finite(momentum_sum)) .or. &
            abs(channel_mass_sum - source%returned_mass) > tol * scale) then
          ierr = ramses_bridge_err_closure
          return
       end if
       do idim = 1, 3
          scale = max(tiny(1.0_stellar_dp), abs(source%momentum(idim)), &
               sum(abs(source%channel_momentum(:,idim))))
          if (.not. ieee_is_finite(scale) .or. &
               abs(momentum_sum(idim) - source%momentum(idim)) > tol * scale) then
             ierr = ramses_bridge_err_closure
             return
          end if
       end do
       bulk_energy = 0.0_stellar_dp
       do channel = 1, size(source%channel_returned_mass)
          channel_mass = source%channel_returned_mass(channel) / scale_mass
          channel_p = source%channel_momentum(channel,:) / scale_momentum
          if (channel_mass > 0.0_stellar_dp) then
             bulk_energy = bulk_energy + 0.5_stellar_dp * &
                  sum((channel_mass * bulk_velocity_code + channel_p)**2) / channel_mass
          else if (any(channel_p /= 0.0_stellar_dp)) then
             ierr = ramses_bridge_err_source
             return
          end if
       end do
    else if (returned_code > 0.0_stellar_dp) then
       ! Aggregate-only callers retain their existing single-component contract.
       bulk_energy = 0.5_stellar_dp * returned_code * &
            sum(bulk_velocity_code**2) + &
            sum(bulk_velocity_code * source_momentum_code) + &
            0.5_stellar_dp * sum(source_momentum_code**2) / returned_code
    else if (maxval(abs(source_momentum_code)) > 0.0_stellar_dp) then
       ierr = ramses_bridge_err_source
       return
    else
       bulk_energy = 0.0_stellar_dp
    end if
    generic_metal_code = generic_metal_ejecta_mass(source%returned_mass, &
         source%ejected_mass) / scale_mass
    snii_code = delayed_cooling_source_mass(source) / scale_mass
    if (.not. ieee_is_finite(returned_code) .or. &
         .not. ieee_is_finite(bulk_energy) .or. &
         .not. ieee_is_finite(generic_metal_code) .or. &
         .not. ieee_is_finite(snii_code) .or. &
         snii_code < -tol * max(1.0_stellar_dp, abs(returned_code))) then
       ierr = ramses_bridge_err_result
       return
    end if
    if (snii_code < 0.0_stellar_dp .or. &
         snii_code > returned_code + tol * max(1.0_stellar_dp, abs(returned_code))) then
       ierr = ramses_bridge_err_closure
       return
    end if
    snii_code = max(0.0_stellar_dp, snii_code)
    generic_metal_code = max(0.0_stellar_dp, generic_metal_code)

    delta(field_map%density_index) = returned_code / volume_code
    do idim = 1, ndim
       delta(field_map%momentum_index(idim)) = &
            (returned_code * bulk_velocity_code(idim) + &
            source_momentum_code(idim)) / volume_code
    end do
    delta(field_map%energy_index) = &
         (source%energy / scale_energy + bulk_energy) / volume_code
    if (field_map%delayed_cooling_index /= 0) then
       delta(field_map%delayed_cooling_index) = snii_code / volume_code
    end if
    if (field_map%total_metal_index /= 0) then
       delta(field_map%total_metal_index) = generic_metal_code / volume_code
    end if
    do element = 1, n_stellar_elements
       if (field_map%element_index(element) == 0) cycle
       delta(field_map%element_index(element)) = &
            source%ejected_mass(element) / scale_mass / volume_code
    end do
    if (.not. all(ieee_is_finite(delta(1:nvar)))) then
       delta = 0.0_stellar_dp
       ierr = ramses_bridge_err_result
    end if
  end subroutine build_stellar_source_unew_delta

  subroutine deposit_snia_budget_to_uold(budget, coupling, bulk_velocity_cm_s, &
       scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, &
       n_cells, cell_volume_code, weights, density_var, energy_var, &
       momentum_var, uold, tolerance, ierr, total_metal_var, element_var)
    type(snia_event_budget_t), intent(in) :: budget
    type(snia_thermal_coupling_t), intent(in) :: coupling
    ! One gas-frame velocity per selected cell is required.  A single velocity
    ! for a distributed stencil would inject a spurious common bulk momentum.
    real(stellar_dp), intent(in) :: bulk_velocity_cm_s(3,n_cells)
    real(stellar_dp), intent(in) :: scale_length_cgs, scale_density_cgs
    real(stellar_dp), intent(in) :: scale_velocity_cgs
    integer, intent(in) :: nvar, n_cells
    real(stellar_dp), intent(in) :: cell_volume_code(n_cells)
    real(stellar_dp), intent(in) :: weights(n_cells)
    integer, intent(in) :: density_var, energy_var
    integer, intent(in) :: momentum_var(3)
    real(stellar_dp), intent(inout) :: uold(nvar,n_cells)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    integer, intent(in), optional :: total_metal_var
    integer, intent(in), optional :: element_var(n_stellar_elements)

    type(snia_cell_increment_t), allocatable :: increments(:)
    real(stellar_dp) :: weight_sum, normalized_weight
    real(stellar_dp) :: volume_cm3, density_scale, momentum_scale
    real(stellar_dp) :: energy_scale
    real(stellar_dp) :: ejected_sum, generic_metal_mass, tol
    integer :: cell, element, allocation_status, increment_ierr

    ierr = ramses_bridge_ok
    tol = max(tolerance, 1.0e-12_stellar_dp)
    if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_stellar_dp .or. &
         .not. all(ieee_is_finite(bulk_velocity_cm_s)) .or. &
         .not. ieee_is_finite(scale_length_cgs) .or. &
         .not. ieee_is_finite(scale_density_cgs) .or. &
         .not. ieee_is_finite(scale_velocity_cgs) .or. &
         scale_length_cgs <= 0.0_stellar_dp .or. &
         scale_density_cgs <= 0.0_stellar_dp .or. &
         scale_velocity_cgs <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (nvar <= 0 .or. n_cells <= 0 .or. density_var < 1 .or. &
         density_var > nvar .or. energy_var < 1 .or. energy_var > nvar .or. &
         any(momentum_var < 1) .or. any(momentum_var > nvar)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    if (energy_var == density_var .or. &
         any(momentum_var == density_var) .or. &
         any(momentum_var == energy_var)) then
       ierr = ramses_bridge_err_index
       return
    end if
    if (present(total_metal_var)) then
       if (total_metal_var < 1 .or. total_metal_var > nvar .or. &
            total_metal_var == density_var .or. total_metal_var == energy_var .or. &
            any(total_metal_var == momentum_var)) then
          ierr = ramses_bridge_err_index
          return
       end if
    end if
    if (present(element_var)) then
       do element = 1, n_stellar_elements
          if (element_var(element) < 0 .or. element_var(element) > nvar) then
             ierr = ramses_bridge_err_index
             return
          end if
          if (element_var(element) == 0) cycle
          if (element_var(element) == density_var .or. &
               element_var(element) == energy_var .or. &
               any(element_var(element) == momentum_var) .or. &
               (present(total_metal_var) .and. &
               element_var(element) == total_metal_var)) then
             ierr = ramses_bridge_err_index
             return
          end if
          if (element > 1) then
             if (any(element_var(element) == element_var(1:element-1))) then
                ierr = ramses_bridge_err_index
                return
             end if
          end if
       end do
    end if
    if (.not. all(ieee_is_finite(cell_volume_code)) .or. &
         .not. all(ieee_is_finite(weights)) .or. &
         any(cell_volume_code <= 0.0_stellar_dp) .or. &
         any(weights < 0.0_stellar_dp)) then
       ierr = ramses_bridge_err_argument
       return
    end if
    weight_sum = sum(weights)
    if (.not. ieee_is_finite(weight_sum) .or. weight_sum <= 0.0_stellar_dp) then
       ierr = ramses_bridge_err_argument
       return
    end if
    ejected_sum = sum(budget%ejected_mass)
    generic_metal_mass = generic_metal_ejecta_mass(budget%returned_mass, &
         budget%ejected_mass)
    if (.not. ieee_is_finite(ejected_sum) .or. &
         .not. ieee_is_finite(generic_metal_mass) .or. &
         ejected_sum > budget%returned_mass + tol * &
         max(1.0_stellar_dp, budget%returned_mass)) then
       ierr = ramses_bridge_err_snia
       return
    end if
    if (.not. present(element_var) .and. ejected_sum > tol * &
         max(1.0_stellar_dp, budget%returned_mass)) then
       ! Never silently discard a selected SNIa chemical payload.
       ierr = ramses_bridge_err_index
       return
    end if
    if (.not. present(total_metal_var) .and. generic_metal_mass > tol * &
         max(1.0_stellar_dp, budget%returned_mass)) then
       ierr = ramses_bridge_err_index
       return
    end if

    allocate(increments(n_cells), stat=allocation_status)
    if (allocation_status /= 0) then
       ierr = ramses_bridge_err_result
       return
    end if

    ! Build every cell increment before modifying uold.  This makes the
    ! bridge transactional when a policy or a budget fails at any cell.
    density_scale = scale_density_cgs
    momentum_scale = scale_density_cgs * scale_velocity_cgs
    energy_scale = scale_density_cgs * scale_velocity_cgs**2
    do cell = 1, n_cells
       volume_cm3 = cell_volume_code(cell) * scale_length_cgs**3
       if (.not. ieee_is_finite(volume_cm3) .or. volume_cm3 <= 0.0_stellar_dp) then
          deallocate(increments)
          ierr = ramses_bridge_err_result
          return
       end if
       call build_snia_cell_increment(budget, volume_cm3, &
            bulk_velocity_cm_s(:,cell), coupling, increments(cell), increment_ierr)
       if (increment_ierr /= snia_deposition_ok) then
          deallocate(increments)
          ierr = ramses_bridge_err_snia
          return
       end if
    end do

    do cell = 1, n_cells
       if (weights(cell) == 0.0_stellar_dp) cycle
       normalized_weight = weights(cell) / weight_sum
       uold(density_var,cell) = uold(density_var,cell) + &
            normalized_weight * increments(cell)%mass_density / density_scale
       uold(energy_var,cell) = uold(energy_var,cell) + &
            normalized_weight * increments(cell)%total_energy_density / energy_scale
       uold(momentum_var(1),cell) = uold(momentum_var(1),cell) + &
            normalized_weight * increments(cell)%momentum_density(1) / momentum_scale
       uold(momentum_var(2),cell) = uold(momentum_var(2),cell) + &
            normalized_weight * increments(cell)%momentum_density(2) / momentum_scale
       uold(momentum_var(3),cell) = uold(momentum_var(3),cell) + &
            normalized_weight * increments(cell)%momentum_density(3) / momentum_scale
       if (present(total_metal_var)) then
          uold(total_metal_var,cell) = uold(total_metal_var,cell) + &
               normalized_weight * increments(cell)%total_metal_density / density_scale
       end if
       if (present(element_var)) then
          do element = 1, n_stellar_elements
             if (element_var(element) == 0) cycle
             uold(element_var(element),cell) = uold(element_var(element),cell) + &
                  normalized_weight * increments(cell)%element_mass_density(element) / &
                  density_scale
          end do
       end if
    end do
    deallocate(increments)
  end subroutine deposit_snia_budget_to_uold

  subroutine build_snia_budget_unew_delta(budget, coupling, bulk_velocity_cm_s, &
       scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, &
       cell_volume_code, density_var, energy_var, momentum_var, delta, &
       tolerance, ierr, total_metal_var, element_var)
    ! Build one SNIa contribution in the RAMSES row-major layout without
    ! touching a production array.  The caller owns the eventual synchronized
    ! commit.  The variable-major scratch is retained only as a unit-conversion
    ! boundary shared with deposit_snia_budget_to_uold.
    type(snia_event_budget_t), intent(in) :: budget
    type(snia_thermal_coupling_t), intent(in) :: coupling
    real(stellar_dp), intent(in) :: bulk_velocity_cm_s(3)
    real(stellar_dp), intent(in) :: scale_length_cgs, scale_density_cgs
    real(stellar_dp), intent(in) :: scale_velocity_cgs
    integer, intent(in) :: nvar
    real(stellar_dp), intent(in) :: cell_volume_code
    integer, intent(in) :: density_var, energy_var, momentum_var(3)
    real(stellar_dp), intent(out) :: delta(:)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    integer, intent(in), optional :: total_metal_var
    integer, intent(in), optional :: element_var(n_stellar_elements)

    real(stellar_dp), allocatable :: scratch_uold(:,:)
    real(stellar_dp) :: velocity_cells(3,1), volume_cells(1), weights(1)
    integer :: allocation_status, scratch_ierr

    ierr = ramses_bridge_ok
    delta = 0.0_stellar_dp
    if (nvar <= 0 .or. size(delta) < nvar) then
       ierr = ramses_bridge_err_argument
       return
    end if
    allocate(scratch_uold(nvar,1), stat=allocation_status)
    if (allocation_status /= 0) then
       ierr = ramses_bridge_err_result
       return
    end if
    scratch_uold = 0.0_stellar_dp
    velocity_cells(:,1) = bulk_velocity_cm_s
    volume_cells(1) = cell_volume_code
    weights(1) = 1.0_stellar_dp

    if (present(total_metal_var)) then
       if (present(element_var)) then
          call deposit_snia_budget_to_uold(budget, coupling, velocity_cells, &
               scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, 1, &
               volume_cells, weights, density_var, energy_var, momentum_var, &
               scratch_uold, tolerance, scratch_ierr, total_metal_var, element_var)
       else
          call deposit_snia_budget_to_uold(budget, coupling, velocity_cells, &
               scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, 1, &
               volume_cells, weights, density_var, energy_var, momentum_var, &
               scratch_uold, tolerance, scratch_ierr, total_metal_var)
       end if
    else if (present(element_var)) then
       call deposit_snia_budget_to_uold(budget, coupling, velocity_cells, &
            scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, 1, &
            volume_cells, weights, density_var, energy_var, momentum_var, &
            scratch_uold, tolerance, scratch_ierr, element_var=element_var)
    else
       call deposit_snia_budget_to_uold(budget, coupling, velocity_cells, &
            scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, 1, &
            volume_cells, weights, density_var, energy_var, momentum_var, &
            scratch_uold, tolerance, scratch_ierr)
    end if
    if (scratch_ierr /= ramses_bridge_ok) then
       deallocate(scratch_uold)
       ierr = scratch_ierr
       return
    end if
    delta(1:nvar) = scratch_uold(:,1)
    deallocate(scratch_uold)
    if (.not. all(ieee_is_finite(delta(1:nvar)))) then
       delta = 0.0_stellar_dp
       ierr = ramses_bridge_err_result
    end if
  end subroutine build_snia_budget_unew_delta

  subroutine deposit_snia_budget_to_unew(budget, coupling, bulk_velocity_cm_s, &
       scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, &
       n_local_cells, n_targets, target_cell, owner_rank, local_rank, &
       cell_volume_code, weights, density_var, energy_var, momentum_var, &
       unew, tolerance, ierr, total_metal_var, element_var)
    ! RAMSES stores conserved state as unew(local_cell,variable), whereas the
    ! unit bridge above uses variable-major scratch storage.  This adapter
    ! validates the AMR/MPI-selected target list, deposits transactionally in
    ! scratch storage, and scatters to addressable local rows.  RAMSES
    ! virtual/reception rows are legal targets and are reconciled later by
    ! make_virtual_reverse_dp.  owner_rank is retained as provenance metadata;
    ! this helper does not claim cross-rank atomicity.
    type(snia_event_budget_t), intent(in) :: budget
    type(snia_thermal_coupling_t), intent(in) :: coupling
    real(stellar_dp), intent(in) :: bulk_velocity_cm_s(3,n_targets)
    real(stellar_dp), intent(in) :: scale_length_cgs, scale_density_cgs
    real(stellar_dp), intent(in) :: scale_velocity_cgs
    integer, intent(in) :: nvar, n_local_cells, n_targets
    integer, intent(in) :: target_cell(n_targets), owner_rank(n_targets)
    integer, intent(in) :: local_rank
    real(stellar_dp), intent(in) :: cell_volume_code(n_targets)
    real(stellar_dp), intent(in) :: weights(n_targets)
    integer, intent(in) :: density_var, energy_var, momentum_var(3)
    real(stellar_dp), intent(inout) :: unew(n_local_cells,nvar)
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    integer, intent(in), optional :: total_metal_var
    integer, intent(in), optional :: element_var(n_stellar_elements)

    real(stellar_dp), allocatable :: selected_uold(:,:)
    integer :: target, previous, scratch_ierr, allocation_status

    ierr = ramses_bridge_ok
    if (n_local_cells <= 0 .or. n_targets <= 0 .or. local_rank < 0) then
       ierr = ramses_bridge_err_target
       return
    end if
    do target = 1, n_targets
       if (target_cell(target) < 1 .or. target_cell(target) > n_local_cells .or. &
            owner_rank(target) < 0) then
          ierr = ramses_bridge_err_target
          return
       end if
       if (target > 1) then
          do previous = 1, target - 1
             if (target_cell(target) == target_cell(previous)) then
                ierr = ramses_bridge_err_target
                return
             end if
          end do
       end if
    end do

    allocate(selected_uold(nvar,n_targets), stat=allocation_status)
    if (allocation_status /= 0) then
       ierr = ramses_bridge_err_result
       return
    end if
    selected_uold = 0.0_stellar_dp
    call deposit_snia_budget_to_uold(budget, coupling, bulk_velocity_cm_s, &
         scale_length_cgs, scale_density_cgs, scale_velocity_cgs, nvar, &
         n_targets, cell_volume_code, weights, density_var, energy_var, &
         momentum_var, selected_uold, tolerance, scratch_ierr, &
         total_metal_var, element_var)
    if (scratch_ierr /= ramses_bridge_ok) then
       deallocate(selected_uold)
       ierr = scratch_ierr
       return
    end if

    do target = 1, n_targets
       if (weights(target) == 0.0_stellar_dp) cycle
       !$omp atomic update
       unew(target_cell(target),density_var) = &
            unew(target_cell(target),density_var) + selected_uold(density_var,target)
       !$omp atomic update
       unew(target_cell(target),energy_var) = &
            unew(target_cell(target),energy_var) + selected_uold(energy_var,target)
       !$omp atomic update
       unew(target_cell(target),momentum_var(1)) = &
            unew(target_cell(target),momentum_var(1)) + selected_uold(momentum_var(1),target)
       !$omp atomic update
       unew(target_cell(target),momentum_var(2)) = &
            unew(target_cell(target),momentum_var(2)) + selected_uold(momentum_var(2),target)
       !$omp atomic update
       unew(target_cell(target),momentum_var(3)) = &
            unew(target_cell(target),momentum_var(3)) + selected_uold(momentum_var(3),target)
       if (present(total_metal_var)) then
          !$omp atomic update
          unew(target_cell(target),total_metal_var) = &
               unew(target_cell(target),total_metal_var) + &
               selected_uold(total_metal_var,target)
       end if
       if (present(element_var)) then
          do previous = 1, n_stellar_elements
             if (element_var(previous) == 0) cycle
             !$omp atomic update
             unew(target_cell(target),element_var(previous)) = &
                  unew(target_cell(target),element_var(previous)) + &
                  selected_uold(element_var(previous),target)
          end do
       end if
    end do
    deallocate(selected_uold)
  end subroutine deposit_snia_budget_to_unew

end module stellar_ramses_bridge

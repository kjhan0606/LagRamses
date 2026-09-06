! Phase 0 interpolation layer for the common stellar-yield table.
!
! Interpolation is performed only when all required grid corners exist.  A
! missing corner is reported to the caller rather than silently replaced by
! a nearest-neighbor or extrapolated value.  This is important for the AGB
! tables, where the metallicity and mass grids are not necessarily identical
! to the massive-star grids.

module stellar_yield_interpolation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels
  use stellar_yield_tables, only: stellar_yield_table_t, &
       yield_mass_assignment_linear, yield_mass_assignment_piecewise_constant
  implicit none

  private

  integer, parameter, public :: interpolation_ok = 0
  integer, parameter, public :: interpolation_err_table = 1
  integer, parameter, public :: interpolation_err_channel = 2
  integer, parameter, public :: interpolation_err_grid = 3
  integer, parameter, public :: interpolation_err_argument = 4
  integer, parameter, public :: interpolation_err_nonfinite = 5
  integer, parameter, public :: interpolation_err_assignment_mode = 6

  public :: interpolate_yield_row

contains

  subroutine interpolate_yield_row(table, channel_id, query_mass, query_z, &
       query_age_gyr, returned_mass, remnant_mass, energy, momentum, &
       ejected_mass, net_yield, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: query_mass, query_z, query_age_gyr
    real(stellar_dp), intent(out) :: returned_mass, remnant_mass, energy
    real(stellar_dp), intent(out) :: momentum(3)
    real(stellar_dp), intent(out) :: ejected_mass(n_stellar_elements)
    real(stellar_dp), intent(out) :: net_yield(n_stellar_elements)
    integer, intent(out) :: ierr

    real(stellar_dp) :: mass_nodes(2), z_nodes(2), age_nodes(2)
    real(stellar_dp) :: mass_weights(2), z_weights(2), age_weights(2)
    real(stellar_dp) :: mass_lo, mass_hi, z_lo, z_hi, age_lo, age_hi
    real(stellar_dp) :: weight
    real(stellar_dp) :: mass_value, z_value, age_value
    real(stellar_dp) :: eps
    integer :: n_mass_nodes, n_z_nodes, n_age_nodes
    integer :: im, iz, ia, row
    logical :: found

    ierr = interpolation_ok
    returned_mass = 0.0_stellar_dp
    remnant_mass = 0.0_stellar_dp
    energy = 0.0_stellar_dp
    momentum = 0.0_stellar_dp
    ejected_mass = 0.0_stellar_dp
    net_yield = 0.0_stellar_dp

    if (.not. table%loaded .or. table%n_rows <= 0) then
       ierr = interpolation_err_table
       return
    end if
    if (channel_id < 1 .or. channel_id > n_stellar_channels) then
       ierr = interpolation_err_channel
       return
    end if
    if (.not. ieee_is_finite(query_mass) .or. &
         .not. ieee_is_finite(query_z) .or. &
         .not. ieee_is_finite(query_age_gyr)) then
       ierr = interpolation_err_argument
       return
    end if
    if (query_mass <= 0.0_stellar_dp .or. query_z < 0.0_stellar_dp .or. &
         query_age_gyr < 0.0_stellar_dp) then
       ierr = interpolation_err_argument
       return
    end if
    if (table%mass_assignment_mode /= yield_mass_assignment_linear .and. &
         table%mass_assignment_mode /= yield_mass_assignment_piecewise_constant) then
       ierr = interpolation_err_assignment_mode
       return
    end if

    call find_bounds(table, channel_id, 1, query_mass, mass_lo, mass_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if
    call find_bounds(table, channel_id, 2, query_z, z_lo, z_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if
    call find_bounds(table, channel_id, 3, query_age_gyr, age_lo, age_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if

    if (table%mass_assignment_mode == yield_mass_assignment_piecewise_constant) then
       call make_piecewise_mass_node(mass_lo, mass_hi, query_mass, mass_nodes, &
            mass_weights, n_mass_nodes)
    else
       call make_nodes(mass_lo, mass_hi, query_mass, mass_nodes, mass_weights, &
            n_mass_nodes)
    end if
    call make_nodes(z_lo, z_hi, query_z, z_nodes, z_weights, n_z_nodes)
    call make_nodes(age_lo, age_hi, query_age_gyr, age_nodes, age_weights, &
         n_age_nodes)

    eps = 1.0e-10_stellar_dp
    do im = 1, n_mass_nodes
       do iz = 1, n_z_nodes
          do ia = 1, n_age_nodes
             mass_value = mass_nodes(im)
             z_value = z_nodes(iz)
             age_value = age_nodes(ia)
             row = find_grid_row(table, channel_id, mass_value, z_value, &
                  age_value)
             if (row < 1) then
                ierr = interpolation_err_grid
                return
             end if

             weight = mass_weights(im) * z_weights(iz) * age_weights(ia)
             if (abs(weight) <= eps) cycle
             returned_mass = returned_mass + weight * table%returned_mass(row)
             remnant_mass = remnant_mass + weight * table%remnant_mass(row)
             energy = energy + weight * table%energy(row)
             momentum = momentum + weight * table%momentum(row,:)
             ejected_mass = ejected_mass + weight * &
                  table%ejected_mass(row,:)
             net_yield = net_yield + weight * table%net_yield(row,:)
          end do
       end do
    end do

    if (.not. finite_result(returned_mass, remnant_mass, energy, momentum, &
         ejected_mass, net_yield)) then
       ierr = interpolation_err_nonfinite
       returned_mass = 0.0_stellar_dp
       remnant_mass = 0.0_stellar_dp
       energy = 0.0_stellar_dp
       momentum = 0.0_stellar_dp
       ejected_mass = 0.0_stellar_dp
       net_yield = 0.0_stellar_dp
    end if
  end subroutine interpolate_yield_row

  subroutine find_bounds(table, channel_id, axis, query, lower, upper, found)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id, axis
    real(stellar_dp), intent(in) :: query
    real(stellar_dp), intent(out) :: lower, upper
    logical, intent(out) :: found

    integer :: i
    real(stellar_dp) :: value
    logical :: has_lower, has_upper

    lower = 0.0_stellar_dp
    upper = 0.0_stellar_dp
    has_lower = .false.
    has_upper = .false.
    found = .false.

    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       value = coordinate_value(table, i, axis)
       found = .true.
       if (value <= query) then
          if (.not. has_lower .or. value > lower) lower = value
          has_lower = .true.
       end if
       if (value >= query) then
          if (.not. has_upper .or. value < upper) upper = value
          has_upper = .true.
       end if
    end do

    ! Out-of-domain requests are hard errors.  In particular, do not clamp to
    ! an endpoint: doing so turns a missing late-time or early-time table row
    ! into a physically different source history.
    if (.not. found .or. .not. has_lower .or. .not. has_upper) then
       found = .false.
       lower = 0.0_stellar_dp
       upper = 0.0_stellar_dp
    end if
  end subroutine find_bounds

  subroutine make_nodes(lower, upper, query, nodes, weights, n_nodes)
    real(stellar_dp), intent(in) :: lower, upper, query
    real(stellar_dp), intent(out) :: nodes(2), weights(2)
    integer, intent(out) :: n_nodes
    real(stellar_dp) :: fraction
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp

    nodes = 0.0_stellar_dp
    weights = 0.0_stellar_dp
    if (abs(upper - lower) <= tolerance * max(1.0_stellar_dp, &
         abs(lower), abs(upper))) then
       n_nodes = 1
       nodes(1) = lower
       weights(1) = 1.0_stellar_dp
       return
    end if

    n_nodes = 2
    nodes(1) = lower
    nodes(2) = upper
    fraction = (query - lower) / (upper - lower)
    fraction = max(0.0_stellar_dp, min(1.0_stellar_dp, fraction))
    weights(1) = 1.0_stellar_dp - fraction
    weights(2) = fraction
  end subroutine make_nodes

  subroutine make_piecewise_mass_node(lower, upper, query, nodes, weights, n_nodes)
    ! A source-node fate is not a quantity that can be linearly blended with
    ! the neighboring node.  Select the left node for an interior half-open
    ! cell and the exact node at a grid edge.  Z and age remain handled by
    ! their ordinary table policy; source-node callers must use exact values
    ! on those axes when their fate semantics are discrete.
    real(stellar_dp), intent(in) :: lower, upper, query
    real(stellar_dp), intent(out) :: nodes(2), weights(2)
    integer, intent(out) :: n_nodes
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp

    nodes = 0.0_stellar_dp
    weights = 0.0_stellar_dp
    n_nodes = 1
    if (abs(upper - lower) <= tolerance * max(1.0_stellar_dp, &
         abs(lower), abs(upper))) then
       nodes(1) = lower
    else if (query >= upper) then
       nodes(1) = upper
    else
       nodes(1) = lower
    endif
    weights(1) = 1.0_stellar_dp
  end subroutine make_piecewise_mass_node

  integer function find_grid_row(table, channel_id, mass, metallicity, age_gyr)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: mass, metallicity, age_gyr
    integer :: i

    find_grid_row = 0
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       if (.not. same_coordinate(table%initial_mass(i), mass)) cycle
       if (.not. same_coordinate(table%birth_metallicity(i), metallicity)) cycle
       if (.not. same_coordinate(table%age_gyr(i), age_gyr)) cycle
       find_grid_row = i
       return
    end do
  end function find_grid_row

  real(stellar_dp) function coordinate_value(table, row, axis)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: row, axis

    select case (axis)
    case (1)
       coordinate_value = table%initial_mass(row)
    case (2)
       coordinate_value = table%birth_metallicity(row)
    case (3)
       coordinate_value = table%age_gyr(row)
    case default
       coordinate_value = 0.0_stellar_dp
    end select
  end function coordinate_value

  logical function same_coordinate(a, b)
    real(stellar_dp), intent(in) :: a, b
    real(stellar_dp) :: scale

    scale = max(1.0_stellar_dp, abs(a), abs(b))
    same_coordinate = abs(a - b) <= 1.0e-10_stellar_dp * scale
  end function same_coordinate

  logical function finite_result(returned_mass, remnant_mass, energy, momentum, &
       ejected_mass, net_yield)
    real(stellar_dp), intent(in) :: returned_mass, remnant_mass, energy
    real(stellar_dp), intent(in) :: momentum(3)
    real(stellar_dp), intent(in) :: ejected_mass(n_stellar_elements)
    real(stellar_dp), intent(in) :: net_yield(n_stellar_elements)
    integer :: i

    finite_result = ieee_is_finite(returned_mass) .and. &
         ieee_is_finite(remnant_mass) .and. ieee_is_finite(energy)
    do i = 1, 3
       finite_result = finite_result .and. ieee_is_finite(momentum(i))
    end do
    do i = 1, n_stellar_elements
       finite_result = finite_result .and. &
            ieee_is_finite(ejected_mass(i)) .and. ieee_is_finite(net_yield(i))
    end do
  end function finite_result

end module stellar_yield_interpolation

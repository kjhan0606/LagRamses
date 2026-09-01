! Contract audit routines for canonical stellar-yield tables.
!
! These checks do not replace a scientific comparison with the source papers.
! They do enforce the numerical contract used by the runtime: finite physical
! values, non-negative cumulative material/energy, mass closure, monotonic
! cumulative histories, unique coordinates, and (when requested) a complete
! Cartesian mass-metallicity-age grid for the required channels.

module stellar_yield_audit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, channel_wind, channel_snii
  use stellar_yield_tables, only: stellar_yield_table_t
  implicit none

  private
  integer, parameter, public :: yield_audit_ok = 0
  integer, parameter, public :: yield_audit_err_table = 1
  integer, parameter, public :: yield_audit_err_value = 2
  integer, parameter, public :: yield_audit_err_mass = 4
  integer, parameter, public :: yield_audit_err_monotonic = 8
  integer, parameter, public :: yield_audit_err_nonfinite = 16
  integer, parameter, public :: yield_audit_err_duplicate = 32
  integer, parameter, public :: yield_audit_err_grid = 64
  integer, parameter, public :: yield_audit_err_energy_monotonic = 128
  integer, parameter, public :: yield_audit_err_remnant_ownership = 256

  public :: audit_yield_table

contains

  subroutine audit_yield_table(table, tolerance, ierr, require_complete, &
       terminal_remnant_owner)
    type(stellar_yield_table_t), intent(in) :: table
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    logical, intent(in), optional :: require_complete
    logical, intent(in), optional :: terminal_remnant_owner(:)

    real(stellar_dp) :: tol, ejected_sum, scale
    logical :: require_grid, row_is_finite, channel_is_bad
    integer :: i, j, channel

    ierr = yield_audit_ok
    require_grid = .false.
    if (present(require_complete)) require_grid = require_complete
    tol = max(tolerance, 1.0e-12_stellar_dp)

    if (.not. table%loaded .or. table%n_rows <= 0) then
       ierr = yield_audit_err_table
       return
    end if
    if (.not. allocated(table%channel) .or. &
         .not. allocated(table%initial_mass) .or. &
         .not. allocated(table%birth_metallicity) .or. &
         .not. allocated(table%age_gyr) .or. &
         .not. allocated(table%returned_mass) .or. &
         .not. allocated(table%remnant_mass) .or. &
         .not. allocated(table%energy) .or. .not. allocated(table%momentum) &
         .or. .not. allocated(table%ejected_mass) .or. &
         .not. allocated(table%net_yield)) then
       ierr = yield_audit_err_table
       return
    end if
    if (table%n_rows > size(table%channel) .or. &
         table%n_rows > size(table%initial_mass) .or. &
         table%n_rows > size(table%birth_metallicity) .or. &
         table%n_rows > size(table%age_gyr) .or. &
         table%n_rows > size(table%returned_mass) .or. &
         table%n_rows > size(table%remnant_mass) .or. &
         table%n_rows > size(table%energy) .or. &
         table%n_rows > size(table%momentum, 1) .or. &
         table%n_rows > size(table%ejected_mass, 1) .or. &
         table%n_rows > size(table%net_yield, 1)) then
       ierr = yield_audit_err_table
       return
    end if
    if (present(terminal_remnant_owner)) then
       if (size(terminal_remnant_owner) /= n_stellar_channels) then
          ierr = yield_audit_err_table
          return
       end if
    end if

    do i = 1, table%n_rows
       row_is_finite = finite_row(table, i)
       if (.not. row_is_finite) then
          ierr = ior(ierr, yield_audit_err_nonfinite)
          cycle
       end if

       if (table%channel(i) < 1 .or. table%channel(i) > n_stellar_channels .or. &
            table%initial_mass(i) <= 0.0_stellar_dp .or. &
            table%birth_metallicity(i) < 0.0_stellar_dp .or. &
            table%age_gyr(i) < 0.0_stellar_dp .or. &
            table%returned_mass(i) < -tol .or. &
            table%remnant_mass(i) < -tol .or. &
            table%energy(i) < -tol .or. &
            minval(table%ejected_mass(i,:)) < -tol) then
          ierr = ior(ierr, yield_audit_err_value)
       end if

       ejected_sum = sum(table%ejected_mass(i,:))
       scale = max(1.0_stellar_dp, abs(table%returned_mass(i)), &
            abs(ejected_sum))
       if (ejected_sum > table%returned_mass(i) + tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if

       scale = max(1.0_stellar_dp, table%initial_mass(i))
       if (table%returned_mass(i) + table%remnant_mass(i) > &
            table%initial_mass(i) + tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if
       if (present(terminal_remnant_owner) .and. &
            table%channel(i) >= 1 .and. table%channel(i) <= n_stellar_channels) then
          if (.not. terminal_remnant_owner(table%channel(i)) .and. &
               table%remnant_mass(i) > tol) then
             ierr = ior(ierr, yield_audit_err_remnant_ownership)
          end if
       end if
    end do

    ! A duplicate coordinate would make the interpolation result depend on
    ! row order.  Reject it even when the duplicate values happen to agree.
    do i = 1, table%n_rows
       do j = i + 1, table%n_rows
          if (table%channel(i) /= table%channel(j)) cycle
          if (.not. same_value(table%initial_mass(i), &
               table%initial_mass(j), tol)) cycle
          if (.not. same_value(table%birth_metallicity(i), &
               table%birth_metallicity(j), tol)) cycle
          if (.not. same_value(table%age_gyr(i), table%age_gyr(j), tol)) cycle
          ierr = ior(ierr, yield_audit_err_duplicate)
       end do
    end do

    ! Compare rows on the same channel, mass, and metallicity grid line.
    ! Actual cumulative material, returned mass, and injected energy must not
    ! decrease with age.  Net yields and momentum are allowed to be signed.
    do i = 1, table%n_rows
       do j = 1, table%n_rows
          if (i == j) cycle
          if (table%channel(i) /= table%channel(j)) cycle
          if (.not. same_value(table%initial_mass(i), &
               table%initial_mass(j), tol)) cycle
          if (.not. same_value(table%birth_metallicity(i), &
               table%birth_metallicity(j), tol)) cycle
          if (table%age_gyr(i) >= table%age_gyr(j) - tol) cycle

          scale = max(1.0_stellar_dp, abs(table%returned_mass(j)))
          if (table%returned_mass(i) > table%returned_mass(j) + tol * scale) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, maxval(abs(table%ejected_mass(j,:))))
          if (any(table%ejected_mass(i,:) > table%ejected_mass(j,:) + &
               tol * scale)) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, &
               abs(table%returned_mass(i) - sum(table%ejected_mass(i,:))), &
               abs(table%returned_mass(j) - sum(table%ejected_mass(j,:))))
          if (table%returned_mass(i) - sum(table%ejected_mass(i,:)) > &
               table%returned_mass(j) - sum(table%ejected_mass(j,:)) + &
               tol * scale) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, abs(table%energy(j)))
          if (table%energy(i) > table%energy(j) + tol * scale) then
             ierr = ior(ierr, yield_audit_err_energy_monotonic)
          end if
       end do
    end do

    if (require_grid) then
       do channel = channel_wind, channel_snii
          call audit_complete_channel(table, channel, tol, channel_is_bad)
          if (channel_is_bad) ierr = ior(ierr, yield_audit_err_grid)
       end do
    end if
  end subroutine audit_yield_table

  subroutine audit_complete_channel(table, channel_id, tolerance, bad)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: tolerance
    logical, intent(out) :: bad

    real(stellar_dp), allocatable :: masses(:), metallicities(:), ages(:)
    integer :: i, j, k, status, n_rows, n_mass, n_z, n_age
    integer(kind=8) :: expected_rows

    bad = .false.
    n_rows = count(table%channel(1:table%n_rows) == channel_id)
    if (n_rows <= 0) then
       bad = .true.
       return
    end if

    allocate(masses(n_rows), metallicities(n_rows), ages(n_rows), stat=status)
    if (status /= 0) then
       bad = .true.
       return
    end if
    masses = 0.0_stellar_dp
    metallicities = 0.0_stellar_dp
    ages = 0.0_stellar_dp
    n_mass = 0
    n_z = 0
    n_age = 0
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       call append_unique(table%initial_mass(i), masses, n_mass, tolerance)
       call append_unique(table%birth_metallicity(i), metallicities, n_z, &
            tolerance)
       call append_unique(table%age_gyr(i), ages, n_age, tolerance)
    end do

    expected_rows = int(n_mass, kind=8) * int(n_z, kind=8) * int(n_age, kind=8)
    if (int(n_rows, kind=8) /= expected_rows) bad = .true.

    do i = 1, n_mass
       do j = 1, n_z
          do k = 1, n_age
             if (.not. grid_row_exists(table, channel_id, masses(i), &
                  metallicities(j), ages(k), tolerance)) bad = .true.
          end do
       end do
    end do
    deallocate(masses, metallicities, ages)
  end subroutine audit_complete_channel

  subroutine append_unique(value, values, n_values, tolerance)
    real(stellar_dp), intent(in) :: value, tolerance
    real(stellar_dp), intent(inout) :: values(:)
    integer, intent(inout) :: n_values
    integer :: i

    do i = 1, n_values
       if (same_value(values(i), value, tolerance)) return
    end do
    n_values = n_values + 1
    values(n_values) = value
  end subroutine append_unique

  logical function grid_row_exists(table, channel_id, mass, metallicity, age, &
       tolerance)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: mass, metallicity, age, tolerance
    integer :: i

    grid_row_exists = .false.
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       if (.not. same_value(table%initial_mass(i), mass, &
            tolerance)) cycle
       if (.not. same_value(table%birth_metallicity(i), metallicity, &
            tolerance)) cycle
       if (.not. same_value(table%age_gyr(i), age, tolerance)) cycle
       grid_row_exists = .true.
       return
    end do
  end function grid_row_exists

  logical function finite_row(table, row)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: row
    integer :: i

    finite_row = ieee_is_finite(table%initial_mass(row)) .and. &
         ieee_is_finite(table%birth_metallicity(row)) .and. &
         ieee_is_finite(table%age_gyr(row)) .and. &
         ieee_is_finite(table%returned_mass(row)) .and. &
         ieee_is_finite(table%remnant_mass(row)) .and. &
         ieee_is_finite(table%energy(row))
    do i = 1, 3
       finite_row = finite_row .and. ieee_is_finite(table%momentum(row,i))
    end do
    do i = 1, n_stellar_elements
       finite_row = finite_row .and. &
            ieee_is_finite(table%ejected_mass(row,i)) .and. &
            ieee_is_finite(table%net_yield(row,i))
    end do
  end function finite_row

  logical function same_value(a, b, tolerance)
    real(stellar_dp), intent(in) :: a, b, tolerance
    real(stellar_dp) :: scale

    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
       same_value = .false.
       return
    end if
    scale = max(1.0_stellar_dp, abs(a), abs(b))
    same_value = abs(a - b) <= tolerance * scale
  end function same_value

end module stellar_yield_audit

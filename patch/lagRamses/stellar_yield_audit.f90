! Phase 0 audit routines for canonical stellar-yield tables.
!
! These checks are data-contract checks, not a replacement for comparing a
! yield set with its source paper.  They reject mass-inconsistent rows and
! cumulative tables whose actual ejecta decrease with age.

module stellar_yield_audit
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements
  use stellar_yield_tables, only: stellar_yield_table_t
  implicit none

  private
  integer, parameter, public :: yield_audit_ok = 0
  integer, parameter, public :: yield_audit_err_table = 1
  integer, parameter, public :: yield_audit_err_value = 2
  integer, parameter, public :: yield_audit_err_mass = 4
  integer, parameter, public :: yield_audit_err_monotonic = 8

  public :: audit_yield_table

contains

  subroutine audit_yield_table(table, tolerance, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr

    real(stellar_dp) :: tol, ejected_sum, scale
    integer :: i, j

    ierr = yield_audit_ok
    tol = max(tolerance, 1.0e-12_stellar_dp)
    if (.not. table%loaded .or. table%n_rows <= 0) then
       ierr = yield_audit_err_table
       return
    end if

    do i = 1, table%n_rows
       if (table%channel(i) < 1 .or. table%channel(i) > 5 .or. &
            table%initial_mass(i) <= 0.0_stellar_dp .or. &
            table%birth_metallicity(i) < 0.0_stellar_dp .or. &
            table%age(i) < 0.0_stellar_dp .or. &
            table%returned_mass(i) < -tol .or. &
            table%remnant_mass(i) < -tol) then
          ierr = ior(ierr, yield_audit_err_value)
       end if

       if (minval(table%ejected_mass(i,:)) < -tol) then
          ierr = ior(ierr, yield_audit_err_value)
       end if

       ejected_sum = sum(table%ejected_mass(i,:))
       scale = max(1.0_stellar_dp, abs(table%returned_mass(i)), &
            abs(ejected_sum))
       if (abs(ejected_sum - table%returned_mass(i)) > tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if

       scale = max(1.0_stellar_dp, table%initial_mass(i))
       if (table%returned_mass(i) + table%remnant_mass(i) > &
            table%initial_mass(i) + tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if
    end do

    ! Compare rows on the same channel, mass, and metallicity grid line.
    ! Actual cumulative ejecta and returned mass must not decrease with age.
    do i = 1, table%n_rows
       do j = 1, table%n_rows
          if (i == j) cycle
          if (table%channel(i) /= table%channel(j)) cycle
          if (.not. same_value(table%initial_mass(i), &
               table%initial_mass(j), tol)) cycle
          if (.not. same_value(table%birth_metallicity(i), &
               table%birth_metallicity(j), tol)) cycle
          if (table%age(i) >= table%age(j) - tol) cycle

          scale = max(1.0_stellar_dp, abs(table%returned_mass(j)))
          if (table%returned_mass(i) > table%returned_mass(j) + tol * scale) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, maxval(abs(table%ejected_mass(j,:))))
          if (any(table%ejected_mass(i,:) > table%ejected_mass(j,:) + &
               tol * scale)) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
       end do
    end do
  end subroutine audit_yield_table

  logical function same_value(a, b, tolerance)
    real(stellar_dp), intent(in) :: a, b, tolerance
    real(stellar_dp) :: scale

    scale = max(1.0_stellar_dp, abs(a), abs(b))
    same_value = abs(a - b) <= tolerance * scale
  end function same_value

end module stellar_yield_audit

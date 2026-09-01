! Phase 0 yield-table backend selector.
!
! The production build may use an embedded generated module by defining
! STELLAR_EMBEDDED_YIELDS.  Development builds can continue to load the
! canonical ASCII table at runtime.  Both paths populate the same table type.

module stellar_yield_backend
  use stellar_yield_tables, only: stellar_yield_table_t, load_yield_table, &
       yield_table_ok
#ifdef STELLAR_EMBEDDED_YIELDS
  use stellar_yield_embedded_data, only: load_embedded_yield_table
#endif
  implicit none

  private
  integer, parameter, public :: backend_ok = 0
  integer, parameter, public :: backend_err_external = 1
  integer, parameter, public :: backend_err_embedded_unavailable = 2

  public :: load_yield_backend

contains

  subroutine load_yield_backend(filename, use_embedded, table, ierr)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: use_embedded
    type(stellar_yield_table_t), intent(inout) :: table
    integer, intent(out) :: ierr
    integer :: table_ierr

    if (use_embedded) then
#ifdef STELLAR_EMBEDDED_YIELDS
       call load_embedded_yield_table(table, table_ierr)
       ierr = table_ierr
#else
       ierr = backend_err_embedded_unavailable
#endif
    else
       call load_yield_table(filename, table, table_ierr)
       if (table_ierr == yield_table_ok) then
          ierr = backend_ok
       else
          ierr = backend_err_external
       end if
    end if
  end subroutine load_yield_backend

end module stellar_yield_backend

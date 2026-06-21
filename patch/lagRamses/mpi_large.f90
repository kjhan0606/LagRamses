module mpi_large
  implicit none

contains

#ifndef WITHOUTMPI
  subroutine mpi_large_isend_dp(buf, ncount, nprops, dest, tag, comm, request, ierr)
    include 'mpif.h'
    real(8), intent(in) :: buf(*)
    integer, intent(in) :: ncount, nprops, dest, tag, comm
    integer, intent(out) :: request, ierr
    integer :: row_type
    integer(8) :: total8

    total8 = int(ncount, 8) * int(nprops, 8)
    if(total8 <= 2147483647_8) then
       call MPI_ISEND(buf(1), int(total8), MPI_DOUBLE_PRECISION, &
            dest, tag, comm, request, ierr)
    else
       call MPI_TYPE_CONTIGUOUS(nprops, MPI_DOUBLE_PRECISION, row_type, ierr)
       call MPI_TYPE_COMMIT(row_type, ierr)
       call MPI_ISEND(buf(1), ncount, row_type, dest, tag, comm, request, ierr)
       call MPI_TYPE_FREE(row_type, ierr)
    end if
  end subroutine mpi_large_isend_dp


  subroutine mpi_large_irecv_dp(buf, ncount, nprops, source, tag, comm, request, ierr)
    include 'mpif.h'
    real(8), intent(out) :: buf(*)
    integer, intent(in) :: ncount, nprops, source, tag, comm
    integer, intent(out) :: request, ierr
    integer :: row_type
    integer(8) :: total8

    total8 = int(ncount, 8) * int(nprops, 8)
    if(total8 <= 2147483647_8) then
       call MPI_IRECV(buf(1), int(total8), MPI_DOUBLE_PRECISION, &
            source, tag, comm, request, ierr)
    else
       call MPI_TYPE_CONTIGUOUS(nprops, MPI_DOUBLE_PRECISION, row_type, ierr)
       call MPI_TYPE_COMMIT(row_type, ierr)
       call MPI_IRECV(buf(1), ncount, row_type, source, tag, comm, request, ierr)
       call MPI_TYPE_FREE(row_type, ierr)
    end if
  end subroutine mpi_large_irecv_dp
#endif

end module mpi_large

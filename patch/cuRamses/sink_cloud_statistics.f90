module sink_cloud_statistics
  use amr_parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  public :: sink_cloud_statistics_reduce

contains

  subroutine sink_cloud_statistics_reduce(nrecord, nsink, ndim, record_sink, &
       record_values, sink_stat, ierr)
    integer, intent(in) :: nrecord, nsink, ndim
    integer, intent(in) :: record_sink(*)
    real(dp), intent(in) :: record_values(:,:)
    real(dp), intent(inout) :: sink_stat(:,:)
    integer, intent(out) :: ierr
    integer, allocatable :: order(:), work(:)
    integer :: i, j, k, left, middle, right, width, sink_id, alloc_error

    ierr = 0
    if (nrecord < 0 .or. nsink < 0 .or. ndim < 1) then
       ierr = 1
       return
    endif
    if (size(sink_stat,1) < nsink .or. size(sink_stat,2) /= 2*ndim+1 .or. &
         size(record_values,1) /= 2*ndim .or. size(record_values,2) < nrecord) then
       ierr = 1
       return
    endif
    do i = 1, nrecord
       if (record_sink(i) < 1 .or. record_sink(i) > nsink) then
          ierr = 2
          return
       endif
    enddo
    if (nrecord == 0) return

    allocate(order(nrecord), work(nrecord), stat=alloc_error)
    if (alloc_error /= 0) then
       ierr = 3
       return
    endif
    do i = 1, nrecord
       order(i) = i
    enddo

    ! Bottom-up mergesort over indices keeps the record arrays compact and
    ! gives a total, input-order-independent key: sink index followed by the
    ! complete coordinate/velocity tuple.
    width = 1
    do while (width < nrecord)
       left = 1
       do while (left <= nrecord)
          middle = min(left + width - 1, nrecord)
          right = min(left + 2*width - 1, nrecord)
          i = left
          j = middle + 1
          k = left
          do while (i <= middle .and. j <= right)
             if (record_less(order(j), order(i), record_sink, record_values, ndim)) then
                work(k) = order(j)
                j = j + 1
             else
                work(k) = order(i)
                i = i + 1
             endif
             k = k + 1
          enddo
          do while (i <= middle)
             work(k) = order(i)
             i = i + 1
             k = k + 1
          enddo
          do while (j <= right)
             work(k) = order(j)
             j = j + 1
             k = k + 1
          enddo
          left = right + 1
       enddo
       order = work
       if (width >= nrecord - width) then
          width = nrecord
       else
          width = 2*width
       endif
    enddo

    do i = 1, nrecord
       j = order(i)
       sink_id = record_sink(j)
       sink_stat(sink_id,2*ndim+1) = sink_stat(sink_id,2*ndim+1) + 1.0_dp
       do k = 1, ndim
          sink_stat(sink_id,k) = sink_stat(sink_id,k) + record_values(k,j)
          sink_stat(sink_id,k+ndim) = sink_stat(sink_id,k+ndim) + &
               record_values(k+ndim,j)
       enddo
    enddo

    deallocate(order, work)
  end subroutine sink_cloud_statistics_reduce


  logical function record_less(a, b, record_sink, record_values, ndim) result(less)
    integer, intent(in) :: a, b, ndim
    integer, intent(in) :: record_sink(*)
    real(dp), intent(in) :: record_values(:,:)
    integer :: field
    real(dp) :: value_a, value_b

    if (record_sink(a) /= record_sink(b)) then
       less = record_sink(a) < record_sink(b)
       return
    endif
    do field = 1, 2*ndim
       value_a = record_values(field,a)
       value_b = record_values(field,b)
       if (real_less(value_a, value_b)) then
          less = .true.
          return
       else if (real_less(value_b, value_a)) then
          less = .false.
          return
       endif
    enddo
    less = .false.
  end function record_less


  logical function real_less(a, b) result(less)
    real(dp), intent(in) :: a, b
    integer(int64) :: bits_a, bits_b
    logical :: nan_a, nan_b

    nan_a = ieee_is_nan(a)
    nan_b = ieee_is_nan(b)
    if (nan_a .neqv. nan_b) then
       less = .not. nan_a
       return
    endif
    if (.not. nan_a) then
       if (a < b) then
          less = .true.
          return
       else if (a > b) then
          less = .false.
          return
       endif
    endif

    ! Numeric equality leaves signed zero ambiguous; NaNs also need a stable
    ! tie-break. Their bit patterns provide a deterministic final key.
    bits_a = transfer(a, bits_a)
    bits_b = transfer(b, bits_b)
    less = bits_a < bits_b
  end function real_less

end module sink_cloud_statistics

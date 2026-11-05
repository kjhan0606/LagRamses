program sink_cloud_statistics_test
  use amr_parameters, only: dp
  use sink_cloud_statistics, only: sink_cloud_statistics_reduce
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none

  integer, parameter :: ndim=3, nsink=4, nrecord=2126
  integer :: sink_a(nrecord), sink_b(nrecord), i, j, field, ierr, failures
  real(dp) :: values_a(2*ndim,nrecord), values_b(2*ndim,nrecord)
  real(dp) :: stats_a(nsink,2*ndim+1), stats_b(nsink,2*ndim+1), value
  integer(int64) :: bits_a, bits_b

  failures = 0
  do i = 1, nrecord
     sink_a(i) = 1 + mod(37*i, nsink)
     do field = 1, 2*ndim
        select case (mod(19*i + 11*field, 9))
        case (0)
           value = 0.125_dp
        case (1)
           value = -0.125_dp
        case (2)
           value = 0.5_dp * epsilon(1.0_dp)
        case (3)
           value = -0.5_dp * epsilon(1.0_dp)
        case default
           value = real(mod(13*i + 7*field, 31)-15,dp) * &
                2.0_dp**(-42-mod(i+field,7))
        end select
        values_a(field,i) = value
     enddo
  enddo

  do i = 1, nrecord
     j = nrecord - i + 1
     sink_b(i) = sink_a(j)
     values_b(:,i) = values_a(:,j)
  enddo

  stats_a = 0.0_dp
  stats_b = 0.0_dp
  call sink_cloud_statistics_reduce(nrecord,nsink,ndim,sink_a,values_a,stats_a,ierr)
  if (ierr /= 0) then
     write(*,'(A,I0)') 'reduction A returned error ', ierr
     failures = failures + 1
  endif
  call sink_cloud_statistics_reduce(nrecord,nsink,ndim,sink_b,values_b,stats_b,ierr)
  if (ierr /= 0) then
     write(*,'(A,I0)') 'reduction B returned error ', ierr
     failures = failures + 1
  endif

  do j = 1, size(stats_a,2)
     do i = 1, size(stats_a,1)
        bits_a = transfer(stats_a(i,j), bits_a)
        bits_b = transfer(stats_b(i,j), bits_b)
        if (bits_a /= bits_b) failures = failures + 1
     enddo
  enddo
  if (sum(stats_a(:,2*ndim+1)) /= real(nrecord,dp)) failures = failures + 1

  if (failures /= 0) then
     write(*,'(A,I0)') 'SINK_CLOUD_STATISTICS_FAIL count=', failures
     stop 1
  endif
  write(*,'(A,I0,A,I0)') 'SINK_CLOUD_STATISTICS_PASS records=', nrecord, &
       ' sinks=', nsink
end program sink_cloud_statistics_test

program test_process_cpu_time
  use omp_lib, only: omp_get_wtime,omp_get_max_threads
  implicit none
  integer,parameter::dp=kind(1.0d0)
  integer::i,nthreads
  integer,parameter::nwork=80000000
  real(dp)::cpu0,cpu1,wall0,wall1,ratio,checksum

  nthreads=omp_get_max_threads()
  if(nthreads/=2)then
     write(*,'(A,I0)') 'FAIL expected two OpenMP threads, got ',nthreads
     error stop 1
  endif

  checksum=0d0
  call cpu_time(cpu0)
  wall0=omp_get_wtime()
!$omp parallel do schedule(static) reduction(+:checksum)
  do i=1,nwork
     checksum=checksum+sqrt(dble(mod(i,97)+1))
  end do
!$omp end parallel do
  wall1=omp_get_wtime()
  call cpu_time(cpu1)

  ratio=(cpu1-cpu0)/max(1d-12,wall1-wall0)
  write(*,'(A,F8.3,A,ES12.3)') &
       'PROCESS_CPU_TIME_RATIO=',ratio,' CHECKSUM=',checksum
  ! The timed load balancer requires CPU_TIME to accumulate CPU consumed by
  ! all OpenMP threads, rather than measuring only the calling thread.  Leave
  ! headroom for scheduler noise while still rejecting thread-local semantics.
  if(ratio<1.5d0.or.ratio>2.3d0)then
     write(*,'(A)') 'FAIL CPU_TIME is not a two-thread process CPU clock'
     error stop 2
  endif
end program test_process_cpu_time

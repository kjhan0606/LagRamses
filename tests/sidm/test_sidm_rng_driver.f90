program test_sidm_rng_driver
  use sidm_test_capture
  implicit none
  integer :: i,j
  do call_number=1,2
     call sidm_scatter(1)
  enddo
  if(any(visits/=1)) error stop 'wrong number of nonempty scatter calls'
  do i=1,4
     do j=1,4
        if(draws(i,1)==draws(j,2)) error stop 'RNG draws recur on the next call'
     enddo
  enddo
  if(any(workers/=1)) error stop 'SIDM scatter used multiple OpenMP workers'
  do i=1,2
     write(*,'(A,I0,4(1X,ES25.17E3))') 'DRAW ',i,draws(:,i)
  enddo
  print '(A)', 'Production SIDM wrapper serial stream: PASS'
end program test_sidm_rng_driver

subroutine sub_sidm_scatter(ilevel,icpu,kgrid,subnump,thread_seeds,thread_id,nthreads)
  use omp_lib,only:omp_get_num_threads
  use amr_commons,only:next
  use random,only:IRandNumSize,ranf
  use sidm_test_capture
  implicit none
  integer,intent(in) :: ilevel,icpu,kgrid,subnump,thread_id,nthreads
  integer,intent(inout) :: thread_seeds(IRandNumSize,0:nthreads-1)
  integer :: i,grid,grid_index
  grid=kgrid
  do grid_index=1,subnump
     ! Grid 1 has no eligible pair; grid 2 consumes four RNG draws. This
     ! reproduces unequal work with the production linked-list partitioner.
     if(grid==2) then
        visits(call_number)=visits(call_number)+1
        workers(call_number)=omp_get_num_threads()
        do i=1,4
           call ranf(thread_seeds(:,thread_id),draws(i,call_number))
        enddo
     endif
     grid=next(grid)
  enddo
end subroutine sub_sidm_scatter

subroutine sidm_init_excited()
  error stop 'unexpected excited-state initialization'
end subroutine sidm_init_excited

! Deterministic dispatcher policy for the prelaunch stream-lease test.
module snrt_runtime_backend
  implicit none
  integer :: test_choice=2,test_threads=1
contains
  subroutine snrt_runtime_mn_policy(choice,sharing,threads,ierr)
    integer,intent(out) :: choice,sharing,threads,ierr
    choice=test_choice;sharing=1;threads=test_threads;ierr=0
  end subroutine
end module

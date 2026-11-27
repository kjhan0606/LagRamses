program test_sidm_collision_subcycle
  implicit none

  integer::invalid,nsub,consistent
  real(kind(1.0d0))::dt_sub,p_remaining,p_sub,remaining_before,last_p
  real(kind(1.0d0))::remaining,interval,elapsed,rate,dt_first,dt_second
  real(kind(1.0d0))::dt_next
  logical::final_fraction_seen
  external sidm_choose_collision_substep
  external sidm_collision_schedule_agrees

  call sidm_collision_schedule_agrees(1.0d0,1.0d0,2.0d0,2.0d0, &
       & 0.1d0,0.1d0,consistent)
  call require(consistent==1,'identical MPI schedule')
  dt_next=nearest(1.0d0,1.0d0)
  call sidm_collision_schedule_agrees(1.0d0,dt_next,2.0d0,2.0d0, &
       & 0.1d0,0.1d0,consistent)
  call require(consistent==0,'one-bit rank endpoint mismatch rejected')

  ! A previously admitted interval must preserve the original one-sweep dt.
  call sidm_choose_collision_substep(4.0d0,5.0d-3,2.0d0,0.1d0, &
       & dt_sub,p_remaining,p_sub,invalid)
  call require(invalid==0,'one-step input accepted')
  call require(abs(dt_sub-4.0d0)<=epsilon(1.0d0),'one-step interval unchanged')
  call require(abs(p_remaining-0.04d0)<1.0d-15,'one-step bound')
  call require(abs(p_sub-p_remaining)<=epsilon(1.0d0), &
       & 'one-step probability unchanged')

  ! The recorded L14 scale of rate x interval must be admitted in bounded
  ! collision-only substeps whose durations add to the original interval.
  interval=2.947497d8
  rate=1.999140d-8
  remaining=interval
  elapsed=0.0d0
  nsub=0
  final_fraction_seen=.false.
  last_p=0.0d0
  do while(remaining>0.0d0)
     nsub=nsub+1
     call require(nsub<=4096,'subcycle cap')
     remaining_before=remaining
     call sidm_choose_collision_substep(remaining,rate,1.0d0,0.1d0, &
          & dt_sub,p_remaining,p_sub,invalid)
     call require(invalid==0,'finite subcycle accepted')
     call require(dt_sub>0.0d0 .and. dt_sub<=remaining,'positive progress')
     call require(p_sub<=0.1d0*(1.0d0-64.0d0*epsilon(1.0d0)), &
          & 'hard probability ceiling')
     last_p=p_sub
     elapsed=elapsed+dt_sub
     if(dt_sub>=remaining_before) then
        final_fraction_seen=(nsub>1)
        remaining=0.0d0
     else
        remaining=remaining-dt_sub
     end if
  end do
  call require(nsub>1,'high-rate interval subcycled')
  call require(final_fraction_seen .and. &
       & last_p<=0.1d0*(1.0d0-64.0d0*epsilon(1.0d0)), &
       & 'final fraction admitted below the hard ceiling')
  call require(abs(elapsed-interval)<=8.0d0*epsilon(interval)*interval, &
       & 'substep duration sum')

  ! A fresh bound after a velocity change must shorten the next substep.
  call sidm_choose_collision_substep(1.0d7,2.0d-8,1.0d0,0.1d0, &
       & dt_first,p_remaining,p_sub,invalid)
  call require(invalid==0,'pre-change rate accepted')
  call sidm_choose_collision_substep(7.5d6,4.0d-8,1.0d0,0.1d0, &
       & dt_second,p_remaining,p_sub,invalid)
  call require(invalid==0,'post-change rate accepted')
  call require(dt_second<dt_first,'increased rate shortens next substep')

  ! Zero rate and zero duration are valid and cannot create a scatter step.
  call sidm_choose_collision_substep(10.0d0,0.0d0,1.0d0,0.1d0, &
       & dt_sub,p_remaining,p_sub,invalid)
  call require(invalid==0 .and. dt_sub>=10.0d0 .and. p_sub<=0.0d0, &
       & 'zero rate')
  call sidm_choose_collision_substep(0.0d0,2.0d-8,1.0d0,0.1d0, &
       & dt_sub,p_remaining,p_sub,invalid)
  call require(invalid==0 .and. dt_sub<=0.0d0 .and. p_sub<=0.0d0, &
       & 'zero duration')

  call sidm_choose_collision_substep(1.0d0,1.0d0,1.0d0,0.0d0, &
       & dt_sub,p_remaining,p_sub,invalid)
  call require(invalid/=0,'invalid probability ceiling rejected')

  write(*,'(A,I0,A,ES12.4)') 'collision-subcycle helper PASS: steps=', &
       & nsub,' interval_s=',interval

contains
  subroutine require(condition,message)
    implicit none
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) then
       write(*,'(A,A)') 'FAIL: ',trim(message)
       error stop 1
    end if
  end subroutine require
end program test_sidm_collision_subcycle

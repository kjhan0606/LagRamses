! All-pair check of the production bound, including odd occupancies and
! unequal particle masses.  Also reject NaN before any ordered comparison.
program test_sidm_rate_bound
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter::nmax=65
  integer::n,i,j,k,invalid,consistent
  real(kind(1d0))::vel(nmax,3),mass(nmax),vmin(3),vmax(3)
  real(kind(1d0))::bound,pair_p,nan,dt,p_remaining,p_sub

  do n=2,nmax
     do i=1,n
        mass(i)=0.3d0+real(mod(17*i,31),kind(1d0))/17d0
        do k=1,3
           vel(i,k)=sin(real(i*k+n,kind(1d0))*0.713d0) &
                +cos(real(i+3*k,kind(1d0))*0.321d0)
        end do
     end do
     vmin=minval(vel(1:n,:),dim=1)
     vmax=maxval(vel(1:n,:),dim=1)
     call sidm_probability_rate_bound(0.13d0,maxval(mass(1:n)),vmin,vmax, &
          4.2d0,n,bound,invalid)
     call require(invalid==0,'valid rate bound')
     do i=1,n-1
        do j=i+1,n
           call sidm_pair_probability(0.13d0,mass(i),mass(j), &
                sqrt(sum((vel(i,:)-vel(j,:))**2)),1d0,4.2d0,n,n/2,pair_p)
           call require(pair_p<=bound,'bound contains every physical pair')
        end do
     end do
  end do

  nan=ieee_value(0d0,ieee_quiet_nan)
  call sidm_probability_rate_bound(nan,1d0,(/0d0,0d0,0d0/), &
       (/1d0,1d0,1d0/),1d0,2,bound,invalid)
  call require(invalid/=0,'NaN cross section rejected')
  call sidm_probability_rate_bound(1d0,1d0,(/0d0,0d0,nan/), &
       (/1d0,1d0,1d0/),1d0,2,bound,invalid)
  call require(invalid/=0,'NaN velocity rejected')
  call sidm_choose_collision_substep(nan,1d0,1d0,0.1d0, &
       dt,p_remaining,p_sub,invalid)
  call require(invalid/=0,'NaN interval rejected')
  call sidm_probability_timestep(1d0,nan,0.1d0,dt,invalid)
  call require(invalid/=0,'NaN time unit rejected')
  call sidm_collision_schedule_agrees(1d0,1d0,nan,1d0,0.1d0,0.1d0,consistent)
  call require(consistent==0,'NaN rank schedule rejected')
  print '(A)', 'All-pair rate bound and nonfinite input: PASS'
contains
  subroutine require(ok,message)
    logical,intent(in)::ok
    character(*),intent(in)::message
    if(.not.ok) then
       write(*,'(A,A)') 'FAIL: ',message
       error stop 1
    end if
  end subroutine require
end program test_sidm_rate_bound

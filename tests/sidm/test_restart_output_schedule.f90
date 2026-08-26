program test_restart_output_schedule
  use amr_parameters, only: dp, restart_output_index
  implicit none

  real(dp)::aout(5),tout(5)

  aout=(/0.20d0,0.25d0,0.33333d0,0.50d0,1.00d0/)
  tout=(/1d0,2d0,3d0,4d0,5d0/)

  call expect('cosmological seed skips completed entry', &
       & restart_output_index(.true.,5,aout,tout,0.20d0,0d0),2)
  call expect('cosmological completed target advances', &
       & restart_output_index(.true.,5,aout,tout,0.33333d0,0d0),4)
  call expect('cosmological in-between target', &
       & restart_output_index(.true.,5,aout,tout,0.40d0,0d0),4)
  call expect('cosmological final target remains bounded', &
       & restart_output_index(.true.,5,aout,tout,1.00d0,0d0),5)
  call expect('non-cosmological schedule', &
       & restart_output_index(.false.,5,aout,tout,0d0,2d0),3)

contains

  subroutine expect(label,actual,expected)
    implicit none
    character(len=*),intent(in)::label
    integer,intent(in)::actual,expected

    if(actual/=expected)then
       write(*,'(A,": expected ",I0,", got ",I0)')trim(label),expected,actual
       error stop 1
    end if
  end subroutine expect

end program test_restart_output_schedule

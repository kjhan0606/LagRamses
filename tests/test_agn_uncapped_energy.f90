program test_agn_uncapped_energy
  use amr_parameters, only: dp
  use agn_feedback_deposition, only: agn_deposit_cell
  implicit none
  real(dp)::row(5),deferred,expected,delta
  integer::i,j,ierr
  do i=-12,12
     do j=1,37
        row=[1d0,0.1d0,0.2d0,0.3d0,1d0]*10d0**i
        delta=real(j,dp)*0.017d0*10d0**i
        expected=row(5)+delta
        call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],delta,1d0, &
             5d0/3d0,1d0,1d6,deferred,ierr)
        if(ierr/=0.or.deferred/=0d0.or.row(5)/=expected) &
             error stop 'uncapped deposit fabricated saved energy'
     enddo
  enddo
  row=[1d0,0d0,0d0,0d0,1d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],10d0,1d0, &
       5d0/3d0,1d0,2d0,deferred,ierr)
  if(ierr/=0.or.deferred<=0d0)error stop 'capped deposit lost saved energy'
  if(abs(row(5)+deferred-11d0)>1d-13)error stop 'capped budget mismatch'
  print *, 'PASS: uncapped exact energy and genuine cap deferral'
end program test_agn_uncapped_energy

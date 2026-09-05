program agn_feedback_deposition_smoke
  use amr_parameters, only: dp
  use agn_feedback_deposition
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  real(dp) :: row(5), before(5), deferred, dm, momentum(3), energy
  real(dp) :: volumes(3), offsets(3), total_mass, total_momentum(3), total_energy
  integer :: ierr, i

  row=[1d0,0d0,0d0,0d0,5d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],3d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==8d0 .and. deferred==0d0,'thermal below cap')
  row=[1d0,0d0,0d0,0d0,5d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],8d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==10d0 .and. deferred==9d0,'cap crossing and nonunit volume')
  call check(abs(row(5)*3d0+deferred-39d0)<1d-12,'gas plus deferred equals old plus input')
  row=[1d0,0d0,0d0,0d0,12d0]
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],4d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==0 .and. row(5)==12d0 .and. deferred==12d0,'already hot gas is not cooled')

  ! A donor with rho=2, v=3, internal energy density 4 loses drho=.5.
  ! average_AGN preserves velocity/internal energy, so its post-removal row
  ! is [1.5,4.5,0,0,10.75]. Return at the captured DONOR velocity, not v_BH.
  row=[1.5d0,4.5d0,0d0,0d0,10.75d0]
  dm=0.5d0
  momentum=[1.5d0,0d0,0d0]
  energy=4d0+0.5d0*dm*3d0**2
  call agn_deposit_cell(row,dm,momentum,energy,2d0,2d0,1d0,100d0,deferred,ierr)
  call check(ierr==0 .and. maxval(abs(row-[2d0,6d0,0d0,0d0,17d0]))<1d-12 .and. &
       deferred==0d0,'moving donor fallback restores mass momentum and kinetic energy')

  ! Synthetic cylinder with unequal cell volumes, including its midplane.
  ! Kernel shape=1 for these axial points; norm=sum(weight*volume)=4.
  volumes=[1d0,2d0,1d0]; offsets=[-1d0,0d0,1d0]
  total_mass=0d0; total_momentum=0d0; total_energy=0d0
  do i=1,3
     call agn_jet_delta(2d0,1d0,sum(volumes),[0d0,0d0,0d0],[1d0,0d0,0d0], &
          sqrt(8d0),offsets(i),dm,momentum,energy)
     if(i==2)call check(all(momentum==0d0) .and. abs(energy-2d0)<1d-12, &
          'midplane retains opposed-stream energy with zero net momentum')
     total_mass=total_mass+dm*volumes(i)
     total_momentum=total_momentum+momentum*volumes(i)
     row=[1d0,0d0,0d0,0d0,1d0]
     call agn_deposit_cell(row,dm,momentum,energy,volumes(i),2d0,1d0,100d0,deferred,ierr)
     call check(ierr==0,'cylinder cell deposition')
     total_energy=total_energy+(row(5)-1d0)*volumes(i)+deferred
  end do
  call check(abs(total_mass-2d0)<1d-12 .and. maxval(abs(total_momentum))<1d-12 .and. &
       abs(total_energy-8d0)<1d-12,'cylinder sums loaded mass and supplied kinetic energy')

  row=[1d0,0d0,0d0,0d0,5d0]; before=row
  call agn_deposit_cell(row,0d0,[0d0,0d0,0d0],ieee_value(0d0,ieee_quiet_nan), &
       3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr/=0 .and. all(row==before) .and. deferred==0d0,'invalid input is non-mutating')
  row=[1d0,3d0,0d0,0d0,1d0]; before=row
  call agn_deposit_cell(row,0.5d0,[1.5d0,0d0,0d0],8d0,3d0,2d0,1d0,10d0,deferred,ierr)
  call check(ierr==agn_deposit_invalid_receiver .and. all(row==before) .and. deferred==0d0, &
       'negative incoming internal energy rejects; do not silently lose returned mass')
  call check(agn_eddington_ratio(0d0,0d0)==0d0,'idle Eddington ratio is finite zero')
  write(*,'(A)') 'AGN_NATIVE_CELL_COUPLING_SMOKE_OK'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       write(*,'(A)') 'FAIL: '//label
       error stop 1
    endif
    write(*,'(A)') 'PASS: '//label
  end subroutine check
end program agn_feedback_deposition_smoke

program test_adm_hpm_closure
  use adm_hpm_mod, only: adm_hpm_pressure, adm_hpm_acceleration, adm_hpm_sound_speed
  implicit none

  real(8)::pressure,accel,cs

  pressure=adm_hpm_pressure(2.0d0,3.0d0,5.0d0/3.0d0)
  if(abs(pressure-4.0d0)>1.0d-14) then
     write(*,*) 'FAIL: HPM ideal-gas pressure',pressure
     stop 1
  end if

  accel=adm_hpm_acceleration(1.0d0,5.0d0,2.0d0,2.0d0)
  if(abs(accel+0.5d0)>1.0d-14) then
     write(*,*) 'FAIL: HPM pressure-gradient sign',accel
     stop 2
  end if
  accel=adm_hpm_acceleration(4.0d0,4.0d0,2.0d0,2.0d0)
  if(abs(accel)>1.0d-14) then
     write(*,*) 'FAIL: HPM uniform-pressure force',accel
     stop 3
  end if

  cs=adm_hpm_sound_speed(3.0d0,5.0d0/3.0d0)
  if(abs(cs-sqrt(10.0d0/3.0d0))>1.0d-14) then
     write(*,*) 'FAIL: HPM sound speed',cs
     stop 4
  end if

  write(*,*) 'ADM HPM closure unit test passed'
end program test_adm_hpm_closure

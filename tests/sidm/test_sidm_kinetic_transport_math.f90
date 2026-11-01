program test_sidm_kinetic_transport_math
  use sidm_kinetic_transport_math
  implicit none
  integer, parameter :: dp=kind(1.0d0),i8=selected_int_kind(18)
  real(dp), parameter :: pi=acos(-1.0d0),box_length=8.0d0*pi
  real(dp), parameter :: dx=box_length/32.0d0,k=0.25d0
  real(dp) :: before1(3),before2(3),after1(3),after2(3)
  real(dp) :: dk1,dk2,p_before(3),p_after(3),e_before,e_after
  real(dp) :: source(2),source_abs,source_k0(2),source_abs_k0
  real(dp) :: flux(2),flux_reverse(2),flux_abs,expected_amplitude
  real(dp) :: raw(12),velocity(3),centered(6),heat,advective,pressure,total
  real(dp) :: volume
  integer :: expected_faces(8),expected_signs(8)

  expected_faces=(/2,0,0,0,0,0,0,0/)
  expected_signs=(/1,0,0,0,0,0,0,0/)
  call check_faces(0.9d0*dx,1.1d0*dx,expected_faces,expected_signs,1, &
       'positive interior crossing')
  expected_signs(1)=-1
  call check_faces(1.1d0*dx,0.9d0*dx,expected_faces,expected_signs,1, &
       'negative interior crossing')
  expected_signs(1)=1
  call check_faces(0.9d0*dx,dx,expected_faces,expected_signs,1, &
       'positive endpoint crossing')
  expected_signs(1)=-1
  call check_faces(dx,0.9d0*dx,expected_faces,expected_signs,1, &
       'negative boundary departure')
  expected_faces=0
  expected_signs=0
  call check_faces(1.1d0*dx,dx,expected_faces,expected_signs,0, &
       'half-open endpoint ownership')
  expected_faces=(/2,3,4,0,0,0,0,0/)
  expected_signs=(/1,1,1,0,0,0,0,0/)
  call check_faces(0.5d0*dx,3.5d0*dx,expected_faces,expected_signs,3, &
       'multi-slab path')
  expected_faces=0
  expected_signs=0
  expected_faces(1)=1
  expected_signs(1)=1
  call check_faces(box_length-0.5d0*dx,box_length+0.5d0*dx, &
       & expected_faces,expected_signs,1,'positive periodic crossing')
  expected_signs(1)=-1
  call check_faces(box_length+0.5d0*dx,box_length-0.5d0*dx, &
       & expected_faces,expected_signs,1,'negative periodic crossing')
  expected_faces=0
  expected_signs=0
  call check_faces(0.3d0*dx,0.3d0*dx,expected_faces,expected_signs,0, &
       'stationary path')

  before1=(/1.0d0,0.0d0,0.0d0/)
  before2=(/0.0d0,0.0d0,0.0d0/)
  after1=(/0.5d0,0.5d0,0.0d0/)
  after2=(/0.5d0,-0.5d0,0.0d0/)
  dk1=0.5d0*(sum(after1**2)-sum(before1**2))
  dk2=0.5d0*(sum(after2**2)-sum(before2**2))
  p_before=before1+before2
  p_after=after1+after2
  e_before=0.5d0*(sum(before1**2)+sum(before2**2))
  e_after=0.5d0*(sum(after1**2)+sum(after2**2))
  call close_value(dk1,-0.25d0,'event particle 1 energy change')
  call close_value(dk2,0.25d0,'event particle 2 energy change')
  call close_value(e_after,e_before,'elastic event energy conservation')
  call close_array(p_after,p_before,'elastic event momentum conservation')
  call sidm_kinetic_event_source(0.0d0,box_length/4.0d0,dk1,dk2, &
       & box_length**3,k,source,source_abs)
  call close_value(source(1),-0.25d0/box_length**3,'finite-k source real')
  call close_value(source(2),-0.25d0/box_length**3,'finite-k source imaginary')
  call close_value(source_abs,0.5d0/box_length**3,'event source absolute scale')
  call sidm_kinetic_event_source(0.0d0,box_length/4.0d0,dk1,dk2, &
       & box_length**3,0.0d0,source_k0,source_abs_k0)
  call close_value(source_k0(1),0.0d0,'elastic zero-mode source')
  call close_value(source_k0(2),0.0d0,'elastic zero-mode imaginary source')

  call sidm_kinetic_mode_flux(0.5d0*dx,2.5d0*dx,1.0d0,0.5d0, &
       & 0.5d0*pi,box_length,k,flux,flux_abs)
  expected_amplitude=pi/(4.0d0*box_length**3)*sin(pi/16.0d0)/(pi/16.0d0)
  call close_value(flux(1),expected_amplitude*cos(3.0d0*pi/32.0d0), &
       'finite-path Fourier flux real')
  call close_value(flux(2),-expected_amplitude*sin(3.0d0*pi/32.0d0), &
       'finite-path Fourier flux imaginary')
  call sidm_kinetic_mode_flux(2.5d0*dx,0.5d0*dx,-1.0d0,0.5d0, &
       & 0.5d0*pi,box_length,k,flux_reverse,flux_abs)
  call close_array(flux_reverse,-flux,'reverse path Fourier flux')
  call sidm_kinetic_mode_flux(0.5d0*dx,0.5d0*dx,0.0d0,0.5d0, &
       & 0.5d0*pi,box_length,k,flux_reverse,flux_abs)
  call close_array(flux_reverse,(/0.0d0,0.0d0/),'stationary Fourier flux')
  call close_value(flux_abs,0.0d0,'stationary absolute Fourier flux')

  volume=2.5d0
  raw=0.0d0
  raw(1)=3.0d0
  raw(2)=6.0d0
  raw(5)=18.0d0
  raw(11)=66.0d0
  raw(12)=9.0d0
  call sidm_kinetic_flux_decomposition(raw,volume,velocity,centered, &
       & heat,advective,pressure,total)
  call close_value(velocity(1),2.0d0,'longitudinal bulk velocity')
  call close_value(centered(1),6.0d0,'longitudinal centered pressure')
  call close_value(heat,3.0d0/volume,'longitudinal heat flux')
  call close_value(advective,18.0d0/volume,'longitudinal advective flux')
  call close_value(pressure,12.0d0/volume,'longitudinal pressure work')
  call close_value(total,33.0d0/volume,'longitudinal total flux')
  call close_value(heat+advective+pressure,total, &
       'longitudinal decomposition closure')

  raw=0.0d0
  raw(1)=3.0d0
  raw(2:4)=6.0d0
  raw(5:7)=18.0d0
  raw(8:10)=18.0d0
  raw(11)=198.0d0
  raw(12)=27.0d0
  call sidm_kinetic_flux_decomposition(raw,volume,velocity,centered, &
       & heat,advective,pressure,total)
  call close_array(velocity,(/2.0d0,2.0d0,2.0d0/),'three-axis bulk velocity')
  call close_array(centered,(/6.0d0,6.0d0,6.0d0,6.0d0,6.0d0,6.0d0/), &
       'three-axis centered pressure tensor')
  call close_value(heat,9.0d0/volume,'three-axis heat flux')
  call close_value(advective,54.0d0/volume,'three-axis advective flux')
  call close_value(pressure,36.0d0/volume,'transverse pressure work')
  call close_value(total,99.0d0/volume,'three-axis total flux')
  call close_value(heat+advective+pressure,total, &
       'three-axis decomposition closure')

  if(.not.sidm_kinetic_bootstrap_allowed(0.0d0,0.0d0,.true., &
       & .false.,.false.,0_i8)) error stop 1
  if(sidm_kinetic_bootstrap_allowed(0.0d0,0.0d0,.true., &
       & .false.,.true.,0_i8)) error stop 2
  if(sidm_kinetic_bootstrap_allowed(0.0d0,1.0d-6,.true., &
       & .false.,.false.,0_i8)) error stop 3
  if(sidm_kinetic_bootstrap_allowed(0.1d0,0.0d0,.true., &
       & .false.,.false.,0_i8)) error stop 4
  if(.not.sidm_kinetic_interval_matches(0.005d0,0.005d0, &
       & 0.005d0,0.005d0,1.0d-12)) error stop 5
  if(.not.sidm_kinetic_interval_matches(0.003d0,0.003d0, &
       & 0.003d0,0.003d0,1.0d-12)) error stop 6
  if(sidm_kinetic_interval_matches(0.005d0,0.005d0, &
       & 0.003d0,0.003d0,1.0d-12)) error stop 7

  write(*,'(A)') 'PASS: 23 production-math fixtures (geometry, Fourier, flux, lifecycle)'

contains

  subroutine check_faces(x0,x1,expected_face,expected_sign,nexpected,name)
    real(dp), intent(in) :: x0,x1
    integer, intent(in) :: expected_face(8),expected_sign(8),nexpected
    character(len=*), intent(in) :: name
    integer(i8) :: first,last,index
    integer :: direction,n,faces(8),signs(8),slot

    faces=0
    signs=0
    n=0
    call sidm_kinetic_crossing_range(x0,x1,dx,first,last,direction)
    if(direction>0 .and. last>=first) then
       do index=first,last
          n=n+1
          slot=1+int(modulo(index,32_i8))
          faces(n)=slot
          signs(n)=1
       end do
    else if(direction<0 .and. first>=last) then
       do index=first,last,-1_i8
          n=n+1
          slot=1+int(modulo(index,32_i8))
          faces(n)=slot
          signs(n)=-1
       end do
    end if
    if(n/=nexpected) then
       write(*,'(A,1X,A,2(1X,I0))') 'FAIL crossing count',trim(name),n,nexpected
       error stop 10
    end if
    if(any(faces/=expected_face) .or. any(signs/=expected_sign)) then
       write(*,'(A,1X,A)') 'FAIL crossing face/sign',trim(name)
       error stop 11
    end if
  end subroutine check_faces

  subroutine close_value(value,expected,name)
    real(dp), intent(in) :: value,expected
    character(len=*), intent(in) :: name
    if(abs(value-expected)>2.0d-13*max(1.0d0,abs(expected))) then
       write(*,'(A,1X,A,2(1X,ES24.16))') &
            'FAIL',trim(name),value,expected
       error stop 12
    end if
  end subroutine close_value

  subroutine close_array(value,expected,name)
    real(dp), intent(in) :: value(:),expected(:)
    character(len=*), intent(in) :: name
    if(size(value)/=size(expected)) error stop 13
    if(any(abs(value-expected)>2.0d-13*max(1.0d0,maxval(abs(expected))))) then
       write(*,'(A,1X,A)') 'FAIL',trim(name)
       error stop 14
    end if
  end subroutine close_array

end program test_sidm_kinetic_transport_math

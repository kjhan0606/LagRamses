! Exercise the real M_N dispatcher with a deliberately unavailable stream.
! No CUDA hardware or RAMSES state is needed for this prelaunch policy test.
program mn_dispatch_policy_test
  use snrt_runtime_backend
  use snrt_moment_transport
  use snrt_moment_dispatch
  implicit none
  type(mn_basis) :: b
  real(mn_dp),allocatable :: angular(:,:),moments(:,:),input(:,:),projected(:,:),cache(:,:),warm(:,:)
  real(mn_dp),allocatable :: harmonic(:,:),flux(:,:)
  real(mn_dp) :: weight(1),direction(3,1),left(4,1,1),right(4,1,1),geometry(5,1)
  real(mn_dp) :: tensor(4,4)
  integer :: status,order

  call mn_initialize(b,5,status)
  if(status/=mn_ok.or.b%nq/=384.or.b%nm/=36)stop 1
  allocate(angular(b%nq,2),moments(b%nm,2),input(b%nm,2),projected(b%nm,2), &
       cache(b%nm+1,2),warm(b%nm-1,2),harmonic(b%nm,1),flux(b%nm,1))
  angular=1;input=0;warm=0;harmonic=0;harmonic(1,1)=1
  weight=1;direction=0;direction(1,1)=1;left=0;right=0
  left(1,1,1)=1;geometry(:,1)=[1d0,1d0,0d0,0d0,1d0]

  call mn_dispatch_initialize(status)
  if(status/=0)stop 2
  call mn_dispatch_closure(b,input,projected,cache,warm,status)
  if(status/=7.or.mn_cpu_closures/=0.or.mn_gpu_closures/=0)stop 3
  call mn_dispatch_project(b,angular,moments,status)
  if(status/=7)stop 4
  call mn_dispatch_flux(harmonic,weight,direction,left,right,geometry,flux,status)
  if(status/=7.or.mn_cpu_faces/=0.or.mn_gpu_faces/=0)stop 5

  test_choice=0
  call mn_dispatch_initialize(status)
  if(status/=0)stop 6
  warm=3
  call mn_dispatch_closure(b,input,projected,cache,warm,status,angular_out=angular)
  if(status/=0.or.mn_cpu_closures/=2.or.any(projected/=0).or.any(cache/=0).or. &
       any(warm/=3).or.any(angular/=0))stop 7
  input(1,2)=1
  warm(:,2)=0
  call mn_dispatch_closure(b,input,projected,cache,warm,status,angular_out=angular)
  if(status/=0.or.mn_cpu_closures/=4.or.any(projected(:,1)/=0).or.any(cache(:,1)/=0).or. &
       any(warm(:,1)/=3).or.any(angular(:,1)/=0).or.abs(projected(1,2)-1d0)>1d-12.or. &
       any(angular(:,2)<=0))stop 11
  input(1,2)=-1
  call mn_dispatch_closure(b,input,projected,cache,warm,status)
  if(status/=mn_bad_input)stop 12
  angular=1
  call mn_dispatch_project(b,angular,moments,status)
  if(status/=0.or.abs(moments(1,1)-1d0)>1d-12)stop 8
  call mn_dispatch_flux(harmonic,weight,direction,left,right,geometry,flux,status)
  if(status/=0.or.mn_cpu_faces/=1.or.abs(flux(1,1)-1d0)>1d-12)stop 9
  deallocate(angular,moments,input,projected,cache,warm,harmonic,flux)
  test_choice=1;test_threads=4
  do order=1,5
    call mn_initialize(b,order,status)
    if(status/=0)stop 13
    allocate(angular(b%nq,64),moments(b%nm,64),input(b%nm,64),projected(b%nm,64), &
         cache(b%nm+1,64),warm(b%nm-1,64))
    input=0;input(1,:)=1;warm=0
    call mn_dispatch_initialize(status)
    if(status/=0)stop 14
    call mn_dispatch_closure(b,input,projected,cache,warm,status,angular_out=angular)
    if(status/=0.or.mn_cpu_closures/=64.or.mn_gpu_closures/=0.or. &
         maxval(abs(projected-input))>1d-10.or.any(angular<=0))stop 15
    call mn_dispatch_project(b,angular,moments,status)
    if(status/=0.or.maxval(abs(moments-input))>1d-10)stop 16
    call mn_stress_tensor(b,input(:,1),tensor,status)
    if(status/=0.or.abs(tensor(1,1)-1d0)>1d-12.or. &
         abs(sum([tensor(2,2),tensor(3,3),tensor(4,4)])-1d0)>1d-12)stop 17
    deallocate(angular,moments,input,projected,cache,warm)
  enddo
  print '(A)','MN_DISPATCH_POLICY_PASS'
end program

function stream_unavailable(bytes,sharers) bind(C,name='snrt_hybrid_try_acquire_c') result(slot)
  use iso_c_binding, only: c_int,c_long_long
  integer(c_long_long),value :: bytes
  integer(c_int),value :: sharers
  integer(c_int) :: slot
  slot=-1
end function

subroutine unused_release(slot) bind(C,name='cuda_release_stream')
  use iso_c_binding, only: c_int
  integer(c_int),value :: slot
  stop 10
end subroutine

function unused_closure(slot,nc,nm,nq,y,w,dir,input,warm,projected,cache) &
    bind(C,name='snrt_mn_cuda_closure_hint_c') result(status)
  use iso_c_binding, only: c_int,c_double
  integer(c_int),value :: slot,nc,nm,nq
  real(c_double) :: y(*),w(*),dir(*),input(*),warm(*),projected(*),cache(*)
  integer(c_int) :: status
  status=99
end function

function unused_project(slot,nc,nm,nq,y,w,angular,moments) &
    bind(C,name='snrt_mn_cuda_project_c') result(status)
  use iso_c_binding, only: c_int,c_double
  integer(c_int),value :: slot,nc,nm,nq
  real(c_double) :: y(*),w(*),angular(*),moments(*)
  integer(c_int) :: status
  status=99
end function

function unused_flux(slot,nf,nm,nq,y,w,dir,left,right,geometry,flux) &
    bind(C,name='snrt_mn_cuda_flux_c') result(status)
  use iso_c_binding, only: c_int,c_double
  integer(c_int),value :: slot,nf,nm,nq
  real(c_double) :: y(*),w(*),dir(*),left(*),right(*),geometry(*),flux(*)
  integer(c_int) :: status
  status=99
end function

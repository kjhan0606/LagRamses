program hlld_device_smoke
  use iso_c_binding
  use ieee_arithmetic
  use omp_lib
  use hydro_parameters
  use mhd_dispatch
  implicit none
  integer,parameter::nf=4096
  real(dp)::ql(nvar,nf),qr(nvar,nf),got(nvar+1,nf),ref(nvar+1,nf),a(nvar),b(nvar)
  real(dp)::x,err,t0,tc,tg
  integer(c_int)::mask(nf),slot
  integer::i,j,rep,caseid
  logical::ok
  interface
     subroutine pool_init(rank,ns) bind(C,name='cuda_pool_init')
       import c_int
       integer(c_int),value::rank,ns
     end subroutine
     integer(c_int) function acquire() bind(C,name='cuda_acquire_stream')
       import c_int
     end function
     subroutine release(slot) bind(C,name='cuda_release_stream')
       import c_int
       integer(c_int),value::slot
     end subroutine
     subroutine finalize() bind(C,name='cuda_pool_finalize')
     end subroutine
  end interface
  call pool_init(0_c_int,1_c_int)
  mhd_gpu_faces=.true.;ischeme=0;iriemann=3
  do caseid=1,2
     gamma=5d0/3d0
     if(caseid==2)gamma=2d0
     do i=1,nf
        x=dble(i)*.037d0
        ql(:,i)=.01d0;qr(:,i)=.02d0
        ql(1,i)=10d0**(2d0*sin(x));qr(1,i)=10d0**(2d0*cos(x))
        ql(2,i)=10d0**(2d0*cos(x*.7d0));qr(2,i)=10d0**(2d0*sin(x*.7d0))
        do j=3,8
           ql(j,i)=.2d0*sin(x*j);qr(j,i)=.2d0*cos(x*j)
        enddo
        if(mod(i,5)==0)then
           ql(4:8:2,i)=0d0;qr(4:8:2,i)=0d0
        endif
        mask(i)=1
        if(mod(i,7)==0)mask(i)=0
        a=ql(:,i);b=qr(:,i)
        call hlld(a,b,ref(:,i))
     enddo
     ok=mhd_try_faces(ql,qr,mask,got,nf)
     if(.not.ok)error stop 'No real GPU work performed'
     err=0d0
     do i=1,nf
        if(mask(i)==0)cycle
        if(any(.not.ieee_is_finite(got(:,i))))error stop 'Nonfinite device flux'
        err=max(err,maxval(abs(got(:,i)-ref(:,i))/(1d0+abs(ref(:,i)))))
     enddo
     write(*,*)'HLLD_DEVICE_PARITY gamma,error=',gamma,err
     if(err>2d-10)error stop 'CPU/device flux mismatch'
  enddo
  slot=acquire()
  if(slot<0)error stop 'Could not reserve stream for busy test'
  ok=mhd_try_faces(ql,qr,mask,got,nf)
  if(ok)error stop 'Busy stream was stolen'
  call release(slot)
  t0=omp_get_wtime()
  do rep=1,20
     do i=1,nf
        if(mask(i)==0)cycle
        a=ql(:,i);b=qr(:,i)
        call hlld(a,b,ref(:,i))
     enddo
  enddo
  tc=omp_get_wtime()-t0
  t0=omp_get_wtime()
  do rep=1,20
     ok=mhd_try_faces(ql,qr,mask,got,nf)
     if(.not.ok)error stop 'Unexpected unavailable stream'
  enddo
  tg=omp_get_wtime()-t0
  write(*,*)'HLLD_BATCH_TIMING CPU,CUDA_seconds=',tc,tg
  call mhd_dispatch_report(0,1)
  call finalize()
  write(*,*)'HLLD_DEVICE_SMOKE_OK'
end program

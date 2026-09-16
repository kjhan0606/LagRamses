! Same compact M_N operator on OMP workers or an available shared CUDA
! stream. GPU failures are transaction failures, never post-launch replay.
module snrt_moment_dispatch
  use snrt_moment_transport
  use snrt_runtime_backend, only: snrt_runtime_mn_policy
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,save::choice=1,sharing=1,threads=1
  integer,public,save::mn_cpu_closures=0,mn_gpu_closures=0,mn_cpu_faces=0,mn_gpu_faces=0
  public::mn_dispatch_initialize,mn_dispatch_closure,mn_dispatch_project,mn_dispatch_flux
#ifdef HYDRO_CUDA
  interface
     function acquire(bytes,sharers) bind(C,name='snrt_hybrid_try_acquire_c') result(slot)
       import c_long_long,c_int
       integer(c_long_long),value::bytes
       integer(c_int),value::sharers
       integer(c_int)::slot
     end function
     subroutine release(slot) bind(C,name='cuda_release_stream')
       import c_int
       integer(c_int),value::slot
     end subroutine
     function gpu_closure(slot,nc,nm,nq,y,w,dir,input,projected,cache) bind(C,name='snrt_mn_cuda_closure_c') result(ierr)
       import c_int,c_double
       integer(c_int),value::slot,nc,nm,nq
       real(c_double),intent(in)::y(*),w(*),dir(*),input(*)
       real(c_double),intent(out)::projected(*),cache(*)
       integer(c_int)::ierr
     end function
     function gpu_project(slot,nc,nm,nq,y,w,angular,moments) bind(C,name='snrt_mn_cuda_project_c') result(ierr)
       import c_int,c_double
       integer(c_int),value::slot,nc,nm,nq
       real(c_double),intent(in)::y(*),w(*),angular(*)
       real(c_double),intent(out)::moments(*)
       integer(c_int)::ierr
     end function
     function gpu_flux(slot,nf,nm,nq,y,w,dir,left,right,geometry,flux) bind(C,name='snrt_mn_cuda_flux_c') result(ierr)
       import c_int,c_double
       integer(c_int),value::slot,nf,nm,nq
       real(c_double),intent(in)::y(*),w(*),dir(*),left(*),right(*),geometry(*)
       real(c_double),intent(out)::flux(*)
       integer(c_int)::ierr
     end function
  end interface
#endif
contains
  subroutine mn_dispatch_initialize(ierr)
    integer,intent(out)::ierr
    call snrt_runtime_mn_policy(choice,sharing,threads,ierr)
    mn_cpu_closures=0;mn_gpu_closures=0;mn_cpu_faces=0;mn_gpu_faces=0
  end subroutine

  subroutine mn_dispatch_closure(b,input,projected,cache,warm,ierr,angular_out)
    type(mn_basis),intent(in)::b
    real(mn_dp),intent(in)::input(:,:)
    real(mn_dp),intent(out)::projected(:,:),cache(:,:)
    real(mn_dp),intent(inout)::warm(:,:)
    integer,intent(out)::ierr
    real(mn_dp),optional,intent(out)::angular_out(:,:)
    integer::first,last,i,q,status,slot,cpu,gpu,local_error,nc
    logical::gpu_supported
    real(mn_dp)::angular(b%nq)
    real(mn_dp)::eta,target(b%nm-1)
    nc=size(input,2);ierr=0;cpu=0;gpu=0
    ! The live CUDA closure currently targets the fixed M5 quadrature.  An
    ! unsupported basis is a normal auto-dispatch choice, not a device error;
    ! force-CUDA must still reject it explicitly.
    gpu_supported=b%nq==384.and.b%nm>=9.and.b%nm<=36
    if(present(angular_out))then
       if(size(angular_out,1)/=b%nq.or.size(angular_out,2)/=nc)then
          ierr=7;return
       endif
    endif
!$omp parallel do schedule(dynamic,1) num_threads(threads) &
!$omp private(first,last,i,q,status,slot,angular,eta,target,local_error) reduction(max:ierr) reduction(+:cpu,gpu)
    do first=1,nc,256
       last=min(nc,first+255);slot=-1;local_error=0
#ifdef HYDRO_CUDA
       if(choice/=1.and.gpu_supported)slot=acquire(33554432_c_long_long+8_c_long_long*(5*b%nq+5*b%nm)*(last-first+1),sharing)
       if(slot>=0)then
          status=gpu_closure(slot,last-first+1,b%nm,b%nq,b%harmonic,b%weight,b%direction,input(:,first:last), &
               projected(:,first:last),cache(:,first:last))
          call release(slot)
          gpu=gpu+last-first+1;local_error=status
          if(status==0)then
             warm(:,first:last)=cache(1:b%nm-1,first:last)
             if(present(angular_out))then
                ! The CUDA closure exports the same dual coefficients, logit
                ! maximum, and partition used by its angular reconstruction.
                ! Replay that closed form here instead of running a second
                ! CPU Newton solve solely to feed the material stage.
                do i=first,last
                   angular_out(:,i)=0d0
                   if(input(1,i)>0d0.and.cache(b%nm+1,i)>0d0)then
                      target=input(2:,i)/input(1,i)
                      do q=1,b%nq
                         eta=dot_product(cache(1:b%nm-1,i),b%harmonic(2:,q)-target)
                         angular_out(q,i)=input(1,i)*exp(eta-cache(b%nm,i))/cache(b%nm+1,i)
                      enddo
                   endif
                enddo
             endif
          endif
       else
#endif
          if(choice==2.and..not.gpu_supported)then
             local_error=7
          else
             do i=first,last
                call mn_reconstruct(b,input(:,i),angular,status,dual_hint=warm(:,i),stream_cache=cache(:,i))
                if(status/=0)call mn_reconstruct(b,input(:,i),angular,status,stream_cache=cache(:,i))
                if(status==0)call mn_project(b,angular,projected(:,i),status)
                if(status==0.and.present(angular_out))angular_out(:,i)=angular
                local_error=max(local_error,status)
             enddo
             cpu=cpu+last-first+1
          endif
#ifdef HYDRO_CUDA
       endif
#endif
       ierr=max(ierr,local_error)
    enddo
!$omp end parallel do
    mn_cpu_closures=mn_cpu_closures+cpu;mn_gpu_closures=mn_gpu_closures+gpu
  end subroutine

  subroutine mn_dispatch_project(b,angular,moments,ierr)
    type(mn_basis),intent(in)::b
    real(mn_dp),intent(in)::angular(:,:)
    real(mn_dp),intent(out)::moments(:,:)
    integer,intent(out)::ierr
    integer::i,status,slot,nc
    logical::gpu_supported
    nc=size(angular,2);ierr=7;slot=-1
    if(size(angular,1)/=b%nq.or.size(moments,1)/=b%nm.or.size(moments,2)/=nc)return
    if(any(.not.ieee_is_finite(angular)).or.any(angular<0d0))then
       ierr=1;return
    endif
#ifdef HYDRO_CUDA
    gpu_supported=b%nq==384.and.b%nm>=9.and.b%nm<=36
    if(choice/=1.and.gpu_supported)slot=acquire(16777216_c_long_long+8_c_long_long*(b%nq+b%nm)*nc,sharing)
    if(slot>=0)then
       status=gpu_project(slot,nc,b%nm,b%nq,b%harmonic,b%weight,angular,moments)
       call release(slot)
       ierr=status
       if(ierr==0.and.any(.not.ieee_is_finite(moments)))ierr=3
       return
    endif
    if(choice==2.and..not.gpu_supported)return
#endif
    ierr=0
!$omp parallel do num_threads(threads) private(i,status) reduction(max:ierr)
    do i=1,nc
       call mn_project(b,angular(:,i),moments(:,i),status)
       ierr=max(ierr,status)
    enddo
!$omp end parallel do
  end subroutine

  subroutine mn_dispatch_flux(y,w,dir,left,right,geometry,flux,ierr)
    real(mn_dp),intent(in)::y(:,:),w(:),dir(:,:),left(:,:,:),right(:,:,:),geometry(:,:)
    real(mn_dp),intent(out)::flux(:,:)
    integer,intent(out)::ierr
    integer::slot,f,j,q,a,nf,nm,nq
    logical::gpu_supported
    real(mn_dp)::sign,mu,il,ir,flow
    nf=size(left,3);nm=size(y,1);nq=size(w);slot=-1;ierr=0;gpu_supported=.false.
    if(nf==0)return
#ifdef HYDRO_CUDA
    gpu_supported=nq<=8
    if(choice/=1.and.gpu_supported)slot=acquire(16777216_c_long_long+8_c_long_long*nf*(8*nq+nm+5),sharing)
    if(slot>=0)then
       ierr=gpu_flux(slot,nf,nm,nq,y,w,dir,left,right,geometry,flux)
       call release(slot)
       mn_gpu_faces=mn_gpu_faces+nf
       return
    endif
#endif
    if(choice==2.and..not.gpu_supported)then
       ierr=7;return
    endif
!$omp parallel do num_threads(threads) private(f,j,q,a,sign,mu,il,ir,flow)
    do f=1,nf
       a=nint(geometry(1,f));sign=geometry(2,f);flux(:,f)=0
       do q=1,nq
          mu=sign*dir(a,q)
          il=left(1,q,f)+sign*.5d0*geometry(3,f)*left(a+1,q,f)
          ir=right(1,q,f)-sign*.5d0*geometry(4,f)*right(a+1,q,f)
          flow=geometry(5,f)*(max(mu,0d0)*il+min(mu,0d0)*ir)
          flux(:,f)=flux(:,f)+y(:,q)*(w(q)*flow)
       enddo
    enddo
!$omp end parallel do
    mn_cpu_faces=mn_cpu_faces+nf
    if(any(.not.ieee_is_finite(flux)))ierr=3
  end subroutine
end module

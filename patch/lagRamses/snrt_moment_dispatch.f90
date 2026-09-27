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
     function gpu_closure(slot,nc,nm,nq,y,w,dir,input,warm,projected,cache) &
          bind(C,name='snrt_mn_cuda_closure_hint_c') result(ierr)
       import c_int,c_double
       integer(c_int),value::slot,nc,nm,nq
       real(c_double),intent(in)::y(*),w(*),dir(*),input(*),warm(*)
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
    integer::first,last,i,q,status,slot,cpu,gpu,local_error,nc,next_column,batch_size
    logical::gpu_supported
    real(mn_dp)::angular(b%nq)
    real(mn_dp)::eta,target(b%nm-1)
    nc=size(input,2);ierr=0;cpu=0;gpu=0
    ! The live CUDA closure supports the default 384-node M1--M5 bases.
    ! A nondefault quadrature is an auto-dispatch CPU choice, while
    ! force-CUDA must still reject it explicitly.
    gpu_supported=b%nq==384.and.b%nm>=4.and.b%nm<=36
    if(present(angular_out))then
       if(size(angular_out,1)/=b%nq.or.size(angular_out,2)/=nc)then
          ierr=7;return
       endif
    endif
    if(choice/=2.and.all(input==0d0))then
       projected=0d0;cache=0d0
       if(present(angular_out))angular_out=0d0
       mn_cpu_closures=mn_cpu_closures+nc
       return
    endif
    next_column=1
!$omp parallel num_threads(threads) &
!$omp private(first,last,i,q,status,slot,angular,eta,target,local_error,batch_size) reduction(max:ierr) reduction(+:cpu,gpu)
    do
       slot=-1;local_error=0;batch_size=256
#ifdef HYDRO_CUDA
       ! The smallest live closure call that beats four CPU workers on A100
       ! for every M1--M5 order is 64 cells. Keep tiny calls on CPU in auto
       ! mode; forced CUDA still probes the device at any valid batch size.
       if(choice/=1.and.gpu_supported.and.(choice==2.or.nc>=64))slot=acquire(33554432_c_long_long+ &
            8_c_long_long*(5*b%nq+5*b%nm)*min(nc,256),sharing)
#endif
       ! All CPU arms, including forced OpenMP and nondefault quadratures,
       ! need work-sharing below 256 cells. GPU leases keep the 256-cell tile.
       if(slot<0.and.choice/=2)batch_size=16
!$omp atomic capture
       first=next_column
       next_column=next_column+batch_size
!$omp end atomic
       if(first>nc)then
#ifdef HYDRO_CUDA
          if(slot>=0)call release(slot)
#endif
          exit
       endif
       last=min(nc,first+batch_size-1)
#ifdef HYDRO_CUDA
       if(choice==0.and.slot>=0.and.last-first+1<64)then
          call release(slot)
          slot=-1
       endif
#endif
       ! Vacuum has an exact closure. Do not allocate/copy/launch a CUDA
       ! batch just to return zeros; forced CUDA still exercises the device.
       if(choice/=2.and.all(input(:,first:last)==0d0))then
#ifdef HYDRO_CUDA
          if(slot>=0)call release(slot)
#endif
          projected(:,first:last)=0d0;cache(:,first:last)=0d0
          if(present(angular_out))angular_out(:,first:last)=0d0
          cpu=cpu+last-first+1
          cycle
       endif
#ifdef HYDRO_CUDA
       if(slot>=0)then
          status=gpu_closure(slot,last-first+1,b%nm,b%nq,b%harmonic,b%weight,b%direction,input(:,first:last), &
               warm(:,first:last), &
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
          ! Forced CUDA must not silently turn into an OpenMP operator when
          ! the basis is unsupported or no stream/memory lease is available.
          if(choice==2)then
             local_error=7
          else
             ! Exact vacuum state: mn_reconstruct returns zero angular/cache
             ! and leaves the warm hint unchanged; projecting zero is zero.
             ! Avoid a 36x384 projection for every dark cell in the live
             ! source-free material sweep. Nonzero/invalid states retain the
             ! full realizability and finiteness checks below.
             if(all(input(:,first:last)==0d0))then
                projected(:,first:last)=0d0
                cache(:,first:last)=0d0
                if(present(angular_out))angular_out(:,first:last)=0d0
             else
                do i=first,last
                   if(all(input(:,i)==0d0))then
                      projected(:,i)=0d0
                      cache(:,i)=0d0
                      if(present(angular_out))angular_out(:,i)=0d0
                      cycle
                   endif
                   call mn_reconstruct(b,input(:,i),angular,status,dual_hint=warm(:,i),stream_cache=cache(:,i))
                   if(status/=0)call mn_reconstruct(b,input(:,i),angular,status,stream_cache=cache(:,i))
                   if(status==0)call mn_project(b,angular,projected(:,i),status)
                   if(status==0.and.present(angular_out))angular_out(:,i)=angular
                   local_error=max(local_error,status)
                enddo
             endif
             cpu=cpu+last-first+1
          endif
#ifdef HYDRO_CUDA
       endif
#endif
       ierr=max(ierr,local_error)
    enddo
!$omp end parallel
    mn_cpu_closures=mn_cpu_closures+cpu;mn_gpu_closures=mn_gpu_closures+gpu
  end subroutine

  subroutine mn_dispatch_project(b,angular,moments,ierr)
    type(mn_basis),intent(in)::b
    real(mn_dp),intent(in)::angular(:,:)
    real(mn_dp),intent(out)::moments(:,:)
    integer,intent(out)::ierr
    integer::i,status,slot,nc,first,last,local_error
    logical::gpu_supported
    nc=size(angular,2);ierr=7;slot=-1
    if(size(angular,1)/=b%nq.or.size(moments,1)/=b%nm.or.size(moments,2)/=nc)return
    if(any(.not.ieee_is_finite(angular)).or.any(angular<0d0))then
       ierr=1;return
    endif
    if(choice/=2.and.all(angular==0d0))then
       moments=0d0;ierr=0;return
    endif
    gpu_supported=b%nq==384.and.b%nm>=4.and.b%nm<=36
    ierr=0
    ! Material projection batches are cheap compared with the device
    ! allocation/copy/sync: A100 1--256-cell measurements favor four CPU
    ! workers for M1--M4, and give no robust M5 device advantage. The
    ! explicit CUDA arm remains available for controlled comparisons.
    if(choice/=2)then
!$omp parallel do schedule(static) num_threads(threads) private(i,status) reduction(max:ierr)
       do i=1,nc
          if(all(angular(:,i)==0d0))then
             moments(:,i)=0d0
          else
             call mn_project(b,angular(:,i),moments(:,i),status)
             ierr=max(ierr,status)
          endif
       enddo
!$omp end parallel do
       return
    endif
#ifdef HYDRO_CUDA
    if(.not.gpu_supported)then
       ierr=7;return
    endif
!$omp parallel do schedule(dynamic,1) num_threads(threads) private(first,last,i,status,slot,local_error) reduction(max:ierr)
    do first=1,nc,256
       last=min(nc,first+255);slot=-1;local_error=0
       slot=acquire(16777216_c_long_long+ &
            8_c_long_long*(b%nq+b%nm)*(last-first+1),sharing)
       if(slot>=0)then
          local_error=gpu_project(slot,last-first+1,b%nm,b%nq,b%harmonic,b%weight, &
               angular(:,first:last),moments(:,first:last))
          call release(slot)
          if(local_error==0.and.any(.not.ieee_is_finite(moments(:,first:last))))local_error=3
       else
          ! No post-launch replay. Forced CUDA rejects a missing lease.
          local_error=7
       endif
       ierr=max(ierr,local_error)
    enddo
!$omp end parallel do
#else
    ierr=7
#endif
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
    ! The live face kernel is bounded to a 32-node angular tile.  Larger
    ! quadrature sets are still split by the caller; this is not a change to
    ! the angular quadrature itself.
    gpu_supported=nq<=32
    ! A100 same-input batches of 1--256 faces are 10--40x faster on CPU:
    ! this 32-angle kernel pays a full allocation/copy/sync per call. In
    ! hybrid mode, reserve shared GPU streams for expensive closures instead.
    ! Explicit CUDA mode still runs the face kernel for parity/comparison.
    if(choice==2.and.gpu_supported)slot=acquire(16777216_c_long_long+8_c_long_long*nf*(8*nq+nm+5),sharing)
    if(slot>=0)then
       ierr=gpu_flux(slot,nf,nm,nq,y,w,dir,left,right,geometry,flux)
       call release(slot)
       mn_gpu_faces=mn_gpu_faces+nf
       return
    endif
#endif
    if(choice==2)then
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

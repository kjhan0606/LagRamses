module mhd_dispatch
  use iso_c_binding, only: c_int,c_double,c_long_long
  use hydro_parameters
  implicit none
  private
  public :: mhd_try_faces,mhd_switch_needed,mhd_dispatch_report
  integer(c_long_long),save :: gpu_faces_done=0,cpu_faces_done=0,busy_batches=0
#ifdef HYDRO_CUDA
  interface
     integer(c_int) function device_faces(ql,qr,mask,flux,n,nv,gam) bind(C,name='mhd_hlld_batch')
       import c_int,c_double
       real(c_double),intent(in) :: ql(*),qr(*)
       integer(c_int),intent(in) :: mask(*)
       real(c_double),intent(out) :: flux(*)
       integer(c_int),value :: n,nv
       real(c_double),value :: gam
     end function
  end interface
#endif
contains
  logical function mhd_switch_needed(ql,qr) result(needed)
    real(dp),intent(in)::ql(nvar),qr(nvar)
    needed=.false.
    if(.not.allow_switch_solver.or.(iriemann/=2.and.iriemann/=3))return
    needed=(qr(4)**2+qr(6)**2+qr(8)**2)/qr(2)>switch_solv_B.or. &
         (ql(4)**2+ql(6)**2+ql(8)**2)/ql(2)>switch_solv_B.or. &
         ql(1)/qr(1)>switch_solv_dens.or.qr(1)/ql(1)>switch_solv_dens.or. &
         min(ql(1),qr(1))<switch_solv_min_dens
  end function

  logical function mhd_try_faces(ql,qr,mask,flux,n) result(done)
    integer,intent(in)::n
    real(dp),intent(in)::ql(nvar,n),qr(nvar,n)
    integer(c_int),intent(in)::mask(n)
    real(dp),intent(out)::flux(nvar+1,n)
    integer::status,ngpu
    done=.false.;ngpu=0;status=0
#ifdef HYDRO_CUDA
    if(mhd_gpu_faces.and.nener==0.and.iriemann==3.and.ischeme/=1.and.any(mask/=0))then
       status=device_faces(ql,qr,mask,flux,int(n,c_int),int(nvar,c_int),gamma)
       if(status<0)then
          write(*,*)'FATAL: CUDA MHD face batch failed; no partial flux committed'
          ! No MPI call from an OpenMP worker. Fail visibly, never consume output.
          error stop 71
       endif
       done=status==1
       if(done)ngpu=count(mask/=0)
       if(status==0)then
          !$omp atomic update
          busy_batches=busy_batches+1
       endif
    endif
#endif
    !$omp atomic update
    gpu_faces_done=gpu_faces_done+ngpu
    !$omp atomic update
    cpu_faces_done=cpu_faces_done+n-ngpu
  end function

  subroutine mhd_dispatch_report(level,rank)
    integer,intent(in)::level,rank
    ! Called outside the parallel region; local cumulative work, no reduction.
    if(mhd_omp.or.mhd_gpu_faces)write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') &
         'MHD_DISPATCH rank=',rank,' level=',level,' gpu_faces=',gpu_faces_done, &
         ' cpu_faces=',cpu_faces_done,' busy_or_unavailable_batches=',busy_batches
  end subroutine
end module

program dust_backend_smoke
  use snrt_dust_ir
  use snrt_runtime_backend
  use mpi_mod
  implicit none
  integer,parameter::nc=17,ng=2,nd=2
  type(dust_ir_table)::table
  type(dust_ir_diagnostics)::ref_diag,trial_diag
  real(dust_dp)::rays(3,nd),weights(nd),rho(nc),heat(nc),cap(nc)
  real(dust_dp)::ref(ng,nd,nc),trial(ng,nd,nc),ref_t(nc),trial_t(nc),ref_p(ng,nc),trial_p(ng,nc)
  real(dust_dp)::ref_e(nc),trial_e(nc),error
  integer::links(6,nc),i,mode,ierr,info
  call MPI_INIT(info)
  call snrt_backend_initialize(ierr)
  if(ierr/=0)then
     write(*,*)'BACKEND_INIT_REJECTED',ierr
     call MPI_ABORT(MPI_COMM_WORLD,2,info)
  endif
  rays(:,1)=[1d0,0d0,0d0];rays(:,2)=[-1d0,0d0,0d0];weights=.5d0;links=0
  do i=1,nc
     rho(i)=.5d0+real(i,dust_dp)/nc
  enddo
  rho(1)=0;cap=1d-24
  do mode=0,1
     if(mode==0)then
        call snrt_dust_ir_initialize(table,[.001d0,.01d0],[.001d0,.01d0],[1d-21,1d-21], &
             [10d0,20d0,50d0,100d0],10d0,ierr)
        ref_e=20d0*cap
     else
        call snrt_dust_ir_initialize(table,[.001d0,.01d0],[.001d0,.01d0],[1d-21,1d-21], &
             [10d0,20d0,50d0,100d0],10d0,ierr,[1d-23,2d-23,8d-23,2d-22])
        ref_e=rho*2d-23
     endif
     if(ierr/=0)stop 3
     heat=1d-30*rho;ref=0;ref_t=20;ref_p=0
     trial=ref;trial_t=ref_t;trial_p=ref_p;trial_e=ref_e
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          ref,ref_t,ref_p,ref_diag,ierr,1d-10,128,ref_e,cap)
     if(ierr/=0)stop 4
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap, &
          material_dispatch=snrt_runtime_dust_material)
     if(ierr/=0)stop 5
     error=max(maxval(abs(ref-trial))/maxval(ref),maxval(abs(ref_t-trial_t))/maxval(ref_t), &
          maxval(abs(ref_p-trial_p))/maxval(ref_p),maxval(abs(ref_e-trial_e))/maxval(ref_e))
     if(error>1d-10)stop 6
     write(*,'(A,I0,A,ES14.6)')'DUST_BACKEND_FORTRAN_PARITY material_u=',mode,' relative=',error
     ! Error after entering the material callback, not just outer shape checks.
     heat=1d99
     ref=trial;ref_e=trial_e;ref_t=trial_t;ref_p=trial_p
     call snrt_dust_ir_advance(table,rays,weights,links,1d12,1d6,1d5,rho,heat, &
          trial,trial_t,trial_p,trial_diag,ierr,1d-10,128,trial_e,cap, &
          material_dispatch=snrt_runtime_dust_material)
     if(ierr==0.or.any(trial/=ref).or.any(trial_e/=ref_e).or.any(trial_p/=ref_p).or.any(trial_t/=ref_t))stop 7
  enddo
  write(*,*)'DUST_BACKEND_PARITY_AND_ROLLBACK_PASS'
  call MPI_FINALIZE(info)
end program

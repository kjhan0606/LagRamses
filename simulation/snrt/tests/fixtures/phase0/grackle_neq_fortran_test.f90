program grackle_neq_fortran_test
  use grackle_neq_runtime
  use iso_c_binding
  use omp_lib
  implicit none
  character(1024)::path
  real(c_double)::rho(4),u(4),s(4,9),z(4),length(4),initial(4)
  integer(c_int)::cool(4)
  integer::ierr,i,tile,failed
  call get_command_argument(1,path)
  call omp_set_num_threads(2)
  call grackle_neq_initialize(trim(path),.false.,.false.,0,.false.,1.6726219d-24,1d0,1d0,1d0,ierr)
  if(ierr/=0)error stop 'initialize'
  failed=0
  !$omp parallel do private(rho,u,s,z,length,initial,cool,ierr,i) reduction(+:failed)
  do tile=1,4
     do i=1,4
        rho(i)=10d0**(i-3);z(i)=1d-3*rho(i);u(i)=1d12;initial(i)=u(i)
        length(i)=3d18;cool(i)=mod(i,2)
        call grackle_neq_seed(.76d0*(rho(i)-z(i)),.24d0*(rho(i)-z(i)),2d-4,2d-6,s(i,:))
     enddo
     call grackle_neq_step(4,rho,u,s,z,length,cool,1.6726219d-24,1d0,1d0,1d0,1d8,ierr)
     if(ierr/=0)failed=failed+1
     do i=1,4
        if(cool(i)==0.and.u(i)/=initial(i))failed=failed+1
     enddo
  enddo
  !$omp end parallel do
  if(failed/=0)error stop 'Fortran C ABI or OpenMP mismatch'
  print *, 'PASS: Intel OpenMP / GNU C Grackle ABI, private owners, nine-species update'
end program

! Only satisfy unused native-startup symbols in this bounded bridge fixture.
subroutine units(l,t,d,v,n,t2)
  real(8)::l,t,d,v,n,t2
  error stop 'unexpected native units call'
end subroutine
subroutine clean_stop
  error stop 'unexpected native clean_stop call'
end subroutine

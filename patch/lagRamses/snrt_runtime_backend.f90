! Runtime choice for the production species+dust primary-photon operator.
! Hydro GPU flags are independent. Host state is authoritative at every call;
! automatic fallback occurs only BEFORE launching a device transaction.
module snrt_runtime_backend
  use iso_c_binding, only: c_int,c_float,c_long_long,c_char
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step_species_dust
  implicit none
  private
  public :: snrt_backend_initialize, snrt_runtime_species_dust_step
  integer,save :: mode=0,init_status=0,sharers=1
  integer(c_long_long),save :: min_cells=256
  logical,save :: initialized=.false.,gpu_ready=.false.
  integer,save :: last_choice=-1,cpu_threads=1
  interface
     function configure_cpu(nrank) bind(C,name='snrt_openmp_configure_c') result(nthreads)
       import c_int
       integer(c_int),value::nrank
       integer(c_int)::nthreads
     end function
     function prepare(rank,uuid) bind(C,name='snrt_cuda_prepare_c') result(ierr)
       import c_int,c_char
       integer(c_int),value :: rank
       character(c_char) :: uuid(*)
       integer(c_int) :: ierr
     end function
     function free_bytes() bind(C,name='snrt_cuda_free_bytes_c') result(n)
       import c_long_long
       integer(c_long_long) :: n
     end function
     function cpu_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt) bind(C,name='snrt_openmp_species_dust_c') result(ierr)
       import c_int,c_float
       integer(c_int),value :: no,nw,nd,ng
       real(c_float),value :: cdt
       real(c_float) :: state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in) :: direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in) :: neighbor(*)
       integer(c_int) :: ierr
     end function
  end interface
contains
  subroutine snrt_backend_initialize(ierr)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    integer,intent(out)::ierr
    integer::status,length,local_rank,local_size,comm,i,info
    character(len=64)::value
    character(c_char)::uuid(33)
    character(c_char),allocatable::uuids(:,:)
    if(initialized)then
       ierr=init_status
       return
    endif
    initialized=.true.
    call get_environment_variable('SNRT_BACKEND',value,length=length,status=status)
    if(status==1.or.length==0)value='auto'
    select case(trim(value))
    case('auto');mode=0
    case('openmp');mode=1
    case('cuda');mode=2
    case default;init_status=1
    end select
    if(status/=0.and.status/=1)init_status=1
    call get_environment_variable('SNRT_GPU_MIN_CELLS',value,length=length,status=status)
    if(length>0)then
       if(status/=0)then
          init_status=1
       else
          read(value,*,iostat=status)min_cells
          if(status/=0.or.min_cells<0)init_status=1
       endif
    endif
    local_rank=0;local_size=1
#ifndef WITHOUTMPI
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,comm,info)
    if(info/=0)then
       init_status=1;ierr=init_status
       return
    endif
    call MPI_COMM_RANK(comm,local_rank,info)
    call MPI_COMM_SIZE(comm,local_size,info)
#endif
    cpu_threads=int(configure_cpu(int(local_size,c_int)))
    uuid=' '
    if(mode/=1)gpu_ready=prepare(int(local_rank,c_int),uuid)==0
    if(.not.gpu_ready)uuid=' '
    allocate(uuids(33,local_size))
    uuids(:,1)=uuid
#ifndef WITHOUTMPI
    call MPI_ALLGATHER(uuid,33,MPI_CHARACTER,uuids,33,MPI_CHARACTER,comm,info)
    if(info/=0)init_status=1
    call MPI_COMM_FREE(comm,info)
#endif
    sharers=0
    do i=1,local_size
       if(all(uuids(:,i)==uuid))sharers=sharers+1
    enddo
    sharers=max(1,sharers)
    if(mode==2.and..not.gpu_ready)init_status=2
    ierr=init_status
  end subroutine

  function snrt_runtime_species_dust_step(state,direction,neighbor,tau,stau,dtau,budget, &
       hhe,dust,returned,raw,absorbed_group,absorbed,no,nw,nd,ng,cdt) result(ierr)
    integer(c_int),intent(in)::no,nw,nd,ng
    real(c_float),intent(in)::cdt,direction(*),tau(*),stau(*),dtau(*)
    integer(c_int),intent(in)::neighbor(*)
    real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
    integer(c_int)::ierr
    integer::status,choice
    integer(c_long_long)::required,free
    call snrt_backend_initialize(status)
    ierr=int(status,c_int)
    if(ierr/=0)return
    ! Exact array budget of the CUDA wrapper, plus 64 MiB headroom. Leave
    ! twenty percent free; divide usable memory between ranks sharing a UUID.
    required=4_c_long_long*(2_c_long_long*nw*nd*ng+3_c_long_long*nd+6_c_long_long*no+ &
         12_c_long_long*no*ng+4_c_long_long*no+1)+67108864_c_long_long
    choice=1
    if(mode/=1.and.gpu_ready)then
       free=free_bytes()
       if(real(required,8)<=0.8d0*real(free,8)/sharers)then
          if(mode==2.or.int(no,c_long_long)>=min_cells)choice=2
       else if(mode==2)then
          ierr=5
          return
       endif
    endif
    if(choice/=last_choice)then
       if(choice==2)then
          write(*,'(A,I0,A,I0)')' SNRT backend=CUDA owned_cells=',no,' device_sharers=',sharers
       else
          write(*,'(A,I0,A,I0)')' SNRT backend=OpenMP owned_cells=',no,' threads=',min(no,cpu_threads)
       endif
       last_choice=choice
    endif
    if(choice==2)then
       ierr=snrt_cuda_multigroup_rt_step_species_dust(state,direction,neighbor,tau,stau,dtau, &
            budget,hhe,dust,returned,raw,absorbed_group,absorbed,no,nw,nd,ng,cdt)
       ! Do not replay on CPU after a device error: the enclosing RAMSES
       ! transaction owns rollback, including partial D2H-copy failures.
    else
       ierr=cpu_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
            absorbed_group,absorbed,no,nw,nd,ng,cdt)
    endif
  end function
end module

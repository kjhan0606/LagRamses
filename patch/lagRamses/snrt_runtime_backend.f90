! Runtime choice for the production species+dust primary-photon operator.
! Hydro GPU flags are independent. Host state is authoritative at every call;
! automatic fallback occurs only BEFORE launching a device transaction.
module snrt_runtime_backend
  use iso_c_binding, only: c_int,c_float,c_double,c_long_long,c_char
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step_species_dust
  implicit none
  private
  public :: snrt_backend_initialize, snrt_runtime_species_dust_step, snrt_runtime_dust_material
  public :: snrt_runtime_ir_transport, snrt_runtime_ir_absorb
  public :: snrt_runtime_isotropic_scatter
  public :: snrt_runtime_dust_exchange
  integer,save :: mode=0,init_status=0,sharers=1
  logical,save :: initialized=.false.,gpu_ready=.false.
  integer,save :: last_choice=-1,cpu_threads=1
  integer,save :: dust_mode=0,last_dust_choice=-1,world_rank=0
  interface
     function dust_exchange(input,table,output,nc,nt,dt,floor_t,choice) &
          bind(C,name='snrt_dust_exchange_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::input(*),table(*)
       real(c_double),intent(inout)::output(*)
       integer(c_int),value::nc,nt,choice
       real(c_double),value::dt,floor_t
       integer(c_int)::ierr
     end function
     function isotropic_scatter(state,tau,weight,nc,ng,nd,choice) &
          bind(C,name='snrt_isotropic_scatter_c') result(ierr)
       import c_float,c_double,c_int
       real(c_float),intent(inout)::state(*)
       real(c_double),intent(in)::tau(*),weight(*)
       integer(c_int),value::nc,ng,nd,choice
       integer(c_int)::ierr
     end function
     function ir_transport(energy,ghosts,neighbor,remote,blocked,density,direction,sigma,transported, &
          transmit,loss,response,nc,ng,nd,nghost,cdt,ratio,choice) bind(C,name='snrt_ir_transport_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::energy(*),ghosts(*),density(*),direction(*),sigma(*)
       integer(c_int),intent(in)::neighbor(*),remote(*),blocked(*)
       real(c_double)::transported(*),transmit(*),loss(*),response(*)
       integer(c_int),value::nc,ng,nd,nghost,choice
       real(c_double),value::cdt,ratio
       integer(c_int)::ierr
     end function
     function ir_absorb(transported,transmit,loss,response,rate,weight,candidate,absorbed, &
          nc,ng,nd,dt,sum_w,choice) bind(C,name='snrt_ir_absorb_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::transported(*),transmit(*),loss(*),response(*),rate(*),weight(*)
       real(c_double)::candidate(*),absorbed(*)
       integer(c_int),value::nc,ng,nd,choice
       real(c_double),value::dt,sum_w
       integer(c_int)::ierr
     end function
     function hybrid_configure(rank,streams,cells,threads,share,gpu) bind(C,name='snrt_hybrid_configure_c') result(ierr)
       import c_int
       integer(c_int),value::rank,streams,cells,threads,share,gpu
       integer(c_int)::ierr
     end function
     function hybrid_dust(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,threads) &
          bind(C,name='snrt_hybrid_dust_material_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::input(*),table(*)
       real(c_double)::output(*)
       integer(c_int),value::nc,ng,nt,use_u,threads
       real(c_double),value::dt,background,bath,tolerance
       integer(c_int)::ierr
     end function
     function hybrid_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt) bind(C,name='snrt_hybrid_species_dust_c') result(ierr)
       import c_int,c_float
       integer(c_int),value::no,nw,nd,ng
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       integer(c_int)::ierr
     end function
     function dust_cpu(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,threads) &
          bind(C,name='snrt_dust_material_openmp_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::input(*),table(*)
       real(c_double)::output(*)
       integer(c_int),value::nc,ng,nt,use_u,threads
       real(c_double),value::dt,background,bath,tolerance
       integer(c_int)::ierr
     end function
     function dust_cuda(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,threads) &
          bind(C,name='snrt_dust_material_cuda_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::input(*),table(*)
       real(c_double)::output(*)
       integer(c_int),value::nc,ng,nt,use_u,threads
       real(c_double),value::dt,background,bath,tolerance
       integer(c_int)::ierr
     end function
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
  subroutine snrt_runtime_dust_exchange(gas,heat_capacity,material,density,n_hydrogen, &
       temperature_grid,material_u,area,accommodation,dt,floor_t,temperature,transfer,ierr)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(c_double),intent(in)::gas(:),heat_capacity(:),density(:),n_hydrogen(:)
    real(c_double),intent(in)::temperature_grid(:),material_u(:),area,accommodation,dt,floor_t
    real(c_double),intent(inout)::material(:),temperature(:),transfer(:)
    integer,intent(out)::ierr
    real(c_double),allocatable::input(:,:),coeff(:,:),output(:,:)
    real(c_double),parameter::kb=1.380649d-16,mp=1.67262192369d-24
    integer::nc,nt,status
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    nc=size(gas);nt=size(temperature_grid)
    if(size(heat_capacity)/=nc.or.size(material)/=nc.or.size(density)/=nc.or.size(n_hydrogen)/=nc)return
    if(size(temperature)/=nc.or.size(transfer)/=nc.or.size(material_u)/=nt.or.nt<2)return
    if(.not.ieee_is_finite(area).or.area<0.or..not.ieee_is_finite(accommodation))return
    if(accommodation<0.or.accommodation>1)return
    if(any(.not.ieee_is_finite(gas)).or.any(gas<0))return
    if(any(.not.ieee_is_finite(heat_capacity)).or.any(heat_capacity<=0))return
    if(any(.not.ieee_is_finite(density)).or.any(density<0))return
    if(any(.not.ieee_is_finite(n_hydrogen)).or.any(n_hydrogen<0))return
    if(any(.not.ieee_is_finite(temperature_grid)).or.any(temperature_grid<=0))return
    ierr=0
    if(nc==0)return
    allocate(input(5,nc),coeff(nt,2),output(4,nc),stat=status)
    if(status/=0)then
       ierr=7
       return
    endif
    input(1,:)=gas;input(2,:)=heat_capacity;input(3,:)=material;input(4,:)=density
    ! Hydrogen-equivalent neutral geometric accommodation law, frozen speed.
    ! Not electron/ion Coulomb collisions or a molecular-species network.
    input(5,:)=2*kb*n_hydrogen*density*area*accommodation* &
         sqrt((8*kb/(acos(-1d0)*mp))*(gas/heat_capacity))
    coeff(:,1)=log(temperature_grid);coeff(:,2)=material_u
    ierr=int(dust_exchange(input,coeff,output,int(nc,c_int),int(nt,c_int),dt,floor_t,int(dust_mode,c_int)))
    if(ierr/=0)return
    material=output(2,:);temperature=output(3,:);transfer=output(4,:)
  end subroutine

  subroutine snrt_runtime_isotropic_scatter(state,weight,density,sigma,cdt,ierr)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(c_float),intent(inout),contiguous::state(:,:,:)
    real(c_double),intent(in)::weight(:),density(:),sigma(:),cdt
    integer,intent(out)::ierr
    real(c_double),allocatable::tau(:,:)
    integer::nc,ng,nd,g,allocation_status
    ! Initialization belongs to the collective startup, never a local subset.
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    nd=size(state,1);ng=size(state,2);nc=size(state,3)
    if(nd/=size(weight).or.ng/=size(sigma).or.nc/=size(density))return
    if(nd<1.or.ng<1)return
    if(.not.ieee_is_finite(cdt).or.cdt<0)return
    if(any(.not.ieee_is_finite(density)).or.any(density<0))return
    if(any(.not.ieee_is_finite(sigma)).or.any(sigma<0))return
    ierr=0
    if(nc==0)return
    allocate(tau(nc,ng),stat=allocation_status)
    if(allocation_status/=0)then
       ierr=7
       return
    endif
    do g=1,ng
       tau(:,g)=density*sigma(g)*cdt
    enddo
    ! The C++ batch transaction publishes only after every CPU/GPU batch passes.
    ierr=int(isotropic_scatter(state,tau,weight,int(nc,c_int),int(ng,c_int), &
         int(nd,c_int),int(mode,c_int)))
  end subroutine

  subroutine snrt_runtime_ir_transport(energy,ghosts,neighbor,remote,blocked,density,direction,sigma, &
       cdt,ratio,transported,transmit,loss,response,ierr)
    real(c_double),intent(in)::energy(:,:,:),ghosts(:,:,:),density(:),direction(:,:),sigma(:),cdt,ratio
    integer,intent(in)::neighbor(:,:),remote(:,:)
    logical,intent(in)::blocked(:,:)
    real(c_double),intent(out)::transported(:,:,:),transmit(:,:),loss(:,:),response(:,:)
    integer,intent(out)::ierr
    integer(c_int),allocatable::flags(:,:)
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    flags=merge(1_c_int,0_c_int,blocked)
    ierr=int(ir_transport(energy,ghosts,int(neighbor,c_int),int(remote,c_int),flags,density,direction,sigma, &
         transported,transmit,loss,response,int(size(density),c_int),int(size(energy,1),c_int), &
         int(size(energy,2),c_int),int(size(ghosts,3),c_int),cdt,ratio,int(dust_mode,c_int)))
  end subroutine

  subroutine snrt_runtime_ir_absorb(transported,transmit,loss,response,rate,weight,dt,sum_w,candidate,absorbed,ierr)
    real(c_double),intent(in)::transported(:,:,:),transmit(:,:),loss(:,:),response(:,:),rate(:,:),weight(:),dt,sum_w
    real(c_double),intent(out)::candidate(:,:,:),absorbed(:)
    integer,intent(out)::ierr
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    ierr=int(ir_absorb(transported,transmit,loss,response,rate,weight,candidate,absorbed, &
         int(size(transported,3),c_int),int(size(transported,1),c_int),int(size(transported,2),c_int), &
         dt,sum_w,int(dust_mode,c_int)))
  end subroutine

  subroutine snrt_backend_initialize(ierr,nstreams)
#ifndef WITHOUTMPI
    use mpi_mod
#endif
    integer,intent(out)::ierr
    integer,intent(in),optional::nstreams
    integer::stream_count,batch_cells
    integer::status,length,local_rank,local_size,comm,i,info
    character(len=64)::value
    character(c_char)::uuid(33)
    character(c_char),allocatable::uuids(:,:)
    if(initialized)then
       ierr=init_status
       return
    endif
    initialized=.true.
    stream_count=1
    if(present(nstreams))stream_count=nstreams
    batch_cells=256
    call get_environment_variable('SNRT_HYBRID_BATCH_CELLS',value,length=length,status=status)
    if(length>0)then
       if(status/=0)then
          init_status=1
       else
          read(value,*,iostat=status)batch_cells
          if(status/=0.or.batch_cells<1.or.batch_cells>1048576)init_status=1
       endif
    endif
    call get_environment_variable('SNRT_BACKEND',value,length=length,status=status)
    if(status==1.or.length==0)value='auto'
    select case(trim(value))
    case('auto');mode=0
    case('openmp');mode=1
    case('cuda');mode=2
    case default;init_status=1
    end select
    if(status/=0.and.status/=1)init_status=1
    dust_mode=mode
    call get_environment_variable('SNRT_DUST_BACKEND',value,length=length,status=status)
    if(status/=0.and.status/=1)init_status=1
    if(length>0.and.status==0)then
       select case(trim(value))
       case('auto');dust_mode=0
       case('openmp');dust_mode=1
       case('cuda');dust_mode=2
       case default;init_status=1
       end select
    endif
    local_rank=0;local_size=1
#ifndef WITHOUTMPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,world_rank,info)
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
    if(mode/=1.or.dust_mode/=1)gpu_ready=prepare(int(local_rank,c_int),uuid)==0
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
    if((mode==2.or.dust_mode==2).and..not.gpu_ready)init_status=2
    if(init_status==0)then
       status=int(hybrid_configure(int(local_rank,c_int),int(stream_count,c_int),int(batch_cells,c_int), &
            int(cpu_threads,c_int),int(sharers,c_int),merge(1_c_int,0_c_int,gpu_ready)))
       if(status/=0)init_status=status
    endif
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
    if(mode==2.and.gpu_ready)then
       free=free_bytes()
       if(real(required,8)<=0.8d0*real(free,8)/sharers)then
          choice=2
       else if(mode==2)then
          ierr=5
          return
       endif
    endif
    if(mode==0)choice=3
    if(choice/=last_choice)then
       if(choice==3)then
          write(*,'(A,I0)')' SNRT backend=hybrid owned_cells=',no
       else if(choice==2)then
          write(*,'(A,I0,A,I0)')' SNRT backend=CUDA owned_cells=',no,' device_sharers=',sharers
       else
          write(*,'(A,I0,A,I0)')' SNRT backend=OpenMP owned_cells=',no,' threads=',min(no,cpu_threads)
       endif
       last_choice=choice
    endif
    if(choice==3)then
       ierr=hybrid_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
            absorbed_group,absorbed,no,nw,nd,ng,cdt)
    else if(choice==2)then
       ierr=snrt_cuda_multigroup_rt_step_species_dust(state,direction,neighbor,tau,stau,dtau, &
            budget,hhe,dust,returned,raw,absorbed_group,absorbed,no,nw,nd,ng,cdt)
       ! Do not replay on CPU after a device error: the enclosing RAMSES
       ! transaction owns rollback, including partial D2H-copy failures.
    else
       ierr=cpu_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
            absorbed_group,absorbed,no,nw,nd,ng,cdt)
    endif
  end function

  subroutine snrt_runtime_dust_material(heating,density,old_energy,capacity,log_t,power,band, &
       material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
       gas_energy,gas_capacity,conductance,gas_transfer)
    real(c_double),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
    real(c_double),intent(in)::material_u(:),dt,background,bath,tolerance
    logical,intent(in)::use_u
    real(c_double),intent(out)::rate(:,:),temperature(:),next_energy(:)
    integer,intent(out)::ierr
    real(c_double),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(c_double),optional,intent(out)::gas_transfer(:)
    real(c_double),allocatable::input(:),coefficients(:),output(:)
    integer::nc,ng,nt,choice,extra,material_mode
    integer(c_long_long)::required,free
    ! Initialization is collective and belongs to read_params, NOT a subset
    ! of ranks that happens to have local dust cells at this AMR level.
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    nc=size(heating);ng=size(band,1);nt=size(log_t)
    if(nc<1.or.ng<1.or.nt<2)return
    extra=0;material_mode=merge(1,0,use_u)
    if(present(gas_energy))then
       if(.not.present(gas_capacity).or..not.present(conductance).or..not.present(gas_transfer))return
       if(.not.use_u)return
       if(size(gas_energy)/=nc.or.size(gas_capacity)/=nc.or.size(conductance)/=nc.or.size(gas_transfer)/=nc)return
       extra=1;material_mode=3 ! Joint BE gas/dust/IR with K(Tgas), frozen chemistry/Cv
    else if(present(gas_capacity).or.present(conductance).or.present(gas_transfer))then
       return
    endif
    required=8_c_long_long*((int(ng,c_long_long)+6)*nc+(int(ng,c_long_long)+3)*nt)+16777216_c_long_long
    required=required+32_c_long_long*extra*nc
    choice=1
    if(dust_mode==2.and.gpu_ready)then
       free=free_bytes()
       if(real(required,8)<=0.8d0*real(free,8)/sharers)then
          choice=2
       else if(dust_mode==2)then
          return
       endif
    endif
    if(dust_mode==0)choice=3
    if(choice/=last_dust_choice)then
       if(choice==3)then
          write(*,'(A,I0,A,I0)')' SNRT dust material backend=hybrid rank=',world_rank,' cells=',nc
       else if(choice==2)then
          write(*,'(A,I0,A,I0,A,I0)')' SNRT dust material backend=CUDA rank=',world_rank, &
               ' cells=',nc,' device_sharers=',sharers
       else
          write(*,'(A,I0,A,I0,A,I0)')' SNRT dust material backend=OpenMP rank=',world_rank, &
               ' cells=',nc,' threads=',min(nc,cpu_threads)
       endif
       last_dust_choice=choice
    endif
    input=[heating,density,old_energy,capacity]
    if(extra==1)input=[input,gas_energy,gas_capacity,conductance]
    coefficients=[log_t,power,material_u,reshape(band,[ng*nt])]
    allocate(output((ng+2+extra)*nc))
    if(choice==3)then
       ierr=int(hybrid_dust(input,coefficients,output,int(nc,c_int),int(ng,c_int),int(nt,c_int), &
            int(material_mode,c_int),dt,background,bath,tolerance,int(cpu_threads,c_int)))
    else if(choice==2)then
       ierr=int(dust_cuda(input,coefficients,output,int(nc,c_int),int(ng,c_int),int(nt,c_int), &
            int(material_mode,c_int),dt,background,bath,tolerance,int(cpu_threads,c_int)))
    else
       ierr=int(dust_cpu(input,coefficients,output,int(nc,c_int),int(ng,c_int),int(nt,c_int), &
            int(material_mode,c_int),dt,background,bath,tolerance,int(cpu_threads,c_int)))
    endif
    if(ierr/=0)return ! no after-launch fallback; outer transaction owns rollback
    rate=reshape(output(1:ng*nc),[ng,nc])
    temperature=output(ng*nc+1:(ng+1)*nc)
    next_energy=output((ng+1)*nc+1:(ng+2)*nc)
    if(extra==1)gas_transfer=output((ng+2)*nc+1:)
  end subroutine
end module

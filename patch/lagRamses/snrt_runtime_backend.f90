! Runtime choice for the production species+dust primary-photon operator.
! Hydro GPU flags are independent. Host state is authoritative at every call;
! automatic fallback occurs only BEFORE launching a device transaction.
module snrt_runtime_backend
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_c_binding, only: c_int,c_float,c_double,c_long_long,c_char,c_ptr,c_loc,c_null_ptr,c_funptr,c_funloc
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step_species_dust
  use snrt_spectral_contract, only: snrt_band_enabled,snrt_group_edges_ev,snrt_node_secondaries_enabled,snrt_band_kind
  use snrt_thermochemistry, only: snrt_secondary_fractions_c
  use snrt_spectral_contract, only: snrt_d03_band_enabled,snrt_fe_band_enabled,snrt_grain_band_bins
  use dust_composition_optics, only: d03_band_ev,d03_band_abs,d03_band_transport
  use dust_iron_optics, only: fe_band_ev,fe_six_band_abs,fe_six_band_transport
  implicit none
  private
  public :: snrt_backend_initialize, snrt_runtime_species_dust_step, snrt_runtime_dust_material
  public :: snrt_runtime_ir_transport, snrt_runtime_ir_absorb
  public :: snrt_runtime_isotropic_scatter,snrt_runtime_ir_scatter
  public :: snrt_runtime_dust_exchange
  public :: snrt_runtime_cpu_material_allowed
  public :: snrt_runtime_energy_admit
  integer,save :: mode=0,init_status=0,sharers=1
  logical,save :: initialized=.false.,gpu_ready=.false.
  integer,save :: last_choice=-1,cpu_threads=1
  integer,save :: dust_mode=0,last_dust_choice=-1,world_rank=0
  interface
     function cpu_band_d03(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt,shift,reference,hhe_e,dust_e,energy_moment,columns,edges, &
          callback,xi,deposition,grains,kabs,ksca,node_ev,weights,nb) bind(C,name='snrt_openmp_band_grains_c') result(ierr)
       import c_int,c_float,c_double,c_funptr
       type(c_funptr),value::callback
       integer(c_int),value::no,nw,nd,ng,nb
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       real(c_double)::shift(*),hhe_e(*),dust_e(*),energy_moment(*),deposition(*)
       real(c_double),intent(in)::reference(*),columns(*),edges(*),xi(*),grains(*),kabs(*),ksca(*),node_ev(*),weights(*)
       integer(c_int)::ierr
     end function
     function cpu_band_secondary(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment,shift,reference,hhe_e,dust_e,energy_moment,columns,edges, &
          callback,xi,deposition) bind(C,name='snrt_openmp_band_secondary_c') result(ierr)
       import c_int,c_float,c_double,c_ptr,c_funptr
       type(c_ptr),value::moment
       type(c_funptr),value::callback
       integer(c_int),value::no,nw,nd,ng
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       real(c_double)::shift(*),hhe_e(*),dust_e(*),energy_moment(*),deposition(*)
       real(c_double),intent(in)::reference(*),columns(*),edges(*),xi(*)
       integer(c_int)::ierr
     end function
     function cpu_band_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment,shift,reference,hhe_e,dust_e,energy_moment,columns,edges) &
          bind(C,name='snrt_openmp_band_energy_c') result(ierr)
       import c_int,c_float,c_double,c_ptr
       type(c_ptr),value::moment
       integer(c_int),value::no,nw,nd,ng
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       real(c_double)::shift(*),hhe_e(*),dust_e(*),energy_moment(*)
       real(c_double),intent(in)::reference(*),columns(*),edges(*)
       integer(c_int)::ierr
     end function
     function cpu_energy_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment,shift,reference,hhe_e,dust_e,energy_moment) &
          bind(C,name='snrt_openmp_species_dust_energy_c') result(ierr)
       import c_int,c_float,c_double,c_ptr
       type(c_ptr),value::moment
       integer(c_int),value::no,nw,nd,ng
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       real(c_double)::shift(*),hhe_e(*),dust_e(*),energy_moment(*)
       real(c_double),intent(in)::reference(*)
       integer(c_int)::ierr
     end function
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
     function ir_scatter(state,tau,weight,nc,ng,nd,choice) bind(C,name='snrt_ir_scatter_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(inout)::state(*)
       real(c_double),intent(in)::tau(*),weight(*)
       integer(c_int),value::nc,ng,nd,choice
       integer(c_int)::ierr
     end function
     function ir_transport(energy,ghosts,neighbor,remote,blocked,density,direction,sigma,transported, &
          transmit,loss,response,nc,ng,nd,nghost,cdt,ratio,choice,cell_sigma) bind(C,name='snrt_ir_transport_c') result(ierr)
       import c_double,c_int
       real(c_double),intent(in)::energy(*),ghosts(*),density(*),direction(*),sigma(*)
       integer(c_int),intent(in)::neighbor(*),remote(*),blocked(*)
       real(c_double)::transported(*),transmit(*),loss(*),response(*)
       integer(c_int),value::nc,ng,nd,nghost,choice,cell_sigma
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
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment) bind(C,name='snrt_hybrid_species_dust_moment_c') result(ierr)
       import c_int,c_float,c_ptr
       type(c_ptr),value::moment
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
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment) bind(C,name='snrt_openmp_species_dust_moment_c') result(ierr)
       import c_int,c_float,c_ptr
       type(c_ptr),value::moment
       integer(c_int),value :: no,nw,nd,ng
       real(c_float),value :: cdt
       real(c_float) :: state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in) :: direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in) :: neighbor(*)
       integer(c_int) :: ierr
     end function
     function cuda_moment_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
          absorbed_group,absorbed,no,nw,nd,ng,cdt,moment) bind(C,name='snrt_cuda_species_dust_moment_c') result(ierr)
       import c_int,c_float,c_ptr
       type(c_ptr),value::moment
       integer(c_int),value::no,nw,nd,ng
       real(c_float),value::cdt
       real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
       real(c_float),intent(in)::direction(*),tau(*),stau(*),dtau(*)
       integer(c_int),intent(in)::neighbor(*)
       integer(c_int)::ierr
     end function
  end interface
contains
  subroutine snrt_runtime_energy_admit(ierr)
    integer,intent(out)::ierr
    call snrt_backend_initialize(ierr)
    ! Paired correction transport is host-only for now. This is admission,
    ! before any halo/transport/device work, including ranks with zero leaves.
    if(mode==2)then
       ierr=8
       if(world_rank==0)write(*,'(A)')' SNRT paired energy transport: forced CUDA unsupported (status 8)'
    endif
  end subroutine

  logical function snrt_runtime_cpu_material_allowed() result(allowed)
    allowed=initialized.and.init_status==0.and.dust_mode/=2
  end function

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

  subroutine snrt_runtime_isotropic_scatter(state,weight,density,sigma,cdt,ierr,cell_sigma)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(c_float),intent(inout),contiguous::state(:,:,:)
    real(c_double),intent(in)::weight(:),density(:),sigma(:),cdt
    integer,intent(out)::ierr
    real(c_double),allocatable::tau(:,:)
    real(c_double),optional,intent(in)::cell_sigma(:,:)
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
    if(present(cell_sigma))then
       if(any(shape(cell_sigma)/=[ng,nc]))return
       if(any(.not.ieee_is_finite(cell_sigma)).or.any(cell_sigma<0))return
    endif
    ierr=0
    if(nc==0)return
    allocate(tau(nc,ng),stat=allocation_status)
    if(allocation_status/=0)then
       ierr=7
       return
    endif
    do g=1,ng
       tau(:,g)=density*sigma(g)*cdt
       if(present(cell_sigma))tau(:,g)=density*cell_sigma(g,:)*cdt
    enddo
    ! The C++ batch transaction publishes only after every CPU/GPU batch passes.
    ierr=int(isotropic_scatter(state,tau,weight,int(nc,c_int),int(ng,c_int), &
         int(nd,c_int),int(mode,c_int)))
  end subroutine

  subroutine snrt_runtime_ir_scatter(state,weight,density,cell_sigma,cdt,ierr)
    real(c_double),intent(inout),contiguous::state(:,:,:)
    real(c_double),intent(in)::weight(:),density(:),cell_sigma(:,:),cdt
    integer,intent(out)::ierr
    real(c_double),allocatable::tau(:,:)
    integer::nc,ng,nd,g,allocation_status
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    ng=size(state,1);nd=size(state,2);nc=size(state,3)
    if(nd/=size(weight).or.nc/=size(density))return
    if(ng<1.or.nd<1.or.any(shape(cell_sigma)/=[ng,nc]))return
    if(.not.ieee_is_finite(cdt).or.cdt<0)return
    if(any(.not.ieee_is_finite(density)).or.any(density<0))return
    if(any(.not.ieee_is_finite(cell_sigma)).or.any(cell_sigma<0))return
    ierr=0
    if(nc==0)return
    allocate(tau(nc,ng),stat=allocation_status)
    if(allocation_status/=0)then
       ierr=7
       return
    endif
    do g=1,ng
       tau(:,g)=density*cell_sigma(g,:)*cdt
    enddo
    ierr=int(ir_scatter(state,tau,weight,int(nc,c_int),int(ng,c_int),int(nd,c_int),int(dust_mode,c_int)))
  end subroutine

  subroutine snrt_runtime_ir_transport(energy,ghosts,neighbor,remote,blocked,density,direction,sigma, &
       cdt,ratio,transported,transmit,loss,response,ierr,cell_sigma)
    real(c_double),intent(in)::energy(:,:,:),ghosts(:,:,:),density(:),direction(:,:),sigma(:),cdt,ratio
    integer,intent(in)::neighbor(:,:),remote(:,:)
    logical,intent(in)::blocked(:,:)
    real(c_double),intent(out)::transported(:,:,:),transmit(:,:),loss(:,:),response(:,:)
    integer,intent(out)::ierr
    integer(c_int),allocatable::flags(:,:)
    real(c_double),optional,intent(in)::cell_sigma(:,:)
    real(c_double),allocatable::alpha(:,:)
    integer::g
    ierr=7
    if(.not.initialized.or.init_status/=0)return
    flags=merge(1_c_int,0_c_int,blocked)
    if(present(cell_sigma))then
       if(any(shape(cell_sigma)/=[size(sigma),size(density)]))return
       if(any(.not.ieee_is_finite(cell_sigma)).or.any(cell_sigma<0))return
       allocate(alpha(size(sigma),size(density)))
       do g=1,size(sigma)
          alpha(g,:)=cell_sigma(g,:)*density
       enddo
       ierr=int(ir_transport(energy,ghosts,int(neighbor,c_int),int(remote,c_int),flags,alpha,direction,sigma, &
            transported,transmit,loss,response,int(size(density),c_int),int(size(energy,1),c_int), &
            int(size(energy,2),c_int),int(size(ghosts,3),c_int),cdt,ratio,int(dust_mode,c_int),1_c_int))
       return
    endif
    ierr=int(ir_transport(energy,ghosts,int(neighbor,c_int),int(remote,c_int),flags,density,direction,sigma, &
         transported,transmit,loss,response,int(size(density),c_int),int(size(energy,1),c_int), &
         int(size(energy,2),c_int),int(size(ghosts,3),c_int),cdt,ratio,int(dust_mode,c_int),0_c_int))
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
    integer::stream_count,batch_cells,band_local,band_min,band_max
    integer::status,length,local_rank,local_size,comm,i,info
    character(len=64)::value
    character(c_char)::uuid(33)
    character(c_char),allocatable::uuids(:,:)
    if(initialized)then
       ierr=init_status
       return
    endif
    initialized=.true.
    band_local=snrt_band_kind();band_min=band_local;band_max=band_local
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(band_local,band_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
    if(info/=0)init_status=1
    call MPI_ALLREDUCE(band_local,band_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)init_status=1
#endif
    if(band_min/=band_max)then
       write(*,*)'SNRT_SPECTRAL_MODEL must agree across MPI ranks'
       init_status=1
    endif
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
       hhe,dust,returned,raw,absorbed_group,absorbed,no,nw,nd,ng,cdt,dust_moment, &
       shift,reference_ev,hhe_energy,dust_energy,dust_energy_moment,species_columns,secondary_xi,band_deposition, &
       grain_columns,angular_weights) result(ierr)
    integer(c_int),intent(in)::no,nw,nd,ng
    real(c_float),intent(in)::cdt,direction(*),tau(*),stau(*),dtau(*)
    integer(c_int),intent(in)::neighbor(*)
    real(c_float)::state(*),budget(*),hhe(*),dust(*),returned(*),raw(*),absorbed_group(*),absorbed(*)
    integer(c_int)::ierr
    integer::status,choice
    integer(c_long_long)::required,free
    real(c_double),optional,intent(inout)::dust_moment(:,:,:)
    real(c_double),optional,intent(inout)::shift(:,:,:),hhe_energy(:,:,:),dust_energy(:,:),dust_energy_moment(:,:,:)
    real(c_double),optional,intent(in)::reference_ev(:)
    real(c_double),optional,intent(in)::species_columns(:,:)
    real(c_double),optional,intent(in)::secondary_xi(:)
    real(c_double),optional,intent(inout)::band_deposition(:,:)
    real(c_double),optional,intent(in)::grain_columns(:,:),angular_weights(:)
    real(c_double),allocatable,target::moment_stage(:,:,:)
    type(c_ptr)::moment_pointer
    moment_pointer=c_null_ptr;ierr=1
    if(present(species_columns).neqv.snrt_band_enabled())return
    if(present(grain_columns).neqv.snrt_d03_band_enabled())return
    if(present(grain_columns).neqv.present(angular_weights))return
    if(present(grain_columns))then
       if(any(shape(grain_columns)/=[no,snrt_grain_band_bins()]).or. &
            size(angular_weights)/=nd.or.present(dust_moment))return
       if(snrt_fe_band_enabled())then
          if(any(fe_band_ev/=d03_band_ev))return
       endif
    endif
    if(present(secondary_xi).neqv.snrt_node_secondaries_enabled())return
    if(present(secondary_xi).neqv.present(band_deposition))return
    if(present(secondary_xi))then
       if(size(secondary_xi)/=no.or.any(shape(band_deposition)/=[no,8]))return
    endif
    if(present(species_columns))then
       if(.not.present(shift))return
       if(any(shape(species_columns)/=[no,3]))return
    endif
    if(present(dust_moment))then
       if(any(shape(dust_moment)/=[no,ng,3]))return
       allocate(moment_stage(no,ng,3));moment_stage=0
       if(no>0.and.ng>0)moment_pointer=c_loc(moment_stage(1,1,1))
    endif
    call snrt_backend_initialize(status)
    if(present(shift))call snrt_runtime_energy_admit(status)
    ierr=int(status,c_int)
    if(ierr/=0)return
    if(present(shift))then
       ierr=1
       if(.not.present(reference_ev).or..not.present(hhe_energy).or..not.present(dust_energy).or. &
            .not.present(dust_energy_moment))return
       if(any(shape(shift)/=[nw,nd,ng]).or.size(reference_ev)/=ng)return
       if(any(shape(hhe_energy)/=[no,ng,3]).or.any(shape(dust_energy)/=[no,ng]).or. &
            any(shape(dust_energy_moment)/=[no,ng,3]))return
       call snrt_runtime_energy_admit(status)
       ierr=int(status,c_int)
       if(ierr/=0)return
       if(last_choice/=4)then
          write(*,'(A,I0)')' SNRT paired energy transport backend=OpenMP owned_cells=',no
          last_choice=4
       endif
       if(present(species_columns))then
          if(present(grain_columns))then
             if(snrt_fe_band_enabled())then
                ierr=cpu_band_d03(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
                     absorbed_group,absorbed,no,nw,nd,ng,cdt,shift,reference_ev,hhe_energy,dust_energy,dust_energy_moment, &
                     species_columns,snrt_group_edges_ev,c_funloc(snrt_secondary_fractions_c),secondary_xi,band_deposition, &
                     grain_columns,fe_six_band_abs,fe_six_band_transport,fe_band_ev,angular_weights,6_c_int)
             else
             ierr=cpu_band_d03(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
                  absorbed_group,absorbed,no,nw,nd,ng,cdt,shift,reference_ev,hhe_energy,dust_energy,dust_energy_moment, &
                  species_columns,snrt_group_edges_ev,c_funloc(snrt_secondary_fractions_c),secondary_xi,band_deposition, &
                  grain_columns,d03_band_abs,d03_band_transport,d03_band_ev,angular_weights,4_c_int)
             endif
          else if(present(secondary_xi))then
             ierr=cpu_band_secondary(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
                  absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer,shift,reference_ev, &
                  hhe_energy,dust_energy,dust_energy_moment,species_columns,snrt_group_edges_ev, &
                  c_funloc(snrt_secondary_fractions_c),secondary_xi,band_deposition)
          else
          ierr=cpu_band_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
               absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer,shift,reference_ev, &
               hhe_energy,dust_energy,dust_energy_moment,species_columns,snrt_group_edges_ev)
          endif
       else
       ierr=cpu_energy_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
            absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer,shift,reference_ev, &
            hhe_energy,dust_energy,dust_energy_moment)
       endif
       if(present(dust_moment).and.ierr==0)dust_moment=moment_stage
       return
    else if(present(reference_ev).or.present(hhe_energy).or.present(dust_energy).or.present(dust_energy_moment))then
       ierr=1;return
    endif
    ! Exact array budget of the CUDA wrapper, plus 64 MiB headroom. Leave
    ! twenty percent free; divide usable memory between ranks sharing a UUID.
    required=4_c_long_long*(2_c_long_long*nw*nd*ng+3_c_long_long*nd+6_c_long_long*no+ &
         12_c_long_long*no*ng+4_c_long_long*no+1)+67108864_c_long_long
    if(present(dust_moment))required=required+24_c_long_long*no*ng
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
            absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer)
    else if(choice==2)then
       ierr=cuda_moment_step(state,direction,neighbor,tau,stau,dtau, &
            budget,hhe,dust,returned,raw,absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer)
       ! Do not replay on CPU after a device error: the enclosing RAMSES
       ! transaction owns rollback, including partial D2H-copy failures.
    else
       ierr=cpu_step(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw, &
            absorbed_group,absorbed,no,nw,nd,ng,cdt,moment_pointer)
    endif
    if(present(dust_moment).and.ierr==0)then
       if(any(.not.ieee_is_finite(moment_stage)))then
          ierr=3 ! Enclosing photon transaction owns rollback on device failure.
       else
          dust_moment=moment_stage
       endif
    endif
  end function

  subroutine snrt_runtime_dust_material(heating,density,old_energy,capacity,log_t,power,band, &
       material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
       gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
    real(c_double),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
    real(c_double),intent(in)::material_u(:),dt,background,bath,tolerance
    logical,intent(in)::use_u
    real(c_double),intent(out)::rate(:,:),temperature(:),next_energy(:)
    integer,intent(out)::ierr
    real(c_double),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(c_double),optional,intent(out)::gas_transfer(:)
    real(c_double),optional,intent(in)::cell_material_u(:,:)
    real(c_double),optional,intent(in)::cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
    real(c_double),allocatable::input(:),coefficients(:),output(:)
    integer::nc,ng,nt,choice,extra,material_mode,s
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
    if(present(cell_material_u))then
       if(.not.use_u.or.any(shape(cell_material_u)/=[nt,nc]))return
       if(any(.not.ieee_is_finite(cell_material_u)).or.any(cell_material_u<=0))return
       if(any(cell_material_u(2:,:)<=cell_material_u(:nt-1,:)))return
       material_mode=material_mode+4
    endif
    required=8_c_long_long*((int(ng,c_long_long)+6)*nc+(int(ng,c_long_long)+3)*nt)+16777216_c_long_long
    required=required+32_c_long_long*extra*nc
    if(present(cell_material_u))required=required+8_c_long_long*nt*nc
    if(present(cell_weights))then
       if(.not.present(basis_power).or..not.present(basis_band).or..not.use_u.or.nt>256)return
       if(any(shape(cell_weights)/=[4,nc]).or.any(shape(basis_power)/=[nt,4]))return
       if(any(shape(basis_band)/=[ng,nt,4]))return
       if(any(.not.ieee_is_finite(cell_weights)).or.any(cell_weights<0))return
       if(any(abs(sum(cell_weights,dim=1)-1)>1d-12))return
       if(any(.not.ieee_is_finite(basis_power)).or.any(basis_power<=0))return
       if(any(basis_power(2:,:)<=basis_power(:nt-1,:)))return
       if(any(.not.ieee_is_finite(basis_band)).or.any(basis_band<0))return
       if(any(basis_band(:,2:,:)<basis_band(:,:nt-1,:)))return
       if(any(abs(sum(basis_band,dim=1)-basis_power)>1d-12*basis_power))return
       material_mode=material_mode+8
       required=required+8_c_long_long*(4_c_long_long*nc+4_c_long_long*(ng+1)*nt)
    else if(present(basis_power).or.present(basis_band))then
       return
    endif
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
    if(present(cell_material_u))input=[input,reshape(transpose(cell_material_u),[nt*nc])]
    if(present(cell_weights))input=[input,reshape(transpose(cell_weights),[4*nc])]
    coefficients=[log_t,power,material_u,reshape(band,[ng*nt])]
    if(present(cell_weights))then
       do s=1,4
          coefficients=[coefficients,basis_power(:,s),reshape(basis_band(:,:,s),[ng*nt])]
       enddo
    endif
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

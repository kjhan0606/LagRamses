! Dark (radiation-transport-off) material evolution for the admitted live
! dust+CHIMES profile.  The RT solvers keep their existing coupled packet
! transaction; this level operator reuses the same CHIMES and dust-material
! kernels with an exactly empty photon field and an analytic CMB bath.
module snrt_dark_matter_runtime
  use iso_c_binding, only: c_int,c_double,c_ptr,c_funptr,c_loc,c_funloc,c_f_pointer,c_null_ptr
  use amr_commons
  use amr_parameters, only: radiation_transport
  use hydro_commons, only: uold,magnetic_energy
  use hydro_parameters, only: gamma,nvar,nvar_all,inener,idust,idust_energy,idust_species,idust_bins,ichimes,ichem
  use dust_mass_physics, only: dust_chimes_enabled,dust_optics_enabled,dust_material_composition_enabled
  use snrt_dust_contract
  use snrt_dust_live, only: snrt_dust_live_prepare,snrt_dust_live_stage,dust_live_coarse_trial
  use snrt_dust_ir, only: dust_ir_diagnostics,snrt_dust_material_temperature
  use dust_composition_material, only: dust_composition_curve,dust_composition_area
  use dust_composition_optics, only: d03_cell_weights
  use snrt_chimes_runtime, only: chimes_live_initialize,chimes_live_capacity, &
       chimes_live_stage,chimes_live_band_stage,chimes_live_cold_stage
  use snrt_chimes, only: chimes_ns
  use snrt_spectral_contract, only: snrt_chimes_band_enabled,snrt_chimes_cold_enabled
  use snrt_atomic_cooling, only: atomic_mh
  use snrt_agn_source, only: snrt_c_cgs
  use omp_lib, only: omp_get_wtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: snrt_dark_matter_advance_level
#ifdef SNRT_CHIMES_CUDA
  ! Each callback writes disjoint staging slots. No AMR publication or MPI is
  ! allowed inside a continuation; the existing level transaction commits.
  type,bind(C)::dark_tile_context
     integer(c_int)::count,width,first,offset,level_queue
     real(c_double)::sd,sv,dt,dx
     type(c_ptr)::cells,chemical,gas_total,hydrogen,td,status
  end type
  interface
     integer(c_int) function dark_cuda_requested() bind(C,name='snrt_chimes_cuda_requested')
       import c_int
     end function
     integer(c_int) function dark_level_requested() bind(C,name='snrt_chimes_rhs_level_requested')
       import c_int
     end function
     integer(c_int) function dark_level_run(n,callback,context) bind(C,name='snrt_chimes_rhs_level_run')
       import c_int,c_funptr,c_ptr
       integer(c_int),value::n
       type(c_funptr),value::callback
       type(c_ptr),value::context
     end function
     integer(c_int) function dark_tile_run(n,callback,context) bind(C,name='snrt_chimes_rhs_tile_run')
       import c_int,c_funptr,c_ptr
       integer(c_int),value::n
       type(c_funptr),value::callback
       type(c_ptr),value::context
     end function
  end interface
#endif
contains
#ifdef SNRT_CHIMES_CUDA
  subroutine dark_tile_cell(index,context) bind(C)
    integer(c_int),value::index
    type(c_ptr),value::context
    type(dark_tile_context),pointer::work
    integer(c_int),pointer::cells(:),status(:)
    real(c_double),pointer::chemical(:,:),gas_total(:),hydrogen(:),td(:)
    integer::j,k,cell,ierr
    real(dp)::row(nvar_all),state(chimes_ns),energy,number(1,9),radiation(1,9)
    real(dp)::next_n(1,9),next_e(1,9),ledger(11),grain_n(9),grain_e(9),cell_td,cell_density
    real(dp)::cell_curve(snrt_dust_contract_number_temperature)
    call c_f_pointer(context,work)
    call c_f_pointer(work%cells,cells,[work%count])
    call c_f_pointer(work%chemical,chemical,[chimes_ns,work%count])
    call c_f_pointer(work%gas_total,gas_total,[work%count])
    call c_f_pointer(work%hydrogen,hydrogen,[work%count])
    if(work%level_queue==0)call c_f_pointer(work%td,td,[work%width])
    call c_f_pointer(work%status,status,[work%width])
    j=work%offset+index;k=work%first+j-1;cell=cells(k)
    row=uold(cell,1:nvar_all);state=row(ichimes:ichimes+chimes_ns-1);energy=row(ndim+2)
    if(work%level_queue/=0)then
       status(j)=1
       if(.not.all(ieee_is_finite(row)).or.row(1)<=0)return
       cell_density=row(idust)*work%sd/snrt_dust_contract_mass_per_h_g
       if(cell_density>0)then
          call dust_composition_curve(snrt_dust_contract_temperature_k(1:size(cell_curve)), &
               row(idust_species:idust_species+1),snrt_dust_contract_mass_per_h_g,cell_curve,ierr)
          if(ierr/=0)return
          call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:size(cell_curve)), &
               cell_curve,row(idust_energy)*work%sv**2/cell_density,cell_td,ierr)
          if(ierr/=0)return
       else
          cell_td=2.727d0/aexp
       endif
    else
       cell_td=td(j)
    endif
    number=0;radiation=0;next_n=0;next_e=0;ledger=0;grain_n=0;grain_e=0
    call chimes_live_cold_stage(cell,work%sd,work%sv,work%dt,work%dx,cell_td,.01d0,1, &
         number,radiation,state,energy,next_n,next_e,ledger,grain_n,grain_e,ierr,staged_row=row)
    status(j)=ierr
    if(ierr/=0)return
    if(any(next_n/=0d0).or.any(next_e/=0d0))then
       status(j)=1;return
    endif
    chemical(:,k)=state;gas_total(k)=energy;hydrogen(k)=max(row(ichem)*work%sd/atomic_mh,0d0)
  end subroutine dark_tile_cell
#endif
  subroutine snrt_dark_matter_advance_level(ilevel)
    integer,intent(in)::ilevel
    integer,parameter::tile_width=128
    integer::i,ind,cell,k,n,first,last,m,j,status,bad,all_bad,info,nt,ni
    integer::iteration_sum,iteration_max,global_iteration_sum,global_iteration_max
    integer,allocatable,target::cells(:)
    integer,allocatable::slots(:),neighbors(:,:)
    real(dp),allocatable,target::chemical(:,:),gas_total(:),hydrogen(:)
    real(dp),allocatable::dust_energy(:)
    real(dp),allocatable::density(:),primary(:),old_material(:),capacity(:),material(:),temperature(:)
    real(dp),allocatable::gas_energy(:),gas_capacity(:),gas_transfer(:),area(:),weights(:,:),curve(:,:)
    real(dp),allocatable::incoming(:,:,:),ir_trial(:,:,:)
    real(dp)::sl,st,sd,sv,snh,st2,dt,dx,volume,energy_unit,chat,physical_cmb
    real(dp)::row(nvar_all),next_state(chimes_ns),nhydrogen
    real(dp)::zero_n(1,9),zero_e(1,9),next_n(1,9),next_e(1,9),ledger(11),grain_n(9),grain_e(9)
    real(dp)::gas_code,nonthermal,td,gas_delta,local_budget(4),global_budget(4),local_balance,global_balance
    real(dp)::chimes_wall,dust_wall,tile_chimes_start,tile_dust_start,tile_chimes_wall,tile_dust_wall
    real(dp)::global_chimes_wall,global_dust_wall
    real(dp)::direction(3,1),angular_weight(1)
    type(dust_ir_diagnostics)::diag
    type(dust_live_coarse_trial)::coarse
    logical,save::entry_reported=.false.
    logical::staged_cuda,level_queue
#ifdef SNRT_CHIMES_CUDA
    integer,parameter::cohort_width=16
    integer::group,dispatch_status
    integer(c_int),target::cell_status(tile_width)
    real(c_double),target::cell_td(tile_width)
    integer(c_int),allocatable,target::level_status(:)
    type(dark_tile_context),target::context
#endif
    include 'mpif.h'

    if(.not.dust_chimes_enabled().or.trim(radiation_transport)/='none')return
    ! The dark orchestrator currently shares the exact, admitted cosmological
    ! two-size/DL01/D03 closure. Other material layouts retain their explicit
    ! admission failures rather than receiving an unvalidated approximation.
    if(.not.cosmo.or..not.dust_optics_enabled().or..not.dust_material_composition_enabled())return
    if(.not.snrt_dust_contract_loaded.or.snrt_dust_contract_version/=4.or. &
         .not.snrt_dust_contract_exchange_enabled)return
    if(.not.ieee_is_finite(aexp).or.aexp<=0)return
    if(ichem<1.or.ichem+10>nvar)return

    call units(sl,st,sd,sv,snh,st2)
    dt=dtnew(ilevel)*st;dx=(boxlen/2d0**ilevel)*sl;volume=dx**3;energy_unit=sd*sv**2
    physical_cmb=2.727d0/aexp
    chat=.01d0*snrt_c_cgs
    if(any(.not.ieee_is_finite([sl,st,sd,sv,snh,st2,dt,dx,volume,energy_unit,physical_cmb])).or. &
         min(sl,st,sd,sv,snh,st2,dt,dx,volume,energy_unit,physical_cmb)<=0)then
       call MPI_ABORT(MPI_COMM_WORLD,31,info)
       return
    endif
    call chimes_live_initialize(status)
    bad=merge(0,1,status==0)
    call snrt_dust_live_prepare(status)
    if(status/=0)bad=1
    call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(all_bad/=0.or.info/=0)then
       if(myid==1)write(*,'(A,I0)')'ERROR: dark dust/CHIMES material initialization failed, status=',status
       call MPI_ABORT(MPI_COMM_WORLD,32,info)
       return
    endif

    nt=snrt_dust_contract_number_temperature;ni=snrt_dust_contract_number_ir
    n=active(ilevel)%ngrid*twotondim
    allocate(cells(n),chemical(chimes_ns,n),gas_total(n),dust_energy(n),hydrogen(n))
    k=0
    do ind=1,twotondim
       do i=1,active(ilevel)%ngrid
          cell=ICELL_OF(active(ilevel)%igrid(i),ind)
          if(son(cell)/=0)cycle
          k=k+1;cells(k)=cell
       enddo
    enddo
    if(.not.entry_reported)then
       if(myid==1)then
          write(*,'(A,I0,A,I0,A,I0,A,I0)')' SNRT dark material begin level=',ilevel, &
               ' local_leaf_cells=',k,' tile_width=',tile_width,' tile_count=',(k+tile_width-1)/tile_width
          flush(6)
       endif
       entry_reported=.true.
    endif
    zero_n=0;zero_e=0;direction(:,1)=[1d0,0d0,0d0];angular_weight=1d0
    allocate(density(tile_width),primary(tile_width),old_material(tile_width),capacity(tile_width), &
         material(tile_width),temperature(tile_width),gas_energy(tile_width),gas_capacity(tile_width), &
         gas_transfer(tile_width),area(tile_width),weights(4,tile_width),curve(nt,tile_width), &
         incoming(ni,1,tile_width))
    allocate(slots(tile_width),neighbors(6,tile_width));slots=0;neighbors=0;incoming=0
    bad=0;local_budget=0;local_balance=0
    chimes_wall=0d0;dust_wall=0d0;iteration_sum=0;iteration_max=0
    staged_cuda=.false.
    level_queue=.false.
#ifdef SNRT_CHIMES_CUDA
    staged_cuda=snrt_chimes_cold_enabled().and.dark_cuda_requested()==1
    level_queue=snrt_chimes_cold_enabled().and.dark_level_requested()==1
    if(level_queue)then
       allocate(level_status(k));level_status=1
       context%count=n;context%width=k;context%first=1;context%offset=1;context%level_queue=1
       context%sd=sd;context%sv=sv;context%dt=dt;context%dx=dx
       context%cells=c_loc(cells);context%chemical=c_loc(chemical)
       context%gas_total=c_loc(gas_total);context%hydrogen=c_loc(hydrogen)
       context%td=c_null_ptr;context%status=c_loc(level_status)
       tile_chimes_start=omp_get_wtime()
       dispatch_status=dark_level_run(int(k,c_int),c_funloc(dark_tile_cell),c_loc(context))
       chimes_wall=omp_get_wtime()-tile_chimes_start
       if(dispatch_status/=0.or.any(level_status/=0))bad=1
       deallocate(level_status)
       call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
       if(all_bad/=0.or.info/=0)then
          call MPI_ABORT(MPI_COMM_WORLD,33,info)
          return
       endif
    endif
#endif

    do first=1,k,tile_width
       last=min(k,first+tile_width-1);m=last-first+1
       ! CHIMES receives exactly zero transported photons. The current dust
       ! temperature is supplied before the independent CMB/material solve.
       tile_chimes_start=omp_get_wtime()
!$omp parallel do schedule(dynamic,1) private(j,cell,row,next_state,zero_n,zero_e,next_n,next_e, &
!$omp& ledger,grain_n,grain_e,gas_code,nonthermal,td,gas_delta,nhydrogen,status) &
!$omp& reduction(max:bad) reduction(+:local_budget)
       do j=1,m
          cell=cells(first+j-1);row=uold(cell,1:nvar_all)
          zero_n=0;zero_e=0
          if(.not.all(ieee_is_finite(row)).or.row(1)<=0)then
             bad=1;cycle
          endif
          ! The cold adapter performs the authoritative state reconciliation
          ! immediately below; doing it here as well duplicated a full 157-
          ! species pass per cell. This admitted bulk-D03 layout contains no
          ! solid hydrogen, so X_H is the exact gas H mass fraction.
          next_state=row(ichimes:ichimes+chimes_ns-1)
          gas_code=row(ndim+2);next_n=0;next_e=0;ledger=0;grain_n=0;grain_e=0
          density(j)=row(idust)*sd/snrt_dust_contract_mass_per_h_g
          nhydrogen=row(ichem)*sd/atomic_mh
          if(density(j)>0)then
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt), &
                  row(idust_species:idust_species+1),snrt_dust_contract_mass_per_h_g,curve(:,j),status)
             if(status==0)call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:nt), &
                  curve(:,j),row(idust_energy)*sv**2/density(j),td,status)
             if(status/=0)then
                bad=1;cycle
             endif
          else
             td=physical_cmb
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:nt),[0d0,0d0], &
                  snrt_dust_contract_mass_per_h_g,curve(:,j),status)
             if(status/=0)then
                bad=1;cycle
             endif
          endif
          call dust_composition_area(row(idust_bins:idust_bins+3), &
               snrt_dust_contract_mass_per_h_g,area(j),status)
          if(status==0)call d03_cell_weights(row(idust_bins:idust_bins+3),weights(:,j),status)
          if(status/=0)then
             bad=1;cycle
          endif
#ifdef SNRT_CHIMES_CUDA
          if(level_queue)then
             gas_delta=(gas_total(first+j-1)-row(ndim+2))*energy_unit*volume
             local_budget(1)=local_budget(1)+gas_delta
             cycle ! Chemistry already finished from the same read-only uold.
          endif
          if(staged_cuda)then
             cell_td(j)=td
             cycle ! Prepare first; independent chemistry runs in cohorts below.
          endif
#endif
          gas_code=row(ndim+2)
          if(snrt_chimes_cold_enabled())then
             call chimes_live_cold_stage(cell,sd,sv,dt,dx,td,.01d0,1,zero_n,zero_e,next_state,gas_code, &
                  next_n,next_e,ledger,grain_n,grain_e,status,staged_row=row)
          else if(snrt_chimes_band_enabled())then
             call chimes_live_band_stage(cell,sd,sv,dt,dx,td,.01d0,1,zero_n,zero_e,next_state,gas_code, &
                  next_n,next_e,ledger(1:9),status,staged_row=row)
          else
             call chimes_live_stage(cell,sd,sv,dt,dx,td,area(j)*density(j)/max(nhydrogen,tiny(1d0)), &
                  .01d0,zero_n(1,:),next_state,gas_code,next_n(1,:),status,staged_row=row)
          endif
          if(status/=0.or.any(next_n/=0d0).or.any(next_e/=0d0))then
             bad=1;cycle
          endif
          chemical(:,first+j-1)=next_state
          gas_total(first+j-1)=gas_code
          gas_delta=(gas_code-row(ndim+2))*energy_unit*volume
          local_budget(1)=local_budget(1)+gas_delta
          hydrogen(first+j-1)=max(nhydrogen,0d0)
       enddo
!$omp end parallel do
#ifdef SNRT_CHIMES_CUDA
       if(staged_cuda.and..not.level_queue.and.bad==0)then
          cell_status=0
!$omp parallel do schedule(dynamic,1) private(group,context,dispatch_status) reduction(max:bad)
          do group=1,m,cohort_width
             context%count=n;context%width=tile_width;context%first=first;context%offset=group;context%level_queue=0
             context%sd=sd;context%sv=sv;context%dt=dt;context%dx=dx
             context%cells=c_loc(cells);context%chemical=c_loc(chemical)
             context%gas_total=c_loc(gas_total);context%hydrogen=c_loc(hydrogen)
             context%td=c_loc(cell_td);context%status=c_loc(cell_status)
             dispatch_status=dark_tile_run(int(min(cohort_width,m-group+1),c_int), &
                  c_funloc(dark_tile_cell),c_loc(context))
             if(dispatch_status/=0)bad=1
          enddo
!$omp end parallel do
          if(any(cell_status(1:m)/=0))bad=1
          if(bad==0)then
             do j=1,m
                cell=cells(first+j-1)
                local_budget(1)=local_budget(1)+(gas_total(first+j-1)-uold(cell,ndim+2))*energy_unit*volume
             enddo
          endif
       endif
#endif
       tile_chimes_wall=omp_get_wtime()-tile_chimes_start
       chimes_wall=chimes_wall+tile_chimes_wall
       if(bad/=0)cycle

       do j=1,m
          cell=cells(first+j-1);row=uold(cell,1:nvar_all)
          nonthermal=.5d0*sum(row(2:ndim+1)**2)/row(1)+magnetic_energy(row)
#if NENER>0
          nonthermal=nonthermal+sum(row(inener:inener+NENER-1))
#endif
          gas_energy(j)=(gas_total(first+j-1)-nonthermal)*energy_unit
          gas_capacity(j)=chimes_live_capacity(chemical(:,first+j-1),sd)
          old_material(j)=row(idust_energy)*energy_unit
          capacity(j)=1d0;primary(j)=0d0
          if(gas_energy(j)<=0.or.gas_capacity(j)<=0.or..not.ieee_is_finite(gas_energy(j)))bad=1
       enddo
       if(bad/=0)cycle
       tile_dust_start=omp_get_wtime()
       call snrt_dust_live_stage(ilevel,cells(first:last),slots(1:m),neighbors(:,1:m),direction,angular_weight, &
            dx,dt,chat,density(1:m),primary(1:m),old_material(1:m),capacity(1:m),ir_trial,material(1:m), &
            temperature(1:m),diag,status,coarse,gas_energy=gas_energy(1:m),gas_capacity=gas_capacity(1:m), &
            n_hydrogen=hydrogen(first:last),gas_transfer=gas_transfer(1:m),cell_material_u=curve(:,1:m), &
            cell_collision_area=area(1:m),cell_weights=weights(:,1:m),incoming_radiation=incoming(:,:,1:m), &
            thin_local_escape=.true.)
       tile_dust_wall=omp_get_wtime()-tile_dust_start
       dust_wall=dust_wall+tile_dust_wall
       iteration_sum=iteration_sum+diag%iterations
       iteration_max=max(iteration_max,diag%iterations)
       if(first==1.and.myid==1)then
          write(*,'(A,I0,A,I0,A,ES12.4,A,ES12.4,A,I0,A,I0)') &
               ' SNRT dark material first_tile level=',ilevel,' cells=',m,' CHIMES_s=',tile_chimes_wall, &
               ' dust_s=',tile_dust_wall,' dust_iterations=',diag%iterations,' status=',status
          flush(6)
       endif
       if(status/=0)then
          bad=1;cycle
       endif
       do j=1,m
          if(.not.ieee_is_finite(material(j)).or..not.ieee_is_finite(gas_transfer(j)))then
             bad=1;cycle
          endif
          dust_energy(first+j-1)=material(j)/energy_unit
          gas_total(first+j-1)=gas_total(first+j-1)-gas_transfer(j)/energy_unit
          cell=cells(first+j-1);row=uold(cell,1:nvar_all)
          nonthermal=.5d0*sum(row(2:ndim+1)**2)/row(1)+magnetic_energy(row)
#if NENER>0
          nonthermal=nonthermal+sum(row(inener:inener+NENER-1))
#endif
          if((gas_total(first+j-1)-nonthermal)<=0d0)bad=1
       enddo
       local_budget(2)=local_budget(2)+diag%escaped_erg
       local_budget(3)=local_budget(3)+diag%background_erg
       local_budget(4)=local_budget(4)+sum(gas_transfer(1:m))*volume
       local_balance=max(local_balance,diag%balance_relative,diag%local_relative)
       if(allocated(ir_trial))deallocate(ir_trial)
    enddo

    call MPI_ALLREDUCE(bad,all_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(all_bad/=0.or.info/=0)then
       if(myid==1)write(*,'(A,I0)')'ERROR: dark dust/CHIMES level transaction rejected, level=',ilevel
       call MPI_ABORT(MPI_COMM_WORLD,33,info)
       return
    endif
    do i=1,k
       cell=cells(i)
       uold(cell,ichimes:ichimes+chimes_ns-1)=chemical(:,i)
       uold(cell,idust_energy)=dust_energy(i)
       uold(cell,ndim+2)=gas_total(i)
    enddo
    call upload_fine(ilevel)
    do i=ichimes,ichimes+chimes_ns-1
       call make_virtual_fine_dp(uold(1,i),ilevel)
    enddo
    call make_virtual_fine_dp(uold(1,idust_energy),ilevel)
    call make_virtual_fine_dp(uold(1,ndim+2),ilevel)
    call MPI_ALLREDUCE(local_budget,global_budget,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(local_balance,global_balance,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(chimes_wall,global_chimes_wall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(dust_wall,global_dust_wall,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(iteration_sum,global_iteration_sum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(iteration_max,global_iteration_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,j)
    if(myid==1)write(*,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES10.3,A,ES12.4,A,ES12.4,A,I0,A,I0)') &
         ' SNRT dark material level=',ilevel,' CHIMES_dE_erg=',global_budget(1), &
         ' IR_escaped_erg=',global_budget(2),' CMB_receipt_erg=',global_budget(3), &
         ' gas_dust_transfer_erg=',global_budget(4),' max_balance=',global_balance, &
         ' CHIMES_wall_s=',global_chimes_wall,' dust_wall_s=',global_dust_wall, &
         ' dust_iteration_sum=',global_iteration_sum,' dust_iteration_max_tile=',global_iteration_max
    deallocate(cells,chemical,gas_total,dust_energy,hydrogen,density,primary,old_material,capacity,material,temperature, &
         gas_energy,gas_capacity,gas_transfer,area,weights,curve,incoming,slots,neighbors)
  end subroutine snrt_dark_matter_advance_level
end module snrt_dark_matter_runtime

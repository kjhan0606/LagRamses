! Live compact M_N operator. All source/radiation/material changes remain
! private until every MPI rank accepts the complete level transaction.
module snrt_moment_ramses
  use amr_commons, only: dp,ndim,ncpu,myid,dtnew,texp,aexp,active,cosmo
  use hydro_commons, only: uold,magnetic_energy
  use hydro_parameters, only: gamma,nvar,nvar_all,inener,idust,idust_energy,idust_species,idust_bins,ichem,dust_relative_motion, &
       idust_iron,idust_pah,idust_momentum,ndust_phase,ichimes
  use pm_commons, only: nsink,xsink,idsink,agn_pending_erg,headp,numbp,nextp,ptypep,PTYPE_STAR,xp,mp0,tpp,zp
  use snrt_state, only: snrt_ngroups,snrt_hydrogen_ii,snrt_helium_ii,snrt_helium_iii,snrt_neutral_fraction
  use snrt_moment_live
  use snrt_moment_transport, only: mn_ok,mn_reconstruct,mn_project
  use snrt_moment_spectrum, only: mn_band_split_angular,mn_band_join
  use snrt_moment_amr
  use snrt_moment_dispatch
  use snrt_moment_material
  use snrt_spectral_contract, only: snrt_band_enabled,snrt_group_energy_fraction,snrt_group_mean_energy_ev, &
       snrt_group_cross_section_cm2,snrt_group_cross_section_hei_cm2,snrt_group_cross_section_heii_cm2, &
       snrt_group_edges_ev,snrt_d03_band_enabled,snrt_fe_band_enabled,snrt_node_secondaries_enabled, &
       snrt_chimes_band_enabled,snrt_chimes_cold_enabled
  use snrt_agn_source, only: snrt_c_cgs,snrt_ev_to_erg,snrt_agn_photon_budget_energy,snrt_agn_source_commit
  use snrt_agn_locator, only: snrt_agn_find_local_leaf
  use snrt_stellar_source, only: stellar_sed_enabled,stellar_photon_interval,stellar_sed_has_energy
  use snrt_rt_transaction, only: snrt_rt_iteration_config
  use snrt_thermochemistry
  use snrt_atomic_cooling, only: atomic_mh,atomic_temperature,atomic_heat_capacity,atomic_advance
  use dust_mass_physics, only: dust_atomic_cooling_enabled,dust_chimes_enabled,dust_gas_elements, &
       dust_optics_enabled,dust_sublimation_rt_enabled,dust_iron_enabled,dust_pah_enabled, &
       dust_pah_nstate,dust_pah_state_mass,dust_pah_inventory,dust_pah_charged,dust_pah_hydrogenated, &
       dust_pah_h2_enabled,dust_pah_atomization,dust_fe_uv_enabled,dust_fe_primary_limit
  use snrt_dust_contract
#ifdef DUST_LIVE
  use snrt_dust_live, only: snrt_dust_live_moment_cell,snrt_dust_live_moment_tile,snrt_dust_live_prepare
  use snrt_dust_ir, only: dust_ir_diagnostics,snrt_dust_material_temperature
  use dust_composition_material, only: dust_material_composition_enabled,dust_composition_curve,dust_composition_area
  use dust_composition_optics, only: d03_cell_weights,d03_opacity_basis,d03_ng,d03_nir, &
       d03_band_abs,d03_band_transport,d03_band_ev
  use dust_iron_optics, only: fe_six_band_abs,fe_six_band_transport,fe_band_ev
  use dust_iron_compare, only: iron_compare_curve,iron_compare_weights,iron_compare_temperature, &
       fe_six_opacity_basis,iron_compare_neutral_area,iron_compare_source_receipt
  use dust_pah_live_model, only: pah_live_prepare,pah_primary_alpha,pah_primary_sigma,pah_primary_supported
  use dust_phase_state, only: dust_phase_read,dust_phase_kinetic,dust_phase_erode
  use snrt_moving_scatter, only: snrt_moving_scatter_cell
#endif
#ifdef SNRT_CHIMES
  use snrt_chimes_runtime, only: chimes_cell_state,chimes_live_capacity,chimes_live_stage,chimes_live_band_stage, &
       chimes_live_cold_stage,chimes_live_grain_scatter,chimes_live_fe_uv_stage
  use snrt_chimes, only: chimes_ns,chimes_boltzmann,chimes_group_binding
#endif
  use iso_c_binding, only: c_int,c_double,c_funptr,c_null_funptr,c_funloc
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use omp_lib, only: omp_get_wtime
#ifndef WITHOUTMPI
  use mpi_mod
#endif
  implicit none
  private
  public :: mn_ramses_advance
  interface
     function band_material(nd,ng,nb,scatter,n,e,edges,columns,budget,an,ae,returned,secondary,xi,deposition, &
          grains,kabs,ksca,nodes,weights) bind(C,name='snrt_mn_band_material_c') result(ierr)
       import c_int,c_double,c_funptr
       integer(c_int),value::nd,ng,nb,scatter
       type(c_funptr),value::secondary
       real(c_double),value::xi
       real(c_double)::n(*),e(*),budget(*),an(*),ae(*),returned(*),deposition(*)
       real(c_double),intent(in)::edges(*),columns(*),grains(*),kabs(*),ksca(*),nodes(*),weights(*)
       integer(c_int)::ierr
     end function
  end interface
#ifdef DUST_LIVE
  type :: material_stage_packet
     ! Pre-CHIMES work is kept in the packet while a tile is assembled.
     ! The packet is also the per-cell write set for the batched IR stage:
     ! no caller-owned material/gas/phase array is exposed until the tile
     ! transaction and its M_N closure checks have succeeded.
     logical :: valid=.false.,has_exchange=.false.,has_curve=.false.,has_weights=.false.
     logical :: has_sublimation=.false.,has_pah=.false.,has_phase=.false.,has_chimes=.false.
     integer :: index=0,cell=0,slot=0
     real(dp) :: rho=0d0,density=0d0,primary=0d0,material=0d0,capacity=0d0
     real(dp) :: material_next=0d0,temperature=0d0,gas_next=0d0,nonthermal=0d0
     real(dp) :: heat=0d0,primary_work=0d0,projection_n=0d0,projection_e=0d0
     real(dp) :: electron_cv=0d0,phase_ke=0d0
     real(dp),allocatable :: nn(:,:),ee(:,:),an(:,:),ae(:,:),old_row(:),next_fraction(:),gasx(:),chemical(:)
     real(dp),allocatable :: gas_energy(:),gas_capacity(:),n_hydrogen(:),transfer(:)
     real(dp),allocatable :: curve(:,:),area(:),weights(:,:)
     real(dp),allocatable :: sub_bins(:,:),sub_next(:,:),pah_number(:,:),pah_spectrum(:,:)
     real(dp),allocatable :: pah_heat(:,:),pah_captures(:,:)
     real(dp),allocatable :: phase_mass(:,:),phase_old(:,:,:),phase_next(:,:,:),phase_work(:,:)
     real(dp),allocatable :: electrons(:),hydrogen(:),molecular_h2(:),carbon(:),carbon_ion(:)
     real(dp),allocatable :: gas_momentum(:),irproj(:,:)
  end type
#endif
contains
  subroutine mn_ramses_advance(lev,chat_factor,config,ierr,step_start)
    integer,intent(in) :: lev
    real(dp),intent(in) :: chat_factor
    type(snrt_rt_iteration_config),intent(in) :: config
    integer,intent(out) :: ierr
    real(dp),optional,intent(in) :: step_start
    type(mn_amr_mesh) :: mesh
    real(dp),allocatable :: number(:,:,:),energy(:,:,:),infrared(:,:,:),oldn(:,:,:),olde(:,:,:),oldir(:,:,:)
    real(dp),allocatable :: rows(:,:),fraction(:,:),source_n(:,:),source_e(:,:)
    real(dp) :: sl,st,sd,sv,snh,st2,dt,subdt,chat,energy_unit,volume,clock0
    real(dp) :: photons(snrt_ngroups),source_energy(snrt_ngroups),lum,age_scale
    real(dp) :: escape(36),projection(36),receipt(12),global_receipt(12),initial_n,final_n,initial_e,final_e
    real(dp) :: pair_escape,pair_projection,lo,hi,source_total_n,source_total_e
    real(dp) :: state_range(4),global_range(4),balance_scale(2),global_scale(2)
    integer :: nm,ng,ni,i,j,g,s,k,nsub,cell,found,grid,part,ip,ig,status,info
    integer :: tile_index,tile_first,tile_last,ntile
    ! Keep the RAMSES-side transaction large enough to amortize the native
    ! material/IR dispatch.  The backend still applies its own memory-aware
    ! GPU batch policy; this is only the caller-owned write-set bound.
    integer,parameter :: material_tile_width=256
    integer,allocatable :: accepted(:),owners(:)
    logical :: atomic_on,band_on,chimes_on
    logical,save :: host_reported=.false.

    ierr=0;clock0=omp_get_wtime();receipt=0
    call mn_dispatch_initialize(ierr)
    if(.not.host_reported)then
       if(myid==1)write(*,*)'SNRT M_N live AMR transport: shared-stream CUDA/OpenMP, MPI, FP64, 16x24 quadrature'
       host_reported=.true.
    endif
    band_on=snrt_band_enabled();chimes_on=dust_chimes_enabled()
    if((band_on.and.chimes_on.and..not.snrt_chimes_band_enabled()).or. &
         (snrt_chimes_band_enabled().and..not.chimes_on))ierr=1
#ifdef SNRT_CHIMES
    if(chimes_on.and..not.snrt_chimes_band_enabled())then
       status=chimes_group_binding(snrt_ngroups,snrt_group_edges_ev,snrt_group_mean_energy_ev)
       ierr=max(ierr,status)
    endif
#endif
    if(config%failure_stage/=0)then
       if(myid==1)write(*,*)'M_N diagnostic fault requested: reject before any source or material publication'
       ierr=1
    endif
    call units(sl,st,sd,sv,snh,st2)
    dt=dtnew(lev)*st;chat=snrt_c_cgs*chat_factor;energy_unit=sd*sv**2
    if(any(.not.ieee_is_finite([sl,st,sd,sv,snh,st2,dt,chat])).or. &
         min(sl,st,sd,sv,snh,st2,dt,chat)<=0)ierr=1
    call mn_collective_status(ierr)
    if(ierr/=0)return
    call mn_amr_prepare(lev,mesh,ierr)
    if(ierr/=0)then
       write(*,*)'M_N topology rejected rank/level/code=',myid,lev,ierr
       return
    endif
    nm=mn_live_basis%nm;ng=snrt_ngroups;ni=0
#ifdef DUST_LIVE
    ni=snrt_dust_contract_number_ir
    if(snrt_dust_contract_version<3.or.ni<1)ierr=1
#endif
    if(nm<9.or.nm>36)ierr=1
    call mn_collective_status(ierr)
    if(ierr/=0)then
       write(*,*)'M_N payload configuration rejected nm/ni/dustversion=',nm,ni,snrt_dust_contract_version
       return
    endif
#ifdef DUST_LIVE
    ! Prepare shared dust/PAH tables and slot storage before the material
    ! loop becomes OpenMP-parallel.  Lazy preparation remains in the stage
    ! for non-M_N callers, but must never be entered concurrently here.
    call snrt_dust_live_prepare(status)
    ierr=max(ierr,status)
    call mn_collective_status(ierr)
    if(ierr/=0)return
#endif
    allocate(number(nm,ng,mesh%nowned),energy(nm,ng,mesh%nowned),infrared(nm,ni,mesh%nowned))
    allocate(rows(nvar_all,mesh%nleaf),fraction(3,mesh%nleaf),source_n(ng,mesh%nleaf),source_e(ng,mesh%nleaf))
    number=0;energy=0;infrared=0;source_n=0;source_e=0
    do i=1,mesh%nowned
       call mn_live_read(mesh%slots(i),number(:,:,i),energy(:,:,i),status)
       ierr=max(ierr,status)
       if(ni>0)then
          call mn_live_ir_read(mesh%slots(i),infrared(:,:,i),status)
          ierr=max(ierr,status)
       endif
    enddo
    do i=1,mesh%nleaf
       rows(:,i)=uold(mesh%cells(i),:)
       s=mesh%slots(i)
       fraction(:,i)=[real(snrt_hydrogen_ii(s),dp),real(snrt_helium_ii(s),dp),real(snrt_helium_iii(s),dp)]
       if(any(.not.ieee_is_finite(rows(:,i))).or.rows(1,i)<=0.or.any(fraction(:,i)<0).or. &
            fraction(1,i)>1.or.sum(fraction(2:3,i))>1)ierr=1
    enddo
    oldn=number;olde=energy;oldir=infrared
    allocate(accepted(max(1,nsink)),owners(max(1,nsink)));accepted=0
    if(nsink>0)then
       if(.not.allocated(agn_pending_erg).or..not.allocated(xsink).or..not.allocated(idsink))then
          ierr=1
       else if(size(agn_pending_erg)<nsink.or.size(xsink,1)<nsink.or.size(idsink)<nsink)then
          ierr=1
       else
          do j=1,nsink
             if(idsink(j)<=0.or.any(idsink(1:j-1)==idsink(j)))ierr=1
             if(.not.ieee_is_finite(agn_pending_erg(j)).or.agn_pending_erg(j)<0)ierr=1
          enddo
       endif
    endif
    call mn_collective_status(ierr)
    if(ierr/=0)then
       write(*,*)'M_N initial payload/material/source rejected rank=',myid
       return
    endif
    volume=mesh%dx**3*sl**3*snh
    do j=1,nsink
       if(agn_pending_erg(j)==0)cycle
       call snrt_agn_find_local_leaf(xsink(j,1:ndim),cell,found)
       if(cell==0.or.found/=lev)cycle
       i=local_leaf(cell)
       if(i==0)then
          ierr=1;cycle
       endif
       do g=1,ng
          call snrt_agn_photon_budget_energy(agn_pending_erg(j),dt,snrt_group_energy_fraction(g), &
               snrt_group_mean_energy_ev(g),lum,photons(g),status)
          ierr=max(ierr,status)
       enddo
       source_n(:,i)=source_n(:,i)+photons/volume
       source_e(:,i)=source_e(:,i)+photons*snrt_group_mean_energy_ev/volume
       accepted(j)=1
    enddo
    owners=accepted
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(accepted,owners,size(accepted),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    if(info/=0)ierr=1
#endif
    if(any(owners>1))ierr=1
    if(stellar_sed_enabled)then
       if(.not.present(step_start))then
          ierr=1
       else
          age_scale=st/aexp**2/31557600d6
          do ig=1,active(lev)%ngrid
             grid=active(lev)%igrid(ig);part=headp(grid)
             do ip=1,numbp(grid)
                if(ptypep(part)==PTYPE_STAR)then
                   if(stellar_sed_has_energy)then
                   call stellar_photon_interval((step_start-tpp(part))*age_scale,(texp-tpp(part))*age_scale, &
                        zp(part),mp0(part)*sd*sl**3/1.98847d33,photons,status,energy_ev=source_energy)
                   else
                   call stellar_photon_interval((step_start-tpp(part))*age_scale,(texp-tpp(part))*age_scale, &
                        zp(part),mp0(part)*sd*sl**3/1.98847d33,photons,status)
                   source_energy=photons*snrt_group_mean_energy_ev
                   endif
                   ierr=max(ierr,status)
                   if(status/=0)write(*,*)'M_N stellar source rejected rank/particle/status=',myid,part,status
                   if(status==0.and.any(photons>0))then
                      call snrt_agn_find_local_leaf(xp(part,1:ndim),cell,found,grid,lev)
                      i=local_leaf(cell)
                      if(i==0.or.found/=lev)then
                         ierr=1
                      else
                         source_n(:,i)=source_n(:,i)+photons/volume
                         source_e(:,i)=source_e(:,i)+source_energy/volume
                      endif
                   endif
                endif
                part=nextp(part)
             enddo
          enddo
       endif
    endif
    call mn_collective_status(ierr)
    if(ierr/=0)return
    source_total_n=sum(source_n)*mesh%dx**3;source_total_e=sum(source_e)*mesh%dx**3
    if(band_on)then
       do i=1,mesh%nleaf
          do g=1,ng
             lo=snrt_group_edges_ev(g);hi=snrt_group_edges_ev(g+1)
             photons(g)=(source_e(g,i)-lo*source_n(g,i))/(hi-lo)
             source_e(g,i)=(hi*source_n(g,i)-source_e(g,i))/(hi-lo)
             source_n(g,i)=photons(g)
          enddo
       enddo
       if(any(source_n<0).or.any(source_e<0))ierr=1
    endif
    call mn_collective_status(ierr)
    if(ierr/=0)return
    initial_n=0;initial_e=0
    do i=1,mesh%nowned
       call inventory(number(:,:,i),energy(:,:,i),lum,volume)
       initial_n=initial_n+lum*mesh%volume(i)
       initial_e=initial_e+volume*mesh%volume(i)
    enddo
    nsub=max(1,ceiling(12*chat*dt/(mesh%dx*sl)));subdt=dt/nsub
    atomic_on=dust_atomic_cooling_enabled()
    ! Source interval is apportioned without consuming replicated AGN fuel.
    ! MC15/SSPRK2 transport plus explicit source/material operator splitting.
    do k=1,nsub
       do i=1,mesh%nleaf
          number(1,:,i)=number(1,:,i)+source_n(:,i)/nsub
          energy(1,:,i)=energy(1,:,i)+source_e(:,i)/nsub
       enddo
       do g=1,ng
          call mn_amr_advance(mn_live_basis,mesh,number(:,g,:),chat*subdt/sl, &
               escape(1:nm),projection(1:nm),ierr)
          if(ierr/=0)then
             if(myid==1)write(*,*)'M_N number transport rejected level/group/substep=',lev,g,k
             return
          endif
          pair_escape=escape(1);pair_projection=projection(1)
          call mn_amr_advance(mn_live_basis,mesh,energy(:,g,:),chat*subdt/sl, &
               escape(1:nm),projection(1:nm),ierr)
          if(ierr/=0)then
             if(myid==1)write(*,*)'M_N energy transport rejected level/group/substep=',lev,g,k
             return
          endif
          if(band_on)then
             lo=snrt_group_edges_ev(g);hi=snrt_group_edges_ev(g+1)
             receipt(3)=receipt(3)+pair_escape+escape(1)
             receipt(4)=receipt(4)+pair_projection+projection(1)
             receipt(10)=receipt(10)+hi*pair_escape+lo*escape(1)
             receipt(11)=receipt(11)+hi*pair_projection+lo*projection(1)
          else
             receipt(3)=receipt(3)+pair_escape;receipt(4)=receipt(4)+pair_projection
             receipt(10)=receipt(10)+escape(1);receipt(11)=receipt(11)+projection(1)
          endif
       enddo
       do g=1,ni
          call mn_amr_advance(mn_live_basis,mesh,infrared(:,g,:),chat*subdt/sl, &
               escape(1:nm),projection(1:nm),ierr)
          if(ierr/=0)return
       enddo
#ifdef DUST_LIVE
       if(snrt_dust_contract_exchange_enabled)then
          ! CHIMES preparation remains per cell, but the expensive material/
          ! IR receiver is committed once per bounded tile.  The tile width
          ! is deliberately fixed here: it is a work-unit bound, not a GPU
          ! stream tuning parameter.
          ntile=(mesh%nleaf+material_tile_width-1)/material_tile_width
!$omp parallel do schedule(dynamic,1) private(tile_index,tile_first,tile_last,status) reduction(max:ierr)
          do tile_index=1,ntile
             tile_first=(tile_index-1)*material_tile_width+1
             tile_last=min(mesh%nleaf,tile_first+material_tile_width-1)
             call material_tile(tile_first,tile_last,subdt,status)
             if(status/=0)then
                write(*,*)'M_N material tile failure rank/level/tile/substep/code=',myid,lev,tile_index,k,status
                ierr=max(ierr,status)
             endif
          enddo
!$omp end parallel do
       else
#endif
!$omp parallel do schedule(dynamic,1) private(i,status) reduction(max:ierr)
          do i=1,mesh%nleaf
             call material_cell(i,subdt,status)
             if(status/=0)then
                write(*,*)'M_N material failure rank/level/cell/substep/code=',myid,lev,mesh%cells(i),k,status
                ierr=max(ierr,status)
             endif
          enddo
!$omp end parallel do
#ifdef DUST_LIVE
       endif
#endif
       call mn_collective_status(ierr)
       if(ierr/=0)return
    enddo
    ! Global number AND actual radiation energy accounting, before any
    ! material publication. Numerical closure receipts are not matter heat.
    final_n=0;final_e=0
    do i=1,mesh%nowned
       call inventory(number(:,:,i),energy(:,:,i),lum,volume)
       final_n=final_n+lum*mesh%volume(i)
       final_e=final_e+volume*mesh%volume(i)
    enddo
    receipt(1)=source_total_n
    receipt(5)=final_n-initial_n-receipt(1)+receipt(2)+receipt(3)-receipt(4)
    receipt(8)=final_n
    receipt(12)=final_e-initial_e-source_total_e+receipt(9)+receipt(10)-receipt(11)
    global_receipt=receipt
    balance_scale=[initial_n+receipt(1)+final_n,initial_e+source_total_e+final_e]
    global_scale=balance_scale
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(receipt,global_receipt,12,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if(info/=0)ierr=1
    call MPI_ALLREDUCE(balance_scale,global_scale,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    if(info/=0)ierr=1
#endif
    if(any(.not.ieee_is_finite(global_receipt)).or. &
         abs(global_receipt(5))>1d-9*max(global_scale(1),tiny(1d0)).or. &
         abs(global_receipt(12))>1d-9*max(global_scale(2),tiny(1d0)))then
       if(myid==1)write(*,'(A,4ES22.13)')'M_N accounting rejected N/E residual/scales=', &
            global_receipt(5),global_receipt(12),global_scale
       ierr=1
    endif
    call mn_collective_status(ierr)
    if(ierr/=0)return
    ! Allocation/realizability is checked collectively BEFORE publication.
    do i=1,mesh%nowned
       call mn_live_write(mesh%slots(i),number(:,:,i),energy(:,:,i),status,reserve_only=.true.)
       ierr=max(ierr,status)
       if(ni>0)then
          call mn_live_ir_write(mesh%slots(i),infrared(:,:,i),status,reserve_only=.true.)
          ierr=max(ierr,status)
       endif
    enddo
    call mn_collective_status(ierr)
    if(ierr/=0)return
    do i=1,mesh%nowned
       call mn_live_write(mesh%slots(i),number(:,:,i),energy(:,:,i),status)
       ierr=max(ierr,status)
       if(ni>0)then
          call mn_live_ir_write(mesh%slots(i),infrared(:,:,i),status)
          ierr=max(ierr,status)
       endif
    enddo
    call mn_collective_status(ierr)
    if(ierr/=0)then
       do i=1,mesh%nowned
          call mn_live_write(mesh%slots(i),oldn(:,:,i),olde(:,:,i),status)
          if(ni>0)call mn_live_ir_write(mesh%slots(i),oldir(:,:,i),status)
       enddo
       return
    endif
    do i=1,mesh%nleaf
       uold(mesh%cells(i),:)=rows(:,i);s=mesh%slots(i)
       snrt_hydrogen_ii(s)=fraction(1,i);snrt_helium_ii(s)=fraction(2,i);snrt_helium_iii(s)=fraction(3,i)
       snrt_neutral_fraction(s)=1-fraction(1,i)
    enddo
    do j=1,nsink
       if(owners(j)==1)call snrt_agn_source_commit(agn_pending_erg(j),.true.)
    enddo
    state_range=0
    if(mesh%nleaf>0)then
       state_range(1)=maxval(fraction(1,:));state_range(2)=maxval(rows(5,:))
#ifdef DUST_LIVE
       state_range(3)=maxval(rows(idust_energy,:))
       state_range(4)=maxval(sum(infrared(1,:,:),dim=1))
#endif
    endif
    global_range=state_range
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(state_range,global_range,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)then
       call MPI_ABORT(MPI_COMM_WORLD,info,status)
    endif
#endif
    if(myid==1)then
       write(*,'(A,3I8,A,F12.3)')' SNRT_MN_COMMIT level/order/substeps=',lev,mn_live_basis%order,nsub, &
            ' wall=',omp_get_wtime()-clock0
       write(*,'(A,8ES22.13)')' SNRT_MN_LEDGER source/absorbed/escaped/projection/residual/gas_heat/dust_heat/N=', &
            global_receipt(1:8)
       write(*,'(A,4ES22.13)')' SNRT_MN_ENERGY absorbed/escaped/projection/residual=',global_receipt(9:12)
       write(*,'(A,4ES22.13)')' SNRT_MN_STATE max_xHII/gasE_code/dustE_code/IR_erg_cm3=',global_range
       write(*,'(A,4I12)')' SNRT_MN_DISPATCH rank1 closureCPU/GPU faceCPU/GPU=', &
            mn_cpu_closures,mn_gpu_closures,mn_cpu_faces,mn_gpu_faces
    endif
  contains
    subroutine inventory(n,e,total_n,total_e)
      real(dp),intent(in)::n(:,:),e(:,:)
      real(dp),intent(out)::total_n,total_e
      if(band_on)then
         total_n=sum(n(1,:)+e(1,:))
         total_e=sum(snrt_group_edges_ev(2:ng+1)*n(1,:)+snrt_group_edges_ev(1:ng)*e(1,:))
      else
         total_n=sum(n(1,:));total_e=sum(e(1,:))
      endif
    end subroutine
    integer function local_leaf(c) result(index)
      integer,intent(in) :: c
      integer :: q
      index=0
      do q=1,mesh%nleaf
         if(mesh%cells(q)==c)then
            index=q;return
         endif
      enddo
    end function

#ifdef DUST_LIVE
    subroutine material_cell(index,cell_dt,status,packet,defer_ir)
#else
    subroutine material_cell(index,cell_dt,status)
#endif
      integer,intent(in) :: index
      integer,intent(out) :: status
      real(dp),intent(in) :: cell_dt
#ifdef DUST_LIVE
      type(material_stage_packet),optional,intent(inout) :: packet
      logical,optional,intent(in) :: defer_ir
#endif
      real(dp) :: nn(nm,ng),ee(nm,ng),start(3),guess(3),next_fraction(3),gasx(11)
      real(dp) :: tau(3,ng),target_tau(3,ng),dust_tau(ng),available(3),an(4,ng),ae(4,ng),af(3,4,ng)
      real(dp) :: returned(ng),unassigned(ng),excess(3,ng),threshold(3),hcode,hecode,rho,internal,temp
      real(dp) :: nonthermal,heat,gas_next,loss,density,primary,capacity,material,temperature,irproj(nm,ni)
      real(dp) :: gas_energy(1),gas_capacity(1),transfer(1),opacity(ng),scatter(ng),factor,mu,material_next
      real(dp) :: rn(mn_live_basis%nq,ng),re(mn_live_basis%nq,ng),pn(mn_live_basis%nq,ng),pe(mn_live_basis%nq,ng)
      real(dp) :: ia_ray(mn_live_basis%nq),ib_ray(mn_live_basis%nq),in_n,in_e,out_n,out_e,projection_n,projection_e
      real(dp) :: columns(3),deposition(8),solid_fe,solid_pah(2),old_temperature,chemical_area,delta_heat
      real(dp) :: primary_spectrum(ng,1),old_row(nvar_all),total_chemical_energy,gas_mass,gas_momentum(3),phase_ke,primary_work
      type(c_funptr)::secondary
#ifdef DUST_LIVE
      real(dp),allocatable :: curve(:,:),area(:),weights(:,:)
      real(dp),allocatable :: sub_bins(:,:),sub_next(:,:),pah_number(:,:),pah_heat(:,:),pah_captures(:,:)
      real(dp),allocatable :: phase_mass(:,:),phase_old(:,:,:),phase_next(:,:,:),phase_work(:,:)
      real(dp),allocatable :: electrons(:),hydrogen(:),molecular_h2(:),carbon(:),carbon_ion(:)
      real(dp),allocatable :: pah_spectrum(:,:)
      real(dp) :: electron_cv,kick(3),velocity(3),impulse(3,ng),phase_gain,phase_fraction(ng,ndust_phase)
      real(dp) :: phase_sigma(ng,ndust_phase),phase_tau(ndust_phase,ng),work(ndust_phase),channel_heat(ng)
      real(dp) :: pa(d03_ng,6),ps(d03_ng,6),pg(d03_ng,6),ia(d03_nir,6),isc(d03_nir,6),isg(d03_nir,6)
      real(dp) :: grain_columns(6),dummy(1),ph_e(9,4),ph_p(3,4),events(2)
      type(dust_ir_diagnostics) :: ir_diag
#endif
#ifdef SNRT_CHIMES
      real(dp)::chemical(chimes_ns),elements(11),next_photons(9),old_photons(9),photo_ledger(11),gn(9),ge(9),fe_ledger(8)
#endif
      type(snrt_thermochemistry_result) :: chemistry
      integer :: iteration,q,a,b,nb,jj
      logical :: converged,restoring
#ifdef DUST_LIVE
      restoring=.false.
      gasx=0
      if(present(packet).and.present(defer_ir).and..not.defer_ir)then
         if(.not.packet%valid.or.packet%index/=index)then
            status=1;return
         endif
         restoring=.true.
         rho=packet%rho;density=packet%density;primary=packet%primary
         material=packet%material;capacity=packet%capacity
         material_next=packet%material_next;temperature=packet%temperature
         gas_next=packet%gas_next;nonthermal=packet%nonthermal;heat=packet%heat
         primary_work=packet%primary_work;projection_n=packet%projection_n
         projection_e=packet%projection_e;phase_ke=packet%phase_ke
         old_row=packet%old_row;next_fraction=packet%next_fraction
         nn=packet%nn;ee=packet%ee;an=packet%an;ae=packet%ae;gasx=packet%gasx
         gas_energy=packet%gas_energy;gas_capacity=packet%gas_capacity
         transfer=packet%transfer
#ifdef SNRT_CHIMES
         if(packet%has_chimes)chemical=packet%chemical
#endif
         if(allocated(packet%sub_bins))then
            allocate(sub_bins(size(packet%sub_bins,1),size(packet%sub_bins,2)));sub_bins=packet%sub_bins
         endif
         if(allocated(packet%sub_next))then
            allocate(sub_next(size(packet%sub_next,1),size(packet%sub_next,2)));sub_next=packet%sub_next
         endif
         if(allocated(packet%pah_number))then
            allocate(pah_number(size(packet%pah_number,1),size(packet%pah_number,2)));pah_number=packet%pah_number
         endif
         if(allocated(packet%phase_mass))then
            allocate(phase_mass(size(packet%phase_mass,1),size(packet%phase_mass,2)));phase_mass=packet%phase_mass
         endif
         if(allocated(packet%phase_old))then
            allocate(phase_old(size(packet%phase_old,1),size(packet%phase_old,2),size(packet%phase_old,3)));phase_old=packet%phase_old
         endif
         if(allocated(packet%phase_next))then
            allocate(phase_next(size(packet%phase_next,1),size(packet%phase_next,2),size(packet%phase_next,3)));phase_next=packet%phase_next
         endif
         if(allocated(packet%electrons))then
            allocate(electrons(size(packet%electrons)));electrons=packet%electrons
         endif
         if(allocated(packet%hydrogen))then
            allocate(hydrogen(size(packet%hydrogen)));hydrogen=packet%hydrogen
         endif
         if(allocated(packet%molecular_h2))then
            allocate(molecular_h2(size(packet%molecular_h2)));molecular_h2=packet%molecular_h2
         endif
         if(allocated(packet%carbon))then
            allocate(carbon(size(packet%carbon)));carbon=packet%carbon
         endif
         if(allocated(packet%carbon_ion))then
            allocate(carbon_ion(size(packet%carbon_ion)));carbon_ion=packet%carbon_ion
         endif
         goto 820
      endif
#else
      restoring=.false.
#endif
      status=1;rho=rows(1,index);start=fraction(:,index);guess=start;solid_fe=0
      old_row=rows(:,index)
      nonthermal=.5d0*sum(rows(2:4,index)**2)/rho+magnetic_energy(rows(:,index))
#ifdef DUST_LIVE
      if(dust_relative_motion)nonthermal=dust_phase_kinetic(rows(1:nvar,index))+magnetic_energy(rows(:,index))
#endif
#if NENER>0
      nonthermal=nonthermal+sum(rows(inener:inener+NENER-1,index))
#endif
      internal=rows(5,index)-nonthermal
      if(.not.ieee_is_finite(internal).or.internal<=0)return
      hcode=rho;hecode=rho*snrt_nhelium_per_hydrogen
      if(atomic_on.or.chimes_on)then
         solid_fe=0;solid_pah=0
         if(dust_iron_enabled())solid_fe=sum(rows(idust_iron:idust_iron+1,index))
         if(dust_pah_enabled())solid_pah=dust_pah_inventory(rows(idust_pah:idust_pah+dust_pah_nstate()-1,index))
         call dust_gas_elements(rows(ichem:ichem+10,index),rows(idust_species:idust_species+1,index),gasx,status, &
              solid_fe,solid_pah)
         if(status/=0)return
         gasx=gasx/rho;hcode=rho*sd*gasx(1)/(atomic_mh*snh);hecode=rho*sd*gasx(2)/(4*atomic_mh*snh)
      endif
      temp=(gamma-1)*internal/rho*st2*snrt_mean_molecular_weight(start(1),start(2),start(3))
      if(atomic_on)temp=atomic_temperature(rho*sd,gasx,start,gamma,internal*energy_unit)
      temp=max(1d0,temp);density=0;opacity=0;scatter=0;dust_tau=0;old_temperature=1;chemical_area=0
      primary_spectrum=0;projection_n=0;projection_e=0;an=0;ae=0;af=0;heat=0;primary_work=0
#ifdef DUST_LIVE
      if(rows(idust,index)<0.or.rows(idust_energy,index)<0)return
      density=rows(idust,index)*sd/snrt_dust_contract_mass_per_h_g
      opacity=snrt_dust_contract_absorption_per_h_cm2(1:ng)
      scatter=snrt_dust_contract_scattering_per_h_cm2(1:ng)
      if(dust_material_composition_enabled().or.dust_iron_enabled())then
         allocate(curve(snrt_dust_contract_number_temperature,1),area(1))
         if(dust_iron_enabled())then
         solid_fe=sum(rows(idust_iron:idust_iron+1,index))
         call iron_compare_curve(snrt_dust_contract_temperature_k(1:size(curve,1)), &
              rows(idust_species:idust_species+1,index),solid_fe,snrt_dust_contract_mass_per_h_g,curve(:,1),status)
         else
         call dust_composition_curve(snrt_dust_contract_temperature_k(1:size(curve,1)), &
              rows(idust_species:idust_species+1,index),snrt_dust_contract_mass_per_h_g,curve(:,1),status)
         endif
         if(status/=0)return
         call dust_composition_area(rows(idust_bins:idust_bins+3,index),snrt_dust_contract_mass_per_h_g,area(1),status)
         if(status/=0)return
      endif
      if(dust_optics_enabled())then
         nb=merge(6,4,dust_iron_enabled());allocate(weights(nb,1));pa=0;ps=0;pg=0
         if(dust_iron_enabled())then
         call iron_compare_weights([rows(idust_bins:idust_bins+3,index),rows(idust_iron:idust_iron+1,index)], &
              weights(:,1),area(1),snrt_dust_contract_mass_per_h_g,status)
         if(status/=0)return
         call fe_six_opacity_basis(snrt_dust_contract_mass_per_h_g,pa,ps,pg,ia,isc,isg,status)
         else
         call d03_cell_weights(rows(idust_bins:idust_bins+3,index),weights(:,1),status)
         if(status/=0)return
         call d03_opacity_basis(snrt_dust_contract_mass_per_h_g,pa(:,1:4),ps(:,1:4),pg(:,1:4), &
              ia(:,1:4),isc(:,1:4),isg(:,1:4),status)
         endif
         if(status/=0)return
         opacity=matmul(pa(:,1:nb),weights(:,1));scatter=matmul(ps(:,1:nb)-pg(:,1:nb),weights(:,1))
         if(dust_fe_uv_enabled())opacity=matmul(pa(:,1:4),weights(1:4,1))
      endif
      dust_tau=density*opacity*chat*cell_dt
      material=rows(idust_energy,index)*energy_unit
      capacity=density*snrt_dust_contract_heat_capacity_per_h_erg_k
      if(snrt_dust_contract_version==4)then
         capacity=1
         if(density>0)then
            if(dust_iron_enabled().or.dust_pah_enabled().or.dust_relative_motion.or.cosmo)then
               ! These receivers use absolute material energy, including
               ! their existing analytic cold branch below the first knot.
               call iron_compare_temperature(rows(idust_species:idust_species+1,index),solid_fe, &
                    rows(idust_energy,index)*sv**2,old_temperature,status)
            else if(allocated(curve))then
               call snrt_dust_material_temperature(snrt_dust_contract_temperature_k(1:size(curve,1)), &
                    curve(:,1),material/density,old_temperature,status)
            else
               call snrt_dust_material_temperature( &
                    snrt_dust_contract_temperature_k(1:snrt_dust_contract_number_temperature), &
                    snrt_dust_contract_internal_energy_per_h_erg(1:snrt_dust_contract_number_temperature), &
                    material/density,old_temperature,status)
            endif
            if(status/=0)return
         endif
      else if(capacity>0)then
         old_temperature=material/capacity
      endif
      if(dust_sublimation_rt_enabled())then
         allocate(sub_bins(4,1),sub_next(4,1));sub_bins(:,1)=rows(idust_bins:idust_bins+3,index)*sd;sub_next=sub_bins
      endif
      if(dust_pah_enabled())then
         call pah_live_prepare(status)
         if(status/=0)return
         allocate(pah_number(dust_pah_nstate(),1))
         do b=1,size(pah_number,1)
            pah_number(b,1)=rows(idust_pah+b-1,index)*sd/dust_pah_state_mass(b)
         enddo
         dust_tau=dust_tau+chat*cell_dt*pah_primary_alpha(pah_number(:,1))
      endif
      if(allocated(area))chemical_area=area(1)*density/(hcode*snh)
      if(dust_iron_enabled().or.dust_pah_charged())then
         call dust_composition_area(rows(idust_bins:idust_bins+3,index),snrt_dust_contract_mass_per_h_g, &
              chemical_area,status)
         if(status/=0)return
         chemical_area=chemical_area*sum(rows(idust_species:idust_species+1,index))*sd/ &
              (hcode*snh*snrt_dust_contract_mass_per_h_g)
      endif
      if(dust_relative_motion)then
         allocate(phase_mass(ndust_phase,1),phase_old(3,ndust_phase,1),phase_work(ndust_phase,1))
         call dust_phase_read(rows(1:nvar,index),phase_mass(:,1),phase_old(:,:,1),gas_mass,gas_momentum,phase_ke,status)
         if(status/=0)return
         phase_mass=phase_mass*sd;phase_old=phase_old*sd*sv;phase_next=phase_old;phase_work=0
         if(dust_pah_enabled())then
            allocate(pah_heat(ng,1),pah_captures(ng,1));pah_heat=0;pah_captures=0
         endif
      endif
#endif
      threshold=[13.60d0,24.59d0,54.42d0];converged=.false.
      do iteration=1,config%max_iterations
         nn=number(:,:,index);ee=energy(:,:,index)
         available=[hcode*(1-start(1)),hecode*(1-start(2)-start(3)),hecode*start(2)]
         call opacity_hhe(guess,hcode,hecode,cell_dt,tau)
         if(chimes_on)then
            tau=0;available=0
         endif
         if(snrt_chimes_band_enabled())then
            an=0;ae=0;af=0
         else if(band_on)then
            call angular_read(nn,ee,rn,re,status)
            if(status/=0)then
               write(*,*)'M_N band angular reconstruction rejected cell/status=',mesh%cells(index),status
               return
            endif
            columns=[hcode*(1-guess(1)),hecode*(1-guess(2)-guess(3)),hecode*guess(2)]*snh*chat*cell_dt
            secondary=c_null_funptr;deposition=0
            if(snrt_node_secondaries_enabled())secondary=c_funloc(snrt_secondary_fractions_c)
#ifdef DUST_LIVE
            grain_columns=0;dummy=0
            if(snrt_d03_band_enabled())then
               grain_columns(1:4)=rows(idust_bins:idust_bins+3,index)*sd*chat*cell_dt
               if(dust_iron_enabled())then
                  grain_columns(5:6)=rows(idust_iron:idust_iron+1,index)*sd*chat*cell_dt
                  status=band_material(mn_live_basis%nq,ng,6,1,rn,re,snrt_group_edges_ev,columns,available, &
                       an,ae,returned,secondary,start(1),deposition,grain_columns,fe_six_band_abs, &
                       fe_six_band_transport,fe_band_ev,mn_live_basis%weight)
               else
                  status=band_material(mn_live_basis%nq,ng,4,1,rn,re,snrt_group_edges_ev,columns,available, &
                       an,ae,returned,secondary,start(1),deposition,grain_columns,d03_band_abs, &
                       d03_band_transport,d03_band_ev,mn_live_basis%weight)
               endif
            else
#endif
               status=band_material(mn_live_basis%nq,ng,0,0,rn,re,snrt_group_edges_ev,columns,available, &
                    an,ae,returned,secondary,start(1),deposition,[0d0],[0d0],[0d0],[0d0],mn_live_basis%weight)
#ifdef DUST_LIVE
            endif
#endif
            if(status/=0)then
               write(*,*)'M_N band spectral material rejected cell/status=',mesh%cells(index),status
               return
            endif
            call angular_write(rn,re,nn,ee,projection_n,projection_e,status)
            if(status/=0)then
               write(*,*)'M_N band angular projection rejected cell/status=',mesh%cells(index),status
               return
            endif
            af=0
         else
         call mn_absorb_groups(mn_live_basis,nn,ee,tau,dust_tau,available,an,ae,af,returned,unassigned,status, &
              validated_input=.true.)
         if(status/=mn_ok)return
         endif
         if(chimes_on)then
            converged=.true.;exit
         endif
         excess=0
         do a=1,3
            do q=1,ng
               if(an(a,q)>0)excess(a,q)=max(0d0,ae(a,q)/an(a,q)-threshold(a))
            enddo
         enddo
         if(snrt_node_secondaries_enabled())then
            call snrt_thermochemistry_advance_cell(hcode*snh,hecode*snh,snh,temp,cell_dt,start(1),start(2),start(3), &
                 an(1:3,:),excess,chemistry,defer_recombination=atomic_on,band_deposition=deposition)
         else
            call snrt_thermochemistry_advance_cell(hcode*snh,hecode*snh,snh,temp,cell_dt,start(1),start(2),start(3), &
                 an(1:3,:),excess,chemistry,defer_recombination=atomic_on)
         endif
         status=chemistry%ierr
         if(status/=0)then
            write(*,*)'M_N HHe chemistry rejected cell/status=',mesh%cells(index),status
            return
         endif
         next_fraction=[chemistry%x_hydrogen_ii,chemistry%x_helium_ii,chemistry%x_helium_iii]
         heat=chemistry%heating_rate_erg_cm3_s*cell_dt
         gas_next=internal*energy_unit+heat
         if(atomic_on)then
            call atomic_advance(rho*sd,gasx,gamma,aexp,cell_dt,gas_next,next_fraction, &
                 gas_energy(1),available,loss,status)
            if(status/=0)return
            gas_next=gas_energy(1);next_fraction=available
         endif
         call opacity_hhe(next_fraction,hcode,hecode,cell_dt,target_tau)
         converged=maxval(abs(next_fraction-guess))<=config%fraction_absolute_tolerance.and. &
              maxval(abs(target_tau-tau)/max(tau,config%tau_floor))<=config%tau_relative_tolerance
         if(converged)exit
         guess=(1-config%relaxation)*guess+config%relaxation*next_fraction
      enddo
      status=1
      if(.not.converged)then
         write(*,*)'M_N HHe Picard interval did not converge cell=',mesh%cells(index)
         status=4;return
      endif
      primary=sum(ae(4,:))*snh*snrt_ev_to_erg
      primary_spectrum(:,1)=ae(4,:)*snh*snrt_ev_to_erg
#ifdef SNRT_CHIMES
      if(chimes_on)then
         if(snrt_chimes_band_enabled())then
            call angular_read(nn,ee,rn,re,status)
            if(status/=0)return
            rn=rn*snh;re=re*snh;pn=rn;pe=re;photo_ledger=0;gn=0;ge=0;ph_e=0;ph_p=0
            if(snrt_chimes_cold_enabled())then
               if(dust_relative_motion)then
                  call chimes_live_cold_stage(mesh%cells(index),sd,sv,cell_dt,mesh%dx*sl,old_temperature,chat_factor, &
                       mn_live_basis%nq,rn,re,chemical,total_chemical_energy,pn,pe,photo_ledger,gn,ge,status,events, &
                       staged_row=rows(:,index),directions=mn_live_basis%direction,phase_energy=ph_e,phase_moment=ph_p)
               else
                  call chimes_live_cold_stage(mesh%cells(index),sd,sv,cell_dt,mesh%dx*sl,old_temperature,chat_factor, &
                       mn_live_basis%nq,rn,re,chemical,total_chemical_energy,pn,pe,photo_ledger,gn,ge,status, &
                       staged_row=rows(:,index))
               endif
               an(4,:)=gn/snh;ae(4,:)=ge/snh
               primary_spectrum(:,1)=ge*snrt_ev_to_erg;primary=sum(primary_spectrum)
            else
               call chimes_live_band_stage(mesh%cells(index),sd,sv,cell_dt,mesh%dx*sl,old_temperature,chat_factor, &
                    mn_live_basis%nq,rn,re,chemical,total_chemical_energy,pn,pe,photo_ledger(1:9),status, &
                    staged_row=rows(:,index))
            endif
            if(status/=0)return
            ! The native chemical receiver owns ALL primary absorption.
            ! Keep its dust channel separate; gas has many species, not just HHe.
            an(1,:)=sum(rn-pn,dim=1)/snh-an(4,:)
            ae(1,:)=sum(re-pe,dim=1)/snh-ae(4,:)
            call angular_write(pn/snh,pe/snh,nn,ee,projection_n,projection_e,status)
            if(status/=0)return
         else
            old_photons=nn(1,:)*snh
            call chimes_live_stage(mesh%cells(index),sd,sv,cell_dt,mesh%dx*sl,old_temperature,chemical_area,chat_factor, &
                 old_photons,chemical,total_chemical_energy,next_photons,status,staged_row=rows(:,index))
            if(status/=0)then
               write(*,*)'M_N fixed-group CHIMES receiver rejected cell/status=',mesh%cells(index),status
               return
            endif
            do q=1,ng
               factor=1
               if(old_photons(q)>0)factor=next_photons(q)/old_photons(q)
               if(factor<0.or.factor>1)then
                  status=1;return
               endif
               an(1,q)=(1-factor)*nn(1,q);ae(1,q)=(1-factor)*ee(1,q)
               ! Frozen group tables use reference energies. Account once
               ! for any actual-minus-reference energy in accepted captures.
               delta_heat=(ae(1,q)-an(1,q)*snrt_group_mean_energy_ev(q))*snh*snrt_ev_to_erg
               total_chemical_energy=total_chemical_energy+delta_heat/energy_unit
               nn(:,q)=factor*nn(:,q);ee(:,q)=factor*ee(:,q)
            enddo
            if(dust_fe_uv_enabled())then
               old_photons=next_photons;fe_ledger=0;material_next=material+primary
               call chimes_live_fe_uv_stage(mesh%cells(index),sd,sv,cell_dt,old_temperature,chat,next_photons, &
                    chemical,total_chemical_energy,material_next,fe_ledger,status,staged_row=rows(:,index))
               if(status/=0)return
               call iron_compare_source_receipt(fe_ledger(3)*snrt_ev_to_erg,material,primary,status)
               if(status/=0)return
               call iron_compare_neutral_area(weights(:,1),snrt_dust_contract_mass_per_h_g, &
                    chemical(2)/(hcode*snh*atomic_mh/sd),area(1),status)
               if(status/=0)return
               do q=1,ng
                  factor=1
                  if(old_photons(q)>0)factor=next_photons(q)/old_photons(q)
                  an(4,q)=an(4,q)+(1-factor)*nn(1,q);ae(4,q)=ae(4,q)+(1-factor)*ee(1,q)
                  nn(:,q)=nn(:,q)*factor;ee(:,q)=ee(:,q)*factor
               enddo
            endif
         endif
         next_fraction(1)=chemical(3)*sd/(hcode*snh*atomic_mh);next_fraction(2:3)=0
         if(hecode>0)next_fraction(2:3)=chemical(6:7)*sd/(hecode*snh*atomic_mh)
         gas_next=(total_chemical_energy-nonthermal)*energy_unit
         heat=gas_next-internal*energy_unit
         rows(ichimes:ichimes+chimes_ns-1,index)=chemical
      endif
#endif
#ifdef DUST_LIVE
      if(dust_pah_enabled())then
         do q=1,ng
            if(sum(pah_number)<=0.or.an(4,q)<=0)cycle
            if(.not.pah_primary_supported(snrt_group_mean_energy_ev(q)))then
               status=1;return
            endif
         enddo
      endif
      if(dust_iron_enabled().and..not.dust_fe_uv_enabled().and..not.band_on)then
         do q=1,ng
            if(an(4,q)>0.and.snrt_group_mean_energy_ev(q)>dust_fe_primary_limit())then
               status=1;return
            endif
         enddo
      endif
      if(dust_relative_motion)then
         phase_gain=0;primary_spectrum=0
         if(snrt_chimes_cold_enabled())then
            do b=1,4
               if(phase_mass(b,1)==0)cycle
               kick=ph_p(:,b)*snrt_ev_to_erg/snrt_c_cgs
               velocity=(phase_old(:,b,1)+.5d0*kick)/phase_mass(b,1)
               channel_heat=ph_e(:,b)*snrt_ev_to_erg
               if(sqrt(sum(velocity**2))>.01d0*snrt_c_cgs.or.dot_product(velocity,kick)>sum(channel_heat))then
                  status=1;return
               endif
               if(sum(channel_heat)>0)channel_heat=channel_heat*(1-dot_product(velocity,kick)/sum(channel_heat))
               primary_spectrum(:,1)=primary_spectrum(:,1)+channel_heat
               phase_next(:,b,1)=phase_old(:,b,1)+kick
            enddo
         else
            phase_sigma=0;phase_fraction=0;phase_tau=0
            do b=1,nb
               phase_sigma(:,b)=pa(:,b)*phase_mass(b,1)/snrt_dust_contract_mass_per_h_g
               phase_tau(b,:)=(ps(:,b)-pg(:,b))*phase_mass(b,1)/snrt_dust_contract_mass_per_h_g*chat*cell_dt
            enddo
            if(dust_pah_enabled())phase_sigma(:,ndust_phase)=sum(pah_number)*pah_primary_sigma
            do q=1,ng
               factor=sum(phase_sigma(q,:))
               if(factor>0)phase_fraction(q,:)=phase_sigma(q,:)/factor
            enddo
            do b=1,ndust_phase
               if(phase_mass(b,1)==0)cycle
               do q=1,ng
                  impulse(:,q)=af(:,4,q)*phase_fraction(q,b)*snh*snrt_ev_to_erg/snrt_c_cgs
               enddo
               kick=sum(impulse,dim=2);velocity=(phase_old(:,b,1)+.5d0*kick)/phase_mass(b,1)
               channel_heat=ae(4,:)*phase_fraction(:,b)*snh*snrt_ev_to_erg-matmul(velocity,impulse)
               if(sqrt(sum(velocity**2))>.01d0*snrt_c_cgs.or.any(channel_heat<0))then
                  status=1;return
               endif
               primary_spectrum(:,1)=primary_spectrum(:,1)+channel_heat
               phase_next(:,b,1)=phase_old(:,b,1)+kick
               if(dust_pah_enabled().and.b==ndust_phase)then
                  pah_heat(:,1)=channel_heat;pah_captures(:,1)=an(4,:)*snh*phase_fraction(:,b)
               endif
            enddo
         endif
         primary=sum(primary_spectrum)
         if(snrt_dust_contract_scattering_enabled)then
            call angular_read(nn,ee,rn,re,status)
            if(status/=0)return
            rn=rn*snh;re=re*snh;work=0
#ifdef SNRT_CHIMES
            if(snrt_chimes_cold_enabled())then
               call chimes_live_grain_scatter(mn_live_basis%nq,rn,re,phase_next(:,:,1),phase_mass(:,1), &
                    mn_live_basis%direction,mn_live_basis%weight,cell_dt,chat_factor,work,status)
            else
#endif
               re=re*snrt_ev_to_erg
               call snrt_moving_scatter_cell(rn,re,phase_next(:,:,1),phase_mass(:,1),phase_tau, &
                    mn_live_basis%direction,mn_live_basis%weight,snrt_c_cgs,work,status)
               re=re/snrt_ev_to_erg
#ifdef SNRT_CHIMES
            endif
#endif
            if(status/=0)return
            primary_work=sum(work)
            call angular_write(rn/snh,re/snh,nn,ee,in_n,in_e,status)
            if(status/=0)return
            projection_n=projection_n+in_n;projection_e=projection_e+in_e
         endif
      endif
      ! Isotropic effective scattering: exactly conserve monopoles; damp all
      ! higher harmonics. This is the existing reduced (sigma_s - sigma_g) model.
      if(snrt_dust_contract_scattering_enabled.and..not.dust_relative_motion.and..not.snrt_d03_band_enabled())then
         do q=1,ng
            factor=exp(-density*scatter(q)*chat*cell_dt)
            nn(2:,q)=nn(2:,q)*factor;ee(2:,q)=ee(2:,q)*factor
         enddo
      endif
      capacity=density*snrt_dust_contract_heat_capacity_per_h_erg_k
      if(snrt_dust_contract_version==4)capacity=1d0
      temperature=0;irproj=0;transfer=0;electron_cv=0
      gas_energy(1)=gas_next
      mu=snrt_mean_molecular_weight(next_fraction(1),next_fraction(2),next_fraction(3))
      gas_capacity(1)=rho*energy_unit/((gamma-1)*st2*mu)
      if(atomic_on)gas_capacity(1)=atomic_heat_capacity(rho*sd,gasx,next_fraction,gamma)
#ifdef SNRT_CHIMES
      if(chimes_on)gas_capacity(1)=chimes_live_capacity(chemical,sd)
      electron_cv=0
      if(dust_pah_charged())then
         electrons=[chemical(1)*sd/atomic_mh];electron_cv=1.5d0*chimes_boltzmann()
         if(dust_pah_hydrogenated())hydrogen=[chemical(2)*sd/atomic_mh]
         if(dust_pah_h2_enabled())molecular_h2=[chemical(138)*sd/atomic_mh]
         if(dust_pah_atomization())then
            carbon=[chemical(8)*sd/atomic_mh];carbon_ion=[chemical(9)*sd/atomic_mh]
         endif
      endif
#endif
      if(dust_pah_enabled())pah_spectrum=primary_spectrum
820   continue
#ifdef DUST_LIVE
      if(present(packet).and.present(defer_ir).and.defer_ir)then
         packet%valid=.false.;packet%has_exchange=snrt_dust_contract_exchange_enabled
         packet%has_curve=allocated(curve);packet%has_weights=allocated(weights)
         packet%has_sublimation=allocated(sub_bins);packet%has_pah=allocated(pah_number)
         packet%has_phase=allocated(phase_mass);packet%has_chimes=chimes_on
         packet%index=index;packet%cell=mesh%cells(index);packet%slot=mesh%slots(index)
         packet%rho=rho;packet%density=density;packet%primary=primary;packet%material=material;packet%capacity=capacity
         packet%gas_next=gas_next;packet%nonthermal=nonthermal;packet%heat=heat
         packet%primary_work=primary_work;packet%projection_n=projection_n;packet%projection_e=projection_e
         packet%electron_cv=electron_cv;packet%phase_ke=phase_ke
         packet%nn=nn;packet%ee=ee;packet%an=an;packet%ae=ae;packet%old_row=old_row
         packet%next_fraction=next_fraction;packet%gasx=gasx;packet%gas_energy=gas_energy
         packet%gas_capacity=gas_capacity;packet%n_hydrogen=[hcode*snh];packet%transfer=transfer
         packet%irproj=irproj
#ifdef SNRT_CHIMES
         if(chimes_on)packet%chemical=chemical
#endif
         if(allocated(curve))packet%curve=curve
         if(allocated(area))packet%area=area
         if(allocated(weights))packet%weights=weights
         if(allocated(sub_bins))packet%sub_bins=sub_bins
         if(allocated(sub_next))packet%sub_next=sub_next
         if(allocated(pah_number))packet%pah_number=pah_number
         if(allocated(pah_spectrum))packet%pah_spectrum=pah_spectrum
         if(allocated(pah_heat))packet%pah_heat=pah_heat
         if(allocated(pah_captures))packet%pah_captures=pah_captures
         if(allocated(phase_mass))packet%phase_mass=phase_mass
         if(allocated(phase_old))packet%phase_old=phase_old
         if(allocated(phase_next))packet%phase_next=phase_next
         if(allocated(phase_work))packet%phase_work=phase_work
         if(allocated(electrons))packet%electrons=electrons
         if(allocated(hydrogen))packet%hydrogen=hydrogen
         if(allocated(molecular_h2))packet%molecular_h2=molecular_h2
         if(allocated(carbon))packet%carbon=carbon
         if(allocated(carbon_ion))packet%carbon_ion=carbon_ion
         packet%valid=.true.;status=0;return
      endif
      if(.not.restoring)then
#endif
      if(snrt_dust_contract_exchange_enabled)then
      call snrt_dust_live_moment_cell(lev,mesh%cells(index),mesh%slots(index),mesh%dx*sl,cell_dt,chat, &
           density,primary,material,capacity,infrared(:,:,index),material_next,temperature, &
           ir_diag,irproj,status,gas_energy=gas_energy,gas_capacity=gas_capacity,n_hydrogen=[hcode*snh], &
           gas_transfer=transfer,cell_material_u=curve,cell_collision_area=area,cell_weights=weights, &
           sublimation_bins=sub_bins,sublimation_next=sub_next,pah_population=pah_number,primary_spectrum=pah_spectrum, &
           primary_pah_heat=pah_heat,primary_pah_captures=pah_captures,phase_density=phase_mass, &
           phase_momentum=phase_next,phase_work=phase_work,gas_electrons=electrons,electron_capacity=electron_cv, &
           gas_atomic_h=hydrogen,gas_molecular_h2=molecular_h2,gas_atomic_c=carbon,gas_carbon_ion=carbon_ion)
      else
      call snrt_dust_live_moment_cell(lev,mesh%cells(index),mesh%slots(index),mesh%dx*sl,cell_dt,chat, &
           density,primary,material,capacity,infrared(:,:,index),material_next,temperature, &
           ir_diag,irproj,status,cell_material_u=curve,cell_collision_area=area,cell_weights=weights)
      endif
#ifdef DUST_LIVE
      endif
#endif
      if(status/=0)then
         write(*,*)'M_N dust/IR receiver rejected cell/status=',mesh%cells(index),status
         return
      endif
      rows(idust_energy,index)=material_next/energy_unit
      ! Positive transfer is gas -> dust; the material receiver already
      ! added it to dust. Subtract it exactly once from the gas reservoir.
      gas_next=gas_next-transfer(1)
      if(dust_relative_motion)then
         rows(2:4,index)=old_row(2:4)+sum(phase_next(:,:,1)-phase_old(:,:,1),dim=2)/(sd*sv)
         do b=1,ndust_phase
            jj=idust_momentum+3*(b-1);rows(jj:jj+2,index)=phase_next(:,b,1)/(sd*sv)
         enddo
         nonthermal=dust_phase_kinetic(rows(1:nvar,index))+magnetic_energy(rows(:,index))
#if NENER>0
         nonthermal=nonthermal+sum(rows(inener:inener+NENER-1,index))
#endif
      endif
      if(dust_pah_enabled())then
         do b=1,dust_pah_nstate()
            rows(idust_pah+b-1,index)=pah_number(b,1)*dust_pah_state_mass(b)/sd
         enddo
      endif
      if(allocated(sub_next))then
         if(dust_relative_motion)then
            phase_ke=dust_phase_kinetic(rows(1:nvar,index))
            gas_momentum=rows(2:4,index)*sd*sv-sum(phase_next(:,:,1),dim=2)
            call dust_phase_erode(sub_bins(:,1),sub_next(:,1),gas_momentum,phase_next(:,1:4,1),status)
            if(status/=0)return
            do b=1,4
               jj=idust_momentum+3*(b-1);rows(jj:jj+2,index)=phase_next(:,b,1)/(sd*sv)
            enddo
         endif
         rows(idust_bins:idust_bins+3,index)=sub_next(:,1)/sd
         rows(idust_species:idust_species+1,index)=[sum(sub_next(1:2,1)),sum(sub_next(3:4,1))]/sd
         rows(idust,index)=sum(sub_next)/sd
         if(dust_relative_motion)then
            delta_heat=(phase_ke-dust_phase_kinetic(rows(1:nvar,index)))*energy_unit
            gas_next=gas_next+delta_heat;nonthermal=nonthermal-delta_heat/energy_unit
         endif
      endif
#ifdef SNRT_CHIMES
      if(chimes_on)then
         if(dust_pah_charged())chemical(1)=electrons(1)*atomic_mh/sd
         if(dust_pah_hydrogenated())chemical(2)=hydrogen(1)*atomic_mh/sd
         if(dust_pah_h2_enabled())chemical(138)=molecular_h2(1)*atomic_mh/sd
         if(dust_pah_atomization())then
            chemical(8)=carbon(1)*atomic_mh/sd;chemical(9)=carbon_ion(1)*atomic_mh/sd
         endif
         rows(ichimes:ichimes+chimes_ns-1,index)=chemical
         if(allocated(sub_next))then
            call chimes_cell_state(mesh%cells(index),rows(idust_species:idust_species+1,index),chemical,elements, &
                 status,staged_row=rows(:,index))
            if(status/=0)return
            rows(ichimes:ichimes+chimes_ns-1,index)=chemical
         endif
         ! PAH H exchange and sublimation can change the gas inventory.
         ! Publish compatibility fractions from the final chemical state,
         ! not the state preceding those accepted material exchanges.
         solid_fe=0;solid_pah=0
         if(dust_iron_enabled())solid_fe=sum(rows(idust_iron:idust_iron+1,index))
         if(dust_pah_enabled())solid_pah=dust_pah_inventory(rows(idust_pah:idust_pah+dust_pah_nstate()-1,index))
         call dust_gas_elements(rows(ichem:ichem+10,index),rows(idust_species:idust_species+1,index), &
              gasx,status,solid_fe,solid_pah)
         if(status/=0)return
         next_fraction(1)=chemical(3)/gasx(1);next_fraction(2:3)=0
         if(gasx(2)>0)next_fraction(2:3)=chemical(6:7)/(gasx(2)/4d0)
      endif
#endif
#endif
      if(.not.ieee_is_finite(gas_next).or.gas_next<=0)then
         status=1;return
      endif
      number(:,:,index)=nn;energy(:,:,index)=ee;fraction(:,index)=next_fraction
      rows(5,index)=nonthermal+gas_next/energy_unit
      receipt(2)=receipt(2)+sum(an)*mesh%volume(index)
      receipt(6)=receipt(6)+heat*mesh%volume(index)*sl**3
      receipt(7)=receipt(7)+primary*mesh%volume(index)*sl**3
      receipt(9)=receipt(9)+(sum(ae)+primary_work/(snh*snrt_ev_to_erg))*mesh%volume(index)
      receipt(4)=receipt(4)+projection_n*mesh%volume(index)
      receipt(11)=receipt(11)+projection_e*mesh%volume(index)
      status=0
    end subroutine

#ifdef DUST_LIVE
    subroutine material_tile(first,last,cell_dt,status)
      ! Keep the expensive CHIMES preparation in the existing per-cell code,
      ! but submit its IR/material write sets as one bounded M_N tile.  This
      ! routine is called by one OpenMP worker per tile; the tile adapter then
      ! owns the only material-stage transaction for all its cells.
      integer,intent(in) :: first,last
      real(dp),intent(in) :: cell_dt
      integer,intent(out) :: status
      type(material_stage_packet),allocatable :: packet(:)
      integer,allocatable :: tile_cells(:),tile_slots(:)
      real(dp),allocatable :: tile_ir(:,:,:),tile_projection(:,:,:),tile_density(:),tile_primary(:)
      real(dp),allocatable :: tile_old(:),tile_capacity(:),tile_material(:),tile_temperature(:)
      real(dp),allocatable :: tile_transfer(:),tile_gas_energy(:),tile_gas_capacity(:),tile_hydrogen(:)
      real(dp),allocatable :: tile_curve(:,:),tile_area(:),tile_weights(:,:)
      real(dp),allocatable :: tile_sub_bins(:,:),tile_sub_next(:,:),tile_pah(:,:),tile_spectrum(:,:)
      real(dp),allocatable :: tile_pah_heat(:,:),tile_pah_captures(:,:)
      real(dp),allocatable :: tile_phase_density(:,:),tile_phase_momentum(:,:,:),tile_phase_work(:,:)
      real(dp),allocatable :: tile_electrons(:),tile_atomic_h(:),tile_molecular_h2(:)
      real(dp),allocatable :: tile_atomic_c(:),tile_carbon_ion(:)
      type(dust_ir_diagnostics) :: tile_diag
      integer :: nc,j,index,status_cell,nb
      logical :: has_exchange,has_curve,has_weights,has_sub,has_pah,has_phase,has_chimes

      status=1;nc=last-first+1
      if(nc<1)return
      allocate(packet(nc),tile_cells(nc),tile_slots(nc))
      do j=1,nc
         index=first+j-1
         call material_cell(index,cell_dt,status_cell,packet(j),.true.)
         if(status_cell/=0)return
      enddo
      has_exchange=packet(1)%has_exchange;has_curve=packet(1)%has_curve
      has_weights=packet(1)%has_weights;has_sub=packet(1)%has_sublimation
      has_pah=packet(1)%has_pah;has_phase=packet(1)%has_phase;has_chimes=packet(1)%has_chimes
      do j=1,nc
         if(.not.packet(j)%valid)return
         if(packet(j)%has_exchange.neqv.has_exchange.or.packet(j)%has_curve.neqv.has_curve.or. &
            packet(j)%has_weights.neqv.has_weights.or.packet(j)%has_sublimation.neqv.has_sub.or. &
            packet(j)%has_pah.neqv.has_pah.or.packet(j)%has_phase.neqv.has_phase.or. &
            packet(j)%has_chimes.neqv.has_chimes)return
      enddo
      if(.not.has_exchange)return

      allocate(tile_ir(nm,ni,nc),tile_projection(nm,ni,nc),tile_density(nc),tile_primary(nc), &
           tile_old(nc),tile_capacity(nc),tile_material(nc),tile_temperature(nc))
      allocate(tile_transfer(nc),tile_gas_energy(nc),tile_gas_capacity(nc),tile_hydrogen(nc))
      tile_ir=0;tile_projection=0;tile_transfer=0;tile_temperature=0
      do j=1,nc
         tile_cells(j)=packet(j)%cell;tile_slots(j)=packet(j)%slot
         tile_ir(:,:,j)=infrared(:,:,packet(j)%index)
         tile_density(j)=packet(j)%density;tile_primary(j)=packet(j)%primary
         tile_old(j)=packet(j)%material;tile_capacity(j)=packet(j)%capacity
         tile_material(j)=packet(j)%material;tile_gas_energy(j)=packet(j)%gas_energy(1)
         tile_gas_capacity(j)=packet(j)%gas_capacity(1);tile_hydrogen(j)=packet(j)%n_hydrogen(1)
      enddo

      if(has_curve)then
         allocate(tile_curve(size(packet(1)%curve,1),nc),tile_area(nc))
         do j=1,nc
            tile_curve(:,j)=packet(j)%curve(:,1);tile_area(j)=packet(j)%area(1)
         enddo
      endif
      if(has_weights)then
         nb=size(packet(1)%weights,1);allocate(tile_weights(nb,nc))
         do j=1,nc;tile_weights(:,j)=packet(j)%weights(:,1);enddo
      endif
      if(has_sub)then
         allocate(tile_sub_bins(size(packet(1)%sub_bins,1),nc),tile_sub_next(size(packet(1)%sub_next,1),nc))
         do j=1,nc
            tile_sub_bins(:,j)=packet(j)%sub_bins(:,1);tile_sub_next(:,j)=packet(j)%sub_next(:,1)
         enddo
      endif
      if(has_pah)then
         allocate(tile_pah(size(packet(1)%pah_number,1),nc),tile_spectrum(size(packet(1)%pah_spectrum,1),nc))
         do j=1,nc
            tile_pah(:,j)=packet(j)%pah_number(:,1);tile_spectrum(:,j)=packet(j)%pah_spectrum(:,1)
         enddo
      endif
      if(has_phase)then
         allocate(tile_phase_density(size(packet(1)%phase_mass,1),nc), &
              tile_phase_momentum(size(packet(1)%phase_next,1),size(packet(1)%phase_next,2),nc), &
              tile_phase_work(size(packet(1)%phase_work,1),nc))
         do j=1,nc
            tile_phase_density(:,j)=packet(j)%phase_mass(:,1)
            tile_phase_momentum(:,:,j)=packet(j)%phase_next(:,:,1)
            tile_phase_work(:,j)=packet(j)%phase_work(:,1)
         enddo
      endif
      if(has_pah.and.allocated(packet(1)%pah_heat))then
         allocate(tile_pah_heat(size(packet(1)%pah_heat,1),nc),tile_pah_captures(size(packet(1)%pah_captures,1),nc))
         do j=1,nc
            tile_pah_heat(:,j)=packet(j)%pah_heat(:,1);tile_pah_captures(:,j)=packet(j)%pah_captures(:,1)
         enddo
      endif
      if(allocated(packet(1)%electrons))then
         allocate(tile_electrons(nc));do j=1,nc;tile_electrons(j)=packet(j)%electrons(1);enddo
      endif
      if(allocated(packet(1)%hydrogen))then
         allocate(tile_atomic_h(nc));do j=1,nc;tile_atomic_h(j)=packet(j)%hydrogen(1);enddo
      endif
      if(allocated(packet(1)%molecular_h2))then
         allocate(tile_molecular_h2(nc));do j=1,nc;tile_molecular_h2(j)=packet(j)%molecular_h2(1);enddo
      endif
      if(allocated(packet(1)%carbon))then
         allocate(tile_atomic_c(nc));do j=1,nc;tile_atomic_c(j)=packet(j)%carbon(1);enddo
      endif
      if(allocated(packet(1)%carbon_ion))then
         allocate(tile_carbon_ion(nc));do j=1,nc;tile_carbon_ion(j)=packet(j)%carbon_ion(1);enddo
      endif

      call snrt_dust_live_moment_tile(lev,tile_cells,tile_slots,mesh%dx*sl,cell_dt,chat, &
           tile_density,tile_primary,tile_old,tile_capacity,tile_ir,tile_material,tile_temperature,tile_diag, &
           tile_projection,status,gas_energy=tile_gas_energy,gas_capacity=tile_gas_capacity,n_hydrogen=tile_hydrogen, &
           gas_transfer=tile_transfer,cell_material_u=tile_curve,cell_collision_area=tile_area,cell_weights=tile_weights, &
           sublimation_bins=tile_sub_bins,sublimation_next=tile_sub_next,pah_population=tile_pah, &
           primary_spectrum=tile_spectrum,primary_pah_heat=tile_pah_heat,primary_pah_captures=tile_pah_captures, &
           phase_density=tile_phase_density,phase_momentum=tile_phase_momentum,phase_work=tile_phase_work, &
           gas_electrons=tile_electrons,electron_capacity=packet(1)%electron_cv,gas_atomic_h=tile_atomic_h, &
           gas_molecular_h2=tile_molecular_h2,gas_atomic_c=tile_atomic_c,gas_carbon_ion=tile_carbon_ion)
      if(status/=0)return

      do j=1,nc
         packet(j)%material_next=tile_material(j);packet(j)%temperature=tile_temperature(j)
         packet(j)%irproj=tile_projection(:,:,j);packet(j)%transfer(1)=tile_transfer(j)
         infrared(:,:,packet(j)%index)=tile_ir(:,:,j)
         if(allocated(packet(j)%sub_next))packet(j)%sub_next(:,1)=tile_sub_next(:,j)
         if(allocated(packet(j)%pah_number))packet(j)%pah_number(:,1)=tile_pah(:,j)
         if(allocated(packet(j)%phase_next))packet(j)%phase_next(:,:,1)=tile_phase_momentum(:,:,j)
         if(allocated(packet(j)%phase_work))packet(j)%phase_work(:,1)=tile_phase_work(:,j)
         if(allocated(packet(j)%electrons))packet(j)%electrons(1)=tile_electrons(j)
         if(allocated(packet(j)%hydrogen))packet(j)%hydrogen(1)=tile_atomic_h(j)
         if(allocated(packet(j)%molecular_h2))packet(j)%molecular_h2(1)=tile_molecular_h2(j)
         if(allocated(packet(j)%carbon))packet(j)%carbon(1)=tile_atomic_c(j)
         if(allocated(packet(j)%carbon_ion))packet(j)%carbon_ion(1)=tile_carbon_ion(j)
         call material_cell(packet(j)%index,cell_dt,status_cell,packet(j),.false.)
         if(status_cell/=0)then;status=status_cell;return;endif
      enddo
      status=0
    end subroutine
#endif

    subroutine angular_read(n,e,rn,re,status)
      real(dp),intent(in)::n(:,:),e(:,:)
      real(dp),intent(out)::rn(:,:),re(:,:)
      integer,intent(out)::status
      real(dp)::a(mn_live_basis%nq),b(mn_live_basis%nq),hint(nm-1)
      integer::q
      do q=1,ng
         hint=0
         call mn_reconstruct(mn_live_basis,n(:,q),a,status,dual_hint=hint)
         if(status/=0)return
         call mn_reconstruct(mn_live_basis,e(:,q),b,status,dual_hint=hint)
         if(status/=0)call mn_reconstruct(mn_live_basis,e(:,q),b,status)
         if(status/=0)return
         if(band_on)then
            call mn_band_join(a,b,snrt_group_edges_ev(q),snrt_group_edges_ev(q+1),rn(:,q),re(:,q))
         else
            rn(:,q)=a;re(:,q)=b
         endif
         rn(:,q)=rn(:,q)*mn_live_basis%weight;re(:,q)=re(:,q)*mn_live_basis%weight
      enddo
    end subroutine

    subroutine angular_write(rn,re,n,e,dn,de,status)
      real(dp),intent(in)::rn(:,:),re(:,:)
      real(dp),intent(out)::n(:,:),e(:,:),dn,de
      integer,intent(out)::status
      real(dp)::a(mn_live_basis%nq),b(mn_live_basis%nq),na,ea
      integer::q
      dn=0;de=0
      do q=1,ng
         if(band_on)then
            call mn_band_split_angular(rn(:,q),re(:,q),snrt_group_edges_ev(q),snrt_group_edges_ev(q+1),a,b,status)
            if(status/=0)then
               write(*,*)'M_N band split rejected group/status/minA/minB=',q,status,minval(a),minval(b)
               return
            endif
         else
            a=rn(:,q);b=re(:,q)
         endif
         call mn_project(mn_live_basis,a/mn_live_basis%weight,n(:,q),status)
         if(status/=0)then
            write(*,*)'M_N first measure projection rejected group/status/min=',q,status,minval(a)
            return
         endif
         call mn_project(mn_live_basis,b/mn_live_basis%weight,e(:,q),status)
         if(status/=0)then
            write(*,*)'M_N second measure projection rejected group/status/min=',q,status,minval(b)
            return
         endif
      enddo
      call inventory(n,e,na,ea);dn=na-sum(rn);de=ea-sum(re)
      call mn_live_validate(n,e,status)
    end subroutine

    subroutine opacity_hhe(x,h,he,cell_dt,tau)
      real(dp),intent(in) :: x(3),h,he,cell_dt
      real(dp),intent(out) :: tau(3,ng)
      tau(1,:)=h*(1-x(1))*snh*chat*cell_dt*snrt_group_cross_section_cm2
      tau(2,:)=he*(1-x(2)-x(3))*snh*chat*cell_dt*snrt_group_cross_section_hei_cm2
      tau(3,:)=he*x(2)*snh*chat*cell_dt*snrt_group_cross_section_heii_cm2
    end subroutine
  end subroutine
end module

! Persistent IR energy uses the primary RT cell-slot map, but a separate
! spectral axis and FP64 physical units (erg/cm3 per normalized direction).
! Trial work never writes persistent radiation or RAMSES material fields.
module snrt_dust_live
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_dust_contract
  use snrt_dust_ir
  use dust_mass_physics, only: dust_optics_enabled,dust_iron_enabled,dust_fe_max_temperature,dust_fe_uv_enabled
  use dust_iron_optics, only: fe_six_opacity_basis,fe_base_binding,fe_radius_cm,fe_density
  use dust_iron_radiation, only: iron_radiative_batch,iron_radiative_cell
  use dust_mass_physics, only: dust_pah_enabled,dust_pah_nstate,dust_pah_molecule_g,dust_pah_charged, &
       dust_pah_hydrogenated,dust_pah_h2_enabled,dust_pah_atomization
  use dust_pah_live_model, only: pah_live_prepare,pah_ir_sigma
  use dust_pah_mixed, only: pah_mixed_advance
  use omp_lib, only: omp_get_max_threads,omp_get_wtime
  use dust_composition_material, only: dust_composition_curve,dust_composition_area
  use dust_sublimation_material, only: dust_radiative_sublimation_evolve
  use dust_composition_optics, only: d03_ng,d03_nir,d03_opacity_basis,d03_cell_weights
  use snrt_runtime_backend, only: snrt_runtime_dust_material, snrt_runtime_ir_transport, snrt_runtime_ir_absorb
  use snrt_runtime_backend, only: snrt_runtime_ir_scatter
  use snrt_runtime_backend, only: snrt_runtime_cpu_material_allowed
  use snrt_state, only: snrt_ndirection, snrt_nslot, snrt_state_get_slot,snrt_state_is_moment
  use snrt_moment_live, only: mn_live_ir_pack,mn_live_ir_unpack,mn_live_ir_expand,mn_live_basis
  use snrt_moment_transport, only: mn_project,mn_reconstruct,mn_ok
  use amr_commons, only: ngridmax,ncoarse,ncpu,myid,headl,next,son,cosmo,aexp
  use snrt_amr_topology, only: snrt_face_kind,snrt_face_cell, &
       snrt_halo_tile_exchange, &
       SNRT_FACE_LOCAL,SNRT_FACE_PHYSICAL,SNRT_FACE_MPI, &
       SNRT_FACE_COARSE_TO_FINE,SNRT_FACE_FINE_TO_COARSE
#ifndef WITHOUTMPI
  use mpi_mod
#endif
#include "amr_index.h"
  implicit none
  private
  ! Stable slabs avoid retaining a full old+new dense IR field on growth.
  ! Pointer components transfer ownership without deep-copying slab data.
  integer, parameter :: ir_slab_slots=128
  type :: ir_slab
     real(dust_dp), pointer :: values(:,:,:)=>null()
  end type
  type(ir_slab), allocatable, save :: radiation(:)
  integer,save::ir_ready_slots=0
  real(dust_dp), allocatable, save :: scattering_basis(:,:)
  real(dust_dp), allocatable, save :: absorption_basis(:,:)
  type(dust_ir_table), save :: table
  logical, save :: initialized=.false.
  type, public :: dust_live_coarse_trial
     integer, allocatable :: slots(:)
     real(dust_dp), allocatable :: energy(:,:,:)
  end type
  public :: snrt_dust_live_stage, snrt_dust_live_commit
  public :: snrt_dust_live_prepare
  public :: snrt_dust_live_pack, snrt_dust_live_restore
  public :: snrt_dust_live_validate_payload
  public :: snrt_dust_live_expand
  public :: snrt_dust_live_moment_cell
  public :: snrt_dust_live_moment_tile
contains
  function ir_slot(slot) result(values)
    integer,intent(in)::slot
    real(dust_dp),pointer::values(:,:)
    values=>radiation((slot-1)/ir_slab_slots+1)%values(:,:,mod(slot-1,ir_slab_slots)+1)
  end function

  subroutine ir_store(slot,values)
    integer,intent(in)::slot
    real(dust_dp),intent(in)::values(:,:)
    real(dust_dp),pointer::destination(:,:)
    destination=>ir_slot(slot)
    destination=values
  end subroutine

  subroutine snrt_dust_live_expand(expansion_ratio,ierr)
    ! Called ONCE at the global clock advance, outside the per-level loop.
    ! This field is physical erg/cm3; fixed-group photon energies therefore
    ! require a^-3 dilution. Spectral redshift is NOT supplied by this hook.
    ! No saved epoch is needed: restart values already belong to restart a.
    real(dust_dp),intent(in)::expansion_ratio
    integer,intent(out)::ierr
    real(dust_dp)::factor
    integer::b
    ierr=dust_err_state
    if(.not.ieee_is_finite(expansion_ratio).or.expansion_ratio<=0)return
    factor=(1d0/expansion_ratio)**3
    if(.not.ieee_is_finite(factor).or.factor<=0)return
    if(snrt_state_is_moment())then
       call mn_live_ir_expand(factor,ierr)
       return
    endif
    if(allocated(radiation))then
       if(factor>1d0)then
          do b=1,size(radiation)
             if(.not.associated(radiation(b)%values))cycle
             if(maxval(radiation(b)%values)>huge(1d0)/factor)return
          enddo
       endif
       do b=1,size(radiation)
          if(.not.associated(radiation(b)%values))cycle
          radiation(b)%values=radiation(b)%values*factor
       enddo
    endif
    ierr=dust_ok
  end subroutine

  subroutine prepare(ierr)
    integer, intent(out) :: ierr
    integer :: ng, nt, old, capacity,nb,b,status
    type(ir_slab), allocatable :: expanded(:)
    real(dust_dp),allocatable::basis(:,:)
    real(dust_dp)::pa(d03_ng,6),ps(d03_ng,6),pg(d03_ng,6),ia(d03_nir,6),isc(d03_nir,6),ig(d03_nir,6)
    ierr=dust_err_config
    if(snrt_dust_contract_version<3.or..not.snrt_dust_contract_runtime_allowed)return
    ! UV signed sensible receipts require the existing implicit, exchanging
    ! material solve, not the absorption-only receiver's precomputed U.
    if(dust_fe_uv_enabled())then
       if(snrt_dust_contract_version/=4.or..not.snrt_dust_contract_exchange_enabled)return
    endif
    ng=snrt_dust_contract_number_ir; nt=snrt_dust_contract_number_temperature
    if(.not.initialized)then
       if(dust_optics_enabled())then
          if(ng/=d03_nir)return
          nb=4
          if(dust_iron_enabled())then
             if(.not.snrt_runtime_cpu_material_allowed())return
             if(snrt_dust_contract_ir_background_k>dust_fe_max_temperature)return
             if(.not.fe_base_binding( &
                  snrt_dust_contract_group_edges_ev(1:snrt_dust_contract_number_groups+1), &
                  snrt_dust_contract_absorption_mean_energy_ev(1:snrt_dust_contract_number_groups), &
                  snrt_dust_contract_ir_energy_ev(1:ng),fe_radius_cm,fe_density(1)))return
             nb=6
             call fe_six_opacity_basis(snrt_dust_contract_mass_per_h_g,pa,ps,pg,ia,isc,ig,ierr)
          else
             call d03_opacity_basis(snrt_dust_contract_mass_per_h_g,pa(:,1:4),ps(:,1:4),pg(:,1:4), &
                  ia(:,1:4),isc(:,1:4),ig(:,1:4),ierr)
          endif
          if(ierr/=0)return
          basis=ia(:,1:nb)
          absorption_basis=ia(:,1:nb)
          scattering_basis=isc(:,1:nb)-ig(:,1:nb)
       endif
       if(snrt_dust_contract_version==4)then
          if(dust_pah_enabled())then
             if(.not.snrt_runtime_cpu_material_allowed())return
             call pah_live_prepare(ierr)
             if(ierr/=0)return
             if(.not.snrt_dust_contract_exchange_enabled)then
                ierr=dust_err_config;return
             endif
             ! The PAH comparison has absolute, initially empty IR radiation.
             ! Its mixed receiver owns both opacity and material components.
             call snrt_dust_ir_initialize(table,snrt_dust_contract_ir_energy_ev(1:ng), &
                  snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
                  snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr)
          else
          call snrt_dust_ir_initialize(table,snrt_dust_contract_ir_energy_ev(1:ng), &
               snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
               snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr, &
               snrt_dust_contract_internal_energy_per_h_erg(1:nt),optical_sigma=basis)
          endif
       else
       call snrt_dust_ir_initialize(table,snrt_dust_contract_ir_energy_ev(1:ng), &
            snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
            snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr)
       endif
       if(ierr/=dust_ok)return
       initialized=.true.
    end if
    if(snrt_state_is_moment())then
       ! Physical table preparation is shared; compact slots must never
       ! acquire a second persistent angular IR field here.
       ierr=dust_ok;return
    endif
    if(snrt_nslot<=ir_ready_slots)then
       ierr=dust_ok;return
    endif
    old=0
    if(allocated(radiation))then
       if(size(radiation)>0)then
       if(associated(radiation(1)%values))then
       if(size(radiation(1)%values,1)/=ng)then
          ierr=dust_err_config
          return
       end if
       endif
       endif
       old=size(radiation)
    end if
    capacity=(snrt_nslot+ir_slab_slots-1)/ir_slab_slots
    if(old<capacity)then
       allocate(expanded(capacity),stat=status)
       if(status/=0)then
          ierr=dust_err_state;return
       endif
       if(old>0)expanded(1:old)=radiation
       call move_alloc(expanded,radiation)
    end if
    do b=1,capacity
       if(associated(radiation(b)%values))cycle
       allocate(radiation(b)%values(ng,snrt_ndirection,ir_slab_slots),stat=status)
       if(status/=0)then
          ierr=dust_err_state;return
       endif
       radiation(b)%values=0
    enddo
    ir_ready_slots=capacity*ir_slab_slots
    ierr=dust_ok
  end subroutine

  subroutine snrt_dust_live_prepare(ierr)
    ! Complete all lazy, shared dust/PAH initialization before the caller
    ! enters a cell-parallel region.  The per-cell stage retains its guard
    ! for serial/restart callers, but must not be the first initializer from
    ! multiple OpenMP workers.
    integer,intent(out)::ierr
    call prepare(ierr)
  end subroutine

  subroutine snrt_dust_live_moment_cell(ilevel,cell,slot,dx,dt,chat,density,primary_energy,old_energy,capacity, &
       energy,material,temperature,diagnostics,projection,ierr, &
       gas_energy,gas_capacity,n_hydrogen,gas_transfer,cell_material_u,cell_collision_area,cell_weights, &
       sublimation_bins,sublimation_next,pah_population,primary_spectrum,primary_pah_heat,primary_pah_captures, &
       phase_density,phase_momentum,phase_work,gas_electrons,electron_capacity,gas_atomic_h,gas_molecular_h2, &
       gas_atomic_c,gas_carbon_ion)
    ! Actual live dust physics, with one-cell angular scratch. The optional
    ! material arrays retain the existing ABI (last extent = one cell).
    ! Neither persistent slots nor RAMSES uold are modified. All mutable
    ! outputs, including PAH/phase/gas state, publish only after final closure.
    integer,intent(in) :: ilevel,cell,slot
    real(dust_dp),intent(in) :: dx,dt,chat,density,primary_energy,old_energy,capacity
    real(dust_dp),intent(inout) :: energy(:,:),material,temperature,projection(:,:)
    type(dust_ir_diagnostics),intent(inout) :: diagnostics
    integer,intent(out) :: ierr
    real(dust_dp),optional,intent(in) :: gas_energy(:),gas_capacity(:),n_hydrogen(:)
    real(dust_dp),optional,intent(inout) :: gas_transfer(:)
    real(dust_dp),optional,intent(in) :: cell_material_u(:,:),cell_collision_area(:),cell_weights(:,:)
    real(dust_dp),optional,intent(in) :: sublimation_bins(:,:)
    real(dust_dp),optional,intent(inout) :: sublimation_next(:,:),pah_population(:,:)
    real(dust_dp),optional,intent(in) :: primary_spectrum(:,:),primary_pah_heat(:,:),primary_pah_captures(:,:)
    real(dust_dp),optional,intent(in) :: phase_density(:,:),electron_capacity
    real(dust_dp),optional,intent(inout) :: phase_momentum(:,:,:),phase_work(:,:),gas_electrons(:)
    real(dust_dp),optional,intent(inout) :: gas_atomic_h(:),gas_molecular_h2(:),gas_atomic_c(:),gas_carbon_ion(:)
    real(dust_dp) :: angular(size(energy,2),mn_live_basis%nq,1),check(mn_live_basis%nq)
    real(dust_dp) :: candidate(size(energy,1),size(energy,2)),receipt(size(energy,1),size(energy,2))
    real(dust_dp) :: next_material(1),next_temperature(1)
    real(dust_dp),allocatable :: next_angular(:,:,:),next_transfer(:),next_sublimation(:,:),next_population(:,:)
    real(dust_dp),allocatable :: next_momentum(:,:,:),next_work(:,:),next_electrons(:),next_h(:),next_h2(:)
    real(dust_dp),allocatable :: next_c(:),next_cp(:)
    type(dust_ir_diagnostics) :: trial
    type(dust_live_coarse_trial) :: coarse
    integer :: g,status,neighbors(6,1)
    ierr=dust_err_config
    if(.not.snrt_state_is_moment().or.mn_live_basis%nm<1)return
    ierr=dust_err_shape
    if(any(shape(energy)/=[mn_live_basis%nm,snrt_dust_contract_number_ir]))return
    if(any(shape(projection)/=shape(energy)))return
    do g=1,size(energy,2)
       call mn_reconstruct(mn_live_basis,energy(:,g),angular(g,:,1),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
       call mn_project(mn_live_basis,angular(g,:,1),candidate(:,g),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
    enddo
    ! Numerical reclosure is a separate density receipt, not material heat.
    receipt=candidate-energy;neighbors=0
    if(present(gas_transfer))allocate(next_transfer(size(gas_transfer)),source=0d0)
    if(present(sublimation_next))allocate(next_sublimation(size(sublimation_next,1),size(sublimation_next,2)),source=0d0)
    if(present(pah_population))next_population=pah_population
    if(present(phase_momentum))next_momentum=phase_momentum
    if(present(phase_work))next_work=phase_work
    if(present(gas_electrons))next_electrons=gas_electrons
    if(present(gas_atomic_h))next_h=gas_atomic_h
    if(present(gas_molecular_h2))next_h2=gas_molecular_h2
    if(present(gas_atomic_c))next_c=gas_atomic_c
    if(present(gas_carbon_ion))next_cp=gas_carbon_ion
    call snrt_dust_live_stage(ilevel,[cell],[slot],neighbors,mn_live_basis%direction,mn_live_basis%weight, &
         dx,dt,chat,[density],[primary_energy],[old_energy],[capacity], &
         next_angular,next_material,next_temperature,trial,ierr,coarse, &
         gas_energy=gas_energy,gas_capacity=gas_capacity,n_hydrogen=n_hydrogen,gas_transfer=next_transfer, &
         cell_material_u=cell_material_u,cell_collision_area=cell_collision_area,cell_weights=cell_weights, &
         sublimation_bins=sublimation_bins,sublimation_next=next_sublimation,pah_population=next_population, &
         primary_spectrum=primary_spectrum,primary_pah_heat=primary_pah_heat,primary_pah_captures=primary_pah_captures, &
         phase_density=phase_density,phase_momentum=next_momentum,phase_work=next_work, &
         gas_electrons=next_electrons,electron_capacity=electron_capacity,gas_atomic_h=next_h,gas_molecular_h2=next_h2, &
         gas_atomic_c=next_c,gas_carbon_ion=next_cp,incoming_radiation=angular)
    if(ierr/=dust_ok)return
    do g=1,size(energy,2)
       call mn_project(mn_live_basis,next_angular(g,:,1),candidate(:,g),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
       call mn_reconstruct(mn_live_basis,candidate(:,g),check,status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
    enddo
    energy=candidate;material=next_material(1);temperature=next_temperature(1);projection=receipt;diagnostics=trial
    if(present(gas_transfer))gas_transfer=next_transfer
    if(present(sublimation_next))sublimation_next=next_sublimation
    if(present(pah_population))pah_population=next_population
    if(present(phase_momentum))phase_momentum=next_momentum
    if(present(phase_work))phase_work=next_work
    if(present(gas_electrons))gas_electrons=next_electrons
    if(present(gas_atomic_h))gas_atomic_h=next_h
    if(present(gas_molecular_h2))gas_molecular_h2=next_h2
    if(present(gas_atomic_c))gas_atomic_c=next_c
    if(present(gas_carbon_ion))gas_carbon_ion=next_cp
  end subroutine

  subroutine snrt_dust_live_moment_tile(ilevel,cells,slots,dx,dt,chat,density,primary_energy,old_energy,capacity, &
       energy,material,temperature,diagnostics,projection,ierr, &
       gas_energy,gas_capacity,n_hydrogen,gas_transfer,cell_material_u,cell_collision_area,cell_weights, &
       sublimation_bins,sublimation_next,pah_population,primary_spectrum,primary_pah_heat,primary_pah_captures, &
       phase_density,phase_momentum,phase_work,gas_electrons,electron_capacity,gas_atomic_h,gas_molecular_h2, &
       gas_atomic_c,gas_carbon_ion)
    ! Batch counterpart of snrt_dust_live_moment_cell.  The incoming M_N
    ! moments are reconstructed for all cells, then one material-only stage
    ! owns the complete tile transaction.  Spatial transport remains outside
    ! this adapter: incoming_radiation is the conservative angular state
    ! already produced by the M_N transport operator.
    !
    ! All optional mutable fields are copied to tile-local next-state arrays.
    ! Thus a rejected cell or material solve cannot partially publish gas,
    ! PAH, sublimation, or phase updates.  The public arrays are assigned only
    ! after every cell has passed projection/reconstruction closure.
    integer,intent(in) :: ilevel,cells(:),slots(:)
    real(dust_dp),intent(in) :: dx,dt,chat,density(:),primary_energy(:),old_energy(:),capacity(:)
    real(dust_dp),intent(inout) :: energy(:,:,:),material(:),temperature(:),projection(:,:,:)
    type(dust_ir_diagnostics),intent(inout) :: diagnostics
    integer,intent(out) :: ierr
    real(dust_dp),optional,intent(in) :: gas_energy(:),gas_capacity(:),n_hydrogen(:)
    real(dust_dp),optional,intent(inout) :: gas_transfer(:)
    real(dust_dp),optional,intent(in) :: cell_material_u(:,:),cell_collision_area(:),cell_weights(:,:)
    real(dust_dp),optional,intent(in) :: sublimation_bins(:,:)
    real(dust_dp),optional,intent(inout) :: sublimation_next(:,:),pah_population(:,:)
    real(dust_dp),optional,intent(in) :: primary_spectrum(:,:),primary_pah_heat(:,:),primary_pah_captures(:,:)
    real(dust_dp),optional,intent(in) :: phase_density(:,:),electron_capacity
    real(dust_dp),optional,intent(inout) :: phase_momentum(:,:,:),phase_work(:,:),gas_electrons(:)
    real(dust_dp),optional,intent(inout) :: gas_atomic_h(:),gas_molecular_h2(:),gas_atomic_c(:),gas_carbon_ion(:)
    integer :: nc,ng,nm,nq,i,g,status
    integer,allocatable :: neighbors(:,:)
    real(dust_dp),allocatable :: angular(:,:,:),candidate(:,:,:),receipt(:,:,:)
    real(dust_dp),allocatable :: next_angular(:,:,:),next_material(:),next_temperature(:)
    real(dust_dp),allocatable :: next_transfer(:),next_sublimation(:,:),next_population(:,:)
    real(dust_dp),allocatable :: next_momentum(:,:,:),next_work(:,:),next_electrons(:)
    real(dust_dp),allocatable :: next_h(:),next_h2(:),next_c(:),next_cp(:)
    real(dust_dp),allocatable :: check(:)
    type(dust_ir_diagnostics) :: trial
    type(dust_live_coarse_trial) :: coarse

    ierr=dust_err_config
    if(.not.snrt_state_is_moment().or.mn_live_basis%nm<1)return
    nc=size(cells);ng=snrt_dust_contract_number_ir;nm=mn_live_basis%nm;nq=mn_live_basis%nq
    if(nc<1)return
    ierr=dust_err_shape
    if(size(slots)/=nc.or.size(density)/=nc.or.size(primary_energy)/=nc.or.size(old_energy).or. &
       size(capacity)/=nc.or.size(material)/=nc.or.size(temperature)/=nc)return
    if(any(shape(energy)/=[nm,ng,nc]).or.any(shape(projection)/=shape(energy)))return
    if(size(mn_live_basis%direction,1)/=3.or.size(mn_live_basis%direction,2)/=nq.or. &
       size(mn_live_basis%weight)/=nq)return
    if(present(gas_energy))then
       if(size(gas_energy)/=nc)return
    endif
    if(present(gas_capacity))then
       if(size(gas_capacity)/=nc)return
    endif
    if(present(n_hydrogen))then
       if(size(n_hydrogen)/=nc)return
    endif
    if(present(gas_transfer))then
       if(size(gas_transfer)/=nc)return
    endif

    allocate(angular(ng,nq,nc),candidate(nm,ng,nc),receipt(nm,ng,nc),check(nq),neighbors(6,nc))
    neighbors=0
    do i=1,nc
       do g=1,ng
          call mn_reconstruct(mn_live_basis,energy(:,g,i),angular(g,:,i),status)
          if(status/=mn_ok)then
             ierr=100+status;return
          endif
          call mn_project(mn_live_basis,angular(g,:,i),candidate(:,g,i),status)
          if(status/=mn_ok)then
             ierr=100+status;return
          endif
       enddo
    enddo
    ! Numerical reclosure is a separate density receipt, not material heat.
    receipt=candidate-energy

    ! These local copies are the tile's transaction write set.  In particular,
    ! do not pass caller-owned inout arrays into the stage itself.
    if(present(gas_transfer))allocate(next_transfer(nc),source=0d0)
    if(present(sublimation_next))then
       allocate(next_sublimation(size(sublimation_next,1),size(sublimation_next,2)),source=0d0)
    endif
    if(present(pah_population))then
       allocate(next_population(size(pah_population,1),size(pah_population,2)))
       next_population=pah_population
    endif
    if(present(phase_momentum))then
       allocate(next_momentum(size(phase_momentum,1),size(phase_momentum,2),size(phase_momentum,3)))
       next_momentum=phase_momentum
    endif
    if(present(phase_work))then
       allocate(next_work(size(phase_work,1),size(phase_work,2)))
       next_work=phase_work
    endif
    if(present(gas_electrons))then
       allocate(next_electrons(nc));next_electrons=gas_electrons
    endif
    if(present(gas_atomic_h))then
       allocate(next_h(nc));next_h=gas_atomic_h
    endif
    if(present(gas_molecular_h2))then
       allocate(next_h2(nc));next_h2=gas_molecular_h2
    endif
    if(present(gas_atomic_c))then
       allocate(next_c(nc));next_c=gas_atomic_c
    endif
    if(present(gas_carbon_ion))then
       allocate(next_cp(nc));next_cp=gas_carbon_ion
    endif
    allocate(next_material(nc),next_temperature(nc))

    call snrt_dust_live_stage(ilevel,cells,slots,neighbors,mn_live_basis%direction,mn_live_basis%weight, &
         dx,dt,chat,density,primary_energy,old_energy,capacity,next_angular,next_material,next_temperature, &
         trial,ierr,coarse,gas_energy=gas_energy,gas_capacity=gas_capacity,n_hydrogen=n_hydrogen, &
         gas_transfer=next_transfer,cell_material_u=cell_material_u,cell_collision_area=cell_collision_area, &
         cell_weights=cell_weights,sublimation_bins=sublimation_bins,sublimation_next=next_sublimation, &
         pah_population=next_population,primary_spectrum=primary_spectrum,primary_pah_heat=primary_pah_heat, &
         primary_pah_captures=primary_pah_captures,phase_density=phase_density,phase_momentum=next_momentum, &
         phase_work=next_work,gas_electrons=next_electrons,electron_capacity=electron_capacity,gas_atomic_h=next_h, &
         gas_molecular_h2=next_h2,gas_atomic_c=next_c,gas_carbon_ion=next_cp,incoming_radiation=angular)
    if(ierr/=dust_ok)return
    do i=1,nc
       do g=1,ng
          call mn_project(mn_live_basis,next_angular(g,:,i),candidate(:,g,i),status)
          if(status/=mn_ok)then
             ierr=100+status;return
          endif
          call mn_reconstruct(mn_live_basis,candidate(:,g,i),check,status)
          if(status/=mn_ok)then
             ierr=100+status;return
          endif
       enddo
    enddo

    energy=candidate;material=next_material;temperature=next_temperature;projection=receipt;diagnostics=trial
    if(present(gas_transfer))gas_transfer=next_transfer
    if(present(sublimation_next))sublimation_next=next_sublimation
    if(present(pah_population))pah_population=next_population
    if(present(phase_momentum))phase_momentum=next_momentum
    if(present(phase_work))phase_work=next_work
    if(present(gas_electrons))gas_electrons=next_electrons
    if(present(gas_atomic_h))gas_atomic_h=next_h
    if(present(gas_molecular_h2))gas_molecular_h2=next_h2
    if(present(gas_atomic_c))gas_atomic_c=next_c
    if(present(gas_carbon_ion))gas_carbon_ion=next_cp
    ierr=dust_ok
  end subroutine

  subroutine snrt_dust_live_stage(ilevel,cells,slots,neighbors,directions,weights,dx,dt,chat, &
       density,primary_energy,old_energy,capacity,trial,material,temperature,diagnostics,ierr,coarse, &
       gas_energy,gas_capacity,n_hydrogen,gas_transfer,cell_material_u,cell_collision_area,cell_weights, &
       sublimation_bins,sublimation_next,pah_population,primary_spectrum,primary_pah_heat,primary_pah_captures, &
       phase_density,phase_momentum,phase_work,gas_electrons,electron_capacity,gas_atomic_h,gas_molecular_h2, &
       gas_atomic_c,gas_carbon_ion,incoming_radiation)
    integer, intent(in) :: ilevel,cells(:),slots(:),neighbors(:,:)
    real(dust_dp), intent(in) :: directions(:,:),weights(:),dx,dt,chat
    real(dust_dp), intent(in) :: density(:),primary_energy(:),old_energy(:),capacity(:)
    real(dust_dp), allocatable, intent(out) :: trial(:,:,:)
    real(dust_dp), intent(out) :: material(:),temperature(:)
    type(dust_ir_diagnostics), intent(out) :: diagnostics
    integer, intent(out) :: ierr
    type(dust_live_coarse_trial), intent(out) :: coarse
    real(dust_dp),optional,intent(in)::gas_energy(:),gas_capacity(:),n_hydrogen(:)
    real(dust_dp),optional,intent(out)::gas_transfer(:)
    real(dust_dp),optional,intent(in)::cell_material_u(:,:),cell_collision_area(:)
    real(dust_dp),optional,intent(in)::cell_weights(:,:)
    real(dust_dp),optional,intent(in)::sublimation_bins(:,:)
    real(dust_dp),optional,intent(out)::sublimation_next(:,:)
    real(dust_dp),optional,intent(inout)::pah_population(:,:)
    real(dust_dp),optional,intent(in)::primary_spectrum(:,:)
    real(dust_dp),optional,intent(in)::primary_pah_heat(:,:),primary_pah_captures(:,:)
    ! Per-phase cgs absolute momenta are private to this stage until every
    ! IR substep and rank succeeds. Bulk order is the optical basis order;
    ! optional PAH fluid follows the 4/6 bulk phases, not 128 excitation bins.
    real(dust_dp),optional,intent(in)::phase_density(:,:)
    real(dust_dp),optional,intent(inout)::phase_momentum(:,:,:),phase_work(:,:)
    real(dust_dp),optional,intent(inout)::gas_electrons(:)
    real(dust_dp),optional,intent(in)::electron_capacity
    real(dust_dp),optional,intent(inout)::gas_atomic_h(:)
    real(dust_dp),optional,intent(inout)::gas_molecular_h2(:)
    real(dust_dp),optional,intent(inout)::gas_atomic_c(:),gas_carbon_ion(:)
    ! Supplying radiation selects a strictly local material stage. Spatial
    ! M_N transport is done by the conservative moment operator, never again
    ! by the legacy IR stencil. This path performs NO MPI collectives, so
    ! the caller may use bounded per-cell/tile scratch on unequal rank loads.
    real(dust_dp),optional,intent(in)::incoming_radiation(:,:,:)
    real(dust_dp),allocatable::electron_work(:),electron_initial(:),pah_capacity(:)
    real(dust_dp),allocatable::primary_pah_population(:,:)
    real(dust_dp),allocatable::h_work(:),h_initial(:)
    real(dust_dp),allocatable::h2_work(:),h2_initial(:)
    real(dust_dp),allocatable::c_work(:),c_initial(:),cp_work(:),cp_initial(:)
    real(dust_dp),allocatable::phase_p_stage(:,:,:),phase_alpha(:,:,:),phase_sca(:,:,:)
    real(dust_dp),allocatable::phase_work_step(:,:),phase_work_sum(:,:)
    real(dust_dp),allocatable::pah_bulk_bins(:,:)
    real(dust_dp),allocatable::pah_heat_step(:,:),pah_captures_step(:,:)
    real(dust_dp),allocatable::sub_bins(:,:),sub_trial(:,:),sub_phase(:),sub_sensible(:)
    real(dust_dp),allocatable::sub_u(:,:),sub_area(:),sub_weights(:,:),sub_density(:)
    real(dust_dp),allocatable::gas_work(:),conductance(:),exchange(:),exchange_sum(:)
    real(dust_dp),allocatable::scatter_sigma(:,:)
    real(dust_dp),parameter::kb=1.380649d-16,mp=1.67262192369d-24
    real(dust_dp), allocatable :: photons(:,:),field(:)
    real(dust_dp), allocatable,target :: ghosts(:,:,:)
    real(dust_dp),pointer :: ghost_arg(:,:,:)
    integer,pointer :: remote_arg(:,:)
    procedure(dust_transport_dispatch),pointer :: transport_arg
    real(dust_dp), allocatable :: halo_field(:,:)
    integer,parameter :: halo_tile=16
    integer :: first_component,ncomponent,component,column
    integer, allocatable,target :: remote(:,:)
    integer, allocatable :: ghost_cells(:),ghost_kind(:),coarse_cells(:)
    logical :: local_material
    logical, allocatable :: blocked(:,:)
    type(dust_ir_diagnostics) :: step
    real(dust_dp) :: cfl,step_dt,mu,q
    real(dust_dp) :: phase_start,halo_wall,solve_wall,scatter_wall
    integer :: nsub,isub,i,ng,face,k,nghost,nfield,g,d,global_nsub,info,has_coarse,b,nbulk
    integer :: grid,child,cell,ncoarse_leaf,axis,nd
    local_material=present(incoming_radiation);nd=size(weights)
    nullify(ghost_arg,remote_arg,transport_arg)
    if(.not.local_material)transport_arg=>snrt_runtime_ir_transport
    call prepare(ierr)
    if(ierr==dust_ok)call validate_stage()
    call stage_error(ierr)
    if(ierr/=dust_ok)return
    if(present(cell_weights))scatter_sigma=matmul(scattering_basis,cell_weights)
    if(present(phase_density))then
       nbulk=size(absorption_basis,2)
       phase_p_stage=phase_momentum
       allocate(phase_alpha(size(phase_density,1),size(absorption_basis,1),size(slots)))
       allocate(phase_sca(size(phase_density,1),size(absorption_basis,1),size(slots)))
       allocate(phase_work_step(size(phase_density,1),size(slots)),phase_work_sum(size(phase_density,1),size(slots)))
       phase_alpha=0;phase_sca=0;phase_work_sum=0;phase_work_step=0
       do i=1,size(slots)
          do b=1,nbulk
             phase_alpha(b,:,i)=absorption_basis(:,b)*phase_density(b,i)/snrt_dust_contract_mass_per_h_g
             phase_sca(b,:,i)=scattering_basis(:,b)*phase_density(b,i)/snrt_dust_contract_mass_per_h_g
          enddo
          if(dust_pah_enabled())phase_alpha(nbulk+1,:,i)= &
               pah_ir_sigma*phase_density(nbulk+1,i)/dust_pah_molecule_g
       enddo
       ! Current neutral PAH model supplies absorption/emission only. No
       ! unprovided PAH scattering cross section is borrowed from bulk dust.
    endif
#ifndef WITHOUTMPI
    if(.not.local_material)then
    call MPI_ALLREDUCE(nsub,global_nsub,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,k)
    nsub=global_nsub
    endif
#endif
    step_dt=dt/nsub
    if(present(primary_pah_heat))then
       pah_heat_step=primary_pah_heat*(step_dt/dt)
       pah_captures_step=primary_pah_captures*(step_dt/dt)
    endif
    if(present(gas_energy))then
       gas_work=gas_energy
       allocate(conductance(size(slots)),exchange(size(slots)),exchange_sum(size(slots)))
       exchange_sum=0
    endif
    if(dust_pah_enabled())then
       allocate(pah_bulk_bins(6,size(slots)));pah_bulk_bins=0
       pah_capacity=gas_capacity
       if(dust_pah_charged())then
          electron_work=gas_electrons;electron_initial=gas_electrons
          if(dust_pah_atomization())then
             c_work=gas_atomic_c;c_initial=gas_atomic_c
             cp_work=gas_carbon_ion;cp_initial=gas_carbon_ion
          endif
          primary_pah_population=pah_population
          if(dust_pah_hydrogenated())then
             h_work=gas_atomic_h;h_initial=gas_atomic_h
             if(dust_pah_h2_enabled())then
                h2_work=gas_molecular_h2;h2_initial=gas_molecular_h2
             endif
          endif
       endif
       do i=1,size(slots)
          pah_bulk_bins(1:size(cell_weights,1),i)=cell_weights(:,i)*density(i)*snrt_dust_contract_mass_per_h_g
       enddo
    endif
    if(present(sublimation_bins))then
       sub_bins=sublimation_bins
       allocate(sub_trial(4,size(slots)),sub_phase(size(slots)),sub_sensible(size(slots)), &
            sub_area(size(slots)),sub_density(size(slots)))
       allocate(sub_u(snrt_dust_contract_number_temperature,size(slots)),sub_weights(4,size(slots)))
    endif
    ng=snrt_dust_contract_number_ir
    has_coarse=0
    if(size(slots)>0.and..not.local_material)then
       if(any(snrt_face_kind==SNRT_FACE_FINE_TO_COARSE))has_coarse=1
    end if
    ! Use the same global decision on empty ranks and on coarse-only owners.
    call stage_error(has_coarse)
    ierr=dust_ok
    if(has_coarse/=0)then
       if(ilevel<=1)ierr=dust_err_config
       if(.not.allocated(headl).or..not.allocated(next).or..not.allocated(son))ierr=dust_err_config
    end if
    call stage_error(ierr)
    if(ierr/=dust_ok)return
    ncoarse_leaf=0
    if(has_coarse/=0)then
       grid=headl(myid,ilevel-1)
       do while(grid>0)
          do child=1,twotondim
             cell=ICELL_OF(grid,child)
             if(son(cell)==0)ncoarse_leaf=ncoarse_leaf+1
          end do
          grid=next(grid)
       end do
    end if
    allocate(coarse_cells(ncoarse_leaf),coarse%slots(ncoarse_leaf))
    allocate(coarse%energy(ng,nd,ncoarse_leaf))
    k=0
    if(has_coarse/=0)then
       grid=headl(myid,ilevel-1)
       do while(grid>0)
          do child=1,twotondim
             cell=ICELL_OF(grid,child)
             if(son(cell)/=0)cycle
             k=k+1; coarse_cells(k)=cell; coarse%slots(k)=snrt_state_get_slot(cell)
          end do
          grid=next(grid)
       end do
    end if
    ierr=dust_ok
    if(any(coarse%slots<1).or.any(coarse%slots>snrt_nslot))ierr=dust_err_state
    call stage_error(ierr)
    if(ierr/=dust_ok)return
    do i=1,ncoarse_leaf
       coarse%energy(:,:,i)=ir_slot(coarse%slots(i))
    enddo
    allocate(trial(ng,nd,size(slots)),photons(ng,size(slots)))
    if(local_material)then
       trial=incoming_radiation
    else
    do i=1,size(slots)
       trial(:,:,i)=ir_slot(slots(i))
    enddo
    endif
    photons=0
    material=old_energy
    temperature=material/capacity
    nghost=0
    allocate(remote(6,size(slots)),blocked(6,size(slots))); remote=0; blocked=.false.
    if(size(slots)>0.and..not.local_material)then
       nghost=count(snrt_face_kind==SNRT_FACE_MPI.or.snrt_face_kind==SNRT_FACE_FINE_TO_COARSE)
       blocked=snrt_face_kind==SNRT_FACE_COARSE_TO_FINE
    end if
    allocate(ghost_cells(nghost),ghost_kind(nghost),ghosts(ng,nd,nghost))
    if(local_material)then
       blocked=.true.
       allocate(field(0))
    else
    allocate(field(nfield))
    ghost_arg=>ghosts;remote_arg=>remote
    if(ncpu>1)allocate(halo_field(nfield,halo_tile))
    k=0
    do i=1,size(slots)
       do face=1,6
          if(snrt_face_kind(face,i)/=SNRT_FACE_MPI.and.snrt_face_kind(face,i)/=SNRT_FACE_FINE_TO_COARSE)cycle
          k=k+1; remote(face,i)=k; ghost_cells(k)=snrt_face_cell(face,i)
          ghost_kind(k)=snrt_face_kind(face,i)
       end do
    end do
    endif
    diagnostics=dust_ir_diagnostics()
    if(has_coarse/=0)then
       ! A missing coarse donor must not silently appear as vacuum. The
       ! owner marker is exchanged using the same mapping as the energy.
       field=0; field(coarse_cells)=1
       if(ncpu>1)call make_virtual_fine_dp(field,ilevel-1)
       ierr=dust_ok
       do k=1,nghost
          if(ghost_kind(k)/=SNRT_FACE_FINE_TO_COARSE)cycle
          if(field(ghost_cells(k))/=1)ierr=dust_err_state
       end do
       call stage_error(ierr)
       if(ierr/=dust_ok)return
    end if
    halo_wall=0;solve_wall=0;scatter_wall=0
    do isub=1,nsub
       phase_start=omp_get_wtime()
       ! Reuse RAMSES' real halo communicator, in FP64. All ranks participate
       ! in every substep, including a rank with no local MPI boundary faces.
       if(ncpu>1.and..not.local_material)then
          do first_component=1,ng*nd,halo_tile
             ncomponent=min(halo_tile,ng*nd-first_component+1)
             halo_field=0d0
             ! Cell-major packing reuses adjacent frequencies in each cache
             ! line instead of rereading the entire strided field per column.
!$omp parallel do private(i,column,component,g,d) schedule(static)
             do i=1,size(cells)
                do column=1,ncomponent
                   component=first_component+column-1
                   g=mod(component-1,ng)+1;d=(component-1)/ng+1
                   halo_field(cells(i),column)=trial(g,d,i)
                enddo
             enddo
!$omp end parallel do
             call snrt_halo_tile_exchange(halo_field(:,1:ncomponent),ilevel,ierr)
             call stage_error(ierr)
             if(ierr/=dust_ok)return
!$omp parallel do private(k,column,component,g,d) schedule(static)
             do k=1,nghost
                if(ghost_kind(k)/=SNRT_FACE_MPI)cycle
                do column=1,ncomponent
                   component=first_component+column-1
                   g=mod(component-1,ng)+1;d=(component-1)/ng+1
                   ghosts(g,d,k)=halo_field(ghost_cells(k),column)
                enddo
             end do
!$omp end parallel do
          end do
       end if
       if(has_coarse/=0)then
          do d=1,nd
             do g=1,ng
                field=0; field(coarse_cells)=coarse%energy(g,d,:)
                if(ncpu>1)call make_virtual_fine_dp(field,ilevel-1)
                do k=1,nghost
                   if(ghost_kind(k)==SNRT_FACE_FINE_TO_COARSE)ghosts(g,d,k)=field(ghost_cells(k))
                end do
                ! Fine-volume/coarse-volume = 1/8 in 3D. Use the OLD
                ! upwind field, exactly as the native transport operator.
                field=0
                do i=1,size(slots)
                   do face=1,6
                      if(snrt_face_kind(face,i)/=SNRT_FACE_FINE_TO_COARSE)cycle
                      k=remote(face,i); cell=ghost_cells(k); axis=(face+1)/2
                      mu=directions(axis,d)
                      if(mod(face,2)==1)mu=-mu
                      q=ghosts(g,d,k)
                      if(mu>=0)q=trial(g,d,i)
                      field(cell)=field(cell)+chat*step_dt/dx*mu*q/real(twotondim,dust_dp)
                   end do
                end do
                if(ncpu>1)call make_virtual_reverse_dp(field,ilevel-1)
                coarse%energy(g,d,:)=coarse%energy(g,d,:)+field(coarse_cells)
             end do
          end do
       end if
       halo_wall=halo_wall+omp_get_wtime()-phase_start
       phase_start=omp_get_wtime()
       step=dust_ir_diagnostics()
       ierr=dust_ok
       if(dust_pah_enabled())then
          if(dust_pah_charged())pah_capacity=gas_capacity+electron_capacity*(electron_work-electron_initial)
          if(dust_pah_hydrogenated())pah_capacity=pah_capacity+electron_capacity*(h_work-h_initial)
          if(dust_pah_h2_enabled())pah_capacity=pah_capacity+electron_capacity*(h2_work-h2_initial)
          if(dust_pah_atomization())pah_capacity=pah_capacity+ &
               electron_capacity*(c_work-c_initial+cp_work-cp_initial)
          conductance=2*kb*n_hydrogen*density*cell_collision_area* &
               snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/pah_capacity))
          if(size(slots)>0)call pah_mixed_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
               pah_bulk_bins,primary_spectrum*(step_dt/dt),trial,pah_population,material,temperature,photons,step,ierr, &
               ghost_arg,remote_arg,blocked,gas_work,pah_capacity,conductance,exchange, &
               primary_pah_heat=pah_heat_step,primary_pah_captures=pah_captures_step, &
               phase_density=phase_density,phase_momentum=phase_p_stage,phase_absorption=phase_alpha, &
               phase_scattering=phase_sca,phase_work=phase_work_step, &
               gas_electrons=electron_work,electron_capacity=electron_capacity,primary_population=primary_pah_population, &
               gas_atomic_h=h_work,gas_molecular_h2=h2_work,gas_atomic_c=c_work,gas_carbon_ion=cp_work, &
               material_only=local_material)
          if(ierr==dust_ok)exchange_sum=exchange_sum+exchange
       else if(allocated(sub_bins))then
          do i=1,size(slots)
             sub_density(i)=sum(sub_bins(:,i))/snrt_dust_contract_mass_per_h_g
             call dust_composition_curve(snrt_dust_contract_temperature_k(1:size(sub_u,1)), &
                  [sum(sub_bins(1:2,i)),sum(sub_bins(3:4,i))],snrt_dust_contract_mass_per_h_g,sub_u(:,i),k)
             if(k/=0)ierr=dust_err_state
             call dust_composition_area(sub_bins(:,i),snrt_dust_contract_mass_per_h_g,sub_area(i),k)
             if(k/=0)ierr=dust_err_state
             call d03_cell_weights(sub_bins(:,i),sub_weights(:,i),k)
             if(k/=0)ierr=dust_err_state
          enddo
          call stage_error(ierr)
          if(ierr/=dust_ok)return
          scatter_sigma=matmul(scattering_basis,sub_weights)
          conductance=2*kb*n_hydrogen*sub_density*sub_area*snrt_dust_contract_accommodation* &
               sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
          if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
               sub_density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
               ghost_arg,remote_arg,blocked,material_dispatch=coupled_material, &
               transport_dispatch=transport_arg,absorb_dispatch=snrt_runtime_ir_absorb, &
               gas_energy=gas_work,gas_capacity=gas_capacity,conductance=conductance,gas_transfer=exchange, &
               cell_material_u=sub_u,cell_weights=sub_weights,thin_reabsorption=.true.,material_only=local_material)
          if(ierr==dust_ok)then
             ! The accepted IR ledger includes phase energy, which is derived
             ! from lost solids rather than stored in the sensible-energy field.
             material=sub_sensible
             sub_bins=sub_trial
             exchange_sum=exchange_sum+exchange
          endif
       else if(cosmo.and.present(gas_energy))then
       conductance=2*kb*n_hydrogen*density*snrt_dust_contract_collision_area_per_h* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(present(cell_collision_area))conductance=2*kb*n_hydrogen*density*cell_collision_area* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghost_arg,remote_arg,blocked,bath_dispatch=cosmological_material, &
            transport_dispatch=transport_arg,absorb_dispatch=snrt_runtime_ir_absorb, &
            gas_energy=gas_work,gas_capacity=gas_capacity,conductance=conductance,gas_transfer=exchange, &
            cell_material_u=cell_material_u,cell_weights=cell_weights,material_only=local_material)
       if(ierr==dust_ok)exchange_sum=exchange_sum+exchange
       else if(present(gas_energy))then
       conductance=2*kb*n_hydrogen*density*snrt_dust_contract_collision_area_per_h* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(present(cell_collision_area))conductance=2*kb*n_hydrogen*density*cell_collision_area* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghost_arg,remote_arg,blocked,material_dispatch=fixed_material, &
            transport_dispatch=transport_arg,absorb_dispatch=snrt_runtime_ir_absorb, &
            gas_energy=gas_work,gas_capacity=gas_capacity,conductance=conductance,gas_transfer=exchange, &
            cell_material_u=cell_material_u,cell_weights=cell_weights,phase_density=phase_density, &
            phase_momentum=phase_p_stage,phase_absorption=phase_alpha,phase_scattering=phase_sca,phase_work=phase_work_step, &
            material_only=local_material)
       if(ierr==dust_ok)exchange_sum=exchange_sum+exchange
       else
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghost_arg,remote_arg,blocked,material_dispatch=fixed_material, &
            transport_dispatch=transport_arg,absorb_dispatch=snrt_runtime_ir_absorb, &
            cell_material_u=cell_material_u,cell_weights=cell_weights,phase_density=phase_density, &
            phase_momentum=phase_p_stage,phase_absorption=phase_alpha,phase_scattering=phase_sca,phase_work=phase_work_step, &
            material_only=local_material)
       endif
       solve_wall=solve_wall+omp_get_wtime()-phase_start
       phase_start=omp_get_wtime()
       ! Conservative delta-isotropic angular relaxation, Lie-split after
       ! absorption/emission. No scattering energy is given to the material.
       if(ierr==dust_ok.and.size(slots)>0.and.allocated(scatter_sigma).and..not.present(phase_density))then
          if(allocated(sub_bins))then
             call snrt_runtime_ir_scatter(trial,weights,sub_density,scatter_sigma,chat*step_dt,ierr)
          else
             call snrt_runtime_ir_scatter(trial,weights,density,scatter_sigma,chat*step_dt,ierr)
          endif
       endif
       scatter_wall=scatter_wall+omp_get_wtime()-phase_start
       if(ierr/=dust_ok.and.size(slots)>0)then
          write(*,'(A,3I6,A,2ES25.16,A,ES14.5)')' SNRT IR rejected state rank/level/error=',myid,ilevel,ierr, &
               ' material_T_range=',minval(material/capacity),maxval(material/capacity), &
               ' min_IR=',minval(trial)
          if(nghost>0)write(*,'(A,ES14.5)')' SNRT IR rejected min_ghost_IR=',minval(ghosts)
       endif
       if(any(.not.ieee_is_finite(coarse%energy)).or.any(coarse%energy<0))ierr=dust_err_state
       call stage_error(ierr)
       if(ierr/=dust_ok)return
       if(present(phase_density))phase_work_sum=phase_work_sum+phase_work_step
       diagnostics%escaped_erg=diagnostics%escaped_erg+step%escaped_erg
       diagnostics%absorbed_erg=diagnostics%absorbed_erg+step%absorbed_erg
       diagnostics%primary_erg=diagnostics%primary_erg+step%primary_erg
       diagnostics%interface_erg=diagnostics%interface_erg+step%interface_erg
       diagnostics%mechanical_erg=diagnostics%mechanical_erg+step%mechanical_erg
       diagnostics%background_erg=diagnostics%background_erg+step%background_erg
       diagnostics%balance_relative=max(diagnostics%balance_relative,step%balance_relative)
       diagnostics%local_relative=max(diagnostics%local_relative,step%local_relative)
       diagnostics%iterations=diagnostics%iterations+step%iterations
    end do
    ! Empty dust cells carry zero material energy, not a fictitious heat bath.
    if(myid==1.and..not.local_material)write(*,'(A,I0,A,I0,3(A,F12.4))')' SNRT_IR timings level=',ilevel,' substeps=',nsub, &
         ' halo=',halo_wall,' solve=',solve_wall,' scatter=',scatter_wall
    ! Their reported temperature is only a harmless diagnostic placeholder.
    do i=1,size(slots)
       if(density(i)==0)temperature(i)=snrt_dust_contract_ir_background_k
       if(cosmo.and.density(i)==0)temperature(i)=2.727d0/aexp
    end do
    if(present(gas_transfer))gas_transfer=exchange_sum
    if(dust_pah_charged())gas_electrons=electron_work
    if(dust_pah_hydrogenated())gas_atomic_h=h_work
    if(dust_pah_h2_enabled())gas_molecular_h2=h2_work
    if(dust_pah_atomization())then
       gas_atomic_c=c_work;gas_carbon_ion=cp_work
    endif
    if(present(sublimation_next))sublimation_next=sub_bins
    if(present(phase_density))then
       phase_momentum=phase_p_stage;phase_work=phase_work_sum
    endif
  contains
    subroutine stage_error(status)
      integer,intent(inout) :: status
      if(.not.local_material)call collective_error(status)
    end subroutine

    subroutine cosmological_material(heating,density,old_energy,log_t,basis_band,cell_weights,dt, &
         gas_energy,gas_capacity,conductance,rate,temperature,next_energy,gas_transfer,background_transfer,ierr)
      ! Optically thin analytic CMB: only positive excess enters transport.
      ! Negative net emission is supplied by the explicitly recorded bath.
      real(dust_dp),intent(in)::heating(:),density(:),old_energy(:),log_t(:),basis_band(:,:,:),cell_weights(:,:),dt
      real(dust_dp),intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
      real(dust_dp),intent(out)::rate(:,:),temperature(:),next_energy(:),gas_transfer(:),background_transfer(:)
      integer,intent(out)::ierr
      real(dust_dp)::bins(6),phase(4),bath
      real(dust_dp),allocatable::bands(:,:,:)
      integer::cell,status,bad
      ierr=dust_err_config
      if(.not.ieee_is_finite(aexp).or.aexp<=0)return
      if(size(cell_weights,1)/=4.or.size(basis_band,3)/=4)return
      if(present(phase_density).or.dust_iron_enabled().or.dust_pah_enabled())return
      bath=2.727d0/aexp
      allocate(bands(size(rate,1),size(log_t),6));bands=0;bands(:,:,1:4)=basis_band
      rate=0;temperature=bath;next_energy=old_energy;gas_transfer=0;background_transfer=0
      bad=0
!$omp parallel do private(cell,bins,phase,status) reduction(max:bad) schedule(static)
      do cell=1,size(density)
         bins=0
         bins(1:4)=cell_weights(:,cell)*density(cell)*snrt_dust_contract_mass_per_h_g
         call iron_radiative_cell(exp(log_t),bath,bins,old_energy(cell),heating(cell),dt, &
              snrt_dust_contract_mass_per_h_g,bands,gas_energy(cell),gas_capacity(cell),conductance(cell), &
              next_energy(cell),temperature(cell),phase,rate(:,cell),gas_transfer(cell),status, &
              photon_ev=snrt_dust_contract_ir_energy_ev(1:size(rate,1)),cmb_exchange=background_transfer(cell))
         bad=max(bad,status)
      enddo
!$omp end parallel do
      ierr=bad
    end subroutine

    subroutine fixed_material(heating,density,old_energy,capacity,log_t,power,band, &
         material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
         gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
      real(dust_dp),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
      real(dust_dp),intent(in)::material_u(:),dt,background,bath,tolerance
      logical,intent(in)::use_u
      real(dust_dp),intent(out)::rate(:,:),temperature(:),next_energy(:)
      integer,intent(out)::ierr
      real(dust_dp),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
      real(dust_dp),optional,intent(out)::gas_transfer(:)
      real(dust_dp),optional,intent(in)::cell_material_u(:,:),cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
      real(dust_dp),allocatable::bins(:,:),phase(:,:)
      real(dust_dp),allocatable::moving_band(:,:,:)
      real(dust_dp)::eg,cv,kappa,qgas
      integer::cell,nb
      if(present(phase_density))then
         ! Moving radiation represents ABSOLUTE energy. Use the existing
         ! absolute C/silicate/Fe receiver, including when Fe masses are zero;
         ! no untracked isotropic bath is subtracted from moving emission.
         ierr=dust_err_config
         if(.not.present(cell_weights).or..not.present(basis_band))return
         nb=size(cell_weights,1)
         if(nb/=4.and.nb/=6)return
         allocate(bins(6,size(density)),phase(4,size(density)),moving_band(size(rate,1),size(log_t),6))
         bins=0;moving_band=0
         bins(1:nb,:)=cell_weights*spread(density*snrt_dust_contract_mass_per_h_g,1,nb)
         moving_band(:,:,1:nb)=basis_band
         rate=0;next_energy=old_energy;temperature=bath;phase=0
         if(present(gas_transfer))gas_transfer=0
         do cell=1,size(density)
            eg=0;cv=1;kappa=0
            if(present(gas_energy))then
               eg=gas_energy(cell);cv=gas_capacity(cell);kappa=conductance(cell)
            endif
            ! Absolute IR has no implicit bath floor. The existing analytic
            ! cold branch requires the first knot and physical band energies.
            call iron_radiative_cell(snrt_dust_contract_temperature_k(1:size(log_t)), &
                 snrt_dust_contract_temperature_k(1),bins(:,cell), &
                 old_energy(cell),heating(cell),dt,snrt_dust_contract_mass_per_h_g,moving_band,eg,cv,kappa, &
                 next_energy(cell),temperature(cell),phase(:,cell),rate(:,cell),qgas,ierr,absolute_emission=.true., &
                 photon_ev=snrt_dust_contract_ir_energy_ev(1:size(rate,1)))
            if(ierr/=0)return
            if(present(gas_transfer))gas_transfer(cell)=qgas
         enddo
         if(dust_iron_enabled().and.any(temperature>dust_fe_max_temperature))ierr=dust_err_range
         return
      endif
      if(.not.dust_iron_enabled())then
         call snrt_runtime_dust_material(heating,density,old_energy,capacity,log_t,power,band, &
              material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
              gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
         return
      endif
      ierr=dust_err_config
      if(.not.present(cell_weights).or..not.present(basis_band))return
      bins=cell_weights*spread(density*snrt_dust_contract_mass_per_h_g,1,6)
      allocate(phase(4,size(density)))
      rate=0;temperature=bath;next_energy=old_energy;phase=0
      if(present(gas_transfer))gas_transfer=0
      call iron_radiative_batch(snrt_dust_contract_temperature_k(1:size(log_t)),bath,bins,old_energy,heating, &
           dt,snrt_dust_contract_mass_per_h_g,basis_band,omp_get_max_threads(), &
           next_energy,temperature,phase,rate,ierr,gas_energy,gas_capacity,conductance,gas_transfer)
      if(ierr/=0)then
         write(*,'(A,2I6)')'Fe material batch rejected rank/native_status=',myid,ierr
         ierr=dust_err_range
      endif
      ! Rejected outputs are trial-only; the enclosing collective transaction
      ! never commits a result outside the explicitly selected comparison.
      if(ierr==0.and.any(temperature>dust_fe_max_temperature))ierr=dust_err_range
    end subroutine

    subroutine coupled_material(heating,density,old_energy,capacity,log_t,power,band, &
         material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
         gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
      real(dust_dp),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
      real(dust_dp),intent(in)::material_u(:),dt,background,bath,tolerance
      logical,intent(in)::use_u
      real(dust_dp),intent(out)::rate(:,:),temperature(:),next_energy(:)
      integer,intent(out)::ierr
      real(dust_dp),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
      real(dust_dp),optional,intent(out)::gas_transfer(:)
      real(dust_dp),optional,intent(in)::cell_material_u(:,:),cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
      integer::j,status,bad
      real(dust_dp)::ed
      ierr=dust_err_config
      if(.not.present(gas_energy).or..not.present(gas_transfer).or..not.present(basis_band))return
      bad=0
      ! Native CPU/OpenMP for this new nonlinear closure; no claim that the
      ! legacy fixed-composition CUDA material kernel implements sublimation.
!$omp parallel do private(j,status,ed) reduction(max:bad)
      do j=1,size(heating)
         call dust_radiative_sublimation_evolve(snrt_dust_contract_temperature_k(1:size(log_t)), &
              bath,sub_bins(:,j),old_energy(j),heating(j),dt, &
              snrt_dust_contract_mass_per_h_g,basis_band,gas_energy(j),gas_capacity(j),conductance(j), &
              sub_trial(:,j),ed,temperature(j),rate(:,j),gas_transfer(j),sub_phase(j),status)
         next_energy(j)=ed+sub_phase(j)
         sub_sensible(j)=ed
         if(status/=0)then
            bad=dust_err_range
            if(j==1)write(*,'(A,2I6,7ES17.8)')'DUST_SUBLIMATION_ROOT_REJECT rank/status/dt/H/Ed/Eg/Cv/K/rhoD=', &
                 myid,status,dt,heating(j),old_energy(j),gas_energy(j),gas_capacity(j),conductance(j),sum(sub_bins(:,j))
         endif
      enddo
!$omp end parallel do
      ierr=bad
    end subroutine

    subroutine validate_stage()
      ierr=dust_err_config
      if(snrt_state_is_moment().and..not.local_material)return
      if(cosmo)then
         if(.not.ieee_is_finite(aexp).or.aexp<=0)return
         if(2.727d0/aexp>snrt_dust_contract_temperature_k(snrt_dust_contract_number_temperature))return
         if(.not.present(gas_energy).or..not.present(cell_weights).or..not.present(cell_material_u))return
         if(present(phase_density).or.present(sublimation_bins).or.dust_iron_enabled().or.dust_pah_enabled())return
         if(.not.snrt_runtime_cpu_material_allowed())return
      endif
      ierr=dust_err_shape
      if(present(phase_density).neqv.present(phase_momentum))return
      if(present(phase_density).neqv.present(phase_work))return
      if(present(phase_density))then
         if(.not.allocated(absorption_basis).or..not.present(cell_weights))return
         ! This local interaction freezes masses. Sublimation with donor
         ! momentum exchange belongs to the separately staged phase update.
         if(present(sublimation_bins))return
         if(any(shape(phase_density)/=[size(absorption_basis,2)+merge(1,0,dust_pah_enabled()),size(slots)]))return
         if(any(shape(phase_momentum)/=[3,size(phase_density,1),size(slots)]))return
         if(any(shape(phase_work)/=shape(phase_density)))return
         if(any(.not.ieee_is_finite(phase_density)).or.any(phase_density<0))return
         if(any(.not.ieee_is_finite(phase_momentum)))return
      endif
      if(present(primary_pah_heat).neqv.present(primary_pah_captures))return
      if(present(primary_pah_heat))then
         if(.not.dust_pah_enabled())return
         if(any(shape(primary_pah_heat)/=[snrt_dust_contract_number_groups,size(slots)]))return
         if(any(shape(primary_pah_captures)/=shape(primary_pah_heat)))return
         if(any(.not.ieee_is_finite(primary_pah_heat)).or.any(primary_pah_heat<0))return
         if(any(.not.ieee_is_finite(primary_pah_captures)).or.any(primary_pah_captures<0))return
      endif
      if(present(pah_population).neqv.dust_pah_enabled())return
      if(present(primary_spectrum).neqv.dust_pah_enabled())return
      if(dust_pah_enabled())then
         if(dust_pah_atomization())then
            if(.not.present(gas_atomic_c).or..not.present(gas_carbon_ion))return
            if(size(gas_atomic_c)/=size(slots).or.size(gas_carbon_ion)/=size(slots))return
            if(any(.not.ieee_is_finite(gas_atomic_c)).or.any(.not.ieee_is_finite(gas_carbon_ion)))return
            if(any(gas_atomic_c<0).or.any(gas_carbon_ion<0))return
         else if(present(gas_atomic_c).or.present(gas_carbon_ion))then
            return
         endif
         if(dust_pah_charged())then
            if(.not.present(gas_electrons).or..not.present(electron_capacity))return
            if(present(phase_density).or.present(primary_pah_heat))return
            if(size(gas_electrons)/=size(slots))return
            if(any(.not.ieee_is_finite(gas_electrons)).or.any(gas_electrons<0))return
            if(dust_pah_hydrogenated())then
               if(.not.present(gas_atomic_h))return
               if(size(gas_atomic_h)/=size(slots))return
               if(any(.not.ieee_is_finite(gas_atomic_h)).or.any(gas_atomic_h<0))return
            endif
            if(dust_pah_h2_enabled())then
               if(.not.present(gas_molecular_h2))return
               if(size(gas_molecular_h2)/=size(slots))return
               if(any(.not.ieee_is_finite(gas_molecular_h2)).or.any(gas_molecular_h2<0))return
            else if(present(gas_molecular_h2))then
               return
            endif
            if(.not.ieee_is_finite(electron_capacity).or.electron_capacity<=0)return
         endif
         if(.not.present(gas_energy).or..not.present(cell_weights).or..not.present(cell_collision_area))return
         if(present(sublimation_bins))return
         if(any(shape(pah_population)/=[dust_pah_nstate(),size(slots)]))return
         if(any(shape(primary_spectrum)/=[snrt_dust_contract_number_groups,size(slots)]))return
         if(any(.not.ieee_is_finite(pah_population)).or.any(pah_population<0))return
         if(any(.not.ieee_is_finite(primary_spectrum)).or.any(primary_spectrum<0))return
         if(any(abs(sum(primary_spectrum,dim=1)-primary_energy)> &
              1d-10*max(primary_energy,tiny(1d0))))return
      endif
      if(present(sublimation_bins).neqv.present(sublimation_next))return
      if(present(sublimation_bins))then
         if(.not.snrt_runtime_cpu_material_allowed())return
         if(.not.present(gas_energy).or..not.present(cell_weights).or..not.present(cell_material_u))return
         if(any(shape(sublimation_bins)/=[4,size(slots)]).or.any(shape(sublimation_next)/=[4,size(slots)]))return
         if(any(.not.ieee_is_finite(sublimation_bins)).or.any(sublimation_bins<0))return
      endif
      if(present(cell_weights).neqv.dust_optics_enabled())return
      if(present(cell_weights))then
         if(any(shape(cell_weights)/=[merge(6,4,dust_iron_enabled()),size(slots)]))return
         if(any(.not.ieee_is_finite(cell_weights)).or.any(cell_weights<0))return
         if(any(abs(sum(cell_weights,dim=1)-1d0)>1d-12))return
      endif
      if(present(cell_material_u).neqv.present(cell_collision_area))return
      if(present(cell_material_u))then
         if(any(shape(cell_material_u)/=[snrt_dust_contract_number_temperature,size(slots)]))return
         if(size(cell_collision_area)/=size(slots))return
         if(any(.not.ieee_is_finite(cell_collision_area)).or.any(cell_collision_area<0))return
         if(any(.not.ieee_is_finite(cell_material_u)).or.any(cell_material_u<=0))return
      endif
      if(snrt_dust_contract_exchange_enabled.neqv.present(gas_energy))return
      if(present(gas_energy))then
         if(.not.present(gas_capacity).or..not.present(n_hydrogen).or..not.present(gas_transfer))return
         if(size(gas_energy)/=size(slots).or.size(gas_capacity)/=size(slots).or. &
              size(n_hydrogen)/=size(slots).or.size(gas_transfer)/=size(slots))return
         ierr=dust_err_state
         if(any(.not.ieee_is_finite(gas_energy)).or.any(gas_energy<0))return
         if(any(.not.ieee_is_finite(gas_capacity)).or.any(gas_capacity<=0))return
         if(any(.not.ieee_is_finite(n_hydrogen)).or.any(n_hydrogen<0))return
      else if(present(gas_capacity).or.present(n_hydrogen).or.present(gas_transfer))then
         return
      endif
      ierr=dust_err_shape
      if(size(material)/=size(slots).or.size(temperature)/=size(slots).or.size(cells)/=size(slots))return
      if(size(old_energy)/=size(slots).or.size(capacity)/=size(slots))return
      if(size(density)/=size(slots).or.size(primary_energy)/=size(slots))return
      if(nd<1.or.any(shape(directions)/=[3,nd]))return
      if(.not.local_material.and.nd/=snrt_ndirection)return
      if(any(shape(neighbors)/=[6,size(slots)]))return
      if(local_material)then
         if(any(shape(incoming_radiation)/=[snrt_dust_contract_number_ir,nd,size(slots)]))return
         if(any(neighbors/=0))return
         ierr=dust_err_state
         if(any(.not.ieee_is_finite(incoming_radiation)).or.any(incoming_radiation<0))return
      endif
      if(any(slots<1).or.any(slots>snrt_nslot))return
      nfield=ICELL_OF(ngridmax,twotondim)
      if(any(cells<1).or.any(cells>nfield))return
      if(size(slots)>0.and..not.local_material)then
         if(.not.allocated(snrt_face_kind).or..not.allocated(snrt_face_cell))return
         if(any(shape(snrt_face_kind)/=[6,size(slots)]))return
         if(any(shape(snrt_face_cell)/=[6,size(slots)]))return
         ierr=dust_err_config
         if(any(snrt_face_kind/=SNRT_FACE_LOCAL.and.snrt_face_kind/=SNRT_FACE_PHYSICAL.and. &
              snrt_face_kind/=SNRT_FACE_MPI.and.snrt_face_kind/=SNRT_FACE_FINE_TO_COARSE.and. &
              snrt_face_kind/=SNRT_FACE_COARSE_TO_FINE))return
         if(any(snrt_face_kind==SNRT_FACE_FINE_TO_COARSE))then
            if(ilevel<=1)return
            if(.not.allocated(headl).or..not.allocated(next).or..not.allocated(son))return
         end if
         do i=1,size(slots)
            do face=1,6
               if(snrt_face_kind(face,i)==SNRT_FACE_MPI.and.ncpu<2)return
               if(snrt_face_kind(face,i)/=SNRT_FACE_MPI.and.snrt_face_kind(face,i)/=SNRT_FACE_FINE_TO_COARSE)cycle
               if(snrt_face_cell(face,i)<1.or.snrt_face_cell(face,i)>nfield)return
            end do
         end do
      end if
      ierr=dust_err_state
      if(any(capacity<=0).or.any(.not.ieee_is_finite(capacity)))return
      ierr=dust_err_config
      if(.not.all(ieee_is_finite([dx,dt,chat])).or.min(dx,dt,chat)<=0)return
      nsub=1
      if(.not.local_material)then
         cfl=chat*dt/dx*maxval(sum(abs(directions),dim=1))
         if(.not.ieee_is_finite(cfl).or.cfl>real(huge(nsub)-1,dust_dp))return
         nsub=max(1,ceiling(cfl))
      endif
      ierr=dust_ok
    end subroutine
  end subroutine

  subroutine collective_error(ierr)
    integer, intent(inout) :: ierr
#ifndef WITHOUTMPI
    integer :: global_error,info,abort_info
    call MPI_ALLREDUCE(ierr,global_error,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,abort_info)
    ierr=global_error
#endif
  end subroutine

  subroutine snrt_dust_live_commit(slots,trial,coarse)
    integer, intent(in) :: slots(:)
    real(dust_dp), intent(in) :: trial(:,:,:)
    type(dust_live_coarse_trial), intent(in) :: coarse
    integer::i
    ! Called only after the primary transaction commits; no allocation or
    ! fallible conversion remains here. Stage has validated the slot window.
    do i=1,size(slots)
       call ir_store(slots(i),trial(:,:,i))
    enddo
    do i=1,size(coarse%slots)
       call ir_store(coarse%slots(i),coarse%energy(:,:,i))
    enddo
  end subroutine

  subroutine snrt_dust_live_pack(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dust_dp), intent(out) :: payload(:)
    integer, intent(out) :: ierr
    integer :: slot
    payload=0
    ierr=dust_err_shape
    if(size(payload)/=snrt_dust_contract_number_ir*snrt_ndirection)return
    slot=snrt_state_get_slot(icell)
    call prepare(ierr)
    if(ierr/=dust_ok.or.slot==0)return
    if(snrt_state_is_moment())then
       call mn_live_ir_pack(slot,snrt_dust_contract_number_ir,snrt_ndirection,payload,ierr)
       return
    endif
    payload=reshape(ir_slot(slot),[size(payload)])
    if(any(.not.ieee_is_finite(payload)).or.any(payload<0))ierr=dust_err_state
  end subroutine

  subroutine snrt_dust_live_restore(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dust_dp), intent(in) :: payload(:)
    integer, intent(out) :: ierr
    integer :: slot
    call snrt_dust_live_validate_payload(payload,ierr)
    if(ierr/=dust_ok)return
    slot=snrt_state_get_slot(icell)
    if(slot==0)then
       ierr=dust_err_state
       if(all(payload==0))ierr=dust_ok
       return
    end if
    call prepare(ierr)
    if(ierr/=dust_ok)return
    if(snrt_state_is_moment())then
       call mn_live_ir_unpack(slot,snrt_dust_contract_number_ir,snrt_ndirection,payload,ierr)
       return
    endif
    call ir_store(slot,reshape(payload,[snrt_dust_contract_number_ir,snrt_ndirection]))
  end subroutine

  subroutine snrt_dust_live_validate_payload(payload,ierr)
    real(dust_dp),intent(in) :: payload(:)
    integer,intent(out) :: ierr
    ierr=dust_err_shape
    if(size(payload)/=snrt_dust_contract_number_ir*snrt_ndirection)return
    ierr=dust_err_state
    if(any(.not.ieee_is_finite(payload)))return
    if(snrt_state_is_moment())then
       if(all(payload==0))then
          ierr=dust_ok;return
       endif
       call mn_live_ir_unpack(0,snrt_dust_contract_number_ir,snrt_ndirection,payload,ierr,validate_only=.true.)
    else
       if(any(payload<0))return
       ierr=dust_ok
    endif
  end subroutine
end module

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
  use omp_lib, only: omp_get_max_threads
  use dust_composition_material, only: dust_composition_curve,dust_composition_area
  use dust_sublimation_material, only: dust_radiative_sublimation_evolve
  use dust_composition_optics, only: d03_ng,d03_nir,d03_opacity_basis,d03_cell_weights
  use snrt_runtime_backend, only: snrt_runtime_dust_material, snrt_runtime_ir_transport, snrt_runtime_ir_absorb
  use snrt_runtime_backend, only: snrt_runtime_ir_scatter
  use snrt_runtime_backend, only: snrt_runtime_cpu_material_allowed
  use snrt_state, only: snrt_ndirection, snrt_nslot, snrt_state_get_slot
  use amr_commons, only: ngridmax,ncoarse,ncpu,myid,headl,next,son
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
  real(dust_dp), allocatable, save :: radiation(:,:,:)
  real(dust_dp), allocatable, save :: scattering_basis(:,:)
  real(dust_dp), allocatable, save :: absorption_basis(:,:)
  type(dust_ir_table), save :: table
  logical, save :: initialized=.false.
  type, public :: dust_live_coarse_trial
     integer, allocatable :: slots(:)
     real(dust_dp), allocatable :: energy(:,:,:)
  end type
  public :: snrt_dust_live_stage, snrt_dust_live_commit
  public :: snrt_dust_live_pack, snrt_dust_live_restore
contains
  subroutine prepare(ierr)
    integer, intent(out) :: ierr
    integer :: ng, nt, old, capacity,nb
    real(dust_dp), allocatable :: expanded(:,:,:)
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
    old=0
    if(allocated(radiation))then
       if(size(radiation,1)/=ng)then
          ierr=dust_err_config
          return
       end if
       old=size(radiation,3)
    end if
    if(old<snrt_nslot)then
       capacity=max(snrt_nslot,max(16,2*old))
       allocate(expanded(ng,snrt_ndirection,capacity)); expanded=0
       if(old>0)expanded(:,:,1:old)=radiation
       call move_alloc(expanded,radiation)
    end if
    ierr=dust_ok
  end subroutine

  subroutine snrt_dust_live_stage(ilevel,cells,slots,neighbors,directions,weights,dx,dt,chat, &
       density,primary_energy,old_energy,capacity,trial,material,temperature,diagnostics,ierr,coarse, &
       gas_energy,gas_capacity,n_hydrogen,gas_transfer,cell_material_u,cell_collision_area,cell_weights, &
       sublimation_bins,sublimation_next,pah_population,primary_spectrum,primary_pah_heat,primary_pah_captures, &
       phase_density,phase_momentum,phase_work,gas_electrons,electron_capacity,gas_atomic_h,gas_molecular_h2, &
       gas_atomic_c,gas_carbon_ion)
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
    real(dust_dp), allocatable :: photons(:,:),ghosts(:,:,:),field(:)
    real(dust_dp), allocatable :: halo_field(:,:)
    integer,parameter :: halo_tile=16
    integer :: first_component,ncomponent,component,column
    integer, allocatable :: remote(:,:),ghost_cells(:),ghost_kind(:),coarse_cells(:)
    logical, allocatable :: blocked(:,:)
    type(dust_ir_diagnostics) :: step
    real(dust_dp) :: cfl,step_dt,mu,q
    integer :: nsub,isub,i,ng,face,k,nghost,nfield,g,d,global_nsub,info,has_coarse,b,nbulk
    integer :: grid,child,cell,ncoarse_leaf,axis
    call prepare(ierr)
    if(ierr==dust_ok)call validate_stage()
    call collective_error(ierr)
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
    call MPI_ALLREDUCE(nsub,global_nsub,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,k)
    nsub=global_nsub
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
    if(size(slots)>0)then
       if(any(snrt_face_kind==SNRT_FACE_FINE_TO_COARSE))has_coarse=1
    end if
    ! Use the same global decision on empty ranks and on coarse-only owners.
    call collective_error(has_coarse)
    ierr=dust_ok
    if(has_coarse/=0)then
       if(ilevel<=1)ierr=dust_err_config
       if(.not.allocated(headl).or..not.allocated(next).or..not.allocated(son))ierr=dust_err_config
    end if
    call collective_error(ierr)
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
    allocate(coarse%energy(ng,snrt_ndirection,ncoarse_leaf))
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
    call collective_error(ierr)
    if(ierr/=dust_ok)return
    if(ncoarse_leaf>0)coarse%energy=radiation(:,:,coarse%slots)
    allocate(trial(ng,snrt_ndirection,size(slots)),photons(ng,size(slots)))
    if(size(slots)>0)trial=radiation(:,:,slots)
    photons=0
    material=old_energy
    temperature=material/capacity
    nghost=0
    allocate(remote(6,size(slots)),blocked(6,size(slots))); remote=0; blocked=.false.
    if(size(slots)>0)then
       nghost=count(snrt_face_kind==SNRT_FACE_MPI.or.snrt_face_kind==SNRT_FACE_FINE_TO_COARSE)
       blocked=snrt_face_kind==SNRT_FACE_COARSE_TO_FINE
    end if
    allocate(ghost_cells(nghost),ghost_kind(nghost),ghosts(ng,snrt_ndirection,nghost),field(nfield))
    if(ncpu>1)allocate(halo_field(nfield,halo_tile))
    k=0
    do i=1,size(slots)
       do face=1,6
          if(snrt_face_kind(face,i)/=SNRT_FACE_MPI.and.snrt_face_kind(face,i)/=SNRT_FACE_FINE_TO_COARSE)cycle
          k=k+1; remote(face,i)=k; ghost_cells(k)=snrt_face_cell(face,i)
          ghost_kind(k)=snrt_face_kind(face,i)
       end do
    end do
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
       call collective_error(ierr)
       if(ierr/=dust_ok)return
    end if
    do isub=1,nsub
       ! Reuse RAMSES' real halo communicator, in FP64. All ranks participate
       ! in every substep, including a rank with no local MPI boundary faces.
       if(ncpu>1)then
          do first_component=1,ng*snrt_ndirection,halo_tile
             ncomponent=min(halo_tile,ng*snrt_ndirection-first_component+1)
             halo_field=0d0
             do column=1,ncomponent
                component=first_component+column-1
                g=mod(component-1,ng)+1;d=(component-1)/ng+1
                halo_field(cells,column)=trial(g,d,:)
             enddo
             call snrt_halo_tile_exchange(halo_field(:,1:ncomponent),ilevel,ierr)
             call collective_error(ierr)
             if(ierr/=dust_ok)return
             do column=1,ncomponent
                component=first_component+column-1
                g=mod(component-1,ng)+1;d=(component-1)/ng+1
                do k=1,nghost
                   if(ghost_kind(k)==SNRT_FACE_MPI)ghosts(g,d,k)=halo_field(ghost_cells(k),column)
                end do
             end do
          end do
       end if
       if(has_coarse/=0)then
          do d=1,snrt_ndirection
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
               ghosts,remote,blocked,gas_work,pah_capacity,conductance,exchange, &
               primary_pah_heat=pah_heat_step,primary_pah_captures=pah_captures_step, &
               phase_density=phase_density,phase_momentum=phase_p_stage,phase_absorption=phase_alpha, &
               phase_scattering=phase_sca,phase_work=phase_work_step, &
               gas_electrons=electron_work,electron_capacity=electron_capacity,primary_population=primary_pah_population, &
               gas_atomic_h=h_work,gas_molecular_h2=h2_work,gas_atomic_c=c_work,gas_carbon_ion=cp_work)
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
          call collective_error(ierr)
          if(ierr/=dust_ok)return
          scatter_sigma=matmul(scattering_basis,sub_weights)
          conductance=2*kb*n_hydrogen*sub_density*sub_area*snrt_dust_contract_accommodation* &
               sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
          if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
               sub_density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
               ghosts,remote,blocked,material_dispatch=coupled_material, &
               transport_dispatch=snrt_runtime_ir_transport,absorb_dispatch=snrt_runtime_ir_absorb, &
               gas_energy=gas_work,gas_capacity=gas_capacity,conductance=conductance,gas_transfer=exchange, &
               cell_material_u=sub_u,cell_weights=sub_weights,thin_reabsorption=.true.)
          if(ierr==dust_ok)then
             ! The accepted IR ledger includes phase energy, which is derived
             ! from lost solids rather than stored in the sensible-energy field.
             material=sub_sensible
             sub_bins=sub_trial
             exchange_sum=exchange_sum+exchange
          endif
       else if(present(gas_energy))then
       conductance=2*kb*n_hydrogen*density*snrt_dust_contract_collision_area_per_h* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(present(cell_collision_area))conductance=2*kb*n_hydrogen*density*cell_collision_area* &
            snrt_dust_contract_accommodation*sqrt((8*kb/(acos(-1d0)*mp))*(gas_work/gas_capacity))
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghosts,remote,blocked,material_dispatch=fixed_material, &
            transport_dispatch=snrt_runtime_ir_transport,absorb_dispatch=snrt_runtime_ir_absorb, &
            gas_energy=gas_work,gas_capacity=gas_capacity,conductance=conductance,gas_transfer=exchange, &
            cell_material_u=cell_material_u,cell_weights=cell_weights,phase_density=phase_density, &
            phase_momentum=phase_p_stage,phase_absorption=phase_alpha,phase_scattering=phase_sca,phase_work=phase_work_step)
       if(ierr==dust_ok)exchange_sum=exchange_sum+exchange
       else
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghosts,remote,blocked,material_dispatch=fixed_material, &
            transport_dispatch=snrt_runtime_ir_transport,absorb_dispatch=snrt_runtime_ir_absorb, &
            cell_material_u=cell_material_u,cell_weights=cell_weights,phase_density=phase_density, &
            phase_momentum=phase_p_stage,phase_absorption=phase_alpha,phase_scattering=phase_sca,phase_work=phase_work_step)
       endif
       ! Conservative delta-isotropic angular relaxation, Lie-split after
       ! absorption/emission. No scattering energy is given to the material.
       if(ierr==dust_ok.and.size(slots)>0.and.allocated(scatter_sigma).and..not.present(phase_density))then
          if(allocated(sub_bins))then
             call snrt_runtime_ir_scatter(trial,weights,sub_density,scatter_sigma,chat*step_dt,ierr)
          else
             call snrt_runtime_ir_scatter(trial,weights,density,scatter_sigma,chat*step_dt,ierr)
          endif
       endif
       if(ierr/=dust_ok.and.size(slots)>0)then
          write(*,'(A,3I6,A,2ES25.16,A,ES14.5)')' SNRT IR rejected state rank/level/error=',myid,ilevel,ierr, &
               ' material_T_range=',minval(material/capacity),maxval(material/capacity), &
               ' min_IR=',minval(trial)
          if(nghost>0)write(*,'(A,ES14.5)')' SNRT IR rejected min_ghost_IR=',minval(ghosts)
       endif
       if(any(.not.ieee_is_finite(coarse%energy)).or.any(coarse%energy<0))ierr=dust_err_state
       call collective_error(ierr)
       if(ierr/=dust_ok)return
       if(present(phase_density))phase_work_sum=phase_work_sum+phase_work_step
       diagnostics%escaped_erg=diagnostics%escaped_erg+step%escaped_erg
       diagnostics%absorbed_erg=diagnostics%absorbed_erg+step%absorbed_erg
       diagnostics%primary_erg=diagnostics%primary_erg+step%primary_erg
       diagnostics%interface_erg=diagnostics%interface_erg+step%interface_erg
       diagnostics%mechanical_erg=diagnostics%mechanical_erg+step%mechanical_erg
       diagnostics%balance_relative=max(diagnostics%balance_relative,step%balance_relative)
       diagnostics%local_relative=max(diagnostics%local_relative,step%local_relative)
       diagnostics%iterations=diagnostics%iterations+step%iterations
    end do
    ! Empty dust cells carry zero material energy, not a fictitious heat bath.
    ! Their reported temperature is only a harmless diagnostic placeholder.
    do i=1,size(slots)
       if(density(i)==0)temperature(i)=snrt_dust_contract_ir_background_k
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
      if(any(shape(directions)/=[3,snrt_ndirection]).or.size(weights)/=snrt_ndirection)return
      if(any(shape(neighbors)/=[6,size(slots)]))return
      if(any(slots<1).or.any(slots>snrt_nslot))return
      nfield=ICELL_OF(ngridmax,twotondim)
      if(any(cells<1).or.any(cells>nfield))return
      if(size(slots)>0)then
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
      cfl=chat*dt/dx*maxval(sum(abs(directions),dim=1))
      if(.not.ieee_is_finite(cfl).or.cfl>real(huge(nsub)-1,dust_dp))return
      nsub=max(1,ceiling(cfl))
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
    ! Called only after the primary transaction commits; no allocation or
    ! fallible conversion remains here. Stage has validated the slot window.
    if(size(slots)>0)radiation(:,:,slots)=trial
    if(size(coarse%slots)>0)radiation(:,:,coarse%slots)=coarse%energy
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
    payload=reshape(radiation(:,:,slot),[size(payload)])
    if(any(.not.ieee_is_finite(payload)).or.any(payload<0))ierr=dust_err_state
  end subroutine

  subroutine snrt_dust_live_restore(icell,payload,ierr)
    integer, intent(in) :: icell
    real(dust_dp), intent(in) :: payload(:)
    integer, intent(out) :: ierr
    integer :: slot
    ierr=dust_err_shape
    if(size(payload)/=snrt_dust_contract_number_ir*snrt_ndirection)return
    ierr=dust_err_state
    if(any(.not.ieee_is_finite(payload)).or.any(payload<0))return
    slot=snrt_state_get_slot(icell)
    if(slot==0)then
       if(all(payload==0))ierr=dust_ok
       return
    end if
    call prepare(ierr)
    if(ierr/=dust_ok)return
    radiation(:,:,slot)=reshape(payload,[snrt_dust_contract_number_ir,snrt_ndirection])
  end subroutine
end module

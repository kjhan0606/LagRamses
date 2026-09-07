! Persistent IR energy uses the primary RT cell-slot map, but a separate
! spectral axis and FP64 physical units (erg/cm3 per normalized direction).
! Trial work never writes persistent radiation or RAMSES material fields.
module snrt_dust_live
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_dust_contract
  use snrt_dust_ir
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
    integer :: ng, nt, old, capacity
    real(dust_dp), allocatable :: expanded(:,:,:)
    ierr=dust_err_config
    if(snrt_dust_contract_version/=3.or..not.snrt_dust_contract_runtime_allowed)return
    ng=snrt_dust_contract_number_ir; nt=snrt_dust_contract_number_temperature
    if(.not.initialized)then
       call snrt_dust_ir_initialize(table,snrt_dust_contract_ir_energy_ev(1:ng), &
            snrt_dust_contract_ir_weight_ev(1:ng),snrt_dust_contract_ir_absorption_per_h_cm2(1:ng), &
            snrt_dust_contract_temperature_k(1:nt),snrt_dust_contract_ir_background_k,ierr)
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
       density,primary_energy,old_energy,capacity,trial,material,temperature,diagnostics,ierr,coarse)
    integer, intent(in) :: ilevel,cells(:),slots(:),neighbors(:,:)
    real(dust_dp), intent(in) :: directions(:,:),weights(:),dx,dt,chat
    real(dust_dp), intent(in) :: density(:),primary_energy(:),old_energy(:),capacity(:)
    real(dust_dp), allocatable, intent(out) :: trial(:,:,:)
    real(dust_dp), intent(out) :: material(:),temperature(:)
    type(dust_ir_diagnostics), intent(out) :: diagnostics
    integer, intent(out) :: ierr
    type(dust_live_coarse_trial), intent(out) :: coarse
    real(dust_dp), allocatable :: photons(:,:),ghosts(:,:,:),field(:)
    real(dust_dp), allocatable :: halo_field(:,:)
    integer,parameter :: halo_tile=16
    integer :: first_component,ncomponent,component,column
    integer, allocatable :: remote(:,:),ghost_cells(:),ghost_kind(:),coarse_cells(:)
    logical, allocatable :: blocked(:,:)
    type(dust_ir_diagnostics) :: step
    real(dust_dp) :: cfl,step_dt,mu,q
    integer :: nsub,isub,i,ng,face,k,nghost,nfield,g,d,global_nsub,info,has_coarse
    integer :: grid,child,cell,ncoarse_leaf,axis
    call prepare(ierr)
    if(ierr==dust_ok)call validate_stage()
    call collective_error(ierr)
    if(ierr/=dust_ok)return
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(nsub,global_nsub,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,k)
    nsub=global_nsub
#endif
    step_dt=dt/nsub
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
       if(size(slots)>0)call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity, &
            ghosts,remote,blocked)
       if(ierr/=dust_ok.and.size(slots)>0)then
          write(*,'(A,3I6,A,2ES25.16,A,ES14.5)')' SNRT IR rejected state rank/level/error=',myid,ilevel,ierr, &
               ' material_T_range=',minval(material/capacity),maxval(material/capacity), &
               ' min_IR=',minval(trial)
          if(nghost>0)write(*,'(A,ES14.5)')' SNRT IR rejected min_ghost_IR=',minval(ghosts)
       endif
       if(any(.not.ieee_is_finite(coarse%energy)).or.any(coarse%energy<0))ierr=dust_err_state
       call collective_error(ierr)
       if(ierr/=dust_ok)return
       diagnostics%escaped_erg=diagnostics%escaped_erg+step%escaped_erg
       diagnostics%absorbed_erg=diagnostics%absorbed_erg+step%absorbed_erg
       diagnostics%primary_erg=diagnostics%primary_erg+step%primary_erg
       diagnostics%interface_erg=diagnostics%interface_erg+step%interface_erg
       diagnostics%balance_relative=max(diagnostics%balance_relative,step%balance_relative)
       diagnostics%local_relative=max(diagnostics%local_relative,step%local_relative)
       diagnostics%iterations=diagnostics%iterations+step%iterations
    end do
    ! Empty dust cells carry zero material energy, not a fictitious heat bath.
    ! Their reported temperature is only a harmless diagnostic placeholder.
    do i=1,size(slots)
       if(density(i)==0)temperature(i)=snrt_dust_contract_ir_background_k
    end do
  contains
    subroutine validate_stage()
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

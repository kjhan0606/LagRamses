! Bounded one-cell representation adapter to the existing IR matter solver.
! Spatial transport is performed separately on compact moments. No persistent
! angular field, no replacement grain physics, no physical-state commit here.
module snrt_moment_dust_ir
  use snrt_moment_transport, only: mn_basis,mn_dp,mn_ok,mn_bad_input,mn_reconstruct,mn_project
  use snrt_dust_ir, only: dust_ir_table,dust_ir_diagnostics,snrt_dust_ir_advance, &
       dust_material_dispatch,dust_population_dispatch,dust_moving_population_dispatch, &
       dust_bath_dispatch,dust_absorb_dispatch
  implicit none
  private
  public :: mn_ir_material_cell
contains
  subroutine mn_ir_material_cell(b,table,dx,dt,c_hat,density,primary,energy,temperature,photons, &
       diagnostics,projection,ierr,tolerance,max_iterations,dust_energy,heat_capacity, &
       material_dispatch,absorb_dispatch,gas_energy,gas_capacity,conductance,gas_transfer, &
       cell_material_u,cell_weights,thin_reabsorption,population,population_dispatch, &
       cell_absorption,phase_density,phase_momentum,phase_absorption,phase_scattering, &
       phase_work,moving_material_dispatch,population_loss,bath_dispatch)
    type(mn_basis),intent(in) :: b
    type(dust_ir_table),intent(in) :: table
    real(mn_dp),intent(in) :: dx,dt,c_hat,density,primary,tolerance
    integer,intent(in) :: max_iterations
    real(mn_dp),intent(inout) :: energy(:,:),temperature,photons(:),projection(:,:)
    type(dust_ir_diagnostics),intent(inout) :: diagnostics
    integer,intent(out) :: ierr
    ! Optional shapes are the existing native IR shapes with ncell=1. Keeping
    ! that ABI permits existing bulk/PAH/moving/bath callbacks without copies
    ! of their scientific implementation. Trial writable arrays ARE copied:
    ! final moment closure failure must also roll back the material callback.
    real(mn_dp),optional,intent(inout) :: dust_energy(:),gas_energy(:),gas_transfer(:),population(:,:)
    real(mn_dp),optional,intent(inout) :: phase_momentum(:,:,:),phase_work(:,:)
    real(mn_dp),optional,intent(in) :: heat_capacity(:),gas_capacity(:),conductance(:)
    real(mn_dp),optional,intent(in) :: cell_material_u(:,:),cell_weights(:,:),cell_absorption(:,:)
    real(mn_dp),optional,intent(in) :: phase_density(:,:),phase_absorption(:,:,:),phase_scattering(:,:,:)
    logical,optional,intent(in) :: thin_reabsorption,population_loss
    procedure(dust_material_dispatch),optional :: material_dispatch
    procedure(dust_absorb_dispatch),optional :: absorb_dispatch
    procedure(dust_population_dispatch),optional :: population_dispatch
    procedure(dust_moving_population_dispatch),optional :: moving_material_dispatch
    procedure(dust_bath_dispatch),optional :: bath_dispatch
    real(mn_dp) :: angular(size(energy,2),b%nq,1),check(b%nq)
    real(mn_dp) :: candidate(size(energy,1),size(energy,2)),receipt(size(energy,1),size(energy,2))
    real(mn_dp) :: next_t(1),emitted(size(photons),1)
    real(mn_dp),allocatable :: next_dust(:),next_gas(:),next_transfer(:),next_population(:,:)
    real(mn_dp),allocatable :: next_momentum(:,:,:),next_work(:,:)
    type(dust_ir_diagnostics) :: trial
    integer :: ng,g,links(6,1),status
    ierr=mn_bad_input;ng=size(energy,2)
    if(b%nm<1.or.size(energy,1)/=b%nm.or.ng<1.or.size(photons)/=ng)return
    if(any(shape(projection)/=shape(energy)))return
    do g=1,ng
       call mn_reconstruct(b,energy(:,g),angular(g,:,1),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
       call mn_project(b,angular(g,:,1),candidate(:,g),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
    enddo
    receipt=candidate-energy
    next_t=temperature;emitted(:,1)=photons;links=0;trial=diagnostics
    if(present(dust_energy))allocate(next_dust,source=dust_energy)
    if(present(gas_energy))allocate(next_gas,source=gas_energy)
    if(present(gas_transfer))allocate(next_transfer,source=gas_transfer)
    if(present(population))allocate(next_population,source=population)
    if(present(phase_momentum))allocate(next_momentum,source=phase_momentum)
    if(present(phase_work))allocate(next_work,source=phase_work)
    call snrt_dust_ir_advance(table,b%direction,b%weight,links,dx,dt,c_hat,[density],[primary], &
         angular,next_t,emitted,trial,status,tolerance,max_iterations, &
         dust_energy=next_dust,heat_capacity=heat_capacity,material_dispatch=material_dispatch, &
         absorb_dispatch=absorb_dispatch,gas_energy=next_gas,gas_capacity=gas_capacity, &
         conductance=conductance,gas_transfer=next_transfer,cell_material_u=cell_material_u, &
         cell_weights=cell_weights,thin_reabsorption=thin_reabsorption,population=next_population, &
         population_dispatch=population_dispatch,cell_absorption=cell_absorption, &
         phase_density=phase_density,phase_momentum=next_momentum,phase_absorption=phase_absorption, &
         phase_scattering=phase_scattering,phase_work=next_work,moving_material_dispatch=moving_material_dispatch, &
         population_loss=population_loss,bath_dispatch=bath_dispatch,material_only=.true.)
    if(status/=0)then
       ierr=200+status;return
    endif
    do g=1,ng
       call mn_project(b,angular(g,:,1),candidate(:,g),status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
       call mn_reconstruct(b,candidate(:,g),check,status)
       if(status/=mn_ok)then
          ierr=100+status;return
       endif
    enddo
    ! Publish only after every group, field and callback trial has succeeded.
    energy=candidate;temperature=next_t(1);photons=emitted(:,1);diagnostics=trial;projection=receipt
    if(present(dust_energy))dust_energy=next_dust
    if(present(gas_energy))gas_energy=next_gas
    if(present(gas_transfer))gas_transfer=next_transfer
    if(present(population))population=next_population
    if(present(phase_momentum))phase_momentum=next_momentum
    if(present(phase_work))phase_work=next_work
    ierr=mn_ok
  end subroutine
end module

! Persistent IR energy uses the primary RT cell-slot map, but a separate
! spectral axis and FP64 physical units (erg/cm3 per normalized direction).
! Trial work never writes persistent radiation or RAMSES material fields.
module snrt_dust_live
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_dust_contract
  use snrt_dust_ir
  use snrt_state, only: snrt_ndirection, snrt_nslot, snrt_state_get_slot
  implicit none
  private
  real(dust_dp), allocatable, save :: radiation(:,:,:)
  type(dust_ir_table), save :: table
  logical, save :: initialized=.false.
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

  subroutine snrt_dust_live_stage(slots,neighbors,directions,weights,dx,dt,chat, &
       density,primary_energy,old_energy,capacity,trial,material,temperature,diagnostics,ierr)
    integer, intent(in) :: slots(:),neighbors(:,:)
    real(dust_dp), intent(in) :: directions(:,:),weights(:),dx,dt,chat
    real(dust_dp), intent(in) :: density(:),primary_energy(:),old_energy(:),capacity(:)
    real(dust_dp), allocatable, intent(out) :: trial(:,:,:)
    real(dust_dp), intent(out) :: material(:),temperature(:)
    type(dust_ir_diagnostics), intent(out) :: diagnostics
    integer, intent(out) :: ierr
    real(dust_dp), allocatable :: photons(:,:)
    type(dust_ir_diagnostics) :: step
    real(dust_dp) :: cfl,step_dt
    integer :: nsub,isub,i,ng
    call prepare(ierr)
    if(ierr/=dust_ok)return
    ierr=dust_err_shape
    if(size(slots)<1.or.size(material)/=size(slots).or.size(temperature)/=size(slots))return
    if(size(old_energy)/=size(slots).or.size(capacity)/=size(slots))return
    if(size(density)/=size(slots).or.size(primary_energy)/=size(slots))return
    if(any(shape(directions)/=[3,snrt_ndirection]).or.size(weights)/=snrt_ndirection)return
    if(any(shape(neighbors)/=[6,size(slots)]))return
    if(any(slots<1).or.any(slots>snrt_nslot))return
    ierr=dust_err_config
    if(.not.all(ieee_is_finite([dx,dt,chat])).or.min(dx,dt,chat)<=0)return
    cfl=chat*dt/dx*maxval(sum(abs(directions),dim=1))
    if(.not.ieee_is_finite(cfl).or.cfl>real(huge(nsub)-1,dust_dp))return
    nsub=max(1,ceiling(cfl)); step_dt=dt/nsub
    ng=snrt_dust_contract_number_ir
    allocate(trial(ng,snrt_ndirection,size(slots)),photons(ng,size(slots)))
    trial=radiation(:,:,slots); photons=0
    material=old_energy
    ierr=dust_err_state
    if(any(capacity<=0).or.any(.not.ieee_is_finite(capacity)))return
    temperature=material/capacity
    diagnostics=dust_ir_diagnostics()
    do isub=1,nsub
       call snrt_dust_ir_advance(table,directions,weights,neighbors,dx,step_dt,chat, &
            density,primary_energy/dt,trial,temperature,photons,step,ierr,1d-9,256,material,capacity)
       if(ierr/=dust_ok)return
       diagnostics%escaped_erg=diagnostics%escaped_erg+step%escaped_erg
       diagnostics%absorbed_erg=diagnostics%absorbed_erg+step%absorbed_erg
       diagnostics%primary_erg=diagnostics%primary_erg+step%primary_erg
       diagnostics%balance_relative=max(diagnostics%balance_relative,step%balance_relative)
       diagnostics%local_relative=max(diagnostics%local_relative,step%local_relative)
       diagnostics%iterations=diagnostics%iterations+step%iterations
    end do
    ! Empty dust cells carry zero material energy, not a fictitious heat bath.
    ! Their reported temperature is only a harmless diagnostic placeholder.
    do i=1,size(slots)
       if(density(i)==0)temperature(i)=snrt_dust_contract_ir_background_k
    end do
  end subroutine

  subroutine snrt_dust_live_commit(slots,trial)
    integer, intent(in) :: slots(:)
    real(dust_dp), intent(in) :: trial(:,:,:)
    ! Called only after the primary transaction commits; no allocation or
    ! fallible conversion remains here. Stage has validated the slot window.
    radiation(:,:,slots)=trial
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

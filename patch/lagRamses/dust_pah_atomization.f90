! Opt-in single-photon C24 complete-atomization yield comparison primitive.
! No opacity, rate fit, ensemble energy pool, source heating or live selector.
! Provenance: pah_single_photon_atomization_plan_2026-09-11.md.
module dust_pah_atomization_physics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  integer,parameter,public::pah_atomization_ok=0,pah_atomization_invalid=1,pah_atomization_exhausted=2
  real(real64),parameter,public::pah_atomization_ev=1.602176634d-12
  real(real64),parameter::kjmol_to_erg=1d10/6.02214076d23
  ! NIST/Roux 2008 gas Hf298=295 +/-11 kJ/mol. Convert with native DL01
  ! C24H12 modes: H298-H0=Uvib(298.15)+4kT=39.51733815680014 kJ/mol.
  ! CODATA graphite/H2 increments=1.050/8.468 kJ/mol. This thermal-model
  ! conversion is NOT a measured 0 K coronene atomization enthalpy.
  real(real64),parameter,public::pah_atomization_hf0_kjmol= &
       295d0-39.51733815680014d0+24*1.050d0+6*8.468d0
  ! Ground-state atomic/ionic references match snrt_chimes_atomization.h
  ! (ATcT v1.130), not graphite's approximate bulk latent heat.
  real(real64),parameter,public::pah_atomization_a12= &
       (24*711.396d0+12*216.034d0-pah_atomization_hf0_kjmol)*kjmol_to_erg
  real(real64),parameter,public::pah_atomization_ip=7.02d0*pah_atomization_ev
  real(real64),parameter,public::pah_atomization_carbon_ip=(1797.849d0-711.396d0)*kjmol_to_erg
  type,public::pah_atomization_receipt
     ! Number densities [cm^-3]. Carbon is neutral C, carbon_ion is C+.
     ! Electrons and H2 do not change in this specified endpoint.
     real(real64)::hydrogen=0,carbon=0,carbon_ion=0,photon_number=0
     ! Energy densities [erg cm^-3]. All positive except pah_level_removed,
     ! which includes the signed existing H-binding offset (H13 is negative).
     real(real64)::photon_energy=0,gas_heat=0,pah_level_removed=0
     real(real64)::binding_increase=0,gas_ionization_increase=0
     ! Identity: photon_energy = gas_heat - pah_level_removed
     !                         + binding_increase + gas_ionization_increase.
     ! Do NOT subtract binding again after using this material residual.
  end type
  public::pah_atomization_bond,pah_atomization_threshold,pah_atomization_batch
contains
  pure real(real64) function pah_atomization_bond(hydrogen) result(bond)
    ! Same H=0..13 offsets as pah_hydrogen_parameters, relative to neutral
    ! C24H12 + free H atoms; same explicit approximation for both charges.
    integer,intent(in)::hydrogen
    integer::h
    bond=ieee_value(0d0,ieee_quiet_nan)
    if(hydrogen<0.or.hydrogen>13)return
    bond=0
    if(hydrogen==13)then
       bond=-3.2d0*pah_atomization_ev
    else
       do h=12,hydrogen+1,-1
          bond=bond+merge(4.8d0,3.2d0,mod(h,2)==0)*pah_atomization_ev
       enddo
    endif
  end function

  pure real(real64) function pah_atomization_threshold(excitation,bond,charge) result(threshold)
    ! Minimum SINGLE photon energy [erg] for this grain's state. The input
    ! excitation excludes H bond/IP offsets. No tolerance promotes a photon
    ! below threshold. Quiet NaN signals invalid thermochemical inputs.
    real(real64),intent(in)::excitation,bond
    integer,intent(in)::charge
    real(real64)::ground
    threshold=ieee_value(0d0,ieee_quiet_nan)
    if(.not.all(ieee_is_finite([excitation,bond])))return
    if(excitation<0.or.charge<0.or.charge>1)return
    ground=pah_atomization_a12-bond+charge*(pah_atomization_carbon_ip-pah_atomization_ip)
    if(.not.ieee_is_finite(ground))return
    if(ground<=0)return
    threshold=max(0d0,ground-excitation)
  end function

  pure subroutine pah_atomization_batch(excitation,bond,hydrogen,charge,photon_energy,captured, &
       old,next,destroyed,receipt,ierr)
    ! Flattened states s, individual monochromatic event streams (group,s).
    ! Energies: erg/event; populations and integrated captures: cm^-3.
    ! Photon energies must be individual material-frame quanta, NOT a mixed
    ! hard/soft spectrum's energy/count ratio. Caller owns optical coverage,
    ! event-to-state assignment, momentum, transport retry and all survivors.
    ! Existing bond(s) may be supplied directly; its ground reference must
    ! match A12, h(s) and the caller's PAH energy ledger.
    !
    ! next/destroyed/receipt are overwritten only on success; on ANY failure
    ! all remain bitwise untouched. No gas state is mutated here. Successful
    ! destroyed(g,s) is removed from competing photo/heating event counts
    ! AND energies exactly once by the caller. Subthreshold events remain.
    ! Parent exhaustion is a distinct rejection, NEVER a yield cap, hidden
    ! energy budget, redistribution to cheaper grains, or silent photon loss.
    real(real64),intent(in)::excitation(:),bond(:),photon_energy(:,:),captured(:,:),old(:)
    integer,intent(in)::hydrogen(:),charge(:)
    real(real64),intent(inout)::next(:),destroyed(:,:)
    type(pah_atomization_receipt),intent(inout)::receipt
    integer,intent(out)::ierr
    real(real64)::trial(size(old)),events(size(captured,1),size(old)),ledger(9),terms(9),delta(9)
    real(real64)::threshold,need,removed,count,photon,excess,level,limit
    type(pah_atomization_receipt)::staged
    integer::ns,ng,s,g,k
    ierr=pah_atomization_invalid
    ns=size(old);ng=size(captured,1)
    if(ns<1.or.ng<1)return
    if(size(excitation)/=ns.or.size(bond)/=ns.or.size(hydrogen)/=ns.or.size(charge)/=ns)return
    if(size(next)/=ns.or.size(captured,2)/=ns)return
    if(any(shape(photon_energy)/=shape(captured)).or.any(shape(destroyed)/=shape(captured)))return
    if(.not.all(ieee_is_finite(old)).or..not.all(ieee_is_finite(captured)))return
    if(.not.all(ieee_is_finite(photon_energy)))return
    if(any(old<0).or.any(captured<0).or.any(photon_energy<0))return
    if(any(captured>0.and.photon_energy==0))return
    if(any(hydrogen<0).or.any(hydrogen>13))return
    trial=old;events=0;ledger=0;limit=huge(1d0)
    do s=1,ns
       threshold=pah_atomization_threshold(excitation(s),bond(s),charge(s))
       if(.not.ieee_is_finite(threshold))return
       if(old(s)==0.and.any(captured(:,s)>0))return
       need=pah_atomization_a12-bond(s)+charge(s)*(pah_atomization_carbon_ip-pah_atomization_ip)-excitation(s)
       level=excitation(s)+bond(s)+charge(s)*pah_atomization_ip
       removed=0
       do g=1,ng
          count=captured(g,s);photon=photon_energy(g,s)
          if(count==0.or.photon<threshold)cycle
          if(count>old(s)-removed)then
             ierr=pah_atomization_exhausted;return
          endif
          removed=removed+count
          if(removed>old(s))then
             ierr=pah_atomization_exhausted;return
          endif
          ! eps - (D-u) avoids summing two large positive energies. Guard
          ! genuinely overflowing excess/products/sums before arithmetic.
          if(need<0)then
             if(photon>limit+need)return
          endif
          excess=photon-need
          terms=[real(hydrogen(s),real64),real(24-charge(s),real64),real(charge(s),real64), &
               1d0,photon,excess,level,pah_atomization_a12,charge(s)*pah_atomization_carbon_ip]
          if(.not.all(ieee_is_finite(terms)))return
          if(count>limit/max(1d0,maxval(abs(terms))))return
          delta=count*terms
          do k=1,size(ledger)
             if(delta(k)>0)then
                if(ledger(k)>limit-delta(k))return
             else if(delta(k)<0)then
                if(ledger(k)<-limit-delta(k))return
             endif
          enddo
          ledger=ledger+delta;events(g,s)=count
       enddo
       trial(s)=old(s)-removed
    enddo
    if(.not.all(ieee_is_finite(ledger)).or..not.all(ieee_is_finite(trial)))return
    staged%hydrogen=ledger(1);staged%carbon=ledger(2);staged%carbon_ion=ledger(3)
    staged%photon_number=ledger(4);staged%photon_energy=ledger(5);staged%gas_heat=ledger(6)
    staged%pah_level_removed=ledger(7);staged%binding_increase=ledger(8);staged%gas_ionization_increase=ledger(9)
    next=trial;destroyed=events;receipt=staged;ierr=pah_atomization_ok
  end subroutine
end module

! Representation adapter, not a new chemistry or grain model. The existing
! H/He + dust partition owns the inventory cap/return approximation.
module snrt_moment_material
  use snrt_moment_transport, only: mn_basis,mn_dp,mn_ok,mn_bad_input,mn_reconstruct
  use snrt_dust_coupling, only: snrt_dust_partition_group,dust_coupling_ok
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: mn_absorb_groups
contains
  subroutine mn_absorb_groups(b,number,energy,tau_hhe,tau_dust,available, &
       absorbed_number,absorbed_energy,absorbed_first,returned,unassigned,ierr,validated_input)
    ! Local direction-independent opacity only. number/energy contain actual
    ! moment densities in consistent units, NOT a signed energy correction.
    ! Keep both fields: energy need not equal E_reference*number. A spectral
    ! band whose effective opacity depends on direction needs its existing
    ! angular material adapter, not this scalar-opacity routine.
    ! Outputs are (H I,He I,He II,dust; group). absorbed_first is the energy
    ! first moment WITHOUT c_hat: the live receiver converts energy units and
    ! divides by physical c exactly once to obtain impulse density.
    ! Returned/unassigned photons stay in radiation, not physical sinks.
    ! No caller state, inventory, or receipt is published on any failure.
    type(mn_basis),intent(in) :: b
    real(mn_dp),intent(inout) :: number(:,:),energy(:,:),available(3)
    real(mn_dp),intent(in) :: tau_hhe(:,:),tau_dust(:)
    real(mn_dp),intent(inout) :: absorbed_number(:,:),absorbed_energy(:,:),absorbed_first(:,:,:)
    real(mn_dp),intent(inout) :: returned(:),unassigned(:)
    integer,intent(out) :: ierr
    logical,optional,intent(in) :: validated_input
    real(mn_dp) :: next_n(size(number,1),size(number,2)),next_e(size(energy,1),size(energy,2))
    real(mn_dp) :: take_n(4,size(number,2)),take_e(4,size(number,2)),take_first(3,4,size(number,2))
    real(mn_dp) :: give_back(size(number,2)),unplaced(size(number,2)),inventory(3),angular(b%nq)
    real(mn_dp) :: tau,fraction,raw,assigned,channel_fraction,total_fraction
    integer :: ng,g,s,status
    logical :: already_valid
    ierr=mn_bad_input;ng=size(number,2)
    already_valid=.false.
    if(present(validated_input))already_valid=validated_input
    if(b%nm<4.or.size(number,1)/=b%nm.or.any(shape(energy)/=shape(number)))return
    if(any(shape(tau_hhe)/=[3,ng]).or.size(tau_dust)/=ng)return
    if(any(.not.ieee_is_finite(number)).or.any(.not.ieee_is_finite(energy)))return
    if(any(shape(absorbed_number)/=[4,ng]).or.any(shape(absorbed_energy)/=[4,ng]))return
    if(any(shape(absorbed_first)/=[3,4,ng]).or.size(returned)/=ng.or.size(unassigned)/=ng)return
    if(any(.not.ieee_is_finite(tau_hhe)).or.any(.not.ieee_is_finite(tau_dust)))return
    if(any(.not.ieee_is_finite(available)).or.any(available<0))return
    if(any(tau_hhe<0).or.any(tau_dust<0))return
    next_n=number;next_e=energy;inventory=available
    take_n=0;take_e=0;take_first=0;give_back=0;unplaced=0
    do g=1,ng
       ! Scalar nonnegative attenuation preserves the realizability cone.
       ! The live driver has just validated the transported input; do not
       ! solve the identical entropy problem again at each chemistry iterate.
       if(.not.already_valid)then
       call mn_reconstruct(b,number(:,g),angular,ierr)
       if(ierr/=mn_ok)return
       call mn_reconstruct(b,energy(:,g),angular,ierr)
       if(ierr/=mn_ok)return
       endif
       ierr=mn_bad_input
       if(number(1,g)<0.or.energy(1,g)<0)return
       if(number(1,g)==0.and.any(energy(:,g)/=0))return
       if(number(1,g)>0.and.energy(1,g)<=0)return
       tau=sum(tau_hhe(:,g))+tau_dust(g)
       if(.not.ieee_is_finite(tau))return
       if(tau==0.or.number(1,g)==0)cycle
       if(tau<1d-6)then
          fraction=tau*(1-tau/2+tau*tau/6-tau*tau*tau/24)
       else
          fraction=1-exp(-tau)
       endif
       raw=number(1,g)*fraction
       call snrt_dust_partition_group(raw,tau_hhe(:,g),tau_dust(g),inventory, &
            take_n(1:3,g),take_n(4,g),give_back(g),unplaced(g),status)
       if(status/=dust_coupling_ok)return
       assigned=sum(take_n(:,g));total_fraction=assigned/number(1,g)
       if(.not.ieee_is_finite(total_fraction).or.total_fraction<0.or.total_fraction>1)return
       next_n(:,g)=number(:,g)*(1-total_fraction)
       next_e(:,g)=energy(:,g)*(1-total_fraction)
       do s=1,4
          channel_fraction=take_n(s,g)/number(1,g)
          take_e(s,g)=energy(1,g)*channel_fraction
          take_first(:,s,g)=[-energy(3,g),-energy(4,g),energy(2,g)]*(channel_fraction/sqrt(3d0))
       enddo
       ! next_n/e are the validated input times a factor in [0,1].
    enddo
    number=next_n;energy=next_e;available=inventory
    absorbed_number=take_n;absorbed_energy=take_e;absorbed_first=take_first
    returned=give_back;unassigned=unplaced;ierr=mn_ok
  end subroutine
end module

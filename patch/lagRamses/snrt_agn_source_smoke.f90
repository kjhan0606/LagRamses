program snrt_agn_source_smoke
  use, intrinsic :: iso_c_binding, only: c_float
  use amr_parameters, only: dp
  use snrt_agn_source, only: snrt_agn_photon_budget, snrt_agn_isotropic_packet, &
       snrt_agn_photons_to_density_code, snrt_agn_deposit_isotropic, &
       snrt_agn_deposit_transaction, snrt_c_cgs, snrt_ev_to_erg, &
       snrt_agn_photon_budget_energy, snrt_agn_source_commit
  use snrt_rt_transaction, only: snrt_rt_transaction_snapshot, &
       snrt_transaction_begin, snrt_transaction_restore, &
       snrt_transaction_commit_level, snrt_transaction_ok
  implicit none

  real(dp) :: luminosity_erg_s, emitted_photons
  real(dp) :: expected_energy_erg, expected_photons
  real(dp) :: angular_weights(4), directional_photons(4)
  real(dp) :: photon_density_code, transaction_deposited
  real(dp) :: emitted_groups(2)
  real(dp) :: pending,reference_photons
  real(c_float) :: state(4,2,1), state_before(4,2,1)
  real(c_float) :: coupled_state(4,2,1), coupled_state_before(4,2,1), &
       coupled_trial(4,2,1), coupled_coarse(4,2,1)
  integer :: coupled_leaf(1)
  real(dp) :: coupled_h(1), coupled_heii(1), coupled_heiii(1), &
       coupled_neutral(1), coupled_thermal(1)
  real(dp) :: coupled_trial_h(1), coupled_trial_heii(1), &
       coupled_trial_heiii(1), coupled_trial_neutral(1), coupled_trial_thermal(1)
  type(snrt_rt_transaction_snapshot) :: coupled_transaction
  integer :: ierr

  ! The first positional argument is supplied inflow, not retained BH mass.
  call snrt_agn_photon_budget(2.0d0, 5.0d33, 4.0d0, 0.1d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)

  expected_energy_erg = 0.1d0 * 2.0d0 * 5.0d33 * snrt_c_cgs**2
  expected_photons = 0.25d0 * expected_energy_erg / (20.0d0 * snrt_ev_to_erg)
  if (abs(luminosity_erg_s - expected_energy_erg / 4.0d0) / &
       (expected_energy_erg / 4.0d0) > 1.0d-13) error stop 1
  if (abs(emitted_photons - expected_photons) / expected_photons > 1.0d-13) error stop 2
  call snrt_agn_photon_budget_energy(expected_energy_erg,4d0,0.25d0,20d0, &
       luminosity_erg_s,emitted_photons,ierr)
  if(ierr/=0.or.abs(emitted_photons/expected_photons-1d0)>1d-13)error stop 15
  call snrt_agn_photon_budget_energy(0d0,4d0,0.25d0,20d0,luminosity_erg_s,emitted_photons,ierr)
  if(ierr/=0.or.luminosity_erg_s/=0d0.or.emitted_photons/=0d0)error stop 16
  call snrt_agn_photon_budget_energy(-1d0,4d0,0.25d0,20d0,luminosity_erg_s,emitted_photons,ierr)
  if(ierr==0)error stop 17

  call snrt_agn_photon_budget(0.0d0, 5.0d33, 4.0d0, 0.1d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)
  if (luminosity_erg_s /= 0.0d0 .or. emitted_photons /= 0.0d0) error stop 3

  call snrt_agn_photon_budget(2.0d0, 5.0d33, 4.0d0, 1.0d0, 0.25d0, 20.0d0, &
       luminosity_erg_s, emitted_photons)
  if (luminosity_erg_s /= 0.0d0 .or. emitted_photons /= 0.0d0) error stop 14

  angular_weights = (/1.0d0, 2.0d0, -1.0d0, 1.0d0/)
  call snrt_agn_isotropic_packet(12.0d0, angular_weights, directional_photons)
  if (abs(sum(directional_photons) - 12.0d0) > 1.0d-13) error stop 4
  if (maxval(abs(directional_photons - (/3.0d0, 6.0d0, 0.0d0, 3.0d0/))) > &
       1.0d-13) error stop 5

  ! Production deposition requires a non-negative quadrature.  Keep the
  ! negative-weight compatibility check above local to the pure splitter.
  angular_weights = (/1.0d0, 2.0d0, 1.0d0, 1.0d0/)

  call snrt_agn_photons_to_density_code(1.0d63, 125.0d0, 1.0d21, 1.0d-3, &
       photon_density_code)
  if (abs(photon_density_code - 8.0d0) > 1.0d-13) error stop 6

  state = 0.0_c_float
  call snrt_agn_deposit_isotropic(state, 1, 1, 1.0d63, 125.0d0, 1.0d21, &
       1.0d-3, angular_weights, photon_density_code, ierr)
  if (ierr /= 0) error stop 7
  if (abs(sum(real(state,dp)) - 8.0d0) > 1.0d-6) error stop 8

  emitted_groups = (/1.0d63, 2.0d63/)
  call snrt_agn_deposit_transaction(state, 1, emitted_groups, 125.0d0, &
       1.0d21, 1.0d-3, angular_weights, transaction_deposited, ierr)
  if (ierr /= 0) error stop 9
  if (abs(transaction_deposited - 24.0d0) > 1.0d-12) error stop 10
  if (abs(sum(real(state,dp)) - 32.0d0) > 1.0d-6) error stop 11

  state_before = state
  emitted_groups = (/1.0d63, -1.0d0/)
  call snrt_agn_deposit_transaction(state, 1, emitted_groups, 125.0d0, &
       1.0d21, 1.0d-3, angular_weights, transaction_deposited, ierr)
  if (ierr == 0) error stop 12
  if (maxval(abs(real(state,dp) - real(state_before,dp))) > 0.0d0) error stop 13

  ! Accepted-energy receipt survives a failed all-group source transaction.
  pending=1d53; state=0
  call snrt_agn_photon_budget_energy(pending,4d0,0.25d0,20d0,luminosity_erg_s,emitted_groups(1),ierr)
  reference_photons=emitted_groups(1)
  emitted_groups(2)=-1d0
  call snrt_agn_deposit_transaction(state,1,emitted_groups,125d0,1d21,1d-3,angular_weights, &
       transaction_deposited,ierr)
  call snrt_agn_source_commit(pending,ierr==0)
  if(ierr==0.or.pending/=1d53.or.any(state/=0))error stop 18
  ! Retry succeeds. Unrepresented bolometric energy is not renormalized.
  emitted_groups(2)=reference_photons
  call snrt_agn_deposit_transaction(state,1,emitted_groups,125d0,1d21,1d-3,angular_weights, &
       transaction_deposited,ierr)
  call snrt_agn_source_commit(pending,ierr==0)
  if(ierr/=0.or.pending/=0d0)error stop 19

  ! The production driver snapshots the persistent RT state before source
  ! injection and consumes pending fuel only after the coupled transaction.
  ! Exercise that ordering directly: a failed coupled step must remove both
  ! the staged photons and no accepted-event fuel; the retry then commits once.
  coupled_state = 0.0_c_float
  coupled_state_before = coupled_state
  coupled_coarse = 0.0_c_float
  coupled_leaf = (/1/)
  coupled_h = 0.10d0
  coupled_heii = 0.05d0
  coupled_heiii = 0.01d0
  coupled_neutral = 0.90d0
  coupled_thermal = 10.0d0
  pending = 1d53
  call snrt_transaction_begin(coupled_transaction, coupled_state, coupled_leaf, &
       coupled_h, coupled_heii, coupled_heiii, coupled_neutral, coupled_thermal, ierr)
  if(ierr/=snrt_transaction_ok.or..not.coupled_transaction%active)error stop 21
  call snrt_agn_photon_budget_energy(pending,4d0,0.25d0,20d0,luminosity_erg_s,emitted_groups(1),ierr)
  emitted_groups(2)=emitted_groups(1)
  call snrt_agn_deposit_transaction(coupled_state,1,emitted_groups,125d0,1d21,1d-3, &
       angular_weights,transaction_deposited,ierr)
  if(ierr/=0.or.all(coupled_state==coupled_state_before))error stop 22
  call snrt_transaction_restore(coupled_transaction,coupled_state,coupled_leaf,coupled_h, &
       coupled_heii,coupled_heiii,coupled_neutral,coupled_thermal,ierr)
  if(ierr/=snrt_transaction_ok.or..not.all(coupled_state==coupled_state_before).or. &
       pending/=1d53)error stop 23
  write(*,'(A)') 'SNRT_AGN_SOURCE_COUPLED_ROLLBACK_PENDING_PASS'

  call snrt_transaction_begin(coupled_transaction, coupled_state, coupled_leaf, coupled_h, &
       coupled_heii, coupled_heiii, coupled_neutral, coupled_thermal, ierr)
  if(ierr/=snrt_transaction_ok.or..not.coupled_transaction%active)error stop 24
  call snrt_agn_photon_budget_energy(pending,4d0,0.25d0,20d0,luminosity_erg_s,emitted_groups(1),ierr)
  emitted_groups(2)=emitted_groups(1)
  call snrt_agn_deposit_transaction(coupled_state,1,emitted_groups,125d0,1d21,1d-3, &
       angular_weights,transaction_deposited,ierr)
  coupled_trial = coupled_state
  coupled_trial_h = coupled_h
  coupled_trial_heii = coupled_heii
  coupled_trial_heiii = coupled_heiii
  coupled_trial_neutral = coupled_neutral
  coupled_trial_thermal = coupled_thermal
  call snrt_transaction_commit_level(coupled_transaction,coupled_state,coupled_leaf, &
       coupled_h,coupled_heii,coupled_heiii,coupled_neutral,coupled_trial,coupled_coarse, &
       coupled_trial_h,coupled_trial_heii,coupled_trial_heiii,coupled_trial_neutral, &
       coupled_thermal,coupled_trial_thermal,ierr)
  if(ierr/=snrt_transaction_ok.or.coupled_transaction%active)error stop 25
  call snrt_agn_source_commit(pending,.true.)
  if(pending/=0d0.or.all(coupled_state==coupled_state_before))error stop 26
  write(*,'(A)') 'SNRT_AGN_SOURCE_COUPLED_COMMIT_PENDING_CLEAR_PASS'

  call paired_source_checks()

  write(*,'(a,es14.6,a,es14.6)') 'SNRT_AGN_SOURCE_OK luminosity=', &
       expected_energy_erg / 4.0d0, ' photons=', expected_photons
contains
  subroutine paired_source_checks()
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    real(c_float)::n(4,2,2),saved_n(4,2,2),legacy(4,2,2)
    real(dp)::shift(4,2,2),saved_shift(4,2,2),reference(2),sources(2),weights(4)
    real(dp)::increment(4,2),expected_shift(4,2),energy_before(2),energy_after(2),added,legacy_added
    real(dp)::h(2),heii(2),heiii(2),neutral(2),heat(1)
    type(snrt_rt_transaction_snapshot)::transaction
    integer::g,k,status,legacy_status
    reference=[20._dp,80._dp];weights=[1._dp,2._dp,0._dp,3._dp]
    n(:,1,1)=[16777216._c_float,.25_c_float,4._c_float,1e-20_c_float]
    n(:,2,1)=2*n(:,1,1);n(:,:,2)=7
    shift(:,:,1)=.125_dp;shift(1,:,1)=-2._dp;shift(:,:,2)=-1
    saved_n=n;saved_shift=shift;legacy=n;sources=[.7_dp,.13_dp]
    do g=1,2
       increment(:,g)=sources(g)*weights/sum(weights)
       energy_before(g)=sum(reference(g)*real(n(:,g,1),dp)+shift(:,g,1))
    enddo
    call snrt_agn_deposit_transaction(legacy,1,sources,1._dp,1._dp,1._dp,weights,legacy_added,legacy_status)
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status, &
         persistent_energy_shift=shift,group_mean_energy_ev=reference)
    if(status/=0.or.legacy_status/=0.or.any(n/=legacy).or.added/=legacy_added)error stop 30
    if(any(n(:,:,2)/=saved_n(:,:,2)).or.any(shift(:,:,2)/=saved_shift(:,:,2)))error stop 31
    do g=1,2
       expected_shift(:,g)=saved_shift(:,g,1)+reference(g)* &
            ((real(saved_n(:,g,1),dp)-real(n(:,g,1),dp))+increment(:,g))
       energy_after(g)=sum(reference(g)*real(n(:,g,1),dp)+shift(:,g,1))
       if(abs(energy_after(g)-energy_before(g)-reference(g)*sources(g))> &
            8*epsilon(1._dp)*energy_before(g))error stop 32
    enddo
    if(maxval(abs(shift(:,:,1)-expected_shift))>1d-13)error stop 33
    ! A zero angular weight leaves both members of that bin unchanged.
    if(any(n(3,:,1)/=saved_n(3,:,1)).or.any(shift(3,:,1)/=saved_shift(3,:,1)))error stop 34
    saved_n=n;saved_shift=shift
    call snrt_agn_deposit_transaction(n,1,[0._dp,0._dp],1._dp,1._dp,1._dp,weights,added,status,shift,reference)
    if(status/=0.or.added/=0.or.any(n/=saved_n).or.any(shift/=saved_shift))error stop 35

    ! Source increments below even the FP64 ulp of oldN must survive in shift.
    n=16777216;shift=0;sources=1d-10;weights=1
    do k=1,100
       call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference)
       if(status/=0.or.any(n/=16777216))error stop 36
    enddo
    do g=1,2
       if(maxval(abs(shift(:,g,1)-100*reference(g)*sources(g)/4))>1d-20)error stop 37
    enddo
    if(any(shift(:,:,2)/=0))error stop 38

    ! Reject missing pairs, dimensions, invalid energy, and late-group failures.
    do k=1,11
       n=1;shift=0;reference=[20._dp,80._dp];sources=1
       select case(k)
       case(4);reference(2)=0
       case(5);reference(2)=ieee_value(1._dp,ieee_quiet_nan)
       case(6);shift(4,2,1)=-81
       case(7);n(4,2,1)=0;shift(4,2,1)=1
       case(8);reference(2)=huge(1._dp)*.1_dp;sources(2)=100
       case(9);n(:,2,1)=0;sources(2)=1d-300
       case(10);sources(2)=-1
       case(11);sources(2)=real(huge(0._c_float),dp)*4
       end select
       saved_n=n;saved_shift=shift
       select case(k)
       case(1)
          call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status, &
               persistent_energy_shift=shift)
       case(2)
          call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status, &
               group_mean_energy_ev=reference)
       case(3)
          call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift(:,:,:1),reference)
       case default
          call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference)
       end select
       if(status==0.or.added/=0.or.any(n/=saved_n).or.any(shift/=saved_shift))error stop 39
    enddo

    ! The enclosing source/transport transaction restores the correction as
    ! well as N, including a source too small to change FP32 photon number.
    n=16777216;shift=0;sources=1d-8;reference=[20._dp,80._dp]
    h=.1_dp;heii=0;heiii=0;neutral=.9_dp;heat=1
    saved_n=n;saved_shift=shift
    do k=1,2
       call snrt_transaction_begin(transaction,n,[1],h,heii,heiii,neutral,heat,status,shift,reference)
       if(status/=snrt_transaction_ok)error stop 40
       call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference)
       if(status/=0.or.any(shift(:,1,1)<=0))error stop 41
       call snrt_transaction_restore(transaction,n,[1],h,heii,heiii,neutral,heat,status,shift)
       if(status/=snrt_transaction_ok.or.any(n/=saved_n).or.any(shift/=saved_shift))error stop 42
    enddo
    n=0;shift=0;sources=[1._dp,1d-300];reference=[20._dp,80._dp];weights=1
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference, &
         quantize_source=.true.)
    if(status/=0.or.any(n(:,2,1)/=0).or.any(shift/=0).or.added/=1)error stop 43
    if(abs(sum(real(n(:,1,1),dp))*20-20)>1d-12)error stop 44
    n=0;shift=0;sources=1d-300
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference, &
         quantize_source=.true.)
    if(status==0.or.any(n/=0).or.any(shift/=0).or.added/=0)error stop 45
    ! Actual stellar mean may lie below/above Eref: encoding stays unchanged.
    n=0;shift=0;sources=[.7_dp,.13_dp];weights=[1._dp,2._dp,0._dp,3._dp]
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference, &
         quantize_source=.true.,source_mean_energy_ev=[15._dp,100._dp])
    if(status/=0.or.any(n(3,:,1)/=0).or.any(shift(3,:,1)/=0))error stop 46
    do g=1,2
       energy_after(g)=sum(reference(g)*real(n(:,g,1),dp)+shift(:,g,1))
    enddo
    if(any(abs(energy_after-sources*[15._dp,100._dp])>8*epsilon(0._c_float)*sum(sources*[15._dp,100._dp]))) &
         error stop 47
    if(any(shift([1,2,4],1,1)>=0).or.any(shift([1,2,4],2,1)<=0))error stop 48
    saved_n=n;saved_shift=shift
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference, &
         source_mean_energy_ev=[15._dp,-1._dp])
    if(status/=6.or.added/=0.or.any(n/=saved_n).or.any(shift/=saved_shift))error stop 49
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status, &
         source_mean_energy_ev=[15._dp,100._dp])
    if(status/=6.or.added/=0.or.any(n/=saved_n))error stop 50
    ! Using Eref for the rounding bound would wrongly accept this lost source.
    n=0;shift=0;weights=1;sources=[1._dp,1d-46];reference=1
    call snrt_agn_deposit_transaction(n,1,sources,1._dp,1._dp,1._dp,weights,added,status,shift,reference, &
         quantize_source=.true.,source_mean_energy_ev=[1._dp,1d46])
    if(status/=6.or.added/=0.or.any(n/=0).or.any(shift/=0))error stop 51
    write(*,'(A)')'SNRT_ACTUAL_STELLAR_ENERGY_INJECTION_AND_ROUNDING_BOUND_PASS'
    write(*,'(A)')'SNRT_BAND_REPRESENTABLE_SOURCE_TAIL_ENERGY_BOUND_PASS'
    write(*,'(A)')'SNRT_PAIRED_SOURCE_FP32_REBASE_EXACT_INCREMENT_LEGACY_PARITY_ROLLBACK_PASS'
  end subroutine paired_source_checks
end program snrt_agn_source_smoke

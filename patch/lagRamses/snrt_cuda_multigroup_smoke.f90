program snrt_cuda_multigroup_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use snrt_cuda_multigroup_interface, only: snrt_cuda_available, &
       snrt_cuda_multigroup_rt_step_species, snrt_cuda_multigroup_rt_step_species_dust
  implicit none

  integer, parameter :: nx = 7, ny = 5, nz = 3, ncell = nx * ny * nz
  ! Exercise the production nine-group CUDA ABI, including the hard-X-ray
  ! group that was absent from the old four-group smoke.
  integer, parameter :: ndirection = 80, ngroup = 9
  integer :: ix, iy, iz, icell, idir, igroup, linear
  integer(c_int) :: ierr, neighbor(6,ncell)
  real(c_float) :: state(ncell,ndirection,ngroup)
  real(c_float) :: direction(3,ndirection), tau(ncell,ngroup)
  real(c_float) :: tau_species(ncell,ngroup,3)
  real(c_float) :: available_species(ncell,3), initial_species(ncell,3)
  real(c_float) :: absorbed(ncell)
  real(c_float) :: absorbed_group(ncell,ngroup)
  real(c_float) :: tau_dust(ncell,ngroup), tau_species_dust(ncell,ngroup,3)
  real(c_float) :: available_species_dust(ncell,3), initial_species_dust(ncell,3)
  real(c_float) :: absorbed_hhe_species(ncell,ngroup,3)
  real(c_float) :: absorbed_dust_group(ncell,ngroup), returned_group(ncell,ngroup)
  real(c_float) :: raw_group(ncell,ngroup), assigned_group(ncell,ngroup)
  real(c_float) :: state_dust(ncell,ndirection,ngroup)
  real(c_float) :: state_zero_reference(ncell,ndirection,ngroup)
  real(c_float) :: available_zero_reference(ncell,3)
  real(c_float) :: absorbed_zero_reference(ncell)
  real(c_float) :: absorbed_group_zero_reference(ncell,ngroup)
  real(c_float) :: tau_zero_reference(ncell,ngroup)
  real(c_float) :: tau_species_zero_reference(ncell,ngroup,3)
  real(c_float) :: tau_dust_zero_reference(ncell,ngroup)
  integer(c_int) :: state_bits(size(state_dust)), state_zero_bits(size(state_zero_reference))
  integer(c_int) :: available_bits(size(available_species_dust))
  integer(c_int) :: available_zero_bits(size(available_zero_reference))
  integer(c_int) :: group_bits(size(assigned_group))
  integer(c_int) :: group_zero_bits(size(absorbed_group_zero_reference))
  integer(c_int) :: scalar_bits(size(absorbed)), scalar_zero_bits(size(absorbed_zero_reference))
  real(c_double) :: before, after, removed, budget_error
  real(c_double) :: ledger_error, inventory_error, zero_error, hhe_zero_error
  real(c_double) :: expected_available(3), expected_hhe_tau, expected_total_tau
  real(c_double) :: expected_eligible, expected_hhe_target, expected_excess
  real(c_double) :: expected_dust, expected_fraction, expected_scale

  if (snrt_cuda_available() <= 0_c_int) then
     write(*,'(a)') 'SNRT_NATIVE_CUDA_UNAVAILABLE'
     error stop 90
  end if

  do iz = 0, nz-1
     do iy = 0, ny-1
        do ix = 0, nx-1
           icell = 1 + ix + nx * (iy + ny * iz)
           neighbor(1,icell) = 1 + modulo(ix-1,nx) + nx * (iy + ny * iz)
           neighbor(2,icell) = 1 + modulo(ix+1,nx) + nx * (iy + ny * iz)
           neighbor(3,icell) = 1 + ix + nx * (modulo(iy-1,ny) + ny * iz)
           neighbor(4,icell) = 1 + ix + nx * (modulo(iy+1,ny) + ny * iz)
           neighbor(5,icell) = 1 + ix + nx * (iy + ny * modulo(iz-1,nz))
           neighbor(6,icell) = 1 + ix + nx * (iy + ny * modulo(iz+1,nz))
        end do
     end do
  end do

  do idir = 1, ndirection
     direction(:,idir) = (/0.30_c_float, -0.20_c_float, 0.40_c_float/)
  end do
  do igroup = 1, ngroup
     do icell = 1, ncell
        tau(icell,igroup) = real(0.1d0 * dble(igroup),c_float)
        state(icell,:,igroup) = real(1 + modulo(11*icell + 17*igroup,23),c_float)
        tau_species(icell,igroup,:) = 0.0_c_float
        if (igroup == 1) then
           tau_species(icell,igroup,1) = tau(icell,igroup)
        else if (igroup == 2) then
           tau_species(icell,igroup,1) = 0.5_c_float*tau(icell,igroup)
           tau_species(icell,igroup,2) = 0.5_c_float*tau(icell,igroup)
        else
           tau_species(icell,igroup,1) = 0.4_c_float*tau(icell,igroup)
           tau_species(icell,igroup,2) = 0.3_c_float*tau(icell,igroup)
           tau_species(icell,igroup,3) = 0.3_c_float*tau(icell,igroup)
        end if
     end do
  end do
  do icell = 1, ncell
     available_species(icell,1) = 0.0_c_float
     available_species(icell,2) = 0.0_c_float
     available_species(icell,3) = real(1 + modulo(3*icell,5),c_float)
  end do
  initial_species = available_species

  before = sum(real(state,c_double))
  ierr = snrt_cuda_multigroup_rt_step_species(state, direction, neighbor, tau, &
       tau_species, available_species, absorbed, absorbed_group, int(ncell,c_int), &
       int(ncell,c_int), int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 1
  after = sum(real(state,c_double))
  removed = sum(real(absorbed,c_double))
  budget_error = abs(before-after-removed) / before
  if (budget_error > 2.0d-6) error stop 2
  if (maxval(real(absorbed,c_double) - &
       sum(real(initial_species,c_double),dim=2)) > 1.0d-5) error stop 3
  if (maxval(abs(real(sum(absorbed_group,dim=2),c_double) - &
       real(absorbed,c_double))) > 2.0d-5) error stop 4
  if (maxval(real(sum(absorbed_group(:,1:2),dim=2),c_double)) > 1.0d-6) error stop 5
  if (minval(available_species) < -1.0e-5_c_float) error stop 6
  if (maxval(abs(real(sum(initial_species,dim=2)-sum(available_species,dim=2),c_double) - &
       real(absorbed,c_double))) > 2.0d-5) error stop 7
  if (minval(state) < 0.0_c_float .or. minval(absorbed) < 0.0_c_float .or. &
       minval(absorbed_group) < 0.0_c_float) error stop 8

  write(*,'(a,es14.6)') 'SNRT_CUDA_MULTIGROUP_OK relative_budget_error=', budget_error

  ! Repeat with all three absorbing species present and deliberately
  ! different inventories.  This exercises opacity-weighted saturation and
  ! redistribution, rather than only the zero-H/He inventory edge case above.
  state = 1.0_c_float
  available_species(:,1) = 0.050_c_float
  available_species(:,2) = 0.200_c_float
  available_species(:,3) = 0.350_c_float
  initial_species = available_species
  before = sum(real(state,c_double))
  ierr = snrt_cuda_multigroup_rt_step_species(state, direction, neighbor, tau, &
       tau_species, available_species, absorbed, absorbed_group, int(ncell,c_int), &
       int(ncell,c_int), int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 9
  after = sum(real(state,c_double))
  removed = sum(real(absorbed,c_double))
  budget_error = abs(before-after-removed) / before
  if (budget_error > 2.0d-6) error stop 10
  if (maxval(real(absorbed,c_double) - &
       sum(real(initial_species,c_double),dim=2)) > 1.0d-5) error stop 11
  if (maxval(abs(real(sum(absorbed_group,dim=2),c_double) - &
       real(absorbed,c_double))) > 2.0d-5) error stop 12
  if (minval(available_species) < -1.0e-5_c_float .or. minval(state) < 0.0_c_float .or. &
       minval(absorbed) < 0.0_c_float .or. minval(absorbed_group) < 0.0_c_float) error stop 13
  if (maxval(abs(real(sum(initial_species,dim=2)-sum(available_species,dim=2),c_double) - &
       real(absorbed,c_double))) > 2.0d-5) error stop 14
  write(*,'(a,es14.6)') 'SNRT_CUDA_MULTIGROUP_SPECIES_MIX_OK relative_budget_error=', &
       budget_error

  ! DUST-7: exercise the separate H/He+dust ABI.  Groups 1--4 are
  ! dust-only production groups below 13.6 eV; groups 5--9 contain both
  ! H/He and dust.  The H/He inventory must not be consumed by the former.
  state_dust = 1.0_c_float
  do igroup = 1, ngroup
     do icell = 1, ncell
        tau_species_dust(icell,igroup,:) = 0.0_c_float
        if (igroup <= 4) then
           tau_dust(icell,igroup) = real(0.05d0*dble(igroup),c_float)
        else
           tau_dust(icell,igroup) = real(0.01d0*dble(igroup),c_float)
           tau_species_dust(icell,igroup,1) = 0.30_c_float
           tau_species_dust(icell,igroup,2) = 0.20_c_float
           tau_species_dust(icell,igroup,3) = 0.50_c_float
        end if
        tau(icell,igroup) = tau_dust(icell,igroup) + &
             sum(tau_species_dust(icell,igroup,:))
     end do
  end do
  do icell = 1, ncell
     available_species_dust(icell,1) = 0.030_c_float
     available_species_dust(icell,2) = 0.100_c_float
     available_species_dust(icell,3) = 0.350_c_float
  end do
  initial_species_dust = available_species_dust
  before = sum(real(state_dust,c_double))
  ierr = snrt_cuda_multigroup_rt_step_species_dust(state_dust, direction, neighbor, tau, &
       tau_species_dust, tau_dust, available_species_dust, absorbed_hhe_species, &
       absorbed_dust_group, returned_group, raw_group, assigned_group, absorbed, &
       int(ncell,c_int), int(ncell,c_int), int(ndirection,c_int), int(ngroup,c_int), &
       0.10_c_float)
  if (ierr /= 0_c_int) error stop 15
  after = sum(real(state_dust,c_double))
  removed = sum(real(absorbed,c_double))
  budget_error = abs(before-after-removed) / before
  if (budget_error > 2.0d-6) error stop 16
  ledger_error = maxval(abs(real(raw_group - sum(absorbed_hhe_species,dim=3) - &
       absorbed_dust_group - returned_group,c_double)))
  if (ledger_error > 3.0d-5) error stop 17
  if (maxval(abs(real(raw_group - returned_group - assigned_group,c_double))) > 3.0d-5) &
       error stop 18
  inventory_error = maxval(abs(real(initial_species_dust - available_species_dust - &
       sum(absorbed_hhe_species,dim=2),c_double)))
  if (inventory_error > 3.0d-5) error stop 19
  if (maxval(absorbed_hhe_species(:,1:4,:)) > 0.0_c_float) error stop 20
  if (maxval(abs(absorbed_dust_group(:,1:4) - raw_group(:,1:4))) /= 0.0_c_float) error stop 21
  if (maxval(abs(returned_group(:,1:4))) /= 0.0_c_float) error stop 22
  if (maxval(absorbed_dust_group(:,1:4)) <= 0.0_c_float) error stop 23
  if (maxval(returned_group) <= 0.0_c_float) error stop 24
  if (maxval(available_species_dust) > 3.0e-5_c_float) error stop 25
  if (maxval(absorbed_hhe_species(:,6:9,:)) > 3.0e-5_c_float) error stop 26
  ! Reconstruct the DUST-6 FP64 split in host precision.  This checks the
  ! proportional, direct-dust, guard-band, and finite-excess paths rather
  ! than only checking that the returned ledgers close.
  do icell = 1, ncell
     expected_available = real(initial_species_dust(icell,:),c_double)
     do igroup = 1, ngroup
        expected_hhe_tau = sum(real(tau_species_dust(icell,igroup,:),c_double))
        expected_total_tau = expected_hhe_tau + real(tau_dust(icell,igroup),c_double)
        if (expected_total_tau > 0.0d0 .and. expected_hhe_tau > 0.0d0) then
           expected_hhe_target = min(real(raw_group(icell,igroup),c_double), &
                real(raw_group(icell,igroup),c_double)*expected_hhe_tau/expected_total_tau)
        else
           expected_hhe_target = 0.0d0
        end if
        expected_eligible = 0.0d0
        do idir = 1, 3
           if (tau_species_dust(icell,igroup,idir) > 0.0_c_float) &
                expected_eligible = expected_eligible + expected_available(idir)
        end do
        expected_excess = max(0.0d0, expected_hhe_target - expected_eligible)
        expected_fraction = 1.0d0 - exp(-real(tau_dust(icell,igroup),c_double))
        expected_dust = max(0.0d0, real(raw_group(icell,igroup),c_double) - &
             expected_hhe_target) + expected_excess*expected_fraction
        expected_scale = max(1.0d0, abs(real(raw_group(icell,igroup),c_double)))
        if (abs(real(absorbed_dust_group(icell,igroup),c_double)-expected_dust) > &
             5.0d-5*expected_scale) error stop 27
        expected_available = expected_available - &
             real(absorbed_hhe_species(icell,igroup,:),c_double)
     end do
  end do
  if (minval(state_dust) < 0.0_c_float .or. minval(available_species_dust) < 0.0_c_float .or. &
       minval(absorbed_hhe_species) < 0.0_c_float .or. minval(absorbed_dust_group) < 0.0_c_float .or. &
       minval(returned_group) < 0.0_c_float .or. minval(raw_group) < 0.0_c_float .or. &
       minval(assigned_group) < 0.0_c_float) error stop 28
  write(*,'(a,es14.6)') 'SNRT_CUDA_MULTIGROUP_SPECIES_DUST_OK relative_budget_error=', &
       budget_error

  ! Non-saturating mixed group: group 5 has ample H/He inventory relative to
  ! its deliberately small photon packet, so the opacity-proportional split
  ! is exercised independently of the finite-reservoir branch.
  state_dust = 0.0_c_float
  state_dust(:,:,5) = 1.0e-3_c_float
  available_species_dust = initial_species_dust
  ierr = snrt_cuda_multigroup_rt_step_species_dust(state_dust, direction, neighbor, tau, &
       tau_species_dust, tau_dust, available_species_dust, absorbed_hhe_species, &
       absorbed_dust_group, returned_group, raw_group, assigned_group, absorbed, &
       int(ncell,c_int), int(ncell,c_int), int(ndirection,c_int), int(ngroup,c_int), &
       0.10_c_float)
  if (ierr /= 0_c_int) error stop 29
  do icell = 1, ncell
     expected_total_tau = sum(real(tau_species_dust(icell,5,:),c_double)) + &
          real(tau_dust(icell,5),c_double)
     expected_hhe_tau = sum(real(tau_species_dust(icell,5,:),c_double))
     expected_hhe_target = real(raw_group(icell,5),c_double)*expected_hhe_tau/expected_total_tau
     expected_dust = real(raw_group(icell,5),c_double) - expected_hhe_target
     expected_scale = max(1.0d0,abs(real(raw_group(icell,5),c_double)))
     if (abs(real(absorbed_dust_group(icell,5),c_double)-expected_dust) > &
          5.0d-5*expected_scale) error stop 30
     if (abs(real(returned_group(icell,5),c_double)) > 5.0d-5*expected_scale) error stop 31
     do idir = 1, 3
        if (abs(real(absorbed_hhe_species(icell,5,idir),c_double) - &
             real(raw_group(icell,5),c_double)*real(tau_species_dust(icell,5,idir),c_double)/ &
             expected_total_tau) > 5.0d-5*expected_scale) error stop 32
     end do
  end do
  write(*,'(a)') 'SNRT_CUDA_MULTIGROUP_SPECIES_DUST_NON_SATURATING_OK'

  ! Zero-dust regression: retain the legacy inputs and compare the state,
  ! scalar/group absorption, and inventories bit-for-bit.  The new ledgers
  ! are allowed to expose only the same FP32 partition rounding.
  state_zero_reference = 1.0_c_float
  state_dust = 1.0_c_float
  tau_species_zero_reference = tau_species
  tau_dust_zero_reference = 0.0_c_float
  tau_zero_reference = sum(tau_species_zero_reference,dim=3)
  available_zero_reference(:,1) = 0.050_c_float
  available_zero_reference(:,2) = 0.200_c_float
  available_zero_reference(:,3) = 0.350_c_float
  available_species_dust = available_zero_reference
  ierr = snrt_cuda_multigroup_rt_step_species(state_zero_reference, direction, neighbor, &
       tau_zero_reference, tau_species_zero_reference, available_zero_reference, &
       absorbed_zero_reference, absorbed_group_zero_reference, int(ncell,c_int), &
       int(ncell,c_int), int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 24
  ierr = snrt_cuda_multigroup_rt_step_species_dust(state_dust, direction, neighbor, &
       tau_zero_reference, tau_species_zero_reference, tau_dust_zero_reference, &
       available_species_dust, absorbed_hhe_species, absorbed_dust_group, returned_group, &
       raw_group, assigned_group, absorbed, int(ncell,c_int), int(ncell,c_int), &
       int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr /= 0_c_int) error stop 25
  if (any(transfer(state_dust,state_bits) /= transfer(state_zero_reference,state_zero_bits))) &
       error stop 26
  if (any(transfer(available_species_dust,available_bits) /= &
       transfer(available_zero_reference,available_zero_bits))) error stop 27
  if (any(transfer(absorbed,scalar_bits) /= transfer(absorbed_zero_reference,scalar_zero_bits))) &
       error stop 28
  if (any(transfer(assigned_group,group_bits) /= &
       transfer(absorbed_group_zero_reference,group_zero_bits))) error stop 29
  if (maxval(absorbed_dust_group) /= 0.0_c_float) error stop 30
  hhe_zero_error = maxval(abs(real(sum(absorbed_hhe_species,dim=3) - &
       absorbed_group_zero_reference,c_double)))
  if (hhe_zero_error > 3.0d-5) error stop 30
  zero_error = maxval(abs(real(raw_group - sum(absorbed_hhe_species,dim=3) - &
       returned_group,c_double)))
  if (zero_error > 3.0d-5) error stop 31
  write(*,'(a,es14.6,a,es14.6)') &
       'SNRT_CUDA_MULTIGROUP_SPECIES_DUST_ZERO_DUST_BITWISE_OK hhe_max_abs=', &
       hhe_zero_error, ' closure_max_abs=', zero_error

  ! Invalid-input transaction: validation must reject before copying any
  ! state, inventory, or output back to the caller.
  state_dust = 1.0_c_float
  available_species_dust = available_zero_reference
  state_zero_reference = state_dust
  available_zero_reference = available_species_dust
  tau_dust_zero_reference = 0.0_c_float
  tau_dust_zero_reference(1,1) = -1.0_c_float
  absorbed_hhe_species = -7.0_c_float
  absorbed_dust_group = -7.0_c_float
  returned_group = -7.0_c_float
  raw_group = -7.0_c_float
  assigned_group = -7.0_c_float
  absorbed = -7.0_c_float
  ierr = snrt_cuda_multigroup_rt_step_species_dust(state_dust, direction, neighbor, &
       tau_zero_reference, tau_species_zero_reference, tau_dust_zero_reference, &
       available_species_dust, absorbed_hhe_species, absorbed_dust_group, returned_group, &
       raw_group, assigned_group, absorbed, int(ncell,c_int), int(ncell,c_int), &
       int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr == 0_c_int) error stop 32
  if (any(transfer(state_dust,state_bits) /= transfer(state_zero_reference,state_zero_bits))) &
       error stop 33
  if (any(transfer(available_species_dust,available_bits) /= &
       transfer(available_zero_reference,available_zero_bits))) error stop 34
  if (maxval(absorbed_hhe_species) /= -7.0_c_float .or. maxval(absorbed_dust_group) /= -7.0_c_float .or. &
       maxval(returned_group) /= -7.0_c_float .or. maxval(raw_group) /= -7.0_c_float .or. &
       maxval(assigned_group) /= -7.0_c_float .or. maxval(absorbed) /= -7.0_c_float) error stop 35
  ! A separate all-nonnegative component mismatch must fail for the
  ! conservation check, not merely because of the sign check above.
  tau_dust_zero_reference = 0.0_c_float
  tau_zero_reference(1,1) = tau_zero_reference(1,1) + 0.10_c_float
  ierr = snrt_cuda_multigroup_rt_step_species_dust(state_dust, direction, neighbor, &
       tau_zero_reference, tau_species_zero_reference, tau_dust_zero_reference, &
       available_species_dust, absorbed_hhe_species, absorbed_dust_group, returned_group, &
       raw_group, assigned_group, absorbed, int(ncell,c_int), int(ncell,c_int), &
       int(ndirection,c_int), int(ngroup,c_int), 0.10_c_float)
  if (ierr == 0_c_int) error stop 36
  if (any(transfer(state_dust,state_bits) /= transfer(state_zero_reference,state_zero_bits))) &
       error stop 37
  if (any(transfer(available_species_dust,available_bits) /= &
       transfer(available_zero_reference,available_zero_bits))) error stop 38
  if (maxval(absorbed_hhe_species) /= -7.0_c_float .or. maxval(absorbed_dust_group) /= -7.0_c_float .or. &
       maxval(returned_group) /= -7.0_c_float .or. maxval(raw_group) /= -7.0_c_float .or. &
       maxval(assigned_group) /= -7.0_c_float .or. maxval(absorbed) /= -7.0_c_float) error stop 39
  write(*,'(a)') 'SNRT_CUDA_MULTIGROUP_SPECIES_DUST_INVALID_INPUT_OK'
end program snrt_cuda_multigroup_smoke

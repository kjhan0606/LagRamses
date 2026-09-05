program snrt_cuda_multigroup_smoke
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use snrt_cuda_multigroup_interface, only: snrt_cuda_multigroup_rt_step_species
  implicit none

  integer, parameter :: nx = 7, ny = 5, nz = 3, ncell = nx * ny * nz
  ! Exercise the production nine-group CUDA ABI, including the hard-X-ray
  ! group that was absent from the old four-group smoke.
  integer, parameter :: ndirection = 80, ngroup = 9
  integer :: ix, iy, iz, icell, idir, igroup
  integer(c_int) :: ierr, neighbor(6,ncell)
  real(c_float) :: state(ncell,ndirection,ngroup)
  real(c_float) :: direction(3,ndirection), tau(ncell,ngroup)
  real(c_float) :: tau_species(ncell,ngroup,3)
  real(c_float) :: available_species(ncell,3), initial_species(ncell,3)
  real(c_float) :: absorbed(ncell)
  real(c_float) :: absorbed_group(ncell,ngroup)
  real(c_double) :: before, after, removed, budget_error

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
end program snrt_cuda_multigroup_smoke

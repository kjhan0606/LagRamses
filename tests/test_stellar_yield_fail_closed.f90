program test_stellar_yield_fail_closed
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, channel_snii
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_yield_interpolation, only: interpolate_yield_row, &
       interpolation_ok, interpolation_err_grid, interpolation_err_argument, &
       interpolation_err_nonfinite
  implicit none

  type(stellar_yield_table_t) :: table
  real(stellar_dp) :: returned_mass, remnant_mass, energy
  real(stellar_dp) :: momentum(3), ejecta(n_stellar_elements)
  real(stellar_dp) :: net_yield(n_stellar_elements), nan_value
  integer :: imass, iz, iage, row, ierr

  table%n_rows = 8
  allocate(table%channel(8), table%initial_mass(8), &
       table%birth_metallicity(8), table%age_gyr(8), &
       table%returned_mass(8), table%remnant_mass(8), table%energy(8), &
       table%momentum(8,3), table%ejected_mass(8,n_stellar_elements), &
       table%net_yield(8,n_stellar_elements))
  row = 0
  do imass = 1, 2
     do iz = 1, 2
        do iage = 1, 2
           row = row + 1
           table%channel(row) = channel_snii
           table%initial_mass(row) = real(imass, stellar_dp)
           table%birth_metallicity(row) = 0.02_stellar_dp * real(iz-1, stellar_dp)
           table%age_gyr(row) = real(iage-1, stellar_dp)
           table%returned_mass(row) = real(row, stellar_dp)
           table%remnant_mass(row) = 0.0_stellar_dp
           table%energy(row) = real(row, stellar_dp)
           table%momentum(row,:) = 0.0_stellar_dp
           table%ejected_mass(row,:) = 0.0_stellar_dp
           table%ejected_mass(row,1) = table%returned_mass(row)
           table%net_yield(row,:) = 0.0_stellar_dp
        end do
     end do
  end do
  table%loaded = .true.

  call evaluate(1.5_stellar_dp, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_ok) error stop 1
  call evaluate(0.5_stellar_dp, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_grid) error stop 2
  call evaluate(2.5_stellar_dp, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_grid) error stop 3
  call evaluate(1.5_stellar_dp, 0.03_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_grid) error stop 4
  call evaluate(1.5_stellar_dp, 0.01_stellar_dp, 1.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_grid) error stop 5
  nan_value = ieee_value(0.0_stellar_dp, ieee_quiet_nan)
  call evaluate(nan_value, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_argument) error stop 6
  call evaluate(1.5_stellar_dp, nan_value, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_argument) error stop 7
  call evaluate(1.5_stellar_dp, 0.01_stellar_dp, nan_value, ierr)
  if (ierr /= interpolation_err_argument) error stop 8
  call evaluate(-1.0_stellar_dp, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_argument) error stop 9
  call evaluate(1.5_stellar_dp, -0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_argument) error stop 10
  call evaluate(1.5_stellar_dp, 0.01_stellar_dp, -0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_argument) error stop 11
  table%energy(1) = nan_value
  call evaluate(1.5_stellar_dp, 0.01_stellar_dp, 0.5_stellar_dp, ierr)
  if (ierr /= interpolation_err_nonfinite) error stop 12

  write(*,'(A)') 'stellar yield fail-closed policy: PASS'

contains

  subroutine evaluate(mass, metallicity, age_gyr, status)
    real(stellar_dp), intent(in) :: mass, metallicity, age_gyr
    integer, intent(out) :: status

    call interpolate_yield_row(table, channel_snii, mass, metallicity, &
         age_gyr, returned_mass, remnant_mass, energy, momentum, ejecta, &
         net_yield, status)
  end subroutine evaluate

end program test_stellar_yield_fail_closed

program test_stellar_yield_audit_contract
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, channel_owns_terminal_remnant
  use stellar_yield_tables, only: stellar_yield_table_t, clear_yield_table
  use stellar_yield_audit, only: audit_yield_table, yield_audit_ok, &
       yield_audit_err_value, yield_audit_err_mass, &
       yield_audit_err_duplicate, yield_audit_err_grid, &
       yield_audit_err_energy_monotonic, &
       yield_audit_err_remnant_ownership
  implicit none

  type(stellar_yield_table_t) :: table
  integer :: ierr

  call make_valid_table(table)
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (ierr /= yield_audit_ok) error stop 1

  ! The canonical contract permits untracked residual ejecta.
  if (sum(table%ejected_mass(2,:)) >= table%returned_mass(2)) error stop 2

  table%ejected_mass(2,1) = 2.0_stellar_dp * table%returned_mass(2)
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_mass) == 0) error stop 3

  call make_valid_table(table)
  call copy_row(table, 18, 19)
  table%n_rows = 19
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_duplicate) == 0) error stop 4

  call make_valid_table(table)
  table%n_rows = 17
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_grid) == 0) error stop 5

  call make_valid_table(table)
  table%returned_mass(1) = 0.01_stellar_dp
  table%ejected_mass(1,1) = 0.005_stellar_dp
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_value) == 0) error stop 6

  call make_valid_table(table)
  table%energy(2) = 2.0_stellar_dp
  table%energy(3) = 1.0_stellar_dp
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_energy_monotonic) == 0) error stop 7

  call make_valid_table(table)
  table%remnant_mass(2) = 0.01_stellar_dp
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_remnant_ownership) == 0) error stop 8

  call make_valid_table(table)
  where (table%channel(1:table%n_rows) == 1)
     table%age_gyr(1:table%n_rows) = table%age_gyr(1:table%n_rows) + &
          1.0_stellar_dp
  end where
  call audit_yield_table(table, 1.0e-10_stellar_dp, ierr, .true., &
       channel_owns_terminal_remnant)
  if (iand(ierr, yield_audit_err_grid) == 0) error stop 9

  call clear_yield_table(table)
  write(*,'(A)') 'stellar yield audit contract: PASS'

contains

  subroutine make_valid_table(table)
    type(stellar_yield_table_t), intent(inout) :: table
    integer :: channel, mass_index, age_index, row
    real(stellar_dp) :: age, returned

    call clear_yield_table(table)
    allocate(table%channel(19), table%initial_mass(19), &
         table%birth_metallicity(19), table%age_gyr(19), &
         table%returned_mass(19), table%remnant_mass(19), table%energy(19), &
         table%momentum(19,3), table%ejected_mass(19,n_stellar_elements), &
         table%net_yield(19,n_stellar_elements))
    table%n_rows = 18
    row = 0
    do channel = 1, 3
       do mass_index = 1, 2
          do age_index = 0, 2
             row = row + 1
             age = real(age_index, stellar_dp)
             returned = age * 0.05_stellar_dp * real(channel, stellar_dp)
             table%channel(row) = channel
             table%initial_mass(row) = real(mass_index, stellar_dp)
             table%birth_metallicity(row) = 0.01_stellar_dp
             table%age_gyr(row) = age
             table%returned_mass(row) = returned
             table%remnant_mass(row) = 0.0_stellar_dp
             if (channel > 1) then
                table%remnant_mass(row) = age * 0.01_stellar_dp
             end if
             table%energy(row) = age * real(channel, stellar_dp)
             table%momentum(row,:) = 0.0_stellar_dp
             table%ejected_mass(row,:) = 0.0_stellar_dp
             table%ejected_mass(row,1) = 0.5_stellar_dp * returned
             table%net_yield(row,:) = 0.0_stellar_dp
          end do
       end do
    end do
    table%loaded = .true.
  end subroutine make_valid_table

  subroutine copy_row(table, source, target)
    type(stellar_yield_table_t), intent(inout) :: table
    integer, intent(in) :: source, target

    table%channel(target) = table%channel(source)
    table%initial_mass(target) = table%initial_mass(source)
    table%birth_metallicity(target) = table%birth_metallicity(source)
    table%age_gyr(target) = table%age_gyr(source)
    table%returned_mass(target) = table%returned_mass(source)
    table%remnant_mass(target) = table%remnant_mass(source)
    table%energy(target) = table%energy(source)
    table%momentum(target,:) = table%momentum(source,:)
    table%ejected_mass(target,:) = table%ejected_mass(source,:)
    table%net_yield(target,:) = table%net_yield(source,:)
  end subroutine copy_row

end program test_stellar_yield_audit_contract

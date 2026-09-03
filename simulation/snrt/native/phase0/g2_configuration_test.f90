program g2_configuration_test
  use stellar_enrichment_config, only: stellar_dp, set_enrichment_defaults, &
       read_enrichment_namelist, default_imf_id, population_model_id, &
       population_binary_ssp, configured_channel_mass_min, &
       yield_source_basis_id, yield_basis_per_star_cumulative, &
       configured_binary_fraction
  implicit none

  integer :: unit, ios, failures
  character(len=256) :: filename

  failures = 0
  filename = 'g2_configuration_test.nml'
  call set_enrichment_defaults()
  open(newunit=unit, file=filename, status='replace', action='write', &
       iostat=ios)
  call expect(ios == 0, 'configuration fixture opens for writing', failures)
  if (ios == 0) then
     write(unit, '(a)') '&stellar_enrichment_params'
     write(unit, '(a)') " feedback_mode='channel_resolved',"
     write(unit, '(a)') ' imf_id=2,'
     write(unit, '(a)') " population_model='binary_ssp',"
     write(unit, '(a)') " yield_source_basis='per_star_cumulative',"
     write(unit, '(a)') ' imf_mass_min_msun=0.08, imf_mass_max_msun=120.0,'
     write(unit, '(a)') ' binary_fraction=0.5,'
     write(unit, '(a)') ' channel_mass_min_msun=1.2, 1.0, 8.0, 3.0, 140.0,'
     write(unit, '(a)') ' channel_mass_max_msun=120.0, 8.0, 40.0, 8.0, 260.0,'
     write(unit, '(a)') '/'
     close(unit)
  end if

  open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
  call expect(ios == 0, 'configuration fixture opens for reading', failures)
  if (ios == 0) then
     call read_enrichment_namelist(unit, ios)
     close(unit)
  end if
  call expect(ios == 0, 'IMF/population/boundary values parse from namelist', &
       failures)
  call expect(default_imf_id == 2, 'runtime IMF is configuration-driven', failures)
  call expect(population_model_id == population_binary_ssp, &
       'binary population model is configuration-driven', failures)
  call expect(yield_source_basis_id == yield_basis_per_star_cumulative, &
       'yield basis is configuration-driven', failures)
  call expect(abs(configured_binary_fraction-0.5_stellar_dp) < 1.0e-12, &
       'binary fraction is configuration-driven', failures)
  call expect(abs(configured_channel_mass_min(1) - 1.2_stellar_dp) < 1.0e-12, &
       'channel mass lower bound is configuration-driven', failures)

  call set_enrichment_defaults()
  open(newunit=unit, file=filename, status='replace', action='write', &
       iostat=ios)
  if (ios == 0) then
     write(unit, '(a)') '&stellar_enrichment_params'
     write(unit, '(a)') ' imf_id=99,'
     write(unit, '(a)') " feedback_mode='channel_resolved',"
     write(unit, '(a)') " population_model='single_star_ssp',"
     write(unit, '(a)') " yield_source_basis='per_star_cumulative',"
     write(unit, '(a)') ' imf_mass_min_msun=0.08, imf_mass_max_msun=120.0,'
     write(unit, '(a)') ' binary_fraction=0.0,'
     write(unit, '(a)') ' channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0,'
     write(unit, '(a)') ' channel_mass_max_msun=120.0,8.0,40.0,8.0,260.0,'
     write(unit, '(a)') '/'
     close(unit)
  end if
  open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
  if (ios == 0) then
     call read_enrichment_namelist(unit, ios)
     close(unit)
  end if
  call expect(ios /= 0, 'invalid IMF identifier is rejected', failures)
  call set_enrichment_defaults()

  if (failures == 0) then
     write(*, '(a)') 'G2_CONFIGURATION_TEST_OK'
  else
     write(*, '(a,i0)') 'G2_CONFIGURATION_TEST_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (condition) then
       write(*, '(a)') 'PASS: ' // trim(label)
    else
       failures = failures + 1
       write(*, '(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine expect

end program g2_configuration_test

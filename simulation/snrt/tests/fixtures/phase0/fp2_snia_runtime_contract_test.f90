program fp2_snia_runtime_contract_test
  use stellar_enrichment_config, only: stellar_dp, stellar_imf_kroupa, stellar_imf_chabrier, &
       population_binary_ssp, population_single_star_ssp
  use stellar_snia_population_contract, only: &
       snia_population_realization_t, read_snia_population_realization_namelist, &
       snia_population_contract_ok, validate_snia_population_binding, validate_snia_population_realization, &
       snia_accounting_strict_wd, snia_accounting_effective_ssp, snia_effective_ssp_approval
  use stellar_snia_physical_contract, only: snia_physical_contract_t, &
       read_snia_physical_contract_namelist, snia_contract_ok
  use stellar_snia_cell_deposition, only: snia_thermal_coupling_t, &
       read_snia_thermal_coupling_namelist, snia_deposition_ok
  implicit none

  type(snia_population_realization_t) :: population
  type(snia_physical_contract_t) :: physical
  type(snia_thermal_coupling_t) :: coupling
  character(len=1024) :: filename
  integer :: unit, status, ierr, failures

  failures = 0
  call get_environment_variable('PHASE0_SNIA_RUNTIME_CONTRACT', filename, &
       status=status)
  if (status /= 0 .or. len_trim(filename) == 0) error stop 1
  open(newunit=unit, file=trim(filename), status='old', action='read', &
       iostat=status)
  if (status /= 0) error stop 2

  call read_snia_population_realization_namelist(unit, population, ierr)
  call expect(ierr == snia_population_contract_ok, &
       'approved population group loads from the ordered handoff', failures)
  call read_snia_physical_contract_namelist(unit, physical, ierr)
  call expect(ierr == snia_contract_ok, &
       'approved physical event group loads from the ordered handoff', failures)
  call read_snia_thermal_coupling_namelist(unit, coupling, ierr)
  close(unit)
  call expect(ierr == snia_deposition_ok, &
       'approved thermal coupling group loads from the ordered handoff', failures)
  call expect(trim(population%source_commit_binding) == &
       trim(physical%source_commit_binding) .and. &
       trim(population%source_commit_binding) == trim(coupling%source_commit_binding), &
       'all runtime groups bind to the same source commit', failures)
  call expect(trim(population%approval_id) == trim(physical%approval_id) .and. &
       trim(population%approval_id) == trim(coupling%approval_id), &
       'all runtime groups bind to the same named approval', failures)
  call expect_close(physical%returned_mass_per_event, &
       1.4004633930489443d0, 'approved returned mass is preserved', failures)
  call expect_close(physical%energy_per_event, &
       1.5063100005966762d51, 'approved event energy is preserved', failures)
  call expect(coupling%mode == 1 .and. coupling%thermal_fraction == &
       1.0_stellar_dp .and. .not. coupling%include_event_momentum_kinetic, &
       'thermal coupling is the approved total-energy policy', failures)

  call validate_snia_population_binding(population, stellar_imf_kroupa, population_binary_ssp, .5d0, ierr)
  call expect(ierr == snia_population_contract_ok, 'approved runtime population matches', failures)
  call validate_snia_population_binding(population, stellar_imf_chabrier, population_binary_ssp, .5d0, ierr)
  call expect(ierr /= snia_population_contract_ok, 'Kroupa DTD cannot silently feed Chabrier particles', failures)
  call validate_snia_population_binding(population, stellar_imf_kroupa, population_single_star_ssp, .5d0, ierr)
  call expect(ierr /= snia_population_contract_ok, 'binary DTD cannot silently feed a single-star population', failures)
  call validate_snia_population_binding(population, stellar_imf_kroupa, population_binary_ssp, .3d0, ierr)
  call expect(ierr /= snia_population_contract_ok, 'baked-in binary metadata still has to match', failures)
  population%imf_conversion_factor=2d0
  call validate_snia_population_binding(population, stellar_imf_chabrier, population_binary_ssp, .5d0, ierr)
  call expect(ierr /= snia_population_contract_ok, 'conversion factor does not override target IMF identity', failures)

  population%mass_accounting=snia_accounting_effective_ssp
  population%accounting_approval_id=''
  call validate_snia_population_realization(population,ierr)
  call expect(ierr /= snia_population_contract_ok, 'effective SSP requires separate accounting approval', failures)
  population%accounting_approval_id=snia_effective_ssp_approval
  call validate_snia_population_realization(population,ierr)
  call expect(ierr == snia_population_contract_ok, 'explicit effective SSP accounting admitted', failures)
  population%mass_accounting='unknown'
  call validate_snia_population_realization(population,ierr)
  call expect(ierr /= snia_population_contract_ok, 'unknown mass accounting rejects', failures)
  population%mass_accounting=snia_accounting_strict_wd
  call validate_snia_population_realization(population,ierr)
  call expect(ierr /= snia_population_contract_ok, 'effective approval cannot masquerade as strict WD', failures)
  population%accounting_approval_id=''
  call validate_snia_population_realization(population,ierr)
  call expect(ierr == snia_population_contract_ok, 'historical strict WD accounting preserved', failures)

  if (failures > 0) then
     write(*, '(a,i0)') 'FP2_SNIa_RUNTIME_CONTRACT_TEST_FAILED failures=', failures
     error stop 3
  end if
  write(*, '(a)') 'FP2_SNIa_RUNTIME_CONTRACT_TEST_OK'

contains

  subroutine expect(condition, label, failure_count)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failure_count

    if (condition) then
       write(*, '(a)') 'PASS: '//trim(label)
    else
       write(*, '(a)') 'FAIL: '//trim(label)
       failure_count = failure_count + 1
    end if
  end subroutine expect

  subroutine expect_close(actual, expected, label, failure_count)
    real(stellar_dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failure_count
    real(stellar_dp) :: scale

    scale = max(1.0e-300_stellar_dp, abs(actual), abs(expected))
    call expect(abs(actual - expected) <= 1.0e-12_stellar_dp * scale, &
         label, failure_count)
  end subroutine expect_close

end program fp2_snia_runtime_contract_test

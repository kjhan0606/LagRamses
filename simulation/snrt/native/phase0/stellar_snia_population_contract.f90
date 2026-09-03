! Approved-baseline SNIa binary-population realization contract.
!
! This module defines the inputs needed to turn the interval DTD kernel into
! an expected event count.  The selected Maoz field DTD baseline is loaded by
! the versioned sidecar/namelist contract and binds its choices to a source
! identifier, immutable source commit, and named approval before a runtime
! caller may use the interface.

module stellar_snia_population_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, population_binary_ssp, &
       stellar_imf_salpeter, stellar_imf_popiii
  use stellar_snia_dtd, only: integrate_snia_dtd_interval, snia_dtd_ok
  implicit none

  private

  integer, parameter, public :: snia_population_contract_ok = 0
  integer, parameter, public :: snia_population_contract_err_argument = 1
  integer, parameter, public :: snia_population_contract_err_unapproved = 2
  integer, parameter, public :: snia_population_contract_err_model = 4
  integer, parameter, public :: snia_population_contract_err_parameter = 8
  integer, parameter, public :: snia_population_contract_err_realization = 16

  integer, parameter, public :: snia_realization_expectation = 1
  integer, parameter, public :: snia_realization_poisson = 2
  integer, parameter, public :: snia_binary_fraction_baked_into_rate = 1
  integer, parameter, public :: snia_binary_fraction_scales_rate = 2
  integer, parameter, public :: snia_metallicity_factor_supplied = 1

  type, public :: snia_population_realization_t
     logical :: approved = .false.
     character(len=128) :: population_source_id = ''
     integer :: population_model_id = -1
     integer :: imf_id = -1
     real(stellar_dp) :: binary_fraction = -1.0_stellar_dp
     real(stellar_dp) :: imf_conversion_factor = -1.0_stellar_dp
     real(stellar_dp) :: minimum_delay_gyr = -1.0_stellar_dp
     real(stellar_dp) :: maximum_delay_gyr = -1.0_stellar_dp
     real(stellar_dp) :: power_law_index = -1.0_stellar_dp
     real(stellar_dp) :: events_per_initial_msun = -1.0_stellar_dp
     integer :: event_realization_policy = 0
     integer :: binary_fraction_policy = 0
     integer :: metallicity_policy = 0
     character(len=128) :: metallicity_factor_source_id = ''
     character(len=128) :: source_commit_binding = ''
     character(len=128) :: approval_id = ''
  end type snia_population_realization_t

  public :: validate_snia_population_realization
  public :: read_snia_population_realization_namelist
  public :: evaluate_snia_interval_events

contains

  subroutine validate_snia_population_realization(realization, ierr)
    type(snia_population_realization_t), intent(in) :: realization
    integer, intent(out) :: ierr

    ierr = snia_population_contract_ok
    if (.not. realization%approved) then
       ierr = snia_population_contract_err_unapproved
       return
    end if
    if (len_trim(realization%population_source_id) == 0 .or. &
         realization%population_model_id /= population_binary_ssp .or. &
         realization%imf_id < stellar_imf_salpeter .or. &
         realization%imf_id > stellar_imf_popiii) then
       ierr = snia_population_contract_err_model
       return
    end if
    if (.not. ieee_is_finite(realization%binary_fraction) .or. &
         .not. ieee_is_finite(realization%imf_conversion_factor) .or. &
         .not. ieee_is_finite(realization%minimum_delay_gyr) .or. &
         .not. ieee_is_finite(realization%maximum_delay_gyr) .or. &
         .not. ieee_is_finite(realization%power_law_index) .or. &
         .not. ieee_is_finite(realization%events_per_initial_msun) .or. &
         realization%binary_fraction < 0.0_stellar_dp .or. &
         realization%binary_fraction > 1.0_stellar_dp .or. &
         realization%imf_conversion_factor <= 0.0_stellar_dp .or. &
         realization%minimum_delay_gyr <= 0.0_stellar_dp .or. &
         realization%maximum_delay_gyr <= realization%minimum_delay_gyr .or. &
         realization%events_per_initial_msun < 0.0_stellar_dp) then
       ierr = snia_population_contract_err_parameter
       return
    end if
    if (realization%event_realization_policy /= snia_realization_expectation .and. &
         realization%event_realization_policy /= snia_realization_poisson) then
       ierr = snia_population_contract_err_realization
       return
    end if
    if (realization%binary_fraction_policy /= snia_binary_fraction_baked_into_rate .and. &
         realization%binary_fraction_policy /= snia_binary_fraction_scales_rate) then
       ierr = snia_population_contract_err_parameter
       return
    end if
    if (realization%metallicity_policy /= snia_metallicity_factor_supplied .or. &
         .not. is_hex_commit(realization%source_commit_binding) .or. &
         len_trim(realization%approval_id) == 0 .or. &
         len_trim(realization%metallicity_factor_source_id) == 0) then
       ierr = snia_population_contract_err_parameter
    end if
  end subroutine validate_snia_population_realization

  subroutine read_snia_population_realization_namelist(iunit, realization, ierr)
    ! JSON remains the review/provenance sidecar.  A production caller must
    ! pass its separately generated, human-auditable namelist through this
    ! loader so the complete record is populated into the Fortran type before
    ! validation; no field is silently defaulted during the handoff.
    integer, intent(in) :: iunit
    type(snia_population_realization_t), intent(out) :: realization
    integer, intent(out) :: ierr

    logical :: approved
    character(len=128) :: population_source_id, metallicity_factor_source_id
    character(len=128) :: source_commit_binding, approval_id
    integer :: population_model_id, imf_id, event_realization_policy
    integer :: binary_fraction_policy, metallicity_policy, read_ierr
    real(stellar_dp) :: binary_fraction, imf_conversion_factor
    real(stellar_dp) :: minimum_delay_gyr, maximum_delay_gyr
    real(stellar_dp) :: power_law_index, events_per_initial_msun

    namelist /snia_population_realization/ approved, population_source_id, &
         population_model_id, imf_id, binary_fraction, binary_fraction_policy, &
         imf_conversion_factor, minimum_delay_gyr, maximum_delay_gyr, &
         power_law_index, events_per_initial_msun, event_realization_policy, &
         metallicity_policy, metallicity_factor_source_id, &
         source_commit_binding, approval_id

    approved = .false.
    population_source_id = ''
    population_model_id = -1
    imf_id = -1
    binary_fraction = -1.0_stellar_dp
    binary_fraction_policy = 0
    imf_conversion_factor = -1.0_stellar_dp
    minimum_delay_gyr = -1.0_stellar_dp
    maximum_delay_gyr = -1.0_stellar_dp
    power_law_index = -1.0_stellar_dp
    events_per_initial_msun = -1.0_stellar_dp
    event_realization_policy = 0
    metallicity_policy = 0
    metallicity_factor_source_id = ''
    source_commit_binding = ''
    approval_id = ''

    read(iunit, nml=snia_population_realization, iostat=read_ierr)
    realization%approved = approved
    realization%population_source_id = population_source_id
    realization%population_model_id = population_model_id
    realization%imf_id = imf_id
    realization%binary_fraction = binary_fraction
    realization%binary_fraction_policy = binary_fraction_policy
    realization%imf_conversion_factor = imf_conversion_factor
    realization%minimum_delay_gyr = minimum_delay_gyr
    realization%maximum_delay_gyr = maximum_delay_gyr
    realization%power_law_index = power_law_index
    realization%events_per_initial_msun = events_per_initial_msun
    realization%event_realization_policy = event_realization_policy
    realization%metallicity_policy = metallicity_policy
    realization%metallicity_factor_source_id = metallicity_factor_source_id
    realization%source_commit_binding = source_commit_binding
    realization%approval_id = approval_id
    if (read_ierr /= 0) then
       ierr = snia_population_contract_err_argument
       return
    end if
    call validate_snia_population_realization(realization, ierr)
  end subroutine read_snia_population_realization_namelist

  subroutine evaluate_snia_interval_events(realization, initial_mass_msun, &
       age_old_gyr, age_new_gyr, metallicity_factor, expected_events, ierr)
    type(snia_population_realization_t), intent(in) :: realization
    real(stellar_dp), intent(in) :: initial_mass_msun
    real(stellar_dp), intent(in) :: age_old_gyr, age_new_gyr
    real(stellar_dp), intent(in) :: metallicity_factor
    real(stellar_dp), intent(out) :: expected_events
    integer, intent(out) :: ierr

    real(stellar_dp) :: events_per_mass
    integer :: contract_ierr, dtd_ierr

    expected_events = 0.0_stellar_dp
    ierr = snia_population_contract_ok
    if (.not. ieee_is_finite(initial_mass_msun) .or. &
         .not. ieee_is_finite(metallicity_factor) .or. &
         initial_mass_msun < 0.0_stellar_dp .or. &
         metallicity_factor < 0.0_stellar_dp) then
       ierr = snia_population_contract_err_argument
       return
    end if
    call validate_snia_population_realization(realization, contract_ierr)
    if (contract_ierr /= snia_population_contract_ok) then
       ierr = contract_ierr
       return
    end if
    if (realization%event_realization_policy /= snia_realization_expectation) then
       ! The interval kernel returns an expectation.  Poisson realization
       ! needs an explicit seeded RNG contract and is not silently sampled.
       ierr = snia_population_contract_err_realization
       return
    end if
    events_per_mass = realization%events_per_initial_msun * &
         realization%imf_conversion_factor * metallicity_factor
    if (realization%binary_fraction_policy == snia_binary_fraction_scales_rate) then
       events_per_mass = events_per_mass * realization%binary_fraction
    end if
    call integrate_snia_dtd_interval(age_old_gyr, age_new_gyr, &
         realization%minimum_delay_gyr, realization%maximum_delay_gyr, &
         realization%power_law_index, events_per_mass, expected_events, dtd_ierr)
    if (dtd_ierr /= snia_dtd_ok) then
       expected_events = 0.0_stellar_dp
       ierr = snia_population_contract_err_parameter
       return
    end if
    expected_events = initial_mass_msun * expected_events
    if (.not. ieee_is_finite(expected_events) .or. expected_events < 0.0_stellar_dp) then
       expected_events = 0.0_stellar_dp
       ierr = snia_population_contract_err_parameter
    end if
  end subroutine evaluate_snia_interval_events

  logical function is_hex_commit(value)
    character(len=*), intent(in) :: value
    integer :: i, code

    is_hex_commit = len_trim(value) == 40
    if (.not. is_hex_commit) return
    do i = 1, 40
       code = iachar(value(i:i))
       if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
            (code >= iachar('a') .and. code <= iachar('f')) .or. &
            (code >= iachar('A') .and. code <= iachar('F')))) then
          is_hex_commit = .false.
          return
       end if
    end do
  end function is_hex_commit

end module stellar_snia_population_contract

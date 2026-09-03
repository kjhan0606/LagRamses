! Review-only SNIa binary-population realization contract.
!
! This module defines the inputs needed to turn the interval DTD kernel into
! an expected event count.  It deliberately contains no project-selected
! population, normalization, or metallicity table.  An approved record must
! bind those choices to a source identifier, an immutable source commit, and a
! named approval before this interface can be used by a runtime caller.

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
     integer :: metallicity_policy = 0
     character(len=128) :: source_commit_binding = ''
     character(len=128) :: approval_id = ''
  end type snia_population_realization_t

  public :: validate_snia_population_realization
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
    if (realization%metallicity_policy /= snia_metallicity_factor_supplied .or. &
         .not. is_hex_commit(realization%source_commit_binding) .or. &
         len_trim(realization%approval_id) == 0) then
       ierr = snia_population_contract_err_parameter
    end if
  end subroutine validate_snia_population_realization

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

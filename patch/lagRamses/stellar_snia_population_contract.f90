! Approved-baseline SNIa binary-population realization contract.
!
! This module defines the inputs needed to turn the interval DTD kernel into
! an expected event count.  The selected Maoz field DTD baseline is loaded by
! the versioned sidecar/namelist contract and binds its choices to a source
! identifier, immutable source commit, and named approval before a runtime
! caller may use the interface.

module stellar_snia_population_contract
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, population_binary_ssp, population_effective_ssp, &
       stellar_imf_salpeter, stellar_imf_miller_scalo, stellar_imf_kroupa
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
  character(len=*), parameter, public :: snia_accounting_strict_wd = 'strict_wd'
  character(len=*), parameter, public :: snia_accounting_effective_ssp = 'effective_ssp'
  character(len=*), parameter, public :: snia_effective_ssp_approval = 'SNIA-EFFECTIVE-SSP-2026-09-07'
  character(len=*), parameter, public :: snia_event_empirical = 'empirical_powerlaw'
  character(len=*), parameter, public :: snia_event_frozen_he = 'tabulated_frozen_he_contact_v1'
  character(len=*), parameter, public :: snia_frozen_he_closure = 'imposed_contact_stable_eta1_frozen135_v1'
  character(len=*), parameter, public :: snia_frozen_he_approval = 'SNIA-FROZEN-HE-HYBRID-2026-09-11'

  type, public :: snia_frozen_he_event_t
     real(stellar_dp) :: birth_z=0, age_yr=0, weight=0, contact_age_yr=0
     real(stellar_dp) :: seed_wd=0, seed_donor=0, donor_background_loss=0
     real(stellar_dp) :: transfer_rate=0, threshold=0, history_end_age=0, initial_binary_mass=0
     real(stellar_dp) :: retention_reference_z=0, retention_mass_ceiling=0
     integer :: history_bin_num=-1
     character(len=64) :: closure_id='', comparison_approval=''
     character(len=40) :: source_commit=''
     ! Raw BPP, BCM, grid, selection, Params, converter, retention code,
     ! N100 JSON and base runtime sidecar, respectively. Contents and every
     ! interpreted row parameter are also bound by the restart identity.
     character(len=64) :: source_sha256(9)=''
  end type snia_frozen_he_event_t

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
     character(len=32) :: mass_accounting = snia_accounting_strict_wd
     character(len=128) :: accounting_approval_id = ''
     character(len=64) :: event_model = snia_event_empirical
     type(snia_frozen_he_event_t) :: event
  end type snia_population_realization_t

  public :: validate_snia_population_realization
  public :: read_snia_population_realization_namelist
  public :: evaluate_snia_interval_events
  public :: validate_snia_population_binding
  public :: snia_event_model_identity, validate_snia_event_source

contains

  subroutine validate_snia_population_binding(realization, imf_id, population_id, binary_fraction, ierr)
    ! The empirical DTD is per INITIAL mass of its declared population.
    ! A valid standalone handoff must also match the actual particle model.
    ! imf_conversion_factor belongs to that handoff's declared IMF; it must
    ! not silently authorize using a different runtime IMF.
    type(snia_population_realization_t), intent(in) :: realization
    integer, intent(in) :: imf_id, population_id
    real(stellar_dp), intent(in) :: binary_fraction
    integer, intent(out) :: ierr
    call validate_snia_population_realization(realization, ierr)
    if (ierr /= snia_population_contract_ok) return
    ierr = snia_population_contract_err_model
    if (imf_id /= realization%imf_id .or. population_id /= realization%population_model_id) return
    if (.not. ieee_is_finite(binary_fraction)) return
    if (abs(binary_fraction-realization%binary_fraction) > &
         32*epsilon(1.0_stellar_dp)*max(1.0_stellar_dp,abs(realization%binary_fraction))) return
    ierr = snia_population_contract_ok
  end subroutine validate_snia_population_binding

  subroutine validate_snia_population_realization(realization, ierr)
    type(snia_population_realization_t), intent(in) :: realization
    integer, intent(out) :: ierr

    ierr = snia_population_contract_ok
    if (.not. realization%approved) then
       ierr = snia_population_contract_err_unapproved
       return
    end if
    select case (trim(realization%mass_accounting))
    case (snia_accounting_strict_wd)
       if (len_trim(realization%accounting_approval_id) /= 0) then
          ierr = snia_population_contract_err_parameter
          return
       end if
    case (snia_accounting_effective_ssp)
       if (trim(realization%accounting_approval_id) /= snia_effective_ssp_approval) then
          ierr = snia_population_contract_err_unapproved
          return
       end if
    case default
       ierr = snia_population_contract_err_model
       return
    end select
    if (len_trim(realization%population_source_id) == 0 .or. &
         (realization%population_model_id /= population_binary_ssp.and. &
          realization%population_model_id /= population_effective_ssp) .or. &
         realization%imf_id < stellar_imf_salpeter .or. &
         realization%imf_id > stellar_imf_miller_scalo) then
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
    if(ierr/=snia_population_contract_ok)return
    if(realization%population_model_id==population_effective_ssp)then
       ! Observed DTD already includes unresolved binary incidence. Never
       ! multiply by the zero placeholder fraction or invent WD ownership.
       if(realization%mass_accounting/=snia_accounting_effective_ssp.or. &
            realization%event_model/=snia_event_empirical.or.realization%binary_fraction/=0d0.or. &
            realization%binary_fraction_policy/=snia_binary_fraction_baked_into_rate)then
          ierr=snia_population_contract_err_model;return
       endif
    endif
    select case(trim(realization%event_model))
    case(snia_event_empirical)
       ! Historical defaults, normalization and restart representation stay unchanged.
    case(snia_event_frozen_he)
       if(realization%mass_accounting/=snia_accounting_effective_ssp.or. &
            realization%event_realization_policy/=snia_realization_expectation.or. &
            realization%binary_fraction_policy/=snia_binary_fraction_baked_into_rate.or. &
            realization%imf_id/=stellar_imf_kroupa.or.realization%binary_fraction/=0.5_stellar_dp.or. &
            realization%imf_conversion_factor/=1.0_stellar_dp)then
          ierr=snia_population_contract_err_model
          return
       endif
       call validate_frozen_he_event(realization%event,ierr)
       if(ierr/=0)return
       if(abs(realization%events_per_initial_msun-realization%event%weight)> &
            64*epsilon(1d0)*realization%event%weight)ierr=snia_population_contract_err_parameter
    case default
       ierr=snia_population_contract_err_model
    end select
  end subroutine validate_snia_population_realization

  subroutine read_snia_population_realization_namelist(iunit, realization, ierr)
    ! JSON remains the review/provenance sidecar.  A production caller must
    ! pass its separately generated, human-auditable namelist through this
    ! loader so the complete record is populated into the Fortran type before
    ! validation. Only the backward-compatible strict mass-accounting choice
    ! defaults; effective SSP accounting requires its separate approval ID.
    integer, intent(in) :: iunit
    type(snia_population_realization_t), intent(out) :: realization
    integer, intent(out) :: ierr

    logical :: approved
    character(len=128) :: population_source_id, metallicity_factor_source_id
    character(len=128) :: source_commit_binding, approval_id
    character(len=32) :: mass_accounting
    character(len=128) :: accounting_approval_id
    integer :: population_model_id, imf_id, event_realization_policy
    integer :: binary_fraction_policy, metallicity_policy, read_ierr
    real(stellar_dp) :: binary_fraction, imf_conversion_factor
    real(stellar_dp) :: minimum_delay_gyr, maximum_delay_gyr
    real(stellar_dp) :: power_law_index, events_per_initial_msun
    character(len=64) :: event_model
    type(snia_frozen_he_event_t) :: row
    namelist /snia_frozen_he_event/ row

    namelist /snia_population_realization/ approved, population_source_id, &
         population_model_id, imf_id, binary_fraction, binary_fraction_policy, &
         imf_conversion_factor, minimum_delay_gyr, maximum_delay_gyr, &
         power_law_index, events_per_initial_msun, event_realization_policy, &
         metallicity_policy, metallicity_factor_source_id, &
         source_commit_binding, approval_id, mass_accounting, accounting_approval_id, event_model

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
    mass_accounting = snia_accounting_strict_wd
    accounting_approval_id = ''
    event_model=snia_event_empirical

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
    realization%mass_accounting = mass_accounting
    realization%accounting_approval_id = accounting_approval_id
    realization%event_model=event_model
    if (read_ierr /= 0) then
       ierr = snia_population_contract_err_argument
       return
    end if
    if(event_model==snia_event_frozen_he)then
       read(iunit,nml=snia_frozen_he_event,iostat=read_ierr)
       if(read_ierr/=0)then
          ierr=snia_population_contract_err_argument
          return
       endif
       realization%event=row
    endif
    call validate_snia_population_realization(realization, ierr)
  end subroutine read_snia_population_realization_namelist

  subroutine evaluate_snia_interval_events(realization, initial_mass_msun, &
       age_old_gyr, age_new_gyr, metallicity_factor, expected_events, ierr, birth_metallicity)
    type(snia_population_realization_t), intent(in) :: realization
    real(stellar_dp), intent(in) :: initial_mass_msun
    real(stellar_dp), intent(in) :: age_old_gyr, age_new_gyr
    real(stellar_dp), intent(in) :: metallicity_factor
    real(stellar_dp), intent(out) :: expected_events
    integer, intent(out) :: ierr
    ! Separate from the legacy multiplicative factor. Required by the table.
    real(stellar_dp), intent(in), optional :: birth_metallicity

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
    if(realization%event_model==snia_event_frozen_he)then
       ierr=snia_population_contract_err_argument
       if(.not.present(birth_metallicity))return
       if(.not.all(ieee_is_finite([age_old_gyr,age_new_gyr,birth_metallicity])))return
       if(age_old_gyr<0.or.age_new_gyr<age_old_gyr.or.metallicity_factor/=1.0_stellar_dp)return
       if(abs(birth_metallicity-realization%event%birth_z)> &
            64*epsilon(1d0)*realization%event%birth_z)return
       ! Exactly the same CDF for interval and prior-return calls. No age
       ! tolerance at the jump: old<event<=new, so an accepted event cannot replay.
       expected_events=initial_mass_msun*(frozen_he_cdf(realization%event,age_new_gyr)- &
            frozen_he_cdf(realization%event,age_old_gyr))
       if(.not.ieee_is_finite(expected_events))then
          expected_events=0
          return
       endif
       ierr=snia_population_contract_ok
       return
    endif
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

  pure real(stellar_dp) function frozen_he_cdf(event,age_gyr) result(count)
    type(snia_frozen_he_event_t),intent(in)::event
    real(stellar_dp),intent(in)::age_gyr
    count=0
    if(age_gyr>=event%age_yr*1d-9)count=event%weight
  end function frozen_he_cdf

  function frozen_he_numbers(event) result(values)
    type(snia_frozen_he_event_t),intent(in)::event
    real(stellar_dp)::values(14)
    values=[event%birth_z,event%age_yr,event%weight,event%contact_age_yr,event%seed_wd,event%seed_donor, &
         event%donor_background_loss,event%transfer_rate,event%threshold,event%history_end_age, &
         event%initial_binary_mass,event%retention_reference_z,event%retention_mass_ceiling, &
         real(event%history_bin_num,stellar_dp)]
  end function frozen_he_numbers

  subroutine validate_frozen_he_event(event,ierr)
    type(snia_frozen_he_event_t),intent(in)::event
    integer,intent(out)::ierr
    real(stellar_dp)::age,mdot_stable,mdot_rg,m
    integer::i,j
    ierr=snia_population_contract_err_parameter
    if(.not.all(ieee_is_finite(frozen_he_numbers(event))))return
    if(event%closure_id/=snia_frozen_he_closure.or.event%comparison_approval/=snia_frozen_he_approval)return
    if(event%birth_z/=0.01_stellar_dp.or.event%retention_reference_z/=0.02_stellar_dp.or. &
         event%retention_mass_ceiling/=1.35_stellar_dp.or.event%transfer_rate/=2d-6)return
    if(event%history_bin_num/=169.or.event%weight<=0.or.event%contact_age_yr<=0.or. &
         event%seed_wd<0.6_stellar_dp.or.event%seed_wd>event%retention_mass_ceiling.or. &
         event%threshold/=1.4004633930489443_stellar_dp.or.event%seed_donor<=0.or. &
         event%donor_background_loss<0.or.event%initial_binary_mass<event%seed_wd+event%seed_donor)return
    if(event%weight*event%initial_binary_mass>1)return
    if(event%seed_donor-event%donor_background_loss<event%threshold-event%seed_wd)return
    age=event%contact_age_yr+(event%threshold-event%seed_wd)/event%transfer_rate
    if(abs(event%age_yr-age)>64*epsilon(1d0)*age.or.event%age_yr>event%history_end_age)return
    ! Both bounds increase on this mass interval. The mass argument is
    ! explicitly held at1.35 through threshold; no unlabelled extrapolation.
    do i=1,2
       m=event%seed_wd
       if(i==2)m=event%retention_mass_ceiling
       mdot_stable=1.46d-6*(-m**3+3.45d0*m**2-2.60d0*m+0.85d0)
       mdot_rg=2.17d-6*(m*m+0.82d0*m-0.38d0)
       if(event%transfer_rate<mdot_stable.or.event%transfer_rate>min(mdot_rg,2.05d-6))return
    enddo
    if(.not.is_hex_commit(event%source_commit))return
    do i=1,size(event%source_sha256)
       if(len_trim(event%source_sha256(i))/=64)return
       do j=1,64
          if(index('0123456789abcdef',event%source_sha256(i)(j:j))==0)return
       enddo
    enddo
    ierr=snia_population_contract_ok
  end subroutine validate_frozen_he_event

  subroutine validate_snia_event_source(realization,returned_mass,wd_debit,terminal_remnant,source_id,source_sha,ierr)
    type(snia_population_realization_t),intent(in)::realization
    real(stellar_dp),intent(in)::returned_mass,wd_debit,terminal_remnant
    character(len=*),intent(in)::source_id,source_sha
    integer,intent(out)::ierr
    call validate_snia_population_realization(realization,ierr)
    if(ierr/=0.or.realization%event_model==snia_event_empirical)return
    ierr=snia_population_contract_err_model
    if(.not.all(ieee_is_finite([returned_mass,wd_debit,terminal_remnant])))return
    if(returned_mass/=realization%event%threshold.or.wd_debit/=returned_mass.or.terminal_remnant/=0)return
    if(source_id/='hesma:yysd4-xap92:n100'.or. &
         source_sha/='bb9325de2e07bf86acd1ef4abf7c950fb44e74e6513b8c9c88336d910f8814c0')return
    ierr=snia_population_contract_ok
  end subroutine validate_snia_event_source

  subroutine snia_event_model_identity(realization,values)
    type(snia_population_realization_t),intent(in)::realization
    real(stellar_dp),allocatable,intent(out)::values(:)
    character(len=:),allocatable::names
    integer::i
    ! Preserve historical strict/effective identities when no table is used.
    allocate(values(0))
    if(realization%event_model==snia_event_empirical)return
    names=realization%event_model//realization%event%closure_id//realization%event%comparison_approval// &
         realization%event%source_commit
    do i=1,size(realization%event%source_sha256)
       names=names//realization%event%source_sha256(i)
    enddo
    values=[4d0,frozen_he_numbers(realization%event),(real(iachar(names(i:i)),stellar_dp),i=1,len(names))]
  end subroutine snia_event_model_identity

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

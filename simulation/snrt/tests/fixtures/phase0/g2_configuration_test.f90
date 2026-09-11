program g2_configuration_test
  use stellar_enrichment_config, only: stellar_dp, set_enrichment_defaults, &
       read_enrichment_namelist, default_imf_id, population_model_id, &
       population_binary_ssp, configured_channel_mass_min, &
       yield_source_basis_id, yield_basis_per_star_cumulative, &
       configured_binary_fraction, high_mass_model, high_mass_max_remnant_adjust_fraction, &
       configured_radioactive_model, configured_radioactive_path, high_mass_history_file, &
       stellar_fate_policy, stellar_feedback_mode, active_element
  implicit none

  integer :: unit, ios, failures, imf
  character(len=256) :: filename

  failures = 0
  filename = 'g2_configuration_test.nml'
  call set_enrichment_defaults()
  call expect(default_imf_id == 2, 'Chabrier is the compiled default', failures)
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
     write(unit, '(a)') ' channel_mass_max_msun=120.0, 8.0, 120.0, 8.0, 260.0,'
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
     write(unit, '(a)') ' channel_mass_max_msun=120.0,8.0,120.0,8.0,260.0,'
     write(unit, '(a)') '/'
     close(unit)
  end if
  open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
  if (ios == 0) then
     call read_enrichment_namelist(unit, ios)
     close(unit)
  end if
  call expect(ios /= 0, 'invalid IMF identifier is rejected', failures)
  do imf = 0, 4
     call check_imf_selection(imf)
  end do
  ! Read with omitted IMF immediately after Miller-Scalo: no inherited state.
  call check_imf_selection(-1)
  call check_high_mass('source_consistent', 0.0_stellar_dp, .true.)
  call check_high_mass('wind_only_collapse', 0.0_stellar_dp, .true.)
  call check_high_mass('mixed_remnant', 0.02_stellar_dp, .true.)
  call check_high_mass('mixed_remnant', 0.0_stellar_dp, .false.)
  call check_high_mass('source_consistent', 0.02_stellar_dp, .false.)
  call check_high_mass('unknown', 0.0_stellar_dp, .false.)
  call set_enrichment_defaults()

  ! Native frontend only: these paths are strings, not scientific fixtures.
  ! Actual LC18 history/companion identity is checked by the runtime consumer.
  block
    integer :: trial, scratch, status, previous_imf
    character(len=64) :: model, fate, previous_model, previous_fate
    character(len=1024) :: companion, history, previous_path, previous_history
    character(len=32) :: mode, previous_mode
    character(len=64) :: label
    logical :: previous_elements(size(active_element))
    real(stellar_dp) :: previous_bounds(size(configured_channel_mass_min))

    call expect(configured_radioactive_model=='none'.and.configured_radioactive_path=='', &
         'radioactive compiled defaults are off', failures)
    do trial=1,9
       previous_model=configured_radioactive_model
       previous_path=configured_radioactive_path
       previous_history=high_mass_history_file
       previous_fate=stellar_fate_policy
       previous_mode=stellar_feedback_mode
       previous_imf=default_imf_id
       previous_elements=active_element
       previous_bounds=configured_channel_mass_min
       model='lc18_al26_fe60_transparent_v1'
       companion='/input/LC18 gas/isotopes.nml'
       history='/input/LC18 gas/history.nml'
       fate='user_selected_model_v1'
       mode='channel_resolved'
       select case(trial)
       case(1); label='paired radioactive option with user-selected history'
       case(2); label='missing companion'; companion=''
       case(3); label='unknown radioactive model'; model='unknown'
       case(4); label='legacy rejects radioactive option'; mode='legacy'
       case(5)
          label='paired option requires history and user-selected fate'
          history=''; fate='review_only_unresolved'
       case(6); label='history requires user-selected fate'; fate='review_only_unresolved'
       case(7); label='user-selected fate requires history'; history=''
       case(8); label='disabled model rejects companion'; model='none'
       case(9); label='omission resets radioactive model and path'
       end select
       open(newunit=scratch,status='scratch',action='readwrite')
       write(scratch,'(a)') '&stellar_enrichment_params'
       write(scratch,'(a)') " feedback_mode='"//trim(mode)//"', population_model='single_star_ssp',"
       write(scratch,'(a)') " yield_source_basis='per_star_cumulative', use_agb=.false.,"
       write(scratch,'(a)') ' imf_mass_min_msun=.08, imf_mass_max_msun=120, binary_fraction=0,'
       write(scratch,'(a)') ' channel_mass_min_msun=13,1,13,3,140, channel_mass_max_msun=120,8,120,8,260,'
       write(scratch,'(a)') " fate_policy='"//trim(fate)//"', high_mass_history_path='"//trim(history)//"',"
       if(trial/=9)then
          write(scratch,'(a)') " radioactive_model='"//trim(model)//"',"
          write(scratch,'(a)') " radioactive_companion_path='"//trim(companion)//"',"
       endif
       if(trial>1.and.trial<9)then
          ! Valid unrelated changes must also remain uncommitted on rejection.
          write(scratch,'(a)') ' imf_id=0, use_fe=.false., channel_mass_min_msun(1)=14,'
       else
          write(scratch,'(a)') ' imf_id=4, use_fe=.true.,'
       endif
       write(scratch,'(a)') '/'
       rewind(scratch)
       call read_enrichment_namelist(scratch,status)
       close(scratch)
       if(trial==1.or.trial==9)then
          call expect(status==0,trim(label),failures)
          if(trial==1)then
             call expect(configured_radioactive_model==model.and.configured_radioactive_path==companion.and. &
                  high_mass_history_file==history.and.stellar_fate_policy==fate, &
                  'radioactive pair and source selection committed together',failures)
          else
             call expect(configured_radioactive_model=='none'.and.configured_radioactive_path=='', &
                  'omitted radioactive fields do not inherit previous active pair',failures)
          endif
       else
          call expect(status/=0,trim(label)//' rejected',failures)
          if(trial==5)call expect(status==1013,'late radioactive history/fate guard reached',failures)
          call expect(configured_radioactive_model==previous_model.and.configured_radioactive_path==previous_path.and. &
               high_mass_history_file==previous_history.and.stellar_fate_policy==previous_fate.and. &
               stellar_feedback_mode==previous_mode.and.default_imf_id==previous_imf.and. &
               all(active_element.eqv.previous_elements).and.all(configured_channel_mass_min==previous_bounds), &
               trim(label)//' preserves old state atomically',failures)
       endif
    enddo
  end block
  call set_enrichment_defaults()

  if (failures == 0) then
     write(*, '(a)') 'G2_CONFIGURATION_TEST_OK'
  else
     write(*, '(a,i0)') 'G2_CONFIGURATION_TEST_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine check_high_mass(preset, limit, accepted)
    character(len=*), intent(in) :: preset
    real(stellar_dp), intent(in) :: limit
    logical, intent(in) :: accepted
    integer :: scratch, status
    character(len=32) :: previous
    real(stellar_dp) :: previous_limit
    previous = high_mass_model
    previous_limit = high_mass_max_remnant_adjust_fraction
    open(newunit=scratch, status='scratch', action='readwrite')
    write(scratch, '(a)') '&stellar_enrichment_params'
    write(scratch, '(a)') " feedback_mode='channel_resolved', population_model='single_star_ssp',"
    write(scratch, '(a)') " yield_source_basis='per_star_cumulative',"
    write(scratch, '(a)') ' imf_mass_min_msun=.08, imf_mass_max_msun=120, binary_fraction=0,'
    write(scratch, '(a)') ' channel_mass_min_msun=.8,1,8,3,140, channel_mass_max_msun=120,8,120,8,260,'
    write(scratch, '(a)') " high_mass_preset='"//preset//"',"
    write(scratch, '(a,es24.16)') ' high_mass_remnant_adjust_max_fraction=',limit
    write(scratch, '(a)') '/'
    rewind(scratch)
    call read_enrichment_namelist(scratch, status)
    close(scratch)
    call expect((status == 0) .eqv. accepted, 'preset combination admission: '//preset, failures)
    if (accepted) then
       call expect(high_mass_model == preset .and. high_mass_max_remnant_adjust_fraction == limit, &
            'selected values propagate: '//preset, failures)
    else
       call expect(high_mass_model == previous .and. high_mass_max_remnant_adjust_fraction == previous_limit, &
            'invalid read preserves previous state', failures)
    end if
  end subroutine check_high_mass

  subroutine check_imf_selection(requested)
    integer, intent(in) :: requested
    integer :: scratch, status, expected

    expected = requested
    if (requested < 0) expected = 2
    open(newunit=scratch, status='scratch', action='readwrite')
    write(scratch, '(a)') '&stellar_enrichment_params'
    write(scratch, '(a)') " feedback_mode='channel_resolved',"
    if (requested >= 0) write(scratch, '(a,i0,a)') ' imf_id=', requested, ','
    write(scratch, '(a)') " population_model='single_star_ssp',"
    write(scratch, '(a)') " yield_source_basis='per_star_cumulative',"
    write(scratch, '(a)') ' imf_mass_min_msun=0.08, imf_mass_max_msun=120.0,'
    write(scratch, '(a)') ' binary_fraction=0.0,'
    write(scratch, '(a)') ' channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0,'
    write(scratch, '(a)') ' channel_mass_max_msun=120.0,8.0,120.0,8.0,260.0,'
    write(scratch, '(a)') '/'
    rewind(scratch)
    call read_enrichment_namelist(scratch, status)
    close(scratch)
    call expect(status == 0 .and. default_imf_id == expected, &
         'IMF explicit selection or deterministic Chabrier omission', failures)
  end subroutine check_imf_selection

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

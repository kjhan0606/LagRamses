program test_stellar_feedback_policy
  use stellar_enrichment_config, only: stellar_dp, channel_agb, channel_snii, &
       channel_snia, stellar_feedback_mode, set_enrichment_defaults, &
       read_enrichment_namelist, use_channel_resolved_feedback, &
       default_imf_id, population_model_id, population_binary_ssp, &
       population_single_star_ssp, &
       yield_source_basis_id, yield_basis_per_star_cumulative, &
       configured_imf_mass_min, configured_imf_mass_max, &
       configured_binary_fraction, &
       configured_channel_mass_min, configured_channel_mass_max, &
       production_source_model_supported, enable_wind
  use stellar_enrichment_contract, only: stellar_source_t, clear_source, &
       delayed_cooling_source_mass
  implicit none

  type(stellar_source_t) :: source
  integer :: unit, status

  call set_enrichment_defaults()
  if (.not. use_channel_resolved_feedback()) error stop 1

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&unrelated_group value=1 /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1005) error stop 2

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') &
       '&stellar_enrichment_params feedback_mode="legacy", use_wind=.false. /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. trim(stellar_feedback_mode) /= 'legacy' .or. &
       enable_wind) error stop 3

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'CHANNEL_RESOLVED',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. .not. use_channel_resolved_feedback()) error stop 4

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params feedback_mode="invalid" /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1001) error stop 5

  call set_enrichment_defaults()
  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  write(unit,'(A)') ' feedback_mode="channel_resolved",'
  write(unit,'(A)') ' imf_id=2,'
  write(unit,'(A)') ' population_model="binary_ssp",'
  write(unit,'(A)') ' yield_source_basis="per_star_cumulative",'
  write(unit,'(A)') ' imf_mass_min_msun=0.08, imf_mass_max_msun=120.0,'
  write(unit,'(A)') ' binary_fraction=0.5,'
  write(unit,'(A)') ' channel_mass_min_msun=1.2,1.0,8.0,3.0,140.0,'
  write(unit,'(A)') ' channel_mass_max_msun=120.0,8.0,40.0,8.0,260.0,'
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. default_imf_id /= 2 .or. &
       population_model_id /= population_binary_ssp) error stop 6
  if (abs(configured_channel_mass_min(1)-1.2_stellar_dp) > 1.0d-12 .or. &
       abs(configured_channel_mass_max(1)-120.0_stellar_dp) > 1.0d-12) error stop 7
  if (yield_source_basis_id /= yield_basis_per_star_cumulative .or. &
       abs(configured_imf_mass_min-0.08_stellar_dp) > 1.0d-12 .or. &
       abs(configured_imf_mass_max-120.0_stellar_dp) > 1.0d-12 .or. &
       abs(configured_binary_fraction-0.5_stellar_dp) > 1.0d-12) error stop 8
  if (production_source_model_supported()) error stop 9

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',99, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1002) error stop 10
  if (default_imf_id /= 2 .or. population_model_id /= population_binary_ssp .or. &
       abs(configured_channel_mass_min(1)-1.2_stellar_dp) > 1.0d-12) error stop 11

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'unknown','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1003) error stop 12
  if (default_imf_id /= 2 .or. population_model_id /= population_binary_ssp .or. &
       abs(configured_channel_mass_min(1)-1.2_stellar_dp) > 1.0d-12) error stop 13

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp,.false.)
  write(unit,'(A)') ' channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0,'
  write(unit,'(A)') ' channel_mass_max_msun=0.7,8.0,40.0,8.0,260.0,'
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1004) error stop 14
  if (default_imf_id /= 2 .or. population_model_id /= population_binary_ssp .or. &
       abs(configured_channel_mass_min(1)-1.2_stellar_dp) > 1.0d-12) error stop 15

  call set_enrichment_defaults()
  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') ' use_snia=.true. /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. production_source_model_supported()) error stop 16

  call set_enrichment_defaults()
  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') ' use_pisn=.true. /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. production_source_model_supported()) error stop 17

  call set_enrichment_defaults()
  if (population_model_id /= population_single_star_ssp) error stop 18

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','ssp_cumulative_per_initial_mass',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. production_source_model_supported()) error stop 19

  call set_enrichment_defaults()
  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1006) error stop 20

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.5_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1008) error stop 21

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'binary_ssp','per_star_cumulative',0.0_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1008) error stop 22

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp,.false.)
  write(unit,'(A)') ' channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0,'
  write(unit,'(A)') ' channel_mass_max_msun=121.0,8.0,40.0,8.0,260.0,'
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1009) error stop 23

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params'
  call write_required_population_contract(unit,'channel_resolved',1, &
       'single_star_ssp','per_star_cumulative',0.0_stellar_dp, &
       imf_lower=0.01_stellar_dp)
  write(unit,'(A)') '/'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1007) error stop 24

  call clear_source(source)
  source%channel_returned_mass(channel_agb) = 2.0_stellar_dp
  source%channel_returned_mass(channel_snii) = 3.0_stellar_dp
  source%channel_returned_mass(channel_snia) = 5.0_stellar_dp
  source%returned_mass = sum(source%channel_returned_mass)
  if (delayed_cooling_source_mass(source) /= 3.0_stellar_dp) error stop 25

  write(*,'(A)') 'stellar feedback policy: PASS'

contains

  subroutine write_required_population_contract(unit,mode,imf,population,basis, &
       binary,write_windows,imf_lower)
    integer, intent(in) :: unit, imf
    character(len=*), intent(in) :: mode, population, basis
    real(stellar_dp), intent(in) :: binary
    logical, intent(in), optional :: write_windows
    real(stellar_dp), intent(in), optional :: imf_lower
    logical :: include_windows
    real(stellar_dp) :: lower

    include_windows=.true.
    if(present(write_windows))include_windows=write_windows
    lower=0.08_stellar_dp
    if(present(imf_lower))lower=imf_lower
    write(unit,'(A)') ' feedback_mode="'//trim(mode)//'",'
    write(unit,'(A,I0,A)') ' imf_id=',imf,','
    write(unit,'(A)') ' population_model="'//trim(population)//'",'
    write(unit,'(A)') ' yield_source_basis="'//trim(basis)//'",'
    write(unit,'(A,ES16.8,A)') ' imf_mass_min_msun=',lower, &
         ', imf_mass_max_msun=120.0,'
    write(unit,'(A,ES16.8,A)') ' binary_fraction=',binary,','
    if(include_windows)then
       write(unit,'(A)') ' channel_mass_min_msun=0.8,1.0,8.0,3.0,140.0,'
       write(unit,'(A)') ' channel_mass_max_msun=120.0,8.0,40.0,8.0,260.0,'
    endif
  end subroutine write_required_population_contract
end program test_stellar_feedback_policy

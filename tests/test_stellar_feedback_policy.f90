program test_stellar_feedback_policy
  use stellar_enrichment_config, only: stellar_dp, channel_agb, channel_snii, &
       channel_snia, stellar_feedback_mode, set_enrichment_defaults, &
       read_enrichment_namelist, use_channel_resolved_feedback
  use stellar_enrichment_contract, only: stellar_source_t, clear_source, &
       delayed_cooling_source_mass
  implicit none

  type(stellar_source_t) :: source
  integer :: unit, status

  call set_enrichment_defaults()
  if (.not. use_channel_resolved_feedback()) error stop 1

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params feedback_mode="legacy" /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. trim(stellar_feedback_mode) /= 'legacy') error stop 2

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params feedback_mode="CHANNEL_RESOLVED" /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 0 .or. .not. use_channel_resolved_feedback()) error stop 3

  open(newunit=unit,status='scratch',action='readwrite')
  write(unit,'(A)') '&stellar_enrichment_params feedback_mode="invalid" /'
  rewind(unit)
  call read_enrichment_namelist(unit,status)
  close(unit)
  if (status /= 1001) error stop 4

  call clear_source(source)
  source%channel_returned_mass(channel_agb) = 2.0_stellar_dp
  source%channel_returned_mass(channel_snii) = 3.0_stellar_dp
  source%channel_returned_mass(channel_snia) = 5.0_stellar_dp
  source%returned_mass = sum(source%channel_returned_mass)
  if (delayed_cooling_source_mass(source) /= 3.0_stellar_dp) error stop 5

  write(*,'(A)') 'stellar feedback policy: PASS'
end program test_stellar_feedback_policy

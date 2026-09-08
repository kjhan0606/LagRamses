program lc18_native_wind_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_enrichment_driver
  use stellar_enrichment_contract
  use stellar_population_ledger
  implicit none
  type(stellar_yield_table_t)::table
  type(stellar_cumulative_t)::wind,terminal
  type(stellar_population_t)::pop
  type(stellar_source_t)::whole,early,late
  type(stellar_population_ledger_t)::ledger
  character(len=1024)::path,history,mode
  real(stellar_dp)::final_age
  logical::lowmass
  integer::ierr,i,row
  call get_command_argument(1,path);call get_command_argument(2,history)
  call get_command_argument(3,mode)
  lowmass=trim(mode)=='sukhbold_low'
  if(len_trim(mode)>0.and..not.lowmass)stop 20
  final_age=.01d0
  call set_enrichment_defaults()
  high_mass_model='wind_only_collapse';enable_agb=.false.
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  if(lowmass)then
     configured_channel_mass_min(1)=9;configured_channel_mass_min(3)=9
     configured_channel_mass_max(1)=13;configured_channel_mass_max(3)=13
     final_age=.1d0
  endif
  call load_yield_table(trim(path),table,ierr)
  if(ierr/=0)stop 1
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0.or..not.table%net_yield_diagnostic_unavailable)stop 2
  do i=1,size(table%hm_mass)
     row=table%hm_wind_row(i)
     call evaluate_channel_cumulative(table,1,table%hm_mass(i),table%hm_z(i),table%hm_age(i),wind,ierr)
     if(ierr/=0.or.wind%returned_mass/=table%returned_mass(row))stop 3
     if(any(wind%ejected_mass/=table%ejected_mass(row,:)))stop 4
     call evaluate_channel_cumulative(table,3,table%hm_mass(i),table%hm_z(i),table%hm_age(i),terminal,ierr)
     if(ierr/=0)stop 5
     if(lowmass)then
        if(terminal%returned_mass<=0.or.terminal%energy/=table%energy(table%hm_terminal_row(i)))stop 15
     else
        if(terminal%returned_mass/=0.or.terminal%energy/=0)stop 5
     endif
     if(abs(wind%returned_mass+terminal%returned_mass+terminal%remnant_mass-table%hm_mass(i))>1d-10)stop 6
     call evaluate_channel_cumulative(table,1,table%hm_mass(i),table%hm_z(i),.5d0*table%hm_age(i),wind,ierr)
     if(ierr/=0.or.abs(wind%returned_mass-.5d0*table%returned_mass(row))>1d-10)stop 7
     if(abs(wind%energy-.5d0*table%energy(row))>1d-12*table%energy(row))stop 8
  enddo
  call evaluate_channel_cumulative(table,1,60d0,0d0,.01d0,wind,ierr)
  if(ierr==0)stop 9
  pop%initial_mass=10000;pop%current_mass=10000;pop%birth_metallicity=.01d0
  if(lowmass)then
     pop%birth_metallicity=.02d0
     if(size(table%hm_mass)/=17)stop 16
     call evaluate_channel_cumulative(table,3,9d0,.019d0,.1d0,terminal,ierr)
     if(ierr==0)stop 17
     call evaluate_channel_cumulative(table,3,8.999d0,.02d0,.1d0,terminal,ierr)
     if(ierr==0)stop 18
     call evaluate_channel_cumulative(table,3,13.001d0,.02d0,.1d0,terminal,ierr)
     if(ierr==0)stop 19
  endif
  pop%imf_id=2;pop%imf_mass_min=.08d0;pop%imf_mass_max=120
  pop%yield_basis_id=yield_basis_per_star_cumulative
  call compute_stellar_source_increment(table,pop,0d0,final_age,configured_channel_mass_min, &
       configured_channel_mass_max,64,whole,ierr,ledger)
  if(ierr/=0.or.whole%returned_mass<=0)stop 10
  call compute_stellar_source_increment(table,pop,0d0,.002d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,early,ierr)
  if(ierr/=0)stop 11
  call compute_stellar_source_increment(table,pop,.002d0,final_age,configured_channel_mass_min, &
       configured_channel_mass_max,64,late,ierr)
  if(ierr/=0.or.abs(whole%returned_mass-early%returned_mass-late%returned_mass)>1d-9)stop 12
  if(abs(whole%energy-early%energy-late%energy)>1d-12*whole%energy)stop 13
  if(abs(ledger%living_mass+ledger%remnant_mass+ledger%returned_mass-10000d0)>1d-9)stop 14
  print *, 'LC18_NATIVE_WIND_TEST_OK',size(table%hm_mass),whole%returned_mass,whole%energy
  if(lowmass)print *, 'SUKHBOLD_LOW_MASS_NATIVE_TEST_OK',whole%channel_returned_mass(3),whole%channel_energy(3)
end program lc18_native_wind_test

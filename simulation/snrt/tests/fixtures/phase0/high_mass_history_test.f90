program high_mass_history_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_enrichment_driver
  use stellar_enrichment_contract
  use stellar_population_ledger
  implicit none
  type(stellar_yield_table_t)::table,raw
  type(stellar_cumulative_t)::state
  type(stellar_population_t)::population
  type(stellar_source_t)::whole,early,late,repeat
  type(stellar_population_ledger_t)::ledger
  character(len=1024)::path,history
  integer::ierr,k
  call get_command_argument(1,path)
  call get_command_argument(2,history)
  call set_enrichment_defaults()
  enable_agb=.false.
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  call load_yield_table(trim(path),raw,ierr)
  if(ierr/=0)stop 1
  do k=1,3
     table=raw
     select case(k)
     case(1)
        high_mass_model='source_consistent';high_mass_max_remnant_adjust_fraction=0
     case(2)
        high_mass_model='wind_only_collapse';high_mass_max_remnant_adjust_fraction=0
     case(3)
        high_mass_model='mixed_remnant';high_mass_max_remnant_adjust_fraction=.02d0
        where(table%channel==3.and.table%age_gyr>=1d-6) &
             table%remnant_mass=table%remnant_mass-.01d0*table%initial_mass
     end select
     call prepare_high_mass_history(table,trim(history),ierr)
     if(ierr/=0.or..not.table%high_mass_ready)stop 2
     call evaluate_channel_cumulative(table,3,60d0,.01d0,.999d-6,state,ierr)
     if(ierr/=0.or.state%returned_mass/=0.or.state%remnant_mass/=0)stop 3
     call evaluate_channel_cumulative(table,3,60d0,.01d0,table%hm_age(2),state,ierr)
     if(ierr/=0)stop 4
     if(k==2)then
        if(state%returned_mass/=0.or.state%remnant_mass/=54d0)stop 5
     else
        if(state%returned_mass/=12d0.or.abs(state%remnant_mass-42d0)>1d-12)stop 6
     endif
     call evaluate_channel_cumulative(table,1,60d0,.01d0,.5d-6,state,ierr)
     if(ierr/=0.or.abs(state%returned_mass-3d0)>1d-12)stop 7
     call evaluate_channel_cumulative(table,1,60d0,nearest(.01d0,1d0),.5d-6,state,ierr)
     if(ierr/=0.or.abs(state%returned_mass-3d0)>1d-12)stop 19
     call evaluate_channel_cumulative(table,1,60d0,.011d0,.5d-6,state,ierr)
     if(ierr==0)stop 8 ! No outcome interpolation across metallicities.
     call evaluate_channel_cumulative(table,1,121d0,.01d0,2d-6,state,ierr)
     if(ierr==0)stop 9
     population%initial_mass=10000;population%current_mass=10000
     population%birth_metallicity=.01d0
     population%imf_id=default_imf_id;population%imf_mass_min=.08d0;population%imf_mass_max=120
     population%yield_basis_id=yield_basis_per_star_cumulative
     call compute_stellar_source_increment(table,population,0d0,2d-6,configured_channel_mass_min, &
          configured_channel_mass_max,64,whole,ierr,ledger)
     if(ierr/=0.or.whole%returned_mass<=0)stop 10
     call compute_stellar_source_increment(table,population,0d0,.5d-6,configured_channel_mass_min, &
          configured_channel_mass_max,64,early,ierr)
     if(ierr/=0)stop 11
     call compute_stellar_source_increment(table,population,.5d-6,2d-6,configured_channel_mass_min, &
          configured_channel_mass_max,64,late,ierr)
     if(ierr/=0)stop 12
     if(abs(early%returned_mass+late%returned_mass-whole%returned_mass)>1d-10)stop 13
     if(maxval(abs(early%ejected_mass+late%ejected_mass-whole%ejected_mass))>1d-10)stop 14
     if(abs(early%energy+late%energy-whole%energy)>1d-12*whole%energy)stop 15
     call compute_stellar_source_increment(table,population,2d-6,3d-6,configured_channel_mass_min, &
          configured_channel_mass_max,64,repeat,ierr)
     if(ierr/=0.or.repeat%returned_mass/=0.or.repeat%energy/=0)stop 16
     if(abs(ledger%living_mass+ledger%returned_mass+ledger%remnant_mass-10000d0)>1d-9)stop 17
  enddo
  table=raw
  table%returned_mass(3)=table%returned_mass(3)+1000
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr==0.or.table%high_mass_ready)stop 18
  print *, 'HIGH_MASS_HISTORY_INCREMENT_TEST_OK'
end program high_mass_history_test

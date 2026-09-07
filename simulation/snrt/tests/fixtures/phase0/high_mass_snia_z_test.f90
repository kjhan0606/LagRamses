program high_mass_snia_z_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_enrichment_driver
  use stellar_enrichment_contract
  use stellar_population_ledger
  use stellar_snia_population_contract
  use stellar_snia_physical_contract
  implicit none
  type(stellar_yield_table_t)::table
  type(stellar_cumulative_t)::a,b,c
  type(stellar_population_t)::pop
  type(stellar_source_t)::whole,early,late
  type(stellar_population_ledger_t)::ledger
  type(snia_population_realization_t)::dtd
  type(snia_physical_contract_t)::event
  type(snia_event_budget_t)::budget,first,last
  character(len=1024)::yields,history,snia
  real(stellar_dp)::n,n1,n2,wd,bh
  integer::ierr,u,j
  call get_command_argument(1,yields);call get_command_argument(2,history);call get_command_argument(3,snia)
  call set_enrichment_defaults()
  default_imf_id=1;population_model_id=1;configured_binary_fraction=.5
  stellar_fate_policy='user_selected_model_v1';high_mass_history_file=history
  enable_snia=.true.;enable_agb=.true.
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  if(.not.user_source_model_requested())stop 1
  call load_yield_table(trim(yields),table,ierr)
  if(ierr/=0)stop 2
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0.or..not.table%high_mass_linear_z)stop 3
  call set_yield_mass_assignment_mode(table,yield_mass_assignment_piecewise_constant,ierr)
  do j=1,3
     if(j==2)cycle
     call evaluate_channel_cumulative(table,j,60d0,0d0,1.5d-6,a,ierr)
     if(ierr/=0)stop 4
     call evaluate_channel_cumulative(table,j,60d0,.02d0,1.5d-6,b,ierr)
     if(ierr/=0)stop 5
     call evaluate_channel_cumulative(table,j,60d0,.01d0,1.5d-6,c,ierr)
     if(ierr/=0)stop 6
     if(abs(c%returned_mass-.5d0*(a%returned_mass+b%returned_mass))>1d-12)stop 7
     if(abs(c%remnant_mass-.5d0*(a%remnant_mass+b%remnant_mass))>1d-12)stop 8
  enddo
  call evaluate_channel_cumulative(table,1,60d0,.021d0,1.5d-6,c,ierr)
  if(ierr==0.or.c%returned_mass/=0)stop 9
  open(newunit=u,file=trim(snia),status='old',action='read')
  call read_snia_population_realization_namelist(u,dtd,ierr)
  if(ierr/=0)stop 10
  call read_snia_physical_contract_namelist(u,event,ierr)
  close(u)
  if(ierr/=0)then
     print *, 'physical event load failed',ierr,event%energy_per_event,event%returned_mass_per_event
     stop 11
  endif
  call validate_snia_population_binding(dtd,default_imf_id,population_model_id,configured_binary_fraction,ierr)
  if(ierr/=0)stop 12
  pop%initial_mass=10000;pop%current_mass=10000;pop%birth_metallicity=.01d0
  pop%imf_id=1;pop%population_id=1;pop%imf_mass_min=.08d0;pop%imf_mass_max=120
  pop%yield_basis_id=yield_basis_per_star_cumulative
  call compute_stellar_source_increment(table,pop,0d0,.1d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,whole,ierr,ledger)
  if(ierr/=0)stop 13
  call compute_stellar_source_increment(table,pop,0d0,1.5d-6,configured_channel_mass_min, &
       configured_channel_mass_max,64,early,ierr)
  if(ierr/=0)stop 14
  call compute_stellar_source_increment(table,pop,1.5d-6,.1d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,late,ierr)
  if(ierr/=0.or.abs(whole%returned_mass-early%returned_mass-late%returned_mass)>1d-9)then
     print *, 'split increment failure',ierr,whole%returned_mass,early%returned_mass,late%returned_mass
     stop 15
  endif
  call evaluate_snia_interval_events(dtd,pop%initial_mass,0d0,.1d0,1d0,n,ierr)
  if(ierr/=0.or.n<=0)stop 16
  call evaluate_snia_interval_events(dtd,pop%initial_mass,0d0,.06d0,1d0,n1,ierr)
  if(ierr/=0)stop 17
  call evaluate_snia_interval_events(dtd,pop%initial_mass,.06d0,.1d0,1d0,n2,ierr)
  if(ierr/=0.or.abs(n-n1-n2)>1d-12)stop 18
  wd=ledger%channel_remnant_mass(channel_agb);bh=ledger%channel_remnant_mass(channel_snii)
  call build_snia_event_budget(event,n,wd,budget,ierr)
  if(ierr/=0.or.budget%returned_mass<=0)stop 19
  call build_snia_event_budget(event,n1,wd,first,ierr)
  if(ierr/=0)stop 20
  call build_snia_event_budget(event,n2,wd-first%wd_reservoir_debit,last,ierr)
  if(ierr/=0.or.abs(budget%returned_mass-first%returned_mass-last%returned_mass)>1d-12)stop 21
  call set_white_dwarf_reservoir(ledger,wd,1d-10,ierr)
  if(ierr/=0)stop 22
  call apply_snia_event_budget(ledger,budget,1d-10,ierr)
  if(ierr/=0.or.ledger%channel_remnant_mass(channel_snii)/=bh)stop 23
  if(abs(ledger%living_mass+ledger%returned_mass+ledger%remnant_mass-pop%initial_mass)>1d-9)stop 24
  call build_snia_event_budget(event,n,0d0,budget,ierr)
  if(ierr==0)stop 25 ! A high-mass remnant is not a WD reservoir.
  print *, 'HIGH_MASS_SNIA_Z_PHYSICAL_EVENT_TEST_OK',n,event%energy_per_event
end program high_mass_snia_z_test

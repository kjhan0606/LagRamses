program kl16_lc18_native_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_enrichment_contract
  use stellar_enrichment_driver
  use stellar_population_ledger
  use stellar_snia_population_contract
  use stellar_snia_physical_contract
  use stellar_snia_runtime_accounting
  implicit none
  type(stellar_yield_table_t)::table
  type(stellar_cumulative_t)::a
  type(stellar_source_t)::whole,early,late
  type(stellar_population_t)::pop
  type(stellar_population_ledger_t)::ledger
  type(snia_population_realization_t)::dtd
  type(snia_physical_contract_t)::event
  type(snia_event_budget_t)::budget
  character(len=1024)::yields,history,snia
  real(stellar_dp)::energy,first_age
  real(stellar_dp)::times(9)=[0d0,.04d0,.05d0,.06889d0,.1d0,.2d0,1d0,5d0,13.7d0]
  ! Common KL16/LC18 Z support, not the broader AGB-only domain.
  real(stellar_dp)::zs(4)=[.007d0,.01d0,.012d0,.01345d0]
  real(stellar_dp)::prior,expected_prior,events,remaining,generic_left,all_return,ia_return
  integer::i,r,ierr,u,iz,it
  call get_command_argument(1,yields);call get_command_argument(2,history);call get_command_argument(3,snia)
  call set_enrichment_defaults()
  default_imf_id=1;population_model_id=1;configured_binary_fraction=.5d0
  enable_agb=.true.;enable_snia=.true.;high_mass_model='wind_only_collapse'
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  configured_channel_mass_max(2)=6
  call load_yield_table(trim(yields),table,ierr)
  if(ierr/=0)stop 1
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0.or..not.allocated(table%agb_terminal_row))stop 2
  if(size(table%agb_terminal_row)/=58.or..not.table%net_yield_diagnostic_unavailable)stop 3
  call set_yield_mass_assignment_mode(table,yield_mass_assignment_piecewise_constant,ierr)
  call audit_yield_table(table,1d-10,ierr,.true.,channel_owns_terminal_remnant, &
       [.true.,.true.,.true.,.false.,.false.])
  if(ierr/=0)stop 4
  first_age=huge(1d0)
  do i=1,size(table%agb_terminal_row)
     r=table%agb_terminal_row(i);first_age=min(first_age,table%age_gyr(r))
     call evaluate_channel_cumulative(table,2,table%initial_mass(r),table%birth_metallicity(r), &
          nearest(table%age_gyr(r),-1d0),a,ierr)
     if(ierr/=0.or.a%returned_mass/=0.or.a%remnant_mass/=0)stop 5
     call evaluate_channel_cumulative(table,2,table%initial_mass(r),table%birth_metallicity(r), &
          table%age_gyr(r),a,ierr)
     if(ierr/=0.or.a%returned_mass/=table%returned_mass(r).or.any(a%ejected_mass/=table%ejected_mass(r,:)))stop 6
     if(abs(a%returned_mass+a%remnant_mass-table%initial_mass(r))>1d-12)stop 7
     energy=.5d0*a%returned_mass*1.98847d33*(15d5)**2
     if(abs(a%energy-energy)>1d-14*energy.or.any(a%momentum/=0))stop 8
  enddo
  pop%initial_mass=10000;pop%current_mass=10000;pop%birth_metallicity=.01d0
  pop%imf_id=1;pop%population_id=1;pop%imf_mass_min=.08d0;pop%imf_mass_max=120
  pop%yield_basis_id=yield_basis_per_star_cumulative
  call compute_stellar_source_increment(table,pop,0d0,.2d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,whole,ierr,ledger)
  if(ierr/=0.or.whole%channel_returned_mass(2)<=0.or.ledger%channel_remnant_mass(2)<=0)stop 9
  call compute_stellar_source_increment(table,pop,0d0,.1d0,configured_channel_mass_min, &
       configured_channel_mass_max,8,early,ierr)
  if(ierr/=0)stop 10
  call compute_stellar_source_increment(table,pop,.1d0,.2d0,configured_channel_mass_min, &
       configured_channel_mass_max,127,late,ierr)
  if(ierr/=0.or.abs(whole%returned_mass-early%returned_mass-late%returned_mass)>1d-9)stop 11
  if(abs(ledger%living_mass+ledger%returned_mass+ledger%remnant_mass-10000d0)>1d-9)stop 12
  open(newunit=u,file=trim(snia),status='old',action='read')
  call read_snia_population_realization_namelist(u,dtd,ierr)
  if(ierr/=0)stop 13
  call read_snia_physical_contract_namelist(u,event,ierr)
  close(u)
  if(ierr/=0)stop 14
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.04d0,1d0,6d0,64,ierr)
  if(ierr/=0)stop 15
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.05d0,1d0,6d0,64,ierr)
  if(ierr/=enrichment_driver_err_ledger)stop 16
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.2d0,1d0,6d0,64,ierr)
  if(ierr/=enrichment_driver_err_ledger)stop 17 ! Large dt cannot conceal the earlier shortfall.
  print *, 'KL16_LC18_REAL_SOURCE_TEST_OK',size(table%agb_terminal_row),whole%returned_mass, &
       whole%channel_returned_mass(2),ledger%channel_remnant_mass(2)
  print *, 'UNCHANGED_40_MYR_DTD_INCOMPATIBLE_WITH_SELECTED_WD_SUPPLY',first_age

  ! Same real sources, now explicitly effective SSP accounting. This checks
  ! aggregate closure over the full declared DTD, NOT a physical WD supplier.
  dtd%mass_accounting=snia_accounting_effective_ssp
  dtd%accounting_approval_id=snia_effective_ssp_approval
  do iz=1,size(zs)
     pop%birth_metallicity=zs(iz);pop%current_mass=pop%initial_mass
     all_return=0d0;ia_return=0d0
     do it=2,size(times)
        call compute_stellar_source_increment(table,pop,times(it-1),times(it),configured_channel_mass_min, &
             configured_channel_mass_max,64,whole,ierr,ledger)
        if(ierr/=0)stop 18
        generic_left=ledger%living_mass+ledger%remnant_mass
        call reconstruct_prior_snia_return(pop%current_mass,whole%returned_mass,generic_left, &
             pop%initial_mass,1d-10,prior,ierr)
        if(ierr/=0)stop 19
        call evaluate_snia_interval_events(dtd,pop%initial_mass,0d0,times(it-1),1d0,expected_prior,ierr)
        if(ierr/=0)stop 20
        call evaluate_snia_interval_events(dtd,pop%initial_mass,times(it-1),times(it),1d0,events,ierr)
        if(ierr/=0)stop 21
        call build_snia_event_budget(event,events,pop%current_mass-whole%returned_mass,budget,ierr)
        if(ierr/=0)stop 22
        call close_effective_snia_return(pop%current_mass,whole%returned_mass,generic_left,prior, &
             expected_prior*event%returned_mass_per_event,budget%returned_mass,pop%initial_mass, &
             1d-10,remaining,ierr)
        if(ierr/=0)stop 23
        all_return=all_return+whole%returned_mass+budget%returned_mass
        ia_return=ia_return+budget%returned_mass
        pop%current_mass=remaining
        if(abs(remaining+all_return-pop%initial_mass)>1d-8)stop 24
        if(abs(budget%energy-events*event%energy_per_event)>1d-12*max(1d0,budget%energy))stop 25
        if(any(abs(budget%ejected_mass-events*event%ejected_mass_per_event)>1d-10))stop 26
     enddo
     if(abs(ia_return-pop%initial_mass*.0013d0*event%returned_mass_per_event)>1d-9)stop 27
     print *, 'EFFECTIVE_SSP_FULL_DTD_MASS_CLOSURE',zs(iz),remaining,all_return,ia_return
  enddo
  print *, 'KL16_LC18_EFFECTIVE_SSP_TEST_OK'
end program kl16_lc18_native_test

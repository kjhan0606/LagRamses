program agb_terminal_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_enrichment_driver
  use stellar_enrichment_contract
  use stellar_ssp_sources
  use stellar_population_ledger
  use stellar_snia_population_contract
  use stellar_snia_physical_contract
  implicit none
  type(stellar_yield_table_t)::raw,table,bad,legacy
  type(stellar_cumulative_t)::a,b
  type(stellar_population_t)::pop
  type(stellar_source_t)::whole,early,late
  type(stellar_population_ledger_t)::ledger
  type(snia_population_realization_t)::dtd
  type(snia_physical_contract_t)::event
  type(snia_event_budget_t)::budget
  character(len=1024)::yields,history,snia,old_history,outfile
  real(stellar_dp)::events
  integer::ierr,i,u,r
  call get_command_argument(1,yields);call get_command_argument(2,history)
  call get_command_argument(3,snia);call get_command_argument(4,old_history)
  call get_command_argument(5,outfile)
  call set_enrichment_defaults()
  default_imf_id=1;population_model_id=1;configured_binary_fraction=.5
  enable_snia=.true.;enable_agb=.true.
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  call load_yield_table(trim(yields),raw,ierr)
  if(ierr/=0)stop 1
  ! Reuse the existing SYNTHETIC material fixture, change only its AGB
  ! timing to a common [0,20 Myr,100 Myr,10 Gyr] grid. Do not label as KL16.
  do i=1,raw%n_rows
     if(raw%channel(i)/=channel_agb)cycle
     if(raw%age_gyr(i)==1000d0*1d-9)raw%age_gyr(i)=.02d0
     if(raw%age_gyr(i)==2000d0*1d-9)raw%age_gyr(i)=.1d0
     if(raw%age_gyr(i)==.02d0.and.(raw%initial_mass(i)==1.or.raw%birth_metallicity(i)==.02d0))then
        raw%returned_mass(i)=0;raw%remnant_mass(i)=0;raw%energy(i)=0
        raw%momentum(i,:)=0;raw%ejected_mass(i,:)=0;raw%net_yield(i,:)=0
     endif
  enddo
  if(len_trim(outfile)>0)then
     open(newunit=u,file=trim(outfile),status='new',action='write')
     write(u,'(A)')'# SYNTHETIC AGB terminal-release timing control, not physical yield data.'
     do i=1,raw%n_rows
        write(u,'(I2,31ES25.17)')raw%channel(i),raw%initial_mass(i),raw%birth_metallicity(i), &
             raw%age_gyr(i)*1d9,raw%returned_mass(i),raw%remnant_mass(i),raw%energy(i), &
             raw%momentum(i,:),raw%ejected_mass(i,:),raw%net_yield(i,:)
     enddo
     close(u)
  endif
  table=raw
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0.or..not.allocated(table%agb_terminal_row))stop 2
  call set_yield_mass_assignment_mode(table,yield_mass_assignment_piecewise_constant,ierr)
  call evaluate_channel_cumulative(table,2,8d0,0d0,nearest(.02d0,-1d0),a,ierr)
  if(ierr/=0.or.a%returned_mass/=0.or.a%remnant_mass/=0.or.a%energy/=0)stop 3
  call evaluate_channel_cumulative(table,2,8d0,0d0,.02d0,a,ierr)
  if(ierr/=0.or.a%returned_mass/=6.4d0.or.a%remnant_mass/=1.6d0)stop 4
  call evaluate_channel_cumulative(table,2,8d0,.01d0,.07d0,a,ierr)
  if(ierr/=0.or.abs(a%returned_mass-3.2d0)>1d-12.or.abs(a%remnant_mass-.8d0)>1d-12)stop 5
  call evaluate_channel_cumulative(table,2,4.5d0,0d0,.07d0,a,ierr)
  if(ierr/=0.or.a%returned_mass/=0)stop 6 ! Lower-node tie, lower node has not died.
  call evaluate_channel_cumulative(table,2,4.6d0,0d0,.07d0,a,ierr)
  if(ierr/=0.or.abs(a%returned_mass-3.68d0)>1d-12)stop 7
  call evaluate_channel_cumulative(table,2,8d0,.021d0,.07d0,a,ierr)
  if(ierr==0.or.a%returned_mass/=0)stop 8
  call evaluate_channel_cumulative(table,2,8.01d0,.01d0,.07d0,a,ierr)
  if(ierr==0.or.a%returned_mass/=0)stop 9
  pop%initial_mass=10000;pop%current_mass=10000;pop%birth_metallicity=.01d0
  pop%imf_id=1;pop%population_id=1;pop%imf_mass_min=.08d0;pop%imf_mass_max=120
  pop%yield_basis_id=yield_basis_per_star_cumulative
  call integrate_ssp_channel(table,pop,2,.07d0,1d0,8d0,8,a,ierr)
  if(ierr/=0)stop 10
  call integrate_ssp_channel(table,pop,2,.07d0,1d0,8d0,127,b,ierr)
  if(ierr/=0.or.abs(a%remnant_mass-b%remnant_mass)>1d-10.or.a%remnant_mass<=0)stop 11
  call compute_stellar_source_increment(table,pop,0d0,.15d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,whole,ierr,ledger)
  if(ierr/=0)stop 12
  call compute_stellar_source_increment(table,pop,0d0,.07d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,early,ierr)
  if(ierr/=0)stop 13
  call compute_stellar_source_increment(table,pop,.07d0,.15d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,late,ierr)
  if(ierr/=0.or.abs(whole%returned_mass-early%returned_mass-late%returned_mass)>1d-9)stop 14
  if(abs(ledger%living_mass+ledger%returned_mass+ledger%remnant_mass-10000d0)>1d-9)stop 15
  open(newunit=u,file=trim(snia),status='old',action='read')
  call read_snia_population_realization_namelist(u,dtd,ierr)
  if(ierr/=0)stop 16
  call read_snia_physical_contract_namelist(u,event,ierr)
  close(u)
  if(ierr/=0)stop 17
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.15d0,1d0,8d0,64,ierr)
  if(ierr/=0)stop 18 ! Some WDs already formed at 20 Myr, before the 40 Myr DTD.
  pop%birth_metallicity=.02d0
  call integrate_ssp_channel(table,pop,2,.15d0,1d0,8d0,64,a,ierr)
  if(ierr/=0)stop 19
  call evaluate_snia_interval_events(dtd,pop%initial_mass,0d0,.15d0,1d0,events,ierr)
  if(ierr/=0)stop 20
  call build_snia_event_budget(event,events,a%remnant_mass,budget,ierr)
  if(ierr/=0)stop 21 ! End-of-step WD supply alone would incorrectly pass.
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.15d0,1d0,8d0,64,ierr)
  if(ierr/=enrichment_driver_err_ledger)stop 22
  call check_agb_wd_causality(table,pop,dtd,event%wd_debit_per_event,0d0,.07d0,1d0,8d0,64,ierr)
  if(ierr/=enrichment_driver_err_ledger)stop 23 ! Both coarse and fine steps reject.
  legacy=raw
  call prepare_high_mass_history(legacy,trim(old_history),ierr)
  if(ierr/=0.or.allocated(legacy%agb_terminal_row))stop 24
  call evaluate_channel_cumulative(legacy,2,8d0,0d0,.01d0,a,ierr)
  if(ierr/=0.or.a%remnant_mass<=0)stop 25 ! Legacy explicit linear policy is preserved.
  bad=raw
  do i=1,bad%n_rows
     if(bad%channel(i)==2.and.bad%age_gyr(i)==.02d0.and.bad%initial_mass(i)==8.and. &
          bad%birth_metallicity(i)==0)then
        bad%returned_mass(i)=.5d0*bad%returned_mass(i);bad%remnant_mass(i)=.5d0*bad%remnant_mass(i)
        bad%energy(i)=.5d0*bad%energy(i);bad%ejected_mass(i,:)=.5d0*bad%ejected_mass(i,:)
     endif
  enddo
  call prepare_high_mass_history(bad,trim(history),ierr)
  if(ierr==0.or.bad%high_mass_ready.or.allocated(bad%agb_terminal_row))stop 26
  bad=raw
  where(bad%channel==2.and.bad%remnant_mass>0)bad%remnant_mass=bad%remnant_mass-.01d0
  call prepare_high_mass_history(bad,trim(history),ierr)
  if(ierr==0.or.bad%high_mass_ready.or.allocated(bad%agb_terminal_row))stop 27
  call clear_yield_table(table)
  if(allocated(table%agb_terminal_row))stop 28
  print *, 'AGB_TERMINAL_TIMING_AND_SNIA_CAUSALITY_TEST_OK',events,budget%returned_mass
end program agb_terminal_test

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
  type(stellar_yield_table_t)::table,co_table
  type(stellar_cumulative_t)::a,b
  type(stellar_source_t)::whole,early,late
  type(stellar_population_t)::pop
  type(stellar_population_ledger_t)::ledger
  type(snia_population_realization_t)::dtd
  type(snia_physical_contract_t)::event
  type(snia_event_budget_t)::budget
  character(len=1024)::yields,history,snia,source_option
  real(stellar_dp)::energy,first_age
  real(stellar_dp)::times(9)=[0d0,.04d0,.05d0,.06889d0,.1d0,.2d0,1d0,5d0,13.7d0]
  ! Common KL16/LC18 Z support, not the broader AGB-only domain.
  real(stellar_dp),allocatable::zs(:)
  real(stellar_dp)::prior,expected_prior,events,remaining,generic_left,all_return,ia_return
  integer::i,r,ierr,u,iz,it,expected_agb_rows,n_exploded
  logical::full_ccsn,lowz7,pulses
  call get_command_argument(1,yields);call get_command_argument(2,history);call get_command_argument(3,snia)
  call get_command_argument(4,source_option)
  expected_agb_rows=58
  pulses=trim(source_option)=='lc18_set_r_agb7_pulses'
  lowz7=trim(source_option)=='lc18_set_r_agb7_lowz_net'.or.pulses
  full_ccsn=trim(source_option)=='lc18_set_r_fishlock'.or.trim(source_option)=='lc18_set_r_phase_fishlock'.or. &
       trim(source_option)=='lc18_set_r_agb7'.or.trim(source_option)=='lc18_set_r_agb7_net'.or.lowz7
  zs=[.007d0,.01d0,.012d0,.01345d0]
  select case(trim(source_option))
  case('')
  case('fishlock2014_raiteri96','lc18_set_r_fishlock','lc18_set_r_phase_fishlock')
     expected_agb_rows=73
     zs=[.001d0,.004d0,.007d0,.01d0,.012d0,.01345d0]
  case('lc18_set_r_agb7','lc18_set_r_agb7_net')
     expected_agb_rows=61
  case('lc18_set_r_agb7_lowz_net','lc18_set_r_agb7_pulses')
     expected_agb_rows=77
     zs=[.001d0,.004d0,.007d0,.01d0,.012d0,.01345d0]
  case default
     error stop 'unknown source option'
  end select
  call set_enrichment_defaults()
  default_imf_id=1;population_model_id=1;configured_binary_fraction=.5d0
  enable_agb=.true.;enable_snia=.true.;high_mass_model='wind_only_collapse'
  configured_channel_mass_min(1)=40;configured_channel_mass_min(3)=40
  if(full_ccsn)then
     configured_channel_mass_min(1)=13;configured_channel_mass_min(3)=13
  endif
  configured_channel_mass_max(2)=6
  if(expected_agb_rows==61.or.lowz7)configured_channel_mass_max(2)=7
  call load_yield_table(trim(yields),table,ierr)
  if(ierr/=0)stop 1
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0.or..not.allocated(table%agb_terminal_row))stop 2
  if(size(table%agb_terminal_row)/=expected_agb_rows.or..not.table%net_yield_diagnostic_unavailable)stop 3
  call set_yield_mass_assignment_mode(table,yield_mass_assignment_piecewise_constant,ierr)
  call audit_yield_table(table,1d-10,ierr,.true.,channel_owns_terminal_remnant, &
       [.true.,.true.,.true.,.false.,.false.])
  if(ierr/=0)stop 4
  if(trim(source_option)=='lc18_set_r_agb7_net'.or.lowz7)then
     if(any(table%net_yield_channel_available.neqv.[.false.,.true.,.false.,.false.,.false.]))stop 49
     if(.not.any(table%net_yield<0).or..not.any(table%net_yield>0))stop 50
     print *, 'AGB_SIGNED_NET_CHANNEL_AVAILABILITY_OK'
  endif
  if(expected_agb_rows==61.or.lowz7)then
     if(count(table%agb_remnant_kind/=0)/=merge(2,1,lowz7))stop 44
     co_table=table;co_table%co_wd_inventory_only=.true.
     call evaluate_channel_cumulative(table,2,7d0,.007d0,.1d0,a,ierr)
     if(ierr/=0.or.a%remnant_mass<=0)stop 45
     call evaluate_channel_cumulative(co_table,2,7d0,.007d0,.1d0,b,ierr)
     if(ierr/=0.or.b%remnant_mass/=0.or.b%returned_mass/=a%returned_mass)stop 46
     call evaluate_channel_cumulative(co_table,2,7d0,.014d0,.1d0,b,ierr)
     if(ierr/=0.or.b%remnant_mass<=0)stop 47
     prior=b%remnant_mass
     call evaluate_channel_cumulative(co_table,2,7d0,.01d0,.1d0,b,ierr)
     if(ierr/=0.or.abs(b%remnant_mass-prior*3d0/7d0)>1d-12)stop 48
     print *, 'AGB_HYBRID_ENVELOPE_INCLUDED_CO_WD_INVENTORY_EXCLUDED'
     if(lowz7)then
        if(count(table%agb_remnant_kind==2)/=1)stop 52
        call evaluate_channel_cumulative(table,2,7d0,.001d0,.1d0,a,ierr)
        if(ierr/=0.or.abs(a%returned_mass-5.855482487283754d0)>1d-12)stop 53
        if(abs(a%remnant_mass-1.144517512716246d0)>1d-12)stop 54
        if(a%net_yield(1)>=0.or.a%net_yield(2)<=0)stop 55
        call evaluate_channel_cumulative(co_table,2,7d0,.001d0,.1d0,b,ierr)
        if(ierr/=0.or.b%remnant_mass/=0.or.b%returned_mass/=a%returned_mass)stop 56
        call evaluate_channel_cumulative(co_table,2,7d0,.004d0,.1d0,b,ierr)
        if(ierr/=0.or.b%remnant_mass/=0.or.b%returned_mass<=0)stop 57
        call evaluate_channel_cumulative(table,2,7d0,.00099d0,.1d0,b,ierr)
        if(ierr==0)stop 58
        call evaluate_channel_cumulative(table,2,7.001d0,.001d0,.1d0,b,ierr)
        if(ierr==0)stop 59
        print *, 'FISHLOCK_ONE_ENVELOPE_NET_AND_CO_INVENTORY_OK',a%returned_mass,a%net_yield(:2)
     endif
  endif
  if(full_ccsn)then
     if(size(table%hm_mass)/=36)stop 29
     n_exploded=0
     do i=1,size(table%hm_mass)
        r=table%hm_terminal_row(i)
        call evaluate_channel_cumulative(table,3,table%hm_mass(i),table%hm_z(i), &
             nearest(table%hm_age(i),-1d0),a,ierr)
        if(ierr/=0.or.a%returned_mass/=0.or.a%energy/=0.or.a%remnant_mass/=0)stop 30
        call evaluate_channel_cumulative(table,3,table%hm_mass(i),table%hm_z(i),table%hm_age(i),a,ierr)
        if(ierr/=0.or.a%returned_mass/=table%returned_mass(r).or.a%energy/=table%energy(r))stop 31
        if(any(a%ejected_mass/=table%ejected_mass(r,:)))stop 32
        if(table%hm_mass(i)<=25)then
           if(a%returned_mass<=0.or.a%energy/=1d51)stop 33
           n_exploded=n_exploded+1
        else
           if(a%returned_mass/=0.or.a%energy/=0)stop 34
        endif
        call evaluate_channel_cumulative(table,1,table%hm_mass(i),table%hm_z(i),table%hm_age(i),b,ierr)
        if(ierr/=0.or.abs(a%returned_mass+a%remnant_mass+b%returned_mass-table%hm_mass(i))>1d-10)stop 35
     enddo
     if(n_exploded/=16)stop 36
     do i=1,2
        call evaluate_channel_cumulative(table,i*2-1,12.999d0,.01d0,.1d0,a,ierr)
        if(ierr==0)stop 37
     enddo
     print *, 'LC18_ORDINARY_CCSN_NATIVE_ENDPOINTS_OK',n_exploded
     if(trim(source_option)=='lc18_set_r_phase_fishlock')then
        ! Independent Table5 checkpoints, nonrotating 20 Msun, solar LC18 Z.
        call evaluate_channel_cumulative(table,1,20d0,.01345d0,68200d-9,a,ierr)
        if(ierr/=0.or.a%returned_mass/=0)stop 41
        call evaluate_channel_cumulative(table,1,20d0,.01345d0,9158200d-9,a,ierr)
        if(ierr/=0)stop 42
        call evaluate_channel_cumulative(table,1,20d0,.01345d0,.1d0,b,ierr)
        if(ierr/=0.or.abs(a%returned_mass/b%returned_mass-.5d0/12.46d0)>1d-12)stop 43
        print *, 'LC18_TABLE5_WIND_RELEASE_OK',a%returned_mass,b%returned_mass
     endif
  endif
  first_age=huge(1d0)
  if(pulses)then
     if(count(table%agb_terminal_jump_fraction<1)/=16)stop 60
     ! Published Mtot at the first/last TP: 6.9846 / 2.1117 Msun.
     ! Last TP aligned to selected Padova age, not a recovered Monash age.
     call evaluate_channel_cumulative(table,2,7d0,.001d0,(53248752.07686006d0-75687.63d0)*1d-9,a,ierr)
     if(ierr/=0.or.abs(a%returned_mass-.0154d0)>1d-9.or.a%remnant_mass/=0)stop 61
     call evaluate_channel_cumulative(table,2,7d0,.001d0,nearest(53248752.07686006d-9,-1d0),a,ierr)
     if(ierr/=0.or.abs(a%returned_mass-4.8883d0)>1d-8.or.a%remnant_mass/=0)stop 62
     print *, 'FISHLOCK_TP_WIND_LEFT_LIMIT_WITHOUT_EARLY_WD',a%returned_mass
  endif
  do i=1,size(table%agb_terminal_row)
     r=table%agb_terminal_row(i);first_age=min(first_age,table%age_gyr(r))
     call evaluate_channel_cumulative(table,2,table%initial_mass(r),table%birth_metallicity(r), &
          nearest(table%age_gyr(r),-1d0),a,ierr)
     if(ierr/=0.or.a%remnant_mass/=0)stop 5
     if(pulses)then
        if(abs(a%returned_mass-(1-table%agb_terminal_jump_fraction(i))*table%returned_mass(r))>1d-8)stop 63
     else
        if(a%returned_mass/=0)stop 5
     endif
     call evaluate_channel_cumulative(table,2,table%initial_mass(r),table%birth_metallicity(r), &
          table%age_gyr(r),a,ierr)
     if(ierr/=0.or.a%returned_mass/=table%returned_mass(r).or.any(a%ejected_mass/=table%ejected_mass(r,:)))stop 6
     if(any(a%net_yield/=table%net_yield(r,:)))stop 51
     if(abs(a%returned_mass+a%remnant_mass-table%initial_mass(r))>1d-12)stop 7
     energy=.5d0*a%returned_mass*1.98847d33*(15d5)**2
     if(abs(a%energy-energy)>1d-14*energy.or.any(a%momentum/=0))stop 8
  enddo
  pop%initial_mass=10000;pop%current_mass=10000;pop%birth_metallicity=.01d0
  if(expected_agb_rows==73.or.lowz7)pop%birth_metallicity=.004d0
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
  if(full_ccsn)then
     if(whole%channel_returned_mass(3)<=0.or.whole%channel_energy(3)<=0)stop 38
     if(abs(whole%energy-early%energy-late%energy)>1d-11*whole%energy)stop 39
     if(any(abs(whole%ejected_mass-early%ejected_mass-late%ejected_mass)>1d-9))stop 40
     print *, 'CCSN_SSP_RETURN_ENERGY',whole%channel_returned_mass(3),whole%channel_energy(3)
  endif
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
  ! The common hull is the intersection of active channel supports; neither
  ! the lower Fishlock edge nor the upper LC18 edge may be extrapolated.
  if(expected_agb_rows==73)then
     do iz=1,2
        pop%current_mass=pop%initial_mass
        if(iz==1)pop%birth_metallicity=.00099d0
        if(iz==2)pop%birth_metallicity=.01346d0
        call compute_stellar_source_increment(table,pop,0d0,.2d0,configured_channel_mass_min, &
             configured_channel_mass_max,64,whole,ierr)
        if(ierr==0)stop 28
     enddo
     print *, 'FISHLOCK_COMMON_Z_HULL_NO_EXTRAPOLATION_OK'
  endif
  print *, 'KL16_LC18_EFFECTIVE_SSP_TEST_OK'
end program kl16_lc18_native_test

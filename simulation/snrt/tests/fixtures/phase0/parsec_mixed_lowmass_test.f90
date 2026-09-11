! Actual-package native consumer test; not a production or AMR evolution test.
program parsec_mixed_lowmass_test
  use stellar_enrichment_config
  use stellar_enrichment_contract
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_provider
  use stellar_ssp_sources
  use stellar_enrichment_driver
  use stellar_population_ledger
  use snrt_spectral_contract
  use snrt_parsec_source
  use snrt_stellar_source
  use dust_mass_physics, only: dust_mass_enabled,dust_mass_model,dust_condensation
  implicit none
  type(stellar_yield_table_t)::table,trial
  type(stellar_population_t)::population
  type(stellar_cumulative_t)::state,states(5)
  type(stellar_source_t)::whole,early,late
  type(stellar_population_ledger_t)::ledger
  character(len=1024)::root,mode,path
  integer::status,u,i,c,k,r,checks=0
  real(stellar_dp)::total,rem,owner(5),fraction,t,stop,death,z,m
  real(stellar_dp)::q(9),e(9),qa(9),ea(9),qb(9),eb(9),ql(9),el(9),qh(9),eh(9)
  real(stellar_dp),allocatable::wa(:),wb(:),identity(:),again(:)
  call get_command_argument(1,root);call get_command_argument(2,mode)
  call set_enrichment_defaults()
  path=trim(root)//'/enrichment.nml'
  open(newunit=u,file=trim(path),status='old',action='read',iostat=status)
  call require(status==0,'open actual comparison config')
  call read_enrichment_namelist(u,status);close(u)
  call require(status==0,'native comparison config')
  call snrt_spectral_contract_load_from_environment(status)
  call require(status==0,'native spectral contract')
  if(mode=='wrong_version')then
     call parsec_sed_load(trim(root)//'/nodes.dat',status,expected_version=4)
     call require(status/=0.and..not.parsec_sed_enabled,'v4 rejects actual mixed V2 nodes')
     call get_command_argument(3,path)
     if(len_trim(path)>0)then
        call parsec_sed_load(trim(path),status,expected_version=5)
        call require(status/=0.and..not.parsec_sed_enabled,'v5 rejects actual legacy V1 nodes')
     endif
  endif
  call stellar_sed_load(status)
  call require(status==0.and.stellar_sed_has_energy,'v5 wrapper loads actual Q/E package')
  call stellar_photon_interval(0d0,1d0,.02d0,1d0,q,status,e)
  call require(status/=0.and.all(q==0).and.all(e==0),'unbound source cannot publish')
  call load_yield_table(trim(root)//'/yields.dat',table,status)
  call require(status==0,'actual material load')
  do k=1,3
     trial=table
     if(k==1)then
        do r=1,trial%n_rows
           if(trial%channel(r)==2.and.trial%remnant_mass(r)>0)exit
        enddo
        trial%remnant_mass(r-1)=trial%remnant_mass(r)
        trial%remnant_mass(r)=0 ! premature remnant, unchanged terminal budget
     else if(k==2)then
        do r=1,trial%n_rows
           if(trial%channel(r)==2.and.trial%returned_mass(r)>0)exit
        enddo
        trial%ejected_mass(r,1)=trial%returned_mass(r)+1d0
     else
        enable_snia=.true.
     endif
     call prepare_high_mass_history(trial,trim(root)//'/history.nml',status)
     call require(status/=0.and..not.trial%high_mass_ready.and..not.allocated(trial%hm_mass), &
          'malformed material or enabled SNIa rejects transactionally')
     enable_snia=.false.
  enddo
  call prepare_high_mass_history(table,trim(root)//'/history.nml',status)
  call require(status==0,'actual v5 history admission')
  call require(table%high_mass_version==5.and.size(table%hm_mass)==979,'all eleven actual Z branches')
  do r=1,table%n_rows
     if(table%initial_mass(r)>=14.or.table%channel(r)/=1)cycle
     m=10d0
     if(table%initial_mass(r)>=9)m=1000d0
     total=.5d0*table%returned_mass(r)*1.98847d33*(m*1d5)**2
     call require(abs(table%energy(r)-total)<=1d-12*max(total,tiny(1d0)), &
          'wind kinetic energy uses scalar-reduced wind mass')
  enddo
  dust_mass_enabled=.true.;dust_mass_model='carbon_olivine_2size_v1'
  dust_condensation=[.1d0,.2d0,.15d0]
  call prepare_dust_yields(table,status)
  call require(status==0,'per-release condensation')
  call parsec_sed_bind(table,status)
  call require(status==0,'shared mixed population binding')
  ! An explicit effective population reuses the SAME ordinary-star history;
  ! enabling Ia in the old single-star population remains rejected above.
  population_model_id=population_effective_ssp;enable_snia=.true.
  call require(user_source_model_requested(),'effective v5 config admits full AGB source window')
  call parsec_sed_bind(table,status)
  call require(status==0,'effective SSP shares unchanged ordinary radiation')
  call load_yield_table(trim(root)//'/yields.dat',trial,status)
  call require(status==0,'effective material reload')
  call prepare_high_mass_history(trial,trim(root)//'/history.nml',status)
  call require(status==0.and.trial%high_mass_ready,'effective SSP admits single-track v5 history')
  population_model_id=population_binary_ssp;configured_binary_fraction=.5d0
  call parsec_sed_bind(table,status)
  call require(status/=0,'resolved binary identity cannot masquerade as v5')
  population_model_id=population_single_star_ssp;configured_binary_fraction=0;enable_snia=.false.
  call parsec_sed_bind(table,status)
  call require(status==0,'original single-star binding restored')
  call stellar_sed_report()
  population%initial_mass=1;population%current_mass=1;population%birth_metallicity=.02d0
  population%imf_id=default_imf_id;population%imf_mass_min=configured_imf_mass_min
  population%imf_mass_max=configured_imf_mass_max;population%population_id=population_model_id
  population%yield_basis_id=yield_source_basis_id;population%pisn_enabled=.true.
  call mixed_source_weights(table,population,32,wa,status)
  call require(status==0,'32-bin shared weights')
  call mixed_source_weights(table,population,64,wb,status)
  call require(status==0.and.maxval(abs(wa-wb))<1d-13,'nearest source cells independent of base bins')
  call calculate_imf_mass_fraction(default_imf_id,.08d0,600d0,2d0,600d0,fraction,status)
  call require(status==0,'full IMF denominator')
  call require(abs(sum(wb*table%hm_mass,mask=table%hm_z==.02d0)-fraction)<1d-12, &
       'supported source fraction is not renormalized to one')
  do i=1,size(table%hm_mass)
     total=0;rem=0;owner=0;m=table%hm_mass(i);z=table%hm_z(i);t=table%hm_age(i)
     do c=1,5
        if(c==4)cycle
        call evaluate_channel_cumulative(table,c,m,z,t,state,status)
        call require(status==0,'actual node cumulative endpoint')
        total=total+state%returned_mass;rem=rem+state%remnant_mass;owner(c)=state%remnant_mass
        call require(state%returned_mass>=sum(state%ejected_mass)-1d-12*m,'gross plus untracked budget')
        if(c/=1)then
           call evaluate_channel_cumulative(table,c,m,z,nearest(t,-1d0),state,status)
           call require(status==0.and.state%returned_mass==0.and.state%remnant_mass==0, &
                'terminal material/remnant not created prematurely')
        endif
     enddo
     call require(abs(total+rem-m)<1d-10*m,'each node lifetime budget closes once')
     if(table%hm_terminal_channel(i)==2)then
        call require(owner(2)>0.and.sum(owner([1,3,4,5]))==0,'AGB remnant belongs only to channel2')
     else
        call require(sum(owner([1,2,4,5]))==0,'massive remnant belongs only to channel3')
     endif
  enddo
  call compute_stellar_source_increment(table,population,0d0,.005148513184858701d0, &
       configured_channel_mass_min,configured_channel_mass_max,64,early,status)
  call require(status==0,'actual young-age native increment')
  print '(A,5ES23.14)','YOUNG_RETURN_PER_INITIAL_MSUN ',early%channel_returned_mass
  print '(A,5ES23.14)','YOUNG_ENERGY_ERG_PER_INITIAL_MSUN ',early%channel_energy
  call require(early%channel_energy(1)>0.and.early%channel_energy(5)>0,'upper wind and pair sources preserved')
  call require(early%channel_energy(3)==0,'no fabricated early ordinary SNII')
  t=maxval(table%hm_age)*1.01d0
  call compute_stellar_source_increment(table,population,0d0,t, &
       configured_channel_mass_min,configured_channel_mass_max,64,whole,status,ledger)
  call require(status==0,'full lifetime native source and population ledger')
  call require(whole%channel_returned_mass(2)>0.and.whole%channel_returned_mass(3)>0, &
       'AGB and ordinary SNII both reach actual native consumer')
  call require(sum(whole%dust_species)>0,'actual native dust subset')
  call require(abs(ledger%returned_mass+ledger%remnant_mass-fraction)<1d-10,'covered mass closes across channels')
  call require(abs(ledger%living_mass-(1-fraction))<1d-10,'below2 stays living with full .08--600 normalization')
  call require(abs(ledger%unresolved_initial_mass-(1-fraction))<1d-10,'below2 explicitly unresolved')
  print '(A,5ES23.14)','LIFETIME_RETURN_PER_INITIAL_MSUN ',whole%channel_returned_mass
  print '(A,5ES23.14)','LIFETIME_REMNANT_PER_INITIAL_MSUN ',ledger%channel_remnant_mass
  print '(A,5ES23.14)','LIFETIME_ENERGY_ERG_PER_INITIAL_MSUN ',whole%channel_energy
  print '(A,4ES23.14)','LEDGER_RETURN_REMNANT_LIVING_UNRESOLVED ', &
       ledger%returned_mass,ledger%remnant_mass,ledger%living_mass,ledger%unresolved_initial_mass
  call compute_stellar_source_increment(table,population,0d0,.02d0, &
       configured_channel_mass_min,configured_channel_mass_max,64,early,status)
  call require(status==0,'first increment')
  call compute_stellar_source_increment(table,population,.02d0,t, &
       configured_channel_mass_min,configured_channel_mass_max,64,late,status)
  call require(status==0,'second increment')
  call require(abs(early%returned_mass+late%returned_mass-whole%returned_mass)<1d-12,'mass telescopes')
  call require(abs((early%energy+late%energy)/whole%energy-1)<1d-12,'energy telescopes')
  call require(maxval(abs(early%dust_species+late%dust_species-whole%dust_species))<1d-12,'dust telescopes')
  call stellar_photon_interval(0d0,t*1000,.02d0,1d0,q,status,e)
  call require(status==0.and.all(q>0).and.all(e>0),'actual full preterminal Q/E')
  call require(all(e>=q*snrt_group_edges_ev(:9)).and.all(e<=q*snrt_group_edges_ev(2:)), 'Q/E band support')
  call stellar_photon_interval(0d0,20d0,.02d0,1d0,qa,status,ea)
  call require(status==0,'early radiation interval')
  call stellar_photon_interval(20d0,t*1000,.02d0,1d0,qb,status,eb)
  call require(status==0,'late radiation interval')
  call require(maxval(abs((qa+qb)/q-1))<1d-12.and.maxval(abs((ea+eb)/e-1))<1d-12,'Q/E telescopes')
  print '(A,2ES23.14)','LIFETIME_Q_EEV_PER_INITIAL_MSUN ',sum(q),sum(e)
  stop=maxval(table%hm_radiation_stop);death=maxval(table%hm_age)
  call stellar_photon_interval(stop*1000*(1+1d-12),death*1000+1d0,.02d0,1d0,q,status,e)
  call require(status==0.and.all(q==0).and.all(e==0),'zero after source tracks, no late SED plateau')
  print '(A,2ES23.14)','GLOBAL_LAST_TRACK_AND_TERMINAL_GYR ',stop,death
  call stellar_photon_interval(0d0,t*1000,.017d0,1d0,ql,status,el)
  call require(status==0,'lower physical radiation Z node')
  call stellar_photon_interval(0d0,t*1000,.02d0,1d0,qh,status,eh)
  call require(status==0,'upper physical radiation Z node')
  call stellar_photon_interval(0d0,t*1000,.0185d0,1d0,q,status,e)
  call require(status==0.and.maxval(abs(.5d0*(ql+qh)/q-1))<1d-12.and. &
       maxval(abs(.5d0*(el+eh)/e-1))<1d-12,'same-age target-Z radiation mixture')
  population%birth_metallicity=.017d0
  call compute_stellar_source_increment(table,population,0d0,t, &
       configured_channel_mass_min,configured_channel_mass_max,64,early,status)
  call require(status==0,'lower physical material Z node')
  population%birth_metallicity=.0185d0
  call compute_stellar_source_increment(table,population,0d0,t, &
       configured_channel_mass_min,configured_channel_mass_max,64,late,status)
  call require(status==0.and.abs(.5d0*(early%returned_mass+whole%returned_mass)-late%returned_mass)<1d-12, &
       'same-age target-Z material mixture')
  call stellar_sed_identity(identity)
  do k=1,3
     trial=table
     select case(k)
     case(1)
        trial%hm_terminal_channel(1)=3
     case(2)
        trial%hm_remnant_kind(1)=3
     case(3)
        trial%hm_radiation_stop(1)=trial%hm_radiation_stop(1)*.99d0
     end select
     call parsec_sed_bind(trial,status)
     call require(status/=0,'changed v5 metadata cannot bind')
     call stellar_photon_interval(0d0,1d0,.02d0,1d0,q,status,e)
     call require(status/=0.and.all(q==0).and.all(e==0),'failed binding cannot publish')
  enddo
  call parsec_sed_bind(table,status)
  call require(status==0,'valid binding restored')
  call stellar_sed_identity(again)
  call require(size(identity)==size(again).and.all(identity==again),'deterministic identity')
  trial=table;trial%co_wd_inventory_only=.true.
  do i=1,size(table%hm_mass)
     if(table%hm_remnant_kind(i)/=3)cycle
     call evaluate_channel_cumulative(trial,2,table%hm_mass(i),table%hm_z(i),table%hm_age(i),state,status)
     call require(status==0.and.state%remnant_mass==0,'unclassified proxy never supplies CO-Ia inventory')
  enddo
  population%birth_metallicity=.001d0
  call compute_stellar_source_increment(table,population,0d0,.005148513184858701d0, &
       configured_channel_mass_min,configured_channel_mass_max,64,early,status)
  call require(status==0,'planned live Z=.001 young native material mixture')
  print '(A,5ES23.14)','Z001_YOUNG_ENERGY_ERG_PER_INITIAL_MSUN ',early%channel_energy
  call compute_stellar_source_increment(table,population,0d0,t, &
       configured_channel_mass_min,configured_channel_mass_max,64,whole,status,ledger)
  call require(status==0.and.abs(ledger%returned_mass+ledger%remnant_mass-fraction)<1d-10, &
       'planned live Z=.001 lifetime material mixture and ledger')
  call stellar_photon_interval(0d0,5.148513184858701d0,.001d0,1d0,q,status,e)
  call require(status==0.and.all(q>0).and.all(e>0),'planned live Z=.001 native Q/E mixture')
  print '(A,2ES23.14)','Z001_YOUNG_Q_EEV_PER_INITIAL_MSUN ',sum(q),sum(e)
  print '(A,I0)','PARSEC_MIXED_ACTUAL_NATIVE_PASS checks=',checks
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       print '(A,A,A,I0)','FAIL ',label,' status=',status
       error stop 1
    endif
    checks=checks+1
  end subroutine
end program

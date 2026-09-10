program parsec_pair_history_test
  use stellar_enrichment_config
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_yield_interpolation, only: source_metallicity_bracket
  use stellar_yield_provider
  use stellar_enrichment_driver
  use stellar_enrichment_contract
  use stellar_population_ledger
  use stellar_source_increment, only: source_condensed_donors
  use cosmic_ray_physics, only: cr_enabled,cr_sn_fraction,cr_source_partition
  implicit none
  type(stellar_yield_table_t)::table,raw,trial,baseline
  type(stellar_cumulative_t)::wind,cc,pair,before,a,b,mixed
  type(stellar_population_t)::population
  type(stellar_source_t)::whole,early,late,refined,repeat
  type(stellar_population_ledger_t)::ledger
  character(len=1024)::path,history,baseline_path,baseline_history
  real(stellar_dp)::mass,z,age,scale
  integer::ierr,i,ch,k,seen(5),row,j
  real(stellar_dp),allocatable::zcheck(:)
  real(stellar_dp)::zl,zh,w
  call get_command_argument(1,path);call get_command_argument(2,history)
  call get_command_argument(3,baseline_path);call get_command_argument(4,baseline_history)
  call set_enrichment_defaults()
  enable_agb=.false.;enable_pisn=.true.
  configured_imf_mass_max=600
  configured_channel_mass_min([1,3,5])=14
  configured_channel_mass_max([1,3,5])=600
  call load_yield_table(trim(path),raw,ierr)
  call require(ierr==0,'load actual PARSEC package')
  call audit_yield_table(raw,1d-10,ierr,exact_source_coordinates=.true.)
  if(ierr/=0)print *,'RAW AUDIT STATUS',ierr
  call require(ierr==0,'raw table invariants')
  table=raw
  call prepare_high_mass_history(table,trim(history),ierr)
  if(ierr/=0)print *,'PREPARE STATUS',ierr
  call require(ierr==0.and.table%high_mass_ready,'prepare v4')
  call require(table%high_mass_version==4,'version')
  zcheck=pack(table%hm_z,[.true.,table%hm_z(2:)/=table%hm_z(:size(table%hm_z)-1)])
  seen=0
  do i=1,size(table%hm_mass)
     mass=table%hm_mass(i);z=table%hm_z(i);age=table%hm_age(i)
     seen(table%hm_fate(i))=seen(table%hm_fate(i))+1
     call evaluate_channel_cumulative(table,1,mass,z,age,wind,ierr)
     call require(ierr==0,'wind endpoint')
     call evaluate_channel_cumulative(table,3,mass,z,age,cc,ierr)
     call require(ierr==0,'CCSN/remnant endpoint')
     call evaluate_channel_cumulative(table,5,mass,z,age,pair,ierr)
     call require(ierr==0.and.pair%remnant_mass==0,'no second remnant')
     call require(abs(wind%returned_mass+cc%returned_mass+pair%returned_mass+cc%remnant_mass-mass)< &
          1d-10*mass,'wind+CC+pair+baryonic remnant closure')
     do ch=3,5,2
        call evaluate_channel_cumulative(table,ch,mass,z,nearest(age,-1d0),before,ierr)
        call require(ierr==0.and.before%returned_mass==0.and.before%energy==0.and. &
             before%remnant_mass==0,'no terminal age ramp')
     enddo
     select case(table%hm_fate(i))
     case(1)
        call require(cc%returned_mass>0.and.pair%returned_mass==0,'CCSN channel')
     case(2,5)
        call require(cc%returned_mass==0.and.pair%returned_mass==0,'failed/direct fate')
     case(3)
        call require(cc%returned_mass==0.and.pair%returned_mass>0.and.cc%remnant_mass>0,'PPISN channel')
     case(4)
        call require(cc%returned_mass==0.and.pair%returned_mass>0.and.cc%remnant_mass==0,'PISN channel')
     end select
  enddo
  if(maxval(table%hm_z)>.014d0)then
     do j=1,size(zcheck)-1
        do ch=1,5,2
           call evaluate_channel_cumulative(table,ch,160d0,zcheck(j),.0032d0,a,ierr)
           call require(ierr==0,'metal-rich lower source')
           call evaluate_channel_cumulative(table,ch,160d0,zcheck(j+1),.0032d0,b,ierr)
           call require(ierr==0,'metal-rich upper source')
           call evaluate_channel_cumulative(table,ch,160d0,.5d0*(zcheck(j)+zcheck(j+1)),.0032d0,mixed,ierr)
           call require(ierr==0,'metal-rich same-age mixture')
           call require(abs(mixed%returned_mass-.5d0*(a%returned_mass+b%returned_mass))<1d-11,'five-Z mass mixture')
           call require(maxval(abs(mixed%ejected_mass-.5d0*(a%ejected_mass+b%ejected_mass)))<1d-11, &
                'five-Z elemental mixture')
           call require(abs(mixed%energy-.5d0*(a%energy+b%energy))<1d-12*max(1d0,a%energy,b%energy), &
                'five-Z energy mixture')
        enddo
     enddo
     print *,'PARSEC_N_Z_MATERIAL_MIXTURE_PASS',size(zcheck)
  endif
  call require(all(seen>0),'all five real source fates exercised')
  if(len_trim(baseline_path)>0)then
     call load_yield_table(trim(baseline_path),baseline,ierr)
     call require(ierr==0,'load fixed-wind baseline')
     call prepare_high_mass_history(baseline,trim(baseline_history),ierr)
     call require(ierr==0,'prepare fixed-wind baseline')
     call require(table%high_mass_identity(1)/=baseline%high_mass_identity(1),'different wind model identities')
     do i=1,size(table%hm_mass)
        mass=table%hm_mass(i);z=table%hm_z(i);age=table%hm_age(i)
        do ch=1,5,2
           call evaluate_channel_cumulative(table,ch,mass,z,age,a,ierr)
           call require(ierr==0,'phase-wind endpoint')
           call evaluate_channel_cumulative(baseline,ch,mass,z,age,b,ierr)
           call require(ierr==0,'fixed-wind endpoint')
           call require(a%returned_mass==b%returned_mass.and.a%remnant_mass==b%remnant_mass, &
                'unchanged endpoint mass and remnant')
           call require(all(a%ejected_mass==b%ejected_mass),'unchanged endpoint elements')
           if(ch==1)then
              call require(a%energy>0.and.abs(a%energy/b%energy-1)>1d-6,'actual phase wind changes energy')
           else
              call require(a%energy==b%energy,'unchanged explosion energy')
           endif
        enddo
     enddo
  endif
  z=.011d0
  if(minval(zcheck)<.001d0)z=.001d0
  call source_metallicity_bracket(zcheck,z,zl,zh,w,ierr)
  call require(ierr==0,'actual reference Z bracket')
  do ch=1,5,2
     call evaluate_channel_cumulative(table,ch,160d0,zl,.0032d0,a,ierr)
     call require(ierr==0,'lower Z')
     call evaluate_channel_cumulative(table,ch,160d0,zh,.0032d0,b,ierr)
     call require(ierr==0,'upper Z')
     call evaluate_channel_cumulative(table,ch,160d0,z,.0032d0,mixed,ierr)
     call require(ierr==0,'interior Z')
     call require(abs(mixed%returned_mass-((1-w)*a%returned_mass+w*b%returned_mass))<1d-11,'same-age Z material')
     call require(abs(mixed%energy-((1-w)*a%energy+w*b%energy))<1d-12*max(1d0,a%energy,b%energy),'same-age Z energy')
  enddo
  population%initial_mass=10000;population%current_mass=10000
  population%birth_metallicity=.011d0;population%imf_id=default_imf_id
  population%imf_mass_min=.08d0;population%imf_mass_max=600
  population%yield_basis_id=yield_basis_per_star_cumulative;population%pisn_enabled=.true.
  call compute_stellar_source_increment(table,population,0d0,.02d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,whole,ierr,ledger)
  call require(ierr==0.and.whole%channel_energy(5)>0,'actual IMF-integrated P(P)ISN source')
  call require(abs(ledger%living_mass+ledger%returned_mass+ledger%remnant_mass-10000d0)<1d-9,'SSP closure')
  call compute_stellar_source_increment(table,population,0d0,.0032d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,early,ierr)
  call require(ierr==0,'early source')
  call compute_stellar_source_increment(table,population,.0032d0,.02d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,late,ierr)
  call require(ierr==0,'late source')
  call require(abs(early%returned_mass+late%returned_mass-whole%returned_mass)<1d-9,'time telescoping mass')
  call require(maxval(abs(early%ejected_mass+late%ejected_mass-whole%ejected_mass))<1d-9,'time telescoping elements')
  call require(abs(early%energy+late%energy-whole%energy)<1d-12*whole%energy,'time telescoping energy')
  call compute_stellar_source_increment(table,population,0d0,.02d0,configured_channel_mass_min, &
       configured_channel_mass_max,128,refined,ierr)
  call require(ierr==0,'refined IMF')
  call require(abs(refined%returned_mass-whole%returned_mass)<1d-9,'common IMF cell mass')
  call require(abs(refined%energy-whole%energy)<1d-12*whole%energy,'common IMF cell energy')
  call compute_stellar_source_increment(table,population,.02d0,.03d0,configured_channel_mass_min, &
       configured_channel_mass_max,64,repeat,ierr)
  call require(ierr==0.and.repeat%returned_mass==0.and.repeat%energy==0,'terminal plateau')
  call require(delayed_cooling_source_mass(whole)==sum(whole%channel_returned_mass([3,5])),'SN-class shock mass')
  cr_enabled=.true.;cr_sn_fraction=.1d0
  call cr_source_partition(sum(whole%channel_energy([3,5])),0d0,scale,ierr)
  call require(ierr==0.and.abs(scale-.1d0*sum(whole%channel_energy([3,5])))<1d-14*whole%energy,'SN-class CR partition')
  call clear_source(repeat)
  repeat%channel_returned_mass(3)=1;repeat%channel_returned_mass(5)=10
  repeat%channel_ejected_mass(3,[1,3,11])=[.9d0,.05d0,.05d0]
  repeat%channel_ejected_mass(5,[1,3,11])=[9d0,.5d0,.5d0]
  call source_condensed_donors(repeat,.01d0,.01d0,ierr)
  call require(ierr==0.and.all(repeat%channel_condensed_mass(5,:)==0),'no pair Fe/PAH direct condensation')
  do k=1,4
     trial=raw
     select case(k)
     case(1)
        row=table%hm_pair_row(findloc(table%hm_fate,4,dim=1));trial%remnant_mass(row)=1
     case(2)
        row=table%hm_pair_row(findloc(table%hm_fate,1,dim=1));trial%returned_mass(row)=1
     case(3)
        enable_pisn=.false.
     case(4)
        configured_channel_mass_max(5)=599
     end select
     call prepare_high_mass_history(trial,trim(history),ierr)
     call require(ierr/=0.and..not.trial%high_mass_ready,'invalid admission has no published history')
     enable_pisn=.true.;configured_channel_mass_max(5)=600
  enddo
  call evaluate_channel_cumulative(table,5,160d0,.5d0*minval(zcheck),.02d0,pair,ierr)
  call require(ierr/=0,'no Z extrapolation')
  print *,'PARSEC_PAIR_HISTORY_OK fates=',seen
  print *,'SSP_RETURN_REMNANT=',whole%returned_mass,ledger%remnant_mass
  print *,'SSP_CHANNEL_MASS=',whole%channel_returned_mass
  print *,'SSP_CHANNEL_ENERGY=',whole%channel_energy
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok)then
       print *,'FAIL ',label
       error stop 1
    endif
  end subroutine
end program

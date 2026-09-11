program fp2_snia_frozen_he_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use stellar_enrichment_config, only: stellar_dp, stellar_imf_kroupa, population_binary_ssp
  use stellar_snia_population_contract
  use stellar_snia_physical_contract
  use stellar_snia_runtime_accounting
  use stellar_snia_cell_deposition
  use stellar_native_units, only: solar_mass_cgs
  implicit none
  type(snia_population_realization_t)::pop,bad
  type(snia_physical_contract_t)::physical
  type(snia_thermal_coupling_t)::coupling
  type(snia_event_budget_t)::budget
  type(snia_cell_increment_t)::cell
  real(stellar_dp),allocatable::identity(:),changed(:)
  real(stellar_dp)::t,n,n1,n2,prior,remaining,restarted,old_mass,expected,volume
  real(stellar_dp),parameter::initial=1d5
  character(len=1024)::path
  integer::unit,ierr,i
  call get_command_argument(1,path)
  open(newunit=unit,file=trim(path),status='old',action='read')
  call read_snia_population_realization_namelist(unit,pop,ierr)
  call check(ierr==0,'native frozen event sidecar load')
  call read_snia_physical_contract_namelist(unit,physical,ierr)
  call check(ierr==0,'ordered unchanged N100 load')
  call read_snia_thermal_coupling_namelist(unit,coupling,ierr)
  close(unit)
  call check(ierr==0,'ordered thermal coupling load')
  call validate_snia_population_binding(pop,stellar_imf_kroupa,population_binary_ssp,.5d0,ierr)
  call check(ierr==0,'full-initial SSP population binding')
  call validate_snia_event_source(pop,physical%returned_mass_per_event,physical%wd_debit_per_event, &
       physical%terminal_remnant_per_event,physical%yield_source_id,physical%yield_source_sha256,ierr)
  call check(ierr==0,'N100 mass and source binding')
  call validate_snia_event_source(pop,1.38d0,1.38d0,0d0,physical%yield_source_id,physical%yield_source_sha256,ierr)
  call check(ierr/=0,'unfunded N100 threshold mismatch rejected')
  t=pop%event%age_yr*1d-9
  expected=initial*pop%event%weight
  call evaluate_snia_interval_events(pop,initial,0d0,t,1d0,n,ierr)
  call check(ierr/=0.and.n==0,'table requires actual birth Z')
  call evaluate_snia_interval_events(pop,initial,0d0,t,1d0,n,ierr,birth_metallicity=.02d0)
  call check(ierr/=0.and.n==0,'single-Z table rejects other Z')
  call evaluate_snia_interval_events(pop,initial,0d0,t,.5d0,n,ierr,birth_metallicity=.01d0)
  call check(ierr/=0.and.n==0,'multiplicative factor cannot stand in for birth Z')
  call evaluate_snia_interval_events(pop,initial,0d0,nearest(t,-1d0),1d0,n,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n==0,'no event immediately before threshold')
  call evaluate_snia_interval_events(pop,initial,0d0,t,1d0,n,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n==expected,'event on upper endpoint and no extra binary scaling')
  call evaluate_snia_interval_events(pop,initial,t,nearest(t,1d0),1d0,n2,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n2==0,'no replay at lower endpoint')
  call evaluate_snia_interval_events(pop,initial,t,t,1d0,n2,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n2==0,'zero-width interval')
  call evaluate_snia_interval_events(pop,initial,0d0,t/2,1d0,n1,ierr,birth_metallicity=.01d0)
  call check(ierr==0,'pre-event subinterval')
  call evaluate_snia_interval_events(pop,initial,t/2,13.7d0,1d0,n2,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n1+n2==expected,'split and lifetime cumulative use identical CDF')
  call evaluate_snia_interval_events(pop,initial,t,0d0,1d0,n2,ierr,birth_metallicity=.01d0)
  call check(ierr/=0.and.n2==0,'reverse interval rejects')

  call build_snia_event_budget(physical,n,.75d0*initial,budget,ierr)
  call check(ierr==0,'actual table event builds native N100 budget')
  call close(budget%returned_mass,n*pop%event%threshold,'funded threshold mass')
  call close(budget%energy,n*1.5063100005966762d51,'unchanged N100 event energy')
  call check(all(budget%ejected_mass==n*physical%ejected_mass_per_event),'unchanged N100 eleven elements')
  call close_effective_snia_return(initial,.25d0*initial,.75d0*initial,0d0,0d0, &
       budget%returned_mass,initial,1d-12,remaining,ierr)
  call check(ierr==0,'ordinary plus hybrid Ia closes effective SSP')
  call close(remaining+budget%returned_mass+.25d0*initial,initial,'combined mass conservation')
  volume=1d55
  call build_snia_cell_increment(budget,volume,[0d0,0d0,0d0],coupling,cell,ierr)
  call check(ierr==0,'native physical source-to-cell receiver')
  call close(cell%mass_density*volume/solar_mass_cgs,budget%returned_mass,'receiver mass conservation')
  call close(cell%total_energy_density*volume,budget%energy,'receiver energy conservation')

  old_mass=remaining
  call reconstruct_prior_snia_return(old_mass,0d0,.75d0*initial,initial,1d-12,prior,ierr)
  call check(ierr==0,'persisted particle mass reconstructs prior Ia return')
  call evaluate_snia_interval_events(pop,initial,0d0,t,1d0,n1,ierr,birth_metallicity=.01d0)
  call check(ierr==0,'restart prior exact CDF')
  call evaluate_snia_interval_events(pop,initial,t,.1d0,1d0,n2,ierr,birth_metallicity=.01d0)
  call check(ierr==0.and.n2==0,'restarted next interval has no event')
  call close_effective_snia_return(old_mass,0d0,.75d0*initial,prior,n1*physical%returned_mass_per_event, &
       n2*physical%returned_mass_per_event,initial,1d-12,restarted,ierr)
  call check(ierr==0.and.restarted==old_mass,'restart no duplicate debit')
  call close_effective_snia_return(old_mass,0d0,.75d0*initial,prior,0d0,0d0,initial,1d-12,restarted,ierr)
  call check(ierr/=0,'incorrect prior CDF rejects')
  call build_snia_event_budget(physical,n,budget%returned_mass/2,budget,ierr)
  call check(ierr/=0,'effective remaining mass shortfall rejects instead of clipping')

  call snia_event_model_identity(pop,identity)
  call check(size(identity)>0,'table carries versioned numerical and source identity')
  bad=pop;bad%mass_accounting=snia_accounting_strict_wd;bad%accounting_approval_id=''
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'new event model cannot masquerade as strict WD')
  bad=pop;bad%binary_fraction=.3d0
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'fixed grid weight cannot acquire different binary-fraction metadata')
  bad=pop;bad%event%age_yr=bad%event%age_yr+1
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'event age must follow prescribed growth')
  call snia_event_model_identity(bad,changed)
  call check(any(changed/=identity),'event age changes restart identity')
  bad=pop;bad%event%weight=bad%event%weight*2
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'conflicting event normalizations reject')
  call snia_event_model_identity(bad,changed)
  call check(any(changed/=identity),'weight changes restart identity')
  do i=1,9
     bad=pop
     bad%event%source_sha256(i)(1:1)='f'
     if(pop%event%source_sha256(i)(1:1)=='f')bad%event%source_sha256(i)(1:1)='e'
     call snia_event_model_identity(bad,changed)
     call check(any(changed/=identity),'each source/parameter/converter hash is bound')
  enddo
  bad=pop;bad%event%donor_background_loss=10
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'comparison donor fuel cannot overdraw')
  bad=pop;bad%event%transfer_rate=ieee_value(0d0,ieee_quiet_nan)
  call validate_snia_population_realization(bad,ierr)
  call check(ierr/=0,'nonfinite physical parameter rejects')
  bad=pop;bad%event_model=snia_event_empirical
  call snia_event_model_identity(bad,changed)
  call check(size(changed)==0,'empirical restart extension remains empty')
  write(*,'(A)')'FP2_SNIA_FROZEN_HE_NATIVE_COUPLING_PASS'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(ok)return
    write(*,*)'FAIL: ',label
    error stop 1
  end subroutine check
  subroutine close(a,b,label)
    real(stellar_dp),intent(in)::a,b
    character(len=*),intent(in)::label
    call check(abs(a-b)<=1d-12*max(abs(a),abs(b),tiny(1d0)),label)
  end subroutine close
end program fp2_snia_frozen_he_test

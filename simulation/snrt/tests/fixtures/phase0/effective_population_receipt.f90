! Native per-initial-mass source receipt for comparison with live particles.
program effective_population_receipt
  use stellar_enrichment_config
  use stellar_enrichment_contract
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_enrichment_driver
  use stellar_snia_population_contract
  use stellar_native_units,only:code_time_to_age_gyr
  implicit none
  type(stellar_yield_table_t)::table
  type(stellar_population_t)::p
  type(stellar_source_t)::s
  type(snia_population_realization_t)::ia
  character(len=1024)::arg
  integer::u,ierr
  real(stellar_dp)::age_code,age,events
  call get_command_argument(1,arg)
  open(newunit=u,file=trim(arg),status='old')
  call read_enrichment_namelist(u,ierr);close(u)
  if(ierr/=0)stop 1
  call get_command_argument(2,arg)
  call load_yield_table(trim(arg),table,ierr)
  if(ierr/=0)stop 2
  call prepare_high_mass_history(table,trim(high_mass_history_file),ierr)
  if(ierr/=0)stop 3
  call get_command_argument(3,arg)
  open(newunit=u,file=trim(arg),status='old')
  call read_snia_population_realization_namelist(u,ia,ierr);close(u)
  if(ierr/=0)stop 4
  call validate_snia_population_binding(ia,default_imf_id,population_model_id,configured_binary_fraction,ierr)
  if(ierr/=0)stop 5
  call get_command_argument(4,arg);read(arg,*)age_code
  call code_time_to_age_gyr(age_code,3.1556926d16,1d0,age,ierr)
  if(ierr/=0)stop 6
  p%initial_mass=1;p%current_mass=1;p%imf_id=default_imf_id
  p%population_id=population_model_id;p%binary_fraction=configured_binary_fraction
  p%imf_mass_min=configured_imf_mass_min;p%imf_mass_max=configured_imf_mass_max
  p%yield_basis_id=yield_source_basis_id;p%pisn_enabled=.true.
  call get_command_argument(5,arg);read(arg,*)p%birth_metallicity
  call compute_stellar_source_increment(table,p,0d0,age,configured_channel_mass_min, &
       configured_channel_mass_max,64,s,ierr)
  if(ierr/=0)stop 7
  call evaluate_snia_interval_events(ia,1d0,0d0,age,1d0,events,ierr,birth_metallicity=p%birth_metallicity)
  if(ierr/=0)stop 8
  write(*,'(a,3es26.17)')'RECEIPT ',age,s%returned_mass,events
end program

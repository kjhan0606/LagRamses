program dust_mass_smoke
  use dust_mass_physics
  use dust_composition_material
  use dust_element_cooling, only: cie_n,wss09_curve,wss09_step,wss09_rate
  use stellar_yield_tables, only: stellar_yield_table_t,prepare_dust_yields,clear_yield_table
  use stellar_yield_interpolation, only: interpolate_yield_row
  implicit none
  integer::ierr
  real(real64)::d,e,q,a,b,half,full
  real(real64)::elements(11),gas(11),grains(2),next(2),dust_at_half(2),dust_at_end(2),p(3),net(11),remnant
  type(stellar_yield_table_t)::table
  real(real64)::tau(cie_n),rates(cie_n),tau2(cie_n),rates2(cie_n),rates0(cie_n),delta,dt_cie,actual_t
  real(real64)::bins(4),new_bins(4),split_bins(4),pair(2),reference(2)
  real(real64)::uc(3),us(3),um(3),area_small,area_large
  if(.not.dust_mass_parameters_ok())stop 1
  dust_mass_enabled=.true.
  call dust_condense([1d0,2d0,3d0],[.7d0,1d0,1d0],[.2d0,.5d0,1d0],d,ierr)
  if(ierr/=0.or.abs(d-.25d0)>1d-14)stop 2
  call dust_mass_step(1d0,.1d0,2d0,1d0,0d0,d,ierr)
  if(ierr/=0.or.abs(d-1/(1+9*exp(-2d0)))>1d-14)stop 3
  call dust_mass_step(1d0,.1d0,2d0,0d0,1d0,d,ierr)
  if(ierr/=0.or.abs(d-.1d0*exp(-2d0))>1d-14)stop 4
  call dust_mass_step(1d0,.1d0,2d0,1d0,1d0,d,ierr)
  if(ierr/=0.or.abs(d-.1d0/1.2d0)>1d-14)stop 5
  call dust_mass_step(1d0,.1d0,1d0,.7d0,.3d0,half,ierr)
  call dust_mass_step(1d0,half,1d0,.7d0,.3d0,d,ierr)
  call dust_mass_step(1d0,.1d0,2d0,.7d0,.3d0,full,ierr)
  if(ierr/=0.or.abs(d-full)>1d-14)stop 6
  call dust_mass_step(1d0,.1d0,1d4,1d0,0d0,d,ierr)
  if(ierr/=0.or.d/=1d0)stop 7
  call dust_mass_step(1d0,0d0,1d4,1d0,0d0,d,ierr)
  if(ierr/=0.or.d/=0d0)stop 8
  call dust_mass_step(1d0,1.1d0,1d0,1d0,0d0,d,ierr)
  if(ierr==0)stop 9
  call dust_mass_exchange(.1d0,.2d0,2d0,10d0,e,q,ierr)
  if(ierr/=0.or.e/=4.or.q/=2.or.(10-q+e)/=12)stop 10
  call dust_mass_exchange(.1d0,0d0,2d0,10d0,e,q,ierr)
  if(ierr/=0.or.e/=0.or.q/=-2)stop 11
  call dust_mass_exchange(.1d0,.2d0,2d0,1d0,e,q,ierr)
  if(ierr==0)stop 12
  call dust_mass_rates(1d-24,1d-26,100d0,a,b,ierr)
  if(ierr/=0.or.a<=0.or.b<0)stop 13
  call dust_mass_rates(1d-24,1d-26,2d6,a,b,ierr)
  if(ierr/=0.or.a/=0.or.abs(b-3*3.2d-18*(1d-24/dust_mp)*.5d0/dust_grain_radius_cm)>1d-25)stop 14
  dust_mass_model='carbon_olivine_v1'
  elements=olivine_fraction;elements(3)=.5d0
  call dust_species_condense(elements,3,grains,ierr)
  if(ierr/=0.or.maxval(abs(grains-[.075d0,.15d0]))>1d-14)stop 15
  call dust_gas_elements(elements,grains,gas,ierr)
  if(ierr/=0.or.abs(sum(gas)+sum(grains)-sum(elements))>1d-14)stop 16
  elements=0;elements(3)=.24d0;elements(5)=.16d0
  call dust_species_condense(elements,2,grains,ierr)
  if(ierr/=0.or.abs(grains(1)-.024d0)>1d-14.or.grains(2)/=0)stop 17
  elements=olivine_fraction*.6d0;elements(3)=.2d0;grains=[.02d0,.06d0]
  call dust_species_step(10d0,1d0,elements,grains,100d0,1d-24,1d16,next,ierr)
  if(ierr/=0.or.any(next<=grains).or.any(next>[.2d0,.6d0]))stop 18
  call dust_gas_elements(elements,next,gas,ierr)
  if(ierr/=0.or.abs(sum(gas)+sum(next)-sum(elements))>1d-14)stop 19
  call dust_gas_elements(elements,[.3d0,0d0],gas,ierr)
  if(ierr==0)stop 20
  ! Distinct C-rich and O-rich AGB nodes: averaging C/O first erases C dust.
  table%n_rows=4;table%loaded=.true.
  allocate(table%channel(4),table%initial_mass(4),table%birth_metallicity(4),table%age_gyr(4), &
       table%returned_mass(4),table%remnant_mass(4),table%energy(4),table%momentum(4,3), &
       table%ejected_mass(4,11),table%net_yield(4,11))
  table%channel=2;table%initial_mass=[1d0,1d0,2d0,2d0]
  table%birth_metallicity=.01d0;table%age_gyr=[0d0,1d0,0d0,1d0]
  table%returned_mass=[0d0,1d0,0d0,2d0];table%remnant_mass=0;table%energy=0
  table%momentum=0;table%ejected_mass=0;table%net_yield=0
  table%ejected_mass(2,3)=1;table%ejected_mass(4,5)=2
  call prepare_dust_yields(table,ierr)
  if(ierr/=0)stop 21
  call interpolate_yield_row(table,2,1.5d0,.01d0,.5d0,d,remnant,e,p,elements,net,ierr,dust_at_half)
  if(ierr/=0.or.abs(dust_at_half(1)-.05d0)>1d-14.or.dust_at_half(2)/=0)stop 22
  call interpolate_yield_row(table,2,1.5d0,.01d0,1d0,d,remnant,e,p,elements,net,ierr,dust_at_end)
  if(ierr/=0.or.maxval(abs(2*dust_at_half-dust_at_end))>1d-14)stop 23
  call clear_yield_table(table)
  ! One physical element table, not a scalar-Z curve relabelled as species.
  elements=0;elements(1)=.74d0;elements(2)=.25d0;elements(3)=.01d0
  call wss09_curve(elements,tau,rates,ierr)
  if(ierr/=0)stop 24
  elements(3)=0;elements(11)=.01d0
  call wss09_curve(elements,tau2,rates2,ierr)
  if(ierr/=0.or.abs(rates(150)-rates2(150))<.1d0*abs(rates(150)))stop 25
  ! Pre-depletion and depleted curves share T nodes. Carbon lock removes
  ! its rate linearly, keeping H/He and Fe unchanged (WSS09 equation 3).
  elements(3)=.001d0;elements(11)=.009d0
  call wss09_curve(elements,tau,rates,ierr)
  grains=[.0005d0,0d0]
  call dust_gas_elements(elements,grains,gas,ierr)
  if(ierr/=0)stop 26
  call wss09_curve(gas,tau2,rates2,ierr)
  elements(3)=0
  call wss09_curve(elements,tau2,rates0,ierr)
  if(ierr/=0.or.maxval(abs(rates2-(rates+rates0)/2))>1d-12*maxval(abs(rates)))stop 33
  ! Time splitting through multiple piecewise-linear energy intervals.
  elements(3)=.001d0
  dt_cie=1d13
  call wss09_step(1.66d-27,elements,tau(150),dt_cie,5d0/3,delta,ierr)
  if(ierr/=0.or.delta>=0)stop 27
  full=tau(150)+delta
  call wss09_step(1.66d-27,elements,tau(150),dt_cie/2,5d0/3,delta,ierr)
  if(ierr/=0)stop 28
  half=tau(150)+delta
  call wss09_step(1.66d-27,elements,half,dt_cie/2,5d0/3,delta,ierr)
  if(ierr/=0.or.abs(half+delta-full)>1d-11*full)stop 29
  call wss09_rate(full,elements,d,actual_t,ierr)
  if(ierr/=0.or.d<=0.or.actual_t<100d0)stop 30
  write(*,'(A,3ES24.16)')'WSS09_REFERENCE_T2_INITIAL_FINAL_T ',tau(150),full,actual_t
  call wss09_step(1.66d-27,elements,1d0,1d0,5d0/3,delta,ierr)
  if(ierr==0)stop 31
  ! Signed low-T H/He NET heating is not clipped into cooling.
  elements(3:11)=0
  call wss09_curve(elements,tau,rates,ierr)
  call wss09_step(1.66d-27,elements,tau(1),1d13,5d0/3,delta,ierr)
  if(ierr/=0.or.delta<=0.or.rates(1)>=0)stop 32
  dust_mass_model='carbon_olivine_2size_v1'
  bins=dust_injection_bins([.2d0,.6d0])
  if(any(bins/=[0d0,.2d0,0d0,.6d0]))stop 34
  dust_small_injection_fraction=[.25d0,.75d0]
  bins=dust_injection_bins([.2d0,.6d0])
  if(maxval(abs(bins-[.05d0,.15d0,.45d0,.15d0]))>1d-16)stop 48
  dust_small_injection_fraction=0
  ! Exact 54 Myr donor lifetime at Dlarge/rho=.01, nH=1, density=3.
  dust_size_density=3d0;pair=[0d0,.01d0]
  call dust_size_exchange(1d0,1d0,100d0,1,54*dust_myr,pair,ierr)
  if(ierr/=0.or.maxval(abs(pair-[.005d0,.005d0]))>1d-16)stop 35
  reference=[0d0,.01d0]
  call dust_size_exchange(1d0,1d0,100d0,1,27*dust_myr,reference,ierr)
  call dust_size_exchange(1d0,1d0,100d0,1,27*dust_myr,reference,ierr)
  if(ierr/=0.or.maxval(abs(pair-reference))>1d-16)stop 36
  pair=[.01d0,0d0]
  call dust_size_exchange(1d0,1d3,100d0,2,.27d0*dust_myr,pair,ierr)
  if(ierr/=0.or.maxval(abs(pair-[.005d0,.005d0]))>1d-16)stop 37
  call dust_size_exchange(1d0,1d3,1d4,2,1d20,pair,ierr)
  if(ierr/=0.or.maxval(abs(pair-[.005d0,.005d0]))>1d-16)stop 38
  dust_size_density=[2.2d0,3.3d0];dust_coagulation=.false.;dust_shattering=.false.
  pair=[1d-30,1d0];reference=pair
  call dust_size_exchange(1d0,1d0,100d0,1,1d10,pair,ierr)
  if(ierr/=0.or.any(pair/=reference))stop 49
  dust_sputtering=.false.
  elements=olivine_fraction*.006d0;elements(3)=.002d0;elements(1)=.74d0;elements(2)=.25d0
  bins=[.0001d0,.0001d0,.0002d0,.0002d0]
  call dust_size_step(1d0,.01d0,elements,bins,100d0,1d-22,1d12,new_bins,ierr)
  if(ierr/=0.or.any(new_bins<=bins).or.new_bins(1)<=new_bins(2).or.new_bins(3)<=new_bins(4))stop 39
  call dust_gas_elements(elements,[sum(new_bins(1:2)),sum(new_bins(3:4))],gas,ierr)
  if(ierr/=0)stop 40
  call dust_size_step(1d0,.01d0,elements,bins,100d0,1d-22,5d11,split_bins,ierr)
  call dust_size_step(1d0,.01d0,elements,split_bins,100d0,1d-22,5d11,bins,ierr)
  if(ierr/=0.or.maxval(abs(bins-new_bins)/new_bins)>1d-6)stop 41
  dust_growth=.false.;dust_sputtering=.true.;bins=.0001d0
  call dust_size_step(1d0,.01d0,elements,bins,2d6,1d-24,1d9,new_bins,ierr)
  b=3*3.2d-18*(1d-24/dust_mp)*.5d0
  if(ierr/=0.or.maxval(abs(new_bins(1:2)/bins(1:2)-exp(-b*1d9/dust_size_radius_cm)))>1d-13)stop 42
  dust_growth=.true.;dust_coagulation=.true.;dust_shattering=.true.
  dust_small_injection_fraction(1)=1.01d0
  if(dust_mass_parameters_ok())stop 43
  dust_small_injection_fraction=0
  dust_sn_shocks=.true.;bins=.001d0
  call dust_shock_step(6800d0,1d51,bins,bins/2,new_bins,ierr)
  if(ierr/=0.or.maxval(abs(new_bins/bins-(1+exp(-[2d0,.1d0,3d0,.15d0]))/2))>1d-14)stop 44
  call dust_shock_step(6800d0,1d99,bins,bins,split_bins,ierr)
  if(ierr/=0.or.any(split_bins/=bins))stop 45
  ! At fixed gas mass, splitting the exposure composes exactly when no new
  ! grains are added between intervals. Fresh survival is NOT applied twice.
  call dust_shock_step(1d4,1d51,bins,0*bins,new_bins,ierr)
  call dust_shock_step(1d4,.5d51,bins,0*bins,split_bins,ierr)
  call dust_shock_step(1d4,.5d51,split_bins,0*bins,bins,ierr)
  if(ierr/=0.or.maxval(abs(bins-new_bins))>1d-16)stop 46
  call dust_shock_step(1d4,1d51,bins,bins*1.1d0,new_bins,ierr)
  if(ierr==0)stop 47
  dust_sn_shocks=.false.
  dust_material_model='dl01_composition_v1'
  call dust_composition_curve([10d0,20d0,100d0],[1d0,0d0],1d0,uc,ierr)
  if(ierr/=0)stop 50
  call dust_composition_curve([10d0,20d0,100d0],[0d0,1d0],1d0,us,ierr)
  if(ierr/=0.or.any(abs(uc-us)<1d-6*uc))stop 51
  call dust_composition_curve([10d0,20d0,100d0],[.3d0,.7d0],1d0,um,ierr)
  if(ierr/=0.or.maxval(abs(um-(.3d0*uc+.7d0*us))/um)>1d-14)stop 52
  call dust_composition_area([1d0,0d0,0d0,0d0],1d0,area_small,ierr)
  call dust_composition_area([0d0,1d0,0d0,0d0],1d0,area_large,ierr)
  if(ierr/=0.or.abs(area_small/area_large-20)>1d-13)stop 53
  call dust_composition_curve([1d0,20d0,100d0],[1d0,0d0],1d0,um,ierr)
  if(ierr==0)stop 54
  dust_material_model='fixed_mix'
  dust_mass_model='bulk_v1'
  write(*,*)'DUST_MASS_CONDENSATION_GROWTH_SPUTTERING_ENERGY_PASS'
  write(*,*)'DUST_COMPOSITION_ELEMENT_BUDGET_PREMIX_SOURCE_PASS'
  write(*,*)'WSS09_ELEMENT_DEPLETION_SIGNED_COOLING_TIMESTEP_PASS'
  write(*,*)'DUST_TWO_SIZE_BUDGET_GROWTH_EROSION_COAG_SHATTER_PASS'
  write(*,*)'DUST_SN_AMBIENT_EXPOSURE_FRESH_PROTECTION_PASS'
  write(*,*)'DL01_COMPOSITION_MATERIAL_SIZE_AREA_DOMAIN_PASS'
end program

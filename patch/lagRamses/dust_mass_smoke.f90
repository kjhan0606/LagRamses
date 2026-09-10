program dust_mass_smoke
  use dust_stochastic_physics
  use dust_pah_radiation
  use snrt_dust_ir, only: dust_ir_diagnostics
  use dust_drag_physics
  use dust_multifluid
  use dust_mass_physics
  use dust_composition_material
  use dust_iron_material
  use dust_sublimation_physics, only: dust_sublimation_step,dust_olivine_surface_flux
  use snrt_dust_ir, only: snrt_dust_material_temperature
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
  real(real64)::material_id(1+3*dl01_n)
  if(.not.dust_mass_parameters_ok())stop 1
  dust_pah_condensation=.2d0
  call dust_pah_condense(1d0,1d0,.25d0,d,ierr)
  if(ierr/=0.or.abs(d*dust_pah_hc(2)-.15d0)>1d-15)stop 201
  call dust_pah_condense(.001d0,1d0,.25d0,d,ierr)
  if(ierr/=0.or.abs(d*dust_pah_hc(1)-.001d0)>1d-16)stop 202
  call dust_pah_condense(1d0,1d0,1.1d0,d,ierr)
  if(ierr==0)stop 203
  dust_pah_condensation=0
  call dust_pah_condense(1d0,1d0,0d0,d,ierr)
  if(ierr/=0.or.d/=0)stop 204
  write(*,*)'PAH_NATIVE_NONIA_POST_GRAPHITE_H_C_INJECTION_PASS'
  call check_pah_charge_events()
  call check_pah_hydrogen()
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
  ! Independent 256-node integration of the published bulk mode densities.
  call dust_composition_curve([300d0,1000d0,3000d0],[1d0,0d0],1d0,uc,ierr)
  if(ierr/=0.or.maxval(abs(uc/[1028104816.3705101d0,10823150400.951021d0, &
       50022246461.35608d0]-1))>1d-11)stop 55
  call dust_composition_curve([300d0,1000d0,3000d0],[0d0,1d0],1d0,us,ierr)
  if(ierr/=0.or.maxval(abs(us/[1248119960.6139193d0,7551109914.322145d0, &
       27535848194.481453d0]-1))>1d-11)stop 56
  call dust_composition_curve([300d0,300.000001d0,1000d0],[.3d0,.7d0],1d0,um,ierr)
  if(ierr/=0.or.um(2)<um(1).or.abs(um(2)/um(1)-1)>1d-7)stop 57
  call dust_composition_curve([3000.001d0,20d0,100d0],[1d0,0d0],1d0,um,ierr)
  if(ierr==0)stop 58
  call dust_material_identity(.false.,material_id)
  if(any(material_id/=[1d0,dl01_t,dl01_carbon,dl01_silicate]))stop 59
  call dust_material_identity(.true.,material_id)
  if(material_id(1)/=2.or.any(.not.ieee_is_finite(material_id)))stop 60
  write(*,*)'DL01_HOT_BULK_ENERGY_CONTINUITY_DOMAIN_IDENTITY_PASS'
  dust_material_model='fixed_mix'
  dust_mass_model='bulk_v1'
  write(*,*)'DUST_MASS_CONDENSATION_GROWTH_SPUTTERING_ENERGY_PASS'
  write(*,*)'DUST_COMPOSITION_ELEMENT_BUDGET_PREMIX_SOURCE_PASS'
  write(*,*)'WSS09_ELEMENT_DEPLETION_SIGNED_COOLING_TIMESTEP_PASS'
  write(*,*)'DUST_TWO_SIZE_BUDGET_GROWTH_EROSION_COAG_SHATTER_PASS'
  write(*,*)'DUST_SN_AMBIENT_EXPOSURE_FRESH_PROTECTION_PASS'
  write(*,*)'DL01_COMPOSITION_MATERIAL_SIZE_AREA_DOMAIN_PASS'
  call check_multibin_reference()
  call check_collision_reference()
  call check_sublimation()
  call check_drag()
  call check_phase_kinetics()
  call check_multifluid()
  call check_stochastic()
  call check_iron_material()
  call check_iron_kinetics()
contains
  subroutine check_pah_hydrogen()
    use dust_pah_hydrogen
    integer,parameter::n=129
    real(real64),parameter::ev=1.602176634d-12,kb=1.380649d-16
    real(real64)::down(n,0:13),band(3,n,0:13)
    real(real64)::levels(n),diss(n,0:13),bond(0:13),attach(0:13),neutral(0:13),r(n),r2(n)
    real(real64)::old(n,0:13),next(n,0:13),saved(n,0:13),light(3),saved_light(3),h,heat,db,energy
    real(real64)::h_before,modes(102),test_rates(5),test_levels(5),balance,reference,mean_h,initial_h
    real(real64)::nh_axis(0:13),rates(0:13),ratio,previous,initial_energy
    real(real64)::carrier(3584),hc(2)
    integer::j,nh,status,refinement,steps
    dust_mass_enabled=.true.;dust_pah_model='pah_hydrogen_m13_dl01_v1'
    if(dust_pah_nstate()/=3584.or.dust_pah_charge_size()/=1792)stop 537
    carrier=0;carrier(3201)=dust_pah_state_mass(3201)
    hc=dust_pah_inventory(carrier)
    if(abs(hc(1)/(11*1.66d-24)-1)>1d-14.or.abs(hc(2)/(288*1.66d-24)-1)>1d-14)stop 538
    if(abs(dust_pah_solid_charge(carrier,1.66d-24)/1.66d-24-1)>1d-14)stop 539
    dust_mass_enabled=.false.;dust_pah_model='none'
    call pah_hydrogen_parameters(1,bond,attach,status)
    if(status/=0.or.attach(12)/=1.4d-10.or.attach(11)/=5d-11.or.attach(13)/=0)stop 501
    call pah_hydrogen_parameters(0,bond,neutral,status)
    if(status/=0.or.any(neutral/=0))stop 502
    ! Exact count: two identical harmonic modes have rho(m*q)=m+1.
    test_levels=[0d0,1d0,2d0,3d0,4d0]*ev;test_rates=-1
    call pah_harmonic_dissociation([ev,ev],ev,test_levels,ev,1d0,test_rates,status)
    if(status/=0.or.any(test_rates(:2)/=0).or.abs(test_rates(5)-.8d0)>1d-15)stop 503
    do j=1,n
       levels(j)=.5d0*(j-1)*ev
    enddo
    diss=0
    call dust_pah_modes(24,12,modes,status)
    if(status/=0)stop 504
    r=0;r2=0
    call pah_harmonic_dissociation(modes,.001d0*ev,levels,4.8d0*ev,6.8d17,r,status)
    if(status/=0)stop 505
    call pah_harmonic_dissociation(modes,.0005d0*ev,levels,4.8d0*ev,6.8d17,r2,status)
    if(status/=0)stop 506
    ratio=maxval(abs(r(21:)-r2(21:))/r2(21:))
    write(*,'(A,4ES22.14)')'PAH_DOS refinement_max_10to64eV,k10,k20,k40=',ratio,r2(21),r2(41),r2(81)
    if(ratio>.15d0.or.any(r2(:10)/=0).or.any(r2(21:)<=0))stop 507
    do nh=0,13
       nh_axis(nh)=nh
    enddo
    ! Algebraic competition with analytically known one-level first branch.
    ! Emitted band powers must match the cooling transition energy exactly.
    down=0;band=0
    diss(n,12)=3d0;down(n,12)=7d0
    band(1,n,12)=7d0*(levels(n)-levels(n-1))
    old=0;old(n,12)=1d-6;next=-1;light=-1;h=0;db=-1
    call pah_hydrogen_cool(levels,down,band,diss,bond,old,1d0,next,h,light,db,status)
    if(status/=0.or.abs(h/(3d-6/11)-1)>1d-13)stop 509
    if(abs(next(n,12)/(1d-6/11)-1)>1d-13.or.sum(light)<=0)stop 510
    balance=sum(matmul(levels,next-old))+sum(light)+db
    if(abs(balance)>1d-12*sum(matmul(levels,old)))stop 511
    if(abs(h+dot_product(nh_axis,sum(next,dim=1))-12d-6)>1d-18)stop 512
    saved=next;saved_light=light;h_before=h;reference=db;diss(1,12)=1
    call pah_hydrogen_cool(levels,down,band,diss,bond,old,1d0,next,h,light,db,status)
    if(status==0.or.any(next/=saved).or.any(light/=saved_light).or.h/=h_before.or.db/=reference)stop 513
    diss=0
    ! Finite-H donor exhausted asymptotically, not clipped after reaction.
    old=0;old(1,11)=1d-6;next=old;h=2d-8;h_before=h;heat=0;db=0
    call pah_hydrogen_attach(levels,bond,attach,old,1d20,kb*300,next,h,heat,db,status)
    if(status/=0.or.h<0.or.h>=h_before*.01d0.or.heat>=0.or.db>=0)stop 514
    balance=sum(matmul(levels,next-old))+heat+db
    if(abs(balance)>1d-12*abs(db))stop 515
    if(abs(h+dot_product(nh_axis,sum(next,dim=1))-h_before-11d-6)>1d-18)stop 516
    saved=next;h_before=h;reference=db;energy=heat
    call pah_hydrogen_attach(levels,bond,attach,old,-1d0,kb*300,next,h,heat,db,status)
    if(status==0.or.any(next/=saved).or.h/=h_before.or.heat/=energy.or.db/=reference)stop 517
    old=0;old(n,11)=1d-6;next=saved;h=1;h_before=h
    call pah_hydrogen_attach(levels,bond,attach,old,1d10,kb*300,next,h,heat,db,status)
    if(status==0.or.any(next/=saved).or.h/=h_before.or.heat/=energy.or.db/=reference)stop 518
    ! Zero neutral attachment is exactly an identity, including gas inventory.
    call pah_hydrogen_attach(levels,bond,neutral,old,1d20,kb*300,next,h,heat,db,status)
    if(status/=0.or.any(next/=old).or.h/=h_before.or.heat/=0.or.db/=0)stop 519
    ! Finite-time refinement against exact single-step bimolecular capture;
    ! chain disabled after H12 so the reference is analytically independent.
    rates=0;rates(11)=5d-11;previous=huge(1d0)
    do refinement=0,3
       steps=2**refinement;old=0;old(1,11)=1d-6;h=1d-6
       do j=1,steps
          next=old;heat=0;db=0
          call pah_hydrogen_attach(levels,bond,rates,old,1d16/steps,0d0,next,h,heat,db,status)
          if(status/=0)stop 520
          old=next
       enddo
       reference=1d-6/(1+5d-11*1d-6*1d16)
       ratio=abs(h-reference)
       write(*,'(A,I4,ES22.14)')'PAH_H_ATTACH dt_refinement,error=',steps,ratio
       if(ratio>=previous)stop 521
       previous=ratio
    enddo
    write(*,*)'PAH_H_LOSS_IR_COMPETITION_FINITE_H_REATTACH_BOND_ENERGY_ROLLBACK_PASS'
  end subroutine
  subroutine check_iron_kinetics()
    real(real64)::x(11),g(2),fe(2),out(2),scalar(2),half_fe(2),pm(3,0:2),p0(3,0:2),hc(2)
    real(real64)::rad(2),y,heat,k0,k1,rg,gas0(11),gas1(11),a,expected,dt
    integer::status,j,mode
    dust_mass_enabled=.true.;dust_mass_model='carbon_olivine_2size_v1'
    dust_cooling='chimes_neq_v1';dust_material_model='dl01_composition_v1';dust_optics_model='d03_transport_v1'
    dust_size_radius_cm=[1d-6,1d-5];dust_size_density=[2.2d0,3.8d0]
    dust_iron_model='fe_electric_compare_v1';dust_fe_kinetics=.true.;dust_fe_sticking=.3d0
    dust_sn_shocks=.false.;dust_sublimation='none';dust_growth=.true.;dust_sputtering=.true.
    if(.not.dust_mass_parameters_ok())stop 701
    x=olivine_fraction*.6d0;x(3)=.2d0;x(11)=x(11)+.1d0;x(1)=9.1d0
    g=[.02d0,.06d0];fe=[.01d0,.02d0];rad=[1d-6,1d-5];hc=[.001d0,.01d0]
    ! Independently evaluated decimal polynomial at log10(T)=6.
    call dust_fe_sputtering_yield(1d6,y,status)
    if(status/=0.or.abs(y/(10d0**(-6.8845312d0)*1d-4/31557600d0)-1)>1d-10)stop 702
    call dust_fe_sputtering_yield(9999d0,y,status)
    if(status/=0.or.y/=0)stop 703
    call dust_fe_sputtering_yield(1.001d9,y,status)
    if(status==0)stop 704
    call dust_gas_elements(x,g,gas0,status,sum(fe),hc)
    if(status/=0)stop 705
    do mode=1,2 ! cold accretion, hot erosion; same gas and solid inventory
       dt=merge(1d14,1d11,mode==1)
       p0=0;p0(:,0)=[.3d0,-.2d0,.1d0]
       p0(:,1)=[.01d0,.02d0,0d0];p0(:,2)=[-.03d0,0d0,.02d0]
       pm=p0;heat=-1
       rg=10d0-sum(g)-sum(hc)-sum(fe)
       k0=.5d0*sum(p0(:,0)**2)/rg
       do j=1,2
          k0=k0+.5d0*sum(p0(:,j)**2)/fe(j)
       enddo
       call dust_fe_kinetics_step(10d0,x,g,fe,merge(100d0,1d6,mode==1),1d-24,dt,rad,7.87d0, &
            out,status,hc,pm,heat)
       if(status/=0.or.heat<0.or.maxval(abs(sum(pm,dim=2)-sum(p0,dim=2)))>1d-13)then
          write(*,*)'Fe kinetic failure mode/status/heat=',mode,status,heat
          stop 706
       endif
       if(mode==1.and.any(out<=fe))stop 707
       if(mode==2.and.any(out>=fe))stop 708
       call dust_fe_kinetics_step(10d0,x,g,fe,merge(100d0,1d6,mode==1),1d-24,dt,rad,7.87d0,scalar,status,hc)
       if(status/=0.or.any(out/=scalar))stop 709
       call dust_gas_elements(x,g,gas1,status,sum(out),hc)
       if(status/=0.or.abs(sum(gas1)+sum(g)+sum(out)+sum(hc)-sum(x))>1d-13)stop 710
       if(any(gas1(1:10)/=gas0(1:10)))stop 711
       rg=10d0-sum(g)-sum(hc)-sum(out)
       k1=.5d0*sum(pm(:,0)**2)/rg
       do j=1,2
          k1=k1+.5d0*sum(pm(:,j)**2)/out(j)
       enddo
       if(abs(k1+heat-k0)>1d-12*k0)stop 712
    enddo
    ! One occupied bin has a closed-form logistic solution, no seed creation.
    fe=[.01d0,0d0];dt=1d15
    a=.75d0*.3d0*(x(11)-olivine_fraction(11)*g(2))*1d-24* &
         sqrt(8*dust_kb*100/(acos(-1d0)*56*dust_mp))/(7.87d0*rad(1))
    expected=(x(11)-olivine_fraction(11)*g(2))/ &
         (1+((x(11)-olivine_fraction(11)*g(2))/fe(1)-1)*exp(-a*dt))
    call dust_fe_kinetics_step(10d0,x,g,fe,100d0,1d-24,dt,rad,7.87d0,out,status,hc)
    if(status/=0.or.out(2)/=0.or.abs(out(1)/expected-1)>1d-12)stop 713
    fe=[.01d0,.02d0]
    call dust_fe_kinetics_step(10d0,x,g,fe,100d0,1d-24,1d16,rad,7.87d0,out,status,hc)
    if(status/=0)stop 714
    call dust_fe_kinetics_step(10d0,x,g,fe,100d0,1d-24,5d15,rad,7.87d0,half_fe,status,hc)
    if(status/=0)stop 715
    call dust_fe_kinetics_step(10d0,x,g,half_fe,100d0,1d-24,5d15,rad,7.87d0,scalar,status,hc)
    if(status/=0.or.maxval(abs(out-scalar))>1d-4)stop 716
    ! Rate flags off and inactive selector are exact identities. Failure is atomic.
    dust_growth=.false.;dust_sputtering=.false.;pm=p0;heat=-1
    call dust_fe_kinetics_step(10d0,x,g,fe,100d0,1d-24,1d15,rad,7.87d0,out,status,hc,pm,heat)
    if(status/=0.or.any(out/=fe).or.any(pm/=p0).or.heat/=0)stop 717
    dust_sputtering=.true.;pm=p0;heat=-7
    call dust_fe_kinetics_step(10d0,x,g,fe,2d9,1d-24,1d15,rad,7.87d0,out,status,hc,pm,heat)
    if(status==0.or.any(out/=fe).or.any(pm/=p0).or.heat/=-7)stop 718
    dust_sn_shocks=.true.
    if(dust_mass_parameters_ok())stop 719
    dust_sn_shocks=.false.;dust_fe_kinetics=.false.;dust_fe_sticking=0;dust_growth=.true.
    call dust_fe_kinetics_step(10d0,x,g,fe,100d0,1d-24,1d15,rad,7.87d0,out,status,hc)
    if(status/=0.or.any(out/=fe))stop 720
    write(*,*)'FE_KINETICS_SOURCE_UNITS_SEED_GROWTH_THERMAL_EROSION_MOMENTUM_ENERGY_ROLLBACK_PASS'
  end subroutine

  subroutine check_phase_kinetics()
    ! Focused source kinetics checks, independent of check_multifluid.
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
    real(real64)::x(11),reserved(11),bins0(4),out(4),scalar(4),pm(3,0:4),p0(3,0:4)
    real(real64)::heat,k0,k1,r0(0:4),r1(0:4),aa,bb,dt,stick0,density0(2),updated,flows(2),sd
    real(real64)::pair_m(2),pair_n(2),pair_p(3,2),pair_old(3,2),pair_q,v0(3),v1(3),fixed,hc(2)
    real(real64)::expected(3),empty(4),failed(4),saved_heat,boost(3),shift_p(3,0:4),shift_heat
    logical::grow0,sputter0,coag0,shatter0
    integer::status,j,mode
    grow0=dust_growth;sputter0=dust_sputtering;coag0=dust_coagulation;shatter0=dust_shattering
    stick0=dust_sticking;density0=dust_size_density
    dust_coagulation=.false.;dust_shattering=.false.
    dust_size_density=3d0

    ! Analytic equilibrium: each gross flow is b*D*dt despite zero net mass.
    call dust_mass_step(1d0,.5d0,.025d0,2d0,1d0,updated,status)
    if(status/=0.or.abs(updated-.5d0)>1d-15)stop 701
    call dust_phase_gross(1d0,.5d0,.025d0,2d0,1d0,updated,flows,status)
    if(status/=0.or.maxval(abs(flows-.0125d0))>1d-15)stop 702
    pair_m=[1d0,.5d0];pair_n=pair_m
    pair_old(:,1)=[3d0,-1d0,2d0];pair_old(:,2)=.5d0*[-2d0,1d0,-3d0]
    pair_p=pair_old;pair_q=-1
    call dust_phase_pair(pair_m,pair_n,flows,pair_p,pair_q,status)
    if(status/=0.or.pair_q<=0.or.all(pair_p==pair_old))stop 703
    ! Independently verify both BE equations, not just summed momentum.
    v0=pair_p(:,1)/pair_n(1);v1=pair_p(:,2)/pair_n(2)
    if(maxval(abs(pair_p(:,1)-pair_old(:,1)+flows(1)*v0-flows(2)*v1))>2d-14)stop 704
    if(maxval(abs(pair_p(:,2)-pair_old(:,2)-flows(1)*v0+flows(2)*v1))>2d-14)stop 705
    k0=.5d0*(sum(pair_old(:,1)**2)/pair_m(1)+sum(pair_old(:,2)**2)/pair_m(2))
    k1=.5d0*(sum(pair_p(:,1)**2)/pair_n(1)+sum(pair_p(:,2)**2)/pair_n(2))
    if(abs(k1+pair_q-k0)>2d-14)stop 706
    ! Tiny valid accretion used to give negative heat from cancellation of
    ! nearly equal reduced kinetic energies (m1=.1, m2=.01, flow=1e-18).
    pair_p(:,1)=[.3d0,-.1d0,.2d0];pair_p(:,2)=[-.02d0,.01d0,-.03d0]
    call dust_phase_pair([.1d0,.01d0],[.1d0-1d-18,.01d0+1d-18],[1d-18,0d0],pair_p,pair_q,status)
    if(status/=0.or.pair_q<=0.or.abs(pair_q/(.5d-18*54d0)-1)>2d-14)stop 741
    ! Repeated coupled pair updates converge to the exact constant-mass ODE.
    pair_p=pair_old
    do j=1,100
       call dust_phase_pair(pair_m,pair_m,flows/100,pair_p,pair_q,status)
       if(status/=0)stop 707
    enddo
    expected=(pair_old(:,1)/pair_m(1)-pair_old(:,2)/pair_m(2))* &
         exp(-.0125d0*(1/pair_m(1)+1/pair_m(2)))
    if(maxval(abs(pair_p(:,1)/pair_m(1)-pair_p(:,2)/pair_m(2)-expected))>4d-5)stop 708

    ! Actual size-step equilibrium with physical rates tuned to a=2b.
    dust_growth=.true.;dust_sputtering=.true.;dust_sticking=1
    call dust_mass_rates(10d-24,.2d-24,100d0,aa,bb,status)
    if(status/=0.or.min(aa,bb)<=0)stop 709
    dust_sticking=2*bb/aa;dt=.2d0/(bb*dust_grain_radius_cm/dust_size_radius_cm(1))
    x=0;x(1)=9.8d0;x(3)=.2d0;bins0=[.1d0,0d0,0d0,0d0]
    r0=[10-sum(bins0),bins0];p0=0
    p0(:,0)=r0(0)*[2d0,-1d0,3d0];p0(:,1)=bins0(1)*[-3d0,2d0,-1d0]
    pm=p0;heat=-1
    call dust_size_step(10d0,.2d0,x,bins0,100d0,1d-24,dt,scalar,status)
    if(status/=0)stop 710
    call dust_size_step(10d0,.2d0,x,bins0,100d0,1d-24,dt,out,status,pm,heat)
    if(status/=0.or.any(out/=scalar).or.maxval(abs(out-bins0))>2d-15)stop 711
    if(heat<=0.or.all(pm==p0))stop 712
    r1=[10-sum(out),out]
    call phase_kinetics_closure(r0,r1,p0,pm,heat)
    ! Galilean shift: heat and relative evolution are independent of bulk flow.
    boost=[11d0,-7d0,5d0];shift_p=p0
    do j=0,4
       shift_p(:,j)=shift_p(:,j)+r0(j)*boost
    enddo
    call dust_size_step(10d0,.2d0,x,bins0,100d0,1d-24,dt,out,status,shift_p,shift_heat)
    if(status/=0.or.abs(shift_heat-heat)>2d-12)stop 713
    do j=0,4
       if(maxval(abs(shift_p(:,j)-pm(:,j)-r1(j)*boost))>2d-12)stop 714
    enddo

    ! One-way accretion with fixed Fe AND PAH excluded from gas inertia.
    dust_sputtering=.false.;dust_sticking=.3d0
    hc=[.01d0,.02d0];fixed=.04d0+sum(hc)
    x=0;x(1)=1.74d0;x(3)=.22d0;x(11)=.04d0
    reserved=x;reserved(1)=x(1)-hc(1);reserved(3)=x(3)-hc(2);reserved(11)=0
    r0=[2-fixed-sum(bins0),bins0];p0=0
    p0(:,0)=r0(0)*[2d0,-1d0,3d0];p0(:,1)=bins0(1)*[-3d0,2d0,-1d0]
    call dust_mass_rates(2d-24,.2d-24,100d0,aa,bb,status)
    dt=.08d0/(aa*dust_grain_radius_cm/dust_size_radius_cm(1))
    call dust_size_step_reserved_iron(2d0,.26d0,x,bins0,.04d0,100d0,1d-24,dt,scalar,status,hc)
    if(status/=0)stop 715
    pm=p0;heat=-1
    call dust_size_step_reserved_iron(2d0,.26d0,x,bins0,.04d0,100d0,1d-24,dt,out,status,hc,pm,heat)
    if(status/=0.or.any(out/=scalar).or.out(1)<=bins0(1))stop 716
    r1=[2-fixed-sum(out),out]
    call phase_kinetics_closure(r0,r1,p0,pm,heat)
    if(maxval(abs(pm(:,0)/r1(0)-p0(:,0)/r0(0)))>2d-14)stop 717
    expected=p0(:,1)+(out(1)-bins0(1))*p0(:,0)/r0(0)
    if(maxval(abs(pm(:,1)-expected))>2d-14)stop 718
    shift_p=p0;shift_heat=-1
    call dust_size_step(2d0,.2d0,reserved,bins0,100d0,1d-24,dt,scalar,status,shift_p,shift_heat,fixed)
    if(status/=0.or.any(scalar/=out).or.any(shift_p/=pm).or.shift_heat/=heat)stop 719

    ! One-way sputtering preserves the current grain donor velocity.
    dust_growth=.false.;dust_sputtering=.true.
    call dust_mass_rates(2d-24,.2d-24,2d6,aa,bb,status)
    dt=.2d0/(bb*dust_grain_radius_cm/dust_size_radius_cm(1));pm=p0
    call dust_size_step_reserved_iron(2d0,.26d0,x,bins0,.04d0,2d6,1d-24,dt,out,status,hc,pm,heat)
    if(status/=0.or.out(1)>=bins0(1))stop 720
    r1=[2-fixed-sum(out),out]
    call phase_kinetics_closure(r0,r1,p0,pm,heat)
    if(maxval(abs(pm(:,1)/out(1)-p0(:,1)/bins0(1)))>2d-14)stop 721

    ! Both directions of bin exchange, including initially empty receiver.
    dust_growth=.false.;dust_sputtering=.false.
    do mode=1,2
       dust_shattering=mode==1;dust_coagulation=mode==2
       x=0;x(3)=.2d0;x(1)=9.8d0
       sd=merge(1d0,2d3,mode==1)*dust_mp/x(1)
       bins0=[.04d0,.06d0,0d0,0d0]
       if(mode==2)bins0=[.1d0,0d0,0d0,0d0]
       r0=[10-sum(bins0),bins0];p0=0
       p0(:,0)=r0(0)*[2d0,-1d0,3d0]
       p0(:,1)=bins0(1)*[-3d0,2d0,-1d0];p0(:,2)=bins0(2)*[4d0,-2d0,1d0]
       pm=p0
       call dust_size_step(10d0,.2d0,x,bins0,100d0,sd,1d14,scalar,status)
       if(status/=0)stop 722
       call dust_size_step(10d0,.2d0,x,bins0,100d0,sd,1d14,out,status,pm,heat)
       if(status/=0.or.any(out/=scalar).or.all(out==bins0).or.any(pm(:,0)/=p0(:,0)))stop 723
       r1=[10-sum(out),out]
       call phase_kinetics_closure(r0,r1,p0,pm,heat)
       j=merge(2,1,mode==1)
       if(maxval(abs(pm(:,j)/out(j)-p0(:,j)/bins0(j)))>2d-14)stop 724
    enddo

    ! Combined processes exercise both compositions and the shared gas state.
    dust_growth=.true.;dust_sputtering=.true.;dust_coagulation=.true.;dust_shattering=.true.
    x=.06d0*olivine_fraction;x(1)=9.9d0;x(3)=.04d0
    bins0=[.005d0,.01d0,.01d0,.02d0];r0=[10-sum(bins0),bins0]
    do j=0,4
       p0(:,j)=r0(j)*[real(j,real64)-2,real(j*j,real64),1d0-real(j,real64)]
    enddo
    pm=p0
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,1d13,scalar,status)
    if(status/=0)stop 725
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,1d13,out,status,pm,heat)
    if(status/=0.or.any(out/=scalar))stop 726
    r1=[10-sum(out),out];call phase_kinetics_closure(r0,r1,p0,pm,heat)

    ! Failure must not publish any partial momentum, heat, or mass update.
    pm=p0;saved_heat=73;heat=saved_heat
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,0d0,failed,status,pm,heat)
    if(status/=0.or.any(pm/=p0).or.heat/=0.or.any(failed/=bins0))stop 737
    heat=saved_heat
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,-1d0,failed,status,pm,heat)
    if(status==0.or.any(pm/=p0).or.heat/=saved_heat.or.any(failed/=bins0))stop 727
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,0d0,failed,status,pm)
    if(status==0.or.any(pm/=p0).or.any(failed/=bins0))stop 728
    call dust_size_step(10d0,.1d0,x,bins0,100d0,1d-22,0d0,failed,status,mixing_heat=heat)
    if(status==0.or.heat/=saved_heat)stop 729
    empty=0;pm=0;pm(1,1)=1;shift_p=pm
    call dust_size_step(10d0,.1d0,x,empty,100d0,1d-22,0d0,failed,status,pm,heat)
    if(status==0.or.any(pm/=shift_p).or.heat/=saved_heat.or.any(failed/=empty))stop 730
    pm=p0
    call dust_size_step_reserved_iron(10d0,.1d0,x,bins0,.01d0,100d0,1d-22,1d13,failed,status, &
         [.01d0,1d0],pm,heat)
    if(status==0.or.any(pm/=p0).or.heat/=saved_heat.or.any(failed/=bins0))stop 731
    ! Gas exhaustion AFTER one valid accretion substep must roll it back.
    dust_sputtering=.false.;dust_coagulation=.false.;dust_shattering=.false.
    x=0;x(1)=9.8d0;x(3)=.2d0;bins0=[.1d0,0d0,0d0,0d0]
    call dust_mass_rates(10d-24,.2d-24,100d0,aa,bb,status)
    aa=aa*dust_grain_radius_cm/dust_size_radius_cm(1);dt=.05d0/aa
    call dust_mass_step(.2d0,.1d0,dt/2,aa,0d0,updated,status)
    if(status/=0.or.updated<=.1d0)stop 738
    fixed=10-sum(bins0)-1.5d0*(updated-.1d0)
    r0=[10-fixed-sum(bins0),bins0];p0=0
    p0(:,0)=r0(0)*[2d0,-1d0,3d0];p0(:,1)=bins0(1)*[-3d0,2d0,-1d0]
    pm=p0
    call dust_size_step(10d0,.2d0,x,bins0,100d0,1d-24,dt,failed,status,pm,heat,fixed)
    if(status==0.or.any(pm/=p0).or.heat/=saved_heat.or.any(failed/=bins0))stop 739
    ! Nonfinite pair input is rejected without touching the caller's values.
    pair_m=[1d0,.1d0];pair_n=[.9d0,.2d0];flows=[.1d0,0d0]
    pair_p=pair_old;pair_p(1,2)=ieee_value(0d0,ieee_quiet_nan);pair_q=saved_heat
    call dust_phase_pair(pair_m,pair_n,flows,pair_p,pair_q,status)
    if(status==0.or.pair_q/=saved_heat.or..not.ieee_is_nan(pair_p(1,2)))stop 732
    if(any(pair_p(:,1)/=pair_old(:,1)).or.any(pair_p(2:3,2)/=pair_old(2:3,2)))stop 740
    call check_phase_kinetics_identity()
    dust_growth=grow0;dust_sputtering=sputter0;dust_coagulation=coag0;dust_shattering=shatter0
    dust_sticking=stick0;dust_size_density=density0
    write(*,*)'DUST_PHASE_KINETICS_GROSS_BE_MOMENTUM_HEAT_PARITY_ATOMIC_PASS'
  end subroutine

  subroutine check_phase_kinetics_identity()
    ! A live aggregate and the sum of its bins may differ by roundoff.
    ! The admitted element budget is NOT a physical mass transfer at zero rate.
    real(real64)::rho,metal,aggregate,bins(4),iron,hc(2),x(11),reserved(11),gas(11)
    real(real64)::pm(3,0:4),p0(3,0:4),out(4),scalar(4),heat,dt,capacity,updated,flows(2)
    real(real64)::r0(0:4),r1(0:4),fixed,aa,bb,expected
    integer::status,j,mode
    rho=9.99520328050314d0;metal=.100951941211516d0
    aggregate=4.997590777300154d-2;bins=[0d0,0d0,0d0,4.997590777300155d-2]
    iron=5.997109458852153d-2-aggregate;hc=[.001d0,.024d0]
    fixed=iron+sum(hc)
    x=olivine_fraction*aggregate;x(1)=rho-metal;x(3)=hc(2);x(11)=x(11)+iron
    reserved=x;reserved(1)=x(1)-hc(1);reserved(3)=0;reserved(11)=x(11)-iron
    call dust_gas_elements(x,[0d0,bins(4)],gas,status,iron,hc)
    if(status/=0)stop 742
    capacity=huge(1d0)
    do j=1,11
       if(olivine_fraction(j)>0)capacity=min(capacity,reserved(j)/olivine_fraction(j))
    enddo
    if(capacity>=bins(4))stop 743
    ! Reproduce the previous mismatch without weakening the strict helper:
    ! scalar min() clips, but that is not a directed kinetic process.
    call dust_mass_step(capacity,min(bins(4),capacity),1d12,0d0,0d0,updated,status)
    if(status/=0.or.updated/=capacity)stop 744
    call dust_phase_gross(capacity,bins(4),1d12,0d0,0d0,updated,flows,status)
    if(status==0)stop 745
    r0=[rho-fixed-sum(bins),bins];p0=0
    p0(:,0)=r0(0)*[2d0,-1d0,3d0];p0(:,4)=bins(4)*[-3d0,2d0,-1d0]
    dust_coagulation=.false.;dust_shattering=.false.
    do mode=1,3
       dust_growth=mode/=1;dust_sputtering=mode==3
       dt=1d12;if(mode==3)dt=0
       pm=p0;heat=73
       ! Mode 2 has growth enabled but zero physical rate at this hot T.
       call dust_size_step_reserved_iron(rho,metal,x,bins,iron,10410.7028398063d0,1d-24, &
            dt,out,status,hc,pm,heat)
       if(status/=0.or.any(out/=bins).or.any(pm/=p0).or.heat/=0)stop 746
       call phase_kinetics_closure(r0,r0,p0,pm,heat)
       ! Exercise the unwrapped API and fixed Fe/PAH gas-inertia exclusion.
       pm=p0;heat=73
       call dust_size_step(rho,metal-iron-hc(2),reserved,bins,10410.7028398063d0,1d-24, &
            dt,out,status,pm,heat,fixed)
       if(status/=0.or.any(out/=bins).or.any(pm/=p0).or.heat/=0)stop 747
       ! The existing no-momentum clipping arithmetic is intentionally unchanged.
       call dust_size_step_reserved_iron(rho,metal,x,bins,iron,10410.7028398063d0,1d-24, &
            dt,scalar,status,hc)
       if(status/=0.or.scalar(4)/=capacity.or.any(scalar(1:3)/=0))stop 748
    enddo
    ! Nonzero sputtering must erode the actual donor, not a capacity-clipped
    ! surrogate, at the same already-admitted fully depleted boundary.
    dust_growth=.false.;dust_sputtering=.true.
    call dust_mass_rates(rho*1d-24,capacity*1d-24,2d6,aa,bb,status)
    if(status/=0.or.bb<=0)stop 751
    dt=.05d0/(bb*dust_grain_radius_cm/dust_size_radius_cm(1))
    call dust_mass_step(bins(4),bins(4),dt,0d0, &
         bb*dust_grain_radius_cm/dust_size_radius_cm(2),expected,status)
    if(status/=0.or.expected>=bins(4))stop 752
    pm=p0;heat=73
    call dust_size_step_reserved_iron(rho,metal,x,bins,iron,2d6,1d-24,dt,out,status,hc,pm,heat)
    if(status/=0.or.out(4)/=expected.or.any(out(1:3)/=0).or.heat<=0)stop 753
    r1=[rho-fixed-sum(out),out];call phase_kinetics_closure(r0,r1,p0,pm,heat)
    if(maxval(abs(pm(:,4)/out(4)-p0(:,4)/bins(4)))>2d-14)stop 754
    if(maxval(abs(pm(:,0)-p0(:,0)-(bins(4)-out(4))*p0(:,4)/bins(4)))>2d-14)stop 755
    ! The same boundary also permits simultaneous active growth/destruction.
    dust_growth=.true.
    pm=p0;heat=73
    call dust_size_step_reserved_iron(rho,metal,x,bins,iron,100d0,1d-24,dt,out,status,hc,pm,heat)
    if(status/=0.or.heat<0)stop 756
    r1=[rho-fixed-sum(out),out];call phase_kinetics_closure(r0,r1,p0,pm,heat)
    call dust_gas_elements(x,[sum(out(1:2)),sum(out(3:4))],gas,status,iron,hc)
    if(status/=0)stop 757
    write(*,*)'DUST_PHASE_KINETICS_ACTIVE_SATURATION_ACTUAL_DONOR_PASS'
    ! Zero gas/bin rates must not skip an active bin/bin exchange substep.
    dust_growth=.false.;dust_sputtering=.false.;dust_shattering=.true.
    pm=p0;heat=73
    call dust_size_step_reserved_iron(rho,metal,x,bins,iron,10410.7028398063d0,1d-24, &
         1d12,out,status,hc,pm,heat)
    if(status/=0.or.out(3)<=0.or.out(4)>=bins(4).or.any(pm(:,0)/=p0(:,0)))stop 749
    r1=[rho-fixed-sum(out),out];call phase_kinetics_closure(r0,r1,p0,pm,heat)
    ! Still reject a genuine budget violation at zero rate, atomically.
    dust_shattering=.false.;x(5)=x(5)*.99d0;pm=p0;heat=73
    call dust_size_step_reserved_iron(rho,metal,x,bins,iron,10410.7028398063d0,1d-24, &
         1d12,out,status,hc,pm,heat)
    if(status==0.or.any(out/=bins).or.any(pm/=p0).or.heat/=73)stop 750
    dust_sputtering=.true.
    call dust_size_step_reserved_iron(rho,metal,x,bins,iron,2d6,1d-24,dt,out,status,hc,pm,heat)
    if(status==0.or.any(out/=bins).or.any(pm/=p0).or.heat/=73)stop 758
    write(*,*)'DUST_PHASE_KINETICS_ZERO_RATE_ROUNDOFF_IDENTITY_PASS'
  end subroutine

  subroutine phase_kinetics_closure(r0,r1,p0,p1,heat)
    real(real64),intent(in)::r0(0:4),r1(0:4),p0(3,0:4),p1(3,0:4),heat
    real(real64)::before,after
    integer::j
    before=0;after=0
    do j=0,4
       if(r0(j)>0)before=before+.5d0*sum(p0(:,j)**2)/r0(j)
       if(r1(j)>0)after=after+.5d0*sum(p1(:,j)**2)/r1(j)
       if(r1(j)==0.and.any(p1(:,j)/=0))stop 733
    enddo
    if(heat<0.or.abs(after+heat-before)>2d-12*max(before,1d0))stop 734
    if(maxval(abs(sum(p0,dim=2)-sum(p1,dim=2)))>2d-13*max(maxval(abs(p0)),1d0))stop 735
    if(abs(sum(r0)-sum(r1))>2d-14*sum(r0))stop 736
  end subroutine

  subroutine check_multifluid()
    ! Seven dust fluids: C/S two-size, Fe two-size and one PAH population.
    ! CR + a gas nuclear carrier + four grain energies + 128 PAH states.
    ! No RAMSES launch, new diagnostic framework or configurable live mode.
    integer,parameter::nb=7,ncarrier=134,nv=5+4*nb+ncarrier,base=5+4*nb
    type(dust_fv_layout)::layout
    integer::owner(ncarrier),status,i,j,k,b,n,r,step
    integer,allocatable::links(:,:)
    real(real64)::gnt(ncarrier),u(nv),v(nv),flux(nv),speed,face_v,heat
    real(real64)::rho(0:nb),vel(3,0:nb),th,kin,pressure,cs,dt,dx,t,x,pi,err(3),exact,delta
    real(real64)::sum0(nv),scale(nv),relative,ts(nb),gamma_test
    real(real64)::impulse(3,nb),debit(nb),solid_heat(nb),photon_count(2),moment(3,2),fractions(nb,2)
    real(real64)::before_gas(3),after_gas(3),photon_energy(2),expected_impulse(3)
    real(real64),allocatable::state(:,:),next(:,:),initial(:,:),times(:,:),drag_heat(:)
    gamma_test=5d0/3;pi=acos(-1d0)
    owner=7;owner(1:2)=0;owner(3:6)=[1,2,3,4];gnt=0;gnt(1)=4d0/3
    call dust_fv_initialize(nb,owner,gnt,layout,status)
    if(status/=0)stop 401
    call make_multifluid_state(.25d0,u)
    call dust_fv_decode(layout,u,gamma_test,rho,vel,th,kin,pressure,cs,status)
    if(status/=0.or.abs(th-1.5d0)>1d-14.or.abs(pressure-1.1d0)>1d-14)stop 402
    if(abs(rho(0)-1)>1d-14.or.abs(vel(1,0)-.2d0)>1d-14)stop 403
    before_gas=rho(0)*vel(:,0)
    impulse=0;impulse(1,1)=.002d0;impulse(2,7)=-.001d0
    debit=0;debit(1)=.01d0;debit(7)=.02d0;solid_heat=-42;v=u
    call dust_fv_radiation_kick(layout,u,gamma_test,impulse,debit,v,solid_heat,status)
    if(status/=0.or.any(solid_heat<0))stop 426
    if(maxval(abs(v(2:4)-u(2:4)-sum(impulse,dim=2)))>1d-15)stop 427
    if(abs(v(5)-u(5)+sum(solid_heat)-sum(debit))>1d-14)stop 428
    call dust_fv_decode(layout,v,gamma_test,rho,vel,th,kin,pressure,cs,status)
    after_gas=rho(0)*vel(:,0)
    if(status/=0.or.abs(th-1.5d0)>1d-14.or.maxval(abs(after_gas-before_gas))>1d-14)stop 429
    ! Do NOT add the whole absorbed energy to E_hydro and then grain heat.
    ! Underfunded work rejects the entire kick, including an earlier grain.
    impulse(1,7)=1d0;solid_heat=-42;v=u
    call dust_fv_radiation_kick(layout,u,gamma_test,impulse,debit,v,solid_heat,status)
    if(status==0.or.any(v/=u).or.any(solid_heat/=-42))stop 430
    photon_count=[2d0,3d0];moment=0;moment(:,1)=[2d0,0d0,0d0];moment(:,2)=[0d0,-1d0,0d0]
    fractions=0;fractions(1,:)=[.25d0,.75d0];fractions(7,:)=1-fractions(1,:)
    photon_energy=[1d-12,2d-12];solid_heat=-42;v=u
    call dust_fv_absorption_kick(layout,u,gamma_test,photon_count,moment,photon_energy,fractions, &
         1d6,1d-3,1d-13,v,solid_heat,status)
    expected_impulse=matmul(moment,photon_energy)*1d6/(2.99792458d10*1d-13)
    if(status/=0.or.maxval(abs(v(2:4)-u(2:4)-expected_impulse))>1d-15)stop 431
    if(abs(v(5)-u(5)+sum(solid_heat)-sum(photon_count*photon_energy)*1d9)>1d-14)stop 432
    ! Same absorbed photons in an isotropic field: heating but no force.
    moment=0;v=u;solid_heat=-42
    call dust_fv_absorption_kick(layout,u,gamma_test,photon_count,moment,photon_energy,fractions, &
         1d6,1d-3,1d-13,v,solid_heat,status)
    if(status/=0.or.any(v/=u).or.abs(sum(solid_heat)-.008d0)>1d-14)stop 433
    fractions(7,2)=.9d0;v=u;solid_heat=-42
    call dust_fv_absorption_kick(layout,u,gamma_test,photon_count,moment,photon_energy,fractions, &
         1d6,1d-3,1d-13,v,solid_heat,status)
    if(status==0.or.any(v/=u).or.any(solid_heat/=-42))stop 434
    fractions(7,2)=.25d0;moment(1,1)=2.01d0
    call dust_fv_absorption_kick(layout,u,gamma_test,photon_count,moment,photon_energy,fractions, &
         1d6,1d-3,1d-13,v,solid_heat,status)
    if(status==0.or.any(v/=u).or.any(solid_heat/=-42))stop 435
    write(*,*)'DUST_FV_ABSORPTION_PHYSICAL_C_GAS_UNCHANGED_WORK_HEAT_ROLLBACK_PASS'
    flux=0;speed=0;face_v=0
    call dust_fv_face(layout,u,u,1,gamma_test,flux,speed,face_v,status)
    if(status/=0.or.abs(flux(6)/u(6)+.4d0)>1d-14.or.abs(flux(10)/u(10)-.7d0)>1d-14)stop 404
    if(abs(sum(flux(base+7:))-flux(30))>1d-14)stop 405
    ! Zero net barycentric mass flux must NOT erase finite counterflows.
    v=u;v(2)=0;v(5)=u(5)-.5d0*.2d0**2+.5d0*(.2d0-u(2))**2
    call dust_fv_face(layout,v,v,1,gamma_test,flux,speed,face_v,status)
    if(status/=0.or.flux(1)/=0.or.flux(6)==0.or.flux(base+2)==0)stop 422
    ! Stiff drag is stable and energy-funded, with frozen local stopping times.
    ! This alone is NOT an asymptotic-preserving dustywave accuracy claim.
    ts=1d0;v=u;heat=-1
    call dust_fv_drag(layout,u,gamma_test,ts,1d200,v,heat,status)
    if(status/=0.or.heat<=0.or.any(v(1:5)/=u(1:5)).or.any(v(base+1:)/=u(base+1:)))stop 406
    call dust_fv_decode(layout,v,gamma_test,rho,vel,th,kin,pressure,cs,status)
    if(status/=0.or.maxval(abs(vel-spread(v(2:4)/v(1),2,nb+1)))>1d-13)stop 407
    if(abs(th-1.5d0-heat)>1d-13)stop 408
    v=u;v(6)=0;flux=42;speed=42;face_v=42
    call dust_fv_face(layout,u,v,1,gamma_test,flux,speed,face_v,status)
    if(status==0.or.any(flux/=42).or.speed/=42.or.face_v/=42)stop 409
    ! Passive owner validation must not corrupt an already valid layout.
    owner(1)=7
    call dust_fv_initialize(nb,owner,gnt,layout,status)
    if(status==0)stop 410
    owner(1)=0
    do r=1,3
       n=16*2**r;dx=1d0/n;dt=.02d0*dx;t=0
       allocate(state(nv,n),next(nv,n),links(6,n))
       do i=1,n
          call make_multifluid_state((i-.5d0)*dx,state(:,i))
          links(:,i)=i;links(1,i)=modulo(i-2,n)+1;links(2,i)=modulo(i,n)+1
       enddo
       sum0=sum(state,dim=2);scale=max(sum(abs(state),dim=2),tiny(1d0));initial=state
       next=state
       do step=1,5*n/2
          call dust_fv_periodic_step(layout,state,links,dx,dt,gamma_test,next,status)
          if(status/=0)then
             write(*,*)'DUST_FV_FAILED',r,step,status
             stop 411
          endif
          state=next;t=t+dt
       enddo
       relative=maxval(abs(sum(state,dim=2)-sum0)/scale)
       if(relative>2d-12)stop 412
       err(r)=0
       do i=1,n
          x=(i-.5d0)*dx
          exact=.01d0*(1+.2d0*sin(2*pi*(x+.4d0*t)))
          err(r)=err(r)+abs(state(6,i)-exact)/n
          if(abs(sum(state(base+7:,i))-state(30,i))>1d-14)stop 413
          call dust_fv_decode(layout,state(:,i),gamma_test,rho,vel,th,kin,pressure,cs,status)
          if(status/=0.or.maxval(abs(vel(1,:)-[.2d0,-.4d0,.7d0,.1d0,.1d0,.1d0,.1d0,.8d0]))>1d-11)stop 414
       enddo
       write(*,'(A,I5,2ES18.9)')'DUST_FV_COUNTERFLOW_CELLS_L1_CONSERVATION=',n,err(r),relative
       ! A CFL rejection must preserve all positive and signed fields.
       next=initial
       call dust_fv_periodic_step(layout,state,links,dx,1d0,gamma_test,next,status)
       if(status/=dust_fv_cfl.or.any(next/=initial))stop 415
       deallocate(state,next,links,initial)
    enddo
    if(any(err(1:2)/err(2:3)<1.7d0))stop 416
    ! Pressure gradient accelerates gas, not pressureless dust; CR pressure
    ! is included in total pressure and sound speed, not assigned to grains.
    n=8;dx=1d0/n;dt=1d-3
    allocate(state(nv,n),next(nv,n),links(6,n))
    do i=1,n
       call make_multifluid_state((i-.5d0)*dx,state(:,i))
       state(2:4,i)=0
       do b=1,nb
          k=6+4*(b-1);state(k+1:k+3,i)=0
       enddo
       state(5,i)=1.5d0*(1+.1d0*sin(2*pi*(i-.5d0)*dx))+.3d0
       links(:,i)=i;links(1,i)=modulo(i-2,n)+1;links(2,i)=modulo(i,n)+1
    enddo
    next=state
    call dust_fv_periodic_step(layout,state,links,dx,dt,gamma_test,next,status)
    if(status/=0.or.maxval(abs(next(2,:)))<=0)stop 417
    do b=1,nb
       k=6+4*(b-1)
       if(any(next(k+1:k+3,:)/=0))stop 418
    enddo
    ! Uniform CR energy with nonuniform gas velocity: advection plus PdV
    ! changes it by -gamma_cr*e_cr*dt*div(v_g), leaving total E conservative.
    do i=1,n
       call make_multifluid_state((i-.5d0)*dx,state(:,i))
       call dust_fv_decode(layout,state(:,i),gamma_test,rho,vel,th,kin,pressure,cs,status)
       state(2,i)=state(2,i)-.2d0+.1d0*sin(2*pi*(i-.5d0)*dx)
       state(5,i)=state(5,i)-.5d0*.2d0**2+.5d0*(.1d0*sin(2*pi*(i-.5d0)*dx))**2
    enddo
    next=state
    call dust_fv_periodic_step(layout,state,links,dx,dt,gamma_test,next,status)
    if(status/=0)stop 419
    do i=1,n
       delta=.1d0*(sin(2*pi*(links(2,i)-.5d0)*dx)-sin(2*pi*(links(1,i)-.5d0)*dx))/(2*dx)
       if(abs(next(base+1,i)-.3d0+dt*(4d0/3)*.3d0*delta)>1d-14)stop 420
    enddo
    if(abs(sum(next(5,:))-sum(state(5,:)))>1d-13)stop 421
    allocate(times(nb,n),drag_heat(n));times=.1d0
    next=state;drag_heat=0
    call dust_fv_periodic_advance(layout,state,links,dx,dt,gamma_test,times,next,drag_heat,status)
    if(status/=0.or.minval(drag_heat)<=0)stop 423
    if(maxval(abs(sum(next(1:5,:),dim=2)-sum(state(1:5,:),dim=2)))>1d-13)stop 424
    initial=next;times(:,n)=-1;drag_heat=-42
    call dust_fv_periodic_advance(layout,state,links,dx,dt,gamma_test,times,next,drag_heat,status)
    if(status==0.or.any(next/=initial).or.any(drag_heat/=-42))stop 425
    write(*,*)'DUST_FV_COUNTERFLOW_ZERO_NET_FLUX_SPATIAL_DRAG_LATE_ROLLBACK_PASS'
    write(*,*)'DUST_FV_PRESSURE_DRIFT_CR_WORK_128_POPULATION_DRAG_AND_ROLLBACK_PASS'
  end subroutine

    subroutine make_multifluid_state(x,u)
      integer,parameter::nb=7,nv=167,base=33
      real(real64),parameter::pi=acos(-1d0)
      real(real64),intent(in)::x
      real(real64),intent(out)::u(nv)
      real(real64)::mass(nb),speed_b(nb)
      integer::b,k
      mass=.01d0;mass(1)=.01d0*(1+.2d0*sin(2*pi*x));mass(7)=.005d0*(1+.1d0*cos(2*pi*x))
      speed_b=[-.4d0,.7d0,.1d0,.1d0,.1d0,.1d0,.8d0]
      u=0;u(1)=1+sum(mass);u(2)=.2d0+sum(mass*speed_b)
      u(5)=1.5d0+.3d0+.5d0*.2d0**2+.5d0*sum(mass*speed_b**2)
      do b=1,nb
         k=6+4*(b-1);u(k)=mass(b);u(k+1)=mass(b)*speed_b(b)
      enddo
      u(base+1)=.3d0;u(base+2)=.74d0;u(base+3:base+6)=.001d0*mass(1:4)
      u(base+7:)=mass(7)/128
    end subroutine

  subroutine check_iron_material()
    real(real64)::lo,hi,t,frac(4),m(3),energy,identity(iron_identity_n),ref,a,b
    real(real64)::transitions(3),latent(3),f,els(11),gas0(11),gas1(11),g(2),previous
    real(real64)::old_bins(4),new_bins(4),adjusted(11),ref_bins(4)
    logical::saved_sputtering
    integer::j,k,status
    transitions=[1184d0,1665d0,1809d0];latent=[.900d0,.837d0,13.807d0]*1d10/55.847d0
    call iron_material_identity(identity)
    if(identity(1)/=2.or.any(.not.ieee_is_finite(identity)))stop 180
    call iron_enthalpy_bounds(298.15d0,lo,hi,status)
    if(status/=0.or.lo/=hi.or.abs(lo/(4.507d10/55.847d0)-1)>1d-14)stop 181
    ! These are held-out numerical HD17 integral values, not self round trips.
    call iron_enthalpy_bounds(5d0,lo,hi,status)
    if(status/=0.or.abs(lo/(.06660780153458609d7/55.847d0*identity(9))-1)>1d-12)stop 182
    call iron_enthalpy_bounds(100d0,lo,hi,status)
    if(status/=0.or.abs(lo/(449.22575861920205d7/55.847d0*identity(9))-1)>1d-12)stop 183
    write(*,'(A,3ES16.8)')'FE_LOW_T_JANAF scale,U100/raw_JANAF,U200/raw_JANAF=', &
         identity(9),lo/(.423d10/55.847d0),identity(9)*2229.5914378286066d0/2192d0
    call iron_enthalpy_bounds(500d0,lo,hi,status)
    if(status/=0.or.abs(lo/((5.527d0+4.507d0)*1d10/55.847d0)-1)>1d-14)stop 184
    ! Real mixtures exercise energy absorbed in host solids during a plateau.
    do j=1,3
       m=[0d0,0d0,1d0]
       if(j==2)m=[.2d0,.3d0,.5d0]
       if(j==3)m=[.2d-27,.3d-27,.5d-27]
       do k=1,3
          call iron_mixture_enthalpy(transitions(k),m,lo,hi,status)
          if(status/=0.or.abs((hi-lo)/(m(3)*latent(k))-1)>1d-12)stop 185
          f=.37d0;energy=lo+f*(hi-lo)
          call iron_mixture_state(m,energy,t,frac,status)
          if(status/=0.or.t/=transitions(k).or.abs(frac(k+1)-f)>2d-12)stop 186
          if(abs(sum(frac)-1)>1d-14.or.any(frac<0))stop 187
          call iron_mixture_state(m,lo,t,frac,status)
          if(status/=0.or.t/=transitions(k).or.frac(k)/=1)stop 188
          call iron_mixture_state(m,hi,t,frac,status)
          if(status/=0.or.t/=transitions(k).or.frac(k+1)/=1)stop 189
       enddo
       previous=0
       do k=1,301
          ref=5d0+(3000d0-5d0)*real(k-1,real64)/300
          if(k==1)ref=5d0
          if(k==301)ref=3000d0
          call iron_mixture_enthalpy(ref,m,lo,hi,status)
          if(status/=0.or.lo<=previous)then
             write(*,*)'FE_MONOTONE_REJECT mixture/k/status/T/U/previous=',j,k,status,ref,lo,previous
             stop 190
          endif
          previous=lo
          call iron_mixture_state(m,lo,t,frac,status)
          if(status/=0.or.abs(t-ref)>2d-9)stop 191
       enddo
       call iron_mixture_state(m,hi*1.001d0,t,frac,status)
       if(status==0.or.t/=0.or.any(frac/=0))stop 192
    enddo
    ! Curie point must not acquire fictitious latent energy.
    call iron_enthalpy_bounds(1042d0,lo,hi,status)
    if(status/=0.or.lo/=hi)stop 193
    call iron_enthalpy_bounds(298.15d0*(1-1d-9),a,b,status)
    call iron_enthalpy_bounds(298.15d0*(1+1d-9),lo,hi,status)
    if(status/=0.or.a>=lo.or.abs(a/lo-1)>1d-8)stop 194
    call iron_mixture_state([.4d0,.6d0,0d0],0d0,t,frac,status)
    if(status/=0.or.t/=0.or.any(frac/=0))stop 195 ! DL01 bulk U(0)=0
    do k=1,20
       ref=real(k,real64)/4
       call iron_mixture_enthalpy(ref,[.2d0,.3d0,.5d0],lo,hi,status)
       if(status/=0)stop 196
       call iron_mixture_state([.2d0,.3d0,.5d0],lo,t,frac,status)
       if(status/=0.or.abs(t-ref)>1d-12)stop 197
    enddo
    write(*,*)'DL01_FE_COLD_ZERO_ENERGY_INVERSE_PASS'
    call iron_mixture_state([0d0,0d0,1d0],0d0,t,frac,status)
    if(status/=0.or.abs(t)>1d-20.or.frac(1)/=1)stop 196
    call iron_enthalpy_bounds(3000.0001d0,lo,hi,status)
    if(status==0.or.lo/=0.or.hi/=0)stop 197
    ! Separate Fe uses only Fe remaining AFTER olivine, not an extra metal pool.
    els=olivine_fraction*.6d0;els(3)=.2d0;els(11)=els(11)+.1d0;g=[.05d0,.3d0]
    call dust_gas_elements(els,g,gas0,status)
    if(status/=0)stop 198
    call dust_gas_elements(els,g,gas1,status,metallic_iron=.1d0)
    if(status/=0.or.any(gas0(1:10)/=gas1(1:10)).or.abs(gas0(11)-gas1(11)-.1d0)>1d-14)stop 199
    if(abs(sum(gas1)+sum(g)+.1d0-sum(els))>1d-14)stop 200
    call dust_gas_elements(els,g,gas1,status,metallic_iron=gas0(11)*1.001d0)
    if(status==0.or.any(gas1/=0))stop 201
    call dust_gas_elements(els,g,gas1,status,metallic_iron=-.1d0)
    if(status==0.or.any(gas1/=0))stop 202
    ! Existing growth/size exchange must not accrete Fe locked in the new solid.
    saved_sputtering=dust_sputtering;dust_sputtering=.false.
    els=olivine_fraction*.6d0;els(3)=.2d0;els(11)=els(11)+.1d0;els(1)=9.1d0
    old_bins=[.01d0,.01d0,.03d0,.03d0]
    call dust_size_step_reserved_iron(10d0,.9d0,els,old_bins,.25d0,100d0,1d-24,1d16,new_bins,status)
    if(status/=0.or.sum(new_bins)<=sum(old_bins))stop 203
    call dust_gas_elements(els,[sum(new_bins(1:2)),sum(new_bins(3:4))],gas1,status,.25d0)
    if(status/=0.or.minval(gas1)<0.or.abs(sum(gas1)+sum(new_bins)+.25d0-sum(els))>1d-13)stop 204
    adjusted=els;adjusted(11)=adjusted(11)-.25d0
    call dust_size_step(10d0,.9d0-.25d0,adjusted,old_bins,100d0,1d-24,1d16,ref_bins,status)
    if(status/=0.or.any(ref_bins/=new_bins))stop 205
    call dust_size_step_reserved_iron(10d0,.9d0,els,old_bins,.4d0,100d0,1d-24,1d16,new_bins,status)
    if(status==0.or.any(new_bins/=old_bins))stop 206
    call dust_size_step_reserved_iron(10d0,.9d0,els,old_bins,0d0,100d0,1d-24,1d16,new_bins,status)
    call dust_size_step(10d0,.9d0,els,old_bins,100d0,1d-24,1d16,ref_bins,j)
    if(status/=0.or.j/=0.or.any(ref_bins/=new_bins))stop 207
    call dust_size_step_reserved_iron(10d0,.9d0,els,old_bins,.25d0,100d0,1d-24,1d16,new_bins,status, &
         pah_hc=[.01d0,.1d0])
    if(status/=0)stop 242
    call dust_gas_elements(els,[sum(new_bins(1:2)),sum(new_bins(3:4))],gas1,status,.25d0,[.01d0,.1d0])
    if(status/=0.or.abs(sum(gas1)+sum(new_bins)+.25d0+.11d0-sum(els))>1d-13)stop 243
    adjusted=els;adjusted(11)=adjusted(11)-.25d0
    adjusted(1)=adjusted(1)-.01d0;adjusted(3)=adjusted(3)-.1d0
    call dust_size_step(10d0,.9d0-.25d0-.1d0,adjusted,old_bins,100d0,1d-24,1d16,ref_bins,status)
    if(status/=0.or.any(ref_bins/=new_bins))stop 244
    call dust_size_step_reserved_iron(10d0,.9d0,els,old_bins,.25d0,100d0,1d-24,1d16,new_bins,status, &
         pah_hc=[.01d0,.21d0])
    if(status==0.or.any(new_bins/=old_bins))stop 245
    call dust_gas_elements(els,[sum(old_bins(1:2)),sum(old_bins(3:4))],gas1,status,.25d0,[9.2d0,.1d0])
    if(status==0.or.any(gas1/=0))stop 246
    dust_sputtering=saved_sputtering
    write(*,*)'FE_JANAF_PHASE_LATENT_MIXTURE_INVERSE_DOMAIN_ELEMENT_RESERVATION_PASS'
    write(*,*)'FE_RESERVED_CARBON_OLIVINE_GROWTH_BUDGET_OLD_PATH_PARITY_PASS'
    write(*,*)'FE_PAH_H_C_RESERVED_GROWTH_CHEMICAL_INVENTORY_PASS'
  end subroutine

  subroutine check_drag()
    real(real64)::rho(3),pg(3),pd(3,3),ng(3),nd(3,3),ts(3),h,k0,k1,vg(3),vd(3,3),q(3)
    real(real64)::f,df,flo,fhi,t,radii(3),material_rho(3),drift(3,3),ciso,linear_ts
    integer::j,status
    rho=[.2d0,.3d0,0d0];pg=[1d0,2d0,-3d0]
    pd(:,1)=[.8d0,-.2d0,.6d0];pd(:,2)=[-.6d0,.9d0,.3d0];pd(:,3)=0
    ts=[.1d0,10d0,1d0]
    k0=.5d0*sum(pg**2)
    do j=1,2
       k0=k0+.5d0*sum(pd(:,j)**2)/rho(j)
    enddo
    call dust_drag_step(1d0,rho,pg,pd,ts,2d0,ng,nd,h,status)
    if(status/=0.or.h<=0)stop 121
    q=pg+sum(pd,dim=2)
    if(maxval(abs(ng+sum(nd,dim=2)-q))>1d-13*maxval(abs(q)))stop 122
    k1=.5d0*sum(ng**2)
    do j=1,2
       k1=k1+.5d0*sum(nd(:,j)**2)/rho(j)
    enddo
    if(abs(k1+h-k0)>1d-13*k0)stop 123
    call dust_drag_step(1d0,rho,pg,pd,ts,1d200,ng,nd,h,status)
    if(status/=0)stop 124
    vg=ng
    do j=1,2
       vd(:,j)=nd(:,j)/rho(j)
       if(maxval(abs(vd(:,j)-vg))>1d-13*maxval(abs(vg)))stop 125
    enddo
    call dust_drag_step(1d0,rho,pg,pd,ts,0d0,ng,nd,h,status)
    if(status/=0.or.any(ng/=pg).or.any(nd/=pd).or.h/=0)stop 126
    ts(2)=-1
    call dust_drag_step(1d0,rho,pg,pd,ts,1d0,ng,nd,h,status)
    if(status==0.or.any(ng/=pg).or.any(nd/=pd).or.h/=0)stop 127
    radii=[1d-6,1d-5,1d-4];material_rho=3;drift=0
    ciso=sqrt(1.380649d-16*100d0/1.67262192369d-24)
    linear_ts=sqrt(acos(-1d0)/8)*3d0*1d-6/(1d-24*ciso)
    call dust_epstein_stopping_time(1d-24,100d0,1.67262192369d-24,1d10, &
         radii,material_rho,drift,ts,status)
    if(status/=0.or.maxval(abs(ts/(linear_ts*[1d0,10d0,100d0])-1))>1d-13)stop 131
    drift(1,:)=100d0*ciso
    call dust_epstein_stopping_time(1d-24,100d0,1.67262192369d-24,1d10, &
         radii,material_rho,drift,ts,status)
    f=sqrt(1+(9*acos(-1d0)/128)*1d4)
    if(status/=0.or.abs(ts(1)*f/linear_ts-1)>1d-13)stop 132
    call dust_epstein_stopping_time(1d-24,100d0,1.67262192369d-24,1d-8, &
         radii,material_rho,drift,ts,status)
    if(status==0.or.any(ts/=0))stop 133
    write(*,*)'DUST_EPSTEIN_CGS_LINEAR_SUPERSONIC_DOMAIN_PASS'
    t=1700d0
    call dust_olivine_surface_flux(t,f,df,status)
    if(status/=0.or.min(f,df)<=0)stop 128
    call dust_olivine_surface_flux(t-.01d0,flo,h,status)
    call dust_olivine_surface_flux(t+.01d0,fhi,h,status)
    if(abs((fhi-flo)/.02d0/df-1)>1d-6)stop 129
    if(abs(f/(sum(exp([20.06d0,16.29d0,22.23d0]-[84780d0,74420d0,90590d0]/t))/3)-1)>1d-14)stop 130
    write(*,*)'DUST_DRAG_MOMENTUM_KINETIC_HEAT_ROLLBACK_AND_OLIVINE_RATE_PASS'
    call check_mixture_receivers()
  end subroutine

  subroutine check_mixture_receivers()
    ! Seven carriers exercise six large-grain bins plus an independent PAH
    ! population; no radius, composition or physical source is assigned here.
    real(real64)::solid(7),relative(3,7),pd(3,7),pg(3),rg,dk,dk1,total_p(3),next_p(3)
    real(real64)::next_j(3,7),impulse(3,7),debit(7),heat(7),ts(7),q,work,old_p(3),old_j(3,7)
    real(real64)::old_heat(7),kg0,kg1,velocity(3)
    integer::b,status
    solid=[.02d0,.03d0,.04d0,.05d0,.01d0,.015d0,0d0]
    total_p=[2d0,-1d0,.5d0];relative=0
    do b=1,6
       relative(:,b)=solid(b)*[real(b,real64),-.2d0*real(b,real64),.5d0]
    enddo
    rg=0;pg=0;pd=0;dk=0
    call dust_mixture_split(1d0,total_p,solid,relative,rg,pg,pd,dk,status)
    if(status/=0.or.maxval(abs(pg+sum(pd,dim=2)-total_p))>1d-14)stop 220
    kg0=.5d0*sum(pg**2)/rg
    do b=1,6
       kg0=kg0+.5d0*sum(pd(:,b)**2)/solid(b)
    enddo
    if(abs(kg0-.5d0*sum(total_p**2)-dk)>1d-13*kg0)stop 221
    ts=[.01d0,.1d0,1d0,10d0,100d0,1000d0,1d0];next_j=relative;q=-1
    call dust_mixture_drag(1d0,solid,relative,ts,2d0,next_j,q,status)
    if(status/=0.or.q<=0)stop 222
    call dust_mixture_split(1d0,total_p,solid,next_j,rg,pg,pd,dk1,status)
    if(status/=0.or.abs(dk1+q-dk)>1d-12*dk)stop 223
    call dust_mixture_split(1d0,total_p*1d12,solid,next_j,rg,pg,pd,kg1,status)
    if(status/=0.or.kg1/=dk1)stop 224 ! Bulk-independent internal kinetic energy.
    impulse=0;debit=0
    do b=1,6
       impulse(:,b)=solid(b)*[.4d0,-.1d0,.2d0]
       velocity=total_p+relative(:,b)/solid(b)
       work=dot_product(velocity,impulse(:,b))+.5d0*sum(impulse(:,b)**2)/solid(b)
       debit(b)=work+.1d0*solid(b)
    enddo
    next_p=total_p;next_j=relative;heat=-1
    call dust_mixture_radiation_kick(1d0,total_p,solid,relative,impulse,debit,next_p,next_j,heat,status)
    if(status/=0.or.maxval(abs(next_p-total_p-sum(impulse,dim=2)))>1d-14)stop 225
    call dust_mixture_split(1d0,next_p,solid,next_j,rg,pg,pd,dk1,status)
    kg1=.5d0*sum(next_p**2)+dk1
    if(status/=0.or.abs(kg1-kg0+sum(heat)-sum(debit))>1d-13*kg0)stop 226
    if(maxval(abs(heat-.1d0*solid))>1d-14)stop 227
    old_p=next_p;old_j=next_j;old_heat=heat
    debit=0 ! Reject a force which adds kinetic energy without a photon debit.
    call dust_mixture_radiation_kick(1d0,total_p,solid,relative,impulse,debit,next_p,next_j,heat,status)
    if(status==0.or.any(next_p/=old_p).or.any(next_j/=old_j).or.any(heat/=old_heat))stop 228
    solid(7)=1 ! Invalid total solid mass; no half-committed momentum or heat.
    call dust_mixture_radiation_kick(1d0,total_p,solid,relative,impulse,debit,next_p,next_j,heat,status)
    if(status==0.or.any(next_p/=old_p).or.any(next_j/=old_j).or.any(heat/=old_heat))stop 229
    write(*,*)'DUST_MIXTURE_DRAG_RADIATION_MOMENTUM_ENERGY_ATOMIC_PASS'
  end subroutine

  subroutine check_stochastic()
    real(real64)::levels(3),up(3,3),down(3),prob(3),absorbed,emitted
    real(real64)::over_rate(3),over_power(3),represented,photons(2),rates(2)
    real(real64)::modes(102),u(5),cv(5),temperatures(5)
    integer::status,i,j
    levels=[0d0,1d-12,3d-12];up=0;up(2,1)=2;up(3,1)=1;up(3,2)=4
    down=[0d0,3d0,5d0]
    call dust_stochastic_equilibrium(levels,up,down,prob,absorbed,emitted,status)
    if(status/=0.or.maxval(abs(prob-1d0/3))>1d-14)stop 134
    if(abs(absorbed/emitted-1)>1d-14)stop 135
    up=0
    call dust_stochastic_equilibrium(levels,up,down,prob,absorbed,emitted,status)
    if(status/=0.or.any(prob/=[1d0,0d0,0d0]).or.absorbed/=0.or.emitted/=0)stop 136
    up(2,1)=1d200;up(3,2)=1d200;down(2:)=1d-200
    call dust_stochastic_equilibrium(levels,up,down,prob,absorbed,emitted,status)
    if(status/=0.or.any(prob<0).or.prob(3)/=1d0.or.sum(prob)/=1d0)stop 137
    if(abs(absorbed/emitted-1)>1d-11.or.abs(emitted/2d-212-1)>1d-11)stop 139
    up(1,2)=1
    call dust_stochastic_equilibrium(levels,up,down,prob,absorbed,emitted,status)
    if(status==0.or.any(prob/=0).or.absorbed/=0.or.emitted/=0)stop 138
    write(*,*)'DUST_STOCHASTIC_TC_STATIONARY_DARK_STRONG_FIELD_ENERGY_PASS'
    photons=[.5d-12,2d-12];rates=[2d0,3d0]
    call dust_stochastic_photon_rates(levels,photons,rates,up,over_rate,over_power,status)
    if(status/=0.or.any(over_rate/=[0d0,0d0,5d0]))stop 140
    do i=1,3
       represented=over_power(i)
       do j=i+1,3
          represented=represented+up(j,i)*(levels(j)-levels(i))
       enddo
       if(abs(represented/sum(photons*rates)-1)>1d-14)stop 141
    enddo
    ! The upper boundary is not silently cooled, clipped, or deleted.
    if(abs(over_power(3)/7d-12-1)>1d-14.or.any(up(:,3)/=0))stop 142
    photons(2)=-1
    call dust_stochastic_photon_rates(levels,photons,rates,up,over_rate,over_power,status)
    if(status==0.or.any(up/=0).or.any(over_rate/=0).or.any(over_power/=0))stop 143
    write(*,*)'DUST_STOCHASTIC_PHOTON_GROUP_MOMENT_AND_EXPLICIT_OVERFLOW_PASS'
    call dust_pah_modes(24,12,modes,status)
    if(status/=0.or.any(modes<=0))stop 144
    if(abs(modes(1)/(6.62607015d-27*2.99792458d10*600*sqrt(.5d0/22))-1)>1d-14)stop 145
    temperatures=[0d0,199.999d0,200d0,200.001d0,1d8]
    call dust_vibrational_curve(modes,temperatures,u,cv,status)
    if(status/=0.or.u(1)/=0.or.cv(1)/=0.or.any(u(2:)<=u(:4)))stop 146
    if(abs((u(4)-u(2))/.002d0/cv(3)-1)>1d-8)stop 147
    if(abs(cv(5)/(102*1.380649d-16)-1)>1d-9)stop 148
    write(*,*)'DUST_PAH_C24H12_MODE_COUNT_ENERGY_HEAT_CAPACITY_PASS'
    call check_stochastic_time()
    call check_pah_radiation()
  end subroutine

  subroutine check_pah_charge_events()
    real(real64),parameter::ev=1.602176634d-12,ip=7.02d0,kb=1.380649d-16
    real(real64)::pop(5,2),old(5,2),captured(3,2),vib(3,2),level(5),ne,heat,stored,expected,nrec
    real(real64)::saved_ne,saved_heat,saved_stored,dt,delta,a,nmin,nremain,charge,energy
    integer::status,j
    level=[0d0,2d0,8d0,16d0,32d0]*ev
    pop=0;pop(1,1)=.6d0;pop(2,1)=.4d0;old=pop
    ne=0;heat=0;stored=0;captured=0;captured(:,1)=[.05d0,.1d0,.2d0];vib=-1
    call pah_coronene_photoionize([2d0,7.02d0,10d0],captured,pop,ne,vib,heat,stored,status)
    expected=.1d0*.8d0*exp(-.00128d0*(7.02d0-14.89d0)**4) &
         +.2d0*.8d0*exp(-.00128d0*(10d0-14.89d0)**4)
    if(status/=0.or.abs(ne-expected)>1d-15.or.abs(sum(pop)-1)>1d-15)stop 240
    if(abs(sum(pop(:,2))-ne)>1d-15.or.maxval(abs(sum(pop,dim=2)-sum(old,dim=2)))>1d-15)stop 241
    energy=sum((captured(:,1)-vib(:,1))*[2d0,7.02d0,10d0])*ev
    if(abs(energy-heat-stored)>1d-27.or.vib(1,1)/=captured(1,1).or.any(vib<0))stop 242
    old=pop;saved_ne=ne;saved_heat=heat;saved_stored=stored
    call pah_coronene_photoionize([2d0,7.02d0,14d0],captured,pop,ne,vib,heat,stored,status)
    if(status==0.or.any(pop/=old).or.ne/=saved_ne.or.heat/=saved_heat.or.stored/=saved_stored)stop 243
    captured(:,1)=1d5
    call pah_coronene_photoionize([2d0,7.02d0,10d0],captured,pop,ne,vib,heat,stored,status)
    if(status==0.or.any(pop/=old).or.ne/=saved_ne.or.heat/=saved_heat.or.stored/=saved_stored)stop 244
    ! Actual photoelectrons recombine with the cations just produced. Equal
    ! abundances obey n(t)=n(0)/(1+alpha*n(0)*t), not frozen-ne exponential.
    dt=1d5;nrec=ne-nesolution(ne,1d-5*dt)
    call pah_coronene_recombine(level,300d0,dt,pop,ne,heat,stored,status)
    if(status/=0.or.abs(saved_ne-ne-nrec)>1d-15.or.abs(sum(pop(:,2))-ne)>1d-15)stop 245
    energy=dot_product(level,sum(pop-old,dim=2))+heat-saved_heat+stored-saved_stored
    if(abs(energy)>1d-27.or.abs(heat-saved_heat+nrec*kb*300)>1d-27)stop 246
    if(abs(sum(pop)-sum(old))>1d-15.or.abs(stored-ip*ev*sum(pop(:,2)))>1d-27)stop 247
    ! Both signs of n_e-n_PAH+ and a trace step; verify analytic counts and
    ! combined charge while the neutral vibrational reservoir takes IP+kT.
    do j=1,3
       pop=0;pop(1,2)=1d-3;ne=1d-3*real(j,real64)/2;old=pop;saved_ne=ne
       heat=0;stored=ip*ev*sum(pop(:,2));saved_stored=stored
       dt=1d8;if(j==2)dt=1d-4
       delta=abs(ne-.001d0);a=1d-5*dt;nmin=min(ne,.001d0)
       if(delta==0)then
          nremain=nesolution(nmin,a)
       else
          nremain=nmin*exp(-a*delta)/(1+nmin/delta*(1-exp(-a*delta)))
       endif
       charge=.001d0-ne;nrec=nmin-nremain
       call pah_coronene_recombine(level,300d0,dt,pop,ne,heat,stored,status)
       if(status/=0.or.abs(ne-(saved_ne-nrec))>1d-17)stop 248
       if(abs(sum(pop(:,2))-ne-charge)>1d-17.or.abs(sum(pop)-.001d0)>1d-17)stop 249
       energy=dot_product(level,sum(pop-old,dim=2))+heat+stored-saved_stored
       if(abs(energy)>1d-29)stop 250
    enddo
    pop=0;pop(5,2)=1;old=pop;ne=1;heat=0;stored=ip*ev
    call pah_coronene_recombine(level,300d0,1d5,pop,ne,heat,stored,status)
    if(status==0.or.any(pop/=old).or.ne/=1.or.heat/=0.or.stored/=ip*ev)stop 251
    call pah_coronene_recombine(level,300d0,0d0,pop,ne,heat,stored,status)
    if(status/=0.or.any(pop/=old).or.ne/=1.or.heat/=0.or.stored/=ip*ev)stop 252
    write(*,*)'PAH_CHARGE_PHOTOELECTRON_RECOMBINATION_ENERGY_NUMBER_ROLLBACK_PASS'
  end subroutine

  real(real64) function nesolution(ne,a) result(n)
    real(real64),intent(in)::ne,a
    n=ne/(1+a*ne)
  end function

  subroutine check_pah_radiation()
    integer,parameter::nl=128
    real(real64),parameter::ev=1.602176634d-12
    type(pah_radiative_model)::model,physical
    type(dust_ir_diagnostics)::diag,saved_diag
    real(real64)::level(nl),old(nl),p(nl),saved(nl),rate(3),u,t,over,hc(2),old_u,initial
    real(real64)::ir(3,2,2),ir_saved(3,2,2),pop(nl,2),pop_saved(nl,2),emitted(3,2),em_saved(3,2)
    real(real64)::direction(3,2),weight(2),source(1,2),before,after
    real(real64),allocatable::opt_e(:),opt_w(:),opt_c(:),opt_rate(:),physical_ir(:,:,:),physical_em(:,:)
    real(real64),allocatable::saved_ir(:,:,:),saved_em(:,:)
    real(real64)::primary_photons(1,2,2),saved_primary(1,2,2),total_u,physical_balance
    integer::i,status,neighbor(6,2)
    character(len=2048)::path
    level(1)=0
    do i=2,nl
       level(i)=ev*.012d0*(64d0/.012d0)**(real(i-2,real64)/(nl-2))
    enddo
    ! Artificial three-band coefficients exercise the transaction algebra;
    ! the optional original-table case below is the physical spectral check.
    call pah_radiative_prepare(model,24,12,level,[.002d0,.01d0,.1d0], &
         [.004d0,.01d0,.1d0],[1d-17,1d-16,1d-16],status)
    if(status/=0)stop 250
    old=0;old(1)=1d-8;p=-1;rate=-1;u=-1;t=-1
    call pah_absorbed_step(model,[2d0],[2d-3*ev*sum(old)],old,1d0,p,rate,u,t,over,status)
    if(status/=0.or.any(p<0).or.any(rate<0).or.sum(rate)<=0)stop 251
    if(abs(sum(p)/sum(old)-1)>1d-13)stop 252
    if(abs(u+sum(rate)-2d-3*ev*sum(old))>1d-12*2d-3*ev*sum(old))stop 253
    saved=p
    call pah_absorbed_step(model,[.5d0],[2d-3*ev*sum(old)],old,1d0,p,rate,u,t,over,status)
    if(status/=0.or.maxval(abs(p-saved))<1d-8*sum(old))stop 254
    call pah_solid_inventory(model,p,hc,u,status)
    if(status/=0.or.abs(hc(2)/hc(1)-(24*12.011d0)/(12*1.008d0))>1d-13)stop 255
    ! Independent captured count: nominal 2 eV photons depositing only .5 eV
    ! per event must reproduce the .5 eV excitation case, not lose events.
    saved=p
    call pah_absorbed_step(model,[2d0],[2d-3*ev*sum(old)],old,1d0,p,rate,u,t,over,status, &
         captured_photons=[4d-3*sum(old)])
    if(status/=0.or.maxval(abs(p-saved))>1d-13*sum(old))stop 306
    saved=p;initial=u
    call pah_absorbed_step(model,[2d0],[2d-3*ev*sum(old)],old,1d0,p,rate,u,t,over,status, &
         captured_photons=[0d0])
    if(status==0.or.any(p/=saved).or.u/=initial)stop 307
    call pah_absorbed_step(model,[2d0],[0d0],old,1d0,p,rate,u,t,over,status, &
         captured_photons=[4d-3*sum(old)])
    if(status/=0.or.any(p/=old).or.u/=0.or.any(rate/=0))stop 308
    p=saved;u=initial
    saved=p;old_u=u
    call pah_absorbed_step(model,[2d0],[0d0],saved,1d0,p,rate,u,t,over,status)
    if(status/=0.or.u>=old_u.or.abs(u+sum(rate)-old_u)>1d-12*old_u)stop 256
    saved=p;initial=u
    call pah_absorbed_step(model,[128d0],[128*ev*sum(old)],old,1d0,p,rate,u,t,over,status)
    if(status/=5.or.over<=0.or.any(p/=saved).or.u/=initial)stop 257
    ! Two reciprocal periodic cells, real transport and local reabsorption.
    direction=0;direction(1,:)=[1d0,-1d0];weight=.5d0
    neighbor=0;neighbor(1:2,1)=2;neighbor(1:2,2)=1
    pop=0;pop(1,:)=1d5;ir=0;emitted=0
    source=0;source(1,1)=2d-3*ev*sum(pop(:,1))
    diag=dust_ir_diagnostics()
    call pah_radiative_advance(model,direction,weight,neighbor,1d12,1d0,2.99792458d10,[2d0],source, &
         ir,pop,emitted,diag,status,1d-10,100)
    if(status/=0.or.sum(ir)<=0.or.sum(emitted)<=0.or.diag%absorbed_erg<=0)stop 258
    if(diag%balance_relative>1d-10.or.diag%local_relative>1d-10)stop 259
    before=sum(source)
    after=sum(matmul(level,pop))+.5d0*sum(ir)
    if(abs(after/before-1)>1d-10.or.any(abs(sum(pop,dim=1)/1d5-1)>1d-12))stop 260
    ! No hidden relaxation reset, and spatial photons heat the second cell.
    source=0
    call pah_radiative_advance(model,direction,weight,neighbor,1d12,1d0,2.99792458d10,[2d0],source, &
         ir,pop,emitted,diag,status,1d-10,100)
    if(status/=0.or.sum(pop(2:,2))<=0)stop 261
    if(abs((sum(matmul(level,pop))+.5d0*sum(ir))/before-1)>2d-10)stop 262
    ir_saved=ir;pop_saved=pop;em_saved=emitted;saved_diag=diag
    source(1,1)=128*ev*sum(pop(:,1))
    call pah_radiative_advance(model,direction,weight,neighbor,1d12,1d0,2.99792458d10,[128d0],source, &
         ir,pop,emitted,diag,status,1d-10,100)
    if(status==0.or.any(ir/=ir_saved).or.any(pop/=pop_saved).or.any(emitted/=em_saved))stop 263
    if(diag%balance_relative/=saved_diag%balance_relative.or.diag%iterations/=saved_diag%iterations)stop 264
    write(*,*)'DUST_PAH_SPECTRAL_POPULATION_IR_TRANSPORT_REABSORPTION_HC_ROLLBACK_PASS'
    call get_environment_variable('SNRT_PAH_NEUTRAL_TABLE',path,status=status)
    if(status/=0.or.len_trim(path)==0)return
    call pah_neutral_optics(trim(path),24,opt_e,opt_w,opt_c,status)
    if(status/=0)stop 265
    if(size(opt_e)/=1201.or.opt_e(1)>.00125d0.or.opt_e(1201)<1239d0)stop 266
    call pah_radiative_prepare(physical,24,12,level,opt_e,opt_w,opt_c,status)
    if(status/=0)stop 267
    allocate(opt_rate(size(opt_e)));opt_rate=0;p=old;u=0;t=0
    call pah_absorbed_step(physical,[2d0],[2d-3*ev*sum(old)],old,1d0,p,opt_rate,u,t,over,status)
    if(status/=0.or.sum(opt_rate)<=0.or.u<=0)stop 268
    if(abs(u+sum(opt_rate)-2d-3*ev*sum(old))>1d-12*2d-3*ev*sum(old))stop 269
    write(*,'(A,4ES14.6)')'DUST_PAH_ORIGINAL_NEUTRAL_TABLE radiation,U,T,overflow=',sum(opt_rate),u,t,over
    write(*,*)'DUST_PAH_C24H12_ORIGINAL_OPTICS_1201_NODES_ENERGY_PASS'
    allocate(physical_ir(size(opt_e),2,2),physical_em(size(opt_e),2))
    physical_ir=0;physical_em=0;pop=0;pop(1,:)=1d3
    primary_photons=0;primary_photons(1,:,1)=1d3
    before=.5d0*sum(primary_photons)*2*ev
    call pah_photon_advance(physical,direction,weight,neighbor,1d12,1d0,2.99792458d10,[2d0], &
         primary_photons,physical_ir,pop,physical_em,diag,status,1d-10,100)
    if(status/=0.or.sum(physical_ir)<=0.or.sum(primary_photons)>=2d3)stop 270
    total_u=sum(matmul(level,pop))
    after=.5d0*sum(primary_photons)*2*ev+.5d0*sum(physical_ir)+total_u
    physical_balance=abs(after/before-1)
    if(physical_balance>1d-12)stop 271
    ! Data coverage failure must not debit a high-energy group or any state.
    saved_primary=primary_photons;saved_ir=physical_ir;saved_em=physical_em;pop_saved=pop;saved_diag=diag
    call pah_photon_advance(physical,direction,weight,neighbor,1d12,1d0,2.99792458d10,[5000d0], &
         primary_photons,physical_ir,pop,physical_em,diag,status,1d-10,100)
    if(status==0.or.any(primary_photons/=saved_primary).or.any(pop/=pop_saved))stop 272
    if(any(physical_ir/=saved_ir).or.any(physical_em/=saved_em))stop 273
    if(diag%iterations/=saved_diag%iterations)stop 274
    ! A LATE excitation overflow also rolls back the already staged debit.
    call pah_photon_advance(physical,direction,weight,neighbor,1d12,1d0,2.99792458d10,[128d0], &
         primary_photons,physical_ir,pop,physical_em,diag,status,1d-10,100)
    if(status/=5.or.any(primary_photons/=saved_primary).or.any(pop/=pop_saved))stop 275
    if(any(physical_ir/=saved_ir).or.any(physical_em/=saved_em))stop 276
    write(*,'(A,3ES14.6)')'DUST_PAH_PHYSICAL_PRIMARY_IR_U_BALANCE=',physical_balance, &
         saved_diag%balance_relative,saved_diag%local_relative
    write(*,*)'DUST_PAH_ORIGINAL_OPTICS_PRIMARY_DEBIT_IR_POPULATION_ATOMIC_PASS'
    call check_pah_charged_receiver(physical,level,opt_e,opt_w,opt_c)
    call check_pah_hydrogen_physical(physical,level)
  end subroutine

  subroutine check_pah_hydrogen_physical(physical,levels)
    use dust_pah_hydrogen
    type(pah_radiative_model),intent(in)::physical
    real(real64),intent(in)::levels(:)
    type(pah_hydrogen_model)::model
    type(pah_hydrogen_model)::charged_models(2)
    type(pah_radiative_model)::ionized
    real(real64)::old(size(levels),0:13),next(size(levels),0:13),saved(size(levels),0:13)
    real(real64),allocatable::e(:),bands(:,:),down(:),rate(:),saved_rate(:)
    real(real64),allocatable::ie(:),iw(:),ic(:)
    real(real64)::charge_old(size(levels),0:13,2),charge_next(size(levels),0:13,2),captures(1,2)
    real(real64)::ne,initial_ne,initial_h,total_before,total_after
    real(real64)::gas_h,db,stored,energy,initial,heat,h_before,bond(0:13),attach(0:13),axis(0:13)
    real(real64),parameter::ev=1.602176634d-12
    integer::status,nh,j
    character(len=2048)::path
    call pah_hydrogen_prepare(model,physical,0,.0005d0*ev,status)
    if(status/=0)stop 522
    call pah_radiative_coefficients(physical,e,bands,down,status)
    if(status/=0)stop 523
    allocate(rate(size(bands,1)));rate=0
    do nh=0,13
       axis(nh)=nh
    enddo
    old=0;old(1,12)=1d-8;next=old;gas_h=0;db=0
    initial=12d0*1d-11*ev
    call pah_hydrogen_absorbed_step(model,[12d0],[1d-11],old,1d0,next,gas_h,rate,db,status)
    if(status/=0.or.gas_h<=0.or.db<=0.or.sum(rate)<=0)then
       write(*,*)'PAH_H_PHYSICAL_FAILURE=',status,gas_h,db,sum(rate)
       stop 524
    endif
    energy=sum(matmul(levels,next-old))
    if(abs((energy+sum(rate)+db)/initial-1)>1d-10)stop 525
    if(abs(gas_h+dot_product(axis,sum(next,dim=1))-12d-8)>1d-20)stop 526
    write(*,'(A,4ES22.14)')'PAH_H_PHYSICAL balance,H_released,bond,IR=', &
         (energy+sum(rate)+db)/initial-1,gas_h,db,sum(rate)
    saved=next;saved_rate=rate;h_before=gas_h;stored=db
    call pah_hydrogen_absorbed_step(model,[20d0],[1d-11],old,1d0,next,gas_h,rate,db,status)
    if(status==0.or.any(next/=saved).or.any(rate/=saved_rate).or.gas_h/=h_before.or.db/=stored)stop 527
    ! Below-barrier initial excitation, with no captures: no H can be lost.
    old=0;j=minloc(abs(levels-2*ev),dim=1);old(j,12)=1d-8;next=old;gas_h=0;db=0
    call pah_hydrogen_absorbed_step(model,[2d0],[0d0],old,10d0,next,gas_h,rate,db,status)
    if(status/=0.or.gas_h/=0.or.db/=0.or.sum(rate)<=0)stop 528
    ! Dehydrogenated carbon skeletons remain explicitly counted; they are
    ! not turned into atomic C by this hydrogen chemistry receiver.
    old=0;old(j,0)=1d-8;next=old;gas_h=0;db=0
    call pah_hydrogen_absorbed_step(model,[2d0],[0d0],old,10d0,next,gas_h,rate,db,status)
    if(status/=0.or.gas_h/=0.or.db/=0.or.abs(sum(next(:,0))/1d-8-1)>1d-12)stop 529
    write(*,*)'PAH_ORIGINAL_OPTICS_MULTIPHOTON_IR_H_LOSS_HC_BOND_ENERGY_PASS'
    call get_environment_variable('SNRT_PAH_ION_TABLE',path,status=status)
    if(status/=0.or.len_trim(path)==0)return
    call pah_charge_optics(trim(path),24,1,ie,iw,ic,status)
    if(status/=0)stop 530
    call pah_radiative_prepare(ionized,24,12,levels,ie,iw,ic,status)
    if(status/=0)stop 531
    charged_models(1)=model
    call pah_hydrogen_prepare(charged_models(2),ionized,1,.0005d0*ev,status)
    if(status/=0)stop 532
    call pah_hydrogen_parameters(1,bond,attach,status)
    charge_old=0;charge_old(1,12,1)=1d-8;charge_old(1,11,2)=5d-9
    charge_next=charge_old;captures=0;captures(1,1)=1d-11
    ne=1d-4;initial_ne=ne;gas_h=100;initial_h=gas_h;heat=0;rate=0
    total_before=dot_product(bond,sum(sum(charge_old,dim=3),dim=1))+sum(charge_old(:,:,2))*7.02d0*ev
    call pah_hydrogen_charged_step(charged_models,[12d0],captures,charge_old,1d8,300d0, &
         charge_next,ne,gas_h,heat,rate,status)
    if(status/=0)then
       write(*,*)'PAH_H_CHARGE_COUPLED_FAILURE',status
       stop 533
    endif
    total_after=sum(matmul(levels,sum(charge_next,dim=3))) &
         +dot_product(bond,sum(sum(charge_next,dim=3),dim=1))+sum(charge_next(:,:,2))*7.02d0*ev
    energy=total_after-total_before+heat+1d8*sum(rate)-sum(captures)*12*ev
    if(abs(energy)>1d-10*abs(total_before))stop 534
    if(abs(sum(charge_next(:,:,2))-sum(charge_old(:,:,2))-ne+initial_ne)>1d-16)stop 535
    if(abs(gas_h-initial_h+dot_product(axis,sum(sum(charge_next-charge_old,dim=3),dim=1)))>1d-12)stop 536
    write(*,'(A,4ES22.14)')'PAH_H_CHARGE_PHYSICAL balance,gasH_delta,ne_delta,heat=', &
         energy/abs(total_before),gas_h-initial_h,ne-initial_ne,heat
    write(*,*)'PAH_H_PHOTOELECTRON_RECOMBINATION_H_ATTACHMENT_SHARED_IR_PASS'
    block
      type(pah_hydrogen_model)::molecular(2)
      real(real64)::molecules,h2_initial,keep_h,keep_ne,keep_heat,keep_h2
      real(real64)::zero_result(size(levels),0:13,2),keep_rate(size(rate)),nuclear_error
      real(real64),allocatable::old_id(:),new_id(:)
      call pah_hydrogen_prepare(molecular(1),physical,0,.0005d0*ev,status,h2_capture=.true.)
      if(status/=0)stop 600
      call pah_hydrogen_prepare(molecular(2),ionized,1,.0005d0*ev,status,h2_capture=.true.)
      if(status/=0)stop 601
      call pah_hydrogen_identity(charged_models(2),old_id,status)
      call pah_hydrogen_identity(molecular(2),new_id,status)
      if(status/=0.or.size(new_id)/=size(old_id)+15)stop 602
      if(any(new_id(2:size(old_id))/=old_id(2:)))stop 603
      zero_result=charge_next;keep_ne=ne;keep_h=gas_h;keep_heat=heat;keep_rate=rate
      ne=initial_ne;gas_h=initial_h;molecules=0;charge_next=charge_old
      call pah_hydrogen_charged_step(molecular,[12d0],captures,charge_old,1d8,300d0, &
           charge_next,ne,gas_h,heat,rate,status,gas_h2=molecules)
      if(status/=0.or.any(charge_next/=zero_result).or.ne/=keep_ne.or.gas_h/=keep_h.or.heat/=keep_heat)stop 604
      if(any(rate/=keep_rate).or.molecules/=0)stop 605
      charge_old=0;charge_old(1,10,2)=1d-3;charge_next=charge_old;captures=0
      gas_h=.005d0;initial_h=gas_h;ne=0;initial_ne=ne;molecules=.3d0;h2_initial=molecules
      total_before=dot_product(bond,sum(sum(charge_old,dim=3),dim=1))+sum(charge_old(:,:,2))*7.02d0*ev
      call pah_hydrogen_charged_step(molecular,[12d0],captures,charge_old,1d11,300d0, &
           charge_next,ne,gas_h,heat,rate,status,gas_h2=molecules)
      if(status/=0.or.molecules>=h2_initial.or.molecules<0)stop 606
      total_after=sum(matmul(levels,sum(charge_next,dim=3))) &
           +dot_product(bond,sum(sum(charge_next,dim=3),dim=1))+sum(charge_next(:,:,2))*7.02d0*ev
      energy=total_after-total_before+heat+1d11*sum(rate)+(h2_initial-molecules)*pah_h2_binding
      nuclear_error=gas_h-initial_h+2*(molecules-h2_initial)+ &
           dot_product(axis,sum(sum(charge_next-charge_old,dim=3),dim=1))
      if(abs(energy)>1d-10*abs(total_before).or.abs(nuclear_error)>1d-13)stop 607
      if(abs(sum(charge_next(:,:,2))-sum(charge_old(:,:,2))-ne+initial_ne)>1d-13)stop 608
      write(*,'(A,3ES22.14)')'PAH_H2_PHYSICAL energy_relative,H_nuclei,H2_consumed=', &
           energy/abs(total_before),nuclear_error,h2_initial-molecules
      zero_result=charge_next;keep_ne=ne;keep_h=gas_h;keep_heat=heat;keep_rate=rate;keep_h2=molecules
      call pah_hydrogen_charged_step(molecular,[20d0],reshape([1d-5,1d-5],[1,2]),charge_old,1d11,300d0, &
           charge_next,ne,gas_h,heat,rate,status,gas_h2=molecules)
      if(status==0.or.any(charge_next/=zero_result).or.ne/=keep_ne.or.gas_h/=keep_h.or.heat/=keep_heat)stop 609
      if(any(rate/=keep_rate).or.molecules/=keep_h2)stop 610
      call pah_hydrogen_charged_step(molecular,[12d0],captures,charge_old,1d11,300d0, &
           charge_next,ne,gas_h,heat,rate,status)
      if(status==0.or.any(charge_next/=zero_result).or.ne/=keep_ne.or.gas_h/=keep_h)stop 611
      ! Fail AFTER charge/atomic staging: an H2 daughter would exceed the
      ! excitation grid. Every caller-owned inventory must remain untouched.
      charge_old=0;charge_old(size(levels),10,2)=1d-3
      gas_h=0;keep_h=gas_h;ne=0;keep_ne=ne;molecules=.3d0;keep_h2=molecules
      call pah_hydrogen_charged_step(molecular,[12d0],captures,charge_old,1d11,300d0, &
           charge_next,ne,gas_h,heat,rate,status,gas_h2=molecules)
      if(status==0.or.any(charge_next/=zero_result).or.ne/=keep_ne.or.gas_h/=keep_h.or.heat/=keep_heat)stop 619
      if(any(rate/=keep_rate).or.molecules/=keep_h2)stop 620
      write(*,*)'PAH_H2_ZERO_DONOR_PARITY_PHYSICAL_CHARGE_ENERGY_IDENTITY_ROLLBACK_PASS'
    end block
    call check_pah_h2_donor(levels)
  end subroutine

  subroutine check_pah_h2_donor(level)
    use dust_pah_hydrogen
    real(real64),intent(in)::level(:)
    real(real64)::old(size(level),0:13),next(size(level),0:13),saved(size(level),0:13)
    real(real64)::bond(0:13),coef(0:13),atomic(0:13),donor,heat,db,root,events,kinetic,balance
    real(real64)::reference(size(level),0:13),errors(3),h,h2,step,previous,unused
    integer::status,k,j,nsteps
    call pah_hydrogen_parameters(1,bond,atomic,status)
    old=0;old(1,10)=1;next=old;coef=0;coef(10)=1;donor=.25d0;kinetic=1.5d0*1.380649d-16*300
    call pah_hydrogen_attach(level,bond,coef,old,1d0,kinetic,next,donor,heat,db,status, &
         capture_stride=2,donor_binding=pah_h2_binding)
    root=.5d0/(1.75d0+sqrt(1.75d0**2+1d0));events=.25d0-root
    if(status/=0.or.abs(donor-root)>1d-13.or.abs(sum(next(:,12))-events)>1d-13)stop 612
    if(abs(sum(next)-1)>1d-13.or.any(next<0))stop 613
    balance=sum(matmul(level,next-old))+heat+db+events*pah_h2_binding
    if(abs(balance)>1d-12*abs(db))stop 614
    saved=next;old=0;old(size(level),10)=1;donor=.25d0;heat=7;db=9
    call pah_hydrogen_attach(level,bond,coef,old,1d0,kinetic,next,donor,heat,db,status, &
         capture_stride=2,donor_binding=pah_h2_binding)
    if(status==0.or.any(next/=saved).or.donor/=.25d0.or.heat/=7.or.db/=9)stop 615
    coef=0;coef(0:10)=5d-13
    do k=0,3
       nsteps=128
       if(k>0)nsteps=2**(k+1)
       step=1d11/nsteps;old=0;old(1,10)=1;h=.5d0;h2=.3d0
       do j=1,nsteps
          next=old
          call pah_hydrogen_attach(level,bond,atomic,old,step,kinetic,next,h,heat,db,status)
          if(status/=0)stop 616
          old=next
          call pah_hydrogen_attach(level,bond,coef,old,step,kinetic,next,h2,heat,db,status, &
               capture_stride=2,donor_binding=pah_h2_binding)
          if(status/=0)stop 617
          old=next
       enddo
       if(k==0)then
          reference=next
       else
          errors(k)=sum(abs(next-reference))
       endif
    enddo
    if(any(errors(2:)>=errors(:2)))stop 618
    write(*,'(A,3ES22.14)')'PAH_H_H2_SPLIT_REFINEMENT_L1_4_8_16=',errors
    write(*,*)'PAH_H2_ANALYTIC_FINITE_DONOR_ENERGY_OVERFLOW_REFINEMENT_PASS'
  end subroutine

  subroutine check_pah_charged_receiver(neutral_model,levels,neutral_ev,neutral_w,neutral_c)
    type(pah_radiative_model),intent(in)::neutral_model
    real(real64),intent(in)::levels(:),neutral_ev(:),neutral_w(:),neutral_c(:)
    type(pah_radiative_model)::models(2)
    real(real64),allocatable::e(:),w(:),cs(:),rate(:),saved_rate(:)
    real(real64)::pop(size(levels),2),next(size(levels),2),saved(size(levels),2),captures(1,2)
    real(real64)::ne,next_ne,heat,ion,balance,absorbed,keep_ne,keep_heat,keep_ion
    real(real64),parameter::ev=1.602176634d-12
    character(len=2048)::path
    integer::status
    call get_environment_variable('SNRT_PAH_ION_TABLE',path,status=status)
    if(status/=0.or.len_trim(path)==0)return
    call pah_neutral_optics(trim(path),24,e,w,cs,status)
    if(status==0)stop 325
    call pah_charge_optics(trim(path),24,1,e,w,cs,status)
    if(status/=0)stop 326
    if(size(e)/=size(neutral_ev))stop 327
    if(any(e/=neutral_ev).or.any(w/=neutral_w).or.maxval(abs(cs/neutral_c-1))<.1d0)stop 328
    models(1)=neutral_model
    call pah_radiative_prepare(models(2),24,12,levels,e,w,cs,status)
    if(status/=0)stop 329
    allocate(rate(size(e)));rate=-1
    pop=0;pop(1,:)=1d-8;next=-1;ne=.01d0;next_ne=-1;heat=-1;ion=-1
    captures=1d-11;absorbed=sum(captures)*10*ev
    call pah_charged_absorbed_step(models,[10d0],captures,pop,ne,300d0,1d0, &
         next,next_ne,heat,ion,rate,status)
    if(status/=0.or.any(next<0).or.sum(rate)<=0.or.heat<=0)stop 330
    balance=dot_product(levels,sum(next-pop,dim=2))+ion+heat+sum(rate)-absorbed
    if(abs(balance/absorbed)>1d-10.or.abs(sum(next)/sum(pop)-1)>1d-12)stop 331
    if(abs((next_ne-ne)-(sum(next(:,2))-sum(pop(:,2))))>1d-17)stop 332
    if(abs(ion-7.02d0*ev*(sum(next(:,2))-sum(pop(:,2))))>1d-34)stop 333
    saved=next;keep_ne=next_ne;keep_heat=heat;keep_ion=ion;saved_rate=rate
    call pah_charged_absorbed_step(models,[20d0],captures,pop,ne,300d0,1d0, &
         next,next_ne,heat,ion,rate,status)
    if(status==0.or.any(next/=saved).or.next_ne/=keep_ne.or.heat/=keep_heat.or.ion/=keep_ion)stop 334
    if(any(rate/=saved_rate))stop 335
    ! Late failure after photoionization and IR evolution must also roll back.
    call pah_charged_absorbed_step(models,[10d0],captures,pop,ne,1d6,1d0, &
         next,next_ne,heat,ion,rate,status)
    if(status==0.or.any(next/=saved).or.next_ne/=keep_ne.or.heat/=keep_heat.or.ion/=keep_ion)stop 336
    if(any(rate/=saved_rate))stop 337
    write(*,'(A,4ES22.14)')'PAH_CHARGED_LOCAL relative_balance,gas_heat,IP_delta,IR=', &
         balance/absorbed,keep_heat,keep_ion,sum(saved_rate)
    write(*,*)'PAH_CHARGED_ORIGINAL_OPTICS_IR_PHOTOELECTRON_RECOMBINATION_TRANSACTION_PASS'
  end subroutine

  subroutine check_stochastic_time()
    real(real64)::e(3),a(3,3),d(3),p(3),old(3),saved(3),pa,pe,saved_a,saved_e,dt,res(3),truth
    integer::status,k,j,n
    e=[0d0,1d-12,3d-12];a=0;a(2,1)=2;d=[0d0,3d0,5d0]
    old=[1d0,0d0,0d0];p=-1;pa=-1;pe=-1
    call dust_stochastic_evolve(e,a,d,old,.1d0,p,pa,pe,status)
    if(status/=0.or.abs(p(2)-(2d0*.1d0)/(1+5d0*.1d0))>1d-14.or.p(3)/=0)stop 230
    if(abs(dot_product(e,p-old)-pa+pe)>1d-26)stop 231
    saved=p;saved_a=pa;saved_e=pe
    call dust_stochastic_evolve(e,a,d,old*1d-30,.1d0,p,pa,pe,status)
    if(status/=0.or.maxval(abs(p*1d30-saved))>1d-14)stop 239
    if(abs(pa*1d30-saved_a)>1d-26.or.abs(pe*1d30-saved_e)>1d-26)stop 240
    call dust_stochastic_evolve(e,a,d,old*0d0,.1d0,p,pa,pe,status)
    if(status/=0.or.any(p/=0).or.pa/=0.or.pe/=0)stop 241
    truth=.4d0*(1-exp(-5d0))
    do k=1,3
       n=2**(k+2);dt=1d0/n;p=old
       do j=1,n
          saved=p
          call dust_stochastic_evolve(e,a,d,saved,dt,p,pa,pe,status)
          if(status/=0)stop 232
       enddo
       res(k)=abs(p(2)-truth)
    enddo
    if(res(2)>=res(1).or.res(3)>=res(2))stop 233
    a(3,1)=1;a(3,2)=4
    call dust_stochastic_evolve(e,a,d,old,.7d0,p,pa,pe,status)
    if(status/=0)stop 234
    do k=1,3
       truth=sum(a(k,:)*p)-sum(a(:,k))*p(k)-d(k)*p(k)
       if(k<3)truth=truth+d(k+1)*p(k+1)
       if(abs(p(k)-old(k)-.7d0*truth)>1d-14)stop 235
    enddo
    ! At dt*rate >> 1, naive diagonal subtraction loses the identity term.
    call dust_stochastic_evolve(e,a*1d200,d*1d200,old,1d0,p,pa,pe,status)
    if(status/=0.or.maxval(abs(p-1d0/3))>1d-14)stop 236
    saved=p;saved_a=pa;saved_e=pe;a(1,2)=1
    call dust_stochastic_evolve(e,a,d,old,1d0,p,pa,pe,status)
    if(status==0.or.any(p/=saved).or.pa/=saved_a.or.pe/=saved_e)stop 237
    a=0;d=0
    call dust_stochastic_evolve(e,a,d,old,1d200,p,pa,pe,status)
    if(status/=0.or.any(p/=old).or.pa/=0.or.pe/=0)stop 238
    write(*,'(A,3ES14.6)')'DUST_STOCHASTIC_TIME_REFINEMENT errors=',res
    write(*,*)'DUST_STOCHASTIC_FINITE_TIME_STIFF_PROBABILITY_ENERGY_ROLLBACK_PASS'
  end subroutine

  subroutine check_sublimation()
    real(real64)::nodes(65),uc(65),us(65),old(4),new(4),energy,ed,heat,latent,total,td
    real(real64)::curve(65),rate,loss,thermal,phase,eg,dc,ls,gas_elements(11),total_elements(11)
    integer::i,status
    do i=1,65
       nodes(i)=10d0*(300d0)**(real(i-1,real64)/64)
    enddo
    nodes(65)=3000d0
    call dust_composition_curve(nodes,[1d0,0d0],1d0,uc,status)
    if(status/=0)stop 101
    call dust_composition_curve(nodes,[0d0,1d0],1d0,us,status)
    if(status/=0)stop 102
    dust_size_radius_cm=[1d-6,1d-5];dust_size_density=[2.2d0,3.8d0]
    old=[.2d0,.3d0,.1d0,.4d0];energy=sum(old(1:2))*uc(65)+sum(old(3:4))*us(65)
    call dust_sublimation_step(nodes,10d0,old,energy,1d3,new,ed,heat,latent,status)
    if(status/=0.or.any(new<0).or.any(new>old).or.any(new(3:4)/=old(3:4)))stop 103
    if(sum(new(1:2))>=sum(old(1:2)).or.min(heat,latent)<=0)stop 104
    if(abs(ed+heat+latent-energy)>1d-12*energy)stop 105
    if(abs(latent-dust_carbon_latent*sum(old(1:2)-new(1:2)))>1d-13*latent)stop 106
    if(new(1)/old(1)>=new(2)/old(2))stop 107
    curve=sum(new(1:2))*uc+sum(new(3:4))*us
    call snrt_dust_material_temperature(nodes,curve,ed,td,status)
    if(status/=0.or.td<10.or.td>=3000)stop 108
    if(abs(heat-2*dust_kb*td/dust_carbon_atom*sum(old(1:2)-new(1:2)))>1d-11*heat)stop 109
    write(*,'(A,4ES17.8)')'GRAPHITE_SUBLIMATION lost/T/heat/phase=',sum(old-new),td,heat,latent
    call dust_sublimation_step(nodes,10d0,old,energy,1d-8,new,ed,heat,latent,status)
    rate=3*dust_carbon_nu*(dust_carbon_atom/dust_size_density(1))**(1d0/3)* &
         exp(-dust_carbon_binding_k/3000d0)/dust_size_radius_cm(1)
    loss=(old(1)-new(1))/old(1)
    if(status/=0.or.abs(loss/(1d-8*rate)-1)>1d-4)stop 110
    call dust_sublimation_step(nodes,10d0,old,energy,0d0,new,ed,heat,latent,status)
    if(status/=0.or.any(new/=old).or.ed/=energy.or.heat/=0.or.latent/=0)stop 111
    call dust_sublimation_step(nodes,10d0,old,energy,-1d0,new,ed,heat,latent,status)
    if(status==0.or.any(new/=old).or.ed/=energy)stop 112
    ! Common binding reference for any graphite mass-exchange operator:
    ! gas delta=-dust sensible delta+L*delta_solid.
    do i=-1,1,2
       dc=i*.01d0;eg=1d13;phase=dust_carbon_latent*(1-sum(old(1:2)))
       thermal=dc*uc(65)
       total=eg+energy+phase
       if(abs((eg-thermal+dust_carbon_latent*dc)+(energy+thermal)+ &
            (phase-dust_carbon_latent*dc)-total)>1d-14*total)stop 113
    enddo
    write(*,*)'GRAPHITE_SUBLIMATION_RATE_NUCLEI_THERMAL_LATENT_DOMAIN_ROLLBACK_PASS'
    dust_mass_enabled=.true.;dust_mass_model='carbon_olivine_2size_v1'
    dust_material_model='dl01_composition_v1';dust_sublimation='gd89_xu25_olivine_v1'
    call dust_olivine_phase_reference(ls,status)
    if(status/=0.or.ls<2d11.or.ls>3d11)stop 114
    old=[0d0,0d0,.4d0,.6d0];energy=us(65)
    call dust_sublimation_step(nodes,10d0,old,energy,1d3,new,ed,heat,latent,status)
    if(status/=0.or.any(new<0).or.any(new>old).or.sum(new)>=1)stop 115
    if(abs(ed+heat+latent-energy)>1d-12*energy)stop 116
    if(abs(latent-ls*sum(old-new))>1d-13*latent)stop 117
    total_elements=olivine_fraction
    call dust_gas_elements(total_elements,[0d0,sum(new)],gas_elements,status)
    if(status/=0.or.maxval(abs(gas_elements-olivine_fraction*sum(old-new)))>1d-14)stop 118
    if(new(3)/old(3)>=new(4)/old(4))stop 119
    write(*,'(A,3ES17.8)')'OLIVINE_SUBLIMATION loss/latent_specific/energy_balance=',sum(old-new),ls, &
         (ed+heat+latent-energy)/energy
    dust_sublimation='none';dust_material_model='fixed_mix';dust_mass_model='bulk_v1'
  end subroutine

  real(real64) function power_integral(lo,hi,p) result(v)
    real(real64),intent(in)::lo,hi,p
    v=(hi**(p+1)-lo**(p+1))/(p+1)
  end function

  subroutine check_collision_reference()
    real(real64),allocatable::r(:),m(:),k(:,:),birth(:,:,:),leak(:,:,:),old(:),new(:),fine(:),edges(:)
    real(real64)::out(2),outfine(2),density,nh,temp,mach,pcrit,rho,dt,amin,amax,lo,hi,total
    real(real64)::area0,area,areafine,small0,small,smallfine,pair(2),twoarea,exactarea,exactsmall
    real(real64)::time_error(3),areas(3),smalls(3),boundary_mass(3),analytic_error(2),start,finish
    real(real64)::fm(4),fn(4),fb(2),expected,save_density(2),closure,maxclosure,vone,angular_error
    integer::n,i,j,res,phase,material,status,refine
    logical::saved_coag,saved_shatter
    saved_coag=dust_coagulation;saved_shatter=dust_shattering;save_density=dust_size_density
    dust_coagulation=.true.;dust_shattering=.true.;dust_size_density=[2.2d0,3.3d0]
    ! Independently specified one-event ejecta and boundary integral.
    fm=[.001d0,.01d0,.1d0,1d0];fn=0;fb=0
    call dust_collision_fragments(fm,1d0,.5d0,fn,fb)
    expected=.5d0*(.001d0**((4d0-3.3d0)/3)-(1d-8)**((4d0-3.3d0)/3))/ &
         (.01d0**((4d0-3.3d0)/3)-(1d-8)**((4d0-3.3d0)/3))
    if(abs(fb(1)-expected)>1d-14.or.abs(sum(fn*fm)+sum(fb)-1)>1d-14)stop 301
    ! Constant K=1 Smoluchowski: N(t)=N0/(1+K*N0*t/2).
    n=32
    allocate(r(n),m(n),k(n,n),birth(n,n,n),leak(2,n,n),old(n),new(n),fine(n),edges(n+1))
    m=[(real(i,real64),i=1,n)];k=1;birth=0;leak=0
    do j=1,n
       do i=1,j
          call dust_collision_deposit(m,m(i)+m(j),birth(:,i,j),leak(:,i,j))
          birth(:,j,i)=birth(:,i,j);leak(:,j,i)=leak(:,i,j)
       enddo
    enddo
    old=0;old(1)=1
    do refine=1,2
       call dust_collision_evolve(m,k,birth,leak,old,.5d0,.01d0/refine,new,out,status)
       if(status/=0)stop 302
       analytic_error(refine)=abs(sum(new)/.8d0-1)
       if(abs(sum(new*m)+sum(out)-1)>1d-12.or.minval(new)<0)stop 303
    enddo
    if(analytic_error(2)>=analytic_error(1).or.analytic_error(2)>1d-3)stop 304
    write(*,'(A,2ES16.7)')'DUST_COLLISION_CONSTANT_KERNEL_N_ERRORS ',analytic_error
    call dust_collision_evolve(m,k,birth,leak,old,0d0,.02d0,new,out,status)
    if(status/=0.or.any(new/=old).or.any(out/=0))stop 305
    ! Entire merged daughter above the grid: not deposited at last pivot.
    old=0;old(n)=1
    call dust_collision_evolve(m,k,birth,leak,old,.5d0,.005d0,new,out,status)
    if(status/=0.or.out(2)<=0.or.out(1)/=0)stop 314
    if(abs(sum(new*m)+sum(out)-m(n))>1d-12.or.abs(new(n)-2d0/3)>2d-3)stop 315
    ! A corrupted product ledger fails before any state is committed.
    birth(1,1,1)=birth(1,1,1)+.1d0
    call dust_collision_evolve(m,k,birth,leak,old,1d0,.02d0,new,out,status)
    if(status==0.or.any(new/=old).or.any(out/=0))stop 316
    birth(1,1,1)=birth(1,1,1)-.1d0
    k(1,2)=-1
    call dust_collision_evolve(m,k,birth,leak,old,1d0,.02d0,new,out,status)
    if(status==0.or.any(new/=old).or.any(out/=0))stop 306
    k=0;birth=0;leak=0
    call dust_collision_evolve(m,k,birth,leak,old,1d10,.02d0,new,out,status)
    if(status/=0.or.any(new/=old).or.any(out/=0))stop 307
    deallocate(r,m,k,birth,leak,old,new,fine,edges)
    amin=5d-7;amax=2.5d-5
    maxclosure=0
    exactsmall=(sqrt(3d-6)-sqrt(amin))/(sqrt(amax)-sqrt(amin))
    do material=1,2
       density=dust_size_density(material);pcrit=4d10
       if(material==2)pcrit=3d11
       exactarea=3/(4*density)*power_integral(amin,amax,-1.5d0)/power_integral(amin,amax,-.5d0)
       do phase=1,2
          nh=1d3;temp=100;mach=1;dt=3*dust_myr
          if(phase==2)then
             nh=.3d0;temp=8000;mach=3;dt=30*dust_myr
          endif
          rho=1.4d0*dust_mp*nh;total=.01d0*rho
          pair=total*[exactsmall,1-exactsmall]
          call dust_size_exchange(rho,nh,temp,material,dt,pair,status)
          if(status/=0)stop 308
          twoarea=sum(pair*3/(4*density*dust_size_radius_cm))/total
          write(*,'(A,2I3,A,4ES16.7)')'DUST_COLLISION_TWO_SIZE material/phase=',material,phase, &
               ' small0/small/area0/area_cm2_g=',exactsmall,pair(1)/total, &
               sum([exactsmall,1-exactsmall]*3/(4*density*dust_size_radius_cm)),twoarea
          do res=1,3
             n=2**(res+3)
             allocate(r(n),m(n),k(n,n),birth(n,n,n),leak(2,n,n),old(n),new(n),fine(n),edges(n+1))
             do i=1,n
                r(i)=3d-8*(1d-3/3d-8)**(real(i-1,real64)/(n-1))
             enddo
             edges(1)=r(1);edges(n+1)=r(n)
             edges(2:n)=sqrt(r(:n-1)*r(2:n))
             call cpu_time(start)
             call dust_collision_kernel(r,density,nh,temp,mach,pcrit,phase==2,m,k,birth,leak,status)
             if(status/=0)stop 309
             ! Equal-speed isotropic mean relative speed is 4v/3, independent
             ! of the fragment model. Fixed quadrature is accurate to 0.1%.
             vone=1.1d5*mach**1.5d0*sqrt(r(1)/1d-5)*(temp/1d4)**.25d0* &
                  nh**(-.25d0)*sqrt(density/3.5d0)
             angular_error=abs(k(1,1)/(acos(-1d0)*(2*r(1))**2*(4*vone/3))-1)
             if(angular_error>1d-3)stop 317
             old=0
             do i=1,n
                lo=max(edges(i),amin);hi=min(edges(i+1),amax)
                if(hi>lo)old(i)=total*power_integral(lo,hi,-.5d0)/power_integral(amin,amax,-.5d0)/m(i)
             enddo
             area0=sum(old*acos(-1d0)*r*r)/total
             small0=sum(old*m,mask=r<3d-6)/total
             call dust_collision_evolve(m,k,birth,leak,old,dt,.02d0,new,out,status)
             if(status/=0)then
                write(*,*)'COLLISION_EVOLVE_FAILED ',material,phase,n,status;stop 310
             endif
             call dust_collision_evolve(m,k,birth,leak,old,dt,.01d0,fine,outfine,status)
             if(status/=0)stop 311
             call cpu_time(finish)
             if(abs(sum(new*m)+sum(out)-total)>2d-10*total.or.any(new<0))stop 312
             closure=abs(sum(fine*m)+sum(outfine)-total)/total
             maxclosure=max(maxclosure,closure)
             area=sum(new*acos(-1d0)*r*r)/total
             areafine=sum(fine*acos(-1d0)*r*r)/total
             small=sum(new*m,mask=r<3d-6)/total
             smallfine=sum(fine*m,mask=r<3d-6)/total
             time_error(res)=max(abs(area/areafine-1),abs(small-smallfine))
             areas(res)=areafine;smalls(res)=smallfine;boundary_mass(res)=sum(outfine)/total
             write(*,'(A,3I4,A,8ES16.7)')'DUST_COLLISION_MULTI material/phase/n=',material,phase,n, &
                  ' small0/small/area0/area/dt_error/under/over/cpu_s=',small0,smallfine,area0,areafine, &
                  time_error(res),outfine/total,finish-start
             if(time_error(res)>.02d0)stop 313
             deallocate(r,m,k,birth,leak,old,new,fine,edges)
          enddo
          write(*,'(A,2I3,A,5ES16.7)')'DUST_COLLISION_COMPARISON material/phase=',material,phase, &
               ' analytic_initial_area/grid32_64_area/small/two64_area/small=',exactarea, &
               areas(2)/areas(3)-1,smalls(2)-smalls(3),twoarea/areas(3)-1,pair(1)/total-smalls(3)
       enddo
    enddo
    dust_coagulation=saved_coag;dust_shattering=saved_shatter;dust_size_density=save_density
    write(*,'(A,2ES16.7)')'DUST_COLLISION_MAX_MASS_RESIDUAL_AND_ANGLE_ERROR ',maxclosure,angular_error
    write(*,*)'DUST_COLLISION_NATIVE_CONSERVATION_POSITIVITY_REFINEMENT_PASS'
  end subroutine

  subroutine check_multibin_reference()
    real(real64),allocatable::edges(:),numbers(:),slopes(:),nn(:),ss(:)
    real(real64)::amin,amax,shift,factor,scale,lo,m0,raw,gas_delta,corr,m(3),mass,area,removed
    real(real64)::truth_mass,truth_area,truth_n,initial_n,elapsed,start,finish
    real(real64)::mass_error(3),area_error(3),number_error(3),radii(2),fixed_mass,split
    integer::n,k,j,status,it
    ! Same initial truncated MRN and imposed radius loss in every resolution.
    ! No D03 dielectric interpolation or new opacity data is required here.
    amin=1d-7;amax=1d-4;shift=-amin/2;factor=4*acos(-1d0)*3d0/3
    scale=1d-3/(factor*power_integral(amin,amax,-.5d0))
    lo=amin-shift
    truth_n=scale*power_integral(lo,amax,-3.5d0)
    truth_area=scale*(power_integral(lo,amax,-1.5d0)+ &
         2*shift*power_integral(lo,amax,-2.5d0)+shift**2*power_integral(lo,amax,-3.5d0))
    truth_mass=factor*scale*(power_integral(lo,amax,-.5d0)+ &
         3*shift*power_integral(lo,amax,-1.5d0)+3*shift**2*power_integral(lo,amax,-2.5d0)+ &
         shift**3*power_integral(lo,amax,-3.5d0))
    do k=1,3
       n=2**(k+2) ! 8,16,32; 32 checks the direction of convergence.
       allocate(edges(n+1),numbers(n),slopes(n),nn(n),ss(n))
       do j=1,n+1
          edges(j)=amin*(amax/amin)**(real(j-1,real64)/n)
       enddo
       do j=1,n
          raw=scale*power_integral(edges(j),edges(j+1),-3.5d0)
          m0=scale*power_integral(edges(j),edges(j+1),-.5d0)
          call dust_bin_reconstruct(edges(j),edges(j+1),raw,m0,numbers(j),slopes(j),status)
          if(status/=0)stop 60
       enddo
       initial_n=sum(numbers)
       call cpu_time(start)
       do it=1,100
          call dust_multibin_shift(edges,numbers,slopes,shift,3d0,0d0,nn,ss,gas_delta,corr,status,removed)
          if(status/=0)stop 61
       enddo
       call cpu_time(finish);elapsed=(finish-start)/100
       mass=0;area=0
       do j=1,n
          call dust_bin_moments(edges(j),edges(j+1),nn(j),ss(j),edges(j),edges(j+1),0d0,m)
          mass=mass+factor*m(3);area=area+m(2)
          if(nn(j)<0.or.abs(ss(j))*(edges(j+1)-edges(j))/2> &
               nn(j)/(edges(j+1)-edges(j))*(1+32*epsilon(1d0)))stop 62
       enddo
       if(abs(mass+gas_delta-1d-3)>3d-17)stop 63
       if(abs(sum(nn)+removed-corr-initial_n)>1d-13*initial_n)stop 70
       mass_error(k)=abs(mass/truth_mass-1)
       area_error(k)=abs(area/truth_area-1)
       number_error(k)=abs(sum(nn)/truth_n-1)
       write(*,'(A,I3,A,5ES13.5)')'DUST_MULTIBIN_MRN n=',n, &
            ' mass/area/number/limiter/cpu_s=',mass_error(k),area_error(k),number_error(k),corr/initial_n,elapsed
       ! Upper overflow and insufficient gas for growth reject without edits.
       call dust_multibin_shift(edges,numbers,slopes,amin,3d0,1d0,nn,ss,gas_delta,corr,status)
       if(status==0.or.any(nn/=numbers).or.any(ss/=slopes).or.gas_delta/=0)stop 64
       numbers(n)=0;slopes(n)=0
       call dust_multibin_shift(edges,numbers,slopes,amin,3d0,0d0,nn,ss,gas_delta,corr,status)
       if(status==0.or.any(nn/=numbers).or.any(ss/=slopes))stop 65
       call dust_multibin_shift(edges,numbers,slopes,amin,3d0,1d0,nn,ss,gas_delta,corr,status)
       if(status/=0.or.gas_delta>=0)stop 66
       call dust_multibin_shift(edges,numbers,slopes,-2*amax,3d0,0d0,nn,ss,gas_delta,corr,status)
       if(status/=0.or.any(nn/=0).or.any(ss/=0).or.gas_delta<=0)stop 67
       deallocate(edges,numbers,slopes,nn,ss)
    enddo
    if(mass_error(2)>=mass_error(1).or.mass_error(3)>=mass_error(2))stop 68
    ! A signed integral moment can have accidental cancellation at 16 bins;
    ! it need not improve monotonically at every resolution. Both refinements
    ! must improve on 8, and all measured errors remain in the report.
    if(maxval(area_error(2:3))>=area_error(1))stop 69
    ! Native two-size fixed-radius closure, normalized to the same initial
    ! MRN mass. This is a physical closure comparison, not a discretization
    ! of the translated distribution: da changes M at fixed representative a.
    split=3d-6;radii=[5d-7,1d-5]
    fixed_mass=factor*scale*(power_integral(amin,split,-.5d0)*exp(3*shift/radii(1))+ &
         power_integral(split,amax,-.5d0)*exp(3*shift/radii(2)))
    write(*,'(A,2ES14.6)')'DUST_TWO_SIZE_FIXED_VS_RADIUS_SHIFT mass_error,reference_mass=', &
         abs(fixed_mass/truth_mass-1),truth_mass
    write(*,*)'DUST_MULTIBIN_MOMENTS_MASS_GAS_BOUNDARY_ROLLBACK_PASS'
  end subroutine
end program

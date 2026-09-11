program stellar_radioactive_test
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use stellar_enrichment_config
  use stellar_enrichment_contract
  use stellar_yield_tables
  use stellar_yield_audit
  use stellar_ssp_sources
  use stellar_radioactive_decay
  use stellar_radioactive_sources
  use stellar_radioactive_transport
  use stellar_source_increment, only: integrate_ssp_channel_increment,age_radioactive_source
  use stellar_cell_deposition
  use stellar_ramses_field_map
  use stellar_ramses_bridge
  use stellar_ramses_mapped_bridge
  implicit none
  type(stellar_yield_table_t)::table,view
  type(radioactive_companion_t)::companion
  type(stellar_population_t)::population
  type(stellar_cumulative_t)::native
  type(stellar_source_t)::source,before,part
  type(stellar_field_map_t)::map
  character(len=1024)::yields,history,path
  real(stellar_dp)::p(2),s(2),d(2),s1(2),d1(2),s2(2),d2(2),aged(2),lost(2),dt
  real(stellar_dp)::x(11),x0(11),p0(2),total(2),t0,t1,tm,fresh(2),saved
  integer::status,j,n,c,r
  real(stellar_dp)::lo(5),hi(5),delta(21),row(21,1),row0(21,1),rho(1),energy(1),mom(3,1)
  real(stellar_dp)::elem(11,1),parent(2,1),metal(1),center(14),faces(14,2),faces0(14,2),g(16)
  call get_command_argument(1,yields);call get_command_argument(2,history);call get_command_argument(3,path)
  ! The fixture consumes the existing combined Kroupa LC18 prompt selection.
  call set_enrichment_defaults()
  default_imf_id=1;population_model_id=1;configured_binary_fraction=.5d0
  enable_agb=.true.;enable_snia=.true.;high_mass_model='wind_only_collapse'
  configured_channel_mass_min(1)=13;configured_channel_mass_min(3)=13
  configured_channel_mass_max(2)=6
  call load_yield_table(trim(yields),table,status)
  if(status/=0)stop 1
  call prepare_high_mass_history(table,trim(history),status)
  if(status/=0)stop 2
  call set_yield_mass_assignment_mode(table,yield_mass_assignment_piecewise_constant,status)
  call load_radioactive_companion(trim(path),table,companion,status)
  if(status/=0.or..not.companion%loaded)stop 3

  p=[.02d0,.03d0]
  do j=1,2
     call radioactive_decay(p,radioactive_half_life_s(j),s,d,status)
     if(status/=0.or.abs(s(j)/p(j)-.5d0)>2d-15.or.maxval(abs(s+d-p))>1d-17)stop 4
  enddo
  call radioactive_decay(p,0d0,s,d,status)
  if(status/=0.or.any(s/=p).or.any(d/=0))stop 5
  call radioactive_decay(p,huge(1d0)/2,s,d,status)
  if(status/=0.or.any(s/=0).or.any(d/=p))stop 6
  call radioactive_release(p,1d-6,0d0,s,d,status)
  if(status/=0.or.any(d<=0).or.any(s<0))stop 7
  call radioactive_release(p,2*radioactive_half_life_s(1),radioactive_half_life_s(1),s,d,status)
  if(status/=0.or.abs(s(1)/p(1)-.5d0*.75d0/(2*log(2d0)))>2d-15)stop 8
  x=0;x(1)=.7d0;x(2)=.2d0;x(7)=.1d0;x(11)=.05d0;x0=x;p0=p
  call radioactive_gas_decay(.3d0,x,p,radioactive_half_life_s(1),status)
  if(status/=0.or.x(7)<=x0(7).or.x(11)>=x0(11))stop 9
  if(abs(x(7)-x0(7)-(p0(1)-p(1)))>2d-17)stop 10
  if(abs(x(11)-x0(11)+(p0(2)-p(2)))>2d-17)stop 11
  x=x0;p=p0
  call radioactive_gas_decay(.01d0,x,p,1d12,status)
  if(status==0.or.any(x/=x0).or.any(p/=p0))stop 12
  call radioactive_gas_decay(.3d0,x,p,ieee_value(0d0,ieee_quiet_nan),status)
  if(status==0.or.any(x/=x0).or.any(p/=p0))stop 13
  print *, 'RADIOACTIVE_ANALYTIC_POSITIVE_SMALL_LARGE_ZERO_GAS_MG_FE_ATOMIC_PASS'

  t0=0;t1=.05d0;tm=.004d0
  do n=1,size(table%hm_mass)
     do c=1,3,2
        call radioactive_node_interval(table,companion,n,c,t0,t1,s,d,status)
        if(status/=0.or.any(s<0).or.any(d<0))stop 14
        r=table%hm_wind_row(n);if(c==3)r=table%hm_terminal_row(n)
        total=companion%cumulative(:,r)
        if(maxval(abs(s+d-total))>2d-13*max(maxval(total),tiny(1d0)))stop 15
        call radioactive_node_interval(table,companion,n,c,t0,tm,s1,d1,status)
        if(status/=0)stop 16
        call radioactive_node_interval(table,companion,n,c,tm,t1,s2,d2,status)
        if(status/=0)stop 17
        call radioactive_decay(s1,(t1-tm)*1d9*radioactive_year_s,aged,lost,status)
        if(status/=0)stop 18
        if(maxval(abs(s-aged-s2))>2d-13*max(maxval(total),tiny(1d0)))stop 19
        if(maxval(abs(d-d1-lost-d2))>2d-13*max(maxval(total),tiny(1d0)))stop 20
        if(c==3)then
           call radioactive_node_interval(table,companion,n,c,0d0,table%hm_age(n),s,d,status)
           if(status/=0.or.any(s/=total).or.any(d/=0))stop 21
           call radioactive_node_interval(table,companion,n,c,table%hm_age(n),t1,s,d,status)
           if(status/=0.or.any(s/=0).or.any(d/=0))stop 22
        endif
     enddo
  enddo
  print *, 'RADIOACTIVE_REAL_LC18_36_NODES_WIND_KNOT_TERMINAL_CONVOLUTION_SPLIT_PASS'

  population%initial_mass=10000;population%current_mass=10000;population%imf_id=1
  population%birth_metallicity=.005d0;population%imf_mass_min=.08d0;population%imf_mass_max=120
  population%yield_basis_id=yield_basis_per_star_cumulative
  view=table
  ! Independent existing native integrator with the two parent payloads in
  ! scratch element columns: compare the SAME source-cell IMF normalization.
  view%ejected_mass(:,7)=companion%cumulative(1,:)
  view%ejected_mass(:,11)=companion%cumulative(2,:)
  do c=1,3,2
     call integrate_radioactive_channel(table,companion,population,c,t0,t1,13d0,120d0,64,s,d,status)
     if(status/=0)stop 23
     call integrate_ssp_channel(view,population,c,t1,13d0,120d0,64,native,status)
     if(status/=0)stop 24
     fresh=[native%ejected_mass(7),native%ejected_mass(11)]
     if(maxval(abs(s+d-fresh))>2d-12*maxval(fresh))stop 25
     if(any(s>=fresh).or.any(d<=0))stop 26 ! Not fresh injection at step end.
     call integrate_radioactive_channel(table,companion,population,c,t0,tm,13d0,120d0,64,s1,d1,status)
     if(status/=0)stop 27
     call integrate_radioactive_channel(table,companion,population,c,tm,t1,13d0,120d0,64,s2,d2,status)
     if(status/=0)stop 28
     call radioactive_decay(s1,(t1-tm)*1d9*radioactive_year_s,aged,lost,status)
     if(maxval(abs(s-aged-s2))>2d-12*maxval(fresh))stop 29
     if(maxval(abs(d-d1-lost-d2))>2d-12*maxval(fresh))stop 30
     write(*,'(A,I2,6ES18.9)')'RADIOACTIVE_IMF_CHANNEL_FRESH_SURVIVING_DECAYED ',c,fresh,s,d
  enddo
  view=table;view%ejected_mass(1,7)=view%ejected_mass(1,7)+1d-10
  saved=sum(companion%cumulative)
  call load_radioactive_companion(trim(path),view,companion,status)
  if(status==0.or..not.companion%loaded.or.sum(companion%cumulative)/=saved)stop 31
  view=table;view%high_mass_identity(2)='PARSEC'
  call load_radioactive_companion(trim(path),view,companion,status)
  if(status==0.or.sum(companion%cumulative)/=saved)stop 32
  print *, 'RADIOACTIVE_SHARED_IMF_LINEAR_Z_SOURCE_MISMATCH_ATOMIC_PASS'

  call configure_radioactive_source(table,trim(path),status)
  if(status/=0)stop 33
  call clear_source(source)
  lo=13;hi=120
  do c=1,3,2
     call integrate_ssp_channel_increment(table,population,c,t0,tm,lo(c),hi(c),64,part,status)
     if(status/=0)stop 34
     source%returned_mass=source%returned_mass+part%returned_mass
     source%ejected_mass=source%ejected_mass+part%ejected_mass
     source%channel_returned_mass(c)=part%returned_mass
     source%channel_ejected_mass(c,:)=part%ejected_mass
  enddo
  before=source
  call age_radioactive_source(table,population,t0,tm,lo,hi,64,source,status)
  if(status/=0.or.any(source%radioactive_parent<=0))stop 35
  if(source%returned_mass/=before%returned_mass.or.source%energy/=before%energy.or. &
       any(source%momentum/=before%momentum).or.any(source%net_yield/=before%net_yield))stop 36
  if(source%ejected_mass(7)<=before%ejected_mass(7).or.source%ejected_mass(11)>=before%ejected_mass(11))stop 37
  if(.not.radioactive_source_valid(source))stop 38
  before=source
  call age_radioactive_source(table,population,t0,tm,lo,hi,64,source,status)
  if(status==0.or.any(source%ejected_mass/=before%ejected_mass).or. &
       any(source%radioactive_parent/=before%radioactive_parent))stop 39
  map%density_index=1;map%energy_index=5;map%momentum_index=[2,3,4];map%total_metal_index=8
  map%element_index=[(j,j=9,19)];map%radioactive_index=[20,21]
  call build_stellar_source_unew_delta(source,[0d0,0d0,0d0],2d0,1d0,1d0,3d0,21,3,map,delta,1d-12,status)
  if(status/=0.or.any(delta(20:21)/=source%radioactive_parent/2d0/3d0))stop 40
  if(abs(delta(1)-source%returned_mass/6)>1d-14)stop 41
  row=0
  call deposit_source_to_uold_mapped(row,21,1,3,[6d0],[1d0],source,map,status)
  if(status/=0.or.maxval(abs(row(:,1)-delta))>1d-13)stop 42
  row0=row
  map%radioactive_index(1)=19
  call deposit_source_to_uold_mapped(row,21,1,3,[6d0],[1d0],source,map,status)
  if(status==0.or.any(row/=row0))stop 43
  map%radioactive_index=[20,21];row=0
  call deposit_source_to_uold(source,21,1,[6d0],[1d0],1,5,[2,3,4],map%element_index,row, &
       1d-12,status,8,[20,21])
  if(status/=0.or.maxval(abs(row(:,1)-delta))>1d-13)stop 44
  rho=0;elem=0;energy=0;mom=0;metal=0;parent=0
  call deposit_stellar_source(source,1,[6d0],[1d0],rho,elem,energy,mom,1d-12,status,metal,parent)
  if(status/=0.or.maxval(abs(parent(:,1)-delta(20:21)))>1d-17)stop 45
  if(abs(metal(1)-delta(8))>1d-14.or.maxval(abs(elem(:,1)-delta(9:19)))>1d-13)stop 46
  call deposit_stellar_source(source,1,[6d0],[1d0],rho,elem,energy,mom,1d-12,status,metal)
  if(status==0.or.maxval(abs(elem(:,1)-delta(9:19)))>1d-13)stop 47
  call configure_radioactive_source(table,'',status)
  before=source
  call age_radioactive_source(table,population,t0,tm,lo,hi,64,source,status)
  if(status/=0.or.any(source%ejected_mass/=before%ejected_mass))stop 48
  print *, 'RADIOACTIVE_ACTUAL_SOURCE_AGING_CELL_MAPPED_CODE_UNITS_NO_EXTRA_MASS_ATOMIC_PASS'

  center=0;center(1)=.3d0;center(2)=.6d0;center(3)=.1d0
  center(8)=.04d0;center(12)=.06d0;center(13)=.03d0;center(14)=.02d0
  faces(:,1)=center;faces(:,2)=center
  faces(12,1)=.01d0;faces(12,2)=.11d0
  faces(13,1)=.25d0;faces(13,2)=-.19d0
  faces0=faces
  call radioactive_limit_states(center,faces,status)
  if(status/=0.or.all(faces==faces0))stop 49
  do j=1,2
     g=radioactive_constraints(faces(:,j))
     if(minval(g)<-1d-16)stop 50
  enddo
  if(maxval(abs(sum(faces,dim=2)/2-center))>1d-16)stop 51
  ! Same theta applies to host Fe as to both subsets, not a posthoc clip.
  if(abs((faces(12,1)-center(12))/(faces0(12,1)-center(12))- &
       (faces(13,1)-center(13))/(faces0(13,1)-center(13)))>2d-15)stop 52
  faces0=faces;center(14)=.2d0
  call radioactive_limit_states(center,faces,status)
  if(status==0.or.any(faces/=faces0))stop 53
  print *, 'RADIOACTIVE_COMMON_FACE_CHILD_LIMITER_HOST_SUBSETS_MEAN_ATOMIC_PASS'
end program

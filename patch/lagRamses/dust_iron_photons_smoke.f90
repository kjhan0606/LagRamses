program dust_iron_photons_smoke
  use iso_fortran_env, only: real64
  use dust_iron_photons
  use dust_mass_physics
  use dust_iron_compare, only: iron_compare_source_receipt,iron_compare_weights,iron_compare_neutral_area
  use dust_iron_optics, only: fe_primary_ev
  use snrt_dust_receiver, only: snrt_dust_receiver_stage
#ifdef SNRT_CHIMES
  use snrt_chimes, only: chimes_secondary_partition
  use snrt_thermochemistry, only: snrt_secondary_tables_load_from_environment
#endif
  implicit none
  real(real64),parameter::ev=1.602176634d-12,kb=1.380649d-16,mp=1.67262192369d-24
  real(real64)::y,pd,lo,hi,ip,k,e,k0,e0,u,old,primary,weights(6),area,area0
  real(real64)::bins(6),number(9),state(157),gas,solid,ledger(8),n0(9),s0(157),g0,d0,l0(8)
  real(real64)::identity(dust_fe_photon_identity_n),data(fe_photon_data_identity_n)
  integer::status
#ifdef SNRT_CHIMES
  call snrt_secondary_tables_load_from_environment(status)
  call check(status==0,'pinned secondary table load')
#endif
  block
    real(real64)::next_u(1),next_t(1),absorbed(1)
    ! Same live receiver used by the thermal limit, including the hard edge.
    call snrt_dust_receiver_stage(reshape([1d-10],[1,1]),[1d4],1d0,[1d0],[1d-10], &
         [1d-18],[20d0],next_u,next_t,absorbed,status,defer_temperature=.true.)
    call check(status==0.and.abs(absorbed(1)-1d-6*ev)<1d-30.and. &
         abs(next_u(1)-1d-18-absorbed(1))<1d-30,'hard photon full retention receiver')
  end block
  dust_mass_model='carbon_olivine_2size_v1';dust_material_model='dl01_composition_v1'
  dust_optics_model='d03_transport_v1';dust_cooling='chimes_neq_v1'
  dust_size_radius_cm=[1d-6,1d-5];dust_size_density=[2.2d0,3.8d0]
  dust_mass_enabled=.true.;dust_iron_model='fe_electric_compare_v1'
  call check(dust_fe_primary_limit()==4d0,'legacy primary limit')
  dust_iron_model='fe_thermal_limit_v1';identity=dust_fe_photon_identity()
  call check(dust_fe_primary_limit()==1d4.and.identity(2)==2,'thermal limit identity')
  dust_iron_model='fe_uv_cycle_v1';identity=dust_fe_photon_identity();data=fe_photon_data_identity()
  call check(dust_fe_primary_limit()==13.6d0.and.identity(2)==3.and.data(1)==2,'UV identity')
  call check(all(data(36:39)==[-110d0,-8604d0,63d0,632d0]).and.all(data(89:92)==fe_gas_ip_ev), &
       'identity binds charge bounds and chemical reference')
  call check(abs(fe_gas_ip_ev(1)-13.598440002498787d0)<1d-13.and. &
       abs(fe_gas_ip_ev(4)-11.260291860855007d0)<1d-13,'shared ATcT H/C binding energies')
  call fe_uv_photo(0,1d-6,4d0,1d-6,y,pd,lo,hi,ip,status)
  call check(status==0.and.y==0.and.pd==0,'below threshold')
  call fe_uv_photo(-3,1d-6,12d0,1d-6,y,pd,lo,hi,ip,status)
  call check(status==0.and.y>0.and.pd>0.and.lo>0.and.hi>lo,'negative PE and detachment')
  call check(abs(fe_charge_energy(-2,1d-6)-fe_charge_energy(-3,1d-6)-ip)<1d-12,'negative F edge')
  call fe_uv_photo(-10,1d-6,fe_primary_ev(4),1d-6,y,pd,lo,hi,ip,status)
  call check(status==0.and.pd>0.and.fe_primary_ev(4)-ip>9d0,'negative detachment exceeds 9 eV')
  call fe_uv_photo(-50,1d-6,fe_primary_ev(4),1d-6,y,pd,lo,hi,ip,status)
  call check(status==0.and.pd>0.and.fe_primary_ev(4)-ip>13.6d0,'negative detachment can ionize H')
  call fe_uv_photo(3,1d-6,12d0,1d-6,y,pd,lo,hi,ip,status)
  call check(status==0.and.y>0.and.pd==0.and.lo<0.and.hi>0,'positive escape distribution')
  call check(abs(fe_charge_energy(4,1d-6)-fe_charge_energy(3,1d-6)-ip)<1d-12,'positive F edge')
  call fe_oml(0,1d-6,100d0,mp,1,1d0,k0,e0,status)
  call check(status==0.and.abs(e0-2*kb*100/ev)<1d-15,'neutral OML moment')
  u=14.3996454784255d0/100/(kb*100/ev)
  call fe_oml(1,1d-6,100d0,mp,1,1d0,k,e,status)
  call check(status==0.and.abs(k/k0-exp(-u))<1d-15.and.abs(e-kb*100/ev*(u+2))<1d-14,'repulsive OML')
  call fe_oml(-1,1d-6,100d0,mp,1,1d0,k,e,status)
  call check(status==0.and.abs(k/k0-(1+u))<1d-12.and.abs(e-kb*100/ev*(2+u)/(1+u))<1d-14,'attractive OML')
  old=2;primary=3
  call iron_compare_source_receipt(-4d0,old,primary,status)
  call check(status==0.and.old==1.and.primary==0,'funded signed solid source')
  call iron_compare_source_receipt(-2d0,old,primary,status)
  call check(status/=0.and.old==1.and.primary==0,'unfunded source rollback')
  bins=[1d-25,1d-25,1d-25,1d-25,1d-29,1d-29]
  call iron_compare_weights(bins,weights,area,1d-26,status);area0=area
  call iron_compare_neutral_area(weights,1d-26,.9d0,area,status)
  call check(status==0.and.area<area0.and.area>0,'Fe-only neutral area')
  state=0;state(2)=.99d0;state(3)=.01d0;state(5)=.079d0
  state(8)=5d-5;state(9)=5d-5;state(1)=state(3)+state(9)
  number=0;number(3:4)=1d-5
  gas=1.5d0*kb*100d0*100d0*sum(state);solid=1d-18;ledger=-1
  n0=number;s0=state;g0=gas;d0=solid;l0=ledger
  call fe_uv_step(1d12,3d8,100d0,20d0,bins,number,state,gas,solid,secondary,ledger,status)
  write(*,*)'UV fixture status / ledger=',status,ledger
  call check(status==0,'finite-time UV actual coupling')
  call check(sum(number)<sum(n0).and.state(9)<s0(9).and.state(1)/=s0(1),'photons/C+/electron change')
  call check(gas/=g0.and.solid/=d0.and.ledger(1)>0.and.ledger(8)>=20,'gas/solid finite-time receipts')
  block
    real(real64)::nr(9,2),sr(157,2),gr(2),dr(2),lr(8,2)
    integer::jr,errors(2)
    nr=spread(n0,2,2);sr=spread(s0,2,2);gr=g0;dr=d0;lr=-1
!$omp parallel do default(shared) private(jr)
    do jr=1,2
       call fe_uv_step(1d12,3d8,100d0,20d0,bins,nr(:,jr),sr(:,jr),gr(jr),dr(jr),secondary,lr(:,jr),errors(jr))
    enddo
!$omp end parallel do
    call check(all(errors==0).and.all(nr(:,1)==number).and.all(sr(:,2)==state).and. &
         all(gr==gas).and.all(dr==solid).and.all(lr(:,1)==ledger).and.all(lr(:,2)==ledger), &
         'concurrent worker scratch matches serial exactly')
  end block
  e=(gas-g0+solid-d0)/ev+100d0*(fe_gas_ip_ev(1)*(state(3)-s0(3))+fe_gas_ip_ev(4)*(state(9)-s0(9)))+ledger(4)
  call check(abs(e-dot_product(n0-number,fe_primary_ev))<2d-10*g0/ev,'radiation gas solid chemical budget')
  number=n0;state=s0;gas=g0;solid=d0;ledger=l0;number(5)=1
  n0=number
  call fe_uv_step(1d12,3d8,100d0,20d0,bins,number,state,gas,solid,secondary,ledger,status)
  call check(status==2,'hard UV rejection');call unchanged()
  number(5)=0;n0=number
  call fe_uv_step(1d-6,3d8,100d0,20d0,bins,number,state,gas,solid,secondary,ledger,status)
  call check(status==6,'slow charge relaxation rejection');call unchanged()
  bins(5:6)=1d-17
  call fe_uv_step(1d12,3d8,100d0,20d0,bins,number,state,gas,solid,secondary,ledger,status)
  call check(status==5,'nontrace inventory rejection');call unchanged()
#ifdef SNRT_CHIMES
  block
    real(real64)::fractions(5)
    call secondary(14.946275201263603d0,state,fractions,status)
    call check(status==0.and.sum(fractions(2:5))>0.and.abs(sum(fractions)-1)<1d-12, &
         'negative-grain emitted electron uses real nonthermal FS partition')
    bins(5:6)=1d-29;number=n0;state=s0;gas=g0*99;solid=d0;ledger=l0
    call fe_uv_step(1d12,3d8,100d0,20d0,bins,number,state,gas,solid,secondary,ledger,status)
    write(*,*)'UV real-FS warm-gas status / ledger=',status,ledger
    call check(status==0.and.ledger(4)>0,'warm negative-grain cycles have actual excitation loss')
  end block
#endif
  print *,'PASS dust_iron_photons_smoke'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    if(.not.ok)then
       print *,'FAIL ',label;stop 1
    endif
  end subroutine
  subroutine unchanged()
    call check(all(number==n0).and.all(state==s0).and.gas==g0.and.solid==d0.and.all(ledger==l0),'atomic rollback')
  end subroutine
  subroutine secondary(energy,s,f,ierr)
    ! CHIMES-enabled target exercises exactly the live target-limited FS
    ! callback. Disabled builds retain an explicit all-heat algebra fixture.
    real(real64),intent(in)::energy,s(157)
    real(real64),intent(out)::f(5)
    integer,intent(out)::ierr
#ifdef SNRT_CHIMES
    f=0;ierr=chimes_secondary_partition(energy,s,f)
#else
    f=[1d0,0d0,0d0,0d0,0d0];ierr=0
    if(energy<0.or.any(s<0))ierr=1
#endif
  end subroutine
end program

program snrt_thermochemistry_smoke
  use amr_parameters, only: dp
  use snrt_atomic_cooling
#ifdef SNRT_CHIMES
  use snrt_chimes
  use iso_c_binding, only: c_null_char,c_ptr,c_null_ptr,c_associated,c_int
#endif
  use dust_element_cooling, only: wss09_metals_rate
  use snrt_thermochemistry, only: &
       snrt_secondary_tables_load_from_environment, snrt_secondary_tables_loaded, &
       snrt_secondary_source_id, snrt_secondary_upstream_commit, &
       snrt_secondary_manifest_sha256, snrt_secondary_fractions, snrt_secondary_fractions_c, &
       snrt_alpha_hydrogen_case_b, snrt_alpha_helium_ii_case_b, &
       snrt_alpha_helium_ii_radiative_case_b, &
       snrt_alpha_helium_ii_dielectronic_case_b, snrt_alpha_helium_iii_case_b, &
       snrt_partition_absorption, &
       snrt_thermochemistry_advance_cell, snrt_thermochemistry_result, &
       snrt_thermochemistry_ok, snrt_thermochemistry_err_inventory
  implicit none

  integer, parameter :: ngroup = 9
  real(dp), parameter :: expected_fion_200 = 0.09607314283559579d0
  ! The public table row has a 9.41e-6 rounding residual; the native
  ! contract normalizes all five deposition channels before returning them.
  real(dp), parameter :: expected_fheat_200 = 0.8050044271114285d0
  real(dp), parameter :: expected_fexc_200 = 0.09892333418982854d0
  real(dp) :: fheat, fhi, fhei, fheii, fexc, fion
  real(dp) :: fheat_low, fhi_low, fhei_low, fheii_low, fexc_low
  real(dp) :: fheat_high, fhi_high, fhei_high, fheii_high, fexc_high
  real(dp) :: max_delta, max_floor_delta
  real(dp) :: alpha_h, alpha_heii, alpha_heiii, alpha_heii_hot
  real(dp) :: alpha_heii_rad_hot, alpha_heii_dielectronic_hot
  real(dp) :: opacity(3), available(3), partition(3)
  real(dp) :: unassigned
  real(dp) :: absorbed(3,ngroup), excess(3,ngroup)
  real(dp) :: simplex
  type(snrt_thermochemistry_result) :: result
  integer :: ierr, failures

  failures = 0
  call snrt_secondary_tables_load_from_environment(ierr)
  call expect(ierr == snrt_thermochemistry_ok .and. snrt_secondary_tables_loaded, &
       'native FS2010 contract and all fourteen tables load', failures)
  call expect(trim(snrt_secondary_source_id) == &
       'furlanetto_stoever_2010_21cmfast', &
       'native table source identity is pinned', failures)
  call expect(len_trim(snrt_secondary_upstream_commit) == 40 .and. &
       len_trim(snrt_secondary_manifest_sha256) == 64, &
       'native table contract exposes upstream and manifest identities', failures)

  call snrt_secondary_fractions(200.0d0, 0.1d0, fheat, fhi, fhei, fheii, &
       fexc, ierr, fion)
  call expect(ierr == 0, 'native FS2010 interpolation returns successfully', failures)
  call expect(abs(fion-expected_fion_200) < 2.0d-12 .and. &
       abs(fheat-expected_fheat_200) < 2.0d-12 .and. &
       abs(fexc-expected_fexc_200) < 2.0d-12, &
       'native 200 eV, xHII=0.1 interpolation matches pinned reference', failures)
  call expect(abs(fheat+fhi+fhei+fheii+fexc-1.0d0) < 2.0d-15 .and. &
       min(fheat,fhi,fhei,fheii,fexc) >= 0.0d0, &
       'native deposition fractions close and remain non-negative', failures)

  call snrt_secondary_fractions(99.9d0, 0.1d0, fheat_low, fhi_low, fhei_low, &
       fheii_low, fexc_low, ierr)
  call snrt_secondary_fractions(100.1d0, 0.1d0, fheat_high, fhi_high, fhei_high, &
       fheii_high, fexc_high, ierr)
  max_delta = maxval(abs((/fheat_high-fheat_low, fhi_high-fhi_low, &
       fhei_high-fhei_low, fheii_high-fheii_low, fexc_high-fexc_low/)))
  call expect(max_delta < 5.0d-3, 'native 99.9/100.1 eV interpolation is continuous', failures)
  call snrt_secondary_fractions(9.999d0, 0.1d0, fheat_low, fhi_low, fhei_low, &
       fheii_low, fexc_low, ierr)
  call snrt_secondary_fractions(10.001d0, 0.1d0, fheat_high, fhi_high, fhei_high, &
       fheii_high, fexc_high, ierr)
  max_floor_delta = maxval(abs((/fheat_high-fheat_low, fhi_high-fhi_low, &
       fhei_high-fhei_low, fheii_high-fheii_low, fexc_high-fexc_low/)))
  call expect(max_floor_delta < 5.0d-3, 'native 10 eV table-floor transition is bounded', failures)
  call snrt_secondary_fractions(5.0d0, 0.1d0, fheat, fhi, fhei, fheii, fexc, ierr)
  call expect(ierr == 0 .and. fheat == 1.0d0 .and. fhi == 0.0d0 .and. &
       fhei == 0.0d0 .and. fheii == 0.0d0 .and. fexc == 0.0d0, &
       'below-table electron energy is assigned entirely to heat', failures)

  alpha_h = snrt_alpha_hydrogen_case_b(10000.0d0)
  alpha_heii = snrt_alpha_helium_ii_case_b(10000.0d0)
  alpha_heiii = snrt_alpha_helium_iii_case_b(10000.0d0)
  call expect(alpha_h > 0.0d0 .and. alpha_heii > 0.0d0 .and. alpha_heiii > 0.0d0, &
       'native case-B recombination coefficients are positive', failures)
  call expect(abs(alpha_heiii-2.0d0*snrt_alpha_hydrogen_case_b(2500.0d0)) < 1.0d-24, &
       'He III uses exactly 2 alpha_H,B(T/4)', failures)
  call expect(abs(alpha_heii-2.616130035d-13)/2.616130035d-13 < 2.0d-6, &
       'He II case-B coefficient matches the temperature-resolved reference', failures)
  alpha_heii_hot = snrt_alpha_helium_ii_case_b(100000.0d0)
  alpha_heii_rad_hot = snrt_alpha_helium_ii_radiative_case_b(100000.0d0)
  alpha_heii_dielectronic_hot = snrt_alpha_helium_ii_dielectronic_case_b(100000.0d0)
  call expect(abs(alpha_heii_hot-alpha_heii_rad_hot-alpha_heii_dielectronic_hot) / &
       alpha_heii_hot < 1.0d-12 .and. alpha_heii_dielectronic_hot > 0.0d0 .and. &
       alpha_heii_dielectronic_hot/alpha_heii_hot > 1.0d-3, &
       'He II case-B retains the non-negligible dielectronic term at 1e5 K', failures)

  opacity = (/1.0d0, 10.0d0, 1.0d0/)
  available = (/0.10d0, 0.02d0, 0.01d0/)
  call snrt_partition_absorption(0.12d0, opacity, available, partition, ierr)
  call expect(ierr == 0 .and. abs(sum(partition)-0.12d0) < 1.0d-12 .and. &
       minval(available) >= -1.0d-14, &
       'native absorption partition redistributes around species inventory caps', failures)
  opacity = (/1.0d0, 0.0d0, 0.0d0/)
  available = (/0.10d0, 0.02d0, 0.01d0/)
  call snrt_partition_absorption(0.10000001d0, opacity, available, partition, ierr, unassigned)
  call expect(ierr == 0 .and. abs(partition(1)-0.10d0) < 1.0d-12 .and. &
       partition(2) == 0.0d0 .and. partition(3) == 0.0d0 .and. &
       unassigned > 0.0d0, &
       'partition redistribution never assigns a group to an opaque-zero species', failures)
  opacity = (/1.0d0, 0.0d0, 0.0d0/)
  available = (/1.0d-10, 0.0d0, 0.0d0/)
  call snrt_partition_absorption(1.00001d-10, opacity, available, partition, ierr, &
       unassigned, inventory_scale_code=1.0d0)
  call expect(ierr == 0 .and. abs(partition(1)-1.0d-10) < 1.0d-24 .and. &
       unassigned > 0.0d0, &
       'partition tolerance remains tied to the pre-partition cell scale', failures)
  available = (/0.10d0, 0.0d0, 0.0d0/)
  call snrt_partition_absorption(0.11d0, opacity, available, partition, ierr, unassigned, &
       inventory_scale_code=0.10d0)
  call expect(ierr == snrt_thermochemistry_err_inventory .and. &
       abs(available(1)-0.10d0) < 1.0d-15, &
       'above-tolerance unassigned absorption is rejected without mutation', failures)

  absorbed = 0.0d0
  excess = 0.0d0
  absorbed(1,7) = 0.05d0
  absorbed(2,7) = 0.002d0
  absorbed(3,7) = 0.001d0
  excess(:,7) = 200.0d0
  call snrt_thermochemistry_advance_cell(1.0d0, 0.0789474d0, 1.0d0, &
       10000.0d0, 1.0d11, 0.10d0, 0.10d0, 0.0d0, absorbed, excess, result)
  call expect(result%ierr == 0, 'native H/He photo-thermochemistry step succeeds', failures)
  simplex = result%x_helium_ii + result%x_helium_iii
  call expect(result%x_hydrogen_ii >= 0.0d0 .and. result%x_hydrogen_ii <= 1.0d0 .and. &
       result%x_helium_ii >= 0.0d0 .and. result%x_helium_iii >= 0.0d0 .and. &
       simplex <= 1.0d0 + 1.0d-12 .and. result%electron_density_cm3 >= 0.0d0, &
       'native H/He fractions remain on their physical simplex', failures)
  call expect(result%primary_hydrogen_ionizations_cm3 > 0.0d0 .and. &
       result%secondary_hydrogen_ionizations_cm3 >= 0.0d0 .and. &
       result%secondary_helium_i_ionizations_cm3 >= 0.0d0 .and. &
       result%secondary_helium_ii_ionizations_cm3 >= 0.0d0 .and. &
       result%recombination_hydrogen_cm3 >= 0.0d0 .and. &
       result%recombination_helium_ii_cm3 >= 0.0d0 .and. &
       result%recombination_helium_iii_cm3 >= 0.0d0, &
       'native primary/secondary/recombination ledgers are non-negative', failures)
  call expect(abs(result%photoelectron_energy_residual_ev_cm3) < &
       1.0d-11*max(1.0d0,result%photoelectron_energy_ev_cm3), &
       'native photoelectron energy closes into heat, ionization, and excitation', failures)
  call expect(result%heating_rate_erg_cm3_s > 0.0d0 .and. &
       result%absorbed_photon_energy_ev_cm3 >= result%photoelectron_energy_ev_cm3, &
       'only the explicitly deposited gas-heating channel is exposed to RAMSES', failures)

  absorbed = 0.0d0
  excess = 0.0d0
  absorbed(2,7) = 0.002d0
  excess(2,7) = 200.0d0
  call snrt_thermochemistry_advance_cell(1.0d0, 0.0789474d0, 1.0d0, &
       10000.0d0, 1.0d11, 1.0d0, 0.0d0, 0.0d0, absorbed, excess, result)
  call expect(result%ierr == 0 .and. result%secondary_hydrogen_ionizations_cm3 == 0.0d0 .and. &
       result%secondary_heating_energy_ev_cm3 > 0.0d0, &
       'secondary ionization unavailable in a saturated H II cell is routed to heat', failures)

  absorbed=0;excess=0
  call snrt_thermochemistry_advance_cell(1d0,.08d0,1d0,1d4,1d12,.5d0,.2d0,.1d0, &
       absorbed,excess,result,defer_recombination=.true.)
  call expect(result%ierr==0.and.result%x_hydrogen_ii==.5d0.and. &
       result%x_helium_ii==.2d0.and.result%x_helium_iii==.1d0.and. &
       result%recombination_hydrogen_cm3==0,'deferred recombination has exactly one downstream owner',failures)
  call check_atomic(failures)
  call check_band_secondary(failures)
#ifdef SNRT_CHIMES
  call check_chimes_band(failures)
  call check_chimes(failures)
  call check_chimes_photo(failures)
  call check_chimes_hot_split(failures)
  call check_chimes_competing_dust(failures)
  call check_chimes_molecular_photo(failures)
#endif
  if (failures == 0) then
     write(*,'(a)') 'SNRT_NATIVE_THERMOCHEMISTRY_OK'
  else
     write(*,'(a,i0)') 'SNRT_NATIVE_THERMOCHEMISTRY_FAIL count=', failures
     error stop 1
  end if

contains

  subroutine check_band_secondary(failures)
    integer,intent(inout)::failures
    real(dp)::a(3,9),e(3,9),d(8),f(5),g(5),initial(8),expected,delta
    type(snrt_thermochemistry_result)::old,new,broad
    integer::status
    a=0;e=0;d=0
    a(1,7)=1d-6;e(1,7)=200
    status=snrt_secondary_fractions_c(200d0,.1d0,f)
    call expect(status==0.and.abs(sum(f)-1)<1d-14,'C callback uses loaded physical FS2010 fractions',failures)
    d(1)=a(1,7);d(4:8)=a(1,7)*200*f
    call snrt_thermochemistry_advance_cell(1d0,.079d0,1d0,1d4,0d0,.1d0,.1d0,.1d0,a,e,old)
    call snrt_thermochemistry_advance_cell(1d0,.079d0,1d0,1d4,0d0,.1d0,.1d0,.1d0,a,e,new, &
         band_deposition=d)
    call expect(old%ierr==0.and.new%ierr==0.and. &
         abs(new%secondary_heating_energy_ev_cm3-old%secondary_heating_energy_ev_cm3)<1d-18.and. &
         abs(new%x_hydrogen_ii-old%x_hydrogen_ii)<1d-15, &
         'monochromatic node moments reproduce existing chemistry',failures)
    ! Two energies within the HeII-and-above band, chosen across the FS2010
    ! low-energy transition. Same count and mean is not the same deposition.
    a=0;e=0;d=0
    a(3,7)=1d-6;e(3,7)=102.5d0
    status=snrt_secondary_fractions_c(5d0,.1d0,f)
    status=snrt_secondary_fractions_c(200d0,.1d0,g)
    d(3)=a(3,7);d(4:8)=.5d-6*(5*f+200*g)
    call snrt_thermochemistry_advance_cell(1d0,.079d0,1d0,1d4,0d0,.1d0,.1d0,.1d0,a,e,old)
    call snrt_thermochemistry_advance_cell(1d0,.079d0,1d0,1d4,0d0,.1d0,.1d0,.1d0,a,e,broad,band_deposition=d)
    delta=abs(broad%secondary_heating_energy_ev_cm3/old%secondary_heating_energy_ev_cm3-1)
    write(*,*)'FS2010_NODE_VS_MEAN relative heating difference=',delta
    call expect(old%ierr==0.and.broad%ierr==0.and.delta>1d-3.and. &
         abs(broad%photoelectron_energy_residual_ev_cm3)<1d-18, &
         'broad-spectrum integral differs from mean while conserving energy',failures)
    ! Primary reservation happens before any secondary demand, independent
    ! of source group/order. No helium targets: its requested energy heats.
    a=0;e=0;d=0;a(1,7)=.1d0;e(1,7)=200
    d(1)=.1d0;d(4:8)=[1d0,4d0,5d0,5d0,5d0]
    call snrt_thermochemistry_advance_cell(1d0,0d0,1d0,1d4,0d0,.89d0,0d0,0d0,a,e,new,band_deposition=d)
    expected=.01d0
    call expect(new%ierr==0.and.abs(new%secondary_hydrogen_ionizations_cm3-expected)<1d-15.and. &
         new%secondary_helium_i_ionizations_cm3==0.and.new%secondary_helium_ii_ionizations_cm3==0.and. &
         abs(new%secondary_heating_energy_ev_cm3-(15-13.6d0*expected))<1d-14.and. &
         abs(new%photoelectron_energy_residual_ev_cm3)<1d-14, &
         'summed secondary caps preserve energy when H/He targets exhaust',failures)
    initial=d;d(1)=.2d0
    call snrt_thermochemistry_advance_cell(1d0,0d0,1d0,1d4,0d0,.89d0,0d0,0d0,a,e,new,band_deposition=d)
    call expect(new%ierr/=0.and.new%x_hydrogen_ii==.89d0, &
         'inconsistent double primary ledger rejects without state publication',failures)
    d=initial;d(4)=2*d(4)
    call snrt_thermochemistry_advance_cell(1d0,0d0,1d0,1d4,0d0,.89d0,0d0,0d0,a,e,new,band_deposition=d)
    call expect(new%ierr/=0,'inconsistent node energy sum rejects',failures)
    d=initial;d(4)=-1
    call snrt_thermochemistry_advance_cell(1d0,0d0,1d0,1d4,0d0,.89d0,0d0,0d0,a,e,new,band_deposition=d)
    call expect(new%ierr/=0,'negative node-channel energy rejects',failures)
  end subroutine
#ifdef SNRT_CHIMES
  subroutine check_chimes_molecular_photo(failures)
    use snrt_dust_receiver,only:snrt_dust_receiver_stage
    integer,intent(inout)::failures
    character(len=1000)::path
    type(c_ptr)::handle,mol
    integer(c_int)::nr,ns,status
    real(dp)::identity(32),elem(11),measured(11),q,a(chimes_ns),b(chimes_ns),c(chimes_ns),alpha(128,9)
    real(dp)::n(2,9),e(2,9),nn(2,9),ee(2,9),l(10),shield(2),h2_lost,co_lost
    real(dp)::ctl(9),full_ledger(11),temp,initial_u,final_u,states(chimes_ns,3),temps(3),factors(3)
    real(dp)::nwork(2,9),ework(2,9),errs(2),gn(9),ge(9),du(1),td(1),accepted(1)
    real(dp),parameter::edge(10)=[.01d0,1d0,5.6d0,11.2d0,13.6d0,24.59d0,54.42d0,500d0,2000d0,10000d0]
    real(dp)::loss
    integer::j,k,parts
    call get_environment_variable('SNRT_CHIMES_MOLECULAR_TABLE',path)
    if(len_trim(path)==0)return
    mol=c_null_ptr
    status=chimes_molecular_load(trim(path)//c_null_char,mol,identity)
    call expect(status==0,'native molecular bank maps pinned CHIMES reactions',failures)
    if(status/=0)return
    call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
    handle=c_null_ptr;nr=0;ns=0
    status=chimes_band_load(trim(path)//c_null_char,handle,nr,ns,identity)
    if(status/=0)then
       call expect(.false.,'molecular test atomic bank loads',failures)
       call chimes_molecular_free(mol);return
    endif
    elem=0;elem(1)=1;status=chimes_neutral(elem,a);a(2)=0;a(3)=1;a(1)=1;shield=1
    do j=1,9
       n(:,j)=[.002d0,.003d0]
       e(:,j)=n(:,j)*[.75d0*edge(j)+.25d0*edge(j+1),.25d0*edge(j)+.75d0*edge(j+1)]
       alpha(:,j)=j*1d-20
    enddo
    status=chimes_band_photo_molecular_groups(handle,mol,2,1d0,1d11,3d8,alpha,shield,0d0, &
         a,n,e,b,nn,ee,l,gn,ge)
    call expect(status==0,'nine-band grain receiver solves two different incident ray spectra',failures)
    do j=1,9
       loss=1-exp(-.3d0*j)
       call expect(abs(gn(j)-sum(n(:,j))*loss)<2d-9.and. &
            abs(ge(j)-sum(e(:,j))*loss)<2d-7*sum(e(:,j)), &
            'individual grain band matches its own analytic optical depth and energy',failures)
    enddo
    ! A tiny charged molecular tail must be depleted without publishing a
    ! negative CVODE constraint-correction remnant. No abundance floor.
    elem=0;elem(1)=1;status=chimes_neutral(elem,a)
    a(2)=0;a(3)=1;a(1)=1;a(145)=1d-150;alpha=0;shield=1
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d14,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0.and.all(b>=0).and.b(145)<a(145).and.all(nn>=0), &
         'vanishing HCO+ tail remains nonnegative under repeated accepted-step checks',failures)
    call expect(status==0.and.abs(sum(e-ee)-l(7))<1d-8.and.abs(sum(l(1:6))-l(7))<1d-8, &
         'molecular tail recovery preserves the original photon and energy ledgers',failures)
    ! Newly injected trace metals can be far below the original absolute
    ! species tolerance even when their relative conservation matters.
    elem=0;elem(1)=1;elem(2)=.08d0;elem(3)=1d-4;elem(5)=1d-4
    elem(4)=1d-13;elem(6)=1d-13;elem(9)=1d-14;elem(10)=1d-16
    status=chimes_neutral(elem,a);alpha=0;shield=1
    status=chimes_band_photo_molecular_step(handle,mol,2,1d-3,1d13,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0,'strong photo step retains freshly injected trace metal inventories',failures)
    if(status==0)then
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.all(abs(measured-elem)<=1d-8*max(elem,1d-20)), &
            'trace and abundant nuclei obey the same unchanged relative budget',failures)
    endif
    elem=0;elem(1)=1;elem(3)=1d-4;elem(5)=1d-4
    status=chimes_neutral(elem,a)
    a(2)=.6d0;a(138)=.2d0;a(8)=0;a(24)=0;a(149)=1d-4
    n=0;n(:,4)=[.002d0,.003d0];e=0;e(:,4)=n(:,4)*[11.7d0,12.9d0]
    alpha=0;shield=1
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0,'joint native H2/CO molecular spectral evolution',failures)
    h2_lost=a(138)-b(138);co_lost=a(149)-b(149)
    call expect(h2_lost>0.and.co_lost>0.and.abs(b(2)-a(2)-2*h2_lost)<1d-10.and. &
         abs(sum(b(8:9))-co_lost)<1d-10.and.abs(b(24)-co_lost)<1d-10, &
         'photodissociation publishes correct H2/CO product multiplicities',failures)
    ! Released CI can itself photoionize above 11.26 eV in this same solve.
    ! Its heat is additional to (not an error in) H2's direct heat.
    call expect(l(8)>h2_lost+co_lost+b(9).and.l(1)>=h2_lost*6.4d-13/1.602176634d-12-1d-10, &
         'fluorescent captures coexist with dissociation and secondary CI photoionization',failures)
    status=chimes_budget(b,measured,q)
    call expect(status==0.and.maxval(abs(measured-elem))<1d-9.and.abs(q)<1d-9.and. &
         abs(sum(n-nn)-l(8))<1d-9.and.abs(sum(e-ee)-l(7))<1d-8.and. &
         abs(sum(l(1:6))-l(7))<1d-8,'molecular nuclei/charge/photon/energy budgets agree',failures)
    c=b;alpha=1d-17
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0.and.l(7)>0.and.l(9)>0.and.b(138)>c(138), &
         'grain competition reduces molecular dissociation without duplicate photons',failures)
    call expect(abs(sum(n-nn)-l(8)-l(10))<1d-9.and.abs(sum(e-ee)-l(7)-l(9))<1d-8, &
         'molecular plus grain combined budgets close',failures)
    alpha=0;shield=0
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0.and.all(b==a).and.all(nn==n).and.all(ee==e).and.all(l==0), &
         'explicit zero molecular coefficients preserve rays and chemistry exactly',failures)
    shield=-1;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status/=0.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and.all(l==-7), &
         'invalid molecular coefficient multiplier rejects before publication',failures)
    elem=0;elem(1)=1;status=chimes_neutral(elem,a);a(2)=.6d0;a(138)=.2d0;shield=1
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,0d0,a,n,e,b,nn,ee,l)
    call expect(status==0.and.abs(l(1)-(a(138)-b(138))*6.4d-13/1.602176634d-12)<1d-10, &
         'pure-H2 direct heat equals native per-dissociation energy, fluorescence is not thermalized',failures)
    status=chimes_band_photo_molecular_step(handle,mol,2,1d0,1d10,3d8,alpha,shield,1d0,a,n,e,c,nn,ee,l)
    call expect(status==0.and.all(abs(b-c)<1d-10).and.l(1)>(a(138)-b(138))*6.4d-13/1.602176634d-12.and. &
         l(1)<l(7).and.l(6)>=0,'density-quenched H2 pumping uses counted non-dissociating photons',failures)
    ! Weakly ionized molecular gas. The exactly electron-free dark CHIMES
    ! corner can create ~1e-27 H- from its numerical floor; the existing
    ! charge reconciliation correctly rejects that state rather than adding
    ! electrons or silently erasing negative ions to pass this test.
    a(1)=1d-4;a(3)=1d-4;a(2)=a(2)-1d-4
    ctl=[1d3,1d3,20d0,1d9,1d18,.1d0,1d0,0d0,.01d0];alpha=1d-19
    status=chimes_molecular_factors(ctl(2),ctl(1),ctl(5),a,factors)
    call expect(status==0.and.all(factors>=0).and.all(factors<=1).and.factors(1)<1, &
         'native H2/CO column factors and pumping use admitted cold tables',failures)
    temp=-7;b=-7;nn=-7;ee=-7;full_ledger=-7
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,elem,a,alpha,n,e,temp,b,nn,ee,full_ledger,gn,ge)
    if(status/=0)write(*,*)'COLD_SPLIT_REJECT_STATUS=',status
    call expect(status==0.and.temp>=10.and.temp<=1d5.and.b(138)>0, &
         'cold spectral receiver reaches actual molecular CHIMES formation/cooling without deleting H2',failures)
    if(status==0)then
       initial_u=1.5d0*chimes_boltzmann()*ctl(1)*sum(a)*ctl(2)/1.602176634d-12
       final_u=1.5d0*chimes_boltzmann()*ctl(1)*sum(b)*temp/1.602176634d-12
       call expect(abs(final_u-initial_u-full_ledger(1)-full_ledger(11))<1d-7.and. &
            abs(sum(e-ee)-full_ledger(7)-full_ledger(9))<1d-8, &
            'cold split separates gas heat, dark chemistry and grain energy without a second photon debit',failures)
       call expect(abs(sum(gn)-full_ledger(10))<1d-14.and.abs(sum(ge)-full_ledger(9))<1d-14.and. &
            all(gn(:3)==0).and.all(gn(5:)==0).and.gn(4)>0, &
            'grain group counters preserve the photon band and aggregate ledger',failures)
       ! Deliberately use nominal means = 1 eV: the material receiver MUST
       ! use actual accepted spectral energy instead, with no second opacity.
       call snrt_dust_receiver_stage(reshape(gn,[9,1]),spread(1d0,1,9),ctl(4),[ctl(6)],[1d-15], &
            [1d-12],[20d0],du,td,accepted,ierr,defer_temperature=.true., &
            deposited_spectrum_erg_cm3=reshape(ge*1.602176634d-12,[9,1]))
       call expect(ierr==0.and.abs(accepted(1)-full_ledger(9)*1.602176634d-12)<1d-26.and. &
            abs(du(1)-1d-12-accepted(1))<1d-26.and.td(1)==20, &
            'actual grain spectral energy reaches the existing dust material stage once',failures)
    endif
    do j=1,3
       parts=2**(j-1);ctl(4)=1d9/parts;ctl(2)=1d3;c=a;nwork=n;ework=e;temp=ctl(2)
       do k=1,parts
          status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,elem,c,alpha,nwork,ework, &
               temp,b,nn,ee,full_ledger)
          call expect(status==0,'cold molecular timestep refinement step',failures)
          if(status/=0)write(*,*)'COLD_SPLIT_REFINE_REJECT=',j,k,status,ctl(2),ctl(4)
          if(status/=0)exit
          c=b;nwork=nn;ework=ee;ctl(2)=temp
       enddo
       states(:,j)=c;temps(j)=temp
    enddo
    errs=[maxval(abs(states(:,1)-states(:,2))),maxval(abs(states(:,2)-states(:,3)))]
    write(*,*)'CHIMES_COLD_SPLIT_STATE_DIFF=',errs
    write(*,*)'CHIMES_COLD_SPLIT_T_DIFF=',abs(temps(1)-temps(2)),abs(temps(2)-temps(3))
    call expect(errs(2)<errs(1).and.abs(temps(2)-temps(3))<abs(temps(1)-temps(2)), &
         'cold molecular split converges under dt halving',failures)
    ctl(2)=1.001d5;temp=-7;b=-7;nn=-7;ee=-7;full_ledger=-7;gn=-7;ge=-7
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,elem,a,alpha,n,e,temp,b,nn,ee,full_ledger,gn,ge)
    call expect(status/=0.and.temp==-7.and.all(b==-7).and.all(full_ledger==-7).and. &
         all(nn==-7).and.all(ee==-7).and.all(gn==-7).and.all(ge==-7), &
         'unsupported high-temperature molecular state rejects intact instead of atomizing it',failures)
    ctl=[1d0,9d4,20d0,1d8,1d18,0d0,1d0,1d-6,.01d0];alpha=0;n=0;e=0
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,elem,a,alpha,n,e,temp,b,nn,ee,full_ledger,gn,ge)
    write(*,*)'COLD_INTERNAL_DOMAIN_STATUS=',status
    call expect(status==50.and.temp==-7.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and. &
         all(gn==-7).and.all(ge==-7).and.all(full_ledger==-7), &
         'internal dark thermal excursion outside molecular tables rejects every staged output',failures)
    ctl=[1d3,1d3,20d0,1d9,1d18,.1d0,1d0,0d0,.01d0]
    a(1)=0;a(3)=0;a(2)=.6d0
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,elem,a,alpha,n,e,temp,b,nn,ee,full_ledger,gn,ge)
    write(*,*)'COLD_EXACT_NEUTRAL_STATUS=',status
    if(status==0)then
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.abs(q)<1d-10.and.maxval(abs(measured-elem))<1d-8, &
            'exactly neutral dark state preserves nuclei and charge on success',failures)
    else
       call expect(status==46.and.temp==-7.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and. &
            all(gn==-7).and.all(ge==-7).and.all(full_ledger==-7), &
            'electron-free numerical-floor charge failure leaves all caller outputs intact',failures)
    endif
    if(chimes_transition_supported()==1)call check_chimes_transition(handle,mol,failures)
    call chimes_molecular_free(mol);call chimes_band_free(handle)
    elem=0;elem(1)=1;status=chimes_neutral(elem,a)
    a(4)=1d-315
    status=chimes_reconcile(elem,a,b)
    call expect(status==0.and.b(1)==0.and.b(4)==a(4).and.b(2)==a(2), &
         'underflow-only negative electron sum retains all ions and atoms without an electron floor',failures)
    a(4)=1d-300
    status=chimes_reconcile(elem,a,b)
    call expect(status/=0,'normal negative electron requirement still rejects, however small',failures)
  end subroutine

  subroutine check_chimes_transition(handle,mol,failures)
    type(c_ptr),intent(in)::handle,mol
    integer,intent(inout)::failures
    real(dp)::a(157),b(157),c(157),el(11),measured(11),q,qn,cost,t,tmax,elapsed,root_time
    real(dp)::ctl(9),n(2,9),e(2,9),nn(2,9),ee(2,9),alpha(128,9),ledger(11),gn(9),ge(9),events(2)
    real(dp)::before,after,temps(6),photo_state(157),reference_cost
    integer::s,status,j
    real(dp),parameter::ev=1.602176634d-12
    tmax=chimes_molecular_temperature_max()
    do s=138,157
       a=0;a(2)=1;a(s)=.01d0
       status=chimes_budget(a,el,q);a(1)=q
       before=1.5d0*chimes_boltzmann()*sum(a)*1d6
       status=chimes_atomize(1d6,a,b,t,cost)
       call expect(status==0.and.cost>0.and.all(b(138:157)==0).and.b(1)==a(1), &
            'source-backed atomization retains electrons for each of twenty species',failures)
       if(status/=0)cycle
       after=1.5d0*chimes_boltzmann()*sum(b)*t
       status=chimes_budget(b,measured,qn)
       call expect(status==0.and.maxval(abs(measured-el))<1d-12.and.abs(qn)<1d-12.and. &
            abs(after+cost*ev-before)<1d-12*before,'atomization conserves each nucleus charge and thermal plus binding',failures)
    enddo
    el=0;el(1)=1;status=chimes_neutral(el,a);a(2)=.8d0;a(138)=.1d0
    status=chimes_atomize(1.01d5,a,b,t,cost)
    reference_cost=.1d0*432.068d3/(6.02214076d23*1.602176634d-19)
    call expect(status==0.and.abs(cost-reference_cost)<1d-12.and.t<1.01d5*.9d0, &
         'H2 uses ATcT dissociation cost and removes the legacy artificial eleven-percent heat',failures)
    b=-7;t=-7;cost=-7
    status=chimes_atomize(10d0,a,b,t,cost)
    call expect(status/=0.and.all(b==-7).and.t==-7.and.cost==-7,'unaffordable atomization rejects transaction intact',failures)
    temps=[tmax*(1-1d-8),tmax*(1+1d-8),1d5,1d6,1d8,1d9]
    n=0;e=0;alpha=1d-22
    do j=1,size(temps)
       ctl=[1d0,temps(j),20d0,1d-12,1d18,.1d0,1d0,0d0,.01d0]
       before=1.5d0*chimes_boltzmann()*sum(a)*ctl(2)
       status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,el,a,alpha,n,e,t,b,nn,ee,ledger,gn,ge, &
            transition=.true.,event_info=events)
       call expect(status==0,'transition admits cold edge through one billion K without a dead temperature band',failures)
       if(status/=0)cycle
       after=1.5d0*chimes_boltzmann()*sum(b)*t
       call expect(abs(after-before-(ledger(1)+ledger(11))*ev)<1d-10*before.and.all(nn==0).and.all(ee==0), &
            'complete spectral split includes entry dissociation cost without inventing photons',failures)
    enddo
    ctl=[1d0,9d4,20d0,1d8,1d18,0d0,1d0,1d-6,.01d0]
    status=chimes_cell_transition_dark(ctl,el,a,0,t,b,elapsed)
    write(*,*)'TRANSITION_ROOT_STATUS_TIME_T=',status,elapsed,t
    call expect(status==51.and.elapsed>0.and.elapsed<ctl(4).and.abs(t-tmax)<1d-7*tmax, &
         'cold dark evolution returns an accepted CVODE thermal root before invalid molecular rates',failures)
    if(status==51)then
       root_time=elapsed;ctl(4)=root_time/2
       status=chimes_cell_transition_dark(ctl,el,a,0,t,b,elapsed)
       call expect(status==0.and.t<tmax.and.elapsed==ctl(4), &
            'root detection never uses an internal overshoot after the requested endpoint',failures)
       ctl(4)=1d8
       status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,el,a,alpha,n,e,t,b,nn,ee,ledger,gn,ge, &
            transition=.true.,event_info=events)
       call expect(status==0.and.events(1)>=1.and.events(2)>0, &
            'dark thermal root atomizes survivors and finishes the atomic remainder',failures)
    endif
    ! Photo heating across the edge, with competing grains and full rollback.
    ctl=[1d0,9d4,20d0,1d8,1d18,.1d0,1d0,0d0,.01d0]
    n=0;n(:,7)=1d5;e=0;e(:,7)=n(:,7)*400d0;alpha=1d-22
    before=1.5d0*chimes_boltzmann()*sum(a)*ctl(2)
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,el,a,alpha,n,e,t,b,nn,ee,ledger,gn,ge, &
         transition=.true.,event_info=events)
    write(*,*)'TRANSITION_PHOTO_STATUS_T_EVENTS=',status,t,events
    call expect(status==0.and.events(1)>=1.and.sum(ge)>0, &
         'photo heating crosses the boundary while gas and grains share one finite photon budget',failures)
    if(status==0)then
       after=1.5d0*chimes_boltzmann()*sum(b)*t
       call expect(abs(after-before-(ledger(1)+ledger(11))*ev)<1d-8*max(before,after).and. &
            abs(sum(e-ee)-ledger(7)-ledger(9))<1d-7*sum(e), &
            'photo-boundary dissociation and grain absorption close their separate energy ledgers',failures)
    endif
    ctl(2)=1.001d9;b=-7;t=-7;nn=-7;ee=-7;ledger=-7;events=-7
    status=chimes_cell_band_cold_molecular(handle,mol,2,ctl,el,a,alpha,n,e,t,b,nn,ee,ledger,gn,ge, &
         transition=.true.,event_info=events)
    call expect(status/=0.and.all(b==-7).and.t==-7.and.all(events==-7).and.all(nn==-7), &
         'finite upper physical domain still rejects instead of claiming unbounded temperature',failures)
  end subroutine

  subroutine check_chimes_competing_dust(failures)
    integer,intent(inout)::failures
    character(len=1000)::path
    type(c_ptr)::handle
    integer(c_int)::nr,ns,status
    real(dp)::identity(32),elem(11),a(chimes_ns),b(chimes_ns),c(chimes_ns),alpha(128,9)
    real(dp)::n(2,9),e(2,9),nn(2,9),ee(2,9),cn(2,9),ce(2,9),l(10),cl(8),loss
    call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
    if(len_trim(path)==0)return
    handle=c_null_ptr;nr=0;ns=0
    status=chimes_band_load(trim(path)//c_null_char,handle,nr,ns,identity)
    call expect(status==0,'competitive gas/dust bank loads',failures)
    if(status/=0)return
    elem=0;elem(1)=1;status=chimes_neutral(elem,a)
    a(2)=0;a(3)=1;a(1)=1
    n=0;n(:,5)=[.002d0,.003d0];e=20d0*n;alpha=1d-20
    status=chimes_band_photo_dust_step(handle,2,1d0,1d11,3d8,alpha,a,n,e,b,nn,ee,l)
    loss=1-exp(-.3d0)
    call expect(status==0,'dust-only finite photon solve',failures)
    call expect(maxval(abs(nn-n*(1-loss)))<1d-10.and.maxval(abs(ee-e*(1-loss)))<2d-9, &
         'dust-only survival matches analytic exponential for both rays',failures)
    write(*,*)'DUST_ONLY_RESIDUAL=',maxval(abs(b-a)),maxval(abs(l(1:8))), &
         abs(l(9)-sum(e)*loss),abs(l(10)-sum(n)*loss)
    call expect(all(b==a).and.all(l(1:8)==0).and.abs(l(9)-sum(e)*loss)<2d-9.and. &
         abs(l(10)-sum(n)*loss)<1d-10,'grain ledger separate from gas chemistry/heat',failures)
    status=chimes_neutral(elem,a)
    alpha=1d-18
    status=chimes_band_photo_dust_step(handle,2,1d0,1d10,3d8,alpha,a,n,e,b,nn,ee,l)
    call expect(status==0.and.l(7)>0.and.l(9)>0,'neutral gas and dust both capture from one photon population',failures)
    call expect(abs(sum(n-nn)-l(8)-l(10))<1d-9.and.abs(sum(e-ee)-l(7)-l(9))<2d-8.and. &
         abs(sum(l(1:6))-l(7))<2d-8,'competitive gas/grain photon and energy closure',failures)
    call expect(abs(b(3)-a(3)-l(8))<1d-9,'one neutral H ion per gas photon; dust photons do not ionize',failures)
    alpha=0
    status=chimes_band_photo_dust_step(handle,2,1d0,1d10,3d8,alpha,a,n,e,b,nn,ee,l)
    status=status+chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,c,cn,ce,cl)
    call expect(status==0.and.all(b==c).and.all(nn==cn).and.all(ee==ce).and.all(l(1:8)==cl).and. &
         all(l(9:10)==0),'zero-grain path exactly matches existing atomic API',failures)
    alpha=1d-12
    status=chimes_band_photo_dust_step(handle,2,1d0,1d10,3d8,alpha,a,n,e,b,nn,ee,l)
    call expect(status==0.and.sum(nn)<1d-9.and.l(9)>0,'optically thick dust exhaustion preserves nonnegative photons',failures)
    alpha=1d-16
    a(9)=1d-200;a(10)=1d-250;a(1)=1d-200+2d-250
    status=chimes_band_photo_dust_step(handle,2,2d3,1d12,3d8,alpha,a,n,e,b,nn,ee,l)
    call expect(status==0.and.all(nn>=0).and.sum(nn)<1d-9, &
         'dense gas and grains exhaust photons on an accepted nonnegative state',failures)
    alpha(1,1)=-1;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_band_photo_dust_step(handle,2,1d0,1d10,3d8,alpha,a,n,e,b,nn,ee,l)
    call expect(status/=0.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and.all(l==-7), &
         'invalid grain opacity rejects without publishing any state',failures)
    nn=0;nn(:,5)=[1d-200,1d-42];ee=20*nn;cn=nn;ce=ee
    status=chimes_round_subnormal_survivors(1d0,sum(ee),nn,ee)
    call expect(status/=0.and.all(nn==cn).and.all(ee==ce), &
         'storage rounding rejects a tail significant to the entire incoming photon budget',failures)
    status=chimes_round_subnormal_survivors(1d0,1d-3,nn,ee)
    call expect(status==0.and.nn(1,5)==0.and.ee(1,5)==0.and.nn(2,5)>0.and. &
         abs(ee(2,5)/nn(2,5)-20d0)<1d-12.and.sum(abs(ee-ce))<64*epsilon(1d0)*1d-3, &
         'negligible FP32 tail rounds paired photon N/E without a physical abundance floor',failures)
    call chimes_band_free(handle)
  end subroutine

  subroutine check_chimes_hot_split(failures)
    integer,intent(inout)::failures
    character(len=1000)::path
    type(c_ptr)::handle
    integer(c_int)::nr,ns,status
    real(dp)::identity(32),elem(11),a(chimes_ns),b(chimes_ns),c(chimes_ns),states(chimes_ns,3)
    real(dp)::ctl(9),n(1,9),e(1,9),nn(1,9),ee(1,9),n2(1,9),e2(1,9),l(9),l2(9)
    real(dp)::temp,t(3),ntotal,initial_u,final_u,zero(9),pn(9),measured(11),q,errors(2)
    integer::j,k,parts
    call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
    if(len_trim(path)==0)return
    handle=c_null_ptr;nr=0;ns=0;identity=0
    status=chimes_band_load(trim(path)//c_null_char,handle,nr,ns,identity)
    if(status/=0)then
       call expect(.false.,'hot split bank loads',failures)
       return
    endif
    elem=0;elem(1)=1;elem(2)=.08d0;elem(3)=1d-4
    status=chimes_neutral(elem,a)
    a(2)=.9d0;a(3)=.1d0;a(1)=.1d0
    ctl=[.01d0,1d6,20d0,1d8,1d17,0d0,1d0,0d0,.01d0]
    n=0;e=0;n(1,5)=.005d0;e(1,5)=13.6d0*n(1,5)
    temp=-7;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,a,0d0,n,e,temp,b,nn,ee,l)
    call expect(status==0,'spectral photo result reaches actual nonradiative CHIMES chemistry/cooling',failures)
    if(status==0)then
       initial_u=1.5d0*chimes_boltzmann()*ctl(1)*sum(a)*ctl(2)/1.602176634d-12
       final_u=1.5d0*chimes_boltzmann()*ctl(1)*sum(b)*temp/1.602176634d-12
       call expect(abs(final_u-initial_u-l(1)-l(9))<1d-10, &
            'split gas energy counts photo heat and dark thermal change once with updated particle count',failures)
       call expect(abs(sum(e)-sum(ee)-l(7))<1d-8.and.l(9)/=0, &
            'nonradiative CHIMES evolves thermal gas without a second photon debit',failures)
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.maxval(abs(measured-elem))<1d-8.and.abs(q)<1d-8, &
            'combined hot photo/nonradiative state conserves nuclei and charge',failures)
    endif
    ! Fixed monochromatic spectrum: changing dt measures the first-order
    ! physical operator split, not a repeated intragroup projection change.
    do j=1,3
       parts=2**(j-1);ctl(4)=1d8/parts;ctl(2)=1d6;c=a;n2=n;e2=e
       do k=1,parts
          status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,c,0d0,n2,e2,temp,b,nn,ee,l2)
          if(status/=0)exit
          c=b;n2=nn;e2=ee;ctl(2)=temp
       enddo
       call expect(status==0,'hot split dt refinement completes',failures)
       states(:,j)=c;t(j)=ctl(2)
    enddo
    errors(1)=maxval(abs(states(:,1)-states(:,2)))
    errors(2)=maxval(abs(states(:,2)-states(:,3)))
    write(*,*)'CHIMES spectral hot split state dt differences=',errors
    write(*,*)'CHIMES spectral hot split temperature dt differences=',abs(t(1)-t(2)),abs(t(2)-t(3))
    call expect(errors(2)<.8d0*errors(1).and.errors(1)<1d-4, &
         'first-order spectral/nonradiative splitting error decreases under dt refinement',failures)
    ! Exactly dark input must reduce to the already validated CHIMES cell.
    ctl(2)=1d6;ctl(4)=1d8;n=0;e=0;zero=0
    status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,a,0d0,n,e,temp,b,nn,ee,l)
    call expect(status==0,'dark hot split comparison succeeds',failures)
    if(status==0)then
       status=chimes_cell(ctl,elem,a,zero,ntotal,c,pn)
       call expect(status==0.and.maxval(abs(c-b))<1d-12.and.abs(ntotal/temp-1)<1d-12, &
            'zero-radiation hot split equals the existing CHIMES receiver',failures)
    endif
    ctl(2)=1.00001d5;ctl(4)=1d9;n(1,5)=.005d0;e(1,5)=13.6d0*n(1,5)
    temp=-7;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,a,0d0,n,e,temp,b,nn,ee,l)
    call expect(status==8.and.temp==-7.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and.all(l==-7), &
         'post-photo crossing into molecular domain rolls back the whole staged split',failures)
    ctl(2)=1d4;temp=-7;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,a,0d0,n,e,temp,b,nn,ee,l)
    call expect(status/=0.and.temp==-7.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and.all(l==-7), &
         'cold molecular-domain split rejects without partial publication',failures)
    ctl(2)=1d6;a(2)=a(2)-.2d0;a(138)=.1d0
    status=chimes_cell_band_hot_atomic(handle,1,ctl,elem,a,0d0,n,e,temp,b,nn,ee,l)
    call expect(status/=0.and.all(b==-7),'hot split never silently discards preexisting H2',failures)
    call chimes_band_free(handle)
  end subroutine
  subroutine check_chimes_photo(failures)
    integer,intent(inout)::failures
    character(len=1000)::path
    type(c_ptr)::handle
    integer(c_int)::nr,ns,status,ps(2)
    real(dp)::identity(32),elem(11),a(chimes_ns),b(chimes_ns),c(chimes_ns),d(chimes_ns)
    real(dp)::n(2,9),e(2,9),nn(2,9),ee(2,9),n2(2,9),e2(2,9),l(8),l2(8),l3(8)
    real(dp)::pb(chimes_ns,2),pn(2,9,2),pe(2,9,2),pl(8,2),measured(11),q,reference_fs(5)
    integer::k
    call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
    if(len_trim(path)==0)return
    handle=c_null_ptr;nr=0;ns=0;identity=0
    status=chimes_band_load(trim(path)//c_null_char,handle,nr,ns,identity)
    call expect(status==0,'photo operator loads native atomic bank',failures)
    if(status/=0)return
    elem=0;elem(1)=1
    status=chimes_neutral(elem,a)
    n=0;e=0;n(1,5)=.4d0;e(1,5)=.4d0*20d0
    b=-7;nn=-7;ee=-7;l=-7
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'CVODE atomic photon/species integration succeeds',failures)
    if(status==0)then
       call expect(abs((b(3)-a(3))-(sum(n)-sum(nn))-l(2)/13.6d0)<1d-7.and.b(2)>0, &
            'HI primary and secondary ionization use the same finite photon ledger',failures)
       call expect(abs(sum(e)-sum(ee)-l(7))<1d-6.and.abs(sum(l(1:6))-l(7))<1d-6, &
            'actual photon energy equals heat plus secondary/excitation and binding reservoir',failures)
       call expect(l(1)>0.and.abs(l(6)-13.6d0*l(8))<1d-7.and.all(l(3:4)==0), &
            'HI shell binding follows primary captures; absent helium receives no secondary events',failures)
    endif
    ! Monochromatic endpoint: two operator half steps do not add a spectral
    ! reconstruction error, isolating the time integration comparison.
    n=0;e=0;n(1,5)=2d0;e(1,5)=2d0*13.6d0
    status=chimes_band_photo_step(handle,2,1d0,1d9,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'photon-rich finite HI target integration succeeds',failures)
    if(status==0)then
       call expect(sum(n)-sum(nn)<=1d0+1d-7.and.b(2)>=0.and.b(3)<=1d0+1d-8, &
            'finite HI inventory cannot consume more primary photons than atoms',failures)
       status=chimes_band_photo_step(handle,2,1d0,5d8,3d8,a,n,e,c,n2,e2,l2)
       if(status==0)status=chimes_band_photo_step(handle,2,1d0,5d8,3d8,c,n2,e2,d,pn(:,:,1),pe(:,:,1),l3)
       call expect(status==0,'two atomic half steps succeed',failures)
       if(status==0)call expect(maxval(abs(d-b))<2d-7.and.maxval(abs(pn(:,:,1)-nn))<2d-7, &
            'atomic whole/half-step integration converges for unchanged monochromatic closure',failures)
    endif
    n=0;e=0;n(1,5)=.1d0;e(1,5)=13.6d0*n(1,5)
    status=chimes_band_photo_step(handle,2,1d0,1d12,3d8,a,n,e,b,nn,ee,l)
    write(*,*)'CHIMES optically thick photo status=',status
    call expect(status==0,'optically thick photon-starved atomic solve succeeds',failures)
    if(status==0)call expect(sum(nn)<1d-8.and.abs(b(3)-.1d0)<1d-7.and.b(2)>.89d0, &
         'photon exhaustion leaves the correct neutral inventory without over-ionization',failures)
    ! H- photodetachment: opposite ions initially make a neutral gas.
    a=0;a(3)=.5d0;a(4)=.5d0;n=0;e=0;n(1,1)=.1d0;e(1,1)=.1d0
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'H- detachment uses its own threshold and electron multiplicity',failures)
    if(status==0)then
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.abs(q)<1d-8.and.abs(l(6)-.755d0*l(8))<1d-8.and.b(1)>0, &
            'H- detachment conserves charge and retains the electron-affinity energy',failures)
    endif
    ! Photo-only H2 ionization is supported; a later hot CHIMES call may
    ! not erase its molecular products (checked by the split wrapper).
    a=0;a(138)=.5d0;n=0;e=0;n(1,5)=.1d0;e(1,5)=.1d0*24.59d0
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'H2 photoionization retains an H2+ carrier rather than atomizing it',failures)
    if(status==0)then
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.abs(q)<1d-8.and.abs(measured(1)-1)<1d-8.and.b(139)>0, &
            'H2 photoionization conserves nuclei and molecular charge',failures)
    endif
    status=chimes_neutral(elem,a);a(2)=.9d0;a(3)=.1d0;a(1)=.1d0
    n=0;e=0;n(1,9)=.001d0;e(1,9)=10d0
    status=chimes_band_photo_step(handle,2,1d0,1d6,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'small hard-photon perturbation integrates',failures)
    if(status==0)then
       status=chimes_secondary_partition(10000d0-13.6d0,a,reference_fs)
       call expect(status==0.and.maxval(abs(l(1:5)/sum(l(1:5))-reference_fs))<1d-7, &
            'compact raw-grid FS interpolation matches the established native physical partition',failures)
    endif
    elem(2)=.08d0;elem(3)=1d-4;elem(11)=1d-5
    status=chimes_neutral(elem,a)
    n=0;e=0;n(1,9)=.001d0;e(1,9)=10d0
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0,'hard-photon metal/Auger and shell-wise secondary solve succeeds',failures)
    if(status==0)then
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.maxval(abs(measured-elem))<1d-8.and.abs(q)<1d-8, &
            'photo/Auger/secondary reactions conserve CHIMES nuclei and charge',failures)
       call expect(sum(l(2:5))>0.and.l(6)>0.and.abs(sum(l(1:6))-l(7))<1d-6, &
            'hard photons partition secondary energy and retain unresolved binding energy',failures)
       call expect(sum(b(112:137))>0,'hard spectrum populates Fe ionization/Auger products',failures)
!$omp parallel do private(k)
       do k=1,2
          ps(k)=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,pb(:,k),pn(:,:,k),pe(:,:,k),pl(:,k))
       enddo
!$omp end parallel do
       call expect(all(ps==0),'concurrent native CVODE photo solves succeed',failures)
       do k=1,2
          call expect(all(pb(:,k)==b).and.all(pn(:,:,k)==nn).and.all(pe(:,:,k)==ee).and.all(pl(:,k)==l), &
               'private CVODE/photo state is thread-independent bitwise',failures)
       enddo
    endif
    status=chimes_band_photo_step(handle,2,1d0,0d0,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0.and.all(b==a).and.all(nn==n).and.all(ee==e).and.all(l==0), &
         'zero photo timestep is exact identity',failures)
    n=0;e=0
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status==0.and.all(b==a).and.all(l==0),'zero radiation photo operator is identity',failures)
    n(1,5)=-1;b=-7;nn=-7;ee=-7;l=-7
    status=chimes_band_photo_step(handle,2,1d0,1d10,3d8,a,n,e,b,nn,ee,l)
    call expect(status/=0.and.all(b==-7).and.all(nn==-7).and.all(ee==-7).and.all(l==-7), &
         'invalid photo input cannot partially publish species, photons or heat',failures)
    call chimes_band_free(handle)
  end subroutine
  subroutine check_chimes_band(failures)
    integer,intent(inout)::failures
    character(len=1000)::path
    type(c_ptr)::handle
    integer(c_int)::nr,ns,status,mapping(5,311),parallel_status(4)
    integer::hi,j,k
    real(dp)::identity(32),n(2,9),e(2,9),one_n(1,9),one_e(1,9)
    real(dp)::mom(311,3,9),a(311,3,9),b(311,3,9),parallel(311,3,9,4),expected,scale,x
    call get_environment_variable('SNRT_CHIMES_BAND_TABLE',path)
    if(len_trim(path)==0)return
    handle=c_null_ptr;nr=0;ns=0;identity=-1
    status=chimes_band_load(trim(path)//c_null_char,handle,nr,ns,identity)
    call expect(status==0.and.c_associated(handle).and.nr==311.and.ns>nr, &
         'atomic spectral bank loads complete reaction and partial-shell mapping',failures)
    if(status/=0)return
    status=chimes_band_reactions(handle,nr,mapping)
    hi=0
    do j=1,nr
       if(mapping(1,j)==1.and.mapping(3,j)==1)hi=j
    enddo
    call expect(status==0.and.hi>0.and.all(identity>=0).and.all(identity<=255), &
         'atomic bank exposes CHIMES HI mapping and content digest',failures)
    if(hi==0)then
       call chimes_band_free(handle)
       return
    endif
    n=0;e=0;n(:,5)=1;e(:,5)=(/13.6d0,24.59d0/)
    status=chimes_band_moments(handle,2,nr,n,e,mom)
    ! Independent analytic Verner96 HI evaluation at the two delta endpoints.
    expected=0
    do j=1,2
       x=e(j,5)/.4298d0
       expected=expected+1d-18*5.475d4*((x-1)**2)*x**(.5d0*2.963d0-5.5d0)* &
            (1+sqrt(x/32.88d0))**(-2.963d0)
    enddo
    call expect(status==0.and.abs(mom(hi,1,5)/expected-1)<1d-12, &
         'spectral HI opacity matches independent Verner endpoint evaluation',failures)
    call expect(abs(mom(hi,3,5)-(mom(hi,2,5)-13.6d0*mom(hi,1,5)))<1d-28, &
         'HI primary electron energy subtracts the actual shell binding',failures)
    call expect(all(mom>=0).and.all(mom(:,3,:)<=mom(:,2,:)), &
         'all atomic shell moments are positive with primary energy bounded by absorbed energy',failures)
    one_n(1,:)=n(1,:);one_e(1,:)=e(1,:)
    status=chimes_band_moments(handle,1,nr,one_n,one_e,a)
    call expect(status==0,'first separate spectral ray evaluates',failures)
    one_n(1,:)=n(2,:);one_e(1,:)=e(2,:)
    status=chimes_band_moments(handle,1,nr,one_n,one_e,b)
    scale=maxval(mom)
    call expect(status==0.and.maxval(abs(mom-a-b))<1d-13*scale, &
         'mixed rays equal separately reconstructed spectral moments',failures)
    one_n(1,:)=sum(n,dim=1);one_e(1,:)=sum(e,dim=1)
    status=chimes_band_moments(handle,1,nr,one_n,one_e,a)
    call expect(status==0.and.abs(a(hi,1,5)/mom(hi,1,5)-1)>1d-3, &
         'directional spectral hardening is not replaced by a merged mean',failures)
!$omp parallel do private(k)
    do k=1,4
       parallel_status(k)=chimes_band_moments(handle,2,nr,n,e,parallel(:,:,:,k))
    enddo
!$omp end parallel do
    call expect(all(parallel_status==0),'shared immutable spectral bank works in threaded calls',failures)
    do k=1,4
       call expect(all(parallel(:,:,:,k)==mom),'threaded atomic moments equal serial bitwise',failures)
    enddo
    n=0;e=0;n(1,4)=1;e(1,4)=13.6d0
    status=chimes_band_moments(handle,2,nr,n,e,a)
    call expect(status==0.and.all(a(hi,:,:)==0),'upper band left limit excludes threshold double assignment',failures)
    n=0;e=0;n(1,9)=1;e(1,9)=10000
    status=chimes_band_moments(handle,2,nr,n,e,a)
    call expect(status==0.and.sum(a(125:311,1,9))>0.and.all(a(:,3,:)<=a(:,2,:)), &
         'hard photons retain nonzero inner-shell Auger channels',failures)
    n(1,9)=0;a=-7
    status=chimes_band_moments(handle,2,nr,n,e,a)
    call expect(status/=0.and.all(a==-7),'invalid N/E rejects without partial moment publication',failures)
    call chimes_band_free(handle)
  end subroutine
  subroutine check_chimes(failures)
    use dust_pah_radiation, only: pah_coronene_photoionize
    integer,intent(inout)::failures
    character(len=500)::path,group_dir,group_paths(9)
    character(len=2)::group_number
    real(dp)::ctl(9),elem(11),a(chimes_ns),b(chimes_ns),c(chimes_ns),phot(9),pnew(9),temp,q,measured(11)
    real(dp)::target(11),locked(11),first_photons(9),identity(320)
    real(dp)::parallel_state(chimes_ns,4),parallel_photons(9,4),parallel_temperature(4)
    real(dp)::deposit(5),reference(5),lo(5),hi(5),absorptions,thermal_change
    real(dp)::solid_q(4),charged_initial(chimes_ns,4),serial_state(chimes_ns),serial_photons(9),serial_temperature
    real(dp)::pah_pop(2,2),pah_captures(1,2),pah_vib(1,2),ne,pah_heat,pah_ip,initial_energy
    integer::parallel_status(4)
    integer::status,ng,j
    call get_environment_variable('SNRT_CHIMES_MAIN_DATA',path)
    call get_environment_variable('SNRT_CHIMES_GROUP_DIR',group_dir)
    ng=0;group_paths=c_null_char
    if(len_trim(group_dir)>0)then
       ng=9
       do j=1,9
          write(group_number,'(I2.2)')j
          group_paths(j)=trim(group_dir)//'/group_'//group_number//'.hdf5'//c_null_char
       enddo
    endif
    status=chimes_initialize(trim(path)//c_null_char,ng,group_paths,500)
    call expect(status==0,'native CHIMES source-backed reaction/cooling tables load',failures)
    if(status/=0)return
    status=chimes_identity(identity)
    call expect(status==0.and.all(identity>=0).and.all(identity<=255).and.any(identity>0), &
         'native chemistry computes data SHA256 identities',failures)
    elem=0;elem(1)=1;phot=0
    status=chimes_neutral(elem,a)
    a(2)=.999d0;a(3)=.001d0;a(1)=.001d0
    ctl=[100d0,100d0,20d0,1d7,1d18,0d0,1d0,0d0,1d0]
    solid_q=[1d-4,-5d-4,1d-8,0d0]
    do j=1,4
       status=chimes_reconcile_charged(elem,a,solid_q(j),charged_initial(:,j))
       call expect(status==0.and.abs(charged_initial(1,j)-a(1)-solid_q(j))<1d-18, &
            'reconcile includes signed solid charge without altering nuclei',failures)
    enddo
    status=chimes_reconcile_charged(elem,a,-.002d0,b)
    call expect(status/=0.and.all(b==a),'impossible negative solid charge cannot manufacture electrons',failures)
    if(chimes_charge_supported()==1)then
       !$omp parallel do private(j)
       do j=1,4
          parallel_status(j)=chimes_cell_charged(ctl,elem,charged_initial(:,j),solid_q(j),phot, &
               parallel_temperature(j),parallel_state(:,j),parallel_photons(:,j))
       enddo
       !$omp end parallel do
       call expect(all(parallel_status==0),'solid-charge cells evolve through the native network',failures)
       do j=1,4
          status=chimes_budget(parallel_state(:,j),measured,q)
          call expect(status==0.and.abs(q+solid_q(j))<1d-15.and.abs(measured(1)-1)<1d-12, &
               'gas plus solid charge is conserved, not forcibly gas-neutralized',failures)
          status=chimes_cell_charged(ctl,elem,charged_initial(:,j),solid_q(j),phot, &
               serial_temperature,serial_state,serial_photons)
          call expect(status==0.and.all(serial_state==parallel_state(:,j)).and. &
               serial_temperature==parallel_temperature(j).and.all(serial_photons==parallel_photons(:,j)), &
               'different solid charges have identical serial/OpenMP solutions',failures)
       enddo
       status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
       call expect(status==0.and.all(b==parallel_state(:,4)).and.temp==parallel_temperature(4), &
            'zero solid charge exactly preserves the original chemistry path',failures)
       status=chimes_cell_charged(ctl,elem,a,solid_q(1),phot,temp,b,pnew)
       call expect(status/=0.and.all(b==a).and.all(pnew==phot).and.temp==ctl(2), &
            'mismatched gas/solid charge rejects without changing the trial',failures)
       ! Photoionization produces real additional electrons before the
       ! chemistry substep; reconstruct T from energy AND particle count.
       ctl(4)=1d0;initial_energy=1.5d0*chimes_boltzmann()*ctl(1)*sum(a)*ctl(2)
       pah_pop=0;pah_pop(1,1)=1d-3;pah_captures=0;pah_captures(1,1)=1d-6
       ne=a(1)*ctl(1);pah_heat=0;pah_ip=0;pah_vib=0
       call pah_coronene_photoionize([10d0],pah_captures,pah_pop,ne,pah_vib,pah_heat,pah_ip,status)
       call expect(status==0.and.ne>a(1)*ctl(1).and.pah_heat>0, &
            'PAH absorption produces photoelectrons and their matched gas heat',failures)
       a(1)=ne/ctl(1)
       ctl(2)=(initial_energy+pah_heat)/(1.5d0*chimes_boltzmann()*ctl(1)*sum(a))
       status=chimes_cell_charged(ctl,elem,a,sum(pah_pop(:,2))/ctl(1),phot,temp,b,pnew)
       call expect(status==0,'native chemistry accepts the actual PAH photoelectron state',failures)
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.abs(q+sum(pah_pop(:,2))/ctl(1))<1d-16, &
            'CHIMES retains PAH-origin charge through chemical/thermal evolution',failures)
    else
       status=chimes_cell_charged(ctl,elem,charged_initial(:,1),solid_q(1),phot,temp,b,pnew)
       call expect(status==4.and.all(b==charged_initial(:,1)).and.all(pnew==phot), &
            'ABI4 library explicitly rejects charged solids but retains neutral compatibility',failures)
    endif
    elem=[1d0,.08d0,2.46d-4,8.51d-5,4.90d-4,1d-4,3.47d-5,3.47d-5,1.86d-5,2.29d-6,2.82d-5]
    status=chimes_neutral(elem,a)
    ctl=[100d0,100d0,20d0,1d10,1d18,1d0,1d0,0d0,1d0];phot=0
    a(2)=1-1d-4;a(3)=1d-4;a(1)=1d-4
    ctl(4)=1d0
    status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
    write(*,*)'CHIMES one-second zero-field T=',temp
    call expect(status==0.and.abs(temp/ctl(2)-1)<1d-6, &
         'zero-field RT network retains the short-step thermal limit',failures)
    ctl(4)=1d10
    status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
    write(*,*)'CHIMES cold status/T/H2/CO=',status,temp,b(138),b(149)
    call expect(status==0.and.b(138)>0.and.b(149)>0.and.temp>0, &
         'native dust H2 formation and CO chemistry with thermal evolution',failures)
    status=chimes_budget(b,measured,q)
    call expect(status==0.and.maxval(abs(measured/elem-1))<1d-5.and.abs(q)<1d-8, &
         '157-species native network conserves all eleven elements and charge',failures)
    ctl(6)=0
    status=chimes_cell(ctl,elem,a,phot,temp,c,pnew)
    call expect(status==0.and.c(138)<b(138),'surface formation responds to actual dust abundance',failures)
    ctl(2)=1d6;ctl(4)=1d8
    status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
    call expect(status==0.and.b(9)>0.and.b(112)>0.and.b(8)<a(8), &
         'hot gas evolves carbon and iron ion populations, not a CIE table lookup',failures)
    if(ng==9)then
       ctl(2)=1d4;ctl(4)=1d8;ctl(9)=.01d0;phot=1d-3
       status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
       write(*,*)'CHIMES RT status/T/photon removal=',status,temp,sum(phot-pnew)
       call expect(status==0.and.sum(pnew)<sum(phot).and.b(9)>0, &
            'native nine-group RT depletes ionizing photons and evolves carbon',failures)
       call expect(temp<10100d0,'photoheating stays within the supplied photon energy budget',failures)
       status=chimes_neutral(elem,a)
       a(138)=.1d0;a(2)=.8d0;a(149)=1d-4;a(8)=elem(3)-a(149);a(24)=elem(5)-a(149)
       ctl=[100d0,100d0,20d0,1d8,0d0,0d0,1d0,0d0,.01d0]
       phot=0;phot(4)=.01d0
       status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
       first_photons=pnew
       write(*,*)'CHIMES molecular thin status/H2/CO/photons=',status,b(138),b(149),pnew(4)
       call expect(status==0.and.b(138)<a(138).and.b(149)<a(149).and.pnew(4)<phot(4), &
            'H2 and CO dissociation consumes the actual LW photon group',failures)
       ctl(5)=1d20
       status=chimes_cell(ctl,elem,a,phot,temp,c,pnew)
       write(*,*)'CHIMES molecular shielded status/H2/CO/photons=',status,c(138),c(149),pnew(4)
       call expect(status==0.and.c(138)>b(138).and.c(149)>b(149).and.pnew(4)>first_photons(4), &
            'local molecular columns suppress dissociation and its matched photon sink',failures)
       call expect(phot(4)-first_photons(4)>=ctl(1)*max(0d0,a(138)-b(138)+a(149)-b(149)), &
            'molecular destruction does not exceed the removed photon count',failures)
       target=elem;target(3)=target(3)*1.1d0
       status=chimes_reconcile(target,a,b)
       call expect(status==0.and.b(149)==a(149).and.b(8)>a(8), &
            'neutral carbon injection preserves existing CO molecules',failures)
       status=chimes_locked(a,locked)
       call expect(status==0.and.locked(3)==a(149).and.locked(5)==a(149), &
            'CO reserves carbon and oxygen against duplicate grain growth',failures)
       target(3)=locked(3)*.5d0
       status=chimes_reconcile(target,a,b)
       call expect(status/=0.and.all(b==a),'grain growth cannot silently consume locked CO',failures)
       !$omp parallel do private(j)
       do j=1,4
          parallel_status(j)=chimes_cell(ctl,elem,a,phot,parallel_temperature(j),parallel_state(:,j),parallel_photons(:,j))
       enddo
       !$omp end parallel do
       call expect(all(parallel_status==0).and.all(parallel_state(:,1)==parallel_state(:,4)).and. &
            all(parallel_temperature==parallel_temperature(1)), &
            'OpenMP cells use independent native chemistry state',failures)
    endif
    ctl(4)=-1
    status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
    call expect(status/=0.and.all(b==a).and.all(pnew==phot),'native chemistry rejects invalid step without publication',failures)
    elem=0;elem(1)=1;elem(2)=.25d0/(4*.74d0);elem(11)=.01d0/(56*.74d0)
    status=chimes_neutral(elem,a)
    ctl=[.00074d0,1d4,20d0,1d13,1d20,0d0,1d0,0d0,1d0];phot=0
    status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
    call expect(status==0.and.all(b>=0),'zero-abundance elements in the initial live composition are supported',failures)
    status=chimes_nuclear_sums(-a,measured)
    call expect(status==0.and.maxval(abs(measured+elem))<1d-14, &
         'signed face flux carries the same eleven nuclear budgets',failures)
    a=0;a(2)=.9d0;a(3)=.1d0;a(5)=.072d0;a(6)=.008d0;a(1)=.108d0
    status=chimes_secondary_partition(200d0,a,deposit)
    call snrt_secondary_fractions(200d0,.1d0,reference(1),reference(2),reference(3), &
         reference(4),reference(5),j)
    call expect(status==0.and.maxval(abs(deposit-reference))<1d-14, &
         'CHIMES electron partition matches native FS2010 in its atomic reference state',failures)
    status=chimes_secondary_partition(99.9d0,a,lo)
    status=chimes_secondary_partition(100.1d0,a,hi)
    call expect(status==0.and.maxval(abs(hi-lo))<5d-3, &
         'CHIMES electron partition has no artificial 100 eV switch',failures)
    a=0;a(138)=.5d0
    status=chimes_secondary_partition(200d0,a,deposit)
    call expect(status==0.and.deposit(1)==1d0.and.all(deposit(2:5)==0d0), &
         'secondary atomic channels cannot consume molecular or absent helium targets',failures)
    a=0;a(2)=1;a(5)=1d-30
    status=chimes_secondary_partition(200d0,a,deposit)
    call expect(status==0.and.deposit(3)<1d-28.and.deposit(4)==0d0, &
         'secondary helium rates vanish continuously with total helium abundance',failures)
    if(ng==9)then
       elem=0;elem(1)=1;status=chimes_neutral(elem,a)
       a(2)=1d0-1d-4;a(3)=1d-4;a(1)=1d-4
       phot=0;phot(9)=1d3
       ctl=[1d0,100d0,20d0,1d5,0d0,0d0,1d0,0d0,1d0]
       status=chimes_cell(ctl,elem,a,phot,temp,b,pnew)
       absorptions=sum(phot-pnew)
       thermal_change=1.5d0*chimes_boltzmann()*(sum(b)*temp-sum(a)*ctl(2))
       write(*,*)'CHIMES hard photons / additional HII / thermal erg: ',absorptions,b(3)-a(3),thermal_change
       call expect(status==0.and.absorptions>0.and.b(3)-a(3)>2*absorptions, &
            'live native reaction RHS produces multiple ionizations per absorbed hard photon',failures)
       call expect(thermal_change>0.and.thermal_change+(b(3)-a(3))*13.6d0*1.602176634d-12 &
            <absorptions*4023.594574013186d0*1.602176634d-12, &
            'secondary ionization is charged to absorbed energy, not added on top of full photoheating',failures)
       status=chimes_budget(b,measured,q)
       call expect(status==0.and.abs(measured(1)-1d0)<1d-12.and.abs(q)<1d-12, &
            'hard-electron chemistry preserves nuclei and charge',failures)
       !$omp parallel do private(j)
       do j=1,4
          parallel_status(j)=chimes_cell(ctl,elem,a,phot,parallel_temperature(j),parallel_state(:,j),parallel_photons(:,j))
       enddo
       !$omp end parallel do
       call expect(all(parallel_status==0).and.all(parallel_state(:,1)==parallel_state(:,4)).and. &
            all(parallel_temperature==parallel_temperature(1)), &
            'OpenMP hard-photon cells share read-only FS tables without sharing chemical state',failures)
    endif
  end subroutine
#endif

  subroutine check_atomic(failures)
    integer,intent(inout)::failures
    real(dp)::g(11),rho,x(3),y(3),z(3),w(3),e,e1,e2,e3,loss,rate,rate2
    real(dp)::beta(3),alpha(3),rec(3),exc(2),bre,die,expected,t
    integer::status,j
    g=0;g(1:2)=[.76d0,.24d0];rho=atomic_mh/.76d0
    do j=1,2
       t=10d0**(2*j+2)
       call atomic_rates(t,beta,alpha,rec,exc,bre,die,status)
       expected=5.85d-11*sqrt(t)*exp(-157809.1d0/t)/(1+sqrt(t/1d5))
       call expect(status==0.and.abs(beta(1)/expected-1)<1d-14.and. &
            alpha(3)==2*snrt_alpha_hydrogen_case_b(t/4),'atomic rate reference at 1e4/1e6 K',failures)
    enddo
    x=0;e=atomic_heat_capacity(rho,g,x,5d0/3)*1d4
    call atomic_advance(rho,g,5d0/3,1d0,1d11,e,x,e1,y,loss,status)
    call expect(status==0.and.e1==e.and.all(y==x),'neutral electron-free gas: no fabricated collisions',failures)
    x=[.8d0,.3d0,.2d0];e=atomic_heat_capacity(rho,g,x,5d0/3)*1d4
    call atomic_advance(rho,g,5d0/3,1d0,1d11,e,x,e1,y,loss,status)
    call expect(status==0.and.y(1)<x(1).and.e1>0.and.e1<e.and. &
         abs(e1+loss-e)<epsilon(e)*e,'case-B recombination and positive thermal budget',failures)
    x=[.01d0,.01d0,.001d0];e=atomic_heat_capacity(rho,g,x,5d0/3)*1d6
    call atomic_advance(rho,g,5d0/3,1d0,1d10,e,x,e1,y,loss,status)
    call expect(status==0.and.y(1)>x(1).and.all(y>=0).and.sum(y(2:3))<=1.and.e1>0, &
         'hot collisional ionization without photons',failures)
    call atomic_advance(rho,g,5d0/3,1d0,5d9,e,x,e2,z,loss,status)
    call atomic_advance(rho,g,5d0/3,1d0,5d9,e2,z,e3,w,loss,status)
    write(*,*)'ATOMIC_DT relative E, absolute fractions: ',abs(e3/e1-1),maxval(abs(w-y))
    call expect(status==0.and.abs(e3/e1-1)<3d-3.and.maxval(abs(w-y))<3d-3, &
         'atomic full versus two half timesteps',failures)
    call atomic_advance(rho,g,5d0/3,1d0,-1d0,e,x,e1,y,loss,status)
    call expect(status/=0.and.e1==e.and.all(y==x).and.loss==0,'atomic rejection leaves input state intact',failures)
    g(1)=.759d0;g(3)=.001d0
    call wss09_metals_rate(1d5,g,rate,status)
    g(3)=.0005d0
    call wss09_metals_rate(1d5,g,rate2,status)
    call expect(status==0.and.rate>0.and.abs(rate2/rate-.5d0)<1d-14, &
         'metal-only CIE follows the depleted element budget',failures)
    call wss09_metals_rate(99d0,g,rate,status)
    call expect(status/=0,'nonzero metals reject out-of-table temperature',failures)
  end subroutine check_atomic

  subroutine expect(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (condition) then
       write(*,'(a)') 'PASS: ' // trim(label)
    else
       failures = failures + 1
       write(*,'(a)') 'FAIL: ' // trim(label)
    end if
  end subroutine expect

end program snrt_thermochemistry_smoke

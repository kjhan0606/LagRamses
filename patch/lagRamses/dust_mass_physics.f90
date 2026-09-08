! Bulk/two-size dust mass closures and bounded native size-shift reference. Metal fields
! denote total (gas+dust) metals; dust is a subset, not additional gas mass.
module dust_mass_physics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: dust_mass_dp=real64
  logical :: dust_mass_enabled=.false.,dust_growth=.true.,dust_sputtering=.true.
  character(len=32) :: dust_mass_model='bulk_v1',dust_cooling='none'
  character(len=32) :: dust_material_model='fixed_mix'
  character(len=32) :: dust_optics_model='fixed_mix'
  ! Species indices follow H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe in the stellar contract.
  integer,parameter :: dust_nc=2,dust_ne=11
  real(real64),parameter :: olivine_mass=24d0+56d0+28d0+4*16d0
  real(real64),parameter :: olivine_fraction(dust_ne)= &
       [0d0,0d0,0d0,0d0,64d0,0d0,24d0,28d0,0d0,0d0,56d0]/olivine_mass
  real(real64) :: dust_condensation(3)=[0d0,.2d0,.15d0] ! wind, AGB, SNII; no Ia
  real(real64) :: dust_grain_radius_cm=1d-5,dust_grain_density=3d0
  real(real64) :: dust_sticking=.3d0,dust_growth_max_temperature=300d0
  real(real64) :: dust_metal_atom_mass=24d0,dust_injection_temperature=20d0
  real(real64),parameter :: dust_mp=1.67262192369d-24,dust_kb=1.380649d-16
  ! Four transported masses: C-small, C-large, silicate-small, silicate-large.
  ! Dubois+2024 two-size reference; resolved dense gas only for coagulation.
  real(real64) :: dust_size_radius_cm(2)=[5d-7,1d-5],dust_size_density(2)=[2.2d0,3.3d0]
  real(real64) :: dust_small_injection_fraction(2)=0d0
  logical :: dust_coagulation=.true.,dust_shattering=.true.
  logical :: dust_sn_shocks=.false. ! explicit energy-equivalent unresolved comparison
  real(real64),parameter :: dust_shock_delta(2)=[.10d0,.15d0],dust_shock_mass_msun=6800d0
  real(real64),parameter :: dust_myr=315576d8
contains
  subroutine dust_bin_reconstruct(left,right,number,mass_moment,new_number,slope,ierr)
    ! McKinnon et al. 2018, eqs 39--43: linear dn/da around bin midpoint.
    ! mass_moment is integral a^3 dn, i.e. mass/(4*pi*rho_s/3).
    ! Positivity limiting conserves MASS, not necessarily NUMBER. The caller
    ! must record that number change rather than claim both were conserved.
    real(real64),intent(in)::left,right,number,mass_moment
    real(real64),intent(out)::new_number,slope
    integer,intent(out)::ierr
    real(real64)::c,h,a,b,limit
    ierr=1;new_number=number;slope=0
    if(.not.all(ieee_is_finite([left,right,number,mass_moment])))return
    if(left<0.or.right<=left.or.min(number,mass_moment)<0)return
    if(number==0.or.mass_moment==0)then
       if(number/=0.or.mass_moment/=0)return
       ierr=0;return
    endif
    c=(left+right)/2;h=(right-left)/2
    a=c**3+c*h*h;b=2*h**3*(c*c+h*h/5)
    if(b<=0.or..not.ieee_is_finite(b))return
    slope=(mass_moment-number*a)/b;limit=number/(2*h*h)
    if(abs(slope)>limit)then
       new_number=mass_moment/(a+sign(1d0,slope)*b/(2*h*h))
       slope=sign(new_number/(2*h*h),slope)
    endif
    if(.not.all(ieee_is_finite([new_number,slope])).or.new_number<0)return
    ierr=0
  end subroutine

  subroutine dust_bin_moments(left,right,number,slope,lo,hi,shift,moments)
    ! Exact three-node Gaussian integration of a linear number distribution
    ! times (a+shift)^k, k=0,2,3 (polynomial degree <=4). Local coordinates
    ! avoid subtracting nearly equal high powers at narrow bin boundaries.
    real(real64),intent(in)::left,right,number,slope,lo,hi,shift
    real(real64),intent(out)::moments(3)
    real(real64),parameter::node(3)=[-.774596669241483377d0,0d0,.774596669241483377d0]
    real(real64),parameter::weight(3)=[5d0/9,8d0/9,5d0/9]
    real(real64)::x,y,f,h,c
    integer::q
    moments=0
    if(hi<=lo)return
    h=(hi-lo)/2;c=(hi+lo)/2
    do q=1,3
       x=c+h*node(q);y=x+shift
       f=number/(right-left)+slope*(x-(left+right)/2)
       moments=moments+h*weight(q)*f*[1d0,y*y,y*y*y]
    enddo
  end subroutine

  subroutine dust_multibin_shift(edges,number,slope,shift,solid_density,gas_available, &
       new_number,new_slope,gas_transfer,limiter_number_change,ierr,destroyed_number)
    ! Bounded native size-advection reference, NOT a live many-bin selector.
    ! Constant da/dt over this transaction: mantle accretion or sputtering.
    ! Grains below edges(1) are destroyed; their residual mass returns to gas.
    ! Upper overflow is rejected (never silently rebinned/clipped). No new
    ! collision kernel is invented to imitate the two-size transfer closure.
    real(real64),intent(in)::edges(:),number(:),slope(:),shift,solid_density,gas_available
    real(real64),intent(out)::new_number(:),new_slope(:),gas_transfer,limiter_number_change
    integer,intent(out)::ierr
    real(real64),optional,intent(out)::destroyed_number
    real(real64)::nn(size(number)),ss(size(number)),mm(size(number)),m(3),old_mass,factor
    real(real64)::lo,hi,h,tol,exchange,correction,raw_number,removed
    integer::i,j,n,status
    ierr=1;gas_transfer=0;limiter_number_change=0
    if(present(destroyed_number))destroyed_number=0
    n=size(number)
    if(size(new_number)/=n.or.size(new_slope)/=n)return
    new_number=number;new_slope=0
    if(n<1.or.size(edges)/=n+1.or.size(slope)/=n)return
    new_slope=slope
    if(.not.all(ieee_is_finite(edges)).or..not.all(ieee_is_finite(number)))return
    if(.not.all(ieee_is_finite(slope)).or..not.all(ieee_is_finite([shift,solid_density,gas_available])))return
    if(edges(1)<0.or.any(edges(2:)<=edges(:n)).or.any(number<0))return
    if(solid_density<=0.or.gas_available<0)return
    factor=4*acos(-1d0)*solid_density/3
    nn=0;mm=0;old_mass=0;removed=0
    do i=1,n
       h=(edges(i+1)-edges(i))/2
       if(abs(slope(i))*h>number(i)/(2*h)*(1+32*epsilon(1d0)))return
       call dust_bin_moments(edges(i),edges(i+1),number(i),slope(i),edges(i),edges(i+1),0d0,m)
       old_mass=old_mass+factor*m(3)
       hi=min(edges(i+1),edges(1)-shift)
       if(hi>edges(i))then
          call dust_bin_moments(edges(i),edges(i+1),number(i),slope(i),edges(i),hi,0d0,m)
          removed=removed+m(1)
       endif
       if(shift>0.and.number(i)>0)then
          ! A positive linear distribution has support in the bin interior.
          if(edges(i+1)+shift>edges(n+1))return
       endif
       do j=1,n
          lo=max(edges(i),edges(j)-shift);hi=min(edges(i+1),edges(j+1)-shift)
          if(hi<=lo)cycle
          call dust_bin_moments(edges(i),edges(i+1),number(i),slope(i),lo,hi,shift,m)
          nn(j)=nn(j)+m(1);mm(j)=mm(j)+m(3)
       enddo
    enddo
    if(.not.ieee_is_finite(old_mass).or.any(mm<0).or.any(nn<0))return
    exchange=old_mass-factor*sum(mm)
    tol=128*epsilon(1d0)*max(old_mass,factor*sum(mm),gas_available,tiny(1d0))
    if(.not.ieee_is_finite(exchange).or.exchange < -gas_available-tol)return
    correction=0
    do j=1,n
       raw_number=nn(j)
       call dust_bin_reconstruct(edges(j),edges(j+1),raw_number,mm(j),nn(j),ss(j),status)
       if(status/=0)return
       correction=correction+nn(j)-raw_number
    enddo
    if(.not.ieee_is_finite(correction))return
    ! Publish only after all bins and the shared gas-reservoir bound pass.
    new_number=nn;new_slope=ss;gas_transfer=exchange;limiter_number_change=correction;ierr=0
    if(present(destroyed_number))destroyed_number=removed
  end subroutine

  logical function dust_mass_parameters_ok() result(ok)
    real(real64)::v(9)
    v=[dust_condensation,dust_grain_radius_cm,dust_grain_density,dust_sticking, &
       dust_growth_max_temperature,dust_metal_atom_mass,dust_injection_temperature]
    ok=all(ieee_is_finite(v)).and.all(v>=0)
    ok=ok.and.all(dust_condensation<=1).and.dust_sticking<=1
    ok=ok.and.dust_grain_radius_cm>0.and.dust_grain_density>0.and.dust_metal_atom_mass>0
    ok=ok.and.dust_growth_max_temperature>0.and.dust_injection_temperature>0
    ok=ok.and.(trim(dust_mass_model)=='bulk_v1'.or.trim(dust_mass_model)=='carbon_olivine_v1'.or. &
         trim(dust_mass_model)=='carbon_olivine_2size_v1')
    ok=ok.and.(trim(dust_cooling)=='none'.or.trim(dust_cooling)=='depleted_scalar'.or. &
         trim(dust_cooling)=='wss09_cie'.or.trim(dust_cooling)=='snrt_hhe_cie_metals')
    if(trim(dust_cooling)=='wss09_cie'.or.trim(dust_cooling)=='snrt_hhe_cie_metals') &
         ok=ok.and.dust_composition_enabled()
    ok=ok.and.all(ieee_is_finite(dust_size_radius_cm)).and.all(dust_size_radius_cm>0)
    ok=ok.and.dust_size_radius_cm(1)<dust_size_radius_cm(2)
    ok=ok.and.all(ieee_is_finite(dust_size_density)).and.all(dust_size_density>0)
    ok=ok.and.all(ieee_is_finite(dust_small_injection_fraction))
    ok=ok.and.all(dust_small_injection_fraction>=0).and.all(dust_small_injection_fraction<=1)
    if(dust_sn_shocks)ok=ok.and.dust_two_size_enabled()
    ok=ok.and.(trim(dust_material_model)=='fixed_mix'.or.trim(dust_material_model)=='dl01_composition_v1')
    if(trim(dust_material_model)=='dl01_composition_v1')ok=ok.and.dust_two_size_enabled()
    ok=ok.and.(trim(dust_optics_model)=='fixed_mix'.or.trim(dust_optics_model)=='d03_transport_v1')
    if(trim(dust_optics_model)=='d03_transport_v1')then
       ok=ok.and.dust_material_composition_enabled()
       ok=ok.and.all(dust_size_radius_cm==[1d-6,1d-5]).and.all(dust_size_density==[2.2d0,3.8d0])
    endif
    if(.not.dust_mass_enabled)ok=ok.and.trim(dust_mass_model)=='bulk_v1'.and.trim(dust_cooling)=='none'
  end function

  logical function dust_atomic_cooling_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_cooling)=='snrt_hhe_cie_metals'
  end function

  logical function dust_composition_enabled() result(enabled)
    enabled=dust_mass_enabled.and.(trim(dust_mass_model)=='carbon_olivine_v1'.or. &
         trim(dust_mass_model)=='carbon_olivine_2size_v1')
  end function

  logical function dust_two_size_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_mass_model)=='carbon_olivine_2size_v1'
  end function

  logical function dust_material_composition_enabled() result(enabled)
    enabled=dust_two_size_enabled().and.trim(dust_material_model)=='dl01_composition_v1'
  end function

  logical function dust_optics_enabled() result(enabled)
    enabled=dust_material_composition_enabled().and.trim(dust_optics_model)=='d03_transport_v1'
  end function

  function dust_size_identity() result(v)
    real(real64)::v(13)
    v=[1d0,dust_size_radius_cm,dust_size_density,dust_small_injection_fraction, &
         merge(1d0,0d0,dust_coagulation),merge(1d0,0d0,dust_shattering), &
         merge(1d0,0d0,dust_sn_shocks),dust_shock_delta,dust_shock_mass_msun]
  end function

  subroutine dust_shock_step(gas_mass_msun,sn_energy_erg,bins,fresh,next,ierr)
    ! Canonical E_SN=1e51 erg, N_eff=COUPLED SN energy / E_SN; not a count
    ! of physical explosions. Apply after MPI source summation, once per cell.
    ! Only ambient grains are exposed. All same-step fresh ejecta are protected;
    ! reverse-shock survival is already implicit in condensation efficiencies.
    real(real64),intent(in)::gas_mass_msun,sn_energy_erg,bins(4),fresh(4)
    real(real64),intent(out)::next(4)
    integer,intent(out)::ierr
    real(real64)::fraction,tau,loss,ambient,tol(4)
    integer::i,j,k
    ierr=1;next=bins
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([gas_mass_msun,sn_energy_erg,bins,fresh])))return
    if(gas_mass_msun<=0.or.min(sn_energy_erg,minval(bins),minval(fresh))<0)return
    tol=128*epsilon(1d0)*max(bins,fresh,tiny(1d0))
    if(any(fresh>bins+tol))return
    do i=1,2
       do j=1,2
          k=2*i+j-2
          fraction=(1-exp(-dust_shock_delta(i)*1d-5/dust_size_radius_cm(j)))* &
               min(dust_shock_mass_msun/gas_mass_msun,1d0)
          if(fraction<1d-5)then
             tau=fraction*(1+fraction/2+fraction*fraction/3)*(sn_energy_erg/1d51)
          else if(fraction<1)then
             tau=-log(1-fraction)*(sn_energy_erg/1d51)
          else
             tau=huge(1d0);if(sn_energy_erg==0)tau=0
          endif
          if(.not.ieee_is_finite(tau))return
          if(tau<1d-5)then
             loss=tau*(1-tau/2+tau*tau/6)
          else
             loss=1-exp(-tau)
          endif
          ambient=max(0d0,bins(k)-fresh(k))
          next(k)=bins(k)-ambient*loss
       enddo
    enddo
    ierr=0
  end subroutine

  function dust_injection_bins(grains) result(bins)
    real(real64),intent(in)::grains(2)
    real(real64)::bins(4)
    bins(1)=grains(1)*dust_small_injection_fraction(1);bins(2)=grains(1)-bins(1)
    bins(3)=grains(2)*dust_small_injection_fraction(2);bins(4)=grains(2)-bins(3)
  end function

  subroutine dust_size_exchange(rho,nh,temperature,material,dt,bins,ierr)
    ! Exact donor-quadratic transfer, dD_donor/dt=-k D_donor^2.
    ! Collisions within each composition; no creation/destruction of metals.
    real(real64),intent(in)::rho,nh,temperature,dt
    integer,intent(in)::material
    real(real64),intent(inout)::bins(2)
    integer,intent(out)::ierr
    real(real64)::rate,donor,transfer,x,p
    integer::j
    ierr=1
    if(.not.all(ieee_is_finite([rho,nh,temperature,dt,bins])))return
    if(rho<=0.or.min(nh,temperature,dt,minval(bins))<0.or.material<1.or.material>2)return
    rate=0;j=1
    if(dust_shattering.and.nh>0.and.nh<1d3)then
       p=1d0;if(nh>=1)p=1d0/3
       j=2
       rate=(bins(2)/rho/.01d0)*nh**p/(54*dust_myr)* &
            (1d-5/dust_size_radius_cm(2))*(3/dust_size_density(material))
    else if(dust_coagulation.and.nh>=1d3.and.temperature<1d4)then
       ! F=1, sigma_small=0.1 km/s, LOCAL resolved n_H. No hidden cloud boost.
       j=1
       rate=(bins(1)/rho/.01d0)*(nh/1d3)/(.27d0*dust_myr)* &
            (5d-7/dust_size_radius_cm(1))*(3/dust_size_density(material))
    endif
    x=rate*dt
    if(.not.ieee_is_finite(x).or.x<0)return
    if(x==0)then
       ierr=0;return
    endif
    donor=bins(j)/(1+x)
    transfer=bins(j)-donor
    bins(j)=donor;bins(3-j)=bins(3-j)+transfer;ierr=0
  end subroutine

  subroutine dust_size_step(rho,metal,elements,bins,temperature,scale_d,dt,next,ierr)
    ! Symmetric, positivity-preserving splits of shared-reservoir growth,
    ! size-dependent thermal erosion and exact coagulation/shattering flows.
    ! Growth retains the explicitly configured effective accreting atom mass;
    ! sputtering is the shared Tsai-Mathews reference, NOT species Hu fits.
    real(real64),intent(in)::rho,metal,elements(11),bins(4),temperature,scale_d,dt
    real(real64),intent(out)::next(4)
    integer,intent(out)::ierr
    real(real64)::gas(11),grains(2),cap(2),aa(2,2),bb(2,2),a,b,h,available,nh,steps,updated
    integer::i,j,k,m,first,status,ns,seq(3)
    next=bins;ierr=1
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([rho,metal,bins,temperature,scale_d,dt])))return
    if(rho<=0.or.scale_d<=0.or.min(metal,temperature,dt,minval(bins))<0.or.metal>rho)return
    grains=[sum(bins(1:2)),sum(bins(3:4))]
    if(sum(grains)>metal*(1+64*epsilon(1d0)))return
    call dust_gas_elements(elements,grains,gas,status)
    if(status/=0)return
    cap(1)=elements(3);cap(2)=huge(1d0)
    do j=1,11
       if(olivine_fraction(j)>0)cap(2)=min(cap(2),elements(j)/olivine_fraction(j))
    enddo
    if(sum(cap)>metal*(1+64*epsilon(1d0)))return
    nh=elements(1)*scale_d/dust_mp
    do i=1,2
       call dust_mass_rates(rho*scale_d,cap(i)*scale_d,temperature,a,b,status)
       if(status/=0)return
       aa(:,i)=a*dust_grain_radius_cm/dust_size_radius_cm*dust_grain_density/dust_size_density(i)
       bb(:,i)=b*dust_grain_radius_cm/dust_size_radius_cm
    enddo
    ! Resolve only the noncommuting growth/erosion splits. Exact collision
    ! transfer needs no rate-based subcycling when growth/erosion is absent.
    steps=maxval(aa+bb)*dt/.1d0
    if(.not.ieee_is_finite(steps).or.steps>4096)return
    ns=max(1,ceiling(steps));h=dt/ns;seq=[1,2,1]
    do k=1,ns
       do i=1,2
          first=2*i-1
          call dust_size_exchange(rho,nh,temperature,i,h/2,next(first:first+1),status)
          if(status/=0)return
          do m=1,3
             j=seq(m)
             available=max(0d0,cap(i)-next(first+2-j))
             a=0;if(cap(i)>0)a=aa(j,i)*available/cap(i)
             b=h/2;if(m==2)b=h
             call dust_mass_step(available,min(next(first+j-1),available),b,a,bb(j,i),updated,status)
             if(status/=0)return
             next(first+j-1)=updated
          enddo
          call dust_size_exchange(rho,nh,temperature,i,h/2,next(first:first+1),status)
          if(status/=0)return
       enddo
    enddo
    grains=[sum(next(1:2)),sum(next(3:4))]
    call dust_gas_elements(elements,grains,gas,ierr)
  end subroutine

  subroutine dust_species_condense(elements,channel,grains,ierr)
    ! Apply to one source-node RELEASE SEGMENT, before age/Z/IMF mixing.
    ! CO binds the smaller number inventory in AGB ejecta. It is not dust.
    real(real64),intent(in)::elements(dust_ne)
    integer,intent(in)::channel
    real(real64),intent(out)::grains(dust_nc)
    integer,intent(out)::ierr
    real(real64)::available(dust_ne),co,silicate
    integer::i
    grains=0;ierr=1
    if(.not.dust_mass_parameters_ok().or.channel<1.or.channel>5)return
    if(.not.all(ieee_is_finite(elements)).or.any(elements<0))return
    ierr=0
    if(channel>3)return ! No Ia/PISN condensation prescription.
    available=elements
    if(channel==2)then
       co=min(available(3)/12d0,available(5)/16d0)
       available(3)=max(0d0,available(3)-12*co)
       available(5)=max(0d0,available(5)-16*co)
    endif
    silicate=huge(1d0)
    do i=1,dust_ne
       if(olivine_fraction(i)>0)silicate=min(silicate,available(i)/olivine_fraction(i))
    enddo
    grains=dust_condensation(channel)*[available(3),silicate]
    if(.not.all(ieee_is_finite(grains)))ierr=1
  end subroutine

  subroutine dust_gas_elements(elements,grains,gas,ierr)
    real(real64),intent(in)::elements(dust_ne),grains(dust_nc)
    real(real64),intent(out)::gas(dust_ne)
    integer,intent(out)::ierr
    real(real64)::locked(dust_ne),tol(dust_ne)
    gas=0;ierr=1
    if(.not.all(ieee_is_finite(elements)).or..not.all(ieee_is_finite(grains)))return
    if(any(elements<0).or.any(grains<0))return
    locked=olivine_fraction*grains(2);locked(3)=grains(1)
    tol=64*epsilon(1d0)*max(elements,locked,tiny(1d0))
    if(any(locked>elements+tol))return
    gas=max(0d0,elements-locked);ierr=0
  end subroutine

  subroutine dust_species_step(rho,metal,elements,grains,temperature,scale_d,dt,next,ierr)
    real(real64),intent(in)::rho,metal,elements(dust_ne),grains(dust_nc),temperature,scale_d,dt
    real(real64),intent(out)::next(dust_nc)
    integer,intent(out)::ierr
    real(real64)::gas(dust_ne),capacity(dust_nc),a,b
    integer::i,j,status
    next=grains;ierr=1
    if(.not.all(ieee_is_finite([rho,metal,temperature,scale_d,dt])).or.scale_d<=0)return
    if(rho<=0.or.metal<0.or.metal>rho.or.sum(grains)>metal*(1+64*epsilon(1d0)))return
    call dust_gas_elements(elements,grains,gas,status)
    if(status/=0)return
    ! Disjoint C and silicate reservoirs. For silicates the gas-phase limiting
    ! element bounds growth; total capacity includes already locked atoms.
    capacity(1)=elements(3);capacity(2)=huge(1d0)
    do j=1,dust_ne
       if(olivine_fraction(j)>0)capacity(2)=min(capacity(2),elements(j)/olivine_fraction(j))
    enddo
    if(sum(capacity)>metal*(1+64*epsilon(1d0)))return
    do i=1,dust_nc
       call dust_mass_rates(rho*scale_d,capacity(i)*scale_d,temperature,a,b,status)
       if(status/=0)return
       call dust_mass_step(capacity(i),min(grains(i),capacity(i)),dt,a,b,next(i),status)
       if(status/=0)return
    enddo
    call dust_gas_elements(elements,next,gas,ierr)
  end subroutine

  function dust_mass_identity() result(v)
    real(real64)::v(13)
    v=[1d0,merge(1d0,0d0,dust_mass_enabled),merge(1d0,0d0,dust_growth), &
       merge(1d0,0d0,dust_sputtering),dust_condensation,dust_grain_radius_cm, &
       dust_grain_density,dust_sticking,dust_growth_max_temperature, &
       dust_metal_atom_mass,dust_injection_temperature]
  end function

  subroutine dust_condense(returned,hydrogen,helium,dust,ierr)
    real(real64),intent(in)::returned(3),hydrogen(3),helium(3)
    real(real64),intent(out)::dust
    integer,intent(out)::ierr
    real(real64)::metals(3)
    ierr=1;dust=0
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite(returned)).or..not.all(ieee_is_finite(hydrogen)).or. &
       .not.all(ieee_is_finite(helium)))return
    if(any(returned<0).or.any(hydrogen<0).or.any(helium<0))return
    metals=returned-hydrogen-helium
    if(any(metals < -1d-10*max(returned,tiny(1d0))))return
    if(dust_mass_enabled)dust=sum(dust_condensation*max(metals,0d0))
    if(.not.ieee_is_finite(dust))return
    ierr=0
  end subroutine

  subroutine dust_mass_rates(rho,metal_density,temperature,growth_rate,loss_rate,ierr)
    real(real64),intent(in)::rho,metal_density,temperature ! physical cgs
    real(real64),intent(out)::growth_rate,loss_rate
    integer,intent(out)::ierr
    real(real64)::speed,x
    growth_rate=0;loss_rate=0;ierr=1
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([rho,metal_density,temperature])))return
    if(min(rho,metal_density,temperature)<0.or.metal_density>rho)return
    if(dust_growth.and.temperature>0.and.temperature<=dust_growth_max_temperature)then
       speed=sqrt(8*dust_kb*temperature/(acos(-1d0)*dust_metal_atom_mass*dust_mp))
       ! Geometric collisions: pi*a^2*rho_Z*v*S / (4*pi*a^3*rho_s/3).
       ! The separate (1-D/Z) factor debits gas-phase metal availability.
       growth_rate=.75d0*dust_sticking*metal_density*speed/(dust_grain_density*dust_grain_radius_cm)
    endif
    if(dust_sputtering.and.temperature>0)then
       ! Tsai & Mathews erosion, McKinnon+2018 eq50; bulk mass lifetime a/(3|adot|).
       if(temperature<2d6)then
          x=(temperature/2d6)**2.5d0;x=x/(1+x)
       else
          x=1/(1+(2d6/temperature)**2.5d0)
       endif
       loss_rate=3*3.2d-18*(rho/dust_mp)*x/dust_grain_radius_cm
    endif
    if(.not.all(ieee_is_finite([growth_rate,loss_rate])))return
    ierr=0
  end subroutine

  subroutine dust_mass_step(metal,dust,dt,a,b,new_dust,ierr)
    ! Exact frozen-rate dD/dt = (A-B)D - A D^2/Z. Preserves D=0 and 0<=D<=Z.
    real(real64),intent(in)::metal,dust,dt,a,b
    real(real64),intent(out)::new_dust
    integer,intent(out)::ierr
    real(real64)::r,x,phi,decay,denominator
    ierr=1;new_dust=dust
    if(.not.all(ieee_is_finite([metal,dust,dt,a,b])))return
    if(min(metal,dust,dt,a,b)<0.or.dust>metal)return
    if(dust==0.or.dt==0.or.(a==0.and.b==0))then
       ierr=0;return
    endif
    r=a-b;x=r*dt
    if(.not.ieee_is_finite(x).or..not.ieee_is_finite(a*dt))return
    if(abs(x)<1d-5)then
       phi=dt*(1+x/2+x*x/6+x*x*x/24)
       new_dust=dust*(1+r*phi)/(1+a*(dust/metal)*phi)
    else if(r>0)then
       decay=exp(-x)
       denominator=decay+(a/r)*(dust/metal)*(1-decay)
       new_dust=dust/denominator
    else
       decay=exp(x);phi=(decay-1)/r
       new_dust=dust*decay/(1+a*(dust/metal)*phi)
    endif
    if(.not.ieee_is_finite(new_dust))return
    if(new_dust<0.or.new_dust>metal*(1+32*epsilon(1d0)))return
    new_dust=min(new_dust,metal);ierr=0
  end subroutine

  subroutine dust_mass_exchange(old_mass,new_mass,old_dust_energy,gas_thermal,new_energy,transfer,ierr)
    ! Fixed dust temperature during mass exchange; no latent heat model.
    ! transfer is charged to gas thermal energy, not CR/kinetic energy.
    real(real64),intent(in)::old_mass,new_mass,old_dust_energy,gas_thermal
    real(real64),intent(out)::new_energy,transfer
    integer,intent(out)::ierr
    new_energy=old_dust_energy;transfer=0;ierr=1
    if(.not.all(ieee_is_finite([old_mass,new_mass,old_dust_energy,gas_thermal])))return
    if(min(old_mass,new_mass,old_dust_energy,gas_thermal)<0)return
    if(old_mass==0)then
       if(new_mass/=0.or.old_dust_energy/=0)return
    else
       new_energy=old_dust_energy*(new_mass/old_mass)
       transfer=new_energy-old_dust_energy
    endif
    if(.not.ieee_is_finite(new_energy).or.transfer>gas_thermal)return
    ierr=0
  end subroutine
end module

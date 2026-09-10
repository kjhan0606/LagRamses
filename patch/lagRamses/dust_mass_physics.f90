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
  character(len=32) :: dust_sublimation='none'
  character(len=32) :: dust_iron_model='none'
  character(len=32) :: dust_pah_model='none'
  ! Neutral C24H12 comparison. Mass-number convention matches CHIMES.
  integer,parameter :: dust_pah_nbin=128
  real(real64),parameter :: dust_pah_hc(2)=[12d0,288d0]/300d0
  real(real64),parameter :: dust_pah_molecule_g=300d0*1.67262192369d-24
  real(real64) :: dust_pah_condensation=0d0
  ! Explicit electric-only, co-advected comparison; NOT full Fe dust physics.
  real(real64) :: dust_fe_condensation=0d0
  ! Opt-in geometric, unfocused accretion onto existing metallic Fe seeds.
  ! Zero disables accretion; kinetics=false preserves the fixed-Fe model.
  ! This is not a calibrated ion/charge model.
  real(real64) :: dust_fe_sticking=0d0
  logical :: dust_fe_kinetics=.false.
  ! Choban et al. 2026, Table 3, doi:10.1093/mnras/stag1020:
  ! log10(Y/[micron yr^-1 cm^3]) = sum_i c_i log10(T/K)^i.
  ! Fit to Nozawa+2006 Fe, not the C/silicate Tsai-Mathews law.
  real(real64),parameter :: dust_fe_sputter_coeff(0:5)= &
       [-156.88d0,82.110d0,-18.238d0,2.0692d0,-.11933d0,.0027788d0]
  real(real64),parameter :: dust_fe_sputter_min_t=1d4,dust_fe_sputter_max_t=1d9
  real(real64),parameter :: dust_fe_max_temperature=300d0,dust_fe_max_primary_ev=4d0
  ! GD89 graphite BULK vacuum evaporation; not Mg2SiO4 coefficients for our MgFeSiO4.
  real(real64),parameter :: dust_carbon_atom=12.011d0*1.66053906660d-24
  real(real64),parameter :: dust_carbon_nu=2d14,dust_carbon_binding_k=81200d0
  real(real64),parameter :: dust_carbon_latent=1.380649d-16*dust_carbon_binding_k/dust_carbon_atom
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
  subroutine dust_collision_deposit(mass,grain,number,boundary)
    ! Fixed pivots: one daughter, conserving number AND mass inside the grid.
    ! Outside mass is an inert numerical reservoir, NOT gas/PAH destruction.
    real(real64),intent(in)::mass(:),grain
    real(real64),intent(inout)::number(:),boundary(2)
    real(real64)::w
    integer::i,n
    n=size(mass)
    if(grain<=0)return
    if(grain<mass(1))then
       boundary(1)=boundary(1)+grain;return
    else if(grain>mass(n))then
       boundary(2)=boundary(2)+grain;return
    endif
    do i=1,n-1
       if(grain>mass(i+1))cycle
       w=(grain-mass(i))/(mass(i+1)-mass(i))
       number(i)=number(i)+1-w;number(i+1)=number(i+1)+w
       return
    enddo
  end subroutine

  subroutine dust_collision_fragments(mass,parent,eject,number,boundary)
    ! HA19 eqs22--25: alpha_f=3.3, mmax=.02*mej, mmin=1e-6*mmax.
    ! Exact integrated fragment MASS in geometric-pivot cells; number is
    ! mass/pivot, a convergent approximation, not exact fragment number.
    real(real64),intent(in)::mass(:),parent,eject
    real(real64),intent(inout)::number(:),boundary(2)
    real(real64),parameter::p=(4d0-3.3d0)/3
    real(real64)::lo,hi,fmin,fmax,norm,left,right
    integer::i,n
    n=size(mass)
    call dust_collision_deposit(mass,parent-eject,number,boundary)
    if(eject<=0)return
    fmax=.02d0*eject;fmin=1d-6*fmax
    norm=fmax**p-fmin**p
    hi=min(fmax,mass(1))
    if(hi>fmin)boundary(1)=boundary(1)+eject*(hi**p-fmin**p)/norm
    lo=max(fmin,mass(n))
    if(fmax>lo)boundary(2)=boundary(2)+eject*(fmax**p-lo**p)/norm
    left=mass(1)
    do i=1,n
       right=mass(n)
       if(i<n)right=sqrt(mass(i)*mass(i+1))
       lo=max(left,fmin);hi=min(right,fmax)
       if(hi>lo)number(i)=number(i)+eject*(hi**p-lo**p)/(norm*mass(i))
       left=right
    enddo
  end subroutine

  subroutine dust_collision_kernel(radius,density,nh,temperature,mach,critical_pressure, &
       shatter,mass,kernel,birth,boundary,ierr)
    ! Standalone, compact spherical grains, same-material collisions only.
    ! Hirashita & Aoyama 2019 (1810.07962), eqs15--25. No Coulomb focusing,
    ! sticking threshold, kinetic heating, opacity or live hydro coupling.
    ! 8-point Gauss-Legendre averages K*yield over isotropic relative angles;
    ! NEVER evaluate the nonlinear shattering yield at an averaged speed.
    real(real64),intent(in)::radius(:),density,nh,temperature,mach,critical_pressure
    logical,intent(in)::shatter
    real(real64),intent(out)::mass(:),kernel(:,:),birth(:,:,:),boundary(:,:,:)
    integer,intent(out)::ierr
    real(real64),parameter::mu(8)=[-.9602898564975363d0,-.7966664774136267d0, &
         -.5255324099163290d0,-.1834346424956498d0,.1834346424956498d0, &
         .5255324099163290d0,.7966664774136267d0,.9602898564975363d0]
    real(real64),parameter::weight(8)=[.1012285362903763d0,.2223810344533745d0, &
         .3137066458778873d0,.3626837833783620d0,.3626837833783620d0, &
         .3137066458778873d0,.2223810344533745d0,.1012285362903763d0]/2
    real(real64)::velocity(size(radius)),products(size(radius)),leak(2),v,k,e,phi,ej,qd
    integer::n,i,j,q,t
    ierr=1;mass=0;kernel=0;birth=0;boundary=0;n=size(radius)
    if(n<2.or.size(mass)/=n)return
    if(any(shape(kernel)/=[n,n]).or.any(shape(birth)/=[n,n,n]))return
    if(any(shape(boundary)/=[2,n,n]))return
    if(.not.all(ieee_is_finite(radius)))return
    if(.not.all(ieee_is_finite([density,nh,temperature,mach,critical_pressure])))return
    if(minval(radius)<=0.or.min(density,nh,temperature,critical_pressure)<=0.or.mach<0)return
    if(any(radius(2:)<=radius(:n-1)))return
    mass=4*acos(-1d0)/3*density*radius**3
    velocity=1.1d5*mach**1.5d0*sqrt(radius/1d-5)*(temperature/1d4)**.25d0* &
         nh**(-.25d0)*sqrt(density/3.5d0)
    qd=critical_pressure/(2*density)
    do j=1,n
       do i=1,j
          do q=1,8
             v=sqrt((velocity(i)-velocity(j))**2+2*velocity(i)*velocity(j)*(1-mu(q)))
             k=weight(q)*acos(-1d0)*(radius(i)+radius(j))**2*v
             if(k==0)cycle
             products=0;leak=0
             if(shatter)then
                e=.5d0*mass(i)*(mass(j)/(mass(i)+mass(j)))*v*v
                do t=1,2
                   if(t==1)then
                      phi=e/(mass(i)*qd);ej=mass(i)*(phi/(1+phi))
                      call dust_collision_fragments(mass,mass(i),ej,products,leak)
                   else
                      phi=e/(mass(j)*qd);ej=mass(j)*(phi/(1+phi))
                      call dust_collision_fragments(mass,mass(j),ej,products,leak)
                   endif
                enddo
             else
                call dust_collision_deposit(mass,mass(i)+mass(j),products,leak)
             endif
             kernel(i,j)=kernel(i,j)+k
             birth(:,i,j)=birth(:,i,j)+k*products
             boundary(:,i,j)=boundary(:,i,j)+k*leak
          enddo
          kernel(j,i)=kernel(i,j);birth(:,j,i)=birth(:,i,j);boundary(:,j,i)=boundary(:,i,j)
       enddo
    enddo
    if(.not.all(ieee_is_finite(mass)).or.any(mass<=0))return
    if(.not.all(ieee_is_finite(kernel)).or..not.all(ieee_is_finite(birth)))return
    if(.not.all(ieee_is_finite(boundary)))return
    ierr=0
  end subroutine

  subroutine dust_collision_evolve(mass,kernel,birth,boundary,number,dt,fraction,next,escaped,ierr)
    ! Native reference only. birth = <K * daughter count>; boundary = <K *
    ! daughter mass outside grid>. Count unordered events, self factor 1/2.
    ! Explicit positivity-limited Euler, NO clipping/renormalization. Bound
    ! the net donor rate (including same-bin remnants), avoiding the severe
    ! fictitious stiffness of almost-elastic small-projectile collisions.
    ! All state staged; failure returns original number and zero escaped mass.
    real(real64),intent(in)::mass(:),kernel(:,:),birth(:,:,:),boundary(:,:,:),number(:),dt,fraction
    real(real64),intent(out)::next(:),escaped(2)
    integer,intent(out)::ierr
    real(real64)::trial(size(number)),change(size(number)),net(size(number),size(number),size(number))
    real(real64)::out(2),out_rate(2),event,h,remaining,rate,expected,actual,initial
    integer::n,i,j,it
    escaped=0;ierr=1;n=size(number)
    if(size(next)/=n)then
       next=0;return
    endif
    next=number
    if(n<2.or.size(mass)/=n)return
    if(any(shape(kernel)/=[n,n]).or.any(shape(birth)/=[n,n,n]))return
    if(any(shape(boundary)/=[2,n,n]))return
    if(.not.all(ieee_is_finite(mass)).or..not.all(ieee_is_finite(number)))return
    if(.not.all(ieee_is_finite([dt,fraction])))return
    if(minval(mass)<=0.or.minval(number)<0.or.dt<0.or.fraction<=0.or.fraction>.1d0)return
    if(any(mass(2:)<=mass(:n-1)))return
    if(.not.all(ieee_is_finite(kernel)).or.any(kernel<0))return
    if(.not.all(ieee_is_finite(birth)).or.any(birth<0))return
    if(.not.all(ieee_is_finite(boundary)).or.any(boundary<0))return
    net=birth
    do j=1,n
       do i=1,j
          if(kernel(i,j)/=kernel(j,i))return
          if(any(birth(:,i,j)/=birth(:,j,i)).or.any(boundary(:,i,j)/=boundary(:,j,i)))return
          expected=kernel(i,j)*(mass(i)+mass(j))
          actual=sum(birth(:,i,j)*mass)+sum(boundary(:,i,j))
          if(abs(actual-expected)>2d-12*max(expected,tiny(1d0)))return
          net(i,i,j)=net(i,i,j)-kernel(i,j)
          net(j,i,j)=net(j,i,j)-kernel(i,j)
       enddo
    enddo
    trial=number;out=0;remaining=dt;initial=sum(number*mass)
    do it=1,100000
       if(remaining<=0)exit
       change=0;out_rate=0
       do j=1,n
          do i=1,j
             event=trial(i)*trial(j)
             if(i==j)event=event/2
             if(event==0)cycle
             change=change+event*net(:,i,j)
             out_rate=out_rate+event*boundary(:,i,j)
          enddo
       enddo
       if(.not.all(ieee_is_finite(change)).or..not.all(ieee_is_finite(out_rate)))return
       rate=0
       do i=1,n
          if(change(i)>=0)cycle
          if(trial(i)<=0)return
          rate=max(rate,-change(i)/trial(i))
       enddo
       h=remaining
       if(rate>0)h=min(h,fraction/rate)
       if(h<=0.or.remaining-h==remaining)return
       trial=trial+h*change;out=out+h*out_rate;remaining=remaining-h
       if(.not.all(ieee_is_finite(trial)).or.any(trial<0))return
    enddo
    if(remaining>0)return
    if(.not.all(ieee_is_finite(out)))return
    if(abs(sum(trial*mass)+sum(out)-initial)>2d-10*max(initial,tiny(1d0)))return
    next=trial;escaped=out;ierr=0
  end subroutine

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
         trim(dust_cooling)=='wss09_cie'.or.trim(dust_cooling)=='snrt_hhe_cie_metals'.or. &
         trim(dust_cooling)=='chimes_neq_v1')
    if(trim(dust_cooling)=='chimes_neq_v1')ok=ok.and.dust_two_size_enabled().and.dust_material_composition_enabled()
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
    ok=ok.and.(trim(dust_sublimation)=='none'.or.trim(dust_sublimation)=='gd89_graphite_bulk_v1'.or. &
         trim(dust_sublimation)=='gd89_xu25_olivine_v1'.or.trim(dust_sublimation)=='gd89_xu25_olivine_rt_v1')
    if(trim(dust_sublimation)/='none')ok=ok.and.dust_material_composition_enabled()
    if(dust_sublimation_rt_enabled())ok=ok.and.dust_optics_enabled().and.dust_chimes_enabled()
    ok=ok.and.(trim(dust_optics_model)=='fixed_mix'.or.trim(dust_optics_model)=='d03_transport_v1')
    if(trim(dust_optics_model)=='d03_transport_v1')then
       ok=ok.and.dust_material_composition_enabled()
       ok=ok.and.all(dust_size_radius_cm==[1d-6,1d-5]).and.all(dust_size_density==[2.2d0,3.8d0])
    endif
    if(.not.dust_mass_enabled)ok=ok.and.trim(dust_mass_model)=='bulk_v1'.and.trim(dust_cooling)=='none'
    ok=ok.and.(trim(dust_iron_model)=='none'.or.trim(dust_iron_model)=='fe_electric_compare_v1')
    ok=ok.and.ieee_is_finite(dust_fe_condensation).and.dust_fe_condensation>=0.and.dust_fe_condensation<=1
    ok=ok.and.ieee_is_finite(dust_fe_sticking).and.dust_fe_sticking>=0.and.dust_fe_sticking<=1
    if(trim(dust_iron_model)/='none')then
       ok=ok.and.dust_optics_enabled().and.trim(dust_sublimation)=='none'
       if(.not.dust_chimes_enabled())then
          ! Static comparison only: no CHIMES donor/reconciliation path.
          ok=ok.and.trim(dust_cooling)=='none'
          ok=ok.and.dust_fe_condensation==0.and..not.dust_fe_kinetics
          ok=ok.and.all(dust_condensation==0).and..not.dust_growth.and..not.dust_sputtering
          ok=ok.and..not.dust_coagulation.and..not.dust_shattering.and..not.dust_sn_shocks
       endif
       ok=ok.and.dust_injection_temperature<=dust_fe_max_temperature
       ! No Fe unresolved-shock efficiency has been selected. Never leave
       ! Fe immune while applying an unresolved C/silicate SN-shock model.
       if(dust_fe_kinetics)ok=ok.and..not.dust_sn_shocks
    else
       ok=ok.and.dust_fe_condensation==0.and.dust_fe_sticking==0.and..not.dust_fe_kinetics
    endif
    if(.not.dust_fe_kinetics)ok=ok.and.dust_fe_sticking==0
    ok=ok.and.(trim(dust_pah_model)=='none'.or.trim(dust_pah_model)=='pah_neutral_absolute_v1'.or. &
         trim(dust_pah_model)=='pah_charge_fixed_h_v1'.or.trim(dust_pah_model)=='pah_hydrogen_m13_dl01_v1'.or. &
         trim(dust_pah_model)=='pah_h2_rehydrogenation_v1')
    if(dust_pah_charged())ok=ok.and..not.dust_iron_enabled()
    ok=ok.and.ieee_is_finite(dust_pah_condensation).and.dust_pah_condensation>=0.and.dust_pah_condensation<=1
    if(trim(dust_pah_model)/='none')then
       ok=ok.and.dust_chimes_enabled().and.dust_optics_enabled().and.trim(dust_sublimation)=='none'
    else
       ok=ok.and.dust_pah_condensation==0
    endif
  end function

  logical function dust_pah_enabled() result(enabled)
    enabled=dust_mass_enabled.and.(trim(dust_pah_model)=='pah_neutral_absolute_v1'.or.dust_pah_charged())
  end function

  logical function dust_pah_charged() result(enabled)
    enabled=dust_mass_enabled.and.(trim(dust_pah_model)=='pah_charge_fixed_h_v1'.or.dust_pah_hydrogenated())
  end function

  logical function dust_pah_hydrogenated() result(enabled)
    enabled=dust_mass_enabled.and.(trim(dust_pah_model)=='pah_hydrogen_m13_dl01_v1'.or.dust_pah_h2_enabled())
  end function

  logical function dust_pah_h2_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_pah_model)=='pah_h2_rehydrogenation_v1'
  end function

  integer function dust_pah_charge_size() result(n)
    n=dust_pah_nbin*merge(14,1,dust_pah_hydrogenated())
  end function

  integer function dust_pah_nstate() result(n)
    n=dust_pah_charge_size()*merge(2,1,dust_pah_charged())
  end function

  real(real64) function dust_pah_state_mass(k) result(m)
    integer,intent(in)::k
    m=dust_pah_molecule_g
    ! New H-exchange model uses the CHIMES nuclear-mass convention, so one
    ! molecular H exchanged is exactly one gas H atom. Old modes unchanged.
    if(dust_pah_hydrogenated())m=(288d0+mod((k-1)/dust_pah_nbin,14))*1.66d-24
  end function

  function dust_pah_inventory(mass) result(hc)
    real(real64),intent(in)::mass(:)
    real(real64)::hc(2),fraction
    integer::k,nh
    hc=sum(mass)*dust_pah_hc
    if(.not.dust_pah_hydrogenated())return
    hc=0
    do k=1,size(mass)
       nh=mod((k-1)/dust_pah_nbin,14);fraction=real(nh,real64)/(288+nh)
       hc(1)=hc(1)+mass(k)*fraction;hc(2)=hc(2)+mass(k)*(1-fraction)
    enddo
  end function

  real(real64) function dust_pah_solid_charge(mass,carrier_mass_g) result(q)
    ! Signed charge in the caller's m_H*n carrier units. The carrier mass
    ! convention need not equal the proton mass used by the PAH molecule.
    real(real64),intent(in)::mass(:),carrier_mass_g
    integer::k
    q=0
    if(dust_pah_hydrogenated())then
       do k=dust_pah_charge_size()+1,size(mass)
          q=q+mass(k)*(carrier_mass_g/dust_pah_state_mass(k))
       enddo
    else if(dust_pah_charged())then
       q=sum(mass(dust_pah_nbin+1:))*(carrier_mass_g/dust_pah_molecule_g)
    endif
  end function

  subroutine dust_pah_condense(hydrogen,carbon,graphite,mass,ierr)
    ! Same mass units throughout. Input is NON-Ia ejecta; graphite has
    ! already reserved its carbon. PAH hydrogen is not a metal reservoir.
    real(real64),intent(in)::hydrogen,carbon,graphite
    real(real64),intent(out)::mass
    integer,intent(out)::ierr
    mass=0;ierr=1
    if(.not.all(ieee_is_finite([hydrogen,carbon,graphite,dust_pah_condensation])))return
    if(min(hydrogen,carbon,graphite,dust_pah_condensation)<0.or.dust_pah_condensation>1)return
    if(graphite>carbon+128*epsilon(1d0)*max(carbon,tiny(1d0)))return
    mass=min(dust_pah_condensation*max(0d0,carbon-graphite)/dust_pah_hc(2),hydrogen/dust_pah_hc(1))
    if(.not.ieee_is_finite(mass))return
    ierr=0
  end subroutine

  logical function dust_iron_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_iron_model)=='fe_electric_compare_v1'
  end function

  logical function dust_atomic_cooling_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_cooling)=='snrt_hhe_cie_metals'
  end function

  logical function dust_composition_enabled() result(enabled)
    enabled=dust_mass_enabled.and.(trim(dust_mass_model)=='carbon_olivine_v1'.or. &
         trim(dust_mass_model)=='carbon_olivine_2size_v1')
  end function

  logical function dust_chimes_enabled() result(enabled)
    enabled=dust_mass_enabled.and.trim(dust_cooling)=='chimes_neq_v1'
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

  logical function dust_sublimation_enabled() result(enabled)
    enabled=dust_material_composition_enabled().and.(trim(dust_sublimation)=='gd89_graphite_bulk_v1'.or. &
         trim(dust_sublimation)=='gd89_xu25_olivine_v1'.or.trim(dust_sublimation)=='gd89_xu25_olivine_rt_v1')
  end function

  logical function dust_silicate_sublimation_enabled() result(enabled)
    enabled=dust_material_composition_enabled().and.(trim(dust_sublimation)=='gd89_xu25_olivine_v1'.or. &
         trim(dust_sublimation)=='gd89_xu25_olivine_rt_v1')
  end function

  logical function dust_sublimation_rt_enabled() result(enabled)
    enabled=dust_material_composition_enabled().and.trim(dust_sublimation)=='gd89_xu25_olivine_rt_v1'
  end function

  function dust_sublimation_identity() result(v)
    real(real64)::v(6)
    ! Version 1: fixed-radius two-size BE loss, vacuum, bulk B, 2 kT escaping atoms.
    v=[1d0,dust_carbon_nu,dust_carbon_binding_k,dust_carbon_atom,2d0,dust_carbon_latent]
    ! v3: coupled radiation/phase, BE step doubling at relative tolerance 1e-4;
    ! incompatible with the earlier un-subcycled v2 comparison checkpoints.
    if(dust_sublimation_rt_enabled())v(1)=3d0
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

  subroutine dust_size_exchange(rho,nh,temperature,material,dt,bins,ierr,phase_momentum,mixing_heat)
    ! Exact donor-quadratic transfer, dD_donor/dt=-k D_donor^2.
    ! Collisions within each composition; no creation/destruction of metals.
    real(real64),intent(in)::rho,nh,temperature,dt
    integer,intent(in)::material
    real(real64),intent(inout)::bins(2)
    integer,intent(out)::ierr
    real(real64),optional,intent(inout)::phase_momentum(3,2),mixing_heat
    real(real64)::rate,donor,transfer,x,p
    real(real64)::pair_mass(2),pair_next(2),pair_p(3,2),heat,flows(2)
    integer::j,status
    ierr=1
    if(present(phase_momentum).neqv.present(mixing_heat))return
    if(.not.all(ieee_is_finite([rho,nh,temperature,dt,bins])))return
    if(rho<=0.or.min(nh,temperature,dt,minval(bins))<0.or.material<1.or.material>2)return
    if(present(phase_momentum))then
       if(.not.all(ieee_is_finite(phase_momentum)))return
       do j=1,2
          if(bins(j)==0.and.any(phase_momentum(:,j)/=0))return
       enddo
    endif
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
       if(present(mixing_heat))mixing_heat=0
       ierr=0;return
    endif
    donor=bins(j)/(1+x)
    transfer=bins(j)-donor
    if(present(phase_momentum))then
       pair_mass=bins;pair_next=bins;pair_next(j)=donor;pair_next(3-j)=bins(3-j)+transfer
       flows=0;flows(j)=transfer;pair_p=phase_momentum;heat=0
       call dust_phase_pair(pair_mass,pair_next,flows,pair_p,heat,status)
       if(status/=0)return
       phase_momentum=pair_p;mixing_heat=heat
    endif
    bins(j)=donor;bins(3-j)=bins(3-j)+transfer;ierr=0
  end subroutine

  subroutine dust_size_step(rho,metal,elements,bins,temperature,scale_d,dt,next,ierr, &
       phase_momentum,mixing_heat,fixed_solid_density)
    ! Symmetric, positivity-preserving splits of shared-reservoir growth,
    ! size-dependent thermal erosion and exact coagulation/shattering flows.
    ! Growth retains the explicitly configured effective accreting atom mass;
    ! sputtering is the shared Tsai-Mathews reference, NOT species Hu fits.
    real(real64),intent(in)::rho,metal,elements(11),bins(4),temperature,scale_d,dt
    real(real64),intent(out)::next(4)
    integer,intent(out)::ierr
    ! Optional pair: absolute MOMENTUM DENSITIES (xyz, gas=0, bins=1:4),
    ! evolved in place, and this call's dissipated kinetic energy density.
    ! Both are unchanged on failure; next then equals the input bins.
    ! rho and the frozen scalar rate laws retain their existing meaning.
    ! Gas inertia = rho - fixed_solid_density - sum(current bins). Fixed
    ! solids do not exchange mass/momentum here (the reserved wrapper sets it).
    ! Away from admitted element roundoff, mass numerics are unchanged. The
    ! coupled path retains actual donor mass at that boundary and preserves
    ! zero-rate/time identity instead of the scalar capacity-roundoff clip.
    ! Momentum uses gross
    ! directed transfers and conservative BE at EACH existing substep: first
    ! order, not an exact coupled solution. Heat is NOT an extra total-E source.
    real(real64),optional,intent(inout)::phase_momentum(3,0:4),mixing_heat
    real(real64),optional,intent(in)::fixed_solid_density
    real(real64)::gas(11),grains(2),cap(2),aa(2,2),bb(2,2),a,b,h,available,nh,steps,updated
    real(real64)::trial(4),mom(3,0:4),fixed,heat,q,pmass(2),pnext(2),pair_p(3,2),flows(2),rg
    integer::i,j,k,m,first,status,ns,seq(3)
    logical::coupled
    next=bins;ierr=1
    coupled=present(phase_momentum)
    if(coupled.neqv.present(mixing_heat))return
    fixed=0
    if(present(fixed_solid_density))fixed=fixed_solid_density
    if(.not.ieee_is_finite(fixed).or.fixed<0)return
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([rho,metal,bins,temperature,scale_d,dt])))return
    if(rho<=0.or.scale_d<=0.or.min(metal,temperature,dt,minval(bins))<0.or.metal>rho)return
    trial=bins;heat=0
    if(coupled)then
       mom=phase_momentum;rg=rho-fixed-sum(bins)
       if(rg<=0.or..not.all(ieee_is_finite(mom)))return
       do j=1,4
          if(bins(j)==0.and.any(mom(:,j)/=0))return
       enddo
    endif
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
          if(coupled)then
             call dust_size_exchange(rho,nh,temperature,i,h/2,trial(first:first+1),status, &
                  mom(:,first:first+1),q)
             if(status==0)heat=heat+q
          else
             call dust_size_exchange(rho,nh,temperature,i,h/2,trial(first:first+1),status)
          endif
          if(status/=0)return
          do m=1,3
             j=seq(m)
             available=max(0d0,cap(i)-trial(first+2-j))
             if(coupled.and.available<trial(first+j-1))then
                ! Reuse the EXISTING elementwise roundoff admission, on the
                ! current bins after any size exchange/previous rate split.
                ! Only an admitted boundary may retain its actual donor mass
                ! as effective capacity; no new budget tolerance or mass floor.
                grains=[sum(trial(1:2)),sum(trial(3:4))]
                call dust_gas_elements(elements,grains,gas,status)
                if(status/=0)return
                available=trial(first+j-1)
             endif
             a=0;if(cap(i)>0)a=aa(j,i)*available/cap(i)
             b=h/2;if(m==2)b=h
             ! The element budget above admits representational roundoff.
             ! At zero rate/time, min(bin,available) would turn that into a
             ! spurious mass loss with no directed process/donor. Preserve
             ! exact identity here; retain the legacy no-momentum arithmetic.
             if(coupled.and.(b==0.or.(a==0.and.bb(j,i)==0)))cycle
             if(coupled)then
                call dust_mass_step(available,trial(first+j-1),b,a,bb(j,i),updated,status)
             else
                call dust_mass_step(available,min(trial(first+j-1),available),b,a,bb(j,i),updated,status)
             endif
             if(status/=0)return
             if(coupled)then
                ! Do not reinterpret the final mass difference as a donor.
                ! The logistic ODE may have nonzero opposing flows at equilibrium.
                call dust_phase_gross(available,trial(first+j-1),b,a,bb(j,i),updated,flows,status)
                if(status/=0)return
                pmass=[rho-fixed-sum(trial),trial(first+j-1)]
                pnext=[pmass(1)-(updated-pmass(2)),updated]
                if(pnext(1)<=0)return
                pair_p(:,1)=mom(:,0);pair_p(:,2)=mom(:,first+j-1)
                call dust_phase_pair(pmass,pnext,flows,pair_p,q,status)
                if(status/=0)return
                mom(:,0)=pair_p(:,1);mom(:,first+j-1)=pair_p(:,2);heat=heat+q
             endif
             trial(first+j-1)=updated
          enddo
          if(coupled)then
             call dust_size_exchange(rho,nh,temperature,i,h/2,trial(first:first+1),status, &
                  mom(:,first:first+1),q)
             if(status==0)heat=heat+q
          else
             call dust_size_exchange(rho,nh,temperature,i,h/2,trial(first:first+1),status)
          endif
          if(status/=0)return
       enddo
    enddo
    grains=[sum(trial(1:2)),sum(trial(3:4))]
    call dust_gas_elements(elements,grains,gas,status)
    if(status/=0)return
    if(coupled)then
       if(.not.ieee_is_finite(heat))return
       phase_momentum=mom;mixing_heat=heat
    endif
    next=trial;ierr=0
  end subroutine

  subroutine dust_size_step_reserved_iron(rho,metal,elements,bins,metallic_iron,temperature,scale_d,dt,next,ierr, &
       pah_hc,phase_momentum,mixing_heat)
    ! Evolve the existing four C/olivine bins while a separate Fe reservoir
    ! is held fixed. Both metal capacity and elemental Fe must be reserved.
    ! Optional PAH H/C inventories are likewise held fixed, not reaccreted
    ! as carbon grains or reused as gas hydrogen. H is not a metal budget.
    ! This supplies no Fe growth/sputtering law or stellar condensation model.
    real(real64),intent(in)::rho,metal,elements(11),bins(4),metallic_iron,temperature,scale_d,dt
    real(real64),optional,intent(in)::pah_hc(2) ! H/C mass densities, not number ratio or grain radius
    real(real64),intent(out)::next(4)
    integer,intent(out)::ierr
    real(real64),optional,intent(inout)::phase_momentum(3,0:4),mixing_heat
    real(real64)::available(11),gas(11),trial(4),grains(2),available_metal
    real(real64)::mom(3,0:4),heat,fixed
    integer::status
    next=bins;ierr=1
    if(present(phase_momentum).neqv.present(mixing_heat))return
    if(.not.all(ieee_is_finite([metal,metallic_iron])))return
    if(metallic_iron<0.or.metallic_iron>metal)return
    grains=[sum(bins(1:2)),sum(bins(3:4))]
    call dust_gas_elements(elements,grains,gas,status,metallic_iron,pah_hc)
    if(status/=0)return
    available=elements
    ! A representational boundary excess already accepted by dust_gas_elements
    ! is at most its elementwise roundoff tolerance; don't create negative Fe.
    available(11)=max(0d0,elements(11)-metallic_iron)
    available_metal=metal-metallic_iron
    if(present(pah_hc))then
       available(1)=max(0d0,elements(1)-pah_hc(1))
       available(3)=max(0d0,elements(3)-pah_hc(2))
       ! PAH hydrogen is baryonic mass but NOT metallicity. Reserve C from
       ! the metal capacity and H from the gas collision/chemical inventory.
       if(pah_hc(2)>available_metal)return
       available_metal=available_metal-pah_hc(2)
    endif
    if(present(phase_momentum))then
       fixed=metallic_iron
       if(present(pah_hc))fixed=fixed+sum(pah_hc)
       mom=phase_momentum;heat=0
       call dust_size_step(rho,available_metal,available,bins,temperature,scale_d,dt,trial,status, &
            mom,heat,fixed)
    else
       call dust_size_step(rho,available_metal,available,bins,temperature,scale_d,dt,trial,status)
    endif
    if(status/=0)return
    grains=[sum(trial(1:2)),sum(trial(3:4))]
    call dust_gas_elements(elements,grains,gas,status,metallic_iron,pah_hc)
    if(status/=0)return
    if(present(phase_momentum))then
       phase_momentum=mom;mixing_heat=heat
    endif
    next=trial;ierr=0
  end subroutine

  subroutine dust_fe_sputtering_yield(temperature,y,ierr)
    ! C2=1 (resolved density); fixed low-Z projectile mixture as in the
    ! source fit, not an arbitrary-metallicity sputtering calculation.
    ! Below 1e4 K we explicitly neglect thermal erosion, not extrapolate a
    ! high-T polynomial. Above 1e9 K reject (bounded subset of source Fig2).
    real(real64),intent(in)::temperature
    real(real64),intent(out)::y ! positive |da/dt| / n_H, cm^4 / s
    integer,intent(out)::ierr
    real(real64)::x,logy
    integer::j
    y=0;ierr=1
    if(.not.ieee_is_finite(temperature))return
    if(temperature<0.or.temperature>dust_fe_sputter_max_t)return
    if(temperature>=dust_fe_sputter_min_t)then
       x=log10(temperature);logy=dust_fe_sputter_coeff(5)
       do j=4,0,-1
          logy=dust_fe_sputter_coeff(j)+x*logy
       enddo
       y=10d0**logy*1d-4/31557600d0
    endif
    ierr=0
  end subroutine

  subroutine dust_fe_kinetics_step(rho,elements,grains,iron,temperature,scale_d,dt,radius,density,next,ierr, &
       pah_hc,phase_momentum,mixing_heat)
    ! Fixed representative-size closure, as for the existing two-size model.
    ! dD_j/dt = 3*S*v_Fe/(4*rho_s*a_j) * D_j * rho_Fe,gas.
    ! v_Fe=sqrt(8*k*T/(pi*56*m_p)); mass-number convention matches CHIMES.
    ! Geometric accretion plus Fe-specific thermal sputtering. No Coulomb
    ! focusing, charge evolution, nucleation, coagulation or nonthermal erosion.
    ! The caller supplies the SAME radius/density as the active Fe optics.
    ! Olivine and PAH are reserved; both Fe bins share the remaining Fe.
    ! Optional absolute momenta: gas=0, Fe bins=1:2. Other solids are fixed.
    ! Symmetric split mass evolution; donor momentum BE is first order.
    real(real64),intent(in)::rho,elements(11),grains(2),iron(2),temperature,scale_d,dt,radius(2),density
    real(real64),intent(out)::next(2)
    integer,intent(out)::ierr
    real(real64),optional,intent(in)::pah_hc(2)
    real(real64),optional,intent(inout)::phase_momentum(3,0:2),mixing_heat
    real(real64)::gas(11),solid_pah(2),trial(2),cap,speed,rate(2),loss(2),yield,steps,h,available,updated,aa
    real(real64)::mom(3,0:2),pmass(2),pnext(2),pair_p(3,2),flows(2),heat,q,fixed,step_time
    integer::status,ns,k,m,j,seq(3)
    logical::coupled
    next=iron;ierr=1;solid_pah=0;heat=0
    coupled=present(phase_momentum)
    if(coupled.neqv.present(mixing_heat))return
    if(present(pah_hc))solid_pah=pah_hc
    if(.not.dust_mass_parameters_ok())return
    if(.not.all(ieee_is_finite([rho,iron,temperature,scale_d,dt,radius,density])))return
    if(min(rho,scale_d,density,minval(radius))<=0.or.min(temperature,dt,minval(iron))<0)return
    call dust_gas_elements(elements,grains,gas,status,sum(iron),solid_pah)
    if(status/=0)return
    fixed=sum(grains)+sum(solid_pah)
    if(rho-fixed-sum(iron)<=0)return
    if(coupled)then
       mom=phase_momentum
       if(.not.all(ieee_is_finite(mom)))return
       do j=1,2
          if(iron(j)==0.and.any(mom(:,j)/=0))return
       enddo
    endif
    ! Exact identity includes admitted element-budget roundoff.
    if(dt==0.or..not.dust_fe_kinetics.or.sum(iron)==0)then
       if(coupled)mixing_heat=0
       ierr=0;return
    endif
    cap=gas(11)+sum(iron)
    rate=0;loss=0
    if(dust_growth.and.temperature>0.and.temperature<=dust_growth_max_temperature)then
       speed=sqrt(8*dust_kb*temperature/(acos(-1d0)*56*dust_mp))
       rate=.75d0*dust_fe_sticking*cap*scale_d*speed/(density*radius)
    endif
    if(dust_sputtering)then
       call dust_fe_sputtering_yield(temperature,yield,status)
       if(status/=0)return
       loss=3*yield*gas(1)*scale_d/dust_mp/radius
    endif
    steps=maxval(rate+loss)*dt/.1d0
    if(.not.all(ieee_is_finite([rate,loss])).or..not.ieee_is_finite(steps))return
    if(steps>4096)return
    ns=max(1,ceiling(steps));h=dt/ns;trial=iron;seq=[1,2,1]
    do k=1,ns
       do m=1,3
          j=seq(m);available=max(trial(j),cap-trial(3-j))
          step_time=h/2;if(m==2)step_time=h
          aa=rate(j)*available/cap
          if(aa==0.and.loss(j)==0)cycle
          call dust_mass_step(available,trial(j),step_time,aa,loss(j),updated,status)
          if(status/=0)return
          if(coupled)then
             call dust_phase_gross(available,trial(j),step_time,aa,loss(j),updated,flows,status)
             if(status/=0)return
             pmass=[rho-fixed-sum(trial),trial(j)]
             pnext=[pmass(1)-(updated-trial(j)),updated]
             if(pnext(1)<=0)return
             pair_p(:,1)=mom(:,0);pair_p(:,2)=mom(:,j)
             call dust_phase_pair(pmass,pnext,flows,pair_p,q,status)
             if(status/=0)return
             mom(:,0)=pair_p(:,1);mom(:,j)=pair_p(:,2);heat=heat+q
          endif
          trial(j)=updated
       enddo
    enddo
    call dust_gas_elements(elements,grains,gas,status,sum(trial),solid_pah)
    if(status/=0.or..not.ieee_is_finite(heat))return
    if(coupled)then
       phase_momentum=mom;mixing_heat=heat
    endif
    next=trial;ierr=0
  end subroutine

  subroutine dust_phase_gross(capacity,dust,dt,a,b,updated,flows,ierr)
    ! Gross integrated gas->grain and grain->gas transfers for the frozen
    ! logistic substep. Only used inside the existing (a+b)*dt <= 0.1 bound.
    ! Integral D dt = D0*phi*log(1+y)/y, y=a*D0*phi/Z,
    ! phi=(exp((a-b)*dt)-1)/(a-b). Destruction=b*integral; growth follows
    ! from mass balance. This retains turnover when updated == dust.
    real(real64),intent(in)::capacity,dust,dt,a,b,updated
    real(real64),intent(out)::flows(2)
    integer,intent(out)::ierr
    real(real64)::x,phi,y,ratio,tol
    ierr=1;flows=0
    if(.not.all(ieee_is_finite([capacity,dust,dt,a,b,updated])))return
    if(min(capacity,dust,dt,a,b,updated)<0.or.dust>capacity)return
    if((a+b)*dt>.100000000001d0)return
    if(dust==0.or.dt==0.or.(a==0.and.b==0))then
       if(updated/=dust)return
       ierr=0;return
    endif
    if(a==0)then
       flows(2)=dust-updated
    else if(b==0)then
       flows(1)=updated-dust
    else
       x=(a-b)*dt
       if(abs(x)<1d-5)then
          phi=dt*(1+x/2+x*x/6+x*x*x/24)
       else
          phi=(exp(x)-1)/(a-b)
       endif
       y=a*(dust/capacity)*phi
       if(abs(y)<1d-4)then
          ratio=1-y/2+y*y/3-y**3/4+y**4/5
       else
          ratio=log(1+y)/y
       endif
       flows(2)=b*(dust*phi*ratio)
       flows(1)=(updated-dust)+flows(2)
    endif
    tol=128*epsilon(1d0)*max(dust,updated,maxval(abs(flows)),tiny(1d0))
    if(any(.not.ieee_is_finite(flows)).or.any(flows < -tol))return
    flows=max(flows,0d0);ierr=0
  end subroutine

  subroutine dust_phase_pair(mass,next_mass,flows,momentum,mixing_heat,ierr)
    ! Pairwise conservative backward-Euler momentum coupling. flows(1) is
    ! 1->2 and flows(2) is 2->1, not signed NET mass. Solve donor velocities
    ! implicitly with the known final masses. For one-way transfer this
    ! leaves the donor velocity unchanged and mixes its mass into receiver.
    ! Work in the barycentric frame to measure heat without subtracting large
    ! bulk kinetic energies. No material/latent/CR energy is changed here.
    real(real64),intent(in)::mass(2),next_mass(2),flows(2)
    real(real64),intent(inout)::momentum(3,2),mixing_heat
    integer,intent(out)::ierr
    real(real64)::total,vel(3,2),bulk(3),relative(3),fraction,af,bf,heat
    real(real64)::trial(3,2),psum(3),tol
    integer::j
    ierr=1
    if(.not.all(ieee_is_finite([mass,next_mass,flows,reshape(momentum,[6])])))return
    if(min(minval(mass),minval(next_mass),minval(flows))<0)return
    total=sum(mass)
    if(.not.ieee_is_finite(total).or.total<=0)return
    tol=128*epsilon(1d0)*max(total,maxval(flows),tiny(1d0))
    if(abs(sum(next_mass)-total)>tol)return
    if(abs(next_mass(1)-(mass(1)-flows(1)+flows(2)))>tol)return
    vel=0
    do j=1,2
       if(mass(j)==0)then
          if(any(momentum(:,j)/=0).or.flows(j)/=0)return
       else
          vel(:,j)=momentum(:,j)/mass(j)
       endif
    enddo
    if(all(flows==0))then
       if(any(next_mass/=mass))return
       mixing_heat=0;ierr=0;return
    endif
    psum=sum(momentum,dim=2);bulk=psum/total
    heat=0;trial=0
    if(any(mass==0))then
       trial(:,1)=next_mass(1)*bulk;trial(:,2)=next_mass(2)*bulk
    else
       relative=vel(:,1)-vel(:,2)
       fraction=1/(1+flows(1)/mass(2)+flows(2)/mass(1))
       ! BE energy identity: 2Q = sum_i m_i_old |dv_i|^2
       !                        + (flow12+flow21)*|v1_new-v2_new|^2.
       ! Positive terms avoid subtracting nearly equal reduced kinetic
       ! energies, which can falsely reject tiny but valid mass transfers.
       af=(flows(1)/mass(2))*fraction;bf=(flows(2)/mass(1))*fraction
       heat=.5d0*(flows(1)*fraction*(fraction+af)+flows(2)*fraction*(fraction+bf))*sum(relative**2)
       trial(:,1)=next_mass(1)*(bulk+(next_mass(2)/total)*fraction*relative)
       trial(:,2)=next_mass(2)*(bulk-(next_mass(1)/total)*fraction*relative)
    endif
    ! Put only roundoff from the solve into the more massive survivor.
    j=maxloc(next_mass,dim=1);trial(:,j)=psum-trial(:,3-j)
    if(.not.all(ieee_is_finite(trial)).or..not.ieee_is_finite(heat))return
    momentum=trial;mixing_heat=heat;ierr=0
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

  subroutine dust_gas_elements(elements,grains,gas,ierr,metallic_iron,pah_hc)
    real(real64),intent(in)::elements(dust_ne),grains(dust_nc)
    real(real64),intent(out)::gas(dust_ne)
    integer,intent(out)::ierr
    real(real64),optional,intent(in)::metallic_iron
    real(real64),optional,intent(in)::pah_hc(2)
    real(real64)::locked(dust_ne),tol(dust_ne)
    gas=0;ierr=1
    if(.not.all(ieee_is_finite(elements)).or..not.all(ieee_is_finite(grains)))return
    if(any(elements<0).or.any(grains<0))return
    locked=olivine_fraction*grains(2);locked(3)=grains(1)
    ! Separate metallic Fe is additional solid Fe, not a relabelling of the
    ! Fe already in MgFeSiO4. Omitted argument preserves the old arithmetic.
    if(present(metallic_iron))then
       if(.not.ieee_is_finite(metallic_iron))return
       if(metallic_iron<0)return
       locked(11)=locked(11)+metallic_iron
       if(.not.ieee_is_finite(locked(11)))return
    endif
    if(present(pah_hc))then
       if(.not.all(ieee_is_finite(pah_hc)).or.any(pah_hc<0))return
       locked(1)=locked(1)+pah_hc(1);locked(3)=locked(3)+pah_hc(2)
       if(.not.all(ieee_is_finite(locked)))return
    endif
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

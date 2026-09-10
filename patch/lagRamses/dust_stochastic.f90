! Native stochastic grain-energy population primitive. No global workspace,
! file IO, runtime Python, or live selector. The caller supplies physical
! transition rates and must establish energy-grid coverage and relaxation.
module dust_stochastic_physics
  use, intrinsic :: iso_fortran_env, only: real64,int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public::dust_stochastic_equilibrium,dust_stochastic_photon_rates
  public::dust_stochastic_evolve
  public::dust_pah_modes,dust_vibrational_curve
contains
  subroutine dust_stochastic_evolve(energy,up,down,old_population,dt,population,absorbed,emitted,ierr, &
       sink_rate,removed)
    ! Finite-time backward Euler of the thermal-continuous master equation.
    ! Unlike equilibrium(), this preserves the incoming grain population;
    ! there is no assumption that relaxation is faster than the hydro step.
    ! The population may be a normalized probability or transported grain
    ! NUMBER DENSITIES per energy bin. Its total is preserved, not normalized
    ! to one; zero population is a genuine empty cell. Rates are frozen in dt.
    ! absorbed/emitted are integrated erg/grain or erg/cm3, respectively.
    ! The caller still owns photon supply, grid overflow,
    ! molecular destruction/charging and the physical emission spectrum.
    !
    ! Positive-rate elimination of the lower-Hessenberg M-matrix. Store its
    ! loss-to-reservoir separately instead of subtracting large diagonal
    ! terms: the latter loses the identity term at dt*rate >> 1. Each
    ! eliminated state feeds only the next state upwards in the matrix, so
    ! the solve costs O(N^2), without a stiff explicit subcycle loop.
    real(real64),intent(in)::energy(:),up(:,:),down(:),old_population(:),dt
    real(real64),intent(inout)::population(:),absorbed,emitted
    integer,intent(out)::ierr
    ! Optional outward reaction channel. Return its integrated state-resolved
    ! counts so the enclosing H-state network can populate daughters with
    ! their proper excitation and binding energies; never discard molecules.
    real(real64),optional,intent(in)::sink_rate(:)
    real(real64),optional,intent(inout)::removed(:)
    real(real64)::flow(size(energy),size(energy)),leak(size(energy)),rhs(size(energy))
    real(real64)::pivot(size(energy)),p(size(energy)),fraction,pa,pe,change,scale,total
    real(real64)::loss(size(energy)),lost_energy
    integer::n,i,j,k
    logical::downward
    ierr=1;n=size(energy)
    if(n<2.or.size(down)/=n.or.size(old_population)/=n.or.size(population)/=n)return
    if(any(shape(up)/=[n,n]))return
    if(.not.all(ieee_is_finite(energy)).or.any(energy<0))return
    if(any(energy(2:)<=energy(:n-1)))return
    if(.not.all(ieee_is_finite(up)).or.any(up<0))return
    if(.not.all(ieee_is_finite(down)).or.down(1)/=0.or.any(down<0))return
    if(.not.all(ieee_is_finite(old_population)).or.any(old_population<0))return
    total=sum(old_population)
    if(.not.ieee_is_finite(total))return
    if(.not.ieee_is_finite(dt).or.dt<0)return
    if(present(sink_rate).neqv.present(removed))return
    loss=0
    if(present(sink_rate))then
       if(size(sink_rate)/=n.or.size(removed)/=n)return
       if(any(.not.ieee_is_finite(sink_rate)).or.any(sink_rate<0))return
       loss=dt*sink_rate
       if(any(.not.ieee_is_finite(loss)))return
    endif
    do i=1,n
       if(any(up(:i,i)/=0))return
    enddo
    downward=.false.
    if(present(sink_rate))downward=all(up==0)
    if(downward)then
       ! H-state relaxation without incident captures has no upward edge.
       ! Solve the identical triangular BE system in O(N), not dense O(N^2).
       rhs=old_population
       do k=n,1,-1
          fraction=dt*down(k);pivot(k)=1+loss(k)+fraction
          if(.not.ieee_is_finite(pivot(k)))return
          p(k)=rhs(k)/pivot(k)
          if(k>1)rhs(k-1)=rhs(k-1)+(fraction/pivot(k))*rhs(k)
       enddo
    else
    flow=dt*up;leak=1+loss;rhs=old_population
    do i=2,n
       flow(i-1,i)=dt*down(i)
    enddo
    if(any(.not.ieee_is_finite(flow)))return
    do k=1,n-1
       pivot(k)=leak(k)+sum(flow(k+1:n,k))
       if(.not.ieee_is_finite(pivot(k)).or.pivot(k)<=0)return
       ! Do not form diagonal -= offdiag*offdiag/pivot. The positive leak
       ! carries the original probability source even in the stiff limit.
       leak(k+1)=leak(k+1)+(leak(k)/pivot(k))*flow(k,k+1)
       do i=k+1,n
          fraction=flow(i,k)/pivot(k)
          rhs(i)=rhs(i)+fraction*rhs(k)
          if(i>k+1)flow(i,k+1)=flow(i,k+1)+fraction*flow(k,k+1)
       enddo
    enddo
    if(any(.not.ieee_is_finite(leak)).or.any(.not.ieee_is_finite(rhs)))return
    p(n)=rhs(n)/leak(n)
    do k=n-1,1,-1
       p(k)=rhs(k)/pivot(k)+(flow(k,k+1)/pivot(k))*p(k+1)
    enddo
    endif
    if(any(.not.ieee_is_finite(p)).or.any(p<0))return
    loss=loss*p
    if(any(.not.ieee_is_finite(loss)))return
    if(abs(sum(p)+sum(loss)-total)>512*epsilon(1d0)*n*max(total,tiny(1d0)))return
    pa=0;pe=0
    do i=1,n
       if(.not.downward)then
          do j=i+1,n
             pa=pa+(dt*up(j,i))*p(i)*(energy(j)-energy(i))
          enddo
       endif
       if(i>1)pe=pe+(dt*down(i))*p(i)*(energy(i)-energy(i-1))
    enddo
    change=dot_product(energy,p-old_population)
    lost_energy=dot_product(energy,loss)
    scale=max(pa,pe,abs(change),lost_energy,tiny(1d0))
    if(.not.all(ieee_is_finite([pa,pe,change,scale])))return
    if(abs(change-pa+pe+lost_energy)>1024*epsilon(1d0)*n*scale)return
    population=p;absorbed=pa;emitted=pe;ierr=0
    if(present(removed))removed=loss
  end subroutine

  subroutine dust_pah_modes(ncarbon,nhydrogen,mode_energy,ierr)
    ! DL01 eqs. 4-7 and section II.1: generic PAH normal modes, NOT a
    ! particular molecule's measured spectrum. Carbon skeleton and C-H
    ! bending/stretching are distinct. Output hbar*omega in erg per mode.
    ! Atom counts are explicit; do not silently impose an H/C relation on
    ! independently conserved H and C. Use the source's wavenumbers rather
    ! than its rounded conversions to Debye temperatures.
    integer,intent(in)::ncarbon,nhydrogen
    real(real64),intent(out)::mode_energy(:)
    integer,intent(out)::ierr
    real(real64),parameter::hc=6.62607015d-27*2.99792458d10
    real(real64),parameter::cc(2)=[600d0,1740d0],ch(3)=[886d0,1161d0,3030d0]
    real(real64)::beta,shift
    integer::s,j,n,k
    ierr=1;mode_energy=0
    if(ncarbon<3.or.nhydrogen<0)return
    if(size(mode_energy,kind=int64)/=3_int64*(int(ncarbon,int64)+nhydrogen-2))return
    k=0
    do s=1,2
       n=s*(ncarbon-2);beta=0
       if(ncarbon>54.and.ncarbon<=102)beta=((real(ncarbon,real64)-54)/52)/(2d0*n-1)
       if(ncarbon>102)beta=((real(ncarbon,real64)-2)/52*(102d0/ncarbon)**(2d0/3)-1)/(2d0*n-1)
       do j=1,n
          shift=.5d0
          if(j==2.or.j==3)shift=1
          k=k+1;mode_energy(k)=hc*cc(s)*sqrt((1-beta)*(j-shift)/n+beta)
       enddo
    enddo
    do s=1,3
       mode_energy(k+1:k+nhydrogen)=hc*ch(s);k=k+nhydrogen
    enddo
    if(k/=size(mode_energy))return
    ierr=0
  end subroutine

  subroutine dust_vibrational_curve(mode_energy,temperature,energy,capacity,ierr)
    ! Harmonic canonical expectation and derivative, DL01 eqs. 2-3.
    ! Excludes zero-point energy. U [erg/grain], C [erg/K/grain], T [K].
    ! This supplies an excitation curve, not proof of molecular survival
    ! or the low-excitation microcanonical temperature of the exact model.
    real(real64),intent(in)::mode_energy(:),temperature(:)
    real(real64),intent(out)::energy(:),capacity(:)
    integer,intent(out)::ierr
    real(real64),parameter::kb=1.380649d-16
    real(real64)::u(size(temperature)),cv(size(temperature)),x,q,logx
    integer::i,j
    ierr=1;energy=0;capacity=0
    if(size(mode_energy)<1.or.size(energy)/=size(temperature).or.size(capacity)/=size(temperature))return
    if(.not.all(ieee_is_finite(mode_energy)).or.any(mode_energy<=0))return
    if(.not.all(ieee_is_finite(temperature)).or.any(temperature<0))return
    u=0;cv=0
    do i=1,size(temperature)
       if(temperature(i)==0)cycle
       do j=1,size(mode_energy)
          logx=log(mode_energy(j))-log(kb)-log(temperature(i))
          if(logx>log(700d0))cycle
          x=exp(logx)
          if(x<1d-3)then
             u(i)=u(i)+kb*temperature(i)*(1-x/2+x*x/12-x**4/720)
             cv(i)=cv(i)+kb*(1-x*x/12+x**4/240)
          else
             q=exp(-x)
             u(i)=u(i)+mode_energy(j)*q/(1-q)
             cv(i)=cv(i)+kb*x*x*q/(1-q)**2
          endif
       enddo
    enddo
    if(.not.all(ieee_is_finite(u)).or..not.all(ieee_is_finite(cv)))return
    energy=u;capacity=cv;ierr=0
  end subroutine

  subroutine dust_stochastic_photon_rates(energy,photon_energy,photon_rate,up,overflow_rate,overflow_power,ierr)
    ! Monochromatic source-group closure: photon_rate is Cabs*c_hat*n_gamma
    ! per grain per group [s^-1]; all energies are erg, not eV. Distribute a
    ! jump between the bracketing level centers to preserve its MEAN energy.
    ! Self-jumps do not enter the generator. This is not DL01's finite-width
    ! energy-bin integral or a resolution-independent variance of the jumps.
    ! Transitions beyond the grid are returned explicitly for EVERY initial
    ! level. The caller must resolve their probability-weighted contribution
    ! (extend the grid or model destruction), never silently omit these photons.
    real(real64),intent(in)::energy(:),photon_energy(:),photon_rate(:)
    real(real64),intent(out)::up(:,:),overflow_rate(:),overflow_power(:)
    integer,intent(out)::ierr
    real(real64)::matrix(size(energy),size(energy)),nr(size(energy)),power(size(energy)),w
    integer::n,i,g,k
    ierr=1;up=0;overflow_rate=0;overflow_power=0;n=size(energy)
    if(n<2.or.size(photon_energy)/=size(photon_rate))return
    if(any(shape(up)/=[n,n]).or.size(overflow_rate)/=n.or.size(overflow_power)/=n)return
    if(.not.all(ieee_is_finite(energy)).or.any(energy<0))return
    if(any(energy(2:)<=energy(:n-1)))return
    if(.not.all(ieee_is_finite(photon_energy)).or.any(photon_energy<=0))return
    if(.not.all(ieee_is_finite(photon_rate)).or.any(photon_rate<0))return
    matrix=0;nr=0;power=0
    do i=1,n
       do g=1,size(photon_energy)
          if(photon_rate(g)==0)cycle
          if(photon_energy(g)>energy(n)-energy(i))then
             nr(i)=nr(i)+photon_rate(g)
             power(i)=power(i)+photon_rate(g)*photon_energy(g)
             cycle
          endif
          k=i
          do while(k<n-1)
             if(energy(k+1)-energy(i)>photon_energy(g))exit
             k=k+1
          enddo
          ! Differences avoid losing a weak photon's energy when E_i is large.
          w=(photon_energy(g)-(energy(k)-energy(i)))/(energy(k+1)-energy(k))
          if(.not.ieee_is_finite(w).or.w<0.or.w>1)return
          if(k>i)matrix(k,i)=matrix(k,i)+photon_rate(g)*(1-w)
          matrix(k+1,i)=matrix(k+1,i)+photon_rate(g)*w
       enddo
    enddo
    if(.not.all(ieee_is_finite(matrix)).or..not.all(ieee_is_finite(nr)))return
    if(.not.all(ieee_is_finite(power)))return
    up=matrix;overflow_rate=nr;overflow_power=power;ierr=0
  end subroutine

  subroutine dust_stochastic_equilibrium(energy,up,down,probability,absorbed_power,emitted_power,ierr)
    ! Thermal-continuous cooling: arbitrary upward jumps and nearest-lower
    ! cooling. GD89 recursion as Draine & Li 2001 eqs. 53-54, evaluated in log
    ! space so strong heating does not overflow X_j=P_j/P_0. O(N^2).
    ! energy: erg/grain, increasing; up(destination,origin): s^-1;
    ! down(j): j -> j-1 rate in s^-1, with down(1)=0.
    ! Outputs are a normalized stationary probability and erg/s/grain.
    ! This is not the exact-statistical/discrete-emission DL01 treatment.
    ! No photon cutoff repair, spectral renormalization or relaxation-time
    ! assumption is performed here. Finite-grid rate coverage is the caller's
    ! responsibility; a solved truncated matrix is not proof of coverage.
    real(real64),intent(in)::energy(:),up(:,:),down(:)
    real(real64),intent(out)::probability(:),absorbed_power,emitted_power
    integer,intent(out)::ierr
    real(real64)::lb(size(energy),size(energy)),lp(size(energy)),p(size(energy))
    real(real64)::total,term,pa,pe
    logical::reachable(size(energy)),have,has_absorption,has_emission
    integer::n,i,j
    ierr=1;probability=0;absorbed_power=0;emitted_power=0;n=size(energy)
    if(n<2.or.size(down)/=n.or.size(probability)/=n)return
    if(any(shape(up)/=[n,n]))return
    if(.not.all(ieee_is_finite(energy)).or.any(energy<0))return
    if(any(energy(2:)<=energy(:n-1)))return
    if(.not.all(ieee_is_finite(up)).or.any(up<0))return
    if(.not.all(ieee_is_finite(down)))return
    if(down(1)/=0.or.any(down(2:)<=0))return
    do i=1,n
       if(any(up(:i,i)/=0))return
    enddo
    ! B(j,i) sums all heating jumps crossing the cut below level j.
    ! Negative HUGE denotes an absent (zero-rate) contribution only.
    lb=-huge(1d0)
    do i=1,n-1
       have=.false.;total=0
       do j=n,i+1,-1
          if(up(j,i)>0)then
             term=log(up(j,i))
             if(have)then
                total=log_add(total,term)
             else
                total=term;have=.true.
             endif
          endif
          if(have)lb(j,i)=total
       enddo
    enddo
    lp=0;reachable=.false.;reachable(1)=.true.
    do j=2,n
       have=.false.;total=0
       do i=1,j-1
          if(.not.reachable(i).or.lb(j,i)==-huge(1d0))cycle
          term=lp(i)+lb(j,i)
          if(have)then
             total=log_add(total,term)
          else
             total=term;have=.true.
          endif
       enddo
       if(have)then
          lp(j)=total-log(down(j));reachable(j)=.true.
       endif
    enddo
    total=maxval(lp,mask=reachable);p=0
    do j=1,n
       if(reachable(j))p(j)=exp(lp(j)-total)
    enddo
    term=log(sum(p))
    do j=1,n
       if(reachable(j))lp(j)=lp(j)-total-term
    enddo
    p=p/sum(p)
    ! Keep fluxes in log space too: a probability below the representable
    ! range times a large rate can still carry a representable energy flux.
    pa=0;pe=0;has_absorption=.false.;has_emission=.false.
    do i=1,n-1
       do j=i+1,n
          if(up(j,i)==0.or..not.reachable(i))cycle
          term=log(up(j,i))+lp(i)+log(energy(j)-energy(i))
          if(has_absorption)then
             pa=log_add(pa,term)
          else
             pa=term;has_absorption=.true.
          endif
       enddo
       if(.not.reachable(i+1))cycle
       term=log(down(i+1))+lp(i+1)+log(energy(i+1)-energy(i))
       if(has_emission)then
          pe=log_add(pe,term)
       else
          pe=term;has_emission=.true.
       endif
    enddo
    if(has_absorption)pa=exp(pa)
    if(has_emission)pe=exp(pe)
    if(.not.all(ieee_is_finite(p)).or..not.all(ieee_is_finite([pa,pe])))return
    if(any(p<0).or.min(pa,pe)<0)return
    probability=p;absorbed_power=pa;emitted_power=pe;ierr=0
  contains
    real(real64) function log_add(a,b)
      real(real64),intent(in)::a,b
      log_add=max(a,b)+log(1+exp(-abs(a-b)))
    end function
  end subroutine
end module

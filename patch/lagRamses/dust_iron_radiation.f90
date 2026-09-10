! Six-bin C/olivine/metallic-Fe material callback for the native IR receiver.
! Admitted Fe optical coefficients and persistent carriers are still required
! before selecting this callback in a simulation. No surrogate silicate Q_Fe.
module dust_iron_radiation
  use dust_iron_material
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: iso_c_binding, only: c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: iron_radiative_cell,iron_radiative_batch
  interface
     function gas_transfer_native(eg,cv,kdt,td) bind(C,name='snrt_dust_gas_transfer_c') result(q)
       import c_double
       real(c_double),value::eg,cv,kdt,td
       real(c_double)::q
     end function
  end interface
contains
  subroutine iron_radiative_batch(nodes,bath,bins,old_energy,heating,dt,reference_mass,basis_band, &
       workers,energy,temperature,phase,rate,ierr,gas_energy,gas_capacity,conductance,gas_transfer)
    ! Native CPU/OpenMP receiver, with the same six physical components as
    ! iron_radiative_cell. The four-basis CUDA material ABI cannot represent
    ! Fe phase plateaus. No GPU fallback/replay or implicit optical admission.
    ! All outputs are committed together ONLY after every cell has succeeded.
    ! Gas is an input reservoir: the enclosing IR transaction applies -Q.
    real(real64),intent(in)::nodes(:),bath,bins(:,:),old_energy(:),heating(:),dt,reference_mass
    real(real64),intent(in)::basis_band(:,:,:)
    integer,intent(in)::workers
    real(real64),intent(inout)::energy(:),temperature(:),phase(:,:),rate(:,:)
    integer,intent(out)::ierr
    real(real64),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
    real(real64),optional,intent(inout)::gas_transfer(:)
    real(real64),allocatable::work_e(:),work_t(:),work_phase(:,:),work_rate(:,:),work_q(:)
    real(real64)::eg,cv,k
    integer::nc,ng,i,status,bad
    logical::exchange
    ierr=1;nc=size(old_energy);ng=size(basis_band,1)
    if(workers<1.or.size(nodes)<2.or.ng<1)return
    if(any(shape(bins)/=[6,nc]).or.size(heating)/=nc)return
    if(size(energy)/=nc.or.size(temperature)/=nc.or.any(shape(phase)/=[4,nc]))return
    if(any(shape(rate)/=[ng,nc]).or.any(shape(basis_band)/=[ng,size(nodes),6]))return
    exchange=present(gas_energy)
    if(exchange.neqv.present(gas_capacity))return
    if(exchange.neqv.present(conductance))return
    if(exchange.neqv.present(gas_transfer))return
    if(exchange)then
       if(size(gas_energy)/=nc.or.size(gas_capacity)/=nc.or.size(conductance)/=nc)return
       if(size(gas_transfer)/=nc)return
    endif
    if(.not.all(ieee_is_finite([bath,dt,reference_mass])).or.min(dt,reference_mass)<=0)return
    if(any(.not.ieee_is_finite(nodes)).or.any(nodes<=0))return
    if(any(nodes(2:)<=nodes(:size(nodes)-1)).or.bath<nodes(1).or.bath>nodes(size(nodes)))return
    if(any(.not.ieee_is_finite(basis_band)).or.any(basis_band<0))return
    allocate(work_e(nc),work_t(nc),work_phase(4,nc),work_rate(ng,nc),work_q(nc))
    bad=0
!$omp parallel do num_threads(workers) private(i,status,eg,cv,k) reduction(max:bad)
    do i=1,nc
       eg=0;cv=1;k=0
       if(exchange)then
          eg=gas_energy(i);cv=gas_capacity(i);k=conductance(i)
       endif
       call iron_radiative_cell(nodes,bath,bins(:,i),old_energy(i),heating(i),dt,reference_mass,basis_band, &
            eg,cv,k,work_e(i),work_t(i),work_phase(:,i),work_rate(:,i),work_q(i),status)
       if(status/=0.and.bad==0)then
!$omp critical(fe_material_rejection)
          write(*,'(A,2I6,*(ES24.16,1X))')'Fe rejected cell/status, bath,bins,Ed,H,dt,ref,Eg,Cv,K: ', &
               i,status,bath,bins(:,i),old_energy(i),heating(i),dt,reference_mass,eg,cv,k
!$omp end critical(fe_material_rejection)
       endif
       bad=max(bad,status)
    enddo
!$omp end parallel do
    ierr=bad
    if(ierr/=0)return
    energy=work_e;temperature=work_t;phase=work_phase;rate=work_rate
    if(exchange)gas_transfer=work_q
  end subroutine

  subroutine iron_radiative_cell(nodes,bath,bins,old_energy,heating,dt,reference_mass,basis_band, &
       eg,cv,conductance,energy,temperature,phase,rate,transfer,ierr,absolute_emission,photon_ev)
    ! Backward Euler, fixed masses/optical coefficients during this call.
    ! cgs: bins=g/cm3, Ed/Eg=erg/cm3, heating=erg/cm3/s; band power/reference H.
    ! phase is [alpha,gamma,delta,liquid] fraction of metallic Fe. Ed contains
    ! the Fe transition enthalpy ALREADY: never add a second latent ledger.
    ! Collision transfer is the same native helper as the C/silicate receiver.
    ! Vacuum evaporation / mass exchange are not supplied by this operator.
    real(real64),intent(in)::nodes(:),bath,bins(6),old_energy,heating,dt,reference_mass
    real(real64),intent(in)::basis_band(:,:,:),eg,cv,conductance
    real(real64),intent(out)::energy,temperature,phase(4),rate(:),transfer
    integer,intent(out)::ierr
    real(real64)::bands(size(rate),size(nodes)),power(size(nodes)),mass(3),work_rate(size(rate))
    logical,optional,intent(in)::absolute_emission
    ! Supplying physical quadrature energies admits Planck emission below
    ! the first knot for the ABSOLUTE field only. Net-bath callers unchanged.
    real(real64),optional,intent(in)::photon_ev(:)
    logical::absolute
    real(real64)::bath_power,lower,upper,mid,lo_u,hi_u,t,q,target,tol,residual,trial_energy
    real(real64)::trans_t(3),trans_power,fraction,start,finish,width,p(4)
    integer::n,k,j,it,status
    energy=old_energy;temperature=bath;phase=0;rate=0;transfer=0;ierr=1
    absolute=.false.;if(present(absolute_emission))absolute=absolute_emission
    n=size(nodes)
    if(n<2.or.size(rate)<1.or.any(shape(basis_band)/=[size(rate),n,6]))return
    if(any(.not.ieee_is_finite(nodes)).or.any(nodes<=0))return
    if(any(nodes(2:)<=nodes(:n-1)))return
    if(.not.all(ieee_is_finite([bath,bins,old_energy,heating,dt,reference_mass,eg,cv,conductance])))return
    if(minval(bins)<0.or.min(old_energy,heating,eg,conductance)<0.or.min(dt,reference_mass,cv)<=0)return
    if(bath<nodes(1).or.bath>nodes(n))return
    if(present(photon_ev))then
       if(size(photon_ev)/=size(rate))return
       if(any(.not.ieee_is_finite(photon_ev)).or.any(photon_ev<=0))return
       ! Cold extension has no bath: the scalar argument is the first knot.
       ! Reject ambiguous alternate floors rather than leave a 5--bath gap.
       if(absolute.and.bath/=nodes(1))return
    endif
    if(.not.any(bins>0))then
       if(old_energy/=0.or.heating/=0.or.conductance/=0)return
       ierr=0;return
    endif
    mass=[sum(bins(1:2)),sum(bins(3:4)),sum(bins(5:6))]
    call iron_mixture_state(mass,old_energy,t,p,status)
    if(status/=0)return
    call iron_mixture_enthalpy(bath,mass,lo_u,hi_u,status)
    if(status/=0)return
    call iron_mixture_enthalpy(nodes(n),mass,lo_u,hi_u,status)
    if(status/=0)return
    if(any(.not.ieee_is_finite(basis_band)).or.any(basis_band<0))return
    bands=0
    do j=1,6
       bands=bands+basis_band(:,:,j)*(bins(j)/reference_mass)
    enddo
    if(any(.not.ieee_is_finite(bands)).or.any(bands(:,2:)<bands(:,:n-1)))return
    power=sum(bands,dim=1)
    if(any(.not.ieee_is_finite(power)).or.any(power(2:)<=power(:n-1)))return
    bath_power=power_at(bath)
    lower=0;upper=power(n)-bath_power
    if(absolute)then
       bath_power=0;lower=power_at(bath);upper=power(n)
       if(present(photon_ev))then
          ! There is no untracked bath in an absolute-radiation calculation.
          ! Below the first knot solve U(T)+dt*P(T)-Q(T)=Uold+dt*H directly.
          ! P(T) uses each band's Planck ratio, not grey T^4 extrapolation.
          call evaluate(power(1),nodes(1))
          if(status/=0)return
          if(lo_u>=target)then
             call cold_absolute()
             return
          endif
       endif
    endif
    ! Check each possible latent plateau explicitly. A smooth-T bisection
    ! cannot converge when energy is inside a discontinuous enthalpy jump.
    trans_t=[1184d0,1665d0,1809d0]
    do j=1,3
       t=trans_t(j)
       if(t<bath.or.t>nodes(n).or.mass(3)==0)cycle
       trans_power=power_at(t)-bath_power
       call evaluate(trans_power,t)
       if(status/=0)return
       if(hi_u>lo_u.and.target>=lo_u.and.target<=hi_u)then
          trial_energy=target;mid=trans_power
          p=0;fraction=(target-lo_u)/(hi_u-lo_u);p(j)=1-fraction;p(j+1)=fraction
          goto 100
       endif
       if(target<lo_u)then
          upper=trans_power;exit
       endif
       lower=trans_power
    enddo
    ierr=2
    t=temperature_at(lower)
    call evaluate(lower,t)
    if(status/=0.or.hi_u-target>tol)return
    ierr=3
    t=temperature_at(upper)
    call evaluate(upper,t)
    if(status/=0.or.lo_u-target < -tol)return
    ierr=4
    ! Cold material near the radiation bath can have net power tens of
    ! decades below the 1184/3000 K bracket. Retain the energy tolerance;
    ! allow enough bisections for long physical timesteps and dilute grains.
    do it=1,256
       mid=lower+(upper-lower)/2;t=temperature_at(mid)
       call evaluate(mid,t)
       if(status/=0)return
       ! Away from a transition lo=hi. At an exact boundary use the
       ! appropriate energy endpoint; don't smooth the phase jump.
       trial_energy=max(lo_u,min(target,hi_u))
       residual=trial_energy-target
       ! Leave headroom for the subsequent spectral partition and its
       ! independently recomputed energy ledger. Stopping exactly at the
       ! final tolerance can reject a converged dilute cell by one rounding
       ! increment; the final acceptance threshold itself is unchanged.
       if(abs(residual)<=.25d0*tol)exit
       if(residual>0)then
          upper=mid
       else
          lower=mid
       endif
    enddo
    if(abs(residual)>tol)then
       ierr=14;return
    endif
    call iron_mixture_state(mass,trial_energy,t,p,status)
    if(status/=0)return
100 continue
    if(eg-q<0.or..not.all(ieee_is_finite([trial_energy,t,q])))return
    work_rate=0
    if(absolute)work_rate=bands(:,1)
    do k=1,n-1
       start=max(power(k)-bath_power,0d0);finish=max(power(k+1)-bath_power,0d0)
       if(mid<=start)exit
       width=max(min(mid,finish)-start,0d0)
       work_rate=work_rate+width*(bands(:,k+1)-bands(:,k))/(power(k+1)-power(k))
    enddo
    if(any(.not.ieee_is_finite(work_rate)).or.any(work_rate<0))return
    residual=trial_energy-old_energy+dt*(sum(work_rate)-heating)-q
    if(.not.ieee_is_finite(residual).or.abs(residual)>tol)then
       ierr=16;return
    endif
    energy=trial_energy;temperature=t;phase=p;rate=work_rate;transfer=q;ierr=0
  contains
    subroutine cold_absolute()
      real(real64)::left_t,right_t,trial_t,x0,x,ratio,one0,one
      real(real64)::cold_rate(size(rate)),cold_power,err
      integer::g,iteration
      ierr=17;left_t=0;right_t=nodes(1)
      if(old_energy==0.and.heating==0.and.(conductance==0.or.eg==0))then
         energy=0;temperature=0;phase=0;rate=0;transfer=0;ierr=0
         if(mass(3)>0)phase(1)=1
         return
      endif
      do iteration=1,256
         trial_t=left_t+(right_t-left_t)/2
         cold_rate=0
         if(trial_t>0)then
            do g=1,size(rate)
               if(bands(g,1)==0)cycle
               x0=photon_ev(g)/(8.617333262145d-5*nodes(1))
               if(trial_t<photon_ev(g)/(8.617333262145d-5*745d0))cycle
               x=photon_ev(g)/(8.617333262145d-5*trial_t)
               one0=one_minus_exp(x0);one=one_minus_exp(x)
               ratio=exp(x0-x)*one0/one
               cold_rate(g)=bands(g,1)*ratio
            enddo
         endif
         cold_power=sum(cold_rate)
         call evaluate(cold_power,trial_t)
         if(status/=0)return
         err=lo_u-target
         if(abs(err)<=.25d0*tol)exit
         if(err>0)then
            right_t=trial_t
         else
            left_t=trial_t
         endif
      enddo
      if(abs(err)>tol.or.eg-q<0.or.any(.not.ieee_is_finite(cold_rate)))return
      if(any(cold_rate<0))return
      ! Publish only a conservative converged state; no clipped temperature.
      energy=lo_u;temperature=trial_t;phase=0
      if(mass(3)>0)phase(1)=1
      rate=cold_rate;transfer=q;ierr=0
    end subroutine
    real(real64) function one_minus_exp(x) result(y)
      real(real64),intent(in)::x
      if(x<1d-3)then
         y=x*(1-x*(.5d0-x*(1d0/6-x*(1d0/24-x/120))))
      else
         y=1-exp(-x)
      endif
    end function
    real(real64) function power_at(td) result(value)
      real(real64),intent(in)::td
      integer::i
      real(real64)::f
      i=1
      do while(i<n-1)
         if(td<=nodes(i+1))exit
         i=i+1
      enddo
      f=log(td/nodes(i))/log(nodes(i+1)/nodes(i))
      value=(1-f)*power(i)+f*power(i+1)
    end function
    real(real64) function temperature_at(net_power) result(td)
      real(real64),intent(in)::net_power
      integer::i
      real(real64)::f
      i=1
      do while(i<n-1)
         if(bath_power+net_power<=power(i+1))exit
         i=i+1
      enddo
      f=(bath_power+net_power-power(i))/(power(i+1)-power(i))
      td=max(bath,min(nodes(n),exp(log(nodes(i))+f*log(nodes(i+1)/nodes(i)))))
    end function
    subroutine evaluate(net_power,td)
      real(real64),intent(in)::net_power,td
      call iron_mixture_enthalpy(td,mass,lo_u,hi_u,status)
      q=0;target=0;tol=0
      if(status/=0)return
      if(conductance>0)q=gas_transfer_native(eg,cv,dt*conductance,td)
      target=old_energy+dt*(heating-net_power)+q
      tol=2d-12*max(old_energy,abs(target),lo_u,hi_u,dt*heating,dt*net_power,abs(q),tiny(1d0))
      if(.not.all(ieee_is_finite([target,q,tol])))status=1
    end subroutine
  end subroutine
end module

! Native secondary dust-energy operator. Equal-width reciprocal cell sets,
! vacuum exterior; no primary photons, species budgets, live AMR or gas writes.
module snrt_dust_ir
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: dust_dp=real64
  integer, parameter, public :: dust_ok=0, dust_err_table=1, dust_err_state=2
  integer, parameter, public :: dust_err_cfl=3, dust_err_convergence=4, dust_err_range=5
  integer, parameter, public :: dust_err_shape=6, dust_err_config=7
  real(real64), parameter :: ev_erg=1.602176634d-12, light_c=2.99792458d10
  type, public :: dust_ir_table
     private
     logical :: ready=.false.
     real(real64) :: background=0
     real(real64) :: background_temperature=0
     real(real64), allocatable :: energy(:), sigma(:), log_t(:), power(:), band(:,:)
  end type
  type, public :: dust_ir_diagnostics
     real(real64) :: escaped_erg=0, absorbed_erg=0, primary_erg=0
     ! Signed outward MPI boundary flux; unlike physical escape, it cancels
     ! when neighboring ranks' ledgers are summed.
     real(real64) :: interface_erg=0
     real(real64) :: balance_relative=0, local_relative=0
     integer :: iterations=0
  end type
  public :: snrt_dust_ir_initialize, snrt_dust_ir_advance
contains
  subroutine snrt_dust_ir_initialize(table, energy, frequency_weight, sigma, temperature, cmb, ierr)
    type(dust_ir_table), intent(out) :: table
    real(real64), intent(in) :: energy(:), frequency_weight(:), sigma(:), temperature(:), cmb
    integer, intent(out) :: ierr
    integer :: ng, nt, g, t, bath
    real(real64) :: x, occupation, factor
    real(real64), parameter :: h=6.62607015d-27, kb_ev=8.617333262145d-5
    ierr=dust_err_table
    ng=size(energy); nt=size(temperature)
    if (ng<1 .or. nt<2 .or. size(frequency_weight)/=ng .or. size(sigma)/=ng) return
    if (.not.all(ieee_is_finite(energy)) .or. .not.all(ieee_is_finite(frequency_weight))) return
    if (.not.all(ieee_is_finite(sigma)) .or. .not.all(ieee_is_finite(temperature))) return
    if (.not.ieee_is_finite(cmb)) return
    if (any(energy<=0) .or. any(frequency_weight<=0) .or. any(sigma<=0) .or. any(temperature<=0)) return
    if (any(energy(2:)<=energy(:ng-1)) .or. any(temperature(2:)<=temperature(:nt-1))) return
    bath=0
    do t=1,nt
       if (temperature(t)==cmb) bath=t
    end do
    if (bath==0) return
    allocate(table%energy(ng),table%sigma(ng),table%log_t(nt),table%power(nt),table%band(ng,nt))
    table%energy=energy; table%sigma=sigma; table%log_t=log(temperature)
    do t=1,nt
       do g=1,ng
          x=energy(g)/(kb_ev*temperature(t))
          if (x<1d-3) then
             occupation=1/x-.5d0+x/12-x**3/720
          else
             occupation=exp(-x)/(1-exp(-x))
          end if
          factor=8*acos(-1d0)*(energy(g)*ev_erg)**3/(h**3*light_c**2)*ev_erg
          table%band(g,t)=factor*sigma(g)*frequency_weight(g)*occupation
       end do
       table%power(t)=sum(table%band(:,t))
    end do
    if (.not.all(ieee_is_finite(table%band)) .or. .not.all(ieee_is_finite(table%power))) return
    if (any(table%band<0) .or. any(table%power<=0)) return
    if (any(table%power(2:)<=table%power(:nt-1))) return
    if (any(table%band(:,2:)<table%band(:,:nt-1))) return
    table%background=table%power(bath)
    table%background_temperature=cmb
    table%ready=.true.
    ierr=dust_ok
  end subroutine

  subroutine emission(table, heating, density, rate, temperature, ierr)
    type(dust_ir_table), intent(in) :: table
    real(real64), intent(in) :: heating(:), density(:)
    real(real64), intent(out) :: rate(:,:), temperature(:)
    integer, intent(out) :: ierr
    integer :: i, t, nt
    real(real64) :: increment, start, finish, width, target, fraction
    ierr=dust_err_state
    if (.not.all(ieee_is_finite(heating)) .or. any(heating<0)) return
    rate=0; temperature=0; nt=size(table%power)
    do i=1,size(density)
       if (heating(i)==0) cycle
       if (density(i)<=0) return
       increment=heating(i)/density(i)
       if (.not.ieee_is_finite(increment) .or. increment>table%power(nt)-table%background) then
          ierr=dust_err_range
          return
       end if
       do t=1,nt-1
          start=max(table%power(t)-table%background,0d0)
          finish=max(table%power(t+1)-table%background,0d0)
          width=max(min(increment,finish)-start,0d0)
          rate(:,i)=rate(:,i)+width*(table%band(:,t+1)-table%band(:,t)) &
               /(table%power(t+1)-table%power(t))*density(i)
       end do
       target=table%background+increment
       do t=1,nt-1
          if (target<=table%power(t+1)) exit
       end do
       t=min(t,nt-1)
       fraction=(target-table%power(t))/(table%power(t+1)-table%power(t))
       temperature(i)=exp(table%log_t(t)+fraction*(table%log_t(t+1)-table%log_t(t)))
    end do
    ierr=dust_ok
  end subroutine

  subroutine transient_emission(table,heating,density,dt,old_energy,capacity, &
       rate,temperature,next_energy,ierr)
    type(dust_ir_table), intent(in) :: table
    real(real64), intent(in) :: heating(:),density(:),dt,old_energy(:),capacity(:)
    real(real64), intent(out) :: rate(:,:),temperature(:),next_energy(:)
    integer, intent(out) :: ierr
    real(real64) :: lower,upper,mid,target,power,fraction,residual
    real(real64) :: emitted(size(heating)),unused_temperature(size(heating))
    integer :: i,k,iteration,nt
    ierr=dust_err_range
    nt=size(table%power)
    do i=1,size(heating)
       if(density(i)==0)then
          if(heating(i)/=0) return
          next_energy(i)=old_energy(i)
          temperature(i)=old_energy(i)/capacity(i)
          emitted(i)=0
          cycle
       end if
       target=old_energy(i)+dt*heating(i)
       lower=table%background_temperature
       upper=exp(table%log_t(nt))
       if(.not.ieee_is_finite(target))return
       ! Unit conversion of a stored C*T can move a bath-temperature state a
       ! few ulps below the bath. Admit rounding only; material/radiation
       ! closure below still accounts for any resulting energy correction.
       if(target<capacity(i)*lower*(1-64*epsilon(1d0))) return
       if(target>(capacity(i)*upper+dt*density(i)*(table%power(nt)-table%background)) &
            *(1+64*epsilon(1d0)))return
       ! Solve for emitted power, not the tiny temperature displacement of a
       ! stiff grain or the difference of two large material energies.
       lower=0d0
       upper=min(density(i)*(table%power(nt)-table%background), &
            max(heating(i)+(old_energy(i)-capacity(i)*table%background_temperature)/dt,0d0))
       do iteration=1,80
          mid=lower+0.5d0*(upper-lower)
          power=table%background+mid/density(i)
          k=1
          do while(k<nt-1)
             if(power<=table%power(k+1))exit
             k=k+1
          end do
          fraction=(power-table%power(k))/(table%power(k+1)-table%power(k))
          temperature(i)=exp(table%log_t(k)+fraction*(table%log_t(k+1)-table%log_t(k)))
          residual=(capacity(i)*temperature(i)-old_energy(i))/dt+mid-heating(i)
          if(residual>0)then
             upper=mid
          else
             lower=mid
          end if
       end do
       emitted(i)=lower+0.5d0*(upper-lower)
       power=table%background+emitted(i)/density(i)
       k=1
       do while(k<nt-1)
          if(power<=table%power(k+1))exit
          k=k+1
       end do
       fraction=(power-table%power(k))/(table%power(k+1)-table%power(k))
       temperature(i)=exp(table%log_t(k)+fraction*(table%log_t(k+1)-table%log_t(k)))
       next_energy(i)=capacity(i)*temperature(i)
    end do
    call emission(table,emitted,density,rate,unused_temperature,ierr)
  end subroutine

  subroutine snrt_dust_ir_advance(table, direction, weight, neighbor, dx, dt, c_hat, density, primary, &
       energy, temperature, photons, diagnostics, ierr, tolerance, max_iterations, dust_energy, heat_capacity, &
       ghost_energy,ghost_index,blocked_face)
    ! energy(g,d,cell): erg/cm3 per normalized direction; density: nH*relative_dust;
    ! primary: erg/cm3/s. photons(g,cell) accumulates emitted photons/cm3.
    ! Only success commits energy/temperature/photons/diagnostics. All trials
    ! begin at the SAME old field. No persistent hydro state is accessed.
    type(dust_ir_table), intent(in) :: table
    real(real64), intent(in) :: direction(:,:), weight(:), dx, dt, c_hat, density(:), primary(:)
    integer, intent(in) :: neighbor(:,:)
    real(real64), intent(inout) :: energy(:,:,:), temperature(:), photons(:,:)
    type(dust_ir_diagnostics), intent(inout) :: diagnostics
    integer, intent(out) :: ierr
    real(real64), intent(in) :: tolerance
    integer, intent(in) :: max_iterations
    ! Optional finite-capacity material state, both in physical volume units:
    ! energy erg/cm3 and capacity erg/cm3/K. Both or neither must be supplied.
    real(real64), optional, intent(inout) :: dust_energy(:)
    real(real64), optional, intent(in) :: heat_capacity(:)
    real(real64), optional, intent(in) :: ghost_energy(:,:,:)
    integer, optional, intent(in) :: ghost_index(:,:)
    ! A coarse face adjoining finer cells is advanced by the fine owner.
    ! Suppress BOTH inflow and outflow here; this is not a vacuum boundary.
    logical, optional, intent(in) :: blocked_face(:,:)
    logical, allocatable :: blocked(:,:)
    real(real64), allocatable :: transported(:,:,:), candidate(:,:,:), rate(:,:), next_t(:)
    real(real64), allocatable :: guess(:), absorbed(:), transmit(:,:), loss(:,:), response(:,:)
    real(real64), allocatable :: emitted_photons(:,:)
    real(real64), allocatable :: trial_dust_energy(:)
    logical :: transient
    real(real64) :: cfl, volume, factor, tau, source, old_total, new_total, scale, balance, sum_w
    integer :: ng, nd, nc, i, j, g, d, axis, face, outgoing, opposite, iteration, ghost
    integer, allocatable :: remote(:,:)
    type(dust_ir_diagnostics) :: trial
    trial=dust_ir_diagnostics()
    ierr=dust_err_table
    if (.not.table%ready) return
    transient=present(dust_energy)
    ierr=dust_err_config
    if(transient.neqv.present(heat_capacity))return
    if(present(ghost_energy).neqv.present(ghost_index))return
    ng=size(table%energy); nd=size(weight); nc=size(density)
    ierr=dust_err_shape
    if (nd<1 .or. nc<1) return
    if (any(shape(direction)/=[3,nd]) .or. any(shape(neighbor)/=[6,nc])) return
    if (any(shape(energy)/=[ng,nd,nc]) .or. any(shape(photons)/=[ng,nc])) return
    if (size(temperature)/=nc .or. size(primary)/=nc) return
    allocate(remote(6,nc)); remote=0
    allocate(blocked(6,nc)); blocked=.false.
    if(present(ghost_index))then
       if(any(shape(ghost_index)/=[6,nc]))return
       if(size(ghost_energy,1)/=ng.or.size(ghost_energy,2)/=nd)return
       if(any(ghost_index<0).or.any(ghost_index>size(ghost_energy,3)))return
       if(any(ghost_index>0.and.neighbor/=0))return
       ierr=dust_err_state
       if(any(.not.ieee_is_finite(ghost_energy)).or.any(ghost_energy<0))return
       remote=ghost_index
    end if
    ierr=dust_err_shape
    if(present(blocked_face))then
       if(any(shape(blocked_face)/=[6,nc]))return
       blocked=blocked_face
       if(any(blocked.and.(neighbor/=0.or.remote/=0)))return
    end if
    if(transient)then
       if(size(dust_energy)/=nc.or.size(heat_capacity)/=nc)return
       ierr=dust_err_state
       if(any(.not.ieee_is_finite(dust_energy)).or.any(.not.ieee_is_finite(heat_capacity)))return
       if(any(dust_energy<0).or.any(heat_capacity<=0))return
       do i=1,nc
          if(density(i)<=0)cycle
          if(dust_energy(i)/heat_capacity(i)<table%background_temperature*(1-64*epsilon(1d0)).or. &
               dust_energy(i)/heat_capacity(i)>exp(table%log_t(size(table%log_t)))*(1+64*epsilon(1d0)))return
       end do
    end if
    ierr=dust_err_config
    if (.not.all(ieee_is_finite([dx,dt,c_hat,tolerance]))) return
    if (min(dx,dt,c_hat,tolerance)<=0 .or. tolerance>=1 .or. max_iterations<1) return
    if (.not.all(ieee_is_finite(direction)) .or. .not.all(ieee_is_finite(weight))) return
    sum_w=sum(weight)
    if (any(weight<=0) .or. abs(sum_w-1)>1d-12) return
    if (any(abs(sqrt(sum(direction**2,dim=1))-1)>1d-6)) return
    cfl=c_hat*dt/dx*maxval(sum(abs(direction),dim=1))
    ierr=dust_err_cfl
    if (.not.ieee_is_finite(cfl) .or. cfl>1+1d-12) return
    ierr=dust_err_state
    if (.not.all(ieee_is_finite(energy)) .or. .not.all(ieee_is_finite(photons))) return
    if (.not.all(ieee_is_finite(temperature)) .or. .not.all(ieee_is_finite(density))) return
    if (.not.all(ieee_is_finite(primary))) return
    if (any(energy<0) .or. any(photons<0) .or. any(temperature<0) .or. any(density<0) .or. any(primary<0)) return
    if (any(neighbor<0) .or. any(neighbor>nc)) return
    do i=1,nc
       do face=1,6
          j=neighbor(face,i)
          if (j==0) cycle
          opposite=face+1
          if (mod(face,2)==0) opposite=face-1
          if (neighbor(opposite,j)/=i) return
       end do
    end do
    allocate(transported(ng,nd,nc),candidate(ng,nd,nc),rate(ng,nc),next_t(nc))
    allocate(guess(nc),absorbed(nc),transmit(ng,nc),loss(ng,nc),response(ng,nc),emitted_photons(ng,nc))
    if(transient)allocate(trial_dust_energy(nc))
    volume=dx**3
    transported=energy
    old_total=0
    do i=1,nc
       do d=1,nd
          old_total=old_total+sum(energy(:,d,i))*weight(d)*volume
          do axis=1,3
             face=2*axis-1; outgoing=2*axis
             if (direction(axis,d)<0) then
                face=2*axis; outgoing=2*axis-1
             end if
             j=neighbor(face,i)
             factor=c_hat*dt/dx*abs(direction(axis,d))
             if(.not.blocked(outgoing,i)) &
                  transported(:,d,i)=transported(:,d,i)-factor*energy(:,d,i)
             if (j>0) transported(:,d,i)=transported(:,d,i)+factor*energy(:,d,j)
             ghost=remote(face,i)
             if(ghost>0)then
                transported(:,d,i)=transported(:,d,i)+factor*ghost_energy(:,d,ghost)
                trial%interface_erg=trial%interface_erg-sum(ghost_energy(:,d,ghost))*weight(d)*factor*volume
             end if
             if (neighbor(outgoing,i)==0.and..not.blocked(outgoing,i)) then
                if(remote(outgoing,i)>0)then
                   trial%interface_erg=trial%interface_erg+sum(energy(:,d,i))*weight(d)*factor*volume
                else
                   trial%escaped_erg=trial%escaped_erg+sum(energy(:,d,i))*weight(d)*factor*volume
                end if
             end if
          end do
       end do
       do g=1,ng
          tau=c_hat*dt*table%sigma(g)*density(i)
          if (.not.ieee_is_finite(tau)) return
          transmit(g,i)=exp(-tau)
          if (tau<1d-4) then
             loss(g,i)=tau*(1-tau/2+tau*tau/6-tau**3/24)
             response(g,i)=1-tau/2+tau*tau/6
          else
             loss(g,i)=1-transmit(g,i)
             response(g,i)=loss(g,i)/max(tau,tiny(tau))
          end if
       end do
    end do
    trial%primary_erg=sum(primary)*dt*volume
    if (.not.all(ieee_is_finite([volume,old_total,trial%primary_erg,trial%escaped_erg,trial%interface_erg]))) return
    scale=max(trial%primary_erg,tiny(scale))
    if(trial%primary_erg==0)scale=max(old_total,abs(trial%interface_erg),scale)
    if(transient)then
       scale=max(scale,sum(dust_energy)*volume)
       if(.not.ieee_is_finite(scale))return
    end if
    guess=0
    do iteration=0,max_iterations
       if(transient)then
          call transient_emission(table,primary+guess/dt,density,dt,dust_energy,heat_capacity, &
               rate,next_t,trial_dust_energy,ierr)
       else
          call emission(table, primary+guess/dt, density, rate, next_t, ierr)
       end if
       if (ierr/=dust_ok) return
       absorbed=0; new_total=0
       do i=1,nc
          do d=1,nd
             do g=1,ng
                source=dt*rate(g,i)/sum_w
                candidate(g,d,i)=transported(g,d,i)*transmit(g,i)+source*response(g,i)
                absorbed(i)=absorbed(i)+weight(d)*(transported(g,d,i)*loss(g,i)+source*(1-response(g,i)))
             end do
             new_total=new_total+sum(candidate(:,d,i))*weight(d)*volume
          end do
       end do
       ierr=dust_err_state
       if (.not.all(ieee_is_finite(candidate)) .or. any(candidate<0) .or. .not.ieee_is_finite(new_total)) return
       if (.not.all(ieee_is_finite(absorbed)) .or. any(absorbed<0)) return
       balance=new_total-old_total+trial%escaped_erg+trial%interface_erg-trial%primary_erg
       if(transient)balance=balance+sum(trial_dust_energy-dust_energy)*volume
       trial%balance_relative=abs(balance)/scale
       trial%local_relative=maxval(abs(absorbed-guess)/max(primary*dt+absorbed,tiny(scale)))
       if (max(trial%balance_relative,trial%local_relative)<=tolerance) then
          do g=1,ng
             emitted_photons(g,:)=rate(g,:)*dt/(table%energy(g)*ev_erg)
          end do
          if (.not.all(ieee_is_finite(photons+emitted_photons))) return
          trial%iterations=iteration; trial%absorbed_erg=sum(absorbed)*volume
          if (.not.ieee_is_finite(trial%absorbed_erg)) return
          energy=candidate; temperature=next_t; photons=photons+emitted_photons
          if(transient)dust_energy=trial_dust_energy
          diagnostics=trial
          ierr=dust_ok
          return
       end if
       guess=.5d0*(guess+absorbed)
    end do
    ierr=dust_err_convergence
  end subroutine
end module

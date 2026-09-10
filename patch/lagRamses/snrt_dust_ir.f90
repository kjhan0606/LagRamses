! Native secondary dust-energy operator. Equal-width reciprocal cell sets,
! vacuum exterior; no primary photons, species budgets, live AMR or gas writes.
module snrt_dust_ir
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use snrt_moving_scatter, only: snrt_moving_scatter_energy_cell,snrt_scatter_beta_limit
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
     ! Vibrational internal energy per reference H, excluding zero-point
     ! energy. Interpolate U linearly in log(T), with no extrapolation.
     real(real64), allocatable :: material_u(:)
     real(real64), allocatable :: optical_sigma(:,:),optical_power(:,:),optical_band(:,:,:)
  end type
  type, public :: dust_ir_diagnostics
     real(real64) :: escaped_erg=0, absorbed_erg=0, primary_erg=0
     real(real64) :: mechanical_erg=0 ! signed grain kinetic gain, never heat
     ! Signed outward MPI boundary flux; unlike physical escape, it cancels
     ! when neighboring ranks' ledgers are summed.
     real(real64) :: interface_erg=0
     real(real64) :: balance_relative=0, local_relative=0
     integer :: iterations=0
  end type
  public :: snrt_dust_ir_initialize, snrt_dust_ir_advance, snrt_dust_material_temperature
  public :: dust_moving_population_dispatch
  abstract interface
     subroutine dust_moving_population_dispatch(ir_heat,ir_captured,dt,old_population,next_population, &
          phase_rate,next_energy,temperature,ierr)
       ! heat/captures/phase_rate: (phase,group,cell). Heat is COMOVING
       ! deposited erg/cm3; captures are nominal LAB group events/cm3 from
       ! absorbed_E/(group_eV*ev_erg), not heat/(group_eV*ev_erg).
       ! phase_rate is comoving emission erg/cm3/s, separately per phase.
       ! next_energy contains exactly the same material reservoirs as the
       ! enclosing dust_energy (e.g. bulk+PAH+gas in pah_mixed_advance).
       ! Every callback trial MUST start from immutable old material inputs.
       import real64
       real(real64),intent(in)::ir_heat(:,:,:),ir_captured(:,:,:),dt,old_population(:,:)
       real(real64),intent(out)::next_population(:,:),phase_rate(:,:,:),next_energy(:),temperature(:)
       integer,intent(out)::ierr
     end subroutine
     subroutine dust_population_dispatch(ir_absorbed,dt,old_population,next_population,rate,next_energy,temperature,ierr)
       ! ir_absorbed(g,cell) is the actual spectral energy fluence, erg/cm3.
       ! Primary spectral fluence belongs to the caller's immutable context;
       ! its sum MUST equal primary*dt in the enclosing energy ledger.
       ! Each nonlinear trial starts from old_population, never the last trial.
       import real64
       real(real64),intent(in)::ir_absorbed(:,:),dt,old_population(:,:)
       real(real64),intent(out)::next_population(:,:),rate(:,:),next_energy(:),temperature(:)
       integer,intent(out)::ierr
     end subroutine
     subroutine dust_transport_dispatch(energy,ghosts,neighbor,remote,blocked,density,direction,sigma, &
          cdt,ratio,transported,transmit,loss,response,ierr,cell_sigma)
       import real64
       real(real64),intent(in)::energy(:,:,:),ghosts(:,:,:),density(:),direction(:,:),sigma(:),cdt,ratio
       integer,intent(in)::neighbor(:,:),remote(:,:)
       logical,intent(in)::blocked(:,:)
       real(real64),intent(out)::transported(:,:,:),transmit(:,:),loss(:,:),response(:,:)
       integer,intent(out)::ierr
       real(real64),optional,intent(in)::cell_sigma(:,:)
     end subroutine
     subroutine dust_absorb_dispatch(transported,transmit,loss,response,rate,weight,dt,sum_w,candidate,absorbed,ierr)
       import real64
       real(real64),intent(in)::transported(:,:,:),transmit(:,:),loss(:,:),response(:,:),rate(:,:),weight(:),dt,sum_w
       real(real64),intent(out)::candidate(:,:,:),absorbed(:)
       integer,intent(out)::ierr
     end subroutine
     subroutine dust_material_dispatch(heating,density,old_energy,capacity,log_t,power,band, &
          material_u,use_u,dt,background,bath,tolerance,rate,temperature,next_energy,ierr, &
          gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,basis_power,basis_band)
       import real64
       real(real64),intent(in)::heating(:),density(:),old_energy(:),capacity(:),log_t(:),power(:),band(:,:)
       real(real64),intent(in)::material_u(:),dt,background,bath,tolerance
       logical,intent(in)::use_u
       real(real64),intent(out)::rate(:,:),temperature(:),next_energy(:)
       integer,intent(out)::ierr
       real(real64),optional,intent(in)::gas_energy(:),gas_capacity(:),conductance(:)
       real(real64),optional,intent(out)::gas_transfer(:)
       real(real64),optional,intent(in)::cell_material_u(:,:)
       real(real64),optional,intent(in)::cell_weights(:,:),basis_power(:,:),basis_band(:,:,:)
     end subroutine
  end interface
contains
  recursive subroutine snrt_dust_ir_initialize(table, energy, frequency_weight, sigma, temperature, cmb, &
       ierr, material_u, optical_sigma)
    type(dust_ir_table), intent(out) :: table
    real(real64), intent(in) :: energy(:), frequency_weight(:), sigma(:), temperature(:), cmb
    integer, intent(out) :: ierr
    real(real64), optional, intent(in) :: material_u(:)
    real(real64),optional,intent(in)::optical_sigma(:,:)
    type(dust_ir_table)::part
    integer :: ng, nt, g, t, bath,s,nbasis
    real(real64) :: x, occupation, factor
    real(real64), parameter :: h=6.62607015d-27, kb_ev=8.617333262145d-5
    ierr=dust_err_table
    ng=size(energy); nt=size(temperature)
    if (ng<1 .or. nt<2 .or. size(frequency_weight)/=ng .or. size(sigma)/=ng) return
    if (.not.all(ieee_is_finite(energy)) .or. .not.all(ieee_is_finite(frequency_weight))) return
    if (.not.all(ieee_is_finite(sigma)) .or. .not.all(ieee_is_finite(temperature))) return
    if (.not.ieee_is_finite(cmb)) return
    if (any(energy<=0) .or. any(frequency_weight<=0) .or. any(sigma<0) .or. any(temperature<=0)) return
    if (any(energy(2:)<=energy(:ng-1)) .or. any(temperature(2:)<=temperature(:nt-1))) return
    bath=0
    do t=1,nt
       if (temperature(t)==cmb) bath=t
    end do
    if (bath==0) return
    if(present(material_u))then
       if(size(material_u)/=nt)return
       if(any(.not.ieee_is_finite(material_u)).or.any(material_u<=0))return
       if(any(material_u(2:)<=material_u(:nt-1)))return
       table%material_u=material_u
    endif
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
    if(present(optical_sigma))then
       nbasis=size(optical_sigma,2)
       if(size(optical_sigma,1)/=ng.or.nbasis<1)return
       allocate(table%optical_sigma(ng,nbasis),table%optical_power(nt,nbasis),table%optical_band(ng,nt,nbasis))
       table%optical_sigma=optical_sigma
       do s=1,nbasis
          call snrt_dust_ir_initialize(part,energy,frequency_weight,optical_sigma(:,s),temperature,cmb,ierr,material_u)
          if(ierr/=dust_ok)return
          table%optical_power(:,s)=part%power;table%optical_band(:,:,s)=part%band
       enddo
    endif
    table%ready=.true.
    ierr=dust_ok
  end subroutine

  subroutine snrt_dust_material_temperature(nodes,u,energy,temperature,ierr)
    ! Inverse of the v4 material interpolation. Inputs are per reference H;
    ! caller divides volume energy by nH*relative_dust, not by heat capacity.
    real(real64), intent(in) :: nodes(:),u(:),energy
    real(real64), intent(out) :: temperature
    integer, intent(out) :: ierr
    integer :: k,n
    real(real64) :: fraction
    ierr=dust_err_table; temperature=0; n=size(nodes)
    if(n<2.or.size(u)/=n)return
    if(any(.not.ieee_is_finite(nodes)).or.any(.not.ieee_is_finite(u)))return
    if(any(nodes<=0).or.any(u<=0))return
    if(any(nodes(2:)<=nodes(:n-1)).or.any(u(2:)<=u(:n-1)))return
    ierr=dust_err_range
    if(.not.ieee_is_finite(energy))return
    if(energy<u(1)*(1-64*epsilon(1d0)).or.energy>u(n)*(1+64*epsilon(1d0)))return
    k=1
    do while(k<n-1)
       if(energy<=u(k+1))exit
       k=k+1
    enddo
    fraction=max(0d0,min(1d0,(energy-u(k))/(u(k+1)-u(k))))
    temperature=exp(log(nodes(k))+fraction*log(nodes(k+1)/nodes(k)))
    ierr=dust_ok
  end subroutine

  real(real64) function material_energy(table,temperature,density,capacity,cell_u) result(energy)
    type(dust_ir_table), intent(in) :: table
    real(real64), intent(in) :: temperature,density,capacity
    real(real64),optional,intent(in)::cell_u(:)
    integer :: k,n
    real(real64) :: log_t,fraction
    if(.not.allocated(table%material_u))then
       energy=capacity*temperature
       return
    endif
    n=size(table%log_t); log_t=log(temperature); k=1
    do while(k<n-1)
       if(log_t<=table%log_t(k+1))exit
       k=k+1
    enddo
    fraction=max(0d0,min(1d0,(log_t-table%log_t(k))/(table%log_t(k+1)-table%log_t(k))))
    energy=density*(table%material_u(k)+fraction*(table%material_u(k+1)-table%material_u(k)))
    if(present(cell_u))energy=density*(cell_u(k)+fraction*(cell_u(k+1)-cell_u(k)))
  end function

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
       rate,temperature,next_energy,ierr,material_tolerance)
    type(dust_ir_table), intent(in) :: table
    real(real64), intent(in) :: heating(:),density(:),dt,old_energy(:),capacity(:)
    real(real64), intent(out) :: rate(:,:),temperature(:),next_energy(:)
    integer, intent(out) :: ierr
    real(real64), intent(in) :: material_tolerance
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
          if(allocated(table%material_u))then
             if(old_energy(i)/=0)return
             temperature(i)=0
          endif
          emitted(i)=0
          cycle
       end if
       target=old_energy(i)+dt*heating(i)
       lower=table%background_temperature
       upper=exp(table%log_t(nt))
       if(.not.ieee_is_finite(target))return
       ! Hydro advects material mass and energy separately. A tiny deficit
       ! below the bath is admissible only within the declared solve error.
       ! Do NOT replace old_energy: the full floor correction is charged to
       ! material+radiation closure in advance(), not hidden as bath heating.
       if(target<material_energy(table,lower,density(i),capacity(i))*(1-material_tolerance)) return
       if(target>(material_energy(table,upper,density(i),capacity(i))+ &
            dt*density(i)*(table%power(nt)-table%background)) &
            *(1+64*epsilon(1d0)))return
       ! Solve for emitted power, not the tiny temperature displacement of a
       ! stiff grain or the difference of two large material energies.
       lower=0d0
       upper=min(density(i)*(table%power(nt)-table%background), &
            max(heating(i)+(old_energy(i)- &
            material_energy(table,table%background_temperature,density(i),capacity(i)))/dt,0d0))
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
          residual=(material_energy(table,temperature(i),density(i),capacity(i))-old_energy(i))/dt+mid-heating(i)
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
       next_energy(i)=material_energy(table,temperature(i),density(i),capacity(i))
    end do
    call emission(table,emitted,density,rate,unused_temperature,ierr)
  end subroutine

  subroutine snrt_dust_ir_advance(table, direction, weight, neighbor, dx, dt, c_hat, density, primary, &
       energy, temperature, photons, diagnostics, ierr, tolerance, max_iterations, dust_energy, heat_capacity, &
       ghost_energy,ghost_index,blocked_face,material_dispatch,transport_dispatch,absorb_dispatch, &
       gas_energy,gas_capacity,conductance,gas_transfer,cell_material_u,cell_weights,thin_reabsorption, &
       population,population_dispatch,cell_absorption,phase_density,phase_momentum,phase_absorption, &
       phase_scattering,phase_work,moving_material_dispatch)
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
    ! When table%material_u is present, capacity is only a positive ABI
    ! placeholder: physical material energy is density*U(T), never capacity*T.
    real(real64), optional, intent(inout) :: dust_energy(:)
    real(real64), optional, intent(in) :: heat_capacity(:)
    real(real64), optional, intent(in) :: ghost_energy(:,:,:)
    integer, optional, intent(in) :: ghost_index(:,:)
    ! A coarse face adjoining finer cells is advanced by the fine owner.
    ! Suppress BOTH inflow and outflow here; this is not a vacuum boundary.
    logical, optional, intent(in) :: blocked_face(:,:)
    procedure(dust_material_dispatch),optional :: material_dispatch
    procedure(dust_transport_dispatch),optional :: transport_dispatch
    procedure(dust_absorb_dispatch),optional :: absorb_dispatch
    real(real64),optional,intent(inout)::gas_energy(:)
    real(real64),optional,intent(in)::gas_capacity(:),conductance(:)
    real(real64),optional,intent(inout)::gas_transfer(:)
    real(real64),optional,intent(in)::cell_material_u(:,:)
    real(real64),optional,intent(in)::cell_weights(:,:)
    logical,optional,intent(in)::thin_reabsorption
    ! Stochastic material uses transported grain-number populations, not a
    ! common-temperature U(T). This branch transports ABSOLUTE IR radiation;
    ! any background illumination must be supplied as photons by the caller.
    ! Mutually exclusive with the equilibrium/gas/material callbacks.
    real(real64),optional,intent(inout)::population(:,:)
    procedure(dust_population_dispatch),optional::population_dispatch
    ! Explicit mixed-material opacity per density unit. The population is
    ! only ONE component, so density need not equal its grain number.
    real(real64),optional,intent(in)::cell_absorption(:,:)
    ! Moving interactions: all physical cgs, frozen phase masses/opacities.
    ! Radiation still stores energy per normalized direction; quadrature
    ! weights enter exactly once in impulses, material heat and captures.
    real(real64),optional,intent(in)::phase_density(:,:),phase_absorption(:,:,:),phase_scattering(:,:,:)
    real(real64),optional,intent(inout)::phase_momentum(:,:,:),phase_work(:,:)
    procedure(dust_moving_population_dispatch),optional::moving_material_dispatch
    real(real64),allocatable::phase_heat(:,:,:),heat_guess(:,:,:),phase_events(:,:,:),event_guess(:,:,:)
    real(real64),allocatable::phase_rate(:,:,:),beta_guess(:,:,:),beta_next(:,:,:),momentum_next(:,:,:)
    real(real64),allocatable::work_next(:,:),scatter_work(:),scatter_e(:,:),scatter_tau(:,:)
    real(real64)::moving_residual,phase_scale,scatter_balance
    integer::nb,b,moving_status
    logical::moving
    real(real64),allocatable::trial_population(:,:),spectral_guess(:,:),spectral_absorbed(:,:)
    logical,allocatable::direct_update(:)
    real(real64),allocatable::cell_sigma(:,:)
    real(real64),allocatable::exchange(:)
    real(real64),allocatable :: empty_ghost(:,:,:)
    real(real64),allocatable :: dispatch_u(:)
    logical, allocatable :: blocked(:,:)
    real(real64), allocatable :: transported(:,:,:), candidate(:,:,:), rate(:,:), next_t(:)
    real(real64), allocatable :: guess(:), absorbed(:), transmit(:,:), loss(:,:), response(:,:)
    real(real64), allocatable :: emitted_photons(:,:)
    real(real64), allocatable :: trial_dust_energy(:)
    logical :: transient
    real(real64) :: cfl, volume, factor, tau, source, old_total, new_total, scale, balance, sum_w
    real(real64) :: material_tolerance
    integer :: ng, nd, nc, i, j, g, d, axis, face, outgoing, opposite, iteration, ghost
    integer, allocatable :: remote(:,:)
    type(dust_ir_diagnostics) :: trial
    trial=dust_ir_diagnostics()
    ierr=dust_err_table
    if (.not.table%ready) return
    transient=present(dust_energy)
    moving=present(phase_density)
    ierr=dust_err_config
    if(moving.neqv.present(phase_momentum))return
    if(moving.neqv.present(phase_absorption))return
    if(moving.neqv.present(phase_work))return
    if(.not.moving.and.(present(phase_scattering).or.present(moving_material_dispatch)))return
    if(moving.and..not.transient)return
    if(transient.neqv.present(heat_capacity))return
    if(present(population).neqv.(present(population_dispatch).or.present(moving_material_dispatch)))return
    if(present(moving_material_dispatch).and.present(population_dispatch))return
    if(moving.and.present(population).and..not.present(moving_material_dispatch))return
    if(present(population))then
       if(.not.transient)return
       if(present(material_dispatch).or.present(gas_energy).or.present(cell_material_u).or. &
            present(cell_weights).or.present(absorb_dispatch))return
       if(size(population,1)<2.or.size(population,2)/=size(density))return
       if(any(.not.ieee_is_finite(population)).or.any(population<0))return
       ! Here density is grain number/cm3 and sigma is cross section/grain.
       if(.not.present(cell_absorption))then
          if(any(abs(sum(population,dim=1)-density)> &
               512*epsilon(1d0)*size(population,1)*max(density,tiny(1d0))))return
       endif
    endif
    if(present(ghost_energy).neqv.present(ghost_index))return
    ng=size(table%energy); nd=size(weight); nc=size(density)
    if(moving)then
       nb=size(phase_density,1)
       if(nb<1.or.size(phase_density,2)/=nc)return
       if(any(shape(phase_momentum)/=[3,nb,nc]).or.any(shape(phase_work)/=[nb,nc]))return
       if(any(shape(phase_absorption)/=[nb,ng,nc]))return
       if(any(.not.ieee_is_finite(phase_density)).or.any(phase_density<0))return
       if(any(.not.ieee_is_finite(phase_momentum)))return
       if(any(.not.ieee_is_finite(phase_absorption)).or.any(phase_absorption<0))return
       if(present(phase_scattering))then
          if(any(shape(phase_scattering)/=[nb,ng,nc]))return
          if(any(.not.ieee_is_finite(phase_scattering)).or.any(phase_scattering<0))return
       endif
       if(any(weight<=0))return
       allocate(phase_heat(nb,ng,nc),heat_guess(nb,ng,nc),phase_events(nb,ng,nc),event_guess(nb,ng,nc))
       allocate(phase_rate(nb,ng,nc),beta_guess(3,nb,nc),beta_next(3,nb,nc),momentum_next(3,nb,nc))
       allocate(work_next(nb,nc));beta_guess=0;heat_guess=0;event_guess=0;work_next=0
       do i=1,nc
          do b=1,nb
             if(phase_density(b,i)==0)then
                if(any(phase_momentum(:,b,i)/=0).or.any(phase_absorption(b,:,i)/=0))return
                if(present(phase_scattering))then
                   if(any(phase_scattering(b,:,i)/=0))return
                endif
             else
                beta_guess(:,b,i)=phase_momentum(:,b,i)/(phase_density(b,i)*light_c)
                if(norm2(beta_guess(:,b,i))>snrt_scatter_beta_limit)return
             endif
          enddo
       enddo
    endif
    if(present(cell_absorption))then
       if(present(cell_weights).or..not.present(population))return
       if(any(shape(cell_absorption)/=[ng,nc]))return
       if(any(.not.ieee_is_finite(cell_absorption)).or.any(cell_absorption<0))return
       cell_sigma=cell_absorption
    endif
    if(present(cell_weights).neqv.allocated(table%optical_sigma))return
    if(present(cell_weights))then
       if(.not.transient.or..not.present(material_dispatch).or..not.present(cell_material_u))return
       if(any(shape(cell_weights)/=[size(table%optical_sigma,2),nc]))return
       if(any(.not.ieee_is_finite(cell_weights)).or.any(cell_weights<0))return
       if(any(abs(sum(cell_weights,dim=1)-1)>1d-12))return
       cell_sigma=matmul(table%optical_sigma,cell_weights)
    endif
    if(present(cell_material_u))then
       if(.not.transient.or..not.present(material_dispatch).or..not.allocated(table%material_u))return
       if(any(shape(cell_material_u)/=[size(table%log_t),nc]))return
       if(any(.not.ieee_is_finite(cell_material_u)).or.any(cell_material_u<=0))return
       if(any(cell_material_u(2:,:)<=cell_material_u(:size(table%log_t)-1,:)))return
    endif
    if(present(gas_energy))then
       if(.not.transient.or..not.present(material_dispatch).or..not.allocated(table%material_u))return
       if(.not.present(gas_capacity).or..not.present(conductance).or..not.present(gas_transfer))return
       if(size(gas_energy)/=nc.or.size(gas_capacity)/=nc.or.size(conductance)/=nc.or.size(gas_transfer)/=nc)return
       if(any(.not.ieee_is_finite(gas_energy)).or.any(gas_energy<0))return
       allocate(exchange(nc));exchange=0
    else if(present(gas_capacity).or.present(conductance).or.present(gas_transfer))then
       return
    endif
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0.or.tolerance>=1)return
    material_tolerance=max(tolerance,64*epsilon(1d0))
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
          if(present(population))cycle
          if(density(i)<=0)cycle
          if(present(cell_material_u))then
             if(dust_energy(i)<material_energy(table,table%background_temperature,density(i),heat_capacity(i), &
                  cell_material_u(:,i))*(1-material_tolerance).or.dust_energy(i)> &
                  material_energy(table,exp(table%log_t(size(table%log_t))),density(i),heat_capacity(i), &
                  cell_material_u(:,i))*(1+64*epsilon(1d0)))return
             cycle
          endif
          if(dust_energy(i)<material_energy(table,table%background_temperature,density(i),heat_capacity(i)) &
               *(1-material_tolerance).or.dust_energy(i)> &
               material_energy(table,exp(table%log_t(size(table%log_t))),density(i),heat_capacity(i)) &
               *(1+64*epsilon(1d0)))return
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
    if(present(population))then
       allocate(trial_population(size(population,1),nc),spectral_guess(ng,nc),spectral_absorbed(ng,nc))
       spectral_guess=0
    endif
    if(transient.and.present(material_dispatch))then
       allocate(dispatch_u(size(table%log_t)));dispatch_u=0
       if(allocated(table%material_u))dispatch_u=table%material_u
    endif
    volume=dx**3
    transported=energy
    if(present(transport_dispatch))then
       if(present(ghost_energy))then
          call transport_dispatch(energy,ghost_energy,neighbor,remote,blocked,density,direction,table%sigma, &
               c_hat*dt,c_hat*dt/dx,transported,transmit,loss,response,ierr,cell_sigma)
       else
          allocate(empty_ghost(ng,nd,0))
          call transport_dispatch(energy,empty_ghost,neighbor,remote,blocked,density,direction,table%sigma, &
               c_hat*dt,c_hat*dt/dx,transported,transmit,loss,response,ierr,cell_sigma)
       endif
       if(ierr/=dust_ok)return
    endif
    ierr=dust_err_state
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
             if(.not.present(transport_dispatch))then
                if(.not.blocked(outgoing,i)) &
                     transported(:,d,i)=transported(:,d,i)-factor*energy(:,d,i)
                if (j>0) transported(:,d,i)=transported(:,d,i)+factor*energy(:,d,j)
             endif
             ghost=remote(face,i)
             if(ghost>0)then
                if(.not.present(transport_dispatch)) &
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
          if(present(transport_dispatch))exit
          tau=c_hat*dt*table%sigma(g)*density(i)
          if(allocated(cell_sigma))tau=c_hat*dt*cell_sigma(g,i)*density(i)
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
    if(present(population))then
       ! Known absorption of transported OLD radiation belongs in the very
       ! first material trial too, independent of newly emitted photons.
       do i=1,nc
          do d=1,nd
             spectral_guess(:,i)=spectral_guess(:,i)+weight(d)*transported(:,d,i)*loss(:,i)
          enddo
       enddo
    endif
    if (.not.all(ieee_is_finite([volume,old_total,trial%primary_erg,trial%escaped_erg,trial%interface_erg]))) return
    scale=max(trial%primary_erg,tiny(scale))
    if(trial%primary_erg==0)scale=max(old_total,abs(trial%interface_erg),scale)
    if(transient)then
       scale=max(scale,sum(dust_energy)*volume)
       if(.not.ieee_is_finite(scale))return
    end if
    guess=0
    allocate(direct_update(nc));direct_update=.false.
    if(present(thin_reabsorption))then
       ! Emission reabsorbed in this cell is bounded in every band by
       ! max(1-response). Avoid 0.5 damping of a weak feedback map, which
       ! otherwise repeats an expensive adaptive hot transient ~30 times.
       ! This only changes the iteration, never the final acceptance tests.
       if(thin_reabsorption)direct_update=maxval(1-response,dim=1)<.25d0
       if(thin_reabsorption.and.moving)then
          ! Same optical feedback bound, using the maximum directional
          ! depth throughout the admitted beta ball, not a static surrogate.
          do i=1,nc
             tau=c_hat*dt*maxval(sum(phase_absorption(:,:,i),dim=1))*(1+snrt_scatter_beta_limit)
             if(tau<1d-4)then
                direct_update(i)=tau*(.5d0-tau/6+tau*tau/24)<.25d0
             else
                direct_update(i)=1-(1-exp(-tau))/tau<.25d0
             endif
          enddo
       endif
    endif
    do iteration=0,max_iterations
       if(transient)then
          if(present(population_dispatch).or.present(moving_material_dispatch))then
             if(moving)then
                call moving_material_dispatch(heat_guess,event_guess,dt,population,trial_population, &
                     phase_rate,trial_dust_energy,next_t,ierr)
                if(ierr/=dust_ok)return
                rate=sum(phase_rate,dim=1)
                if(any(.not.ieee_is_finite(phase_rate)).or.any(phase_rate<0))ierr=dust_err_state
             else
                call population_dispatch(spectral_guess,dt,population,trial_population,rate,trial_dust_energy,next_t,ierr)
             endif
             if(ierr/=dust_ok)return
             ierr=dust_err_state
             if(any(.not.ieee_is_finite(trial_population)).or.any(trial_population<0))return
             if(any(abs(sum(trial_population,dim=1)-sum(population,dim=1))> &
                  512*epsilon(1d0)*size(population,1)*max(sum(population,dim=1),tiny(1d0))))return
             if(any(.not.ieee_is_finite(rate)).or.any(rate<0))return
             if(any(.not.ieee_is_finite(trial_dust_energy)).or.any(trial_dust_energy<0))return
             if(any(.not.ieee_is_finite(next_t)).or.any(next_t<0))return
             ierr=dust_ok
          else if(present(material_dispatch))then
             if(present(gas_energy))then
             call material_dispatch(primary+guess/dt,density,dust_energy,heat_capacity,table%log_t, &
                  table%power,table%band,dispatch_u,allocated(table%material_u),dt,table%background, &
                  table%background_temperature,material_tolerance,rate,next_t,trial_dust_energy,ierr, &
                  gas_energy,gas_capacity,conductance,exchange,cell_material_u,cell_weights, &
                  table%optical_power,table%optical_band)
             else
             call material_dispatch(primary+guess/dt,density,dust_energy,heat_capacity,table%log_t, &
                  table%power,table%band,dispatch_u,allocated(table%material_u),dt,table%background, &
                  table%background_temperature,material_tolerance,rate,next_t,trial_dust_energy,ierr, &
                  cell_material_u=cell_material_u,cell_weights=cell_weights, &
                  basis_power=table%optical_power,basis_band=table%optical_band)
             endif
          else
             call transient_emission(table,primary+guess/dt,density,dt,dust_energy,heat_capacity, &
                  rate,next_t,trial_dust_energy,ierr,material_tolerance)
          endif
       else
          call emission(table, primary+guess/dt, density, rate, next_t, ierr)
       end if
       if (ierr/=dust_ok) return
       absorbed=0; new_total=0
       if(present(population))spectral_absorbed=0
       if(moving)then
          if(.not.present(moving_material_dispatch))then
             ! Common-temperature bulk emission: all phase spectra share the
             ! same Planck factor within a group, so split by absorption, not
             ! scattering, opacity. This is NOT valid for a PAH spectrum.
             phase_rate=0
             do i=1,nc
                do g=1,ng
                   phase_scale=sum(phase_absorption(:,g,i))
                   if(phase_scale>0)then
                      phase_rate(:,g,i)=rate(g,i)*phase_absorption(:,g,i)/phase_scale
                   else if(rate(g,i)>0)then
                      ierr=dust_err_state;return
                   endif
                enddo
             enddo
          endif
          call moving_candidate(transported,phase_rate,phase_absorption,phase_density,phase_momentum,beta_guess, &
               direction,weight,table%energy,dt,c_hat,candidate,phase_heat,phase_events,momentum_next,beta_next, &
               work_next,absorbed,ierr)
          if(ierr/=dust_ok)return
          do i=1,nc
             do d=1,nd
                new_total=new_total+sum(candidate(:,d,i))*weight(d)*volume
             enddo
          enddo
          moving_residual=0
          do i=1,nc
             phase_scale=max(primary(i)*dt+sum(phase_heat(:,:,i)),tiny(1d0))
             moving_residual=max(moving_residual,sum(abs(phase_heat(:,:,i)-heat_guess(:,:,i)))/phase_scale)
             do g=1,ng
                moving_residual=max(moving_residual, &
                     sum(abs(phase_events(:,g,i)-event_guess(:,g,i)))*table%energy(g)*ev_erg/phase_scale)
             enddo
             moving_residual=max(moving_residual,maxval(abs(beta_next(:,:,i)-beta_guess(:,:,i)))/ &
                  snrt_scatter_beta_limit)
          enddo
       else
       if(present(absorb_dispatch))then
          call absorb_dispatch(transported,transmit,loss,response,rate,weight,dt,sum_w,candidate,absorbed,ierr)
          if(ierr/=dust_ok)return
       endif
       do i=1,nc
          do d=1,nd
             do g=1,ng
                if(present(absorb_dispatch))exit
                source=dt*rate(g,i)/sum_w
                candidate(g,d,i)=transported(g,d,i)*transmit(g,i)+source*response(g,i)
                absorbed(i)=absorbed(i)+weight(d)*(transported(g,d,i)*loss(g,i)+source*(1-response(g,i)))
                if(present(population))spectral_absorbed(g,i)=spectral_absorbed(g,i)+ &
                     weight(d)*(transported(g,d,i)*loss(g,i)+source*(1-response(g,i)))
             end do
             new_total=new_total+sum(candidate(:,d,i))*weight(d)*volume
          end do
       end do
       endif ! moving versus static interaction
       ierr=dust_err_state
       if (.not.all(ieee_is_finite(candidate)) .or. any(candidate<0) .or. .not.ieee_is_finite(new_total)) return
       if (.not.all(ieee_is_finite(absorbed)) .or. any(absorbed<0)) return
       balance=new_total-old_total+trial%escaped_erg+trial%interface_erg-trial%primary_erg
       if(transient)balance=balance+sum(trial_dust_energy-dust_energy)*volume
       if(moving)balance=balance+sum(work_next)*volume
       if(present(gas_energy))then
          if(any(.not.ieee_is_finite(exchange)).or.any(gas_energy-exchange<0))return
          balance=balance-sum(exchange)*volume
          scale=max(scale,sum(abs(exchange))*volume)
       endif
       trial%balance_relative=abs(balance)/scale
       trial%local_relative=maxval(abs(absorbed-guess)/max(primary*dt+absorbed,tiny(scale)))
       if(moving)trial%local_relative=moving_residual
       if(present(population).and..not.moving)then
          ! Equal total heating does not imply equal PAH excitation: require
          ! spectral fixed-point convergence, not just cancellation in sum(g).
          do i=1,nc
             trial%local_relative=max(trial%local_relative, &
                  sum(abs(spectral_absorbed(:,i)-spectral_guess(:,i)))/ &
                  max(primary(i)*dt+absorbed(i),tiny(scale)))
          enddo
       endif
       if (max(trial%balance_relative,trial%local_relative)<=tolerance) then
          if(moving.and.present(phase_scattering))then
             ! First-order split AFTER absorption/emission, wholly inside
             ! this accepted transaction. Only local interaction energy is
             ! used for impulses; transported/escaped energy is never force.
             allocate(scatter_e(nd,ng),scatter_tau(nb,ng),scatter_work(nb))
             do i=1,nc
                do d=1,nd
                   scatter_e(d,:)=candidate(:,d,i)*weight(d)
                enddo
                scatter_tau=phase_scattering(:,:,i)*(c_hat*dt)
                scatter_work=0
                call snrt_moving_scatter_energy_cell(scatter_e,momentum_next(:,:,i),phase_density(:,i), &
                     scatter_tau,direction,weight,light_c,scatter_work,moving_status)
                if(moving_status/=0)then
                   ierr=dust_err_convergence;return
                endif
                work_next(:,i)=work_next(:,i)+scatter_work
                do d=1,nd
                   candidate(:,d,i)=scatter_e(d,:)/weight(d)
                enddo
             enddo
             new_total=0
             do i=1,nc
                do d=1,nd
                   new_total=new_total+sum(candidate(:,d,i))*weight(d)*volume
                enddo
             enddo
             scatter_balance=new_total-old_total+trial%escaped_erg+trial%interface_erg-trial%primary_erg+ &
                  (sum(trial_dust_energy-dust_energy)+sum(work_next))*volume
             if(present(gas_energy))scatter_balance=scatter_balance-sum(exchange)*volume
             trial%balance_relative=abs(scatter_balance)/scale
             if(trial%balance_relative>tolerance)then
                ierr=dust_err_convergence;return
             endif
          endif
          do g=1,ng
             emitted_photons(g,:)=rate(g,:)*dt/(table%energy(g)*ev_erg)
          end do
          if (.not.all(ieee_is_finite(photons+emitted_photons))) return
          trial%iterations=iteration; trial%absorbed_erg=sum(absorbed)*volume
          if (.not.ieee_is_finite(trial%absorbed_erg)) return
          energy=candidate; temperature=next_t; photons=photons+emitted_photons
          if(transient)dust_energy=trial_dust_energy
          if(moving)then
             phase_momentum=momentum_next;phase_work=work_next
             trial%mechanical_erg=sum(work_next)*volume
          endif
          if(present(population))population=trial_population
          if(present(gas_energy))then
             gas_energy=gas_energy-exchange
             gas_transfer=exchange
          endif
          diagnostics=trial
          ierr=dust_ok
          return
       end if
       if(moving)then
          do i=1,nc
             if(direct_update(i))then
                heat_guess(:,:,i)=phase_heat(:,:,i)
                event_guess(:,:,i)=phase_events(:,:,i)
             else
                heat_guess(:,:,i)=.5d0*(heat_guess(:,:,i)+phase_heat(:,:,i))
                event_guess(:,:,i)=.5d0*(event_guess(:,:,i)+phase_events(:,:,i))
             endif
          enddo
          beta_guess=.5d0*(beta_guess+beta_next)
          guess=sum(sum(heat_guess,dim=1),dim=1)
       else
       where(direct_update)
          guess=absorbed
       elsewhere
          guess=.5d0*(guess+absorbed)
       endwhere
       if(present(population))then
          do i=1,nc
             if(direct_update(i))then
                spectral_guess(:,i)=spectral_absorbed(:,i)
             else
                spectral_guess(:,i)=.5d0*(spectral_guess(:,i)+spectral_absorbed(:,i))
             endif
          enddo
       endif
       endif ! moving iteration
    end do
    ierr=dust_err_convergence
  end subroutine

  subroutine moving_candidate(incident,rate,alpha,rho,pold,beta,direction,weight,group_ev,dt,chat, &
       candidate,heat,events,pnext,beta_next,work,absorbed,ierr)
    ! Local interaction only; incident has ALREADY been spatially transported.
    ! Frozen-opacity exact attenuation with constant trial emission, solved
    ! simultaneously with material energy and midpoint phase velocities by
    ! the enclosing fixed-point iteration. No arbitrary angular work debit.
    real(real64),intent(in)::incident(:,:,:),rate(:,:,:),alpha(:,:,:),rho(:,:),pold(:,:,:),beta(:,:,:)
    real(real64),intent(in)::direction(:,:),weight(:),group_ev(:),dt,chat
    real(real64),intent(out)::candidate(:,:,:),heat(:,:,:),events(:,:,:),pnext(:,:,:),beta_next(:,:,:)
    real(real64),intent(out)::work(:,:),absorbed(:)
    integer,intent(out)::ierr
    real(real64)::q(size(rho,1),size(weight)),a(size(rho,1),size(weight)),depth(size(rho,1))
    real(real64)::emit(size(rho,1)),capture(size(rho,1)),dp(3,size(rho,1))
    real(real64)::total_depth,keep,loss,response,reabsorbed,source,old,removed,beta_actual(3),vnext(3)
    integer::i,b,g,d,nb
    ierr=dust_err_state;nb=size(rho,1);heat=0;events=0;work=0;absorbed=0
    pnext=pold;beta_next=0
    do i=1,size(rho,2)
       dp=0
       do b=1,nb
          if(norm2(beta(:,b,i))>snrt_scatter_beta_limit)return
          q(b,:)=1-matmul(beta(:,b,i),direction)
          a(b,:)=weight/q(b,:)**2;a(b,:)=a(b,:)/sum(a(b,:))
          if(rho(b,i)==0.and.any(rate(b,:,i)/=0))return
       enddo
       do g=1,size(group_ev)
          do d=1,size(weight)
             old=incident(g,d,i)*weight(d)
             depth=alpha(:,g,i)*(chat*dt)*q(:,d);total_depth=sum(depth)
             if(.not.ieee_is_finite(total_depth).or.total_depth<0)return
             keep=exp(-total_depth)
             if(total_depth<1d-4)then
                loss=total_depth*(1-total_depth/2+total_depth**2/6-total_depth**3/24)
                reabsorbed=total_depth*(.5d0-total_depth/6+total_depth**2/24-total_depth**3/120)
                response=1-reabsorbed
             else
                loss=1-keep;response=loss/total_depth;reabsorbed=1-response
             endif
             ! rate is the emitted comoving energy, not lab energy: sum(q*emit)
             ! equals dt*rate for each phase. Emission recoil is not omitted.
             emit=dt*rate(:,g,i)*a(:,d)/q(:,d);source=sum(emit)
             candidate(g,d,i)=(keep*old+source*response)/weight(d)
             removed=old*loss+source*reabsorbed
             capture=0
             if(total_depth>0)capture=removed*(depth/total_depth)
             heat(:,g,i)=heat(:,g,i)+q(:,d)*capture
             events(:,g,i)=events(:,g,i)+capture/(group_ev(g)*ev_erg)
             absorbed(i)=absorbed(i)+sum(capture)
             do b=1,nb
                dp(:,b)=dp(:,b)+(capture(b)-emit(b))*direction(:,d)/light_c
             enddo
          enddo
       enddo
       do b=1,nb
          if(rho(b,i)==0)then
             if(any(dp(:,b)/=0))return
             cycle
          endif
          pnext(:,b,i)=pold(:,b,i)+dp(:,b)
          vnext=pnext(:,b,i)/(rho(b,i)*light_c)
          if(.not.all(ieee_is_finite(vnext)))return
          if(norm2(vnext)>snrt_scatter_beta_limit)return
          beta_actual=.5d0*(pold(:,b,i)/(rho(b,i)*light_c)+vnext)
          beta_next(:,b,i)=beta_actual
          work(b,i)=dot_product(beta_actual,dp(:,b))*light_c
       enddo
    enddo
    if(any(.not.ieee_is_finite(heat)).or.any(.not.ieee_is_finite(events)))return
    if(any(.not.ieee_is_finite(work)).or.any(.not.ieee_is_finite(pnext)))return
    ierr=dust_ok
  end subroutine moving_candidate
end module

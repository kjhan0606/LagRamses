! Contract audit routines for canonical stellar-yield tables.
!
! These checks do not replace a scientific comparison with the source papers.
! They do enforce the numerical contract used by the runtime: finite physical
! values, non-negative cumulative material/energy, mass closure, monotonic
! cumulative histories, unique coordinates, and (when requested) a complete
! Cartesian mass-metallicity-age grid for the required channels.

module stellar_yield_audit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, channel_wind, channel_snii, channel_agb, valid_high_mass_choice
  use stellar_yield_tables, only: stellar_yield_table_t
  use stellar_enrichment_config, only: default_imf_id, population_model_id, configured_binary_fraction, &
       configured_imf_mass_min, configured_imf_mass_max, high_mass_model, high_mass_max_remnant_adjust_fraction
  implicit none

  private
  integer, parameter, public :: yield_audit_ok = 0
  integer, parameter, public :: yield_audit_err_table = 1
  integer, parameter, public :: yield_audit_err_value = 2
  integer, parameter, public :: yield_audit_err_mass = 4
  integer, parameter, public :: yield_audit_err_monotonic = 8
  integer, parameter, public :: yield_audit_err_nonfinite = 16
  integer, parameter, public :: yield_audit_err_duplicate = 32
  integer, parameter, public :: yield_audit_err_grid = 64
  integer, parameter, public :: yield_audit_err_energy_monotonic = 128
  integer, parameter, public :: yield_audit_err_remnant_ownership = 256

  public :: audit_yield_table
  public :: resolve_high_mass_endpoint
  public :: prepare_high_mass_history

  type, public :: high_mass_endpoint_t
     real(stellar_dp) :: initial_mass = 0.0_stellar_dp
     real(stellar_dp) :: wind_mass = 0.0_stellar_dp
     real(stellar_dp) :: terminal_mass = 0.0_stellar_dp
     ! Baryonic residual-mass bookkeeping, NOT a gravitational BH mass.
     real(stellar_dp) :: remnant_mass = 0.0_stellar_dp
     real(stellar_dp) :: wind_elements(n_stellar_elements) = 0.0_stellar_dp
     real(stellar_dp) :: terminal_elements(n_stellar_elements) = 0.0_stellar_dp
     real(stellar_dp) :: wind_energy = 0.0_stellar_dp
     real(stellar_dp) :: terminal_energy = 0.0_stellar_dp
     real(stellar_dp) :: wind_momentum(3) = 0.0_stellar_dp
     real(stellar_dp) :: terminal_momentum(3) = 0.0_stellar_dp
  end type high_mass_endpoint_t

contains

  subroutine prepare_high_mass_history(table, filename, ierr)
    type(stellar_yield_table_t), intent(inout) :: table
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ierr
    integer, parameter :: max_nodes=512
    integer :: version, node_count, input_imf_id, input_population_id
    integer :: terminal_outcome(max_nodes), i, j, k, w, s, unit, status
    real(stellar_dp) :: mass_msun(max_nodes), metallicity(max_nodes), terminal_age_yr(max_nodes)
    real(stellar_dp) :: input_binary_fraction, input_imf_min, input_imf_max, delta, history_min, history_max
    real(stellar_dp) :: mass_domain_min, mass_domain_max
    integer :: agb_non_co_count, agb_non_co_kind(max_nodes), found
    real(stellar_dp) :: agb_non_co_mass(max_nodes), agb_non_co_z(max_nodes)
    integer :: agb_wind_jump_count
    real(stellar_dp) :: agb_wind_jump_mass(max_nodes),agb_wind_jump_z(max_nodes),agb_wind_jump_fraction(max_nodes)
    logical :: wind_history,declared(max_nodes)
    logical :: net_channel_available(n_stellar_channels)
    character(len=128) :: model_id, wind_source_id, terminal_source_id, model_coordinates
    character(len=32) :: timing_policy
    character(len=32) :: metallicity_policy
    character(len=32) :: net_yield_policy, agb_release_policy
    type(stellar_yield_table_t) :: trial
    type(high_mass_endpoint_t) :: raw, resolved
    namelist /stellar_high_mass_history/ version, node_count, model_id, wind_source_id, &
         terminal_source_id, model_coordinates, timing_policy, input_imf_id, input_population_id, &
         input_binary_fraction, input_imf_min, input_imf_max, mass_msun, metallicity, &
         terminal_age_yr, terminal_outcome, metallicity_policy, net_yield_policy, agb_release_policy, &
         mass_domain_min, mass_domain_max, agb_non_co_count, agb_non_co_mass, agb_non_co_z, agb_non_co_kind, &
         net_channel_available,agb_wind_jump_count,agb_wind_jump_mass,agb_wind_jump_z,agb_wind_jump_fraction

    ierr=yield_audit_err_table
    if (allocated(table%hm_mass)) return ! A second application would double-transform a model.
    call audit_yield_table(table,1d-10,status,exact_source_coordinates=.true.)
    if(status/=0)return
    version=0; node_count=0; input_imf_id=-1; input_population_id=-1
    input_binary_fraction=-1; input_imf_min=-1; input_imf_max=-1
    model_id=''; wind_source_id=''; terminal_source_id=''; model_coordinates=''; timing_policy=''
    metallicity_policy='exact_nodes'
    net_yield_policy='supplied'
    agb_release_policy='cumulative_linear'
    mass_domain_min=-1d0;mass_domain_max=-1d0
    agb_non_co_count=0;agb_non_co_mass=-1;agb_non_co_z=-1;agb_non_co_kind=-1
    agb_wind_jump_count=0;agb_wind_jump_mass=-1;agb_wind_jump_z=-1;agb_wind_jump_fraction=-1
    net_channel_available=.false.
    mass_msun=-1; metallicity=-1; terminal_age_yr=-1; terminal_outcome=-1
    open(newunit=unit,file=filename,status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,nml=stellar_high_mass_history,iostat=status)
    close(unit)
    if(status/=0.or.version<1.or.version>3.or.node_count<2.or.node_count>max_nodes)return
    ! v1 remains the original [40,120] contract; v2 explicitly includes the
    ! source-supported ordinary CCSN branch, never an 8--13 extrapolation.
    history_min=40d0;history_max=120d0
    if(version==2)history_min=13d0
    if(version==3)then
       ! A source-specific comparison may cover only part of the CCSN
       ! domain. Admission never extrapolates it to the rest of the IMF.
       if(.not.all(ieee_is_finite([mass_domain_min,mass_domain_max])))return
       if(mass_domain_min<8d0.or.mass_domain_max>120d0.or.mass_domain_max<=mass_domain_min)return
       history_min=mass_domain_min;history_max=mass_domain_max
    else
       if(mass_domain_min/=-1d0.or.mass_domain_max/=-1d0)return
    endif
    if(min(len_trim(model_id),len_trim(wind_source_id),len_trim(terminal_source_id), &
         len_trim(model_coordinates))==0)return
    if(timing_policy/='wind_linear_terminal_step')return
    if(metallicity_policy/='exact_nodes'.and.metallicity_policy/='linear_Z_cumulative_mixture')return
    if(net_yield_policy/='supplied'.and.net_yield_policy/='unavailable_diagnostic_zero'.and. &
         net_yield_policy/='channel_mask')return
    wind_history=agb_release_policy=='wind_history_terminal_remnant'
    if(agb_release_policy/='cumulative_linear'.and.agb_release_policy/='terminal_step'.and..not.wind_history)return
    if(agb_wind_jump_count<0.or.agb_wind_jump_count>max_nodes)return
    if(wind_history.neqv.(agb_wind_jump_count>0))return
    if(net_yield_policy=='unavailable_diagnostic_zero')then
       if(any(table%net_yield/=0))return
    endif
    if(net_yield_policy=='channel_mask')then
       if(.not.any(net_channel_available))return
       do i=1,table%n_rows
          if(.not.net_channel_available(table%channel(i)).and.any(table%net_yield(i,:)/=0))return
       enddo
    endif
    if(input_imf_id/=default_imf_id.or.input_population_id/=population_model_id)return
    if(.not.all(ieee_is_finite([input_binary_fraction,input_imf_min,input_imf_max])))return
    if(input_binary_fraction/=configured_binary_fraction.or.input_imf_min/=configured_imf_mass_min.or. &
         input_imf_max/=configured_imf_mass_max)return
    if(.not.all(ieee_is_finite(mass_msun(:node_count))).or. &
         .not.all(ieee_is_finite(metallicity(:node_count))).or. &
         .not.all(ieee_is_finite(terminal_age_yr(:node_count))))return
    if(any(mass_msun(:node_count)<history_min).or.any(mass_msun(:node_count)>history_max).or. &
         any(metallicity(:node_count)<0).or.any(terminal_age_yr(:node_count)<=0))return
    if(any(terminal_outcome(:node_count)<0).or.any(terminal_outcome(:node_count)>1))return
    ! Each exact-Z branch must cover the declared domain in source-node order.
    do i=1,node_count
       if(i==1)then
          if(mass_msun(i)/=history_min)return
       else if(metallicity(i)==metallicity(i-1))then
          if(mass_msun(i)<=mass_msun(i-1))return
       else
          if(metallicity(i)<=metallicity(i-1).or.mass_msun(i-1)/=history_max.or.mass_msun(i)/=history_min)return
       endif
    enddo
    if(mass_msun(node_count)/=history_max)return
    trial=table
    allocate(trial%hm_mass(node_count),trial%hm_z(node_count),trial%hm_age(node_count), &
         trial%hm_remnant(node_count),trial%hm_adjustment(node_count), &
         trial%hm_wind_row(node_count),trial%hm_terminal_row(node_count))
    trial%hm_mass=mass_msun(:node_count); trial%hm_z=metallicity(:node_count)
    trial%hm_age=terminal_age_yr(:node_count)*1d-9
    trial%high_mass_identity=[model_id,wind_source_id,terminal_source_id,model_coordinates]
    trial%high_mass_wind_only=high_mass_model=='wind_only_collapse'
    trial%high_mass_linear_z=metallicity_policy=='linear_Z_cumulative_mixture'
    trial%net_yield_channel_available=net_yield_policy=='supplied'
    if(net_yield_policy=='channel_mask')trial%net_yield_channel_available=net_channel_available
    trial%net_yield_diagnostic_unavailable=.not.all(trial%net_yield_channel_available)
    if(trial%high_mass_linear_z.and.maxval(trial%hm_z)<=minval(trial%hm_z))return
    do i=1,node_count
       w=0; s=0
       do j=1,table%n_rows
          if(table%initial_mass(j)/=mass_msun(i).or.table%birth_metallicity(j)/=metallicity(i))cycle
          if(table%age_gyr(j)/=trial%hm_age(i))cycle
          if(table%channel(j)==channel_wind)w=j
          if(table%channel(j)==channel_snii)s=j
       enddo
       if(w==0.or.s==0)return
       raw=high_mass_endpoint_t()
       raw%initial_mass=mass_msun(i)
       raw%wind_mass=table%returned_mass(w); raw%terminal_mass=table%returned_mass(s)
       raw%remnant_mass=table%remnant_mass(s)
       raw%wind_elements=table%ejected_mass(w,:); raw%terminal_elements=table%ejected_mass(s,:)
       raw%wind_energy=table%energy(w); raw%terminal_energy=table%energy(s)
       raw%wind_momentum=table%momentum(w,:); raw%terminal_momentum=table%momentum(s,:)
       if(terminal_outcome(i)==0.and.(raw%terminal_mass/=0.or.raw%terminal_energy/=0.or. &
            any(raw%terminal_momentum/=0)))return
       if(terminal_outcome(i)==1.and.raw%terminal_mass<=0)return
       if(mass_msun(i)<40d0)then
          ! High-mass alternatives must not suppress ordinary CCSN. These
          ! endpoints must close without any mixed-source remnant repair.
          call resolve_high_mass_endpoint(raw,'source_consistent',wind_source_id==terminal_source_id, &
               0d0,resolved,delta,status,allow_ordinary=.true.)
       else
          call resolve_high_mass_endpoint(raw,high_mass_model,wind_source_id==terminal_source_id, &
               high_mass_max_remnant_adjust_fraction,resolved,delta,status)
       endif
       if(status/=0)then
          ierr=status
          return
       endif
       trial%hm_wind_row(i)=w; trial%hm_terminal_row(i)=s
       trial%hm_remnant(i)=resolved%remnant_mass; trial%hm_adjustment(i)=delta
       ! Input wind starts at zero; both channels have a complete endpoint.
       ! No terminal leakage before the specified event and no new release
       ! after it. The earlier generic linear SN ramp is explicitly rejected.
       k=0
       do j=1,table%n_rows
          if(table%initial_mass(j)/=mass_msun(i).or.table%birth_metallicity(j)/=metallicity(i))cycle
          if(table%channel(j)==channel_wind)then
             if(table%remnant_mass(j)/=0)return
             if(table%age_gyr(j)==0)k=k+1
             if(table%age_gyr(j)>trial%hm_age(i))then
                if(.not.same_payload(table,j,w))return
             endif
          else if(table%channel(j)==channel_snii)then
             if(table%age_gyr(j)<trial%hm_age(i))then
                if(table%returned_mass(j)/=0.or.table%remnant_mass(j)/=0.or.table%energy(j)/=0.or. &
                     any(table%momentum(j,:)/=0).or.any(table%ejected_mass(j,:)/=0).or. &
                     any(table%net_yield(j,:)/=0))return
             else
                if(.not.same_payload(table,j,s))return
             endif
          endif
       enddo
       if(k/=1)return
    enddo
    ! A canonical row above the seam cannot escape the declared node map.
    do j=1,table%n_rows
       if(table%channel(j)/=channel_wind.and.table%channel(j)/=channel_snii)cycle
       if(table%initial_mass(j)<history_min)cycle
       if(.not.any(trial%hm_mass==table%initial_mass(j).and.trial%hm_z==table%birth_metallicity(j)))return
    enddo
    if(agb_release_policy=='terminal_step'.or.wind_history)then
       call prepare_agb_terminal_rows(trial,status,wind_history)
       if(status/=0)return
    endif
    if(wind_history)then
       if(size(trial%agb_terminal_row)>max_nodes)return
       declared=.false.
       do i=1,agb_wind_jump_count
          if(.not.all(ieee_is_finite([agb_wind_jump_mass(i),agb_wind_jump_z(i),agb_wind_jump_fraction(i)])))return
          if(agb_wind_jump_fraction(i)<0.or.agb_wind_jump_fraction(i)>1)return
          found=0
          do j=1,size(trial%agb_terminal_row)
             k=trial%agb_terminal_row(j)
             if(trial%initial_mass(k)/=agb_wind_jump_mass(i).or.trial%birth_metallicity(k)/=agb_wind_jump_z(i))cycle
             if(declared(j))return
             trial%agb_terminal_jump_fraction(j)=agb_wind_jump_fraction(i);declared(j)=.true.;found=found+1
          enddo
          if(found/=1)return
       enddo
       ! No cumulative wind knot may exceed the LEFT limit of the terminal
       ! jump. Otherwise interpolation would subtract material/energy.
       do j=1,size(trial%agb_terminal_row)
          k=trial%agb_terminal_row(j);delta=1-trial%agb_terminal_jump_fraction(j)
          do i=1,trial%n_rows
             if(trial%channel(i)/=channel_agb.or.trial%initial_mass(i)/=trial%initial_mass(k).or. &
                  trial%birth_metallicity(i)/=trial%birth_metallicity(k))cycle
             if(trial%age_gyr(i)>=trial%age_gyr(k))cycle
             if(trial%returned_mass(i)>delta*trial%returned_mass(k)+1d-10*trial%initial_mass(k))return
             if(trial%energy(i)>delta*trial%energy(k)+1d-10*max(1d0,trial%energy(k)))return
             if(any(trial%ejected_mass(i,:)>delta*trial%ejected_mass(k,:)+1d-10*trial%initial_mass(k)))return
             if(trial%returned_mass(i)-sum(trial%ejected_mass(i,:))> &
                  delta*(trial%returned_mass(k)-sum(trial%ejected_mass(k,:)))+1d-10*trial%initial_mass(k))return
          enddo
       enddo
    endif
    if(agb_non_co_count<0.or.agb_non_co_count>max_nodes)return
    if(agb_non_co_count>0)then
       if(.not.allocated(trial%agb_terminal_row))return
       if(.not.all(ieee_is_finite(agb_non_co_mass(:agb_non_co_count))).or. &
            .not.all(ieee_is_finite(agb_non_co_z(:agb_non_co_count))))return
       do i=1,agb_non_co_count
          if(agb_non_co_kind(i)<1.or.agb_non_co_kind(i)>2)return
          found=0
          do j=1,size(trial%agb_terminal_row)
             k=trial%agb_terminal_row(j)
             if(trial%initial_mass(k)/=agb_non_co_mass(i).or.trial%birth_metallicity(k)/=agb_non_co_z(i))cycle
             if(trial%agb_remnant_kind(j)/=0)return ! Duplicate declaration.
             trial%agb_remnant_kind(j)=agb_non_co_kind(i);found=found+1
          enddo
          if(found/=1)return
       enddo
    endif
    trial%high_mass_ready=.true.
    table=trial
    ierr=yield_audit_ok
  end subroutine prepare_high_mass_history

  subroutine prepare_agb_terminal_rows(table,ierr,allow_wind)
    ! The first remnant row defines the supplied terminal age when wind
    ! histories are selected; legacy terminal-only inputs must remain zero
    ! before that age. Never infer lifetimes or duplicate wind ownership.
    type(stellar_yield_table_t),intent(inout)::table
    integer,intent(out)::ierr
    logical,intent(in),optional::allow_wind
    logical::wind
    integer::i,j,k,n,first,zero
    integer::rows(table%n_rows)
    real(stellar_dp)::m,z
    ierr=yield_audit_err_table;n=0
    wind=.false.
    if(present(allow_wind))wind=allow_wind
    do i=1,table%n_rows
       if(table%channel(i)/=channel_agb)cycle
       m=table%initial_mass(i);z=table%birth_metallicity(i)
       if(m>8d0)return ! Ordinary AGB/WD route only, no SAGB fate inference.
       if(n>0)then
          if(any(table%initial_mass(rows(:n))==m.and.table%birth_metallicity(rows(:n))==z))cycle
       endif
       first=0;zero=0
       do j=1,table%n_rows
          if(table%initial_mass(j)/=m.or.table%birth_metallicity(j)/=z)cycle
          if(table%channel(j)==channel_wind)then
             if(.not.zero_payload(table,j))return
          endif
          if(table%channel(j)/=channel_agb)cycle
          if(table%age_gyr(j)==0d0.and.zero_payload(table,j))zero=j
          if(zero_payload(table,j))cycle
          if(wind.and.table%remnant_mass(j)==0)cycle
          if(first==0)then
             first=j
          else if(table%age_gyr(j)<table%age_gyr(first))then
             first=j
          endif
       enddo
       if(first==0.or.zero==0)return
       if(table%age_gyr(first)<=0.or.table%returned_mass(first)<=0.or.table%remnant_mass(first)<=0)return
       if(abs(table%returned_mass(first)+table%remnant_mass(first)-m)>1d-10*max(1d0,m))return
       do k=1,table%n_rows
          if(table%channel(k)/=channel_agb.or.table%initial_mass(k)/=m.or.table%birth_metallicity(k)/=z)cycle
          if(table%age_gyr(k)<table%age_gyr(first))then
             if(wind)then
                if(table%remnant_mass(k)/=0)return
             else
                if(.not.zero_payload(table,k))return
             endif
          else
             if(.not.same_payload(table,k,first))return
          endif
       enddo
       n=n+1;rows(n)=first
    enddo
    if(n==0)return
    table%agb_terminal_row=rows(:n)
    allocate(table%agb_remnant_kind(n));table%agb_remnant_kind=0
    allocate(table%agb_terminal_jump_fraction(n));table%agb_terminal_jump_fraction=1d0
    ierr=yield_audit_ok
  end subroutine prepare_agb_terminal_rows

  logical function zero_payload(table,i)
    type(stellar_yield_table_t),intent(in)::table
    integer,intent(in)::i
    zero_payload=table%returned_mass(i)==0.and.table%remnant_mass(i)==0.and.table%energy(i)==0.and. &
         all(table%momentum(i,:)==0).and.all(table%ejected_mass(i,:)==0).and.all(table%net_yield(i,:)==0)
  end function zero_payload

  logical function same_payload(table,i,j)
    type(stellar_yield_table_t), intent(in) :: table
    integer,intent(in)::i,j
    same_payload=table%returned_mass(i)==table%returned_mass(j).and. &
         table%remnant_mass(i)==table%remnant_mass(j).and.table%energy(i)==table%energy(j).and. &
         all(table%momentum(i,:)==table%momentum(j,:)).and. &
         all(table%ejected_mass(i,:)==table%ejected_mass(j,:)).and. &
         all(table%net_yield(i,:)==table%net_yield(j,:))
  end function same_payload

  subroutine resolve_high_mass_endpoint(raw, preset, same_source, adjustment_limit, &
       resolved, remnant_adjustment, ierr, allow_ordinary)
    ! Resolve a COMPLETE per-initial-star endpoint, before IMF integration.
    ! Never use this on a cumulative row at an intermediate age. Timing,
    ! isotope projection, source identities and fate classification belong
    ! to the source adapter; a zero Wind table is not filled by this routine.
    type(high_mass_endpoint_t), intent(in) :: raw
    character(len=*), intent(in) :: preset
    logical, intent(in) :: same_source
    real(stellar_dp), intent(in) :: adjustment_limit
    type(high_mass_endpoint_t), intent(out) :: resolved
    real(stellar_dp), intent(out) :: remnant_adjustment
    integer, intent(out) :: ierr
    logical, intent(in), optional :: allow_ordinary
    type(high_mass_endpoint_t) :: trial
    real(stellar_dp), parameter :: tolerance = 1.0e-10_stellar_dp
    real(stellar_dp) :: residual, scale, minimum_mass

    ! Failure has no publishable material or energy, including the correction.
    resolved = high_mass_endpoint_t()
    remnant_adjustment = 0.0_stellar_dp
    ierr = yield_audit_err_value
    if (.not. valid_high_mass_choice(preset, adjustment_limit)) return
    if (.not. all(ieee_is_finite((/raw%initial_mass, raw%wind_mass, &
         raw%terminal_mass, raw%remnant_mass, raw%wind_energy, raw%terminal_energy/)))) return
    if (.not. all(ieee_is_finite(raw%wind_elements)) .or. &
         .not. all(ieee_is_finite(raw%terminal_elements)) .or. &
         .not. all(ieee_is_finite(raw%wind_momentum)) .or. &
         .not. all(ieee_is_finite(raw%terminal_momentum))) return
    minimum_mass=40d0
    if(present(allow_ordinary))then
       if(allow_ordinary)then
          if(preset/='source_consistent')return
          minimum_mass=8d0
       endif
    endif
    if (raw%initial_mass < minimum_mass .or. raw%initial_mass > 120.0_stellar_dp) return
    if (min(raw%wind_mass, raw%terminal_mass, raw%remnant_mass, &
         raw%wind_energy, raw%terminal_energy) < 0.0_stellar_dp) return
    if (minval(raw%wind_elements) < 0.0_stellar_dp .or. &
         minval(raw%terminal_elements) < 0.0_stellar_dp) return
    if (raw%remnant_mass > raw%initial_mass) return
    scale = raw%initial_mass
    if (sum(raw%wind_elements) > raw%wind_mass + tolerance*scale .or. &
         sum(raw%terminal_elements) > raw%terminal_mass + tolerance*scale) return
    ! These energies are non-bulk source energy; the existing deposition
    ! bridge separately accounts for each channel's directed kinetic energy.
    if (raw%wind_mass == 0.0_stellar_dp .and. &
         (any(raw%wind_momentum /= 0.0_stellar_dp) .or. raw%wind_energy /= 0.0_stellar_dp)) return
    if (raw%terminal_mass == 0.0_stellar_dp .and. &
         (any(raw%terminal_momentum /= 0.0_stellar_dp) .or. raw%terminal_energy /= 0.0_stellar_dp)) return
    trial = raw
    ierr = yield_audit_err_mass
    select case (trim(preset))
    case ('source_consistent')
       if (.not. same_source) return
    case ('wind_only_collapse')
       ! Explicit alternative, not an assertion about all massive stars:
       ! retain supplied wind and its composition/energy, suppress terminal
       ! ejecta/energy, put all remaining baryonic mass in the remnant.
       trial%terminal_mass = 0.0_stellar_dp
       trial%terminal_elements = 0.0_stellar_dp
       trial%terminal_energy = 0.0_stellar_dp
       trial%terminal_momentum = 0.0_stellar_dp
       trial%remnant_mass = raw%initial_mass - raw%wind_mass
    case ('mixed_remnant')
       trial%remnant_mass = raw%initial_mass - raw%wind_mass - raw%terminal_mass
       if (abs(trial%remnant_mass-raw%remnant_mass) > adjustment_limit*scale) return
    end select
    if (trial%remnant_mass < 0.0_stellar_dp) return
    residual = trial%initial_mass-trial%wind_mass-trial%terminal_mass-trial%remnant_mass
    if (abs(residual) > tolerance*scale) return
    resolved = trial
    remnant_adjustment = trial%remnant_mass-raw%remnant_mass
    ierr = yield_audit_ok
  end subroutine resolve_high_mass_endpoint

  subroutine audit_yield_table(table, tolerance, ierr, require_complete, &
       terminal_remnant_owner, required_channels, exact_source_coordinates)
    type(stellar_yield_table_t), intent(in) :: table
    real(stellar_dp), intent(in) :: tolerance
    integer, intent(out) :: ierr
    logical, intent(in), optional :: require_complete
    logical, intent(in), optional :: terminal_remnant_owner(:)
    logical, intent(in), optional :: required_channels(n_stellar_channels)
    logical, intent(in), optional :: exact_source_coordinates

    real(stellar_dp) :: tol, coordinate_tol, ejected_sum, scale
    logical :: require_grid, row_is_finite, channel_is_bad
    integer :: i, j, channel

    ierr = yield_audit_ok
    require_grid = .false.
    if (present(require_complete)) require_grid = require_complete
    tol = max(tolerance, 1.0e-12_stellar_dp)
    coordinate_tol=tol
    ! Source-node histories use exact table coordinates in the evaluator.
    ! A physical budget tolerance is not a 0.1-year resolution limit: late
    ! nuclear phases can be separated by hours. Preserve the legacy generic
    ! interpolation contract until the explicit history route is requested.
    if(table%high_mass_ready)coordinate_tol=0d0
    if(present(exact_source_coordinates))then
       if(exact_source_coordinates)coordinate_tol=0d0
    endif

    if (.not. table%loaded .or. table%n_rows <= 0) then
       ierr = yield_audit_err_table
       return
    end if
    if (.not. allocated(table%channel) .or. &
         .not. allocated(table%initial_mass) .or. &
         .not. allocated(table%birth_metallicity) .or. &
         .not. allocated(table%age_gyr) .or. &
         .not. allocated(table%returned_mass) .or. &
         .not. allocated(table%remnant_mass) .or. &
         .not. allocated(table%energy) .or. .not. allocated(table%momentum) &
         .or. .not. allocated(table%ejected_mass) .or. &
         .not. allocated(table%net_yield)) then
       ierr = yield_audit_err_table
       return
    end if
    if (table%n_rows > size(table%channel) .or. &
         table%n_rows > size(table%initial_mass) .or. &
         table%n_rows > size(table%birth_metallicity) .or. &
         table%n_rows > size(table%age_gyr) .or. &
         table%n_rows > size(table%returned_mass) .or. &
         table%n_rows > size(table%remnant_mass) .or. &
         table%n_rows > size(table%energy) .or. &
         table%n_rows > size(table%momentum, 1) .or. &
         table%n_rows > size(table%ejected_mass, 1) .or. &
         table%n_rows > size(table%net_yield, 1)) then
       ierr = yield_audit_err_table
       return
    end if
    if (present(terminal_remnant_owner)) then
       if (size(terminal_remnant_owner) /= n_stellar_channels) then
          ierr = yield_audit_err_table
          return
       end if
    end if

    do i = 1, table%n_rows
       row_is_finite = finite_row(table, i)
       if (.not. row_is_finite) then
          ierr = ior(ierr, yield_audit_err_nonfinite)
          cycle
       end if

       if (table%channel(i) < 1 .or. table%channel(i) > n_stellar_channels .or. &
            table%initial_mass(i) <= 0.0_stellar_dp .or. &
            table%birth_metallicity(i) < 0.0_stellar_dp .or. &
            table%age_gyr(i) < 0.0_stellar_dp .or. &
            table%returned_mass(i) < -tol .or. &
            table%remnant_mass(i) < -tol .or. &
            table%energy(i) < -tol .or. &
            minval(table%ejected_mass(i,:)) < -tol) then
          ierr = ior(ierr, yield_audit_err_value)
       end if
       if (abs(table%age_gyr(i)) <= tol) then
          scale = max(1.0_stellar_dp, abs(table%returned_mass(i)), &
               abs(table%remnant_mass(i)), abs(table%energy(i)), &
               maxval(abs(table%momentum(i,:))), &
               maxval(abs(table%ejected_mass(i,:))), &
               maxval(abs(table%net_yield(i,:))))
          if (abs(table%returned_mass(i)) > tol*scale .or. &
               abs(table%remnant_mass(i)) > tol*scale .or. &
               abs(table%energy(i)) > tol*scale .or. &
               maxval(abs(table%momentum(i,:))) > tol*scale .or. &
               maxval(abs(table%ejected_mass(i,:))) > tol*scale .or. &
               maxval(abs(table%net_yield(i,:))) > tol*scale) then
             ierr = ior(ierr, yield_audit_err_value)
          end if
       end if

       ejected_sum = sum(table%ejected_mass(i,:))
       scale = max(1.0_stellar_dp, abs(table%returned_mass(i)), &
            abs(ejected_sum))
       if (ejected_sum > table%returned_mass(i) + tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if

       scale = max(1.0_stellar_dp, table%initial_mass(i))
       if (table%returned_mass(i) + table%remnant_mass(i) > &
            table%initial_mass(i) + tol * scale) then
          ierr = ior(ierr, yield_audit_err_mass)
       end if
       if (present(terminal_remnant_owner) .and. &
            table%channel(i) >= 1 .and. table%channel(i) <= n_stellar_channels) then
          if (.not. terminal_remnant_owner(table%channel(i)) .and. &
               table%remnant_mass(i) > tol) then
             ierr = ior(ierr, yield_audit_err_remnant_ownership)
          end if
       end if
    end do

    ! A duplicate coordinate would make the interpolation result depend on
    ! row order.  Reject it even when the duplicate values happen to agree.
    do i = 1, table%n_rows
       do j = i + 1, table%n_rows
          if (table%channel(i) /= table%channel(j)) cycle
          if (.not. same_value(table%initial_mass(i), &
               table%initial_mass(j), coordinate_tol)) cycle
          if (.not. same_value(table%birth_metallicity(i), &
               table%birth_metallicity(j), coordinate_tol)) cycle
          if (.not. same_value(table%age_gyr(i), table%age_gyr(j), coordinate_tol)) cycle
          ierr = ior(ierr, yield_audit_err_duplicate)
       end do
    end do

    ! Compare rows on the same channel, mass, and metallicity grid line.
    ! Actual cumulative material, returned mass, and injected energy must not
    ! decrease with age.  Net yields and momentum are allowed to be signed.
    do i = 1, table%n_rows
       do j = 1, table%n_rows
          if (i == j) cycle
          if (table%channel(i) /= table%channel(j)) cycle
          if (.not. same_value(table%initial_mass(i), &
               table%initial_mass(j), coordinate_tol)) cycle
          if (.not. same_value(table%birth_metallicity(i), &
               table%birth_metallicity(j), coordinate_tol)) cycle
          if (table%age_gyr(i) >= table%age_gyr(j) - coordinate_tol) cycle

          scale = max(1.0_stellar_dp, abs(table%returned_mass(j)))
          if (table%returned_mass(i) > table%returned_mass(j) + tol * scale) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, maxval(abs(table%ejected_mass(j,:))))
          if (any(table%ejected_mass(i,:) > table%ejected_mass(j,:) + &
               tol * scale)) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, &
               abs(table%returned_mass(i) - sum(table%ejected_mass(i,:))), &
               abs(table%returned_mass(j) - sum(table%ejected_mass(j,:))))
          if (table%returned_mass(i) - sum(table%ejected_mass(i,:)) > &
               table%returned_mass(j) - sum(table%ejected_mass(j,:)) + &
               tol * scale) then
             ierr = ior(ierr, yield_audit_err_monotonic)
          end if
          scale = max(1.0_stellar_dp, abs(table%energy(j)))
          if (table%energy(i) > table%energy(j) + tol * scale) then
             ierr = ior(ierr, yield_audit_err_energy_monotonic)
          end if
       end do
    end do

    if (require_grid) then
       do channel = channel_wind, channel_snii
          if(present(required_channels))then
             if(.not.required_channels(channel))cycle
          endif
          ! A validated terminal AGB map is sparse by design: each M,Z has
          ! its own event age and mass grid. prepare_agb_terminal_rows already
          ! checks every row, the zero origin, endpoint closure and plateau.
          ! Do not manufacture rectangular interpolation corners for it.
          if(channel==channel_agb.and.allocated(table%agb_terminal_row))cycle
          ! The validated source-node evaluator likewise uses each wind's
          ! actual phase knots and the terminal event directly, never a
          ! Cartesian cross-star age interpolation. Its complete endpoint
          ! map and every history row were checked by prepare_high_mass_history.
          if(table%high_mass_ready.and.(channel==channel_wind.or.channel==channel_snii))cycle
          call audit_complete_channel(table, channel, coordinate_tol, channel_is_bad)
          if (channel_is_bad) ierr = ior(ierr, yield_audit_err_grid)
       end do
    end if
  end subroutine audit_yield_table

  subroutine audit_complete_channel(table, channel_id, tolerance, bad)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: tolerance
    logical, intent(out) :: bad

    real(stellar_dp), allocatable :: masses(:), metallicities(:), ages(:)
    integer :: i, j, k, status, n_rows, n_mass, n_z, n_age
    integer(kind=8) :: expected_rows

    bad = .false.
    n_rows = count(table%channel(1:table%n_rows) == channel_id)
    if (n_rows <= 0) then
       bad = .true.
       return
    end if

    allocate(masses(n_rows), metallicities(n_rows), ages(n_rows), stat=status)
    if (status /= 0) then
       bad = .true.
       return
    end if
    masses = 0.0_stellar_dp
    metallicities = 0.0_stellar_dp
    ages = 0.0_stellar_dp
    n_mass = 0
    n_z = 0
    n_age = 0
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       call append_unique(table%initial_mass(i), masses, n_mass, tolerance)
       call append_unique(table%birth_metallicity(i), metallicities, n_z, &
            tolerance)
       call append_unique(table%age_gyr(i), ages, n_age, tolerance)
    end do

    expected_rows = int(n_mass, kind=8) * int(n_z, kind=8) * int(n_age, kind=8)
    if (int(n_rows, kind=8) /= expected_rows) bad = .true.
    if (.not. any(abs(ages(1:n_age)) <= tolerance)) bad = .true.

    do i = 1, n_mass
       do j = 1, n_z
          do k = 1, n_age
             if (.not. grid_row_exists(table, channel_id, masses(i), &
                  metallicities(j), ages(k), tolerance)) bad = .true.
          end do
       end do
    end do
    deallocate(masses, metallicities, ages)
  end subroutine audit_complete_channel

  subroutine append_unique(value, values, n_values, tolerance)
    real(stellar_dp), intent(in) :: value, tolerance
    real(stellar_dp), intent(inout) :: values(:)
    integer, intent(inout) :: n_values
    integer :: i

    do i = 1, n_values
       if (same_value(values(i), value, tolerance)) return
    end do
    n_values = n_values + 1
    values(n_values) = value
  end subroutine append_unique

  logical function grid_row_exists(table, channel_id, mass, metallicity, age, &
       tolerance)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: mass, metallicity, age, tolerance
    integer :: i

    grid_row_exists = .false.
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       if (.not. same_value(table%initial_mass(i), mass, &
            tolerance)) cycle
       if (.not. same_value(table%birth_metallicity(i), metallicity, &
            tolerance)) cycle
       if (.not. same_value(table%age_gyr(i), age, tolerance)) cycle
       grid_row_exists = .true.
       return
    end do
  end function grid_row_exists

  logical function finite_row(table, row)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: row
    integer :: i

    finite_row = ieee_is_finite(table%initial_mass(row)) .and. &
         ieee_is_finite(table%birth_metallicity(row)) .and. &
         ieee_is_finite(table%age_gyr(row)) .and. &
         ieee_is_finite(table%returned_mass(row)) .and. &
         ieee_is_finite(table%remnant_mass(row)) .and. &
         ieee_is_finite(table%energy(row))
    do i = 1, 3
       finite_row = finite_row .and. ieee_is_finite(table%momentum(row,i))
    end do
    do i = 1, n_stellar_elements
       finite_row = finite_row .and. &
            ieee_is_finite(table%ejected_mass(row,i)) .and. &
            ieee_is_finite(table%net_yield(row,i))
    end do
  end function finite_row

  logical function same_value(a, b, tolerance)
    real(stellar_dp), intent(in) :: a, b, tolerance
    real(stellar_dp) :: scale

    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
       same_value = .false.
       return
    end if
    scale = max(1.0_stellar_dp, abs(a), abs(b))
    same_value = abs(a - b) <= tolerance * scale
  end function same_value

end module stellar_yield_audit

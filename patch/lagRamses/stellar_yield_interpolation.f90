! Phase 0 interpolation layer for the common stellar-yield table.
!
! Interpolation is performed only when all required grid corners exist.  A
! missing corner is reported to the caller rather than silently replaced by
! a nearest-neighbor or extrapolated value.  This is important for the AGB
! tables, where the metallicity and mass grids are not necessarily identical
! to the massive-star grids.

module stellar_yield_interpolation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels, channel_wind, channel_snii, channel_agb
  use stellar_yield_tables, only: stellar_yield_table_t, &
       yield_mass_assignment_linear, yield_mass_assignment_piecewise_constant
  implicit none

  private

  integer, parameter, public :: interpolation_ok = 0
  integer, parameter, public :: interpolation_err_table = 1
  integer, parameter, public :: interpolation_err_channel = 2
  integer, parameter, public :: interpolation_err_grid = 3
  integer, parameter, public :: interpolation_err_argument = 4
  integer, parameter, public :: interpolation_err_nonfinite = 5
  integer, parameter, public :: interpolation_err_assignment_mode = 6

  public :: interpolate_yield_row

contains

  subroutine interpolate_yield_row(table, channel_id, query_mass, query_z, &
       query_age_gyr, returned_mass, remnant_mass, energy, momentum, &
       ejected_mass, net_yield, ierr)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: query_mass, query_z, query_age_gyr
    real(stellar_dp), intent(out) :: returned_mass, remnant_mass, energy
    real(stellar_dp), intent(out) :: momentum(3)
    real(stellar_dp), intent(out) :: ejected_mass(n_stellar_elements)
    real(stellar_dp), intent(out) :: net_yield(n_stellar_elements)
    integer, intent(out) :: ierr

    real(stellar_dp) :: mass_nodes(2), z_nodes(2), age_nodes(2)
    real(stellar_dp) :: mass_weights(2), z_weights(2), age_weights(2)
    real(stellar_dp) :: mass_lo, mass_hi, z_lo, z_hi, age_lo, age_hi
    real(stellar_dp) :: weight
    real(stellar_dp) :: mass_value, z_value, age_value
    real(stellar_dp) :: eps
    integer :: n_mass_nodes, n_z_nodes, n_age_nodes
    integer :: im, iz, ia, row, first_row, last_row, age_count
    logical :: found, constant_age_payload

    ierr = interpolation_ok
    returned_mass = 0.0_stellar_dp
    remnant_mass = 0.0_stellar_dp
    energy = 0.0_stellar_dp
    momentum = 0.0_stellar_dp
    ejected_mass = 0.0_stellar_dp
    net_yield = 0.0_stellar_dp

    if (.not. table%loaded .or. table%n_rows <= 0) then
       ierr = interpolation_err_table
       return
    end if
    if (channel_id < 1 .or. channel_id > n_stellar_channels) then
       ierr = interpolation_err_channel
       return
    end if
    if (.not. ieee_is_finite(query_mass) .or. &
         .not. ieee_is_finite(query_z) .or. &
         .not. ieee_is_finite(query_age_gyr)) then
       ierr = interpolation_err_argument
       return
    end if
    if (query_mass <= 0.0_stellar_dp .or. query_z < 0.0_stellar_dp .or. &
         query_age_gyr < 0.0_stellar_dp) then
       ierr = interpolation_err_argument
       return
    end if
    if (table%mass_assignment_mode /= yield_mass_assignment_linear .and. &
         table%mass_assignment_mode /= yield_mass_assignment_piecewise_constant) then
       ierr = interpolation_err_assignment_mode
       return
    end if

    if(table%high_mass_ready.and.allocated(table%agb_terminal_row).and.channel_id==channel_agb)then
       call agb_terminal_value(table,query_mass,query_z,query_age_gyr, &
            returned_mass,remnant_mass,energy,momentum,ejected_mass,net_yield,ierr)
       return
    endif
    if(table%high_mass_ready.and.query_mass>=40.and. &
         (channel_id==channel_wind.or.channel_id==channel_snii))then
       call high_mass_history_value(table,channel_id,query_mass,query_z,query_age_gyr, &
            returned_mass,remnant_mass,energy,momentum,ejected_mass,net_yield,ierr)
       return
    endif

    call find_bounds(table, channel_id, 1, query_mass, mass_lo, mass_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if
    call find_bounds(table, channel_id, 2, query_z, z_lo, z_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if
    call find_bounds(table, channel_id, 3, query_age_gyr, age_lo, age_hi, found)
    if (.not. found) then
       ierr = interpolation_err_grid
       return
    end if

    if (table%mass_assignment_mode == yield_mass_assignment_piecewise_constant) then
       call make_piecewise_mass_node(mass_lo, mass_hi, query_mass, mass_nodes, &
            mass_weights, n_mass_nodes)
    else
       call make_nodes(mass_lo, mass_hi, query_mass, mass_nodes, mass_weights, &
            n_mass_nodes)
    end if
    call make_nodes(z_lo, z_hi, query_z, z_nodes, z_weights, n_z_nodes)
    call make_nodes(age_lo, age_hi, query_age_gyr, age_nodes, age_weights, &
         n_age_nodes)

    eps = 1.0e-10_stellar_dp
    do im = 1, n_mass_nodes
       do iz = 1, n_z_nodes
          age_count=n_age_nodes;constant_age_payload=.false.
          if(table%high_mass_ready.and.n_age_nodes==2)then
             first_row=find_grid_row(table,channel_id,mass_nodes(im),z_nodes(iz),age_nodes(1))
             last_row=find_grid_row(table,channel_id,mass_nodes(im),z_nodes(iz),age_nodes(2))
             if(first_row>0.and.last_row>0)then
                constant_age_payload=table%returned_mass(first_row)==table%returned_mass(last_row).and. &
                     table%remnant_mass(first_row)==table%remnant_mass(last_row).and. &
                     table%energy(first_row)==table%energy(last_row).and. &
                     all(table%momentum(first_row,:)==table%momentum(last_row,:)).and. &
                     all(table%ejected_mass(first_row,:)==table%ejected_mass(last_row,:)).and. &
                     all(table%net_yield(first_row,:)==table%net_yield(last_row,:))
                if(constant_age_payload)age_count=1
             endif
          endif
          do ia = 1, age_count
             mass_value = mass_nodes(im)
             z_value = z_nodes(iz)
             age_value = age_nodes(ia)
             row = find_grid_row(table, channel_id, mass_value, z_value, &
                  age_value)
             if (row < 1) then
                ierr = interpolation_err_grid
                return
             end if

             weight = mass_weights(im) * z_weights(iz) * age_weights(ia)
             if(constant_age_payload)weight=mass_weights(im)*z_weights(iz)
             ! The selected history model uses fractions per initial mass in
             ! source cells, including the low-mass AGB WD supplier. Preserve
             ! old table semantics outside this explicitly selected route.
             if(table%high_mass_ready.and.table%mass_assignment_mode==yield_mass_assignment_piecewise_constant) &
                  weight=weight*query_mass/mass_value
             if (abs(weight) <= eps) cycle
             returned_mass = returned_mass + weight * table%returned_mass(row)
             remnant_mass = remnant_mass + weight * table%remnant_mass(row)
             energy = energy + weight * table%energy(row)
             momentum = momentum + weight * table%momentum(row,:)
             ejected_mass = ejected_mass + weight * &
                  table%ejected_mass(row,:)
             net_yield = net_yield + weight * table%net_yield(row,:)
          end do
       end do
    end do

    if (.not. finite_result(returned_mass, remnant_mass, energy, momentum, &
         ejected_mass, net_yield)) then
       ierr = interpolation_err_nonfinite
       returned_mass = 0.0_stellar_dp
       remnant_mass = 0.0_stellar_dp
       energy = 0.0_stellar_dp
       momentum = 0.0_stellar_dp
       ejected_mass = 0.0_stellar_dp
       net_yield = 0.0_stellar_dp
    end if
  end subroutine interpolate_yield_row

  subroutine agb_terminal_value(table,mass,z,age,returned,remnant,energy,p,elements,net,ierr)
    type(stellar_yield_table_t),intent(in)::table
    real(stellar_dp),intent(in)::mass,z,age
    real(stellar_dp),intent(out)::returned,remnant,energy,p(3),elements(n_stellar_elements),net(n_stellar_elements)
    integer,intent(out)::ierr
    real(stellar_dp)::zl,zh,zs(2),weights(2),best,distance,factor,ml,mh
    integer::i,j,r,node,nz
    returned=0;remnant=0;energy=0;p=0;elements=0;net=0
    ierr=interpolation_err_grid;zl=-huge(1d0);zh=huge(1d0)
    do i=1,size(table%agb_terminal_row)
       r=table%agb_terminal_row(i)
       if(abs(z-table%birth_metallicity(r))<=32*epsilon(1d0)* &
            max(abs(z),abs(table%birth_metallicity(r)),tiny(1d0)))then
          zl=table%birth_metallicity(r);zh=zl;exit
       endif
       if(table%birth_metallicity(r)<z)zl=max(zl,table%birth_metallicity(r))
       if(table%birth_metallicity(r)>z)zh=min(zh,table%birth_metallicity(r))
    enddo
    if(zl<0.or.zh==huge(1d0))return
    nz=1;zs=[zl,zh];weights=[1d0,0d0]
    if(zl/=zh)then
       if(.not.table%high_mass_linear_z)return
       nz=2;weights(2)=(z-zl)/(zh-zl);weights(1)=1-weights(2)
    endif
    do j=1,nz
       best=huge(1d0);node=0;ml=huge(1d0);mh=0
       do i=1,size(table%agb_terminal_row)
          r=table%agb_terminal_row(i)
          if(table%birth_metallicity(r)/=zs(j))cycle
          ml=min(ml,table%initial_mass(r));mh=max(mh,table%initial_mass(r))
          distance=abs(mass-table%initial_mass(r))
          if(distance<best)then
             best=distance;node=r
          else if(distance==best.and.node>0)then
             if(table%initial_mass(r)<table%initial_mass(node))node=r
          endif
       enddo
       if(node==0.or.mass<ml.or.mass>mh)then
          returned=0;remnant=0;energy=0;p=0;elements=0;net=0
          return
       endif
       if(age<table%age_gyr(node))cycle
       ! Same-age mixture, nearest source-mass cells, fractional budgets.
       ! No early envelope return, early WD creation, or interpolated lifetime.
       factor=weights(j)*mass/table%initial_mass(node)
       returned=returned+factor*table%returned_mass(node)
       remnant=remnant+factor*table%remnant_mass(node)
       energy=energy+factor*table%energy(node)
       p=p+factor*table%momentum(node,:)
       elements=elements+factor*table%ejected_mass(node,:)
       net=net+factor*table%net_yield(node,:)
    enddo
    ierr=interpolation_ok
  end subroutine agb_terminal_value

  subroutine high_mass_history_value(table,channel,mass,z,age,returned,remnant,energy,p,elements,net,ierr)
    type(stellar_yield_table_t),intent(in)::table
    integer,intent(in)::channel
    real(stellar_dp),intent(in)::mass,z,age
    real(stellar_dp),intent(out)::returned,remnant,energy,p(3),elements(n_stellar_elements),net(n_stellar_elements)
    integer,intent(out)::ierr
    real(stellar_dp)::zl,zh,w,a(6+2*n_stellar_elements),b(6+2*n_stellar_elements),v(6+2*n_stellar_elements)
    integer::i,status
    returned=0;remnant=0;energy=0;p=0;elements=0;net=0
    if(.not.table%high_mass_linear_z)then
       call high_mass_node_value(table,channel,mass,z,age,returned,remnant,energy,p,elements,net,ierr)
       return
    endif
    ierr=interpolation_err_grid
    zl=-huge(1d0);zh=huge(1d0)
    do i=1,size(table%hm_z)
       if(abs(z-table%hm_z(i))<=32*epsilon(1d0)*max(abs(z),abs(table%hm_z(i)),tiny(1d0)))then
          call high_mass_node_value(table,channel,mass,table%hm_z(i),age,returned,remnant,energy,p,elements,net,ierr)
          return
       endif
       if(table%hm_z(i)<z)zl=max(zl,table%hm_z(i))
       if(table%hm_z(i)>z)zh=min(zh,table%hm_z(i))
    enddo
    if(zl<0.or.zh==huge(1d0))return ! Never extrapolate or clamp physical Z.
    call high_mass_node_value(table,channel,mass,zl,age,a(1),a(2),a(3),a(4:6), &
         a(7:6+n_stellar_elements),a(7+n_stellar_elements:),status)
    if(status/=0)return
    call high_mass_node_value(table,channel,mass,zh,age,b(1),b(2),b(3),b(4:6), &
         b(7:6+n_stellar_elements),b(7+n_stellar_elements:),status)
    if(status/=0)return
    ! Interpolate already resolved cumulative budgets at the SAME age. Fixed
    ! birth-Z weights preserve monotonic return and timestep telescoping even
    ! across different lifetimes/outcomes. This is an SSP expectation mixture,
    ! not interpolation of an individual star's fate or a new explosion time.
    w=(z-zl)/(zh-zl);v=(1-w)*a+w*b
    returned=v(1);remnant=v(2);energy=v(3);p=v(4:6)
    elements=v(7:6+n_stellar_elements);net=v(7+n_stellar_elements:)
    ierr=interpolation_ok
  end subroutine high_mass_history_value

  subroutine high_mass_node_value(table,channel,mass,z,age,returned,remnant,energy,p,elements,net,ierr)
    type(stellar_yield_table_t),intent(in)::table
    integer,intent(in)::channel
    real(stellar_dp),intent(in)::mass,z,age
    real(stellar_dp),intent(out)::returned,remnant,energy,p(3),elements(n_stellar_elements),net(n_stellar_elements)
    integer,intent(out)::ierr
    integer::i,node,lo,hi,row
    real(stellar_dp)::distance,best,tq,fraction,factor
    returned=0;remnant=0;energy=0;p=0;elements=0;net=0
    ierr=interpolation_err_grid
    if(mass<40.or.mass>120)return
    ! No interpolation of discrete fates across Z, rotation or source engines.
    ! Mass cells use nearest-node midpoints, with ties assigned to the lower
    ! node. Budgets scale by M_query/M_node, keeping mass fractions conserved.
    best=huge(1d0);node=0
    do i=1,size(table%hm_mass)
       ! Birth-metallicity division/advection can differ by a few ULPs.
       ! This is relative to Z itself (not max(1,Z)), so primordial and
       ! ultra-low-Z branches cannot acquire a spurious absolute tolerance.
       if(abs(z-table%hm_z(i))>32*epsilon(1d0)*max(abs(z),abs(table%hm_z(i)),tiny(1d0)))cycle
       distance=abs(mass-table%hm_mass(i))
       if(distance<best)then
          best=distance;node=i
       endif
    enddo
    if(node==0)return
    factor=mass/table%hm_mass(node)
    if(channel==channel_snii)then
       ierr=interpolation_ok
       if(age<table%hm_age(node))return
       remnant=factor*table%hm_remnant(node)
       if(table%high_mass_wind_only)return
       row=table%hm_terminal_row(node)
       returned=factor*table%returned_mass(row); energy=factor*table%energy(row)
       p=factor*table%momentum(row,:);elements=factor*table%ejected_mass(row,:);net=factor*table%net_yield(row,:)
       return
    endif
    tq=min(age,table%hm_age(node));lo=0;hi=0
    do i=1,table%n_rows
       if(table%channel(i)/=channel_wind.or.table%initial_mass(i)/=table%hm_mass(node).or. &
            table%birth_metallicity(i)/=table%hm_z(node))cycle
       if(table%age_gyr(i)<=tq)then
          if(lo==0)then
             lo=i
          else if(table%age_gyr(i)>table%age_gyr(lo))then
             lo=i
          endif
       endif
       if(table%age_gyr(i)>=tq)then
          if(hi==0)then
             hi=i
          else if(table%age_gyr(i)<table%age_gyr(hi))then
             hi=i
          endif
       endif
    enddo
    if(lo==0.or.hi==0)return
    fraction=0
    if(lo/=hi)fraction=(tq-table%age_gyr(lo))/(table%age_gyr(hi)-table%age_gyr(lo))
    returned=factor*((1-fraction)*table%returned_mass(lo)+fraction*table%returned_mass(hi))
    energy=factor*((1-fraction)*table%energy(lo)+fraction*table%energy(hi))
    p=factor*((1-fraction)*table%momentum(lo,:)+fraction*table%momentum(hi,:))
    elements=factor*((1-fraction)*table%ejected_mass(lo,:)+fraction*table%ejected_mass(hi,:))
    net=factor*((1-fraction)*table%net_yield(lo,:)+fraction*table%net_yield(hi,:))
    ierr=interpolation_ok
  end subroutine high_mass_node_value

  subroutine find_bounds(table, channel_id, axis, query, lower, upper, found)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id, axis
    real(stellar_dp), intent(in) :: query
    real(stellar_dp), intent(out) :: lower, upper
    logical, intent(out) :: found

    integer :: i
    real(stellar_dp) :: value
    logical :: has_lower, has_upper

    lower = 0.0_stellar_dp
    upper = 0.0_stellar_dp
    has_lower = .false.
    has_upper = .false.
    found = .false.

    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       value = coordinate_value(table, i, axis)
       found = .true.
       if (value <= query) then
          if (.not. has_lower .or. value > lower) lower = value
          has_lower = .true.
       end if
       if (value >= query) then
          if (.not. has_upper .or. value < upper) upper = value
          has_upper = .true.
       end if
    end do

    if (.not. found .or. .not. has_lower .or. .not. has_upper) then
       found = .false.
       lower = 0.0_stellar_dp
       upper = 0.0_stellar_dp
    end if
  end subroutine find_bounds

  subroutine make_nodes(lower, upper, query, nodes, weights, n_nodes)
    real(stellar_dp), intent(in) :: lower, upper, query
    real(stellar_dp), intent(out) :: nodes(2), weights(2)
    integer, intent(out) :: n_nodes
    real(stellar_dp) :: fraction
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp

    nodes = 0.0_stellar_dp
    weights = 0.0_stellar_dp
    if (abs(upper - lower) <= tolerance * max(1.0_stellar_dp, &
         abs(lower), abs(upper))) then
       n_nodes = 1
       nodes(1) = lower
       weights(1) = 1.0_stellar_dp
       return
    end if

    n_nodes = 2
    nodes(1) = lower
    nodes(2) = upper
    fraction = (query - lower) / (upper - lower)
    fraction = max(0.0_stellar_dp, min(1.0_stellar_dp, fraction))
    weights(1) = 1.0_stellar_dp - fraction
    weights(2) = fraction
  end subroutine make_nodes

  subroutine make_piecewise_mass_node(lower, upper, query, nodes, weights, n_nodes)
    ! A source-node fate is not a quantity that can be linearly blended with
    ! the neighboring node.  Select the left node for an interior half-open
    ! cell and the exact node at a grid edge.  Z and age remain handled by
    ! their ordinary table policy; source-node callers must use exact values
    ! on those axes when their fate semantics are discrete.
    real(stellar_dp), intent(in) :: lower, upper, query
    real(stellar_dp), intent(out) :: nodes(2), weights(2)
    integer, intent(out) :: n_nodes
    real(stellar_dp), parameter :: tolerance = 1.0e-12_stellar_dp

    nodes = 0.0_stellar_dp
    weights = 0.0_stellar_dp
    n_nodes = 1
    if (abs(upper - lower) <= tolerance * max(1.0_stellar_dp, &
         abs(lower), abs(upper))) then
       nodes(1) = lower
    else if (query >= upper) then
       nodes(1) = upper
    else
       nodes(1) = lower
    endif
    weights(1) = 1.0_stellar_dp
  end subroutine make_piecewise_mass_node

  integer function find_grid_row(table, channel_id, mass, metallicity, age)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: channel_id
    real(stellar_dp), intent(in) :: mass, metallicity, age
    integer :: i

    find_grid_row = 0
    do i = 1, table%n_rows
       if (table%channel(i) /= channel_id) cycle
       if (.not. same_coordinate(table%initial_mass(i), mass)) cycle
       if (.not. same_coordinate(table%birth_metallicity(i), metallicity)) cycle
       if (.not. same_coordinate(table%age_gyr(i), age)) cycle
       find_grid_row = i
       return
    end do
  end function find_grid_row

  real(stellar_dp) function coordinate_value(table, row, axis)
    type(stellar_yield_table_t), intent(in) :: table
    integer, intent(in) :: row, axis

    select case (axis)
    case (1)
       coordinate_value = table%initial_mass(row)
    case (2)
       coordinate_value = table%birth_metallicity(row)
    case (3)
       coordinate_value = table%age_gyr(row)
    case default
       coordinate_value = 0.0_stellar_dp
    end select
  end function coordinate_value

  logical function same_coordinate(a, b)
    real(stellar_dp), intent(in) :: a, b
    real(stellar_dp) :: scale

    scale = max(1.0_stellar_dp, abs(a), abs(b))
    same_coordinate = abs(a - b) <= 1.0e-10_stellar_dp * scale
  end function same_coordinate

  logical function finite_result(returned_mass, remnant_mass, energy, momentum, &
       ejected_mass, net_yield)
    real(stellar_dp), intent(in) :: returned_mass, remnant_mass, energy
    real(stellar_dp), intent(in) :: momentum(3)
    real(stellar_dp), intent(in) :: ejected_mass(n_stellar_elements)
    real(stellar_dp), intent(in) :: net_yield(n_stellar_elements)
    integer :: i

    finite_result = ieee_is_finite(returned_mass) .and. &
         ieee_is_finite(remnant_mass) .and. ieee_is_finite(energy)
    do i = 1, 3
       finite_result = finite_result .and. ieee_is_finite(momentum(i))
    end do
    do i = 1, n_stellar_elements
       finite_result = finite_result .and. &
            ieee_is_finite(ejected_mass(i)) .and. ieee_is_finite(net_yield(i))
    end do
  end function finite_result

end module stellar_yield_interpolation

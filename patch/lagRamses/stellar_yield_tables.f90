! Phase 0 common stellar-yield table reader.
!
! The table is an ASCII, one-row-per-channel/per-initial-star table.  Values
! are read but not interpolated here; interpolation and IMF integration belong
! to the source engine.  This keeps file I/O, table provenance, and physics
! evaluation separate.
!
! Data-row format (free format, after optional '#' comment lines):
!
! channel initial_mass birth_metallicity age_yr returned_mass remnant_mass energy
! momentum_x momentum_y momentum_z ejecta[H..Fe] net_yield[H..Fe]
!
! Units:
!   mass values       : Msun per initial star
!   age_yr            : yr on disk; converted once to age_gyr in memory
!   energy            : erg per initial star
!   momentum          : g cm/s per initial star
!   elemental arrays  : Msun per initial star

module stellar_yield_tables
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dust_mass_physics, only: dust_species_condense
  use stellar_enrichment_config, only: stellar_dp, n_stellar_elements, &
       n_stellar_channels
  implicit none

  private

  integer, parameter, public :: yield_table_ok = 0
  integer, parameter, public :: yield_table_err_open = 1
  integer, parameter, public :: yield_table_err_read = 2
  integer, parameter, public :: yield_table_err_empty = 3
  integer, parameter, public :: yield_table_err_alloc = 4
  integer, parameter, public :: yield_table_err_format = 5
  integer, parameter, public :: yield_table_err_nonfinite = 6
  integer, parameter, public :: yield_table_err_assignment_mode = 7
  integer, parameter, public :: yield_mass_assignment_linear = 0
  integer, parameter, public :: yield_mass_assignment_piecewise_constant = 1

  type, public :: stellar_yield_table_t
     logical :: loaded = .false.
     integer :: n_rows = 0
     integer :: mass_assignment_mode = yield_mass_assignment_linear
     ! Per-source-node high-mass adapter, admitted only after full validation.
     logical :: high_mass_ready = .false.
     logical :: high_mass_wind_only = .false.
     logical :: high_mass_linear_z = .false.
     logical :: net_yield_diagnostic_unavailable = .false.
     logical :: net_yield_channel_available(n_stellar_channels) = .true.
     character(len=128) :: high_mass_identity(4) = ''
     real(stellar_dp), allocatable :: hm_mass(:), hm_z(:), hm_age(:), hm_remnant(:), hm_adjustment(:)
     integer, allocatable :: hm_wind_row(:), hm_terminal_row(:)
     ! Optional single terminal envelope/WD event, indexed by source M,Z.
     integer, allocatable :: agb_terminal_row(:)
     ! CO=0, hybrid CO(Ne)=1, ONe=2. Only CO can supply the strict Ia ledger.
     integer, allocatable :: agb_remnant_kind(:)
     ! Fraction of the final cumulative envelope released as a terminal jump.
     ! The rest follows source wind knots; remnant creation is always a step.
     real(stellar_dp), allocatable :: agb_terminal_jump_fraction(:)
     logical :: co_wd_inventory_only = .false. ! Internal read-only inventory view.
     integer, allocatable :: channel(:)
     real(stellar_dp), allocatable :: initial_mass(:)
     real(stellar_dp), allocatable :: birth_metallicity(:)
     real(stellar_dp), allocatable :: age_gyr(:)
     real(stellar_dp), allocatable :: returned_mass(:)
     real(stellar_dp), allocatable :: remnant_mass(:)
     real(stellar_dp), allocatable :: energy(:)
     real(stellar_dp), allocatable :: momentum(:,:)
     real(stellar_dp), allocatable :: ejected_mass(:,:)
     real(stellar_dp), allocatable :: net_yield(:,:)
     real(stellar_dp), allocatable :: dust_ejected(:,:),dust_jump(:,:)
  end type stellar_yield_table_t

  public :: clear_yield_table
  public :: load_yield_table
  public :: set_yield_mass_assignment_mode
  public :: prepare_dust_yields

contains

  subroutine clear_yield_table(table)
    type(stellar_yield_table_t), intent(inout) :: table

    if (allocated(table%channel)) deallocate(table%channel)
    if (allocated(table%initial_mass)) deallocate(table%initial_mass)
    if (allocated(table%birth_metallicity)) deallocate(table%birth_metallicity)
    if (allocated(table%age_gyr)) deallocate(table%age_gyr)
    if (allocated(table%returned_mass)) deallocate(table%returned_mass)
    if (allocated(table%remnant_mass)) deallocate(table%remnant_mass)
    if (allocated(table%energy)) deallocate(table%energy)
    if (allocated(table%momentum)) deallocate(table%momentum)
    if (allocated(table%ejected_mass)) deallocate(table%ejected_mass)
    if (allocated(table%net_yield)) deallocate(table%net_yield)
    if (allocated(table%dust_ejected)) deallocate(table%dust_ejected,table%dust_jump)
    if (allocated(table%hm_mass)) deallocate(table%hm_mass, table%hm_z, table%hm_age, &
         table%hm_remnant, table%hm_adjustment, table%hm_wind_row, table%hm_terminal_row)
    if (allocated(table%agb_terminal_row)) deallocate(table%agb_terminal_row)
    if (allocated(table%agb_remnant_kind)) deallocate(table%agb_remnant_kind)
    if (allocated(table%agb_terminal_jump_fraction)) deallocate(table%agb_terminal_jump_fraction)
    table%co_wd_inventory_only = .false.
    table%high_mass_ready = .false.
    table%high_mass_wind_only = .false.
    table%high_mass_linear_z = .false.
    table%net_yield_diagnostic_unavailable = .false.
    table%net_yield_channel_available = .true.
    table%high_mass_identity = ''

    table%loaded = .false.
    table%n_rows = 0
    table%mass_assignment_mode = yield_mass_assignment_linear
  end subroutine clear_yield_table

  subroutine prepare_dust_yields(table,ierr)
    ! Integrate condensation over physical release segments, BEFORE the
    ! existing age/Z/IMF mixture. In particular no SSP-averaged C/O switch.
    type(stellar_yield_table_t),intent(inout)::table
    integer,intent(out)::ierr
    integer::r,j,k,previous,selected,pass,status
    logical,allocatable::done(:)
    real(stellar_dp)::age,delta(n_stellar_elements),jump(n_stellar_elements),part(2),jp(2),before(2),tol
    ierr=1
    if(.not.table%loaded.or.n_stellar_elements/=11)return
    if(allocated(table%dust_ejected))deallocate(table%dust_ejected,table%dust_jump)
    allocate(table%dust_ejected(table%n_rows,2),table%dust_jump(table%n_rows,2),done(table%n_rows))
    table%dust_ejected=0;table%dust_jump=0;done=.false.
    do pass=1,table%n_rows
       selected=0;age=huge(1d0)
       do r=1,table%n_rows
          if(done(r))cycle
          if(table%age_gyr(r)<age)then
             age=table%age_gyr(r);selected=r
          endif
       enddo
       if(selected==0)return
       r=selected;done(r)=.true.
       if(table%channel(r)>3)cycle
       previous=0
       do j=1,table%n_rows
          if(table%channel(j)/=table%channel(r).or.table%initial_mass(j)/=table%initial_mass(r).or. &
               table%birth_metallicity(j)/=table%birth_metallicity(r).or.table%age_gyr(j)>=age)cycle
          if(previous==0)then
             previous=j
          else if(table%age_gyr(j)>table%age_gyr(previous))then
             previous=j
          endif
       enddo
       delta=table%ejected_mass(r,:);before=0;jump=0
       if(previous>0)then
          delta=delta-table%ejected_mass(previous,:);before=table%dust_ejected(previous,:)
       endif
       if(allocated(table%agb_terminal_jump_fraction))then
          do k=1,size(table%agb_terminal_row)
             if(table%agb_terminal_row(k)==r)jump=table%agb_terminal_jump_fraction(k)*table%ejected_mass(r,:)
          enddo
       endif
       delta=delta-jump
       tol=64*epsilon(1d0)*max(table%returned_mass(r),tiny(1d0))
       if(any(delta < -tol))return
       call dust_species_condense(max(delta,0d0),table%channel(r),part,status)
       if(status/=0)return
       call dust_species_condense(jump,table%channel(r),jp,status)
       if(status/=0)return
       table%dust_ejected(r,:)=before+part+jp;table%dust_jump(r,:)=jp
       if(sum(table%dust_ejected(r,:))>table%returned_mass(r)- &
            sum(table%ejected_mass(r,1:2))+tol)return
    enddo
    ierr=0
  end subroutine prepare_dust_yields

  subroutine set_yield_mass_assignment_mode(table, mode, ierr)
    type(stellar_yield_table_t), intent(inout) :: table
    integer, intent(in) :: mode
    integer, intent(out) :: ierr

    ierr = yield_table_ok
    if (mode /= yield_mass_assignment_linear .and. &
         mode /= yield_mass_assignment_piecewise_constant) then
       ierr = yield_table_err_assignment_mode
       return
    end if
    table%mass_assignment_mode = mode
  end subroutine set_yield_mass_assignment_mode

  subroutine load_yield_table(filename, table, ierr)
    character(len=*), intent(in) :: filename
    type(stellar_yield_table_t), intent(inout) :: table
    integer, intent(out) :: ierr

    integer :: unit, ios, row, n_rows
    integer :: channel_id, ios_row
    real(stellar_dp) :: initial_mass, birth_metallicity, age_yr, age_gyr
    real(stellar_dp) :: returned_mass, remnant_mass, energy
    real(stellar_dp) :: momentum(3)
    real(stellar_dp) :: ejected_mass(n_stellar_elements)
    real(stellar_dp) :: net_yield(n_stellar_elements)
    character(len=4096) :: line

    ierr = yield_table_ok
    call clear_yield_table(table)

    open(newunit=unit, file=filename, status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
       ierr = yield_table_err_open
       return
    end if

    ! First pass counts valid data rows.  Comments and blank lines are ignored.
    n_rows = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios < 0) exit
       if (ios /= 0) then
          close(unit)
          ierr = yield_table_err_read
          return
       end if
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle

       read(line, *, iostat=ios_row) channel_id, initial_mass, &
            birth_metallicity, age_yr
       if (ios_row == 0) n_rows = n_rows + 1
    end do

    if (n_rows == 0) then
       close(unit)
       ierr = yield_table_err_empty
       return
    end if

    allocate(table%channel(n_rows), table%initial_mass(n_rows), &
         table%birth_metallicity(n_rows), table%age_gyr(n_rows), &
         table%returned_mass(n_rows), table%remnant_mass(n_rows), &
         table%energy(n_rows), table%momentum(n_rows,3), &
         table%ejected_mass(n_rows,n_stellar_elements), &
         table%net_yield(n_rows,n_stellar_elements), stat=ios)
    if (ios /= 0) then
       close(unit)
       call clear_yield_table(table)
       ierr = yield_table_err_alloc
       return
    end if

    rewind(unit, iostat=ios)
    if (ios /= 0) then
       close(unit)
       call clear_yield_table(table)
       ierr = yield_table_err_read
       return
    end if

    row = 0
    do
       read(unit, '(A)', iostat=ios) line
       if (ios < 0) exit
       if (ios /= 0) then
          close(unit)
          call clear_yield_table(table)
          ierr = yield_table_err_read
          return
       end if
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle

       read(line, *, iostat=ios_row) channel_id, initial_mass, &
            birth_metallicity, age_yr, returned_mass, remnant_mass, energy, &
            momentum(1), momentum(2), momentum(3), ejected_mass, net_yield
       if (ios_row /= 0) then
          ! A non-data line is allowed only if it was also ignored in pass 1.
          read(line, *, iostat=ios_row) channel_id, initial_mass, &
               birth_metallicity, age_yr
          if (ios_row == 0) then
             close(unit)
             call clear_yield_table(table)
             ierr = yield_table_err_format
             return
          end if
          cycle
       end if

       if (.not. finite_row_values(initial_mass, birth_metallicity, age_yr, &
            returned_mass, remnant_mass, energy, momentum, ejected_mass, &
            net_yield)) then
          close(unit)
          call clear_yield_table(table)
          ierr = yield_table_err_nonfinite
          return
       end if

       if (channel_id < 1 .or. channel_id > n_stellar_channels .or. &
            initial_mass <= 0.0_stellar_dp .or. birth_metallicity < 0.0_stellar_dp .or. &
            age_yr < 0.0_stellar_dp .or. returned_mass < 0.0_stellar_dp .or. &
            remnant_mass < 0.0_stellar_dp .or. energy < 0.0_stellar_dp .or. &
            minval(ejected_mass) < 0.0_stellar_dp) then
          close(unit)
          call clear_yield_table(table)
          ierr = yield_table_err_format
          return
       end if

       age_gyr = age_yr * 1.0e-9_stellar_dp
       if (.not. ieee_is_finite(age_gyr)) then
          close(unit)
          call clear_yield_table(table)
          ierr = yield_table_err_nonfinite
          return
       end if

       row = row + 1
       table%channel(row) = channel_id
       table%initial_mass(row) = initial_mass
       table%birth_metallicity(row) = birth_metallicity
       table%age_gyr(row) = age_gyr
       table%returned_mass(row) = returned_mass
       table%remnant_mass(row) = remnant_mass
       table%energy(row) = energy
       table%momentum(row,:) = momentum
       table%ejected_mass(row,:) = ejected_mass
       table%net_yield(row,:) = net_yield
    end do

    close(unit)

    if (row /= n_rows) then
       call clear_yield_table(table)
       ierr = yield_table_err_read
       return
    end if

    table%n_rows = row
    table%loaded = .true.
  end subroutine load_yield_table

  logical function finite_row_values(initial_mass, birth_metallicity, age_yr, &
       returned_mass, remnant_mass, energy, momentum, ejected_mass, net_yield)
    real(stellar_dp), intent(in) :: initial_mass, birth_metallicity, age_yr
    real(stellar_dp), intent(in) :: returned_mass, remnant_mass, energy
    real(stellar_dp), intent(in) :: momentum(3)
    real(stellar_dp), intent(in) :: ejected_mass(n_stellar_elements)
    real(stellar_dp), intent(in) :: net_yield(n_stellar_elements)
    integer :: i

    finite_row_values = ieee_is_finite(initial_mass) .and. &
         ieee_is_finite(birth_metallicity) .and. ieee_is_finite(age_yr) .and. &
         ieee_is_finite(returned_mass) .and. ieee_is_finite(remnant_mass) .and. &
         ieee_is_finite(energy)
    do i = 1, 3
       finite_row_values = finite_row_values .and. ieee_is_finite(momentum(i))
    end do
    do i = 1, n_stellar_elements
       finite_row_values = finite_row_values .and. &
            ieee_is_finite(ejected_mass(i)) .and. ieee_is_finite(net_yield(i))
    end do
  end function finite_row_values

end module stellar_yield_tables

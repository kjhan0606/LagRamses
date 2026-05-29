! Patch changes:
! - added mp0 variable for particles
module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)::msink,r2sink,v2sink,c2sink,oksink_new,oksink_all,tsink
  real(dp),allocatable,dimension(:)::msink_new,msink_all,r2k,v2sink_new,c2sink_new,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)::v2sink_all,c2sink_all
  real(dp),allocatable,dimension(:)::dMBHoverdt,dMEdoverdt,wdens,wvol,wc2
  real(dp),allocatable,dimension(:)::wdens_new,wvol_new,wc2_new,total_volume
  real(dp),allocatable,dimension(:,:)::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_c2
  real(dp),allocatable,dimension(:,:)::jsink,jsink_new,jsink_all
  real(dp),allocatable,dimension(:)::dMBH_coarse,dMEd_coarse,dMsmbh,dMBH_coarse_new
  real(dp),allocatable,dimension(:)::dMEd_coarse_new,dMsmbh_new,dMBH_coarse_all,dMEd_coarse_all,dMsmbh_all
  real(dp),allocatable,dimension(:)::Esave,Esave_new,Esave_all
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:,:,:)::sink_stat,sink_stat_all
  real(dp),allocatable,dimension(:)::c_avgptr,v_avgptr,d_avgptr
  real(dp),allocatable,dimension(:)::spinmag,spinmag_new,spinmag_all
  real(dp),allocatable,dimension(:,:)::bhspin,bhspin_new,bhspin_all
  real(dp),allocatable,dimension(:)::eps_sink
  integer ,allocatable,dimension(:)::idsink,idsink_new,idsink_all
  integer::nindsink=0

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch (pure data; type in ptypep)
  real(dp),allocatable,dimension(:,:)::weightp  ! weight of cloud parts for sink accretion only
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  real(dp),allocatable,dimension(:)  ::edp      ! Dark internal energy (aDM)
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  integer(kind=1),allocatable,dimension(:)::ptypep  ! Particle type code (see PTYPE_* below)

  ! Particle type codes — authoritative species tag for each particle
  integer(kind=1),parameter::PTYPE_DM        =  0_1  ! cold DM (ground state)
  integer(kind=1),parameter::PTYPE_STAR      =  1_1  ! stellar particle
  integer(kind=1),parameter::PTYPE_SINK      =  2_1  ! BH sink cloud particle (idp<0)
  integer(kind=1),parameter::PTYPE_ISIDM_EX1 = 10_1  ! iSIDM excited state 1
  integer(kind=1),parameter::PTYPE_ISIDM_EX2 = 11_1  ! iSIDM excited state 2 (multi-state)
  integer(kind=1),parameter::PTYPE_ADM       = 20_1  ! atomic DM

  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1

  !for chemo components
  real(dp),allocatable,dimension(:)  ::tpp, mp0, indtab
!  real(dp),allocatable,dimension(:,:)::cep
  character(len=3),allocatable,dimension(:) ::elem_list

  type yield_table
     integer::na
     integer::nz
     real,dimension(:)      ,pointer::astar
     real,dimension(:)      ,pointer::zstar
     real,dimension(:,:,:)  ,pointer::Eeject
     real,dimension(:,:)    ,pointer::Zeject
     real,dimension(:,:)    ,pointer::Meject
     real,dimension(:,:)    ,pointer::NSN1  
     real,dimension(:,:)    ,pointer::NSN2 
  end type yield_table
  type(yield_table)::yieldtab


  contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  ! iSIDM state index <-> ptype mapping
  function isidm_state_to_ptype(istate) result(pt)
    integer, intent(in) :: istate
    integer(kind=1) :: pt
    select case(istate)
    case(0);  pt = PTYPE_DM
    case(1);  pt = PTYPE_ISIDM_EX1
    case(2);  pt = PTYPE_ISIDM_EX2
    case default; pt = PTYPE_DM
    end select
  end function isidm_state_to_ptype

  function ptype_to_isidm_state(pt) result(istate)
    integer(kind=1), intent(in) :: pt
    integer :: istate
    select case(pt)
    case(PTYPE_DM);        istate = 0
    case(PTYPE_ISIDM_EX1); istate = 1
    case(PTYPE_ISIDM_EX2); istate = 2
    case default;          istate = -1
    end select
  end function ptype_to_isidm_state

end module pm_commons

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
  real(dp),allocatable,dimension(:)  ::xh2p     ! Dark-H2 nucleus fraction (aDM, non-eq)
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  ! This patch uses compact ptypep in place of the baseline typep/FAM_UNDEF
  ! structure; its empty-slot/default value is PTYPE_DM.
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
  ! Set by read_params when npartmax was supplied as the automatic-capacity
  ! sentinel.  Explicit positive npartmax values retain the fixed-capacity
  ! behaviour.
  logical::npartmax_auto=.false.
  ! The free list is built by init_tree after IC/restart loading.  Growth can
  ! therefore happen before it exists, but must append new slots once it does.
  logical::particle_free_list_ready=.false.
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

  subroutine grow_particle_bundle(new_npartmax)
    ! Grow every particle-sized array as one logical bundle.  The arrays are
    ! resized one at a time: move_alloc briefly retains the old and new copy
    ! of only the current array, bounding transient extra memory by one array
    ! rather than the whole particle bundle.
    integer,intent(in)::new_npartmax
    integer::old_npartmax,target_npartmax,needed,headroom,chunk
    integer::old_tail,nnew,ip
    real(dp),allocatable::new_xp(:,:),new_vp(:,:),new_mp(:)
    real(dp),allocatable::new_tp(:),new_zp(:),new_edp(:),new_xh2p(:)
    real(dp),allocatable::new_tpp(:),new_mp0(:),new_indtab(:)
    real(dp),allocatable::new_ptcl_phi(:)
    real(dp),allocatable::new_weightp(:,:)
    integer,allocatable::new_nextp(:),new_prevp(:),new_levelp(:)
    integer(i8b),allocatable::new_idp(:)
    integer(kind=1),allocatable::new_ptypep(:)

    if(new_npartmax<=npartmax)return

    old_npartmax=npartmax
    needed=new_npartmax-old_npartmax
    headroom=max(1,needed/4)
    chunk=max(1,old_npartmax/2)
    target_npartmax=max(new_npartmax,old_npartmax+chunk, &
         old_npartmax+needed+headroom)

    ! The old bundle remains untouched until each replacement has been fully
    ! allocated, copied, and initialized.  No capacity scalar is changed
    ! until the final array below has been moved into place.
    allocate(new_xp(target_npartmax,ndim))
    if(old_npartmax>0)new_xp(1:old_npartmax,:)=xp
    if(target_npartmax>old_npartmax)new_xp(old_npartmax+1:target_npartmax,:)=0d0
    call move_alloc(new_xp,xp)

    allocate(new_vp(target_npartmax,ndim))
    if(old_npartmax>0)new_vp(1:old_npartmax,:)=vp
    if(target_npartmax>old_npartmax)new_vp(old_npartmax+1:target_npartmax,:)=0d0
    call move_alloc(new_vp,vp)

    allocate(new_mp(target_npartmax))
    if(old_npartmax>0)new_mp(1:old_npartmax)=mp
    new_mp(old_npartmax+1:target_npartmax)=0d0
    call move_alloc(new_mp,mp)

    allocate(new_nextp(target_npartmax))
    if(old_npartmax>0)new_nextp(1:old_npartmax)=nextp
    new_nextp(old_npartmax+1:target_npartmax)=0
    call move_alloc(new_nextp,nextp)

    allocate(new_prevp(target_npartmax))
    if(old_npartmax>0)new_prevp(1:old_npartmax)=prevp
    new_prevp(old_npartmax+1:target_npartmax)=0
    call move_alloc(new_prevp,prevp)

    allocate(new_levelp(target_npartmax))
    if(old_npartmax>0)new_levelp(1:old_npartmax)=levelp
    new_levelp(old_npartmax+1:target_npartmax)=0
    call move_alloc(new_levelp,levelp)

    allocate(new_idp(target_npartmax))
    if(old_npartmax>0)new_idp(1:old_npartmax)=idp
    new_idp(old_npartmax+1:target_npartmax)=0_i8b
    call move_alloc(new_idp,idp)

    allocate(new_ptypep(target_npartmax))
    if(old_npartmax>0)new_ptypep(1:old_npartmax)=ptypep
    new_ptypep(old_npartmax+1:target_npartmax)=PTYPE_DM
    call move_alloc(new_ptypep,ptypep)

#ifdef OUTPUT_PARTICLE_POTENTIAL
    if(allocated(ptcl_phi))then
       allocate(new_ptcl_phi(target_npartmax))
       if(old_npartmax>0)new_ptcl_phi(1:old_npartmax)=ptcl_phi
       new_ptcl_phi(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_ptcl_phi,ptcl_phi)
    endif
#endif

    if(allocated(tp))then
       allocate(new_tp(target_npartmax))
       if(old_npartmax>0)new_tp(1:old_npartmax)=tp
       new_tp(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_tp,tp)
    endif

    if(allocated(zp))then
       allocate(new_zp(target_npartmax))
       if(old_npartmax>0)new_zp(1:old_npartmax)=zp
       new_zp(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_zp,zp)
    endif

    if(allocated(edp))then
       allocate(new_edp(target_npartmax))
       if(old_npartmax>0)new_edp(1:old_npartmax)=edp
       new_edp(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_edp,edp)
    endif

    if(allocated(xh2p))then
       allocate(new_xh2p(target_npartmax))
       if(old_npartmax>0)new_xh2p(1:old_npartmax)=xh2p
       new_xh2p(old_npartmax+1:target_npartmax)=adm_fH2
       call move_alloc(new_xh2p,xh2p)
    endif

    if(allocated(tpp))then
       allocate(new_tpp(target_npartmax))
       if(old_npartmax>0)new_tpp(1:old_npartmax)=tpp
       new_tpp(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_tpp,tpp)
    endif

    if(allocated(mp0))then
       allocate(new_mp0(target_npartmax))
       if(old_npartmax>0)new_mp0(1:old_npartmax)=mp0
       new_mp0(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_mp0,mp0)
    endif

    if(allocated(indtab))then
       allocate(new_indtab(target_npartmax))
       if(old_npartmax>0)new_indtab(1:old_npartmax)=indtab
       new_indtab(old_npartmax+1:target_npartmax)=0d0
       call move_alloc(new_indtab,indtab)
    endif

    ! weightp is a particle-sized extension in pm_commons but is not allocated
    ! by this branch.  Preserve and grow it if a sink backend has allocated it.
    if(allocated(weightp))then
       allocate(new_weightp(target_npartmax,size(weightp,2)))
       if(old_npartmax>0)new_weightp(1:old_npartmax,:)=weightp
       new_weightp(old_npartmax+1:target_npartmax,:)=0d0
       call move_alloc(new_weightp,weightp)
    endif

    if(particle_free_list_ready)then
       nnew=target_npartmax-old_npartmax
       old_tail=tailp_free
       if(numbp_free>0)then
          nextp(old_tail)=old_npartmax+1
          prevp(old_npartmax+1)=old_tail
       else
          headp_free=old_npartmax+1
          prevp(headp_free)=0
       endif
       do ip=old_npartmax+1,target_npartmax-1
          nextp(ip)=ip+1
          prevp(ip+1)=ip
       enddo
       tailp_free=target_npartmax
       nextp(tailp_free)=0
       numbp_free=numbp_free+nnew
    endif

    ! This is deliberately the last capacity update: all allocated bundle
    ! members have the same extent before npartmax changes.
    npartmax=target_npartmax
    if(particle_free_list_ready)npart=npartmax-numbp_free
  end subroutine grow_particle_bundle

  ! Count particles in the actual child cell that contains them.  numbp is
  ! a grid total and must not be charged to every leaf cell during domain
  ! decomposition.  The SIDM count follows sidm_scatter's definition of a
  ! DM-like particle so the pair proxy represents the work really attempted.
  subroutine count_particles_by_leaf(igrid,npart_leaf,ndm_leaf)
    use amr_commons, only: xg
    implicit none
    integer,intent(in)::igrid
    integer,dimension(1:twotondim),intent(out)::npart_leaf,ndm_leaf
    integer::ipart,jpart,ind,ix,iy,iz

    npart_leaf=0
    ndm_leaf=0
    if(igrid<=0 .or. .not.allocated(headp) .or. .not.allocated(numbp) .or. &
         .not.allocated(nextp) .or. .not.allocated(xp))return

    ipart=headp(igrid)
    do jpart=1,numbp(igrid)
       if(ipart<=0)exit
       ix=0
       iy=0
       iz=0
       if(xp(ipart,1)>xg(igrid,1))ix=1
#if NDIM>1
       if(xp(ipart,2)>xg(igrid,2))iy=1
#endif
#if NDIM>2
       if(xp(ipart,3)>xg(igrid,3))iz=1
#endif
       ind=1+ix+2*iy+4*iz
       npart_leaf(ind)=npart_leaf(ind)+1
       if(allocated(idp).and.allocated(ptypep))then
          if(idp(ipart)>0 .and. ptypep(ipart)/=PTYPE_STAR .and. &
               ptypep(ipart)/=PTYPE_SINK) ndm_leaf(ind)=ndm_leaf(ind)+1
       end if
       ipart=nextp(ipart)
    end do
  end subroutine count_particles_by_leaf

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

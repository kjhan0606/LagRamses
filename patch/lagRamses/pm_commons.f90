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
  ! Physical erg accepted during gas accretion, awaiting SNRT source commit.
  ! Persisted with the HDF5 sink payload; MPI source ownership is separate.
  real(dp),allocatable,dimension(:)::agn_pending_erg
  logical::agn_checkpoint_restored=.false.
  ! Reference model only: (heat erg, jet erg, retained loading mass in code
  ! units, deferred erg) by sink slot. Legacy/MAD leaves all four zero.
  real(dp),allocatable,dimension(:,:)::agn_mechanical_pending
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
  ! Particle slot of the canonical sink for each sink identity.  Cloud
  ! particles intentionally share PTYPE_SINK, so position-based selection is
  ! not stable while a sink moves or after a restart.
  integer,allocatable,dimension(:)::canonical_sink_part

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
  ! Grow particle storage at runtime when the current capacity is exhausted.
  ! This is the default; AMR_PARAMS can disable it explicitly.
  logical::npartmax_auto=.true.
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
    use amr_commons, only: myid
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
    write(*,'(A,I0,A,I0,A,I0)') &
         '[RESIZABLE] PARTICLE_GROW rank=',myid,' old=',old_npartmax, &
         ' new=',target_npartmax
    npartmax=target_npartmax
    if(particle_free_list_ready)npart=npartmax-numbp_free
  end subroutine grow_particle_bundle

  ! [RESIZABLE] Grow every grid- and cell-capacity array as one bundle.
  subroutine grow_grid_capacity(new_ngridmax)
    use amr_commons, only: myid, ncoarse, ngridmax, twotondim, amr_block_size, &
         xg, nbor, father, next, prev, son, flag1, flag2, cpu_map, cpu_map2, &
         varcpu_grid_file_idx, hilbert_key, headf, tailf, numbf, used_mem
    use hydro_commons, only: uold, unew, divu, enew
    use poisson_commons, only: lookup_mg, rho, rho_star, phi, phi_old, f, &
         scalar_gr, scalar_gr_old, psi_re, psi_im, rho_top
    use radiation_commons, only: Erad, Srad
    use morton_hash, only: grid_level
    implicit none
    ! Grow the AMR grid/cell storage one array at a time.  Each replacement
    ! retains the old allocation until its replacement has been copied and
    ! moved into place, so the transient extra memory is bounded by one array.
    ! Callers must enter this routine at a serial AMR safe point.
    integer,intent(in)::new_ngridmax
    integer::old_ngridmax,target_ngridmax,old_ncell,new_ncell,nnew,igrid
    integer::old_tail
    real(dp),allocatable::new_xg(:,:)
    integer,allocatable::new_nbor(:,:),new_father(:),new_next(:),new_prev(:)
    integer,allocatable::new_son(:),new_flag1(:),new_flag2(:)
    integer,allocatable::new_cpu_map(:),new_cpu_map2(:)
    integer,allocatable::new_headp(:),new_tailp(:),new_numbp(:)
    integer,allocatable::new_lookup_mg(:),new_grid_level(:)
    integer,allocatable::new_varcpu_grid_file_idx(:)
    real(qdp),allocatable::new_hilbert_key(:)
    real(dp),allocatable::new_uold(:,:),new_unew(:,:),new_divu(:),new_enew(:)
    real(dp),allocatable::new_rho(:),new_rho_star(:),new_phi(:),new_phi_old(:)
    real(dp),allocatable::new_f(:,:),new_rho_top(:)
    real(dp),allocatable::new_scalar_gr(:),new_scalar_gr_old(:)
    real(dp),allocatable::new_psi_re(:),new_psi_im(:)
    real(dp),allocatable::new_Erad(:),new_Srad(:)

    if(new_ngridmax<=ngridmax)return
    if(amr_block_size<=0)stop 'grow_grid_capacity: invalid amr_block_size'
    if(mod(ngridmax,amr_block_size)/=0)stop &
         'grow_grid_capacity: ngridmax is not a whole number of blocks'

    old_ngridmax=ngridmax
    target_ngridmax=max(new_ngridmax,old_ngridmax)
    if(mod(target_ngridmax,amr_block_size)/=0)then
       target_ngridmax=((target_ngridmax/amr_block_size)+1)*amr_block_size
    endif
    if(target_ngridmax<=old_ngridmax)return

    old_ncell=ncoarse+twotondim*old_ngridmax
    new_ncell=ncoarse+twotondim*target_ngridmax
    nnew=target_ngridmax-old_ngridmax

    ! Grid-indexed AMR arrays.
    if(allocated(xg))then
       allocate(new_xg(target_ngridmax,size(xg,2)))
       new_xg=0.0d0
       if(old_ngridmax>0)new_xg(1:old_ngridmax,:)=xg
       call move_alloc(new_xg,xg)
    endif

    if(allocated(father))then
       allocate(new_father(target_ngridmax))
       new_father=0
       if(old_ngridmax>0)new_father(1:old_ngridmax)=father
       call move_alloc(new_father,father)
    endif

    if(allocated(nbor))then
       allocate(new_nbor(target_ngridmax,size(nbor,2)))
       new_nbor=0
       if(old_ngridmax>0)new_nbor(1:old_ngridmax,:)=nbor
       call move_alloc(new_nbor,nbor)
    endif

    if(allocated(next))then
       allocate(new_next(target_ngridmax))
       new_next=0
       if(old_ngridmax>0)new_next(1:old_ngridmax)=next
       call move_alloc(new_next,next)
    endif

    if(allocated(prev))then
       allocate(new_prev(target_ngridmax))
       new_prev=0
       if(old_ngridmax>0)new_prev(1:old_ngridmax)=prev
       call move_alloc(new_prev,prev)
    endif

    if(allocated(headp))then
       allocate(new_headp(target_ngridmax))
       new_headp=0
       if(old_ngridmax>0)new_headp(1:old_ngridmax)=headp
       call move_alloc(new_headp,headp)
    endif

    if(allocated(tailp))then
       allocate(new_tailp(target_ngridmax))
       new_tailp=0
       if(old_ngridmax>0)new_tailp(1:old_ngridmax)=tailp
       call move_alloc(new_tailp,tailp)
    endif

    if(allocated(numbp))then
       allocate(new_numbp(target_ngridmax))
       new_numbp=0
       if(old_ngridmax>0)new_numbp(1:old_ngridmax)=numbp
       call move_alloc(new_numbp,numbp)
    endif

    if(allocated(lookup_mg))then
       allocate(new_lookup_mg(target_ngridmax))
       new_lookup_mg=0
       if(old_ngridmax>0)new_lookup_mg(1:old_ngridmax)=lookup_mg
       call move_alloc(new_lookup_mg,lookup_mg)
    endif

    if(allocated(grid_level))then
       allocate(new_grid_level(target_ngridmax))
       new_grid_level=0
       if(old_ngridmax>0)new_grid_level(1:old_ngridmax)=grid_level
       call move_alloc(new_grid_level,grid_level)
    endif

    ! This HDF5 restart-only mapping is normally released before refinement,
    ! but it is capacity-indexed if a caller grows while it is still live.
    if(allocated(varcpu_grid_file_idx))then
       allocate(new_varcpu_grid_file_idx(target_ngridmax))
       new_varcpu_grid_file_idx=0
       if(old_ngridmax>0)new_varcpu_grid_file_idx(1:old_ngridmax)= &
            varcpu_grid_file_idx
       call move_alloc(new_varcpu_grid_file_idx,varcpu_grid_file_idx)
    endif

    ! hilbert_key(1:1) is a k-section scratch allocation, not a cell array.
   if(allocated(hilbert_key).and.size(hilbert_key)>1.and.&
      size(hilbert_key)==old_ncell)then
       allocate(new_hilbert_key(new_ncell))
       new_hilbert_key=0.0_qdp
       if(old_ncell>0)new_hilbert_key(1:old_ncell)=hilbert_key
       call move_alloc(new_hilbert_key,hilbert_key)
    endif

    ! Cell-indexed AMR topology and domain-decomposition arrays.  flag1 and
    ! flag2 retain their index-zero sentinel; all newly appended cells start
    ! with the same zero state as a fresh init_amr allocation.
    if(allocated(son))then
       allocate(new_son(new_ncell))
       new_son=0
       if(old_ncell>0)new_son(1:old_ncell)=son
       call move_alloc(new_son,son)
    endif

    if(allocated(flag1))then
       allocate(new_flag1(0:new_ncell))
       new_flag1=0
       if(old_ncell>=0)new_flag1(0:old_ncell)=flag1
       call move_alloc(new_flag1,flag1)
    endif

    if(allocated(flag2))then
       allocate(new_flag2(0:new_ncell))
       new_flag2=0
       if(old_ncell>=0)new_flag2(0:old_ncell)=flag2
       call move_alloc(new_flag2,flag2)
    endif

    if(allocated(cpu_map))then
       allocate(new_cpu_map(new_ncell))
       new_cpu_map=0
       if(old_ncell>0)new_cpu_map(1:old_ncell)=cpu_map
       call move_alloc(new_cpu_map,cpu_map)
    endif

    if(allocated(cpu_map2))then
       allocate(new_cpu_map2(new_ncell))
       new_cpu_map2=0
       if(old_ncell>0)new_cpu_map2(1:old_ncell)=cpu_map2
       call move_alloc(new_cpu_map2,cpu_map2)
    endif

    ! Hydro fields are allocated only when the corresponding configuration is
    ! active.  Their first dimension is the AMR cell capacity.
    if(allocated(uold))then
       allocate(new_uold(new_ncell,size(uold,2)))
       new_uold=0.0d0
       if(old_ncell>0)new_uold(1:old_ncell,:)=uold
       call move_alloc(new_uold,uold)
    endif

    if(allocated(unew))then
       allocate(new_unew(new_ncell,size(unew,2)))
       new_unew=0.0d0
       if(old_ncell>0)new_unew(1:old_ncell,:)=unew
       call move_alloc(new_unew,unew)
    endif

    if(allocated(divu))then
       allocate(new_divu(new_ncell))
       new_divu=0.0d0
       if(old_ncell>0)new_divu(1:old_ncell)=divu
       call move_alloc(new_divu,divu)
    endif

    if(allocated(enew))then
       allocate(new_enew(new_ncell))
       new_enew=0.0d0
       if(old_ncell>0)new_enew(1:old_ncell)=enew
       call move_alloc(new_enew,enew)
    endif

    ! Poisson fields, including the conditional modified-gravity, FDM, and
    ! CIC arrays, are grown only when init_poisson actually allocated them.
    if(allocated(rho))then
       allocate(new_rho(new_ncell))
       new_rho=0.0d0
       if(old_ncell>0)new_rho(1:old_ncell)=rho
       call move_alloc(new_rho,rho)
    endif

    if(allocated(rho_star))then
       allocate(new_rho_star(new_ncell))
       new_rho_star=0.0d0
       if(old_ncell>0)new_rho_star(1:old_ncell)=rho_star
       call move_alloc(new_rho_star,rho_star)
    endif

    if(allocated(phi))then
       allocate(new_phi(new_ncell))
       new_phi=0.0d0
       if(old_ncell>0)new_phi(1:old_ncell)=phi
       call move_alloc(new_phi,phi)
    endif

    if(allocated(phi_old))then
       allocate(new_phi_old(new_ncell))
       new_phi_old=0.0d0
       if(old_ncell>0)new_phi_old(1:old_ncell)=phi_old
       call move_alloc(new_phi_old,phi_old)
    endif

    if(allocated(f))then
       allocate(new_f(new_ncell,size(f,2)))
       new_f=0.0d0
       if(old_ncell>0)new_f(1:old_ncell,:)=f
       call move_alloc(new_f,f)
    endif

    if(allocated(scalar_gr))then
       allocate(new_scalar_gr(new_ncell))
       new_scalar_gr=0.0d0
       if(old_ncell>0)new_scalar_gr(1:old_ncell)=scalar_gr
       call move_alloc(new_scalar_gr,scalar_gr)
    endif

    if(allocated(scalar_gr_old))then
       allocate(new_scalar_gr_old(new_ncell))
       new_scalar_gr_old=0.0d0
       if(old_ncell>0)new_scalar_gr_old(1:old_ncell)=scalar_gr_old
       call move_alloc(new_scalar_gr_old,scalar_gr_old)
    endif

    if(allocated(psi_re))then
       allocate(new_psi_re(new_ncell))
       new_psi_re=0.0d0
       if(old_ncell>0)new_psi_re(1:old_ncell)=psi_re
       call move_alloc(new_psi_re,psi_re)
    endif

    if(allocated(psi_im))then
       allocate(new_psi_im(new_ncell))
       new_psi_im=0.0d0
       if(old_ncell>0)new_psi_im(1:old_ncell)=psi_im
       call move_alloc(new_psi_im,psi_im)
    endif

    if(allocated(rho_top))then
       allocate(new_rho_top(new_ncell))
       new_rho_top=0.0d0
       if(old_ncell>0)new_rho_top(1:old_ncell)=rho_top
       call move_alloc(new_rho_top,rho_top)
    endif

    ! Radiation fields are linked into every build but remain unallocated
    ! unless the ATON/radiation path is active.
    if(allocated(Erad))then
       allocate(new_Erad(new_ncell))
       new_Erad=0.0d0
       if(old_ncell>0)new_Erad(1:old_ncell)=Erad
       call move_alloc(new_Erad,Erad)
    endif

    if(allocated(Srad))then
       allocate(new_Srad(new_ncell))
       new_Srad=0.0d0
       if(old_ncell>0)new_Srad(1:old_ncell)=Srad
       call move_alloc(new_Srad,Srad)
    endif

    ! Append the new grid indices after the existing free-list tail.  The
    ! appended chain has exactly the ascending order used by init_amr.
    old_tail=tailf
    if(numbf>0)then
       next(old_tail)=old_ngridmax+1
       prev(old_ngridmax+1)=old_tail
    else
       headf=old_ngridmax+1
       prev(headf)=0
    endif
    do igrid=old_ngridmax+1,target_ngridmax-1
       next(igrid)=igrid+1
       prev(igrid+1)=igrid
    enddo
    tailf=target_ngridmax
    next(tailf)=0
    numbf=numbf+nnew
    used_mem=target_ngridmax-numbf

    ! Capacity is deliberately updated last, after every capacity-indexed
    ! array and the free-grid list have reached the new extent.
    if(myid==1) write(*,'(A,I0,A,I0,A,I0)') &
         '[RESIZABLE] GRID_GROW rank=',myid,' old=',old_ngridmax, &
         ' new=',target_ngridmax
    ngridmax=target_ngridmax
  end subroutine grow_grid_capacity

  ! [RESIZABLE] Choose one capacity collectively, then resize every rank at
  ! the common safe point selected by the caller.
  subroutine ensure_grid_capacity_collective(required_local,context)
    use amr_commons, only: myid, ngridmax, amr_block_size, ngridmax_auto
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer,intent(in)::required_local
    character(len=*),intent(in)::context
    integer::target_local,target_global,growth_chunk
    integer::capacity_min,capacity_max,info

    if(amr_block_size<=0)stop &
         'ensure_grid_capacity_collective: invalid amr_block_size'

    target_local=ngridmax
    if(required_local>ngridmax)then
       growth_chunk=max(amr_block_size,max(1,ngridmax/4))
       target_local=max(required_local,ngridmax+growth_chunk)
       if(mod(target_local,amr_block_size)/=0) &
            target_local=((target_local/amr_block_size)+1)*amr_block_size
    endif

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(target_local,target_global,1,MPI_INTEGER,MPI_MAX, &
         MPI_COMM_WORLD,info)
#else
    target_global=target_local
#endif

    if(target_global>ngridmax)then
       if(.not.ngridmax_auto)then
          if(required_local>ngridmax) write(*,'(A,I0,A,A,A,I0,A,I0)') &
               '[RESIZABLE] GRID_CAPACITY_FIXED rank=',myid, &
               ' context=',trim(context),' required=',required_local, &
               ' capacity=',ngridmax
#ifndef WITHOUTMPI
          call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
          stop 'Increase ngridmax'
#endif
       else
          call grow_grid_capacity(target_global)
       endif
    endif

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(ngridmax,capacity_min,1,MPI_INTEGER,MPI_MIN, &
         MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(ngridmax,capacity_max,1,MPI_INTEGER,MPI_MAX, &
         MPI_COMM_WORLD,info)
#else
    capacity_min=ngridmax
    capacity_max=ngridmax
#endif
    if(capacity_min/=capacity_max)then
       write(*,'(A,I0,A,A,A,I0,A,I0)') &
            '[RESIZABLE] GRID_CAPACITY_MISMATCH rank=',myid, &
            ' context=',trim(context),' min=',capacity_min,' max=',capacity_max
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
       stop 'grid capacity mismatch'
#endif
    endif
  end subroutine ensure_grid_capacity_collective

  ! Count particles in the actual child cell that contains them.  numbp is
  ! a grid total and must not be charged to every leaf cell during domain
  ! decomposition.  The SIDM count follows sidm_scatter's definition of a
  ! DM-like particle so the pair proxy represents the work really attempted.
  subroutine count_particles_by_leaf(igrid,npart_leaf,ndm_leaf,count_dm)
    use amr_commons, only: xg
    implicit none
    integer,intent(in)::igrid
    integer,dimension(1:twotondim),intent(out)::npart_leaf,ndm_leaf
    logical,intent(in),optional::count_dm
    integer::ipart,jpart,ind,ix,iy,iz,np,ivec
    integer,dimension(1:nvector)::ind_part,ind_leaf
    logical,dimension(1:nvector)::dm_like
    logical::need_dm

    npart_leaf=0
    ndm_leaf=0
    if(igrid<=0 .or. .not.allocated(headp) .or. .not.allocated(numbp) .or. &
         .not.allocated(nextp) .or. .not.allocated(xp))return

    need_dm=.true.
    if(present(count_dm))need_dm=count_dm
    need_dm=need_dm.and.allocated(idp).and.allocated(ptypep)

    ! Linked-list chasing is inherently serial, but gather particle indices in
    ! NVECTOR chunks so the coordinate/type classification below can be
    ! vectorized and does not interleave random nextp and xp/idp loads.
    ipart=headp(igrid)
    jpart=0
    do while(jpart<numbp(igrid).and.ipart>0)
       np=0
       do while(np<nvector.and.jpart<numbp(igrid).and.ipart>0)
          np=np+1
          jpart=jpart+1
          ind_part(np)=ipart
          ipart=nextp(ipart)
       end do
       !$OMP SIMD PRIVATE(ix,iy,iz)
       do ivec=1,np
          ix=0
          iy=0
          iz=0
          if(xp(ind_part(ivec),1)>xg(igrid,1))ix=1
#if NDIM>1
          if(xp(ind_part(ivec),2)>xg(igrid,2))iy=1
#endif
#if NDIM>2
          if(xp(ind_part(ivec),3)>xg(igrid,3))iz=1
#endif
          ind_leaf(ivec)=1+ix+2*iy+4*iz
          dm_like(ivec)=.false.
          if(need_dm)then
             dm_like(ivec)=idp(ind_part(ivec))>0 .and. &
                  ptypep(ind_part(ivec))/=PTYPE_STAR .and. &
                  ptypep(ind_part(ivec))/=PTYPE_SINK
          end if
       end do
       do ivec=1,np
          ind=ind_leaf(ivec)
          npart_leaf(ind)=npart_leaf(ind)+1
          if(dm_like(ivec))ndm_leaf(ind)=ndm_leaf(ind)+1
       end do
    end do
  end subroutine count_particles_by_leaf

  ! Spread an exact aggregate count over the selected direct leaf children.
  ! The quotient/remainder rule is deterministic and conserves npart_direct;
  ! it deliberately approximates only the within-grid particle positions.
  subroutine distribute_particle_total_by_leaf(npart_direct,is_leaf,npart_leaf)
    implicit none
    integer(kind=8),intent(in)::npart_direct
    logical,dimension(1:twotondim),intent(in)::is_leaf
    integer,dimension(1:twotondim),intent(out)::npart_leaf
    integer::ind,nleaf,ileaf,nbase,nextra

    npart_leaf=0
    nleaf=count(is_leaf)
    if(npart_direct<=0_8.or.nleaf<=0)return
    if(npart_direct>int(huge(nbase),kind=8))then
       write(*,*)'distribute_particle_total_by_leaf: integer overflow ', &
            npart_direct
       error stop 1
    end if
    nbase=int(npart_direct/int(nleaf,kind=8))
    nextra=int(mod(npart_direct,int(nleaf,kind=8)))
    ileaf=0
    do ind=1,twotondim
       if(is_leaf(ind))then
          ileaf=ileaf+1
          npart_leaf(ind)=nbase
          if(ileaf<=nextra)npart_leaf(ind)=npart_leaf(ind)+1
       end if
    end do
  end subroutine distribute_particle_total_by_leaf

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

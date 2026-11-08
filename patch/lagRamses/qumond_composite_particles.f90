! Controlled serial, fixed-mesh particle coupling to the native leaf solver.
! A fixed physical cloud is integrated over leaves. The identical weights
! deposit particle mass and gather cell acceleration; no uniform field bridge.
module qumond_composite_particles
  use amr_commons
  use pm_commons
  use qumond_composite, only: composite_mesh_t
  use qumond_composite_clouds, only: native_cloud_map
  use poisson_commons, only: qh_sink_sources
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: native_particles_source,native_particles_force,native_particles_cached
  public :: native_particles_begin,native_particles_drift,native_particles_finish,native_particles_export
  public :: native_stellar_begin,native_stellar_reconcile,native_stellar_identity
  integer,allocatable :: slots(:),ptr(:),leaf(:),slot_index(:)
  integer(i8b),allocatable :: identities(:)
  integer(kind=1),allocatable :: species(:)
  real(dp),allocatable :: masses(:),weights(:),source_x(:,:),accel(:,:)
  real(dp),allocatable :: start_x(:,:),start_v(:,:),start_a(:,:),half_v(:,:)
  real(dp) :: cloud_width=0d0,domain_width=0d0,displacement_limit=0d0
  integer :: np=0,nweight=0,live_count=0
  logical :: initialized=.false.,step_open=.false.,drifted=.false.,stellar_open=.false.
  integer(i8b) :: next_stellar_identity=0_i8b
  character(len=32) :: source_role=''
contains
  subroutine fail(reason)
    character(len=*),intent(in) :: reason
    write(*,*) 'Native composite particles failed: ',trim(reason)
    call clean_stop
  end subroutine

  subroutine initialize_roster(mesh,is_mond,canonical)
    type(composite_mesh_t),intent(in) :: mesh
    logical,intent(in) :: is_mond
    logical,intent(in),optional :: canonical
    integer :: lev,g,grid,p,j,count,status,ios,last_level,total_live
    character(len=128) :: width_env
    integer(i8b),allocatable :: sorted(:)
    if(.not.allocated(idp).or..not.allocated(xp).or..not.allocated(vp).or..not.allocated(mp)) &
         call fail('particle arrays are not allocated')
    np=npart;live_count=npart
    if(np<0.or.np>size(idp))call fail('invalid live particle count')
    allocate(slots(np),identities(np),species(np),masses(np),slot_index(size(idp)))
    slot_index=0;count=0;total_live=0
    ! Called after the root kill_tree, before recursion: the level lists
    ! partition the live roster. Free slots and merged child aliases are not
    ! inferred from a default particle species or an array prefix.
    last_level=nlevelmax
    if(present(canonical))then
       if(canonical)last_level=levelmin
    endif
    do lev=levelmin,last_level
       do g=1,active(lev)%ngrid
          grid=active(lev)%igrid(g);p=headp(grid)
          do j=1,numbp(grid)
             if(p<1.or.p>size(idp))call fail('invalid live linked-list slot')
             if(slot_index(p)/=0.or.total_live>=live_count)call fail('duplicate particle or count mismatch')
             total_live=total_live+1
             if(ptypep(p)==PTYPE_SINK)then
                if(.not.qh_sink_sources.or.idp(p)>=0.or.idp(p)<-nsink)call fail('foreign sink cloud')
                slot_index(p)=-1
             else
                count=count+1;slots(count)=p;slot_index(p)=count
                identities(count)=idp(p);species(count)=ptypep(p);masses(count)=mp(p)
             endif
             p=nextp(p)
          enddo
          if(p/=0)call fail('particle list exceeds declared length')
       enddo
    enddo
    if(total_live/=live_count)call fail('live tree roster does not match npart')
    np=count;slots=slots(:np);identities=identities(:np);species=species(:np);masses=masses(:np)
    allocate(source_x(np,3),accel(np,3),start_x(np,3),start_v(np,3),start_a(np,3),half_v(np,3))
    if(any(identities<=0))call fail('nonpositive live particle identity')
    if(any(species/=PTYPE_DM.and.species/=PTYPE_STAR))call fail('unsupported particle species')
    allocate(sorted(np));sorted=identities
    if(np>1)call sort_ids(sorted,1,np)
    do j=2,np
       if(sorted(j)==sorted(j-1))call fail('duplicate live particle identity')
    enddo
    if(.not.initialized)then
       source_role=''
       call get_environment_variable('RAR_NATIVE_PARTICLE_ROLE',source_role,status=status)
       if(status/=0.or.len_trim(source_role)==0)then
          source_role='baryon'
          if(any(species==PTYPE_DM))source_role='cdm'
       endif
    endif
    select case(trim(source_role))
    case('baryon')
       if(any(species/=PTYPE_STAR))call fail('baryon role requires stellar particle species')
    case('cdm')
       if(is_mond)call fail('CDM particles are excluded from the QUMOND science source')
    case('diagnostic')
       ! Manufactured gravitating particles; no galaxy-model interpretation.
    case default
       call fail('unknown particle source role')
    end select
    if(.not.initialized)then
    cloud_width=mesh%box*0.5d0**levelmin
    width_env=''
    call get_environment_variable('RAR_NATIVE_PARTICLE_CLOUD_WIDTH',width_env,status=status)
    if(status==-1)call fail('particle cloud width setting is truncated')
    if(status==0.and.len_trim(width_env)>0)then
       read(width_env,*,iostat=ios)cloud_width
       if(ios/=0)call fail('invalid particle cloud width setting')
    endif
    if(.not.ieee_is_finite(cloud_width))call fail('nonfinite particle cloud width')
    if(cloud_width<mesh%box*0.5d0**levelmin.or.cloud_width>mesh%box) &
         call fail('particle cloud width must span at least a base cell and at most the domain')
    domain_width=mesh%box
    displacement_limit=0.5d0*minval(mesh%width)
    endif
    source_x=-huge(1d0);accel=0d0
    initialized=.true.
    call validate_roster()
    write(*,*) 'RAR_NATIVE_PARTICLE_MODEL_V1',np,trim(source_role),cloud_width, &
         count_species(PTYPE_DM),count_species(PTYPE_STAR)
  end subroutine

  integer function count_species(tag) result(n)
    integer(kind=1),intent(in) :: tag
    n=count(species==tag)
  end function

  recursive subroutine sort_ids(a,lo,hi,order)
    integer(i8b),intent(inout) :: a(:)
    integer,intent(in) :: lo,hi
    integer,intent(inout),optional :: order(:)
    integer :: i,j,itmp
    integer(i8b) :: pivot,tmp
    i=lo;j=hi;pivot=a((lo+hi)/2)
    do
       do while(a(i)<pivot)
          i=i+1
       enddo
       do while(a(j)>pivot)
          j=j-1
       enddo
       if(i>j)exit
       tmp=a(i);a(i)=a(j);a(j)=tmp
       if(present(order))then
          itmp=order(i);order(i)=order(j);order(j)=itmp
       endif
       i=i+1;j=j-1
       if(i>j)exit
    enddo
    if(lo<j)call sort_ids(a,lo,j,order)
    if(i<hi)call sort_ids(a,i,hi,order)
  end subroutine

  subroutine native_stellar_begin()
    call validate_roster()
    if(step_open.or.stellar_open)call fail('stellar transaction overlaps a particle step')
    if(.not.star.or..not.allocated(tpp).or..not.allocated(mp0).or..not.allocated(indtab)) &
         call fail('stellar transaction requires initialized SSP particle arrays')
    if(.not.native_particles_cached())call fail('stellar transaction begins with a stale source')
    start_x=xp(slots,:);start_v=vp(slots,:)
    next_stellar_identity=0_i8b
    if(np>0)next_stellar_identity=maxval(identities)
    stellar_open=.true.
  end subroutine

  integer(i8b) function native_stellar_identity() result(identity)
    if(.not.stellar_open)call fail('stellar identity requested outside a source transaction')
    if(next_stellar_identity==huge(next_stellar_identity))call fail('stellar identity range exhausted')
    next_stellar_identity=next_stellar_identity+1_i8b
    identity=next_stellar_identity
  end function

  subroutine native_stellar_reconcile(mesh,is_mond)
    type(composite_mesh_t),intent(in) :: mesh
    logical,intent(in) :: is_mond
    integer,allocatable :: current_slots(:),order(:)
    integer(i8b),allocatable :: keys(:)
    logical,allocatable :: retained(:),seen_slots(:)
    integer :: g,grid,j,p,i,k,count_live,total_live,lo,hi,mid,found,old_count
    real(dp) :: old_mass,new_mass
    if(.not.stellar_open.or.step_open)call fail('invalid stellar roster reconciliation phase')
    if(npart<np)call fail('stellar transaction deleted a live particle')
    allocate(current_slots(npart),keys(npart),order(npart),retained(npart))
    allocate(seen_slots(size(idp)));seen_slots=.false.
    count_live=0;total_live=0;retained=.false.
    ! Only the canonical base lists are authoritative after the final merge.
    do g=1,active(levelmin)%ngrid
       grid=active(levelmin)%igrid(g);p=headp(grid)
       do j=1,numbp(grid)
          if(p<1.or.p>size(idp).or.count_live>=npart)call fail('invalid post-source canonical list')
          if(seen_slots(p))call fail('duplicate post-source live slot')
          seen_slots(p)=.true.;total_live=total_live+1
          if(ptypep(p)==PTYPE_SINK)then
             if(.not.qh_sink_sources.or.idp(p)>=0.or.idp(p)<-nsink)call fail('foreign post-source cloud')
          else
             count_live=count_live+1;current_slots(count_live)=p
             keys(count_live)=idp(p);order(count_live)=count_live
          endif
          p=nextp(p)
       enddo
       if(p/=0)call fail('post-source list exceeds declared length')
    enddo
    if(total_live/=npart.or.count_live<np.or.any(keys(:count_live)<=0)) &
         call fail('invalid post-source live particle coverage')
    if(count_live>1)call sort_ids(keys,1,count_live,order)
    do j=2,count_live
       if(keys(j)==keys(j-1))call fail('duplicate post-source particle identity')
    enddo
    do i=1,np
       lo=1;hi=count_live;found=0
       do while(lo<=hi)
          mid=(lo+hi)/2
          if(keys(mid)==identities(i))then
             found=mid;exit
          else if(keys(mid)<identities(i))then
             lo=mid+1
          else
             hi=mid-1
          endif
       enddo
       if(found==0)call fail('stellar transaction lost an existing particle identity')
       k=order(found);p=current_slots(k);retained(k)=.true.
       if(ptypep(p)/=species(i))call fail('stellar transaction changed an existing source role')
       if(any(xp(p,:)/=start_x(i,:)).or.any(vp(p,:)/=start_v(i,:))) &
            call fail('isotropic stellar return changed surviving particle position or velocity')
       if(.not.ieee_is_finite(mp(p)).or.mp(p)<0d0.or.mp(p)>masses(i)) &
            call fail('stellar return must leave finite nonnegative nonincreasing stellar mass')
       if(species(i)==PTYPE_DM.and.mp(p)/=masses(i))call fail('stellar transaction changed CDM mass')
    enddo
    do j=1,count_live
       if(retained(j))cycle
       p=current_slots(j)
       if(ptypep(p)/=PTYPE_STAR.or.mp(p)<=0d0.or.mp0(p)/=mp(p).or.tpp(p)/=texp.or.indtab(p)/=0d0) &
            call fail('newborn SSP mass, birth time, role or progress is inconsistent')
    enddo
    old_count=np;old_mass=sum(masses);new_mass=sum(mp(current_slots(:count_live)))
    deallocate(slots,ptr,leaf,slot_index,identities,species,masses,weights,source_x, &
         accel,start_x,start_v,start_a,half_v)
    ! Keep the originally selected source role, physical cloud width and domain.
    call initialize_roster(mesh,is_mond,canonical=.true.)
    call check_base_lists()
    stellar_open=.false.
    write(*,*) 'RAR_NATIVE_STELLAR_ROSTER_V1',nstep_coarse,t,old_count,np,old_mass,new_mass
  end subroutine

  subroutine validate_roster()
    if(.not.initialized)call fail('roster not initialized')
    if(npart/=live_count)call fail('particle/cloud creation outside source transaction')
    if(any(idp(slots)/=identities).or.any(ptypep(slots)/=species).or.any(mp(slots)/=masses)) &
         call fail('particle identity, source role or mass changed')
    if(.not.all(ieee_is_finite(xp(slots,:))).or..not.all(ieee_is_finite(vp(slots,:))).or. &
         .not.all(ieee_is_finite(masses)).or.any(masses<0d0))call fail('invalid particle state')
  end subroutine

  logical function native_particles_cached() result(ok)
    ok=.false.
    if(.not.initialized)return
    call validate_roster()
    ok=all(xp(slots,:)==source_x)
  end function

  subroutine native_particles_source(mesh,density,is_mond)
    type(composite_mesh_t),intent(in) :: mesh
    real(dp),intent(inout) :: density(:)
    logical,intent(in) :: is_mond
    real(dp),allocatable :: particle_rho(:)
    real(dp) :: max_error,integral,expected
    integer :: p,k,i
    if(.not.initialized)call initialize_roster(mesh,is_mond)
    call validate_roster()
    if(mesh%box/=domain_width)call fail('particle assignment domain changed')
    allocate(particle_rho(mesh%n));particle_rho=0d0
    source_x=xp(slots,:)
    call native_cloud_map(mesh,source_x,cloud_width,ptr,leaf,weights,max_error)
    nweight=size(weights)
    do p=1,np
       do k=ptr(p),ptr(p+1)-1
          i=leaf(k)
          particle_rho(i)=particle_rho(i)+masses(p)*weights(k)/mesh%volume(i)
       enddo
    enddo
    integral=sum(particle_rho*mesh%volume);expected=sum(masses)
    if(abs(integral-expected)>1d-12*max(expected,tiny(1d0)))call fail('particle deposit mass mismatch')
    density=density+particle_rho
    write(*,*) 'RAR_NATIVE_PARTICLE_SOURCE_V1',nstep_coarse,t,np,nweight,expected,integral,max_error
  end subroutine

  subroutine native_particles_force(mesh,force)
    type(composite_mesh_t),intent(in) :: mesh
    real(dp),intent(in) :: force(:,:)
    integer :: p,k
    call validate_roster()
    if(.not.native_particles_cached())call fail('force gather has stale source weights')
    accel=0d0
    do p=1,np
       do k=ptr(p),ptr(p+1)-1
          accel(p,:)=accel(p,:)+weights(k)*force(mesh%cell(leaf(k)),:)
       enddo
    enddo
    if(.not.all(ieee_is_finite(accel)))call fail('nonfinite gathered acceleration')
  end subroutine

  subroutine native_particles_begin(dt)
    real(dp),intent(inout) :: dt
    real(dp) :: speed,acc,bound
    integer :: p,d
    call validate_roster()
    if(step_open.or.stellar_open)call fail('particle step overlaps an unfinished transaction')
    start_x=xp(slots,:);start_v=vp(slots,:);start_a=accel
    do p=1,np
       do d=1,3
          speed=abs(start_v(p,d));acc=abs(start_a(p,d))
          if(speed==0d0.and.acc==0d0)cycle
          ! dt*(|v|+dt*|a|/2)<=distance bounds the prospective drift.
          bound=2d0*displacement_limit/(speed+sqrt(speed**2+2d0*acc*displacement_limit))
          dt=min(dt,bound)
       enddo
    enddo
    step_open=.true.;drifted=.false.
  end subroutine

  subroutine check_base_lists()
    logical,allocatable :: seen(:)
    integer :: g,grid,p,j,i,previous,count_live,count_body
    allocate(seen(size(idp)));seen=.false.;count_live=0;count_body=0
    do g=1,active(levelmin)%ngrid
       grid=active(levelmin)%igrid(g);p=headp(grid);previous=0
       do j=1,numbp(grid)
          if(p<1.or.p>size(slot_index))call fail('invalid canonical particle slot')
          i=slot_index(p)
          if(i==0)call fail('canonical list contains a particle outside roster')
          if(i<0)then
             if(.not.qh_sink_sources.or.ptypep(p)/=PTYPE_SINK.or.idp(p)>=0.or.idp(p)<-nsink) &
                  call fail('invalid canonical cloud')
          else
             count_body=count_body+1
          endif
          if(seen(p).or.prevp(p)/=previous)call fail('canonical particle list is duplicated or broken')
          seen(p)=.true.;count_live=count_live+1;previous=p;p=nextp(p)
       enddo
       if(p/=0.or.tailp(grid)/=previous)call fail('canonical list tail/count mismatch')
    enddo
    if(count_live/=live_count.or.count_body/=np.or..not.all(seen(slots))) &
         call fail('canonical lists do not cover roster')
  end subroutine

  subroutine native_particles_drift(dt)
    real(dp),intent(in) :: dt
    call validate_roster()
    if(.not.step_open.or.drifted)call fail('invalid particle drift phase')
    if(any(xp(slots,:)/=start_x).or.any(vp(slots,:)/=start_v)) &
         call fail('legacy consumer changed particles within the common step')
    call check_base_lists()
    half_v=start_v+0.5d0*dt*start_a
    if(np>0)then
       if(maxval(abs(dt*half_v))>displacement_limit*(1d0+1d-12))call fail('particle displacement bound')
    endif
    vp(slots,:)=half_v;xp(slots,:)=start_x+dt*half_v
    ! Root merge makes these canonical base lists. Keep coordinates unwrapped
    ! for make_tree's adjacent-grid lookup; it applies periodic wrapping.
    call make_tree_fine(levelmin)
    call check_base_lists()
    if(any(xp(slots,:)<0d0).or.any(xp(slots,:)>=boxlen))call fail('particle periodic relink failed')
    drifted=.true.
  end subroutine

  subroutine native_particles_finish(dt)
    real(dp),intent(in) :: dt
    integer :: p,grid,lev,child,d,unit,ios,status
    real(dp) :: q(3)
    character(len=8) :: env
    character(len=80) :: filename
    if(.not.step_open.or..not.drifted)call fail('invalid endpoint particle kick')
    call validate_roster()
    vp(slots,:)=half_v+0.5d0*dt*accel
    do p=1,np
       q=xp(slots(p),:)/boxlen;grid=son(1)
       do lev=1,nlevelmax
          child=1
          do d=1,3
             if(q(d)>=xg(grid,d))child=child+2**(d-1)
          enddo
          if(son(ICELL_OF(grid,child))==0)exit
          grid=son(ICELL_OF(grid,child))
       enddo
       if(lev>nlevelmax)call fail('particle centre has no containing leaf')
       levelp(slots(p))=lev
    enddo
    call validate_roster()
    env='';call get_environment_variable('RAR_QUMOND_GRAVITY_KICKS',env,status=status)
    if(status==0.and.trim(env)=='1')then
       write(filename,'("particle_step_",I5.5,".tsv")')nstep_coarse
       open(newunit=unit,file=trim(filename),status='new',action='write',iostat=ios)
       if(ios/=0)call fail('particle transaction export already exists or cannot be created')
       write(unit,'(A)') '# slot id species mass dt x0[3] v0[3] a0[3] vhalf[3] x1[3] a1[3] v1[3]'
       do p=1,np
          write(unit,'(3(I0,1X),23(ES24.16,1X))')slots(p),identities(p),species(p),masses(p),dt, &
               start_x(p,:),start_v(p,:),start_a(p,:),half_v(p,:),xp(slots(p),:),accel(p,:),vp(slots(p),:)
       enddo
       close(unit)
    endif
    write(*,*) 'RAR_NATIVE_PARTICLE_END_V1',nstep_coarse,t,np,sum(masses)
    step_open=.false.;drifted=.false.
  end subroutine

  subroutine native_particles_export(epoch)
    integer,intent(in) :: epoch
    integer :: unit,ios,p
    character(len=80) :: filename
    write(filename,'("particles_epoch_",I5.5,".tsv")')epoch
    open(newunit=unit,file=trim(filename),status='new',action='write',iostat=ios)
    if(ios/=0)call fail('particle source export already exists or cannot be created')
    write(unit,'(A,1X,ES24.16,1X,A)') '# cloud_width',cloud_width,trim(source_role)
    write(unit,'(A,1X,I0,1X,ES24.16)') '# epoch',nstep_coarse,t
    write(unit,'(A)') '# slot id species mass x[3] v[3] acceleration[3]; endpoint source velocities precede final kick'
    do p=1,np
       write(unit,'(3(I0,1X),10(ES24.16,1X))')slots(p),identities(p),species(p),masses(p), &
            xp(slots(p),:),vp(slots(p),:),accel(p,:)
    enddo
    close(unit)
  end subroutine
end module qumond_composite_particles

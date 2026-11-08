! Native serial periodic fixed-mesh gravity publication and common-time evolution.
! Stellar endpoint transactions are separately enabled; sink gravity remains excluded.
module qumond_composite_runtime
  use amr_commons
  use poisson_commons
  use qumond_composite_particles
  use hydro_commons, only: uold,unew,smallr,nvar,courant_factor
#ifdef PHASE0_STELLAR_ENRICHMENT
  use stellar_enrichment_config, only: use_channel_resolved_feedback
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: qumond_publish_fixed_source,qumond_restore_published
  public :: qumond_evolution_begin,qumond_evolution_end
  public :: qumond_stellar_endpoint
  real(dp),allocatable :: saved_phi(:),saved_force(:,:),saved_rho(:)
  real(dp),allocatable :: saved_gas(:)
  real(dp) :: saved_time=-huge(1d0),physical_rho_max=0d0,evolution_start=0d0,evolution_dt=0d0
  logical :: evolution_open=.false.
  integer :: export_count=0,evolution_epoch=-1
  integer :: published_epoch=-1
  integer :: source_generation=0,published_generation=-1
contains
  subroutine qumond_publish_fixed_source()
    real(dp),allocatable :: density(:),rhs(:),pn(:),pq(:),phantom(:),gradient(:,:),draft_n(:)
    real(dp)::xx(nvector,3),dd(nvector),mean_density,mean_source,fourpi,a0,external(3)
    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,res_n,res_q,balance
    integer::i,j,k,lev,ids(nvector),ierr,it_n,it_q,cell,export_status,unit
    character(len=8)::export_env
    character(len=80)::leaf_file
    logical :: same_time_refresh
    if(.not.qh_publish)call clean_stop
    if(qh_composite_mesh%n==0)then
       call composite_native_build(qh_composite_mesh,ierr)
       if(ierr/=0)call fail('native geometry')
    endif
    if(nstep_coarse==0)call continuum_check()
    associate(mesh=>qh_composite_mesh)
    same_time_refresh=qh_evolve.and.saved_time==t.and.allocated(saved_gas).and. &
         source_generation/=published_generation
    if(same_time_refresh)then
       if(.not.qh_stellar_sources.or.evolution_open.or.source_generation/=published_generation+1) &
            call fail('unauthorized same-time physical source replacement')
    endif
    if(qh_evolve.and.saved_time==t.and.allocated(saved_gas).and..not.same_time_refresh)then
       if(any(uold(mesh%cell,1)/=saved_gas))call fail('completed endpoint gas source changed')
       if(pic)then
          if(.not.native_particles_cached())call fail('completed endpoint particle source changed')
       endif
       published_epoch=nstep_coarse
       qh_epoch(:,levelmin:nlevelmax)=nstep_coarse
       call qumond_restore_published()
       return
    endif
    if(qh_evolve)then
       if(.not.all(ieee_is_finite(uold(mesh%cell,:))).or.any(uold(mesh%cell,1)<smallr)) &
            call fail('invalid raw evolving gas before density floor')
       if(any(uold(mesh%cell,5)-0.5d0*sum(uold(mesh%cell,2:4)**2,dim=2)/uold(mesh%cell,1)<=0d0)) &
            call fail('nonpositive raw evolving gas internal energy')
    endif
    allocate(density(mesh%n),rhs(mesh%n),pn(mesh%n),pq(mesh%n),phantom(mesh%n),gradient(3,mesh%n))
    density=max(uold(mesh%cell,1),smallr)
    if(gravity_type<0)then
       do lev=levelmin,nlevelmax
          k=0
          do i=1,mesh%n
             if(mesh%level(i)/=lev)cycle
             k=k+1;ids(k)=i;xx(k,:)=mesh%center(:,i)
             if(k==nvector)then
                call add_analytic(k,boxlen*0.5d0**lev);k=0
             endif
          enddo
          if(k>0)call add_analytic(k,boxlen*0.5d0**lev)
       enddo
    endif
    if(pic)call native_particles_source(mesh,density,use_mond)
    if(.not.all(ieee_is_finite(density)))call fail('source is nonfinite')
    physical_rho_max=maxval(abs(density))
    fourpi=4d0*acos(-1d0)
    mean_density=sum(density*mesh%volume)/sum(mesh%volume)
    rhs=fourpi*(density-mean_density)
    pn=0d0
    if(all(qh_valid(1,levelmin:nlevelmax)))pn=qh_current(mesh%cell,1)
    call composite_solve(mesh,rhs,pn,1d-9,1200,40,it_n,res_n,ierr)
    if(ierr/=0)call fail('Newtonian solve did not converge')
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    a0=a0_mond*scale_t**2/scale_l;external=g_ext_mond*scale_t**2/scale_l
    if(.not.use_mond)then
       a0=0d0;external=0d0
    endif
    if(use_mond)then
       call composite_phantom(mesh,pn,a0,fourpi,external,phantom,ierr)
       if(ierr/=0)call fail('phantom source')
    else
       phantom=0d0
    endif
    balance=abs(sum(phantom*mesh%volume))/max(sum(abs(phantom)*mesh%volume),tiny(1d0))
    if(balance>1d-12)call fail('phantom flux is not conservative')
    rhs=rhs+fourpi*phantom
    mean_source=sum(rhs*mesh%volume)/sum(mesh%volume)
    rhs=rhs-mean_source
    pq=pn
    if(all(qh_valid(2,levelmin:nlevelmax)))pq=qh_current(mesh%cell,2)
    if(use_mond)then
       call composite_solve(mesh,rhs,pq,1d-9,1200,40,it_q,res_q,ierr)
       if(ierr/=0)call fail('QUMOND solve did not converge')
    else
       ! A matched Newtonian control uses precisely the same native operator.
       ! The second history slot mirrors N; it is not a second physics solve.
       pq=pn;it_q=0;res_q=res_n
    endif
    call composite_slopes(mesh,pq,gradient,ierr)
    if(ierr/=0)call fail('force reconstruction')
    if(.not.allocated(saved_phi))then
       allocate(saved_phi(size(phi)),saved_force(size(phi),3),saved_rho(size(phi)))
    endif
    allocate(draft_n(size(phi)))
    saved_phi=0d0;saved_force=0d0;saved_rho=0d0;draft_n=0d0
    do i=1,mesh%n
       cell=mesh%cell(i)
       draft_n(cell)=pn(i);saved_phi(cell)=pq(i)
       saved_force(cell,:)=-gradient(:,i);saved_rho(cell)=density(i)+phantom(i)
    enddo
    do lev=nlevelmax-1,1,-1
       do j=1,active(lev)%ngrid
          do k=1,twotondim
             cell=ICELL_OF(active(lev)%igrid(j),k)
             if(son(cell)>0)call restrict_cell(cell)
          enddo
       enddo
    enddo
    do cell=1,ncoarse
       if(son(cell)>0)call restrict_cell(cell)
    enddo
    if(.not.all(ieee_is_finite(saved_phi)).or..not.all(ieee_is_finite(saved_force)).or. &
         .not.all(ieee_is_finite(saved_rho)))call fail('nonfinite draft fields')
    ! Publish a complete hierarchy only after both stages converge. Particle
    ! gathers use the same fixed-cloud weights as the physical source deposit.
    if(pic)call native_particles_force(mesh,saved_force)
    do lev=levelmin,nlevelmax
       do i=1,size(phi)
          if(qh_cell_level(i)/=lev)cycle
          if(qh_valid(1,lev))then
             if(.not.same_time_refresh)qh_previous(i,1)=qh_current(i,1)
          else
             qh_previous(i,1)=draft_n(i)
          endif
          if(qh_valid(2,lev))then
             if(.not.same_time_refresh)qh_previous(i,2)=qh_current(i,2)
          else
             qh_previous(i,2)=saved_phi(i)
          endif
          qh_current(i,1)=draft_n(i);qh_current(i,2)=saved_phi(i)
          phi_old(i)=qh_previous(i,2)
       enddo
       do k=1,2
          if(qh_valid(k,lev))then
             if(.not.same_time_refresh)qh_old_time(k,lev)=qh_time(k,lev)
          else
             qh_old_time(k,lev)=t
          endif
          write(*,*) 'RAR_QUMOND_HISTORY_V1',k,lev,t,qh_old_time(k,lev),.not.qh_valid(k,lev)
          qh_time(k,lev)=t;qh_valid(k,lev)=.true.;qh_epoch(k,lev)=nstep_coarse
       enddo
    enddo
    published_epoch=nstep_coarse
    published_generation=source_generation
    saved_time=t
    if(qh_evolve)then
       if(.not.allocated(saved_gas))allocate(saved_gas(mesh%n))
       saved_gas=uold(mesh%cell,1)
    endif
    call qumond_restore_published()
    write(*,*) 'RAR_COMPOSITE_MODEL_V1',use_mond,gravity_type
    if(qh_stellar_sources)write(*,*) 'RAR_QUMOND_SOURCE_GENERATION_V1', &
         nstep_coarse,t,source_generation,same_time_refresh
    write(*,*) 'RAR_QUMOND_PUBLICATION_V1',nstep_coarse,t,mesh%n,it_n,it_q,res_n,res_q,balance, &
         sum(saved_phi(mesh%cell)*mesh%volume),sqrt(sum(mesh%volume*sum(gradient**2,dim=1)))
    export_env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_EXPORT',export_env,status=export_status)
    if(export_status==0.and.trim(export_env)=='1'.and.(nstep_coarse==0.or.qh_evolve))then
       leaf_file='composite_leaf_state.tsv'
       if(qh_evolve)write(leaf_file,'("composite_epoch_",I5.5,"_leaf.tsv")')export_count
       open(newunit=unit,file=trim(leaf_file),status='new',action='write',iostat=ierr)
       if(ierr/=0)call fail('leaf export exists or cannot be created')
       write(unit,'(A,5(1X,ES24.16))') '# constants',fourpi,a0,external
       if(qh_evolve)write(unit,'(A,1X,I0,1X,ES24.16)') '# epoch',nstep_coarse,t
       do i=1,mesh%n
          write(unit,'(2(I0,1X),11(ES24.16,1X))') mesh%cell(i),mesh%level(i),mesh%center(:,i), &
               mesh%width(i),density(i),pn(i),pq(i),phantom(i),-gradient(:,i)
       enddo
       close(unit)
       if(.not.qh_evolve.or.export_count==0)then
       open(newunit=unit,file='composite_faces.tsv',status='new',action='write',iostat=ierr)
       if(ierr/=0)call fail('face export exists or cannot be created')
       do i=1,mesh%nface
          write(unit,'(3(I0,1X),2(ES24.16,1X))') mesh%left(i),mesh%right(i),mesh%axis(i), &
               mesh%area(i),mesh%distance(i)
       enddo
       close(unit)
       endif
       if(pic)call native_particles_export(export_count)
       export_count=export_count+1
    endif
    end associate
  contains
    subroutine add_analytic(n,dx)
      integer,intent(in)::n
      real(dp),intent(in)::dx
      integer::q
      call rho_ana(xx,dd,dx,n)
      do q=1,n
         density(ids(q))=density(ids(q))+dd(q)
      enddo
    end subroutine
    subroutine restrict_cell(parent)
      integer,intent(in)::parent
      integer::c,child,grid
      grid=son(parent)
      draft_n(parent)=0d0;saved_phi(parent)=0d0;saved_force(parent,:)=0d0;saved_rho(parent)=0d0
      do c=1,twotondim
         child=ICELL_OF(grid,c)
         draft_n(parent)=draft_n(parent)+draft_n(child)/real(twotondim,dp)
         saved_phi(parent)=saved_phi(parent)+saved_phi(child)/real(twotondim,dp)
         saved_force(parent,:)=saved_force(parent,:)+saved_force(child,:)/real(twotondim,dp)
         saved_rho(parent)=saved_rho(parent)+saved_rho(child)/real(twotondim,dp)
      enddo
    end subroutine
  end subroutine qumond_publish_fixed_source

  subroutine continuum_check()
    ! Smooth periodic continuum solution: phi_N=cos(k*x)/k, g_N=sin(k*x).
    ! A constant transverse external field keeps nu smooth at internal zeros.
    ! The internal QUMOND force is then (nu*sin(k*x),0,0), with zero mean.
    ! Sources come from the analytic Laplacian, not composite_apply(exact).
    real(dp),allocatable::exact(:),rhs(:),pn(:),pq(:),ph(:),grad(:,:),gx(:),fq(:)
    real(dp)::k,resn,resq,en,egn,egq,external(3),removed
    integer::status,ierr,itn,itq
    character(len=8)::env
    env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_CONTINUUM',env,status=status)
    if(status/=0.or.trim(env)/='1')return
    associate(mesh=>qh_composite_mesh)
    allocate(exact(mesh%n),rhs(mesh%n),pn(mesh%n),pq(mesh%n),ph(mesh%n),grad(3,mesh%n),gx(mesh%n),fq(mesh%n))
    k=2d0*acos(-1d0)/mesh%box
    exact=cos(k*mesh%center(1,:))/k
    rhs=-k*k*exact
    removed=sum(rhs*mesh%volume)/sum(mesh%volume)
    rhs=rhs-removed
    exact=exact-sum(exact*mesh%volume)/sum(mesh%volume)
    pn=0d0
    call composite_solve(mesh,rhs,pn,1d-9,1200,40,itn,resn,ierr)
    if(ierr/=0)call fail('continuum Newtonian solve')
    en=sqrt(sum(mesh%volume*(pn-exact)**2)/sum(mesh%volume*exact**2))
    gx=sin(k*mesh%center(1,:))
    call composite_slopes(mesh,pn,grad,ierr)
    if(ierr/=0)call fail('continuum Newtonian force')
    grad(1,:)=grad(1,:)+gx
    egn=sqrt(sum(mesh%volume*sum(grad**2,dim=1))/sum(mesh%volume*gx**2))
    external=[0d0,1d0,0d0]
    call composite_phantom(mesh,pn,1d0,4d0*acos(-1d0),external,ph,ierr)
    if(ierr/=0)call fail('continuum phantom source')
    rhs=rhs+4d0*acos(-1d0)*ph
    rhs=rhs-sum(rhs*mesh%volume)/sum(mesh%volume)
    pq=pn
    call composite_solve(mesh,rhs,pq,1d-9,1200,40,itq,resq,ierr)
    if(ierr/=0)call fail('continuum QUMOND solve')
    fq=0.5d0*(1d0+sqrt(1d0+4d0/sqrt(1d0+gx**2)))*gx
    call composite_slopes(mesh,pq,grad,ierr)
    if(ierr/=0)call fail('continuum QUMOND force')
    grad(1,:)=grad(1,:)+fq
    egq=sqrt(sum(mesh%volume*sum(grad**2,dim=1))/sum(mesh%volume*fq**2))
    write(*,*) 'RAR_QUMOND_CONTINUUM_V1',mesh%n,itn,itq,resn,resq,en,egn,egq,removed
    end associate
  end subroutine continuum_check

  subroutine qumond_evolution_begin()
    real(dp)::dt,amax,old_dt
    integer::lev,status,ios
    character(len=80)::env,filename
    if(.not.qh_evolve.or.evolution_open)call fail('invalid evolution begin phase')
#ifdef PHASE0_STELLAR_ENRICHMENT
    if(qh_stellar_sources.and..not.use_channel_resolved_feedback()) &
         call fail('native stellar transaction requires channel-resolved feedback')
#endif
    call restrict_gas_and_ghosts()
    call qumond_publish_fixed_source()
    call refresh_gas_and_gravity_ghosts()
    old_dt=dtnew(levelmin)
    if(any(dtnew(levelmin:nlevelmax)/=old_dt))call fail('previous hierarchy timestep differs')
    rho_max(levelmin:nlevelmax)=max(physical_rho_max,tiny(1d0))
    do lev=levelmin,nlevelmax
       call newdt_fine(lev)
    enddo
    dt=minval(dtnew(levelmin:nlevelmax))
    amax=maxval(sqrt(sum(saved_force(qh_composite_mesh%cell,:)**2,dim=2)))
    if(amax>0d0)dt=min(dt,courant_factor*sqrt(minval(qh_composite_mesh%width)/amax))
    env=''
    call get_environment_variable('RAR_QUMOND_EVOLVE_DT',env,status=status)
    if(status==0.and.len_trim(env)>0)then
       read(env,*,iostat=ios)old_dt
       if(ios/=0)call fail('invalid evolution timestep cap')
       if(.not.ieee_is_finite(old_dt).or.old_dt<=0d0)call fail('invalid evolution timestep cap')
       dt=min(dt,old_dt)
    endif
    if(pic)call native_particles_begin(dt)
    if(.not.ieee_is_finite(dt).or.dt<=0d0.or.t+dt<=t)call fail('invalid common evolution timestep')
    dtnew(levelmin:nlevelmax)=dt
    evolution_start=t;evolution_dt=dt;evolution_epoch=nstep_coarse;evolution_open=.true.
    env=''
    call get_environment_variable('RAR_QUMOND_GRAVITY_KICKS',env,status=status)
    if(status==0.and.trim(env)=='1')then
       write(filename,'("gravity_kicks_",I5.5,".tsv")')nstep_coarse
       open(newunit=qh_kick_unit,file=trim(filename),status='new',action='write',iostat=ios)
       if(ios/=0)call fail('gravity kick diagnostic exists or cannot be created')
       write(qh_kick_unit,'(A)') '# phase cell level rho_source rho_receiver halfdt force[3] oldp[3] oldE newp[3] newE'
    endif
    write(*,*) 'RAR_QUMOND_EVOLVE_BEGIN_V1',nstep_coarse,t,dt,saved_time,qh_composite_mesh%n
  end subroutine qumond_evolution_begin

  subroutine qumond_evolution_end()
    integer::i,j,n,ids(nvector),lev
    real(dp)::oldp(nvector,3),olde(nvector),clock_error
    if(.not.evolution_open.or.evolution_epoch/=nstep_coarse)call fail('invalid evolution end phase')
    clock_error=abs(t-(evolution_start+evolution_dt))
    if(clock_error>16d0*spacing(max(abs(t),abs(evolution_start),evolution_dt))) &
         call fail('hierarchy clock did not advance exactly once')
    if(any(dtnew(levelmin:nlevelmax)/=evolution_dt))call fail('recursive hierarchy timestep changed')
    if(pic)call native_particles_drift(evolution_dt)
    call restrict_gas_and_ghosts()
    call qumond_history_enter(1,levelmin)
    call qumond_publish_fixed_source()
    if(pic)call native_particles_finish(evolution_dt)
    associate(mesh=>qh_composite_mesh)
    do i=1,mesh%n,nvector
       n=min(nvector,mesh%n-i+1);ids=0;ids(1:n)=mesh%cell(i:i+n-1)
       oldp(1:n,:)=uold(ids(1:n),2:4);olde(1:n)=uold(ids(1:n),5)
       call synchydrofine1(ids,n,0.5d0*evolution_dt)
       if(qh_kick_unit/=0)then
          do j=1,n
             lev=mesh%level(i+j-1)
             call qumond_record_kick(1,ids(j),lev,uold(ids(j),1),uold(ids(j),1), &
                  0.5d0*evolution_dt,f(ids(j),:),oldp(j,:),olde(j),uold(ids(j),2:4),uold(ids(j),5))
          enddo
       endif
    enddo
    if(qh_kick_unit/=0)then
       close(qh_kick_unit);qh_kick_unit=0
    endif
    if(.not.all(ieee_is_finite(uold(mesh%cell,:))))call fail('nonfinite gas after endpoint kick')
    write(*,*) 'RAR_QUMOND_EVOLVE_END_V1',nstep_coarse,evolution_start,t,evolution_dt, &
         sum(uold(mesh%cell,1)*mesh%volume)
    end associate
    call restrict_gas_and_ghosts()
    evolution_open=.false.
  end subroutine qumond_evolution_end

  subroutine qumond_stellar_endpoint()
    integer :: lev,ivar,i
    if(.not.qh_stellar_sources.or.evolution_open.or.saved_time/=t.or.published_generation/=source_generation) &
         call fail('stellar source transaction requires a synchronized completed endpoint')
    call native_stellar_begin()
    ! The old roster has already completed KDK. Repartition for the existing
    ! level-local source routines; no newborn belongs to the completed drift.
    do lev=levelmin,nlevelmax
       call kill_tree_fine(lev)
    enddo
    call diag_rar_global_inventory('before_native_return')
    ! Initialize EVERY level before any deposit; a later set_unew must not
    ! overwrite a deposit made near a refinement boundary.
    do lev=levelmin,nlevelmax
       call set_unew(lev)
    enddo
    do lev=levelmin,nlevelmax
       call thermal_feedback(lev)
    enddo
    do lev=nlevelmax,levelmin,-1
       do ivar=1,nvar
          call make_virtual_reverse_dp(unew(1,ivar),lev)
       enddo
    enddo
    associate(mesh=>qh_composite_mesh)
    do lev=levelmin,nlevelmax
       do i=1,size(qh_cell_level)
          if(qh_cell_level(i)/=lev)cycle
          if(son(i)>0.and.any(unew(i,:)/=uold(i,:))) &
               call fail('stellar source deposited into a covered cell')
       enddo
    enddo
    if(.not.all(ieee_is_finite(unew(mesh%cell,:))).or.any(unew(mesh%cell,1)<smallr)) &
         call fail('invalid stellar receiver state')
    ! set_uold also adds a gravity half-source and must not be used here.
    uold(mesh%cell,:)=unew(mesh%cell,:)
    call restrict_gas_and_ghosts()
    call diag_rar_global_inventory('after_native_return')
    if(cooling)then
       do lev=nlevelmax,levelmin,-1
          call diag_rar_cooling_inventory('before_cooling',lev)
          call cooling_fine(lev)
          call diag_rar_cooling_inventory('after_cooling',lev)
       enddo
    endif
    call restrict_gas_and_ghosts()
    call diag_rar_global_inventory('before_native_birth')
    do lev=nlevelmax,levelmin,-1
       call star_formation(lev)
       call restrict_gas_and_ghosts()
    enddo
    call diag_rar_global_inventory('after_native_birth')
    do lev=nlevelmax,levelmin,-1
       call merge_tree_fine(lev)
    enddo
    call native_stellar_reconcile(mesh,use_mond)
    if(source_generation==huge(source_generation))call fail('physical source generation overflow')
    source_generation=source_generation+1
    call qumond_publish_fixed_source()
    ! Refresh acceleration for the next interval only. There is no third
    ! gravity kick, newborn drift, or rotation of distinct-time histories.
    call refresh_gas_and_gravity_ghosts()
    write(*,*) 'RAR_NATIVE_STELLAR_ENDPOINT_V1',nstep_coarse,t,source_generation, &
         sum(uold(mesh%cell,1)*mesh%volume)
    end associate
  end subroutine qumond_stellar_endpoint

  subroutine restrict_gas_and_ghosts()
    integer::lev
    do lev=nlevelmax-1,1,-1
       call upload_fine(lev)
    enddo
    call refresh_gas_and_gravity_ghosts()
  end subroutine restrict_gas_and_ghosts

  subroutine refresh_gas_and_gravity_ghosts()
    integer::lev,v,d
    do lev=1,nlevelmax
       if(numbtot(1,lev)==0)cycle
       do v=1,nvar
          call make_virtual_fine_dp(uold(1,v),lev)
       enddo
       call make_virtual_fine_dp(phi(1),lev)
       do d=1,3
          call make_virtual_fine_dp(f(1,d),lev)
       enddo
    enddo
  end subroutine refresh_gas_and_gravity_ghosts

  subroutine qumond_restore_published()
    if(.not.allocated(saved_phi).or.published_epoch/=nstep_coarse)call fail('stale publication')
    phi=saved_phi;f=saved_force;rho=saved_rho
  end subroutine

  subroutine fail(reason)
    character(len=*),intent(in)::reason
    write(*,*) 'QUMOND composite publication failed: ',trim(reason)
    call clean_stop
  end subroutine
end module qumond_composite_runtime

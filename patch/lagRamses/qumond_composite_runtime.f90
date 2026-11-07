! First publication path: frozen gas/analytic source, serial periodic fixed mesh.
! The caller rejects particles, evolving hydro/subgrid physics and other gravity.
module qumond_composite_runtime
  use amr_commons
  use poisson_commons
  use hydro_commons, only: uold,smallr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: qumond_publish_fixed_source,qumond_restore_published
  real(dp),allocatable :: saved_phi(:),saved_force(:,:),saved_rho(:)
  integer :: published_epoch=-1
contains
  subroutine qumond_publish_fixed_source()
    real(dp),allocatable :: density(:),rhs(:),pn(:),pq(:),phantom(:),gradient(:,:),draft_n(:)
    real(dp)::xx(nvector,3),dd(nvector),mean_density,mean_source,fourpi,a0,external(3)
    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,res_n,res_q,balance
    integer::i,j,k,lev,ids(nvector),ierr,it_n,it_q,cell,export_status,unit
    character(len=8)::export_env
    if(.not.qh_publish)call clean_stop
    if(qh_composite_mesh%n==0)then
       call composite_native_build(qh_composite_mesh,ierr)
       if(ierr/=0)call fail('native geometry')
    endif
    if(nstep_coarse==0)call continuum_check()
    associate(mesh=>qh_composite_mesh)
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
    if(.not.all(ieee_is_finite(density)))call fail('source is nonfinite')
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
    ! No force consumer has run in this explicitly frozen-source path. Update
    ! stage histories and publish the complete hierarchy only after both solves.
    do lev=levelmin,nlevelmax
       do i=1,size(phi)
          if(qh_cell_level(i)/=lev)cycle
          if(qh_valid(1,lev))then
             qh_previous(i,1)=qh_current(i,1)
          else
             qh_previous(i,1)=draft_n(i)
          endif
          if(qh_valid(2,lev))then
             qh_previous(i,2)=qh_current(i,2)
          else
             qh_previous(i,2)=saved_phi(i)
          endif
          qh_current(i,1)=draft_n(i);qh_current(i,2)=saved_phi(i)
          phi_old(i)=qh_previous(i,2)
       enddo
       do k=1,2
          if(qh_valid(k,lev))then
             qh_old_time(k,lev)=qh_time(k,lev)
          else
             qh_old_time(k,lev)=t
          endif
          write(*,*) 'RAR_QUMOND_HISTORY_V1',k,lev,t,qh_old_time(k,lev),.not.qh_valid(k,lev)
          qh_time(k,lev)=t;qh_valid(k,lev)=.true.;qh_epoch(k,lev)=nstep_coarse
       enddo
    enddo
    published_epoch=nstep_coarse
    call qumond_restore_published()
    write(*,*) 'RAR_COMPOSITE_MODEL_V1',use_mond,gravity_type
    write(*,*) 'RAR_QUMOND_PUBLICATION_V1',nstep_coarse,t,mesh%n,it_n,it_q,res_n,res_q,balance, &
         sum(saved_phi(mesh%cell)*mesh%volume),sqrt(sum(mesh%volume*sum(gradient**2,dim=1)))
    export_env=''
    call get_environment_variable('RAR_QUMOND_COMPOSITE_EXPORT',export_env,status=export_status)
    if(export_status==0.and.trim(export_env)=='1'.and.nstep_coarse==0)then
       open(newunit=unit,file='composite_leaf_state.tsv',status='new',action='write',iostat=ierr)
       if(ierr/=0)call fail('leaf export exists or cannot be created')
       write(unit,'(A,5(1X,ES24.16))') '# constants',fourpi,a0,external
       do i=1,mesh%n
          write(unit,'(2(I0,1X),11(ES24.16,1X))') mesh%cell(i),mesh%level(i),mesh%center(:,i), &
               mesh%width(i),density(i),pn(i),pq(i),phantom(i),-gradient(:,i)
       enddo
       close(unit)
       open(newunit=unit,file='composite_faces.tsv',status='new',action='write',iostat=ierr)
       if(ierr/=0)call fail('face export exists or cannot be created')
       do i=1,mesh%nface
          write(unit,'(3(I0,1X),2(ES24.16,1X))') mesh%left(i),mesh%right(i),mesh%axis(i), &
               mesh%area(i),mesh%distance(i)
       enddo
       close(unit)
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

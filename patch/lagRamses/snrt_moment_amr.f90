! Conservative M_N transport on the actual RAMSES leaf/halo maps.
! Persistent state remains compact; angular work is one small tile.
module snrt_moment_amr
  use amr_commons
  use snrt_moment_transport
  use snrt_moment_dispatch
  use snrt_amr_topology
  use snrt_state, only: snrt_state_get_slot,snrt_state_sync_level
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifndef WITHOUTMPI
  use mpi_mod
#endif
#include "amr_index.h"
  implicit none
  private
  type,public :: mn_amr_mesh
     integer :: level=0,nleaf=0,nowned=0,nall=0,has_coarse=0,nfield=0
     integer,allocatable :: cells(:),slots(:),kind(:,:),neighbor(:,:),owner(:,:)
     real(dp) :: dx=0
     real(dp),allocatable :: width(:),volume(:)
  end type
  public :: mn_amr_prepare,mn_amr_advance,mn_collective_status
contains
  subroutine mn_collective_status(status)
    integer,intent(inout) :: status
#ifndef WITHOUTMPI
    integer :: global,info,ignored
    call MPI_ALLREDUCE(status,global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
    if(info/=0)call MPI_ABORT(MPI_COMM_WORLD,info,ignored)
    status=global
#endif
  end subroutine

  subroutine mn_amr_prepare(lev,mesh,ierr)
    integer,intent(in) :: lev
    type(mn_amr_mesh),intent(out) :: mesh
    integer,intent(out) :: ierr
    integer,allocatable :: cells(:),slots(:),neighbors(:,:),map(:),allcells(:)
    integer :: interfaces,nc,i,f,j,n,grid,child,cell,extra,unused,newslots,owner
    mesh%level=lev
    mesh%dx=boxlen/dble(icoarse_max-icoarse_min+1)*.5d0**lev
    call snrt_amr_build_same_level_neighbors(lev,cells,slots,neighbors,nc,interfaces)
    mesh%nleaf=nc;ierr=0
    allocate(mesh%kind(6,nc),mesh%neighbor(6,nc),mesh%owner(6,nc))
    if(nc>0)mesh%kind=snrt_face_kind
    if(nc>0)then
       if(any(mesh%kind==SNRT_FACE_FINE_TO_COARSE))mesh%has_coarse=1
       if(any(mesh%kind==SNRT_FACE_UNMAPPED))then
          write(*,*)'M_N unmapped faces rank/count=',myid,count(mesh%kind==SNRT_FACE_UNMAPPED)
          ierr=1
       endif
    endif
    call mn_collective_status(mesh%has_coarse)
    if(mesh%has_coarse/=0.and.lev<=1)ierr=1
    call mn_collective_status(ierr)
    if(ierr/=0)return
    extra=0
    if(mesh%has_coarse/=0)then
       call snrt_state_sync_level(lev-1,unused,newslots)
       grid=headl(myid,lev-1)
       do while(grid>0)
          do child=1,twotondim
             if(son(ICELL_OF(grid,child))==0)extra=extra+1
          enddo
          grid=next(grid)
       enddo
    endif
    mesh%nfield=ICELL_OF(ngridmax,twotondim)
    allocate(map(mesh%nfield),allcells(nc+extra+6*nc));map=0
    n=nc;allcells(1:nc)=cells
    do i=1,nc
       map(cells(i))=i
    enddo
    if(mesh%has_coarse/=0)then
       grid=headl(myid,lev-1)
       do while(grid>0)
          do child=1,twotondim
             cell=ICELL_OF(grid,child)
             if(son(cell)/=0)cycle
             n=n+1;allcells(n)=cell;map(cell)=n
          enddo
          grid=next(grid)
       enddo
    endif
    mesh%nowned=n;mesh%neighbor=0;mesh%owner=0
    do i=1,nc
       do f=1,6
          select case(mesh%kind(f,i))
          case(SNRT_FACE_LOCAL,SNRT_FACE_MPI,SNRT_FACE_FINE_TO_COARSE)
             cell=snrt_face_cell(f,i)
             ! Local faces are encoded by the leaf-index stencil, whereas
             ! snrt_face_cell carries only interface/remote cell identities.
             if(mesh%kind(f,i)==SNRT_FACE_LOCAL)cell=cells(neighbors(f,i))
             if(cell<1.or.cell>mesh%nfield)then
                ierr=1;cycle
             endif
             owner=snrt_cell_grid_owner(cell)
             ! Only owned entries precede the appended ghost entries.
             if(map(cell)>0.and.map(cell)<=mesh%nowned)owner=myid
             if(owner<1.or.owner>ncpu)then
                write(*,*)'M_N missing cell owner rank/cell/kind=',myid,cell,mesh%kind(f,i)
                ierr=1;cycle
             endif
             mesh%owner(f,i)=owner
             j=map(cell)
             if(j==0.and.owner==myid)then
                ierr=1;cycle
             endif
             if(j==0)then
                n=n+1;j=n;allcells(n)=cell;map(cell)=n
             endif
             mesh%neighbor(f,i)=j
          case(SNRT_FACE_PHYSICAL,SNRT_FACE_COARSE_TO_FINE)
          case default
             ierr=1
          end select
       enddo
    enddo
    call mn_collective_status(ierr)
    if(ierr/=0)return
    mesh%nall=n;mesh%cells=allcells(1:n)
    allocate(mesh%slots(mesh%nowned),mesh%width(n),mesh%volume(mesh%nowned))
    mesh%width=mesh%dx
    mesh%width(nc+1:mesh%nowned)=2*mesh%dx
    do i=1,nc
       do f=1,6
          if(mesh%kind(f,i)==SNRT_FACE_FINE_TO_COARSE) &
               mesh%width(mesh%neighbor(f,i))=2*mesh%dx
       enddo
    enddo
    mesh%volume=mesh%width(1:mesh%nowned)**3
    do i=1,mesh%nowned
       mesh%slots(i)=snrt_state_get_slot(mesh%cells(i))
       if(mesh%slots(i)<1)ierr=1
    enddo
    call mn_collective_status(ierr)
  end subroutine

  subroutine mn_amr_advance(b,mesh,u,cdt,escape,projection,ierr)
    type(mn_basis),intent(in) :: b
    type(mn_amr_mesh),intent(in) :: mesh
    real(dp),intent(inout) :: u(:,:)
    real(dp),intent(in) :: cdt
    real(dp),intent(out) :: escape(:),projection(:)
    integer,intent(out) :: ierr
    real(dp),allocatable :: base(:,:),stage(:,:),candidate(:,:),second_base(:,:),warm(:,:),final_cache(:,:)
    real(dp) :: e1(b%nm),e2(b%nm),p1(b%nm),p2(b%nm),check(b%nq),residual
    integer :: i,active_group,diagnostic(2)
    ierr=mn_bad_input;escape=0;projection=0
    if(any(shape(u)/=[b%nm,mesh%nowned]).or.size(escape)/=b%nm.or.size(projection)/=b%nm)return
    if(.not.ieee_is_finite(cdt).or.cdt<0.or.cdt/mesh%dx>1d0/12+epsilon(1d0))return
    ierr=0
    active_group=0
    if(any(u/=0))active_group=1
    call mn_collective_status(active_group)
    if(active_group==0.or.cdt==0)return
    allocate(base(b%nm,mesh%nowned),stage(b%nm,mesh%nowned),candidate(b%nm,mesh%nowned),second_base(b%nm,mesh%nowned))
    allocate(warm(b%nm-1,mesh%nall));warm=0
    call euler(u,stage,base,e1,p1,ierr)
    if(ierr/=0)return
    call euler(stage,candidate,second_base,e2,p2,ierr)
    if(ierr/=0)return
    candidate=.5d0*(base+candidate)
    allocate(final_cache(b%nm+1,mesh%nowned))
    call mn_dispatch_closure(b,candidate,stage,final_cache,warm(:,1:mesh%nowned),ierr)
    call mn_collective_status(ierr)
    if(ierr/=0)return
    u=candidate;escape=.5d0*(e1+e2);projection=p1+.5d0*p2
  contains
    subroutine exchange(values,reverse,status)
      real(dp),intent(inout) :: values(:,:)
      logical,intent(in) :: reverse
      integer,intent(out) :: status
      status=0
      if(ncpu<=1)return
      call snrt_halo_tile_exchange(values,mesh%level,status,reverse)
      if(status/=0)return
      if(mesh%has_coarse/=0)call snrt_halo_tile_exchange(values,mesh%level-1,status,reverse)
    end subroutine

    subroutine euler(input,output,projected_base,boundary,receipt,status)
      real(dp),intent(in) :: input(:,:)
      real(dp),intent(out) :: output(:,:),projected_base(:,:),boundary(:),receipt(:)
      integer,intent(out) :: status
      integer,parameter :: angle_tile=8,field_tile=32
      real(dp),allocatable :: field(:,:),extended(:,:),cache(:,:),theta(:),samples(:),slopes(:,:)
      real(dp),allocatable :: donor(:,:,:),delta(:,:)
      real(dp),allocatable :: all_projected(:,:),left(:,:,:),right(:,:,:),geometry(:,:),face_flux(:,:)
      integer,allocatable :: faces(:,:)
      real(dp) :: angular(b%nq),moments(b%nm),flux(b%nm),mu,il,ir,signface,flow,dl,dr,s
      integer :: first,last,col,j,i,q,f,a,r,k,nq,component,ncol,nface,ff,fl,fi,fs
      allocate(field(mesh%nfield,field_tile),extended(b%nm,mesh%nall),cache(b%nm+1,mesh%nall))
      allocate(theta(mesh%nowned),samples(mesh%nall),slopes(3,mesh%nowned))
      allocate(donor(4,angle_tile,mesh%nall),delta(b%nm,mesh%nall))
      allocate(all_projected(b%nm,mesh%nall),faces(3,6*mesh%nleaf))
      allocate(left(4,angle_tile,256),right(4,angle_tile,256),geometry(5,256),face_flux(b%nm,256))
      nface=0
      do i=1,mesh%nleaf
         do f=1,6
            r=mesh%neighbor(f,i)
            select case(mesh%kind(f,i))
            case(SNRT_FACE_COARSE_TO_FINE)
               cycle
            case(SNRT_FACE_LOCAL)
               if(i>r)cycle
            case(SNRT_FACE_MPI)
               if(myid>mesh%owner(f,i))cycle
            end select
            nface=nface+1;faces(:,nface)=[i,r,f]
         enddo
      enddo
      extended=0
      do first=1,b%nm,field_tile
         last=min(b%nm,first+field_tile-1);ncol=last-first+1;field=0
         do i=1,mesh%nowned
            field(mesh%cells(i),1:ncol)=input(first:last,i)
         enddo
         call exchange(field(:,1:ncol),.false.,status)
         if(status/=0)return
         do i=1,mesh%nall
            extended(first:last,i)=field(mesh%cells(i),1:ncol)
         enddo
      enddo
      status=0;receipt=0;boundary=0;theta=1;delta=0
      call mn_dispatch_closure(b,extended,all_projected,cache,warm,status)
      call mn_collective_status(status)
      if(status/=0)return
      projected_base=all_projected(:,1:mesh%nowned)
      do i=1,mesh%nowned
         receipt=receipt+mesh%volume(i)*(projected_base(:,i)-input(:,i))
      enddo
      ! The same theta applies to all quadrature nodes of a cell. Coarse
      ! donors use conservative injection, as in the existing AMR interface.
      do q=1,b%nq
         call sample_slopes(b,mesh,extended,cache,q,samples,slopes)
         do i=1,mesh%nleaf
            do a=1,3
               s=abs(slopes(a,i))*mesh%dx
               if(s>0)theta(i)=min(theta(i),samples(i)/s)
            enddo
         enddo
      enddo
      do first=1,b%nq,angle_tile
         nq=min(angle_tile,b%nq-first+1);donor=0;field=0
         do k=1,nq
            q=first+k-1
            call sample_slopes(b,mesh,extended,cache,q,samples,slopes)
            do i=1,mesh%nowned
               donor(1,k,i)=samples(i)
               donor(2:4,k,i)=theta(i)*slopes(:,i)
               field(mesh%cells(i),4*(k-1)+1:4*k)=donor(:,k,i)
            enddo
         enddo
         call exchange(field(:,1:4*nq),.false.,status)
         if(status/=0)return
         do i=mesh%nowned+1,mesh%nall
            do k=1,nq
               donor(:,k,i)=field(mesh%cells(i),4*(k-1)+1:4*k)
            enddo
         enddo
         do fs=1,nface,256
            fl=min(nface,fs+255);right=0;geometry=0
            do fi=1,fl-fs+1
               ff=fs+fi-1;i=faces(1,ff);r=faces(2,ff);f=faces(3,ff)
               a=(f+1)/2;signface=1
               if(mod(f,2)==1)signface=-1
               left(:,1:nq,fi)=donor(:,1:nq,i)
               geometry(:,fi)=[real(a,dp),signface,mesh%width(i),0d0,cdt*mesh%dx**2]
               if(r>0)then
                  right(:,1:nq,fi)=donor(:,1:nq,r);geometry(4,fi)=mesh%width(r)
               endif
            enddo
            call mn_dispatch_flux(b%harmonic(:,first:first+nq-1),b%weight(first:first+nq-1), &
                 b%direction(:,first:first+nq-1),left(:,1:nq,1:fl-fs+1),right(:,1:nq,1:fl-fs+1), &
                 geometry(:,1:fl-fs+1),face_flux(:,1:fl-fs+1),status)
            if(status/=0)exit
            do fi=1,fl-fs+1
               ff=fs+fi-1;i=faces(1,ff);r=faces(2,ff)
               delta(:,i)=delta(:,i)-face_flux(:,fi)
               if(r>0)then
                  delta(:,r)=delta(:,r)+face_flux(:,fi)
               else
                  boundary=boundary+face_flux(:,fi)
               endif
            enddo
         enddo
         call mn_collective_status(status)
         if(status/=0)return
      enddo
      ! Reverse SUM sends the face owner's already integrated flux to the
      ! remote cell owner, for same-level and coarse receivers separately.
      do first=1,b%nm,field_tile
         last=min(b%nm,first+field_tile-1);ncol=last-first+1;field=0
         do i=1,mesh%nall
            field(mesh%cells(i),1:ncol)=delta(first:last,i)
         enddo
         call exchange(field(:,1:ncol),.true.,status)
         if(status/=0)return
         do i=1,mesh%nowned
            output(first:last,i)=projected_base(first:last,i)+field(mesh%cells(i),1:ncol)/mesh%volume(i)
         enddo
      enddo
      status=0
      if(any(.not.ieee_is_finite(output)))status=mn_bad_step
      call mn_collective_status(status)
    end subroutine
  end subroutine

      subroutine sample_slopes(b,mesh,extended,cache,q,samples,slopes)
        type(mn_basis),intent(in) :: b
        type(mn_amr_mesh),intent(in) :: mesh
        real(dp),intent(in) :: extended(:,:),cache(:,:)
        real(dp),intent(out) :: samples(:),slopes(:,:)
        integer,intent(in) :: q
        integer :: i,a,l,r
        real(dp) :: low,high
        samples=0;slopes=0
        do i=1,mesh%nall
           if(extended(1,i)>0)samples(i)=extended(1,i)* &
                exp(dot_product(cache(1:b%nm-1,i),b%harmonic(2:,q)- &
                extended(2:,i)/extended(1,i))-cache(b%nm,i))/cache(b%nm+1,i)
        enddo
        do i=1,mesh%nleaf
           do a=1,3
              l=mesh%neighbor(2*a-1,i);r=mesh%neighbor(2*a,i)
              if(l==0.or.r==0)cycle
              low=(samples(i)-samples(l))/(.5d0*(mesh%width(i)+mesh%width(l)))
              high=(samples(r)-samples(i))/(.5d0*(mesh%width(i)+mesh%width(r)))
              slopes(a,i)=mn_limited_slope(low,high,b%spatial_limiter)
           enddo
        enddo
      end subroutine
end module

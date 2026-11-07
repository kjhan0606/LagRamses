! Native periodic leaf/subface geometry and conservative scalar operators.
! No gravity publication or composite elliptic convergence is implied here.
module qumond_composite
  use amr_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: composite_mesh_t, composite_native_build, composite_apply, composite_phantom
  public :: composite_slopes, composite_face_gradient, composite_solve
  type composite_mesh_t
     integer :: n=0,nface=0
     real(dp) :: box=0d0
     integer,allocatable :: cell(:),level(:),left(:),right(:),axis(:)
     real(dp),allocatable :: center(:,:),width(:),volume(:),area(:),distance(:),inverse(:,:,:)
  end type
contains
  integer function locate_leaf(point) result(cell)
    real(dp),intent(in)::point(3)
    real(dp)::q(3)
    integer::grid,lev,child,d
    q=modulo(point,boxlen)/boxlen
    grid=son(1);cell=0
    if(grid<=0)return
    do lev=1,nlevelmax
       child=1
       do d=1,3
          if(q(d)>=xg(grid,d))child=child+2**(d-1)
       enddo
       cell=ICELL_OF(grid,child)
       if(son(cell)==0)return
       grid=son(cell)
    enddo
    cell=0
  end function

  subroutine composite_native_build(mesh,ierr)
    type(composite_mesh_t),intent(out)::mesh
    integer,intent(out)::ierr
    integer,allocatable::lookup(:)
    real(dp),allocatable::coverage(:,:,:),matrix(:,:,:)
    integer::lev,g,grid,c,cell,i,j,a,b,d,sgn1,sgn2,tangent(2),k,f,l,r
    real(dp)::point(3),delta(3),w,det,mat(3,3),scale
    ierr=1
    if(ndim/=3.or.ncpu/=1.or.nx/=1.or.ny/=1.or.nz/=1.or.simple_boundary.or. &
         nlevelmax/=levelmin+1.or.levelmin<2.or.boxlen<=0d0)return
    mesh%box=boxlen
    mesh%n=0
    do lev=levelmin,nlevelmax
       do g=1,active(lev)%ngrid
          do c=1,twotondim
             cell=ICELL_OF(active(lev)%igrid(g),c)
             if(son(cell)==0)mesh%n=mesh%n+1
          enddo
       enddo
    enddo
    if(mesh%n==0)return
    allocate(mesh%cell(mesh%n),mesh%level(mesh%n),mesh%center(3,mesh%n))
    allocate(mesh%width(mesh%n),mesh%volume(mesh%n),lookup(size(son)))
    lookup=0;i=0
    do lev=levelmin,nlevelmax
       w=boxlen*0.5d0**lev
       do g=1,active(lev)%ngrid
          grid=active(lev)%igrid(g)
          do c=1,twotondim
             cell=ICELL_OF(grid,c)
             if(son(cell)/=0)cycle
             i=i+1;mesh%cell(i)=cell;mesh%level(i)=lev;lookup(cell)=i
             mesh%width(i)=w;mesh%volume(i)=w**3
             do d=1,3
                mesh%center(d,i)=xg(grid,d)*boxlen+ &
                     (real(ibits(c-1,d-1,1),dp)-0.5d0)*w
             enddo
          enddo
       enddo
    enddo
    if(abs(sum(mesh%volume)-boxlen**3)>1d-12*boxlen**3)return
    allocate(mesh%left(12*mesh%n),mesh%right(12*mesh%n),mesh%axis(12*mesh%n))
    allocate(mesh%area(12*mesh%n),mesh%distance(12*mesh%n))
    do i=1,mesh%n
       do a=1,3
          point=mesh%center(:,i);point(a)=point(a)+mesh%width(i)*(0.5d0+1d-7)
          cell=locate_leaf(point)
          if(cell<=0)return
          j=lookup(cell)
          if(j<=0)return
          if(mesh%width(j)<mesh%width(i))then
             if(mesh%width(i)/=2d0*mesh%width(j))return
             k=0
             do d=1,3
                if(d==a)cycle
                k=k+1;tangent(k)=d
             enddo
             do sgn1=-1,1,2
                do sgn2=-1,1,2
                   point=mesh%center(:,i);point(a)=point(a)+mesh%width(i)*(0.5d0+1d-7)
                   point(tangent(1))=point(tangent(1))+sgn1*mesh%width(i)*0.25d0
                   point(tangent(2))=point(tangent(2))+sgn2*mesh%width(i)*0.25d0
                   cell=locate_leaf(point)
                   if(cell<=0)return
                   j=lookup(cell)
                   if(j<=0)return
                   if(mesh%width(i)/=2d0*mesh%width(j))return
                   call append_face(i,j,a)
                enddo
             enddo
          else
             if(mesh%width(j)>2d0*mesh%width(i))return
             call append_face(i,j,a)
          endif
       enddo
    enddo
    allocate(coverage(3,2,mesh%n),matrix(3,3,mesh%n),mesh%inverse(3,3,mesh%n))
    coverage=0d0;matrix=0d0
    do f=1,mesh%nface
       l=mesh%left(f);r=mesh%right(f);a=mesh%axis(f)
       coverage(a,1,l)=coverage(a,1,l)+mesh%area(f)
       coverage(a,2,r)=coverage(a,2,r)+mesh%area(f)
       delta=mesh%center(:,r)-mesh%center(:,l)
       delta=delta-mesh%box*anint(delta/mesh%box)
       w=mesh%area(f)/mesh%distance(f)
       do a=1,3
          do b=1,3
             matrix(a,b,l)=matrix(a,b,l)+w*delta(a)*delta(b)
             matrix(a,b,r)=matrix(a,b,r)+w*delta(a)*delta(b)
          enddo
       enddo
    enddo
    do i=1,mesh%n
       if(any(abs(coverage(:,:,i)-mesh%width(i)**2)>1d-12*mesh%width(i)**2))return
       mat=matrix(:,:,i);scale=maxval(abs(mat))
       det=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
          -mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
          +mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))
       if(det<=1d-12*scale**3)return
       mesh%inverse(1,1,i)=(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))/det
       mesh%inverse(1,2,i)=(mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3))/det
       mesh%inverse(1,3,i)=(mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2))/det
       mesh%inverse(2,1,i)=(mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3))/det
       mesh%inverse(2,2,i)=(mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1))/det
       mesh%inverse(2,3,i)=(mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3))/det
       mesh%inverse(3,1,i)=(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))/det
       mesh%inverse(3,2,i)=(mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2))/det
       mesh%inverse(3,3,i)=(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))/det
    enddo
    ierr=0
  contains
    subroutine append_face(l,r,a)
      integer,intent(in)::l,r,a
      mesh%nface=mesh%nface+1
      f=mesh%nface
      mesh%left(f)=l;mesh%right(f)=r;mesh%axis(f)=a
      mesh%area(f)=min(mesh%width(l),mesh%width(r))**2
      mesh%distance(f)=0.5d0*(mesh%width(l)+mesh%width(r))
    end subroutine
  end subroutine

  subroutine composite_slopes(mesh,p,gradient,ierr)
    type(composite_mesh_t),intent(in)::mesh
    real(dp),intent(in)::p(:)
    real(dp),intent(out)::gradient(:,:)
    integer,intent(out)::ierr
    integer::f,l,r,i
    real(dp)::delta(3),term(3)
    ierr=1;gradient=0d0
    if(size(p)/=mesh%n.or.size(gradient,1)/=3.or.size(gradient,2)/=mesh%n)return
    if(.not.all(ieee_is_finite(p)))return
    do f=1,mesh%nface
       l=mesh%left(f);r=mesh%right(f)
       delta=mesh%center(:,r)-mesh%center(:,l)
       delta=delta-mesh%box*anint(delta/mesh%box)
       term=mesh%area(f)/mesh%distance(f)*delta*(p(r)-p(l))
       gradient(:,l)=gradient(:,l)+term;gradient(:,r)=gradient(:,r)+term
    enddo
    do i=1,mesh%n
       gradient(:,i)=matmul(mesh%inverse(:,:,i),gradient(:,i))
    enddo
    if(.not.all(ieee_is_finite(gradient)))return
    ierr=0
  end subroutine

  subroutine composite_face_gradient(mesh,p,gradient,f,vector)
    type(composite_mesh_t),intent(in)::mesh
    real(dp),intent(in)::p(:),gradient(:,:)
    integer,intent(in)::f
    real(dp),intent(out)::vector(3)
    integer::l,r,a
    real(dp)::pl,pr,delta(3)
    l=mesh%left(f);r=mesh%right(f);a=mesh%axis(f)
    pl=p(l);pr=p(r)
    delta=mesh%center(:,r)-mesh%center(:,l)
    delta=delta-mesh%box*anint(delta/mesh%box);delta(a)=0d0
    if(mesh%width(l)>mesh%width(r))pl=pl+dot_product(gradient(:,l),delta)
    if(mesh%width(r)>mesh%width(l))pr=pr-dot_product(gradient(:,r),delta)
    vector=0.5d0*(gradient(:,l)+gradient(:,r))
    vector(a)=(pr-pl)/mesh%distance(f)
  end subroutine

  subroutine composite_apply(mesh,p,laplacian,ierr)
    type(composite_mesh_t),intent(in)::mesh
    real(dp),intent(in)::p(:)
    real(dp),intent(out)::laplacian(:)
    integer,intent(out)::ierr
    real(dp),allocatable::gradient(:,:)
    real(dp)::vector(3),term
    integer::f,l,r
    ierr=1;laplacian=0d0
    if(size(laplacian)/=mesh%n)return
    allocate(gradient(3,mesh%n));call composite_slopes(mesh,p,gradient,ierr)
    if(ierr/=0)return
    do f=1,mesh%nface
       call composite_face_gradient(mesh,p,gradient,f,vector)
       term=mesh%area(f)*vector(mesh%axis(f))
       l=mesh%left(f);r=mesh%right(f)
       laplacian(l)=laplacian(l)+term;laplacian(r)=laplacian(r)-term
    enddo
    laplacian=laplacian/mesh%volume
    if(.not.all(ieee_is_finite(laplacian)))then
       ierr=2;return
    endif
    ierr=0
  end subroutine

  subroutine composite_phantom(mesh,p,a0,fourpi,external,rho,ierr)
    type(composite_mesh_t),intent(in)::mesh
    real(dp),intent(in)::p(:),a0,fourpi,external(3)
    real(dp),intent(out)::rho(:)
    integer,intent(out)::ierr
    real(dp),allocatable::gradient(:,:)
    real(dp)::vector(3),mag,boost,term
    integer::f,l,r
    ierr=1;rho=0d0
    if(size(rho)/=mesh%n.or.a0<=0d0.or.fourpi<=0d0)return
    if(.not.all(ieee_is_finite([a0,fourpi,external])))return
    allocate(gradient(3,mesh%n));call composite_slopes(mesh,p,gradient,ierr)
    if(ierr/=0)return
    do f=1,mesh%nface
       call composite_face_gradient(mesh,p,gradient,f,vector)
       vector=-vector+external;mag=sqrt(sum(vector**2))
       if(mag==0d0)cycle
       boost=2d0*(a0/mag)/(sqrt(1d0+4d0*(a0/mag))+1d0)
       term=mesh%area(f)*boost*vector(mesh%axis(f))/fourpi
       l=mesh%left(f);r=mesh%right(f)
       rho(l)=rho(l)-term;rho(r)=rho(r)+term
    enddo
    rho=rho/mesh%volume
    if(.not.all(ieee_is_finite(rho)))then
       ierr=2;return
    endif
    ierr=0
  end subroutine

  ! Restarted flexible GMRES for L phi = rhs on the volume-weighted zero-mean
  ! space. Face-Jacobi is an initial preconditioner, not a claimed AMR MG solve.
  ! Success requires the independently recomputed, unpreconditioned residual.
  subroutine composite_solve(mesh,rhs,x,rtol,maxiter,restart,iterations,relative_residual,ierr)
    type(composite_mesh_t),intent(in)::mesh
    real(dp),intent(in)::rhs(:),rtol
    real(dp),intent(inout)::x(:)
    integer,intent(in)::maxiter,restart
    integer,intent(out)::iterations,ierr
    real(dp),intent(out)::relative_residual
    real(dp),allocatable::v(:,:),z(:,:),h(:,:),cs(:),sn(:),g(:),y(:),b(:),r(:),w(:),diagonal(:)
    real(dp)::bnorm,beta,hn,entry,tmp,radius,volume_sum,mean_rhs
    integer::n,m,f,l,rr,i,j,k,pass,status
    logical::happy
    ierr=1;iterations=0;relative_residual=huge(1d0)
    n=mesh%n;m=min(restart,maxiter)
    if(n<=0.or.size(rhs)/=n.or.size(x)/=n.or.m<1.or.rtol<=0d0.or.rtol>=1d0)return
    if(.not.all(ieee_is_finite(rhs)).or..not.all(ieee_is_finite(x)))return
    if(.not.ieee_is_finite(rtol))return
    allocate(v(n,m+1),z(n,m),h(m+1,m),cs(m),sn(m),g(m+1),y(m),b(n),r(n),w(n),diagonal(n))
    volume_sum=sum(mesh%volume)
    b=-rhs;bnorm=weighted_norm(b)
    if(.not.ieee_is_finite(bnorm))return
    mean_rhs=sum(b*mesh%volume)/volume_sum
    if(abs(mean_rhs)*sqrt(volume_sum)>1d-12*bnorm)then
       ierr=4;return ! Incompatible periodic source must be fixed by caller.
    endif
    call project(b);bnorm=weighted_norm(b)
    if(bnorm==0d0)then
       x=0d0;relative_residual=0d0;ierr=0;return
    endif
    diagonal=0d0
    do f=1,mesh%nface
       l=mesh%left(f);rr=mesh%right(f)
       entry=mesh%area(f)/mesh%distance(f)
       diagonal(l)=diagonal(l)+entry;diagonal(rr)=diagonal(rr)+entry
    enddo
    diagonal=diagonal/mesh%volume
    if(any(diagonal<=0d0).or..not.all(ieee_is_finite(diagonal)))return
    call project(x)
    call composite_apply(mesh,x,w,status)
    if(status/=0)return
    r=b+w
    do
       beta=weighted_norm(r);relative_residual=beta/bnorm
       if(.not.ieee_is_finite(relative_residual))return
       if(relative_residual<=rtol)then
          ierr=0;return
       endif
       if(iterations>=maxiter)then
          ierr=2;return
       endif
       call project(r);beta=weighted_norm(r)
       if(beta<=tiny(1d0))then
          ierr=3;return
       endif
       h=0d0;g=0d0;g(1)=beta;v(:,1)=r/beta;k=0
       do j=1,min(m,maxiter-iterations)
          z(:,j)=v(:,j)/diagonal;call project(z(:,j))
          call composite_apply(mesh,z(:,j),w,status)
          if(status/=0)return
          w=-w;call project(w)
          ! Two modified Gram-Schmidt passes control loss of orthogonality.
          do pass=1,2
             do i=1,j
                entry=sum(mesh%volume*v(:,i)*w)
                h(i,j)=h(i,j)+entry;w=w-entry*v(:,i)
             enddo
          enddo
          hn=weighted_norm(w);h(j+1,j)=hn;happy=hn<=tiny(1d0)
          if(.not.happy)v(:,j+1)=w/hn
          do i=1,j-1
             tmp=cs(i)*h(i,j)+sn(i)*h(i+1,j)
             h(i+1,j)=-sn(i)*h(i,j)+cs(i)*h(i+1,j);h(i,j)=tmp
          enddo
          radius=sqrt(h(j,j)**2+h(j+1,j)**2)
          if(radius<=tiny(1d0).or..not.ieee_is_finite(radius))then
             ierr=3;return
          endif
          cs(j)=h(j,j)/radius;sn(j)=h(j+1,j)/radius
          h(j,j)=radius;h(j+1,j)=0d0
          g(j+1)=-sn(j)*g(j);g(j)=cs(j)*g(j)
          iterations=iterations+1;k=j
          if(abs(g(j+1))/bnorm<=rtol.or.happy)exit
       enddo
       y=0d0
       do i=k,1,-1
          if(abs(h(i,i))<=tiny(1d0))then
             ierr=3;return
          endif
          y(i)=(g(i)-dot_product(h(i,i+1:k),y(i+1:k)))/h(i,i)
       enddo
       x=x+matmul(z(:,1:k),y(1:k));call project(x)
       if(.not.all(ieee_is_finite(x)))return
       call composite_apply(mesh,x,w,status)
       if(status/=0)return
       r=b+w
    enddo
  contains
    subroutine project(q)
      real(dp),intent(inout)::q(:)
      q=q-sum(q*mesh%volume)/volume_sum
    end subroutine
    real(dp) function weighted_norm(q) result(norm)
      real(dp),intent(in)::q(:)
      norm=sqrt(sum(mesh%volume*q*q))
    end function
  end subroutine composite_solve

end module qumond_composite

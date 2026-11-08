! Shared volume-integrated physical cloud assignment on the native octree.
module qumond_composite_clouds
  use amr_commons
  use qumond_composite, only: composite_mesh_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#include "amr_index.h"
  implicit none
  private
  public :: native_cloud_map
contains
  subroutine fail(reason)
    character(len=*),intent(in) :: reason
    write(*,*) 'Native composite cloud failed: ',trim(reason)
    call clean_stop
  end subroutine

  subroutine native_cloud_map(mesh,positions,width,ptr,leaf,weights,max_error)
    type(composite_mesh_t),intent(in) :: mesh
    real(dp),intent(in) :: positions(:,:),width
    integer,allocatable,intent(out) :: ptr(:),leaf(:)
    real(dp),allocatable,intent(out) :: weights(:)
    real(dp),intent(out) :: max_error
    integer,allocatable :: lookup(:)
    integer :: p,k,sx,sy,sz,nweight
    real(dp) :: lo(3),hi(3),position(3),normalization
    if(size(positions,2)/=3.or..not.all(ieee_is_finite(positions)))call fail('invalid positions')
    if(.not.ieee_is_finite(width).or.width<=0d0.or.width>mesh%box)call fail('invalid physical width')
    allocate(ptr(size(positions,1)+1),leaf(max(1,64*size(positions,1))))
    allocate(weights(size(leaf)),lookup(size(son)))
    lookup=0;nweight=0;max_error=0d0
    do k=1,mesh%n
       lookup(mesh%cell(k))=k
    enddo
    do p=1,size(positions,1)
       ptr(p)=nweight+1
       position=modulo(positions(p,:),mesh%box)
       do sz=-1,1
          do sy=-1,1
             do sx=-1,1
                lo=max(0d0,position+mesh%box*real([sx,sy,sz],dp)-0.5d0*width)
                hi=min(mesh%box,position+mesh%box*real([sx,sy,sz],dp)+0.5d0*width)
                if(any(hi<=lo))cycle
                call visit(son(1),1,lo,hi)
             enddo
          enddo
       enddo
       if(nweight<ptr(p))call fail('cloud has no leaf recipients')
       normalization=sum(weights(ptr(p):nweight))
       max_error=max(max_error,abs(normalization-1d0))
       if(abs(normalization-1d0)>1d-12)call fail('cloud does not cover unit weight')
       weights(ptr(p):nweight)=weights(ptr(p):nweight)/normalization
    enddo
    ptr(size(positions,1)+1)=nweight+1
    leaf=leaf(:nweight);weights=weights(:nweight)
  contains
    recursive subroutine visit(grid,lev,low,high)
      integer,intent(in) :: grid,lev
      real(dp),intent(in) :: low(3),high(3)
      integer :: child,cell,d,index
      real(dp) :: center(3),dx,volume
      integer,allocatable :: new_leaf(:)
      real(dp),allocatable :: new_weights(:)
      if(grid<=0.or.lev>nlevelmax)call fail('invalid cloud traversal')
      dx=mesh%box*0.5d0**lev
      do child=1,twotondim
         do d=1,3
            center(d)=xg(grid,d)*mesh%box+(real(ibits(child-1,d-1,1),dp)-0.5d0)*dx
         enddo
         volume=product(max(0d0,min(high,center+0.5d0*dx)-max(low,center-0.5d0*dx)))
         if(volume==0d0)cycle
         cell=ICELL_OF(grid,child)
         if(son(cell)>0)then
            call visit(son(cell),lev+1,low,high)
         else
            index=lookup(cell)
            if(index<=0)call fail('cloud reached a leaf outside the native mesh')
            if(nweight==size(leaf))then
               allocate(new_leaf(2*size(leaf)),new_weights(2*size(weights)))
               new_leaf(1:nweight)=leaf;new_weights(1:nweight)=weights
               call move_alloc(new_leaf,leaf);call move_alloc(new_weights,weights)
            endif
            nweight=nweight+1;leaf(nweight)=index;weights(nweight)=volume/width**3
         endif
      enddo
    end subroutine
  end subroutine native_cloud_map
end module qumond_composite_clouds

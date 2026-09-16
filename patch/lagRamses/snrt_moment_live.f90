! Compact radiation payloads indexed by the EXISTING RAMSES SNRT slot map.
! No cell ownership, chemistry, MPI, or new radiation physics is owned here.
module snrt_moment_live
  use snrt_moment_transport
  use snrt_spectral_contract, only: snrt_band_enabled
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer,parameter :: slab_size=128
  type :: radiation_slab
     ! Fixed groups: N/E. Spectral bands: native positive endpoint weights
     ! A/B. Never recover a weak component by subtracting stored N and E.
     real(mn_dp),allocatable :: primary(:,:,:),ir(:,:,:) ! moment, group/channel, slot
  end type
  type(radiation_slab),allocatable,save :: slabs(:)
  type(mn_basis),public,save :: mn_live_basis
  integer,save :: primary_groups=0
  public :: mn_live_initialize,mn_live_reserve,mn_live_clear,mn_live_read,mn_live_write
  public :: mn_live_pack,mn_live_unpack,mn_live_validate,mn_live_ir_read,mn_live_ir_write
  public :: mn_live_ir_pack,mn_live_ir_unpack,mn_live_ir_expand
contains
  subroutine mn_live_initialize(order,ng,ierr)
    integer,intent(in) :: order,ng
    integer,intent(out) :: ierr
    ierr=mn_bad_input
    if(ng<1)return
    if(mn_live_basis%nm>0)then
       if(mn_live_basis%order==order.and.primary_groups==ng)ierr=mn_ok
       return
    endif
    call mn_initialize(mn_live_basis,order,ierr)
    if(ierr/=mn_ok)return
    call mn_configure_spatial(mn_live_basis,'mc15',ierr)
    if(ierr/=mn_ok)return
    primary_groups=ng
  end subroutine

  subroutine mn_live_reserve(required,ierr)
    integer,intent(in) :: required
    integer,intent(out) :: ierr
    type(radiation_slab),allocatable :: next(:)
    integer :: old,n,k,status
    ierr=mn_bad_input
    if(required<0.or.mn_live_basis%nm==0)return
    n=(required+slab_size-1)/slab_size;old=0
    if(allocated(slabs))old=size(slabs)
    if(n<=old)then
       ierr=mn_ok;return
    endif
    allocate(next(n),stat=status)
    if(status/=0)return
    do k=1,old
       call move_alloc(slabs(k)%primary,next(k)%primary)
       call move_alloc(slabs(k)%ir,next(k)%ir)
    enddo
    call move_alloc(next,slabs)
    ierr=mn_ok
  end subroutine

  subroutine location(slot,block,offset,ierr)
    integer,intent(in) :: slot
    integer,intent(out) :: block,offset,ierr
    ierr=mn_bad_input;block=0;offset=0
    if(slot<1.or..not.allocated(slabs))return
    block=(slot-1)/slab_size+1;offset=mod(slot-1,slab_size)+1
    if(block>size(slabs))return
    ierr=mn_ok
  end subroutine

  subroutine mn_live_clear(slot)
    integer,intent(in) :: slot
    integer :: block,offset,status
    call location(slot,block,offset,status)
    if(status/=mn_ok)return
    if(allocated(slabs(block)%primary))slabs(block)%primary(:,:,offset)=0
    if(allocated(slabs(block)%ir))slabs(block)%ir(:,:,offset)=0
  end subroutine

  subroutine mn_live_validate(number,energy,ierr)
    real(mn_dp),intent(in) :: number(:,:),energy(:,:)
    integer,intent(out) :: ierr
    real(mn_dp) :: angular(mn_live_basis%nq),hint(mn_live_basis%nm-1)
    integer :: g
    ierr=mn_bad_input
    if(any(shape(number)/=[mn_live_basis%nm,primary_groups]))return
    if(any(shape(energy)/=shape(number)).or.mn_live_basis%nm==0)return
    do g=1,primary_groups
       hint=0
       call mn_reconstruct(mn_live_basis,number(:,g),angular,ierr,dual_hint=hint)
       if(ierr/=mn_ok)return
       call mn_reconstruct(mn_live_basis,energy(:,g),angular,ierr,dual_hint=hint)
       ! A warm hint is only an accelerator; it cannot change admission.
       if(ierr/=mn_ok)call mn_reconstruct(mn_live_basis,energy(:,g),angular,ierr)
       if(ierr/=mn_ok)return
       if(snrt_band_enabled())cycle
       ierr=mn_bad_input
       if(number(1,g)==0.and.any(energy(:,g)/=0))return
       if(number(1,g)>0.and.energy(1,g)<=0)return
    enddo
    ierr=mn_ok
  end subroutine

  subroutine mn_live_read(slot,number,energy,ierr)
    integer,intent(in) :: slot
    real(mn_dp),intent(out) :: number(:,:),energy(:,:)
    integer,intent(out) :: ierr
    integer :: block,offset
    number=0;energy=0;ierr=mn_bad_input
    if(any(shape(number)/=[mn_live_basis%nm,primary_groups]).or.any(shape(energy)/=shape(number)))return
    call location(slot,block,offset,ierr)
    if(ierr/=mn_ok)return
    if(.not.allocated(slabs(block)%primary))return
    number=slabs(block)%primary(:,1:primary_groups,offset)
    energy=slabs(block)%primary(:,primary_groups+1:,offset)
  end subroutine

  subroutine mn_live_write(slot,number,energy,ierr,validate_only,reserve_only)
    integer,intent(in) :: slot
    real(mn_dp),intent(in) :: number(:,:),energy(:,:)
    integer,intent(out) :: ierr
    logical,optional,intent(in) :: validate_only,reserve_only
    integer :: block,offset,status
    call mn_live_validate(number,energy,ierr)
    if(ierr/=mn_ok)return
    if(present(validate_only))then
       if(validate_only)return
    endif
    call location(slot,block,offset,ierr)
    if(ierr/=mn_ok)return
    if(.not.allocated(slabs(block)%primary))then
       if(all(number==0).and.all(energy==0))return
       allocate(slabs(block)%primary(mn_live_basis%nm,2*primary_groups,slab_size),stat=status)
       if(status/=0)then
          ierr=mn_bad_input;return
       endif
       slabs(block)%primary=0
    endif
    if(present(reserve_only))then
       if(reserve_only)return
    endif
    slabs(block)%primary(:,1:primary_groups,offset)=number
    slabs(block)%primary(:,primary_groups+1:,offset)=energy
  end subroutine

  subroutine mn_live_pack(slot,stride,number,energy,ierr)
    integer,intent(in) :: slot,stride
    real(mn_dp),intent(out) :: number(:),energy(:)
    integer,intent(out) :: ierr
    real(mn_dp) :: n(mn_live_basis%nm,primary_groups),e(mn_live_basis%nm,primary_groups)
    real(mn_dp) :: padded(stride,primary_groups)
    ierr=mn_bad_input;number=0;energy=0
    if(stride<mn_live_basis%nm.or.size(number)/=stride*primary_groups.or.size(energy)/=size(number))return
    call mn_live_read(slot,n,e,ierr)
    if(ierr/=mn_ok)return
    ! Private payloads enter ONLY through a validated write (or exact zero
    ! initialization). Packing is an exact FP64 copy, not a new closure.
    ierr=mn_bad_input
    if(any(.not.ieee_is_finite(n)).or.any(.not.ieee_is_finite(e)))return
    padded=0;padded(1:mn_live_basis%nm,:)=n;number=reshape(padded,[size(number)])
    padded=0;padded(1:mn_live_basis%nm,:)=e;energy=reshape(padded,[size(energy)])
    ierr=mn_ok
  end subroutine

  subroutine mn_live_unpack(slot,stride,number,energy,ierr,validate_only)
    integer,intent(in) :: slot,stride
    real(mn_dp),intent(in) :: number(:),energy(:)
    integer,intent(out) :: ierr
    logical,optional,intent(in) :: validate_only
    real(mn_dp) :: n(stride,primary_groups),e(stride,primary_groups)
    ierr=mn_bad_input
    if(stride<mn_live_basis%nm.or.size(number)/=stride*primary_groups.or.size(energy)/=size(number))return
    if(any(.not.ieee_is_finite(number)).or.any(.not.ieee_is_finite(energy)))return
    n=reshape(number,shape(n));e=reshape(energy,shape(e))
    if(any(n(mn_live_basis%nm+1:,:)/=0).or.any(e(mn_live_basis%nm+1:,:)/=0))return
    call mn_live_write(slot,n(1:mn_live_basis%nm,:),e(1:mn_live_basis%nm,:),ierr,validate_only)
  end subroutine

  subroutine mn_live_ir_read(slot,energy,ierr)
    integer,intent(in) :: slot
    real(mn_dp),intent(out) :: energy(:,:)
    integer,intent(out) :: ierr
    integer :: block,offset
    energy=0;ierr=mn_bad_input
    if(size(energy,1)/=mn_live_basis%nm.or.size(energy,2)<1)return
    call location(slot,block,offset,ierr)
    if(ierr/=mn_ok)return
    if(.not.allocated(slabs(block)%ir))return
    ierr=mn_bad_input
    if(size(energy,2)/=size(slabs(block)%ir,2))return
    energy=slabs(block)%ir(:,:,offset);ierr=mn_ok
  end subroutine

  subroutine mn_live_ir_write(slot,energy,ierr,validate_only,reserve_only)
    integer,intent(in) :: slot
    real(mn_dp),intent(in) :: energy(:,:)
    integer,intent(out) :: ierr
    logical,optional,intent(in) :: validate_only,reserve_only
    real(mn_dp) :: angular(mn_live_basis%nq)
    integer :: block,offset,g,status
    ierr=mn_bad_input
    if(size(energy,1)/=mn_live_basis%nm.or.size(energy,2)<1)return
    do g=1,size(energy,2)
       call mn_reconstruct(mn_live_basis,energy(:,g),angular,ierr)
       if(ierr/=mn_ok)return
    enddo
    if(present(validate_only))then
       if(validate_only)return
    endif
    call location(slot,block,offset,ierr)
    if(ierr/=mn_ok)return
    if(.not.allocated(slabs(block)%ir))then
       if(all(energy==0))return
       allocate(slabs(block)%ir(mn_live_basis%nm,size(energy,2),slab_size),stat=status)
       if(status/=0)then
          ierr=mn_bad_input;return
       endif
       slabs(block)%ir=0
    endif
    ierr=mn_bad_input
    if(size(energy,2)/=size(slabs(block)%ir,2))return
    ierr=mn_ok
    if(present(reserve_only))then
       if(reserve_only)return
    endif
    slabs(block)%ir(:,:,offset)=energy;ierr=mn_ok
  end subroutine

  subroutine mn_live_ir_pack(slot,ng,stride,payload,ierr)
    integer,intent(in) :: slot,ng,stride
    real(mn_dp),intent(out) :: payload(:)
    integer,intent(out) :: ierr
    real(mn_dp) :: e(mn_live_basis%nm,ng),padded(ng,stride)
    ierr=mn_bad_input;payload=0
    if(stride<mn_live_basis%nm.or.size(payload)/=ng*stride)return
    call mn_live_ir_read(slot,e,ierr)
    if(ierr/=mn_ok)return
    ierr=mn_bad_input
    if(any(.not.ieee_is_finite(e)))return
    padded=0;padded(:,1:mn_live_basis%nm)=transpose(e);payload=reshape(padded,[size(payload)])
    ierr=mn_ok
  end subroutine

  subroutine mn_live_ir_unpack(slot,ng,stride,payload,ierr,validate_only)
    integer,intent(in) :: slot,ng,stride
    real(mn_dp),intent(in) :: payload(:)
    integer,intent(out) :: ierr
    logical,optional,intent(in) :: validate_only
    real(mn_dp) :: padded(ng,stride)
    ierr=mn_bad_input
    if(stride<mn_live_basis%nm.or.size(payload)/=ng*stride)return
    if(any(.not.ieee_is_finite(payload)))return
    padded=reshape(payload,shape(padded))
    if(any(padded(:,mn_live_basis%nm+1:)/=0))return
    call mn_live_ir_write(slot,transpose(padded(:,1:mn_live_basis%nm)),ierr,validate_only)
  end subroutine

  subroutine mn_live_ir_expand(factor,ierr)
    real(mn_dp),intent(in) :: factor
    integer,intent(out) :: ierr
    integer :: block
    ierr=mn_bad_input
    if(.not.ieee_is_finite(factor).or.factor<=0)return
    if(allocated(slabs))then
       if(factor>1)then
          do block=1,size(slabs)
             if(.not.allocated(slabs(block)%ir))cycle
             if(any(abs(slabs(block)%ir)>huge(1d0)/factor))return
          enddo
       endif
       do block=1,size(slabs)
          if(allocated(slabs(block)%ir))slabs(block)%ir=slabs(block)%ir*factor
       enddo
    endif
    ierr=mn_ok
  end subroutine
end module

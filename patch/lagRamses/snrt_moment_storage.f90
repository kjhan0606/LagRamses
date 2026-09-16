! Compact M_N state/transfer ABI. This does not alias or reinterpret SN arrays.
module snrt_moment_storage
  use snrt_moment_transport
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter :: version=1
  character(16), parameter :: magic='SNRT_MN_MB_V1'
  type, public :: mn_state
     type(mn_basis) :: basis
     real(mn_dp), allocatable :: value(:,:,:) ! moment, spectral group, cell
  end type
  public :: mn_state_initialize,mn_state_write,mn_state_read,mn_state_restrict,mn_state_inject
contains
  subroutine mn_state_initialize(state,order,ng,ncell,ierr)
    type(mn_state), intent(out) :: state
    integer, intent(in) :: order,ng,ncell
    integer, intent(out) :: ierr
    integer :: status
    ierr=mn_bad_input
    if(ng<1.or.ncell<0)return
    call mn_initialize(state%basis,order,ierr)
    if(ierr/=mn_ok)return
    allocate(state%value(state%basis%nm,ng,ncell),stat=status)
    if(status/=0)then
       ierr=mn_bad_input;return
    endif
    state%value=0
  end subroutine

  subroutine mn_state_write(state,path,ierr)
    ! Dedicated typed format; never overwrite an existing file. Existing SN
    ! restart files must not be interpreted as compact moments, or vice versa.
    type(mn_state), intent(in) :: state
    character(*), intent(in) :: path
    integer, intent(out) :: ierr
    integer :: unit,ios,closing,g,cell
    real(mn_dp) :: angular(state%basis%nq)
    ierr=mn_bad_input
    if(.not.allocated(state%value))return
    do cell=1,size(state%value,3)
       do g=1,size(state%value,2)
          call mn_reconstruct(state%basis,state%value(:,g,cell),angular,ierr)
          if(ierr/=mn_ok)return
       enddo
    enddo
    ierr=mn_bad_input
    open(newunit=unit,file=path,access='stream',form='unformatted',status='new',action='write',iostat=ios)
    if(ios/=0)return
    write(unit,iostat=ios)magic,version,state%basis%order,shape(state%value),state%value
    close(unit,iostat=closing)
    if(ios==0.and.closing==0)ierr=mn_ok
  end subroutine

  subroutine mn_state_read(state,path,ierr)
    ! Caller declares expected dimensions and closure BEFORE reading. Read
    ! into a trial, validate every cell, then publish; never resize from file.
    type(mn_state), intent(inout) :: state
    character(*), intent(in) :: path
    integer, intent(out) :: ierr
    character(16) :: signature
    integer :: unit,ios,closing,revision,order,dims(3),g,cell
    integer(int64) :: bytes,expected
    real(mn_dp), allocatable :: candidate(:,:,:)
    real(mn_dp) :: angular(state%basis%nq)
    ierr=mn_bad_input
    if(.not.allocated(state%value))return
    inquire(file=path,size=bytes,iostat=ios)
    if(ios/=0)return
    ! Storage sizes, rather than an assumed byte width for default integers.
    expected=len(magic)+5_int64*storage_size(revision)/8+ &
         size(state%value,kind=int64)*storage_size(1d0)/8
    if(bytes/=expected)return
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,iostat=ios)signature,revision,order,dims
    if(ios/=0)then
       close(unit);return
    endif
    if(signature/=magic.or.revision/=version.or.order/=state%basis%order.or.any(dims/=shape(state%value)))then
       close(unit);return
    endif
    allocate(candidate(size(state%value,1),size(state%value,2),size(state%value,3)))
    read(unit,iostat=ios)candidate
    close(unit,iostat=closing)
    if(ios/=0.or.closing/=0)return
    do cell=1,size(candidate,3)
       do g=1,size(candidate,2)
          call mn_reconstruct(state%basis,candidate(:,g,cell),angular,ierr)
          if(ierr/=mn_ok)return
       enddo
    enddo
    state%value=candidate;ierr=mn_ok
  end subroutine

  subroutine mn_state_restrict(b,child,volume,parent,ierr)
    ! Convex volume restriction preserves all moments and realizability.
    ! Prolongation without a slope is simply injection of the parent moments.
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(in) :: child(:,:),volume(:)
    real(mn_dp), intent(inout) :: parent(:)
    integer, intent(out) :: ierr
    real(mn_dp) :: candidate(b%nm),angular(b%nq),weights(size(volume))
    integer :: i
    ierr=mn_bad_input
    if(size(parent)/=b%nm.or.size(child,1)/=b%nm.or.size(child,2)/=size(volume).or.size(volume)==0)return
    if(any(.not.ieee_is_finite(volume)))return
    if(any(volume<=0))return
    do i=1,size(volume)
       call mn_reconstruct(b,child(:,i),angular,ierr)
       if(ierr/=mn_ok)return
    enddo
    weights=volume/maxval(volume);weights=weights/sum(weights)
    candidate=matmul(child,weights)
    call mn_reconstruct(b,candidate,angular,ierr)
    if(ierr/=mn_ok)return
    parent=candidate
  end subroutine

  subroutine mn_state_inject(b,u,source,ierr)
    type(mn_basis), intent(in) :: b
    real(mn_dp), intent(inout) :: u(:)
    real(mn_dp), intent(in) :: source(:)
    integer, intent(out) :: ierr
    real(mn_dp) :: candidate(b%nm),angular(b%nq)
    ierr=mn_bad_input
    if(size(u)/=b%nm.or.size(source)/=b%nm)return
    call mn_reconstruct(b,u,angular,ierr)
    if(ierr/=mn_ok)return
    call mn_reconstruct(b,source,angular,ierr)
    if(ierr/=mn_ok)return
    candidate=u+source
    call mn_reconstruct(b,candidate,angular,ierr)
    if(ierr/=mn_ok)return
    u=candidate
  end subroutine
end module

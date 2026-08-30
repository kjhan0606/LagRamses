! Tensor-Core angular-block interface for the lagRamses S_N RT backend.
! All arrays are contiguous single-precision buffers in row-major logical order:
! directional(nrow,ndirection), projection(ndirection,nbin), binned(nrow,nbin).
module snrt_cuda_interface
  use iso_c_binding, only: c_float, c_int, c_ptr, c_loc
  implicit none
  private

  public :: snrt_cuda_available
  public :: snrt_cuda_angular_reduce_tf32

  interface
     function snrt_cuda_available_c() bind(C,name='snrt_cuda_available_c') result(device_count)
       import :: c_int
       integer(c_int) :: device_count
     end function snrt_cuda_available_c

     function snrt_cuda_angular_reduce_tf32_c(directional,projection,binned, &
          & nrow,ndirection,nbin) bind(C,name='snrt_cuda_angular_reduce_tf32_c') result(ierr)
       import :: c_ptr, c_int
       type(c_ptr),value :: directional,projection,binned
       integer(c_int),value :: nrow,ndirection,nbin
       integer(c_int) :: ierr
     end function snrt_cuda_angular_reduce_tf32_c
  end interface

contains

  function snrt_cuda_available() result(device_count)
    integer :: device_count
    device_count=int(snrt_cuda_available_c())
  end function snrt_cuda_available

  subroutine snrt_cuda_angular_reduce_tf32(directional,projection,binned, &
       & nrow,ndirection,nbin,ierr)
    real(c_float),target,intent(in) :: directional(:),projection(:)
    real(c_float),target,intent(out) :: binned(:)
    integer,intent(in) :: nrow,ndirection,nbin
    integer,intent(out) :: ierr

    if(nrow<1 .or. ndirection<1 .or. nbin<1) then
       ierr=-1
       return
    endif
    if(size(directional)/=nrow*ndirection .or. size(projection)/=ndirection*nbin .or. &
         & size(binned)/=nrow*nbin) then
       ierr=-2
       return
    endif

    ierr=int(snrt_cuda_angular_reduce_tf32_c(c_loc(directional(1)),c_loc(projection(1)), &
         & c_loc(binned(1)),int(nrow,c_int),int(ndirection,c_int),int(nbin,c_int)))
  end subroutine snrt_cuda_angular_reduce_tf32

end module snrt_cuda_interface

program snrt_cuda_smoke
  use iso_c_binding, only: c_float
  use snrt_cuda_interface, only: snrt_cuda_available, &
       snrt_cuda_angular_reduce_tf32
  use snrt_cuda_ledger_interface, only: snrt_cuda_weighted_sum_fp32
  implicit none

  integer, parameter :: nrow = 64
  integer, parameter :: ndirection = 80
  integer, parameter :: nbin = 16
  integer :: idir, ibin, ierr, irow
  real(c_float) :: max_tensor_error, max_ledger_error, reference
  real(c_float), allocatable :: directional(:), weights(:), projection(:), &
       binned(:), scalar(:)

  if (.not. snrt_cuda_available()) error stop 'SNRT CUDA smoke: no CUDA device'
  allocate(directional(nrow * ndirection), weights(ndirection), &
       projection(ndirection * nbin), binned(nrow * nbin), scalar(nrow))
  directional = 1.0_c_float
  weights = 1.0_c_float / real(ndirection, c_float)
  projection = 0.0_c_float
  do idir = 1, ndirection
     ibin = 1 + mod(idir - 1, nbin)
     projection((idir - 1) * nbin + ibin) = weights(idir)
  enddo

  call snrt_cuda_angular_reduce_tf32(directional, projection, binned, nrow, &
       ndirection, nbin, ierr)
  if (ierr /= 0) error stop 'SNRT CUDA smoke: TF32 angular reduction failed'
  call snrt_cuda_weighted_sum_fp32(directional, weights, scalar, nrow, &
       ndirection, ierr)
  if (ierr /= 0) error stop 'SNRT CUDA smoke: FP32 photon ledger failed'

  reference = sum(weights)
  max_tensor_error = 0.0_c_float
  max_ledger_error = 0.0_c_float
  do irow = 1, nrow
     max_tensor_error = max(max_tensor_error, &
          abs(sum(binned((irow - 1) * nbin + 1:irow * nbin)) - reference))
     max_ledger_error = max(max_ledger_error, abs(scalar(irow) - reference))
  enddo
  if (max_ledger_error > 5.0e-6_c_float) then
     error stop 'SNRT CUDA smoke: FP32 photon ledger mismatch'
  endif
  write(*,'(A,ES12.4,A,ES12.4)') 'SNRT_CUDA_SMOKE_OK tensor_abs=', &
       max_tensor_error, ' ledger_abs=', max_ledger_error
end program snrt_cuda_smoke

program snrt_rt_transaction_mpi_smoke
  use, intrinsic :: iso_c_binding, only: c_float
  use mpi
  use amr_parameters, only: dp
  use snrt_rt_transaction
  implicit none

  integer :: ierr_mpi, rank, nrank, ierr, global_failed, global_converged
  real(dp) :: global_residual
  integer :: leaf_slot(0)
  real(c_float) :: persistent(2,2,0), coarse_flux(2,2,0)
  real(dp) :: hydrogen(0), helium_ii(0), helium_iii(0), neutral(0), thermal(0)
  real(dp) :: residual
  type(snrt_rt_transaction_snapshot) :: transaction

  call MPI_Init(ierr_mpi)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr_mpi)
  call MPI_Comm_size(MPI_COMM_WORLD, nrank, ierr_mpi)
  if (nrank < 2) then
     if (rank == 0) write(*,'(A)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_SKIP: need two ranks'
     call MPI_Finalize(ierr_mpi)
     stop 0
  end if

  call snrt_transaction_begin(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  if (ierr /= snrt_transaction_ok .or. .not. transaction%active) then
     write(*,'(A,I0)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_FAIL begin rank=', rank
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr_mpi)
  end if
  call snrt_transaction_restore(transaction, persistent, leaf_slot, hydrogen, &
       helium_ii, helium_iii, neutral, thermal, ierr)
  if (ierr /= snrt_transaction_ok) then
     write(*,'(A,I0)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_FAIL restore rank=', rank
     call MPI_Abort(MPI_COMM_WORLD, 2, ierr_mpi)
  end if

  residual = real(rank + 1, dp)
  call snrt_transaction_reduce_decision(merge(0,snrt_failure_partition,rank == 0), &
       1, residual, global_failed, global_converged, global_residual, ierr)
  if (ierr /= snrt_transaction_ok .or. global_failed /= snrt_failure_partition .or. &
       global_converged /= 1 .or. global_residual /= real(nrank,dp)) then
     write(*,'(A,I0)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_FAIL reduce rank=', rank
     call MPI_Abort(MPI_COMM_WORLD, 3, ierr_mpi)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
  if (rank == 0) write(*,'(A,I0)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS ranks=', nrank
  call MPI_Finalize(ierr_mpi)
end program snrt_rt_transaction_mpi_smoke

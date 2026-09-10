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
  call dry_commit_cases()
  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
  if (rank == 0) write(*,'(A,I0)') 'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_PASS ranks=', nrank
  call MPI_Finalize(ierr_mpi)
contains
  subroutine dry_commit_cases()
    integer :: slots(1), attempt, local_failed
    real(c_float) :: photons(2,2,2), trial(2,2,1), flux(2,2,2), staged(2,2,2)
    real(dp) :: shift(2,2,2), trial_shift(2,2,1), flux_shift(2,2,2), staged_shift(2,2,2)
    real(dp) :: h(2), he2(2), he3(2), hi(2), heat(1), reference(2)
    real(dp) :: th(1), the2(1), the3(1), thi(1), theat(1)
    slots=1;reference=[10.0_dp,20.0_dp]
    trial=3.0_c_float;flux=0.25_c_float
    trial_shift=-2.0_dp
    th=0.4_dp;the2=0.2_dp;the3=0.1_dp;thi=0.6_dp;theat=12.0_dp
    do attempt=1,2
       photons=2.0_c_float;shift=-1.0_dp
       h=0.1_dp;he2=0.05_dp;he3=0.01_dp;hi=0.9_dp;heat=10.0_dp
       call snrt_transaction_begin(transaction,photons,slots,h,he2,he3,hi,heat,ierr,shift,reference)
       call check(ierr==0.and.transaction%active,'shift begin')
       ! Model source deposition after the snapshot, before the dry commit.
       photons(:,:,1)=4.0_c_float;shift(:,:,1)=-0.5_dp
       staged=photons;staged_shift=shift
       flux_shift=0.125_dp
       ! Only rank 1 rejects, and only after checking the non-leaf coarse slot.
       if(attempt==1.and.rank==1)flux_shift(1,1,2)=-100.0_dp
       call snrt_transaction_commit_level(transaction,photons,slots,h,he2,he3,hi,trial,flux, &
            th,the2,the3,thi,heat,theat,ierr,shift,trial_shift,flux_shift,validate_only=.true.)
       call check((ierr/=0).eqv.(attempt==1.and.rank==1),'dry commit local decision')
       call check(transaction%active.and.all(photons==staged).and.all(shift==staged_shift).and. &
            all(h==0.1_dp).and.all(he2==0.05_dp).and.all(he3==0.01_dp).and. &
            all(hi==0.9_dp).and.all(heat==10.0_dp),'dry commit must not publish or clear snapshot')
       local_failed=snrt_failure_none
       if(ierr/=0)local_failed=snrt_failure_receiver
       call snrt_transaction_reduce_decision(local_failed,1,0.0_dp,global_failed, &
            global_converged,global_residual,ierr)
       call check(ierr==0,'dry commit collective')
       if(attempt==1)then
          call check(global_failed==snrt_failure_receiver,'peer rejection must reach all ranks')
          call snrt_transaction_restore(transaction,photons,slots,h,he2,he3,hi,heat,ierr,shift)
          call check(ierr==0.and..not.transaction%active.and.all(photons==2.0_c_float).and. &
               all(shift==-1.0_dp).and.all(h==0.1_dp).and.all(he2==0.05_dp).and. &
               all(he3==0.01_dp).and.all(hi==0.9_dp).and.all(heat==10.0_dp),'collective rollback')
       else
          call check(global_failed==snrt_failure_none,'valid dry commit must pass globally')
          call snrt_transaction_commit_level(transaction,photons,slots,h,he2,he3,hi,trial,flux, &
               th,the2,the3,thi,heat,theat,ierr,shift,trial_shift,flux_shift,validate_only=.false.)
          call check(ierr==0.and..not.transaction%active.and.all(photons(:,:,1)==3.25_c_float).and. &
               all(photons(:,:,2)==2.25_c_float).and.all(shift(:,:,1)==-1.875_dp).and. &
               all(shift(:,:,2)==-0.875_dp).and.h(1)==th(1).and.he2(1)==the2(1).and. &
               he3(1)==the3(1).and.hi(1)==thi(1).and.all(heat==theat).and. &
               h(2)==0.1_dp.and.he2(2)==0.05_dp.and.he3(2)==0.01_dp.and.hi(2)==0.9_dp, &
               'publish after collective acceptance')
       endif
    enddo
    if(rank==0)write(*,'(A)')'SNRT_NATIVE_RT_TRANSACTION_MPI_DRY_COMMIT_PASS'
  end subroutine dry_commit_cases

  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(condition)return
    write(*,'(A,I0,2A)')'SNRT_NATIVE_RT_TRANSACTION_MPI_SMOKE_FAIL rank=',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,4,ierr_mpi)
  end subroutine check
end program snrt_rt_transaction_mpi_smoke

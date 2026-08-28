! Adaptive record exchange shared by AMR grids and particles.
!
! The caller sorts records by destination once and supplies record counts and
! displacements.  A single collective decision then selects one of:
!   * sparse point-to-point for a sparse communication graph,
!   * MPI_Alltoallv for dense, bounded batches, or
!   * an incremental k-section destination schedule for the intermediate case.
!
! Keeping the policy here is important: FDM, DE, CDM and SIDM must not grow
! separate communication implementations with different memory behaviour.
module dynamic_exchange
  use amr_parameters, only: dp, i8b, exchange_method
  use amr_commons, only: myid, ncpu
  use mpi_large, only: mpi_large_isend_dp, mpi_large_irecv_dp, &
       mpi_large_dp_needed
  implicit none

  integer, parameter :: EXCHANGE_SPARSE_P2P = 1
  integer, parameter :: EXCHANGE_KSECTION   = 2
  integer, parameter :: EXCHANGE_ALLTOALLV  = 3

  real(dp), parameter :: P2P_DENSITY_LIMIT = 0.08_dp
  real(dp), parameter :: ALLTOALL_DENSITY_MIN = 0.35_dp
  integer(i8b), parameter :: ALLTOALL_MAX_RANK_BYTES = 256_i8b*1024_i8b*1024_i8b
  integer, parameter :: RECV_GUARD_COLUMNS = 8
  real(dp), parameter :: RECV_GUARD_VALUE = -huge(1.0_dp)

contains

  subroutine exchange_dp_records(sendbuf, nitem, dest_cpu, nprops, &
       recvbuf, nrecv, backend, label)
    real(dp), intent(in) :: sendbuf(:,:)
    integer, intent(in) :: nitem, dest_cpu(:), nprops
    real(dp), allocatable, intent(out) :: recvbuf(:,:)
    integer, intent(out) :: nrecv, backend
    character(len=*), intent(in), optional :: label
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ns(ncpu), nr(ncpu), sd(ncpu), rd(ncpu), offsets(ncpu)
    integer :: nse(ncpu), nre(ncpu), sde(ncpu), rde(ncpu)
    integer :: i, icpu, ierr, nq, req(2*ncpu), unsafe_local, unsafe_global
    real(dp), allocatable :: sorted(:,:)

    ns=0
    do i=1,nitem
       ns(dest_cpu(i))=ns(dest_cpu(i))+1
    end do
#ifndef WITHOUTMPI
    call MPI_ALLTOALL(ns,1,MPI_INTEGER,nr,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
#else
    nr=ns
#endif
    sd(1)=0; rd(1)=0
    do icpu=2,ncpu
       sd(icpu)=sd(icpu-1)+ns(icpu-1)
       rd(icpu)=rd(icpu-1)+nr(icpu-1)
    end do
    nrecv=sum(nr)

    if(present(label)) then
       call choose_exchange_backend(ns,nr,8*nprops,backend,label)
    else
       call choose_exchange_backend(ns,nr,8*nprops,backend)
    end if

    ! Alltoallv and the flat sectioned schedule use default-integer element
    ! counts and displacements.  A record count can be perfectly valid while
    ! count*nprops overflows 32 bits, especially after k-section aggregation.
    ! Direct P2P retains record counts and lets mpi_large construct a record
    ! datatype, so force that path before doing any default-integer product.
    unsafe_local=0
    do icpu=1,ncpu
       if(mpi_large_dp_needed(ns(icpu),nprops) .or. &
            mpi_large_dp_needed(nr(icpu),nprops) .or. &
            mpi_large_dp_needed(sd(icpu),nprops) .or. &
            mpi_large_dp_needed(rd(icpu),nprops)) unsafe_local=1
    end do
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(unsafe_local,unsafe_global,1,MPI_INTEGER,MPI_MAX, &
         MPI_COMM_WORLD,ierr)
#else
    unsafe_global=unsafe_local
#endif
    if(unsafe_global/=0 .and. backend/=EXCHANGE_SPARSE_P2P)then
       backend=EXCHANGE_SPARSE_P2P
       if(myid==1) write(*,*) &
            ' Dynamic exchange: forcing large-count P2P (32-bit element limit)'
    end if
    allocate(sorted(nprops,max(nitem,1)))
    offsets=sd
    do i=1,nitem
       icpu=dest_cpu(i)
       offsets(icpu)=offsets(icpu)+1
       sorted(:,offsets(icpu))=sendbuf(:,i)
    end do
    ! Guard columns turn a bad MPI count/displacement into an immediate,
    ! attributable failure instead of delayed allocator corruption.
    allocate(recvbuf(nprops,max(nrecv,1)+RECV_GUARD_COLUMNS))
    recvbuf(:,max(nrecv,1)+1:max(nrecv,1)+RECV_GUARD_COLUMNS)=RECV_GUARD_VALUE
#ifndef WITHOUTMPI
    if(backend==EXCHANGE_ALLTOALLV) then
       nse=ns*nprops; nre=nr*nprops
       sde=sd*nprops; rde=rd*nprops
       call MPI_ALLTOALLV(sorted,nse,sde,MPI_DOUBLE_PRECISION, &
            recvbuf,nre,rde,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    else if(backend==EXCHANGE_KSECTION) then
       nse=ns*nprops; nre=nr*nprops
       sde=sd*nprops; rde=rd*nprops
       call exchange_sectioned_dp(sorted,recvbuf,nse,nre,sde,rde,701)
    else
       nq=0
       do icpu=1,ncpu
          if(icpu/=myid.and.nr(icpu)>0) then
             nq=nq+1
             call mpi_large_irecv_dp(recvbuf(1,rd(icpu)+1),nr(icpu), &
                  nprops,icpu-1,700,MPI_COMM_WORLD,req(nq),ierr)
          end if
       end do
       do icpu=1,ncpu
          if(icpu/=myid.and.ns(icpu)>0) then
             nq=nq+1
             call mpi_large_isend_dp(sorted(1,sd(icpu)+1),ns(icpu), &
                  nprops,icpu-1,700,MPI_COMM_WORLD,req(nq),ierr)
          end if
       end do
       if(ns(myid)>0) recvbuf(:,rd(myid)+1:rd(myid)+nr(myid)) = &
            sorted(:,sd(myid)+1:sd(myid)+ns(myid))
       if(nq>0) call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
    end if
#else
    if(nitem>0) recvbuf(:,1:nitem)=sorted(:,1:nitem)
#endif
    if(any(recvbuf(:,max(nrecv,1)+1:max(nrecv,1)+RECV_GUARD_COLUMNS) /= &
         RECV_GUARD_VALUE)) then
       write(*,*) 'FATAL: dynamic exchange receive guard overwritten',myid,nrecv,nprops
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD,91,ierr)
#else
       stop 91
#endif
    end if
    deallocate(sorted)
  end subroutine exchange_dp_records

  subroutine choose_exchange_backend(nsend, nrecv, record_bytes, backend, label, &
       p2p_payload_limit)
    integer, intent(in) :: nsend(:), nrecv(:), record_bytes
    integer, intent(out) :: backend
    character(len=*), intent(in), optional :: label
    integer(i8b), intent(in), optional :: p2p_payload_limit
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: icpu, ierr, local_partners, max_partners
    integer(i8b) :: local_edges, global_edges, local_bytes, max_rank_bytes
    real(dp) :: density
    character(len=32) :: backend_name

    local_partners = 0
    local_edges = 0_i8b
    do icpu = 1, ncpu
       if(icpu == myid) cycle
       if(nsend(icpu) > 0) then
          local_edges = local_edges + 1_i8b
          local_partners = local_partners + 1
       end if
       if(nrecv(icpu) > 0) local_partners = max(local_partners, 1)
    end do
    local_partners = max(local_partners, count(nrecv > 0) - merge(1, 0, nrecv(myid) > 0))
    local_bytes = int(sum(nsend), i8b) * int(record_bytes, i8b)

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(local_edges, global_edges, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(local_partners, max_partners, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(local_bytes, max_rank_bytes, 1, MPI_INTEGER8, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
    global_edges = local_edges
    max_partners = local_partners
    max_rank_bytes = local_bytes
#endif

    if(trim(exchange_method)=='p2p')then
       density=0d0
       backend=EXCHANGE_SPARSE_P2P
    else if(trim(exchange_method)=='ksection')then
       density=0d0
       backend=EXCHANGE_KSECTION
    else if(trim(exchange_method)=='alltoallv')then
       density=1d0
       backend=EXCHANGE_ALLTOALLV
    else if(ncpu <= 1) then
       density = 0.0_dp
       backend = EXCHANGE_SPARSE_P2P
    else if(present(p2p_payload_limit) .and. &
         max_rank_bytes <= p2p_payload_limit) then
       ! Repeated exchanges with pre-existing rank buffers should not pay a
       ! collective pack/unpack penalty for tiny payloads, even when their
       ! rank graph is topologically dense.
       density = real(global_edges, dp) / real(int(ncpu, i8b)*int(ncpu-1, i8b), dp)
       backend = EXCHANGE_SPARSE_P2P
    else
       density = real(global_edges, dp) / real(int(ncpu, i8b)*int(ncpu-1, i8b), dp)
       if(max_partners <= max(4, ncpu/16) .or. density <= P2P_DENSITY_LIMIT) then
          backend = EXCHANGE_SPARSE_P2P
       else if(density >= ALLTOALL_DENSITY_MIN .and. &
            max_rank_bytes <= ALLTOALL_MAX_RANK_BYTES) then
          backend = EXCHANGE_ALLTOALLV
       else
          backend = EXCHANGE_KSECTION
       end if
    end if

    select case(backend)
    case(EXCHANGE_SPARSE_P2P); backend_name = 'sparse-p2p'
    case(EXCHANGE_KSECTION);   backend_name = 'k-section'
    case default;              backend_name = 'mpi-alltoallv'
    end select
    if(myid == 1 .and. present(label)) then
          write(*,'(A,A,A,A,A,F6.3,A,I6,A,F10.2,A)') ' Dynamic exchange [', trim(label), &
               ']: ', trim(backend_name), ' density=', density, ' max_partners=', max_partners, &
               ' max_rank_payload=', real(max_rank_bytes,dp)/(1024.0_dp*1024.0_dp), ' MiB'
    end if
  end subroutine choose_exchange_backend

  subroutine exchange_dp_sorted(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    real(dp), intent(in) :: sbuf(:)
    real(dp), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    call exchange_native_dp(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
  end subroutine exchange_dp_sorted

  subroutine exchange_i8_sorted(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(i8b), intent(in) :: sbuf(:)
    integer(i8b), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    real(dp), allocatable :: send_dp(:), recv_dp(:)
    integer :: i

    if(backend /= EXCHANGE_KSECTION) then
       call exchange_native_i8(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
       return
    end if
    allocate(send_dp(max(sum(ns),1)), recv_dp(max(sum(nr),1)))
    do i = 1, sum(ns)
       send_dp(i) = transfer(sbuf(i), send_dp(i))
    end do
    call exchange_dp_sorted(send_dp, recv_dp, ns, nr, sd, rd, backend, tag)
    do i = 1, sum(nr)
       rbuf(i) = transfer(recv_dp(i), rbuf(i))
    end do
    deallocate(send_dp, recv_dp)
  end subroutine exchange_i8_sorted

  subroutine exchange_int_sorted(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer, intent(in) :: sbuf(:)
    integer, intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    real(dp), allocatable :: send_dp(:), recv_dp(:)
    integer :: i

    if(backend /= EXCHANGE_KSECTION) then
       call exchange_native_int(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
       return
    end if
    allocate(send_dp(max(sum(ns),1)), recv_dp(max(sum(nr),1)))
    do i = 1, sum(ns); send_dp(i) = real(sbuf(i), dp); end do
    call exchange_dp_sorted(send_dp, recv_dp, ns, nr, sd, rd, backend, tag)
    do i = 1, sum(nr); rbuf(i) = nint(recv_dp(i)); end do
    deallocate(send_dp, recv_dp)
  end subroutine exchange_int_sorted

  subroutine exchange_i1_sorted(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(kind=1), intent(in) :: sbuf(:)
    integer(kind=1), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    real(dp), allocatable :: send_dp(:), recv_dp(:)
    integer :: i

    if(backend /= EXCHANGE_KSECTION) then
       call exchange_native_i1(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
       return
    end if
    allocate(send_dp(max(sum(ns),1)), recv_dp(max(sum(nr),1)))
    do i = 1, sum(ns); send_dp(i) = real(sbuf(i), dp); end do
    call exchange_dp_sorted(send_dp, recv_dp, ns, nr, sd, rd, backend, tag)
    do i = 1, sum(nr); rbuf(i) = int(nint(recv_dp(i)), kind=1); end do
    deallocate(send_dp, recv_dp)
  end subroutine exchange_i1_sorted

  subroutine exchange_native_dp(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    real(dp), intent(in) :: sbuf(:)
    real(dp), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    call exchange_native_impl_dp(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
  end subroutine exchange_native_dp

  subroutine exchange_native_i8(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(i8b), intent(in) :: sbuf(:)
    integer(i8b), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    call exchange_native_impl_i8(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
  end subroutine exchange_native_i8

  subroutine exchange_native_int(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer, intent(in) :: sbuf(:)
    integer, intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    call exchange_native_impl_int(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
  end subroutine exchange_native_int

  subroutine exchange_native_i1(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(kind=1), intent(in) :: sbuf(:)
    integer(kind=1), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:), nr(:), sd(:), rd(:), backend, tag
    call exchange_native_impl_i1(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
  end subroutine exchange_native_i1

  ! The following four routines deliberately remain type-specific.  This keeps
  ! compatibility with the legacy mpif.h interface used throughout RAMSES.
  subroutine exchange_native_impl_dp(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    real(dp), intent(in) :: sbuf(:); real(dp), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:),nr(:),sd(:),rd(:),backend,tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ierr
#ifndef WITHOUTMPI
    if(backend == EXCHANGE_ALLTOALLV) then
       call MPI_ALLTOALLV(sbuf,ns,sd,MPI_DOUBLE_PRECISION,rbuf,nr,rd,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    else if(backend == EXCHANGE_KSECTION) then
       call exchange_sectioned_dp(sbuf,rbuf,ns,nr,sd,rd,tag)
    else
       call p2p_dp(sbuf,rbuf,ns,nr,sd,rd,tag)
    end if
#else
    if(sum(ns)>0) rbuf(1:sum(ns))=sbuf(1:sum(ns))
#endif
  end subroutine exchange_native_impl_dp

  subroutine exchange_sectioned_dp(sbuf,rbuf,ns,nr,sd,rd,tag)
    real(dp),intent(in)::sbuf(*)
    real(dp),intent(out)::rbuf(*)
    integer,intent(in)::ns(:),nr(:),sd(:),rd(:),tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer::block_lo,block_hi,block_width,icpu,nq,ierr
    integer::req(2*ncpu)

    ! Incremental k-section schedule.  A round activates only one contiguous
    ! destination section.  All payload remains in the caller's sorted grid
    ! or particle arrays; no Hilbert/Morton topology is reconstructed here.
    block_width=max(1,int(sqrt(real(ncpu,dp))))
#ifndef WITHOUTMPI
    do block_lo=1,ncpu,block_width
       block_hi=min(ncpu,block_lo+block_width-1)
       nq=0
       if(myid>=block_lo.and.myid<=block_hi)then
          do icpu=1,ncpu
             if(icpu/=myid.and.nr(icpu)>0)then
                nq=nq+1
                call MPI_IRECV(rbuf(rd(icpu)+1),nr(icpu),MPI_DOUBLE_PRECISION, &
                     icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
             endif
          end do
       endif
       do icpu=block_lo,block_hi
          if(icpu/=myid.and.ns(icpu)>0)then
             nq=nq+1
             call MPI_ISEND(sbuf(sd(icpu)+1),ns(icpu),MPI_DOUBLE_PRECISION, &
                  icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
          endif
       end do
       if(myid>=block_lo.and.myid<=block_hi.and.ns(myid)>0) &
            rbuf(rd(myid)+1:rd(myid)+nr(myid))= &
            sbuf(sd(myid)+1:sd(myid)+ns(myid))
       if(nq>0)call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
    end do
#else
    if(sum(ns)>0)rbuf(1:sum(ns))=sbuf(1:sum(ns))
#endif
  end subroutine exchange_sectioned_dp

  subroutine exchange_native_impl_i8(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(i8b), intent(in) :: sbuf(:); integer(i8b), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:),nr(:),sd(:),rd(:),backend,tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ierr
#ifndef WITHOUTMPI
    if(backend == EXCHANGE_ALLTOALLV) then
       call MPI_ALLTOALLV(sbuf,ns,sd,MPI_INTEGER8,rbuf,nr,rd,MPI_INTEGER8,MPI_COMM_WORLD,ierr)
    else
       call p2p_i8(sbuf,rbuf,ns,nr,sd,rd,tag)
    end if
#else
    if(sum(ns)>0) rbuf(1:sum(ns))=sbuf(1:sum(ns))
#endif
  end subroutine exchange_native_impl_i8

  subroutine exchange_native_impl_int(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer, intent(in) :: sbuf(:); integer, intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:),nr(:),sd(:),rd(:),backend,tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ierr
#ifndef WITHOUTMPI
    if(backend == EXCHANGE_ALLTOALLV) then
       call MPI_ALLTOALLV(sbuf,ns,sd,MPI_INTEGER,rbuf,nr,rd,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    else
       call p2p_int(sbuf,rbuf,ns,nr,sd,rd,tag)
    end if
#else
    if(sum(ns)>0) rbuf(1:sum(ns))=sbuf(1:sum(ns))
#endif
  end subroutine exchange_native_impl_int

  subroutine exchange_native_impl_i1(sbuf, rbuf, ns, nr, sd, rd, backend, tag)
    integer(kind=1), intent(in) :: sbuf(:); integer(kind=1), intent(out) :: rbuf(:)
    integer, intent(in) :: ns(:),nr(:),sd(:),rd(:),backend,tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ierr
#ifndef WITHOUTMPI
    if(backend == EXCHANGE_ALLTOALLV) then
       call MPI_ALLTOALLV(sbuf,ns,sd,MPI_INTEGER1,rbuf,nr,rd,MPI_INTEGER1,MPI_COMM_WORLD,ierr)
    else
       call p2p_i1(sbuf,rbuf,ns,nr,sd,rd,tag)
    end if
#else
    if(sum(ns)>0) rbuf(1:sum(ns))=sbuf(1:sum(ns))
#endif
  end subroutine exchange_native_impl_i1

#ifndef WITHOUTMPI
  subroutine p2p_dp(sbuf,rbuf,ns,nr,sd,rd,tag)
    include 'mpif.h'
    real(dp),intent(in)::sbuf(:); real(dp),intent(out)::rbuf(:)
    integer,intent(in)::ns(:),nr(:),sd(:),rd(:),tag
    integer::icpu,nq,ierr,req(2*ncpu)
    nq=0
    do icpu=1,ncpu
       if(icpu/=myid.and.nr(icpu)>0)then
          nq=nq+1
          call MPI_IRECV(rbuf(rd(icpu)+1),nr(icpu),MPI_DOUBLE_PRECISION, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    do icpu=1,ncpu
       if(icpu/=myid.and.ns(icpu)>0)then
          nq=nq+1
          call MPI_ISEND(sbuf(sd(icpu)+1),ns(icpu),MPI_DOUBLE_PRECISION, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    if(ns(myid)>0) rbuf(rd(myid)+1:rd(myid)+nr(myid))=sbuf(sd(myid)+1:sd(myid)+ns(myid))
    if(nq>0) call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
  end subroutine p2p_dp
  subroutine p2p_i8(sbuf,rbuf,ns,nr,sd,rd,tag)
    include 'mpif.h'
    integer(i8b),intent(in)::sbuf(:); integer(i8b),intent(out)::rbuf(:)
    integer,intent(in)::ns(:),nr(:),sd(:),rd(:),tag
    integer::icpu,nq,ierr,req(2*ncpu)
    nq=0
    do icpu=1,ncpu
       if(icpu/=myid.and.nr(icpu)>0)then
          nq=nq+1
          call MPI_IRECV(rbuf(rd(icpu)+1),nr(icpu),MPI_INTEGER8, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    do icpu=1,ncpu
       if(icpu/=myid.and.ns(icpu)>0)then
          nq=nq+1
          call MPI_ISEND(sbuf(sd(icpu)+1),ns(icpu),MPI_INTEGER8, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    if(ns(myid)>0) rbuf(rd(myid)+1:rd(myid)+nr(myid))=sbuf(sd(myid)+1:sd(myid)+ns(myid))
    if(nq>0) call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
  end subroutine p2p_i8
  subroutine p2p_int(sbuf,rbuf,ns,nr,sd,rd,tag)
    include 'mpif.h'
    integer,intent(in)::sbuf(:); integer,intent(out)::rbuf(:)
    integer,intent(in)::ns(:),nr(:),sd(:),rd(:),tag
    integer::icpu,nq,ierr,req(2*ncpu)
    nq=0
    do icpu=1,ncpu
       if(icpu/=myid.and.nr(icpu)>0)then
          nq=nq+1
          call MPI_IRECV(rbuf(rd(icpu)+1),nr(icpu),MPI_INTEGER, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    do icpu=1,ncpu
       if(icpu/=myid.and.ns(icpu)>0)then
          nq=nq+1
          call MPI_ISEND(sbuf(sd(icpu)+1),ns(icpu),MPI_INTEGER, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    if(ns(myid)>0) rbuf(rd(myid)+1:rd(myid)+nr(myid))=sbuf(sd(myid)+1:sd(myid)+ns(myid))
    if(nq>0) call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
  end subroutine p2p_int
  subroutine p2p_i1(sbuf,rbuf,ns,nr,sd,rd,tag)
    include 'mpif.h'
    integer(kind=1),intent(in)::sbuf(:); integer(kind=1),intent(out)::rbuf(:)
    integer,intent(in)::ns(:),nr(:),sd(:),rd(:),tag
    integer::icpu,nq,ierr,req(2*ncpu)
    nq=0
    do icpu=1,ncpu
       if(icpu/=myid.and.nr(icpu)>0)then
          nq=nq+1
          call MPI_IRECV(rbuf(rd(icpu)+1),nr(icpu),MPI_INTEGER1, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    do icpu=1,ncpu
       if(icpu/=myid.and.ns(icpu)>0)then
          nq=nq+1
          call MPI_ISEND(sbuf(sd(icpu)+1),ns(icpu),MPI_INTEGER1, &
               icpu-1,tag,MPI_COMM_WORLD,req(nq),ierr)
       endif
    enddo
    if(ns(myid)>0) rbuf(rd(myid)+1:rd(myid)+nr(myid))=sbuf(sd(myid)+1:sd(myid)+ns(myid))
    if(nq>0) call MPI_WAITALL(nq,req,MPI_STATUSES_IGNORE,ierr)
  end subroutine p2p_i1
#endif

  subroutine check_receive_count(actual, expected, tag)
    integer, intent(in) :: actual, expected, tag
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    integer :: ierr
    if(actual /= expected) then
       write(*,*) 'Dynamic exchange count mismatch:', myid, tag, actual, expected
#ifndef WITHOUTMPI
       call MPI_ABORT(MPI_COMM_WORLD, 91, ierr)
#else
       stop 91
#endif
    end if
  end subroutine check_receive_count

end module dynamic_exchange

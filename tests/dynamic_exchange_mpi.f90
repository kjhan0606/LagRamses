program test_dynamic_exchange_mpi
  use amr_parameters, only: dp, i8b
  use amr_commons, only: myid, ncpu, ksec_root
  use ksection, only: init_ksection_comm_tree
  use dynamic_exchange
  use mpi_large, only: mpi_large_dp_needed
  implicit none
  include 'mpif.h'

  integer :: ierr, nfail_local, nfail_global

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
  myid = myid + 1

  ksec_root = 0
  call init_ksection_comm_tree()

  nfail_local = 0
  if(mpi_large_dp_needed(23598721,91)) then
     write(*,*) 'large-count threshold false positive',myid
     nfail_local=nfail_local+1
  end if
  if(.not.mpi_large_dp_needed(23598722,91)) then
     write(*,*) 'large-count threshold false negative',myid
     nfail_local=nfail_local+1
  end if
  call run_pattern('sparse', 1, EXCHANGE_SPARSE_P2P, 1000, nfail_local)
  call run_pattern('medium', min(6,ncpu-1), EXCHANGE_KSECTION, 1100, nfail_local)
  call run_pattern('dense', ncpu-1, EXCHANGE_ALLTOALLV, 1200, nfail_local)
  call run_record_pattern('record-sparse', 1, EXCHANGE_SPARSE_P2P, nfail_local)
  call run_record_pattern('record-medium', min(6,ncpu-1), EXCHANGE_KSECTION, nfail_local)
  call run_record_pattern('record-dense', ncpu-1, EXCHANGE_ALLTOALLV, nfail_local)

  call MPI_ALLREDUCE(nfail_local, nfail_global, 1, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD, ierr)
  if(myid == 1) then
     if(nfail_global == 0) then
        write(*,*) 'DYNAMIC_EXCHANGE_TEST PASS'
     else
        write(*,*) 'DYNAMIC_EXCHANGE_TEST FAIL count=', nfail_global
     end if
  end if
  call MPI_FINALIZE(ierr)
  if(nfail_global /= 0) stop 2

contains

  subroutine run_pattern(name, npartner, expected_backend, tag0, nfail)
    character(len=*), intent(in) :: name
    integer, intent(in) :: npartner, expected_backend, tag0
    integer, intent(inout) :: nfail
    integer :: ns(ncpu), nr(ncpu), sd(ncpu), rd(ncpu), offsets(ncpu)
    integer :: backend, i, k, dest, ntot_send, ntot_recv, code_int
    integer(i8b) :: code
    integer, allocatable :: dest_for_item(:), send_int(:), recv_int(:)
    integer(i8b), allocatable :: send_i8(:), recv_i8(:)
    integer(kind=1), allocatable :: send_i1(:), recv_i1(:)
    real(dp), allocatable :: send_dp(:), recv_dp(:)

    ns = 0
    do k = 1, npartner
       dest = mod(myid-1+k, ncpu) + 1
       ns(dest) = 2
    end do
    call MPI_ALLTOALL(ns, 1, MPI_INTEGER, nr, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)

    sd(1)=0; rd(1)=0
    do i=2,ncpu
       sd(i)=sd(i-1)+ns(i-1)
       rd(i)=rd(i-1)+nr(i-1)
    end do
    ntot_send=sum(ns); ntot_recv=sum(nr)
    allocate(dest_for_item(max(ntot_send,1)), send_dp(max(ntot_send,1)), &
         recv_dp(max(ntot_recv,1)), send_i8(max(ntot_send,1)), &
         recv_i8(max(ntot_recv,1)), send_int(max(ntot_send,1)), &
         recv_int(max(ntot_recv,1)), send_i1(max(ntot_send,1)), &
         recv_i1(max(ntot_recv,1)))

    offsets=sd
    do k=1,npartner
       dest=mod(myid-1+k,ncpu)+1
       do i=1,2
          offsets(dest)=offsets(dest)+1
          code=int(myid,i8b)*1000000_i8b+int(dest,i8b)*1000_i8b+int(i,i8b)
          dest_for_item(offsets(dest))=dest
          send_i8(offsets(dest))=code
          send_dp(offsets(dest))=real(code,dp)
          code_int=int(mod(code,100000_i8b))
          send_int(offsets(dest))=code_int
          send_i1(offsets(dest))=int(mod(code,101_i8b)-50_i8b,kind=1)
       end do
    end do

    call choose_exchange_backend(ns,nr,21,backend,name)
    if(backend /= expected_backend) then
       write(*,*) 'backend mismatch',myid,trim(name),backend,expected_backend
       nfail=nfail+1
    end if
    call exchange_dp_sorted(send_dp,recv_dp,ns,nr,sd,rd,backend,tag0)
    call exchange_i8_sorted(send_i8,recv_i8,ns,nr,sd,rd,backend,tag0+1)
    call exchange_int_sorted(send_int,recv_int,ns,nr,sd,rd,backend,tag0+2)
    call exchange_i1_sorted(send_i1,recv_i1,ns,nr,sd,rd,backend,tag0+3)

    do i=1,ntot_recv
       code=recv_i8(i)
       dest=int(mod(code/1000_i8b,1000_i8b))
       if(dest /= myid) nfail=nfail+1
       if(recv_dp(i) /= real(code,dp)) nfail=nfail+1
       if(recv_int(i) /= int(mod(code,100000_i8b))) nfail=nfail+1
       if(recv_i1(i) /= int(mod(code,101_i8b)-50_i8b,kind=1)) nfail=nfail+1
    end do

    deallocate(dest_for_item,send_dp,recv_dp,send_i8,recv_i8, &
         send_int,recv_int,send_i1,recv_i1)
  end subroutine run_pattern

  subroutine run_record_pattern(name, npartner, expected_backend, nfail)
    character(len=*), intent(in) :: name
    integer, intent(in) :: npartner, expected_backend
    integer, intent(inout) :: nfail
    integer :: i, k, dest, nsend, nrecv, backend, source, serial
    integer, allocatable :: dest_cpu(:)
    logical :: seen(ncpu,2)
    real(dp), allocatable :: send_record(:,:), recv_record(:,:)

    nsend = 2*npartner
    allocate(dest_cpu(max(nsend,1)), send_record(3,max(nsend,1)))
    i = 0
    do k = 1, npartner
       dest = mod(myid-1+k,ncpu)+1
       i = i+1
       dest_cpu(i)=dest
       send_record(:,i)=(/real(myid,dp),real(dest,dp),1.0_dp/)
       i = i+1
       dest_cpu(i)=dest
       send_record(:,i)=(/real(myid,dp),real(dest,dp),2.0_dp/)
    end do

    call exchange_dp_records(send_record,nsend,dest_cpu,3,recv_record, &
         nrecv,backend,name)
    if(backend /= expected_backend) then
       write(*,*) 'record backend mismatch',myid,trim(name),backend,expected_backend
       nfail=nfail+1
    end if
    if(nrecv /= 2*npartner) then
       write(*,*) 'record receive count mismatch',myid,trim(name),nrecv,2*npartner
       nfail=nfail+1
    end if
    seen=.false.
    do i=1,nrecv
       source=nint(recv_record(1,i))
       dest=nint(recv_record(2,i))
       serial=nint(recv_record(3,i))
       if(source<1 .or. source>ncpu .or. dest/=myid .or. &
            serial<1 .or. serial>2) then
          write(*,*) 'invalid record',myid,trim(name),recv_record(:,i)
          nfail=nfail+1
       else if(seen(source,serial)) then
          write(*,*) 'duplicate record',myid,trim(name),source,serial
          nfail=nfail+1
       else
          seen(source,serial)=.true.
       end if
    end do
    if(count(seen)/=2*npartner) then
       write(*,*) 'record set mismatch',myid,trim(name),count(seen),2*npartner
       nfail=nfail+1
    end if
    deallocate(dest_cpu,send_record,recv_record)
  end subroutine run_record_pattern

end program test_dynamic_exchange_mpi

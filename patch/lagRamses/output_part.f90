! Fail-closed lagRamses wrapper around the native binary particle writer.
! The payload implementation is included unchanged below, under a private
! name, so valid checkpoints retain their established record layout.
subroutine backup_part(filename)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=80),intent(in)::filename
  integer::npart_active,info

  ! Every payload field below is packed by this exact predicate.  Verify its
  ! cardinality before opening the output or allocating a buffer sized npart;
  ! otherwise xdp/ii8/ll/l1 writes run past their bounds and can overwrite
  ! live AMR topology before the checkpoint is complete.
  npart_active=count(levelp(1:npartmax)>0)
  if(npart<0.or.npart>npartmax.or.npart_active/=npart)then
     write(*,*)'BACKUP_PART_CARDINALITY_INVALID', &
          ' myid=',myid,' npart=',npart,' active=',npart_active, &
          ' npartmax=',npartmax,' numbp_free=',numbp_free
     call flush(6)
#ifndef WITHOUTMPI
     ! This check precedes the rank-serialized output token.  Abort directly:
     ! the selected clean_stop finalizes collectively and can deadlock when
     ! only one rank detects a local cardinality mismatch.
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop 1
#endif
  endif

  call backup_part_payload(filename)
end subroutine backup_part

#define backup_part backup_part_payload
#include "../cuRamses/output_part.f90"
#undef backup_part

program ramses
#ifdef PHASE0_STELLAR_ENRICHMENT
  use stellar_enrichment_config, only: use_channel_resolved_feedback
  use stellar_ramses_runtime, only: phase0_initialize
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
#ifdef PHASE0_STELLAR_ENRICHMENT
  integer :: phase0_ierr, phase0_mpi_ierr
#endif

  ! Read run parameters
  call read_params

#ifdef PHASE0_STELLAR_ENRICHMENT
  ! Validate the selected feedback source before allocating ICs or advancing
  ! the hydro state.  The cached initialization is reused by feedback calls.
  if(use_channel_resolved_feedback())then
     call phase0_initialize(phase0_ierr)
     if(phase0_ierr/=0)then
        write(*,*)'Phase 0 startup preflight failed: ',phase0_ierr
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,phase0_ierr,phase0_mpi_ierr)
#else
        error stop 1
#endif
     endif
  endif
#endif

  ! Set signal handler
  call set_signal_handler

  ! Start time integration
  call adaptive_loop

end program ramses

! sets the hook to catch signal 10
subroutine set_signal_handler
  implicit none
  external output_signal
  integer::jsigact,jhandle
  integer::istatus=-1

#ifndef NOSYSTEM
  call SIGNAL(10,output_signal,istatus)
#else
  call PXFSTRUCTCREATE("sigaction",jsigact,istatus)
  call PXFGETSUBHANDLE(output_signal,jhandle,istatus)
  call PXFINTSET(jsigact,"sa_handler",jhandle,istatus)
  call PXFSIGACTION(10,jsigact,0,istatus)
#endif
end subroutine set_signal_handler

! signal handler subroutine
subroutine output_signal
  use amr_commons
  implicit none

  if (myid==1) write (*,*) 'SIGNAL 10: Output will be written to disk during next main step.'

  ! output will be written to disk at next main step
  output_now = .true.

end subroutine output_signal

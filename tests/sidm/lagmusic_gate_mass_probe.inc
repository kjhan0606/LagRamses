! Test-only subroutine, appended to an isolated init_part source snapshot.
! No particle state is changed and no output snapshot is requested.
subroutine lagmusic_gate_mass_probe
  use amr_commons
  use pm_commons
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'mpif.h'
  real(dp) :: local_mass, global_mass, tier_mass
  integer(kind=8) :: local_count, global_count
  integer :: ierr, ilevel, bad_local, bad_global
  local_mass=sum(mp(1:npart))
  call MPI_ALLREDUCE(local_mass,global_mass,1,MPI_DOUBLE_PRECISION, &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  local_count=int(npart,kind=8)
  call MPI_ALLREDUCE(local_count,global_count,1,MPI_INTEGER8, &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  bad_local=count(.not.ieee_is_finite(mp(1:npart)))
  bad_local=bad_local+count(.not.ieee_is_finite(xp(1:npart,1:ndim)))
  bad_local=bad_local+count(.not.ieee_is_finite(vp(1:npart,1:ndim)))
  call MPI_ALLREDUCE(bad_local,bad_global,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myid==1)write(*,'(A,ES24.16,A,I0,A,I0)') &
       'LM_GATE_MASS: mass=',global_mass,' count=',global_count,' bad=',bad_global
  do ilevel=levelmin,nlevelmax
     tier_mass=0.5d0**(3*ilevel)
     local_count=int(count(mp(1:npart)==tier_mass),kind=8)
     call MPI_ALLREDUCE(local_count,global_count,1,MPI_INTEGER8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
     if(myid==1)write(*,'(A,I0,A,ES24.16,A,I0)') &
          'LM_GATE_TIER: level=',ilevel,' mass=',tier_mass,' count=',global_count
  end do
end subroutine lagmusic_gate_mass_probe

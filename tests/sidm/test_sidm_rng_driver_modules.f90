! Fixture modules for the actual SIDM wrapper and the production RNG.
module amr_parameters
  implicit none
  integer,parameter :: twotondim=8,amr_block_size=64
  logical :: sidm_estimator_diagnostics=.false.
end module amr_parameters

module amr_commons
  use omp_lib
  use amr_parameters
  implicit none
  integer,parameter :: dp=kind(1.0d0)
  integer,parameter :: MPI_IN_PLACE=0,MPI_INTEGER=1,MPI_DOUBLE_PRECISION=2
  integer,parameter :: MPI_SUM=1,MPI_MAX=2,MPI_COMM_WORLD=0
  integer :: ncpu=1,myid=1,numbtot(1,1)=2,headl(1,1)=1,numbl(1,1)=2,next(2)=(/2,0/)
  logical :: sidm=.true.,verbose=.false.,sidm_inelastic=.false.
  real(dp) :: sidm_Pmax(1)=0.0d0,sidm_fdiss=0.0d0,sidm_vrel_max=0.0d0
end module amr_commons

module pm_parameters
  implicit none
  integer :: iseed=12345
end module pm_parameters

module pm_commons
  implicit none
  integer :: localseed(4)=-1
end module pm_commons

module sidm_test_capture
  implicit none
  integer :: call_number=0,workers(2)=0,visits(2)=0
  real(kind(1.0d0)) :: draws(4,2)=0.0d0
end module sidm_test_capture

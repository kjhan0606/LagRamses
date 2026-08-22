module poisson_parameters
  use amr_parameters

  ! Convergence criterion for Poisson solvers
  real(dp)::epsilon=1.0D-4

  ! Maximum number of fine-level multigrid V-cycles.  Keep the historical
  ! default, but allow difficult production meshes to request more cycles.
  integer :: maxiter_fine=10

  ! Restarts use the standard predictor by default.  A valid checkpoint
  ! marker is necessary but not sufficient for the restored-phi warm start;
  ! this explicit diagnostic opt-in is also required.
  logical :: restart_phi_warm_start=.false.

  ! Production jobs can fail closed when fine MG reaches maxiter_fine while
  ! still above epsilon.  The default remains backward compatible.
  logical :: abort_on_mg_nonconvergence=.false.

  ! Type of force computation
  integer ::gravity_type=0

  ! Gravity parameters
  real(dp),dimension(1:10)::gravity_params=0.0

  ! Maximum level for CIC dark matter interpolation
  integer :: cic_levelmax=0

  ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  integer :: cg_levelmin=999

  ! Gauss-Seidel smoothing sweeps for fine multigrid
  integer, parameter :: ngs_fine   = 2
  integer, parameter :: ngs_coarse = 2

  ! Number of multigrid cycles for coarse levels *in safe mode*
  !   1 is the fastest,
  !   2 is slower but can give much better convergence in some cases
  integer, parameter :: ncycles_coarse_safe = 1

end module poisson_parameters

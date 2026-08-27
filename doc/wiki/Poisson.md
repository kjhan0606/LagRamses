

# Poisson Parameters

This namelist, `&POISSON_PARAMS`, is used to specify runtime parameters for the Poisson solver. It is used only if `poisson=.true.` or `pic=.true.`


The lagRamses project baseline uses the FFTW3 MPI/OpenMP direct solver on a
fully uniform base level (`use_fftw=.true.`).  It solves the periodic discrete
Poisson equation without an iterative residual tolerance.  MG remains
available by explicitly setting `use_fftw=.false.` for diagnostic comparisons,
and it is still used on non-uniform refined AMR levels.

On refined levels, conjugate gradient (CG) and multigrid (MG) are available.
Unlike CG, MG has an initialization overhead cost at every call, but is more
efficient on large levels with few holes.

MG can be used on refined levels in conjunction with CG. After the FFTW
uniform-level dispatch, `cg_levelmin` selects the iterative solver as follows:

* Fully uniform levels are solved directly with FFTW when `use_fftw=.true.`
* Remaining coarse or non-uniform levels are solved with MG
* Refined levels with *l* < `cg_levelmin` are solved with MG
* Refined levels with *l* >=  `cg_levelmin` are solved with CG

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `gravity_type`      | `int`  | 0     | Type of gravity force. Possible choices are: `gravity_type=0` self-gravity (Poisson solver); `gravity_type>0` analytical gravity vector; `gravity_type<0` self-gravity plus additional analytical density profile
| `epsilon`           | `real`  | 1e-4  | Stopping criterion for the iterative Poisson solver: residual 2-norm should be lower than `epsilon` times the right hand side 2-norm.
| `gravity_params`    | `real array`  | 0.0, | Parameters used to define the analytical gravity field (routine `gravana.f90`) or the analytical mass density field (routine `rho_ana.f90`).
| `cg_levelmin`       | `integer`  | 999 | Minimum level from which the Conjugate Gradient solver is used in place of the Multigrid solver.
| `cic_levelmax`      |	`integer`  | 0	 | Maximum level for CIC dark matter interpolation (default `cic_levelmax=nlevelmax`)
| `use_fftw`          | `logical`  | `.true.` | Use the FFTW3 direct Poisson solver on fully uniform levels. Requires the default FFTW-enabled build; set `.false.` only to request MG explicitly.

For the Paper-Ib $256^3$ DMO restart on `grammar-debug` (4 MPI ranks), OpenMP
CIC linked-list partitioning reduced the four-thread density-deposition phase
from 2.192 s to 1.275 s.  Coarse-step times were 60.500, 15.620, and 9.455 s
with 1, 4, and 8 threads per rank, respectively.  The corresponding parallel
efficiencies were 96.83% and 79.98%; going from 4 to 8 threads reduced the
wall time by 39.47%.  At eight threads, atomic mesh contention increased the
density-deposition phase to 1.558 s, identifying the current scaling limit.
The reported final conservation state was identical in all three runs.

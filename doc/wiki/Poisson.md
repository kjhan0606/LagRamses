

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

The default OpenMP CIC path uses atomic-free spatially colored scatter.  A
particle grid can touch only its $3^3$ neighboring-grid stencil.  Grid
coordinates modulo four therefore define 64 phases in which no two source
grids can update the same target cell, including across the periodic
power-of-two mesh.  One-thread and sparse-level paths bypass the color setup.
The former atomic implementation remains available for controlled regression
tests with `make CIC_ATOMIC_SCATTER=1`.

An exclusive-node Paper-Ib $256^3$ DMO restart benchmark on `grammar-debug`
compared the two implementations with identical executables apart from that
build switch.  At 8 MPI ranks x 8 OpenMP threads, the colored path reduced
the final timer total from 5.903 to 4.470 s (1.32x) and the `rho` phase from
0.827 to 0.310 s (2.67x).  At 16x4 the corresponding values were
4.094 to 4.046 s and 0.325 to 0.277 s; at 32x2 they were 4.437 to
4.444 s and 0.267 to 0.252 s.  Thus 16x4 was the fastest tested 64-core
layout, while 32x2 had the highest colored-path OpenMP efficiencies among
those layouts (89.3% for the full timer and 69.2% for `rho`).  A separate
fixed-four-rank test measured full-timer efficiencies of 87.5% and 76.9% at
four and eight threads, with the eight-thread run still reducing wall time by
43.1% relative to four threads.  All runs printed the same final conservation
state.  `Time elapsed since last coarse step` includes restart preparation;
these comparisons use the final timer `TOTAL` and report `rho` separately.

# Void-only V-web refinement

`void_web_refine` adds two minimum mesh levels for cosmological void zoom
simulations.  It is disabled by default.  When disabled, no state array is
allocated, no V-web calculation is run, and no additional MPI exchange is
performed; the standard RAMSES refinement path is unchanged.

The feature is separate from the legacy geometry-based `void_refine` option.
It does not require a void centre, radius, or wall trajectory.

The isolation also covers hydro IC ingestion.  Only with
`void_web_refine=.true.` does the reader require a velocity field for the
initial V-web state.  It prefers the two-fluid baryon files
`ic_velb{x,y,z}`.  If they are absent but `ic_velc{x,y,z}` exist, the latter
are used for the initial gas velocity and an explicit single-fluid IC warning
is printed.  With the feature disabled, the historical lagRamses IC file
selection remains unchanged.

## Refinement rule

At `void_web_env_level` (currently required to equal `levelmin`), the code
constructs the symmetric velocity-shear tensor

```text
Sigma_ij = -(d v_i / d x_j + d v_j / d x_i) / (2 a H / h)
```

from density-weighted gas velocities and centred differences.  RAMSES
cosmological units reduce the denominator to `aexp*hexp` after differentiating
the code velocity with respect to the code coordinate.  The largest symmetric
eigenvalue is the local directional-compression statistic.  A cell enters the
wall state at `void_web_lambda_on` and leaves it below
`void_web_lambda_off`; the two thresholds prevent one-step mesh flicker.

The environment state is inherited by all AMR descendants through their
`father` links.  Scoped cells receive `void_web_base_level`; compressive cells
receive `void_web_wall_level`.  These are floors, not caps.  Ordinary mass,
Jeans, shock, and galaxy refinement criteria remain OR-active above the wall
level, so a genuine halo is not suppressed merely because it lies in a wall.

## Scope choices

- `void_web_scope_ivar=-1`: whole box.  This is useful only for small tests and
  is usually too expensive for production.
- `void_web_scope_ivar=0`: the static GRAFIC `ic_refmap` represented by
  `cpu_map2`.
- `void_web_scope_ivar>0`: a passive-scalar mass fraction.  The cell is in the
  target scope when `uold(ivar)/rho > void_web_scope_cut`.

The V4 VoidSim GRAFIC hierarchy writes the target mask to both `ic_refmap` and
`ic_pvar_00001`.  With the current `NVAR=11` hydro build, `ic_pvar_00001` is
hydro variable 6.  Using `void_web_scope_ivar=6` is preferred because the mask
is advected with the gas instead of remaining fixed in Eulerian coordinates.

## V4 pilot configuration

The following block implements an L12 void floor and an L13 compressive-wall
floor.  The lambda thresholds are provisional pilot values and must be checked
against the offline V-web/T-web diagnostics before a production run.

```fortran
&REFINE_PARAMS
  void_web_refine=.true.
  void_web_env_level=9
  void_web_base_level=12
  void_web_wall_level=13
  void_web_scope_ivar=6
  void_web_scope_cut=0.5
  void_web_lambda_on=0.6
  void_web_lambda_off=0.4
  void_web_update_interval=4
/
```

For the 512 Mpc/h V4 box, level 9 cells are 1 Mpc/h wide and the centred
gradient spans 2 Mpc/h.  The current implementation has no additional Gaussian
smoothing.  The pilot must therefore record the flagged volume, wall
persistence, shock-leading behaviour, and cost before the thresholds are
adopted for production.

The feature requires a three-dimensional, periodic, cosmological hydro run.
It is intentionally rejected for DMO and non-cosmological simulations.  One
integer cell-state array is allocated only while enabled, costing approximately
`4*(ncoarse + 8*ngridmax)` bytes per rank in 3D.  Each state update uses one
integer ghost exchange and two packed MPI reductions.

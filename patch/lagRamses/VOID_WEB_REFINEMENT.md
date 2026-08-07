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
Jeans, and galaxy refinement criteria remain OR-active above the wall level,
so a genuine halo is not suppressed merely because it lies in a wall.

When the V-web option is active, the hydro gradient, velocity-jump, and
kinetic-flux triggers are restricted to the advected passive-scalar scope and
can open levels only through `void_web_hydro_max_level`.  A negative input
selects `min(levelmax, void_web_wall_level+1)`.  The finite ceiling prevents
cold, supersonic void gas from refining to the galaxy scale.  It does not cap
the mass or Jeans criteria.  Enabling a hydro trigger in V-web mode therefore
requires `void_web_scope_ivar` to select a passive scalar.  With
`void_web_refine=.false.`, these restrictions are bypassed and the historical
all-domain hydro refinement path is unchanged.

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

The following block implements an L12 void floor, an L13 compressive-wall
floor, and scoped hydro triggers through L14.  The thresholds are calibration
values rather than final production choices.

```fortran
&REFINE_PARAMS
  void_web_refine=.true.
  void_web_env_level=9
  void_web_base_level=12
  void_web_wall_level=13
  void_web_hydro_max_level=14
  void_web_scope_ivar=6
  void_web_scope_cut=0.5
  void_web_lambda_on=1.0
  void_web_lambda_off=0.8
  void_web_update_interval=4
  err_grad_d=0.5
  err_jump_u=32.0
  ekin_flux_refine=-1.0
  d_keflux_max=0.2
/
```

For the 512 Mpc/h V4 box, level 9 cells are 1 Mpc/h wide and the centred
gradient spans 2 Mpc/h.  The current implementation has no additional Gaussian
smoothing.  The pilot must therefore record the flagged volume, wall
persistence, shock-leading behaviour, and cost before the thresholds are
adopted for production.

The high-redshift V4 audit found that the earlier
`err_jump_u=0.5, ekin_flux_refine=4` example selected 97.9 per cent of the
level-9 scope.  At level 13, even `err_jump_u=2, ekin_flux_refine=64` selected
99.34 per cent of active scoped cells.  The bulk-Mach trigger remains active
in smooth cold flow and is therefore disabled in the staging template.  A
level-13 audit of `err_grad_d=0.3, err_jump_u=16` selected 5.79 per cent of
active scoped cells.  AMR expansion raised the realized level-14 grid count
to 31.4 per cent of the audited active scope, so that choice was rejected.
Raising the jump threshold to 32 reduced the raw union to 2.49 per cent.  A
five-step online test completed without a grid-pool, NaN, negative-energy, or
Poisson-convergence error.  Its level-14 count nevertheless grew from 833,212
to 1,259,622 and 1,956,627 grids, and the final count used 65.2 per cent of the
grid pool.

A matched density-gradient sweep selected `err_grad_d=0.5` for the next
preflight.  At level 13, density or jump 32 then selected 104,112 cells, or
1.385 per cent of the active scope.  Only 349 cells passed both conditions.
The matched online run contained 667,957, 831,812, and 1,230,549 level-14
grids at main steps 3 through 5.  The final count was 37.1 per cent below the
`err_grad_d=0.3` result and used 41.0 per cent of the grid pool.  Launcher
time decreased by 17.3 per cent and the maximum RAMSES memory diagnostic
decreased by 4.6 percentage points.  The run completed with all 54 explicit
NaN checks at zero and without a negative-energy, Poisson, or grid-pool error.
The combination `err_grad_d=0.5, err_jump_u=32` is therefore the preferred
high-redshift preflight setting.  It is not a production calibration.

A matched density-only control set `err_jump_u=-1` and changed no other
parameter.  Its level-14 counts were 123,663, 656,781, and 942,903 at main
steps 3 through 5.  The final count was 23.4 per cent below the jump-32 result
and used 31.4 per cent of the grid pool.  Launcher time decreased by 5.1 per
cent and the maximum RAMSES memory diagnostic decreased by 3.8 percentage
points.  The central slice placed 352 of 376 level-14 pixels inside the
advected scope.  The large late increase shows that the density-gradient
condition drives continuing mesh growth even without the jump trigger.  The
jump trigger advances the onset and adds complementary cells, so its physical
value must be judged from low-redshift shock, thermal, and velocity profiles.
If the added cells do not improve those profiles, production should use
`err_jump_u=-1`.  If they do, the trigger needs an additional convergent-flow
or pressure-jump gate.  The five-point PPM prototype cannot provide the AMR
comparison yet because it intentionally rejects AMR and self-gravity.

The feature requires a three-dimensional, periodic, cosmological hydro run.
It is intentionally rejected for DMO and non-cosmological simulations.  One
integer cell-state array is allocated only while enabled, costing approximately
`4*(ncoarse + 8*ngridmax)` bytes per rank in 3D.  Each state update uses one
integer ghost exchange and two packed MPI reductions.

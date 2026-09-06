# Native IR coarse–fine transport and populated-AMR execution

This bundle connects IR transport across refinement interfaces, then verifies
it in actual CPU-hydro/CUDA-SNRT runs. It supersedes the missing coarse–fine
connection in `snrt_level_synchronization_2026-09-07.md`. It does **not** declare
the whole dust/feedback model production-ready.

## Implementation

- Fine cells own the coarse–fine face flux. Their incoming IR uses the current
  coarse trial, exchanged through the level-minus-one RAMSES halo. Opposite
  directional energy-density corrections are scaled by fine/coarse volume
  (1/8 in 3D), summed by the reverse halo, and retained in the owner trial.
- Coarse trial donors are refreshed after every IR CFL substep. Coarse faces
  adjoining finer cells block both inflow and outflow in the coarse solve;
  they are neither vacuum boundaries nor a second application of the flux.
- Fine IR, coarse IR and dust material commit only after the enclosing RT
  transaction succeeds. Missing coarse donors, nonfinite/negative coarse
  energy, invalid topology or a failed material solve reject the transaction.
- Ranks with no local coarse face still participate in the globally selected
  exchange protocol. The retained test has fine/coarse cells on both ranks;
  this does not establish every completely empty-rank configuration.
- Populated-grid coarsening exposed a second missing connection: `defrag`
  changed hydro cell IDs but not the SNRT slot map. `snrt_regrid_defrag` now
  applies that exact old/new grid map to primary and IR identities together.
  The large radiation arrays stay in their slots. Retired slots are cleared
  and parked on free cells, preserving unique nonzero checkpoint identities
  and subsequent pool reuse. Duplicate destination IDs reject before rebinding.

VPATH matters here: the actual defrag source is
`patch/lagRamses/load_balance.kjhan.f90`, not its `patch/cuRamses` fallback.
`make -n -W` verified the selected source. Its new module dependency is inside
`ifeq ($(SNRT),1)`; a legacy `SNRT=0` dry run does not pull SNRT objects.

Version-3 AMR admission remains **reference-control only**, requiring the
existing explicit nonproduction opt-in. Noncosmological operation, `nremap=0`,
HDF5 output/build and HDF5 restart input remain required. `SNRT_RT_LEVEL` must
be unset for AMR because all participating levels must advance. No simulation
namelist schema was added; existing namelist generator schemas need no change.

## Retained native evidence

Root: `/gpfs/kjhan/LRD_JWST/.dust-ir-coarse.ppElta`.
Build: `SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1`, Intel MPI/ifx, NVAR=30.
Targets `snrt_regrid_native ramses`; default executable left untouched.

Initial simulation SHA256:
`675e9da71c801bafc32de73ce341e850f4c7ce1a22440727ad68823c2fcda483`.
Final simulation `rebound/ramses3d`:
`86cf354bd03aa3ed1e93f647796e598994ae59e9b59e8874acc4622812e591aa`.
Final native test `rebound/ramses_regrid`:
`dc9898e3ebc3622befeb992502e82571292fd99e776ea7945d7fd5953fb37e9b`.
The final post-build Makefile change only conditions the new dependency on SNRT.

The native test, linked against the production objects, reports:

- `SNRT_NATIVE_IR_COARSE_FINE_PASS`: four fine faces sharing a coarse face,
  both propagation signs, three CFL substeps, exact analytic upwind values,
  volume-weighted closure within 2e-14 in normalized units, persistent states
  unchanged before commit and after rejection, no double coarse flux.
- `SNRT_NATIVE_DEFRAG_PASS`: cyclic overlapping grid renumbering preserves
  primary/IR payloads exactly; duplicate identity rejection is nonmutating.
- Existing `SNRT_NATIVE_REGRID_PASS` and `SNRT_NATIVE_LEVEL_UPLOAD_PASS`.

`JAX_PLATFORMS=cpu simulation/snrt/.venv/bin/python
simulation/snrt/tests/dust_ir_transport.py --native` passed GNU and Intel
native/reference differentials (largest energy error 1.11e-15, accumulated
closure 7.407e-10), including previous rollback and halo checks.

## Actual executions

Every case retains `effective.nml` and `run.log`; runners record the executable,
contracts and MPI environment. No preexisting output was overwritten.
The synthetic v3 opacity/thermal fixture is NOT a physical-data approval.

| Cases | Evolution and output | Result |
| --- | --- | --- |
| `pattern1`, `pattern2` | Prior nonuniform primary/IR HDF5 restart; fixed 8³; 2 further steps; one 4.2 MiB dump each | Coordinate-matched primary/IR, gas energy, dust mass/energy at levels 2/3 bitwise equal to the previous implementation and between 1/2 ranks |
| `amr_cpu1`, `amr_cpu2` | Fresh level 3/4; 2 coarse steps, two fine subcycles; one ~17 MiB dump each | 296 coarse leaves + 1728 fine leaves; volume and gas mass exactly 1; nonzero IR; 1/2-rank fields bitwise equal |
| `rebound/regrid1`, `rebound/regrid2` | Copy of the same `amr_cpu2` checkpoint; 2 further coarse steps; refinement diameter 0.5→0.25; one 7.8 MiB dump each | 152 grids actually deleted; 448 coarse leaves + 512 fine leaves; no zeroed live IR cells; volume/gas mass exactly 1; 1/2-rank fields bitwise equal |

Output policy was audited on the copied files: `noutput=1,aout=2,tout=1e30`,
scheduled epochs not reached; `fbackup=1000000`. Fresh AMR uses `foutput=2`,
restarts use `foutput=4,nrestart=1,nstepmax=4`. Free space before launch was
173 TiB. Hydro, periodic boundaries, gamma=5/3 and units were explicit; no
stars/sinks/AGN, gravity, gas cooling, cosmology or load balancing. The accepted
AMR runs use CPU hydro and CUDA primary transport, with MPI-local-rank GPU
binding in `rank_exec.sh`. IR itself is the native FP64 CPU operator.

Physical closure sums **leaves only**, with volume `2^(-3*level)` and the native
normalized S8 angular weights. IR records use group-fast reshape `(cell,80,2)`;
dust energy is uold_18 times `1.66e-24*(3.0857e21/3.1556926e13)^2`.
Initial mean material+IR energy is 2e-23 erg/cm³. After fresh AMR it is
1.9999999998400766e-23 (relative error 7.9962e-11). After restart and actual
coarsening it is 1.9999999996504576e-23 on both rank counts: relative change
9.48095e-11 from the input checkpoint.

## Failures retained, not counted as passes

- `amr1`: initial fixture incorrectly used `ivar_refine=0`, which requested
  an external `ic_refmap`. Startup stopped; corrected fixture uses internal
  density/geometry refinement with `ivar_refine=-1,m_refine=0`.
- `amr1b`: CUDA hydro completed, but a nominally uniform state gained
  0.2925% gas mass and changed velocities. Total dust+IR changed 0.1458%.
  CPU hydro preserves exactly the uniform gas mass/state in the same setup.
  This isolates a CUDA-hydro AMR discrepancy; its underlying cause is not
  diagnosed here. Do not qualify this CUDA-hydro AMR profile from the RT pass.
- `amr2`: concurrent/unbound GPU use exhausted memory. The later accepted
  CPU-hydro cases explicitly bind MPI local ranks to separate idle GPUs.
- Initial `regrid1`, `regrid2` completed but their snapshots lost IR and
  disagreed across ranks because defrag did not rebind radiation identities.
  They are superseded by the final `rebound/regrid*` evidence, not erased.

## Remaining scope

AMR **reference** execution, conservative IR faces, actual coarsening and
defrag are demonstrated. Moving/nonuniform gas, additional populated refinement
creation histories, physical boundaries, general empty-rank behavior,
load-balance migration, cosmology, production performance and CUDA-hydro AMR
conservation remain unqualified. Actual production opacity/heat-capacity/SED
inputs and live AGN pending-energy MPI/restart ownership remain separate open
requirements. Continue high-level physics work; do not interpret this bundle
as completion of all feedback/dust readiness.

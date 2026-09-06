# Live IR MPI halo connection on a fixed mesh

This extends `dust_live_ir_execution_2026-09-07.md` from one rank to distributed
fixed meshes. It does not remove the fixed-level, noncosmological, `nremap=0`,
HDF5 build/output/restart requirements or approve physical reference inputs.

## Implementation

The IR operator accepts paired optional ghost energies and a face-to-ghost
index. A ghost face must have no local neighbor; indices, dimensions, finite
values and nonnegative energies are validated before state mutation. Local
faces retain the reciprocal-neighbor check. MPI outgoing minus incoming
energy is a signed `interface_erg` ledger, distinct from `escaped_erg` at a
physical boundary. It is included in local energy closure and cancels across
matching rank faces. The source-present convergence normalization is unchanged;
source-free receiving domains include incoming boundary energy in their scale.

The live adapter reuses `make_virtual_fine_dp`, the same RAMSES halo communicator
used by the existing primary adapter, without downcasting IR to FP32. Each
IR substep exchanges the **current trial** field, not the original persistent
field. No all-gather of the global radiation field is introduced. Ghost buffers
scale with MPI boundary faces; one cell-indexed scalar scratch array is reused
for each IR group/direction.

All ranks perform admission-error reduction before entering halo exchange.
Substep counts are made collective. Every rank participates in each exchange,
including one without local MPI faces; errors after local emission solves are
also reduced before another substep can start. The adapter handles empty local
arrays, but this bundle does not claim to qualify the entire primary RT/hydro
path with empty ranks. Persistent IR and material still commit only after the
primary transaction succeeds. Coarse/fine and unmapped interfaces are rejected
rather than treated as vacuum. The single-rank startup restriction is removed;
other restrictions are retained.

## Validation

The existing native IR smoke now compares a two-cell periodic transport problem
with two separate one-cell domains fed the same frozen ghost states. Split and
whole updates agree; interface ledgers cancel with nonzero individual fluxes,
physical escape is zero, and a negative ghost payload is rejected without a
state update. The existing GNU/Intel native differential, nonlinear rollback,
weak/stiff emission and background-roundoff tests still pass:

`JAX_PLATFORMS=cpu simulation/snrt/.venv/bin/python simulation/snrt/tests/dust_ir_transport.py --native`

An initial edit broadened the source-present normalization to old radiation
energy. The existing accumulated-balance test caught that change; the original
source-present normalization was restored before execution qualification. The
tolerance was not relaxed to obtain a pass.

Retained real-run root: `/gpfs/kjhan/LRD_JWST/.dust-ir-mpi.HJvK2V`.
Binary: `ramses3d`, SHA256
`acac362ebbb8daf10f668f8c0625c026e039fb6ed3722bc72c42529c4d79f2da`.
Built with the already-active Intel MPI/ifx live profile:
`make -s -C bin -j4 SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1 EXEC=/gpfs/kjhan/LRD_JWST/.dust-ir-mpi.HJvK2V/ramses`.

The root's `run_case.sh` records environment and launcher. It uses one OpenMP
thread, GPUs 1/2, shared-memory Intel MPI, the canonical reference primary and
secondary contracts, and `dust_native_reference_control_v3.nml` with explicit
NONPRODUCTION opt-in. Initial conditions and physics are the tracked
`simulation/snrt/config/dust_live_ir_smoke.nml`: stationary gas, 20 K dust,
512 periodic cells; no stars/AGN, gravity or gas cooling.

- `fresh2`: fresh two-rank run, two steps, 256 active cells per rank, successful
  RT/IR commit and one 4.2 MiB HDF5 dump. Local closure errors <9e-11.
- `pattern2`: copy of that dump, advanced two further steps on two ranks.
  Primary group-1 directional bins are set to FP32 `1e-11*f`, where
  `f=1+0.75*x_grid+0.15*y_grid+0.10*z_grid+0.01*child_index` (child index 0..7).
  IR is multiplied by `f*(0.5+0.01*direction_index)` (direction index 0..79).
  Thus both space and direction vary; f spans 1.125..1.945.
- `pattern1`: identical copied input, two-rank checkpoint restored and evolved
  on one rank. Both runs report HDF5 restore, RT/IR commits and completion.

After sorting grid coordinates, all 452608 primary/IR entries, all 512 dust
mass/energy entries and all 512 gas-energy entries match **bitwise** between
`pattern1/output_00002` and `pattern2/output_00002`. Gas energy and dust mass
also match the input exactly. Summed lost primary energy is
`4.116128178964145e-20 erg/cm3`; summed IR+material gain is
`4.1161281238488174e-20`, discrepancy `1.3390090121236643e-8` relative to the
loss (<2e-6 allowance for the FP32 primary state). Both rank counts yield the
same ledger. Per-step two-rank local IR closure errors are <8.7e-10.

Effective namelists are preserved in each case. All keep `noutput=1`,
`aout=2`, `tout=1e30`, `fbackup=1000000`; the scheduled times are not reached.
`fresh2` uses `nstepmax=2,foutput=2`; the pattern cases use
`nrestart=1,nstepmax=4,foutput=4,informat='hdf5'`. Each writes only one new
4.2 MiB dump. Preflight free storage was 173 TiB. Only copied test checkpoints
were seeded; original outputs are preserved.

## Remaining scope

This demonstrates same-level IR MPI transport and two-rank/one-rank HDF5
continuation, not arbitrary rank changes: the previously observed 1-to-2
hydro remap failure is not repaired by this bundle. Refinement/coarsening,
coarse-fine IR fluxes, runtime slot migration, cosmological radiation changes
and approved physical dust inputs remain open. Group-by-direction scalar halo
exchange follows the existing communicator and is correct at this scale; its
large-production performance has not been established. No new simulation
namelist fields were introduced, so generator schemas are unchanged.

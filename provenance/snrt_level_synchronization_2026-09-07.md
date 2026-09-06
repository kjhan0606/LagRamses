# Radiation level synchronization before coarse-fine coupling

The actual `amr_step` advances finer levels before the current level's RT
solve. The primary coarse-flux correction iterates existing state slots; a
coarse leaf without a slot cannot receive that correction. Also, before this
bundle, the new radiation restriction ran only when destroying a grid, leaving
the retained coarse representation stale while children evolved.

Two runtime connections now address these prerequisites:

- `snrt_state_sync_level` runs before recursion, preparing coarse leaf receiver
  slots before a fine step can deposit flux in them.
- `snrt_regrid_upload(ilevel)` runs after the level RT call, restricting the
  completed fine level to `ilevel-1` without clearing children. Photon/IR
  densities use volume averages; ion fractions use density weights. The merge
  implementation is shared with grid-deletion restriction.

Each fine-grid owner computes its oct average. RAMSES' reverse halo exchange
delivers it to the owning parent rank. Every oct must have exactly one receiver:
global source/receiver counts and per-parent multiplicity are checked before
state writes. Received payloads are validated collectively, and child arrays
are retained. Scalar scratch buffers are reused per payload component; this is
not a global radiation all-gather. Large-production synchronization cost still
needs measurement.

Be careful with the two APIs: hydro `upload_fine(ilevel)` refreshes split cells
**on** ilevel from ilevel+1, whereas radiation upload takes the completed child
level and updates ilevel-1. Comments at the call site now state this explicitly.

The Makefile dependency for `amr_step.jaehyun.o` includes the regrid module.
VPATH/compile-path verification selected
`patch/lagRamses/amr_step.jaehyun.f90`, not a fallback implementation.

## Evidence

Retained root: `/gpfs/kjhan/LRD_JWST/.snrt-level-sync.YZE1Hh`.

Build profile: `SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1`, Intel MPI/ifx. Targets
`snrt_regrid_native ramses`, with `EXEC` pointing at the retained root's `ramses`.
The default executable was not overwritten. Native test SHA256:
`147f8e9be6a08d1eb53cd369c780cc16892fca263de1d422e8e29223028ed37a`.
Simulation SHA256:
`97b6f9c5ef25931ab1fb11da6ab5d4983c1b54635a63505aebdb165baf994beb`.
Only an explanatory call-site comment changed after these builds.

The extended native test prints `SNRT_NATIVE_LEVEL_UPLOAD_PASS` and retains
the previous `SNRT_NATIVE_REGRID_PASS`. A populated level-2 oct restricts into
its level-1 parent with density-weighted ions and averaged IR, while every
child's original radiation is checked unchanged. Existing regrid, invalid
input and capacity-growth checks also pass.

`pattern1` and `pattern2` run the prior nonuniform primary/IR checkpoint on one
and two ranks for two further steps. Their effective namelists are retained:
fixed 8^3 hydro mesh; `nrestart=1,nstepmax=4,nremap=0`; no stars/AGN, gravity or
gas cooling. Output policy is `noutput=1,aout=2,tout=1e30,foutput=4,
fbackup=1000000` with HDF5 input/output. Scheduled times are not reached. Each
case writes one approximately 4.2 MiB dump; free space before launch was
173 TiB. Only checkpoint copies were used.

After coordinate matching, one/two-rank results are bitwise identical for:

- 56,576 level-2 parent record values (64 cells times 884 entries);
- 452,608 level-3 leaf record values (512 cells times 884 entries).

Every level-2 parent has present state. For all 64 parents, the primary payload
equals the mean of its actual eight children rounded to the stored FP32
representation, and IR equals the FP64 child mean exactly. The level-3
primary/IR records, gas energy, dust mass and dust energy are bitwise unchanged
from the previous implementation's two-rank result. Physical totals must still
sum leaves, not the extra coarse representations as well.

## Not yet complete

This establishes actual per-step coarse representations and pre-recursion
receiver allocation, not IR coarse-fine transport/reflux or an evolving-AMR
run. The v3 fixed-level restriction remains. Coarse-fine IR face flux ownership,
its conservative commit alongside the fine step, live regridding qualification,
load-balance migration, cosmology and physical input approvals remain open.
No simulation namelist schema changed; generator updates were not required.

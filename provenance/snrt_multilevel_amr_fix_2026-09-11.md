# Kind7 actual-source multilevel AMR repair

Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`.
The preceding uniform-source bundle was committed/pushed as `9b2c289` before
this operator-requested bounded AMR repair. No new physics model or gate chain.

## Causes and bounded changes

1. `snrt_dust_prepare_cell_optical_depth` and `snrt_dust_receiver_stage`
   rejected a consistent zero-cell batch. A fully covered coarse level, or an
   MPI rank without leaves, is not an invalid hydro state. Remove only the
   `nc<1` rejection; retain group, shape, finite, positivity and timestep checks.
   Do not return early from the driver: MPI reductions/exchanges, source
   bookkeeping, transport/IR transactions and subsequent restriction remain.
2. The partial-refinement run then exposed wrong face ownership. A leaf's
   `cpu_map(cell)` selects ownership of a prospective child oct, whereas its
   current hydro grid is owned by `cpu_map(father(IGRID_OF(cell)))`. These can
   differ across a Hilbert boundary. The old classifier treated a legitimate
   reception-grid leaf as an unmapped local leaf (rank1 cell131, neighbor10,
   local hash lookup absent, slot0). Use validated current-grid ownership in
   both classifier paths and their diagnostics. Keep local hash matches and
   coarse-to-fine classification precedence. Invalid ownership is not vacuum.

This is a correction to existing topology semantics, not a new communicator,
particle-tree rewrite, hash-size increase or physical coupling prescription.
Makefile/VPATH, namelist schema and both generators are unchanged by this repair.

## Reproduction and execution identity

Evidence root: `.chimes-amr.vHVvTY`; source-build root:
`.chimes-sources.r1TaD8`. Existing Makefile invocation:

```
make -f ../bin/Makefile -j1 SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1 \
 USE_CUDA=0 USE_FFTW=0 NENER=0 NVECTOR=32 \
 CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.chimes-transition.YdVRvD/chimes \
 SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install \
 EXEC=ramses_amr_owner ramses
```

Receiver-only binary `ramses_amr_fix3d` SHA256:
`2d456eedf17cd69cdef47dd5ecb8f2d2708b26ef7b1ecf4b7261a22a82c1ac38`.
Both-fixes binary `ramses_amr_owner3d` SHA256:
`627f050e46cab327e2b655f75cd6612262bfd0ef0b7d83c30dc4bef71a7fd4ae`.
Build logs: `build.log`, `build-owner.log`.

Effective inputs under the evidence root: `refined/physical.nml` (4 steps),
`mixed/physical.nml` (failed before owner fix), `mixed-owner/physical.nml`
(same 2-step input after owner fix), with individual `ic_sink` and `run.log`.
Source the preceding bundle's `environment.sh`, **unset `SNRT_RT_LEVEL`**, then
MPI2/OMP2. Kind7, actual PARSEC stars plus reference non-MAD Bondi/AGN,
noncosmological periodic NENER0 hydro and evolving C/silicate dust; no fake
photon seeding. `nsubcycle=1`, `nremap=0`, initial normal load balancing enabled.
Lmin2/Lmax4 with `m_refine=0,-1` gives L3 refinement and no L4 cells, preserving
the L/4 sink-cloud radius. Partial case limits refinement to a small corner
box (`r_refine=.26`, centers `.125`, exponent10, `nexpand=0`).

Outputs were audited before each launch: `noutput=1,aout=2,tout=1e30` outside
the short run, `foutput=2,fbackup=1000000`; free storage about98TiB. The full
case wrote two HDF5 files, about108MiB each (slightly above the initial100MiB
per-file estimate but within the0.3GiB run allowance). No production outputs.

## Measurements

- Existing native receiver smoke, with empty-batch cases: before fix exits20;
  after fix passes. Bad shapes and negative timesteps still reject empty
  batches. Compiled with `ifx -g -check all`; `receiver-before.log` and
  `receiver-after.log` retained.
- Fully covered case: 512 L3 leaves, zero L2 leaves; MPI grid counts32/32 at
  L3. Four steps complete, eight RT commits including four empty L2 commits.
  Maximum IR balance residual2.0481e-10. Snapshot gas+solid element error
  <=8.912e-16, charge<=5.427e-16, global gas+stars+BH mass including radiated
  rest-mass loss2.147e-16 relative. Final1495 stars, positive actual source
  gas/grain absorption, finite positive temperatures1605.2--1.187e7K.
- Mixed before owner fix: L3 commit succeeds with only one rank owning its
  single oct; L2 then rolls back on unmapped faces. Exit0 from `clean_stop`
  is **not** a successful run. This failure is retained, not overwritten.
- Mixed after owner fix: two steps complete, four RT/IR commits, no unmapped
  faces or rollback. During both steps L3 has one oct on rank1 and no local
  leaves on rank2; L2 has63 leaves across3/5 owned grids. Both levels receive
  actual stellar radiation and gas/grain absorption in the second step.
  The terminal mesh derefines to64 L2 leaves (no stored L3); do not describe
  that final snapshot as a persistent mixed mesh. Its volume-weighted budgets
  therefore also cover restriction/derefinement. Element error8.471e-16,
  charge4.912e-16, global mass2.147e-16, IR balance4.5432e-10 maximum;
  141 stars, thermal-mode BH feedback, temperature1605.6--19151.2K. Single
  HDF5 file61,895,296bytes; full metrics in `mixed-results.json`.

The one-off evaluator reuses the preceding bundle's element/charge equations,
selects `son_flag==0` across all stored levels and weights mass by `2**(-3L)`.
Covered parents are not double-counted. No new generic validation framework.
Detailed metrics and input/log/build identities are retained in the evidence
root. This short case does not claim cosmological, GPU, changed-MPI restart,
subcycling or resolved-jet convergence coverage.

## Driver disposition and retention

PASS for the recorded zero-leaf AMR failure and the newly reproduced mixed
AMR face-ownership failure. No blanket claim for every physics/AMR option.
No new state format, output field or namelist option; no extra restart matrix
is introduced for this bounded fix. Prior uniform-source restart evidence
remains scoped to that preceding bundle. No additional routine external audit.

Three evaluated output directories were removed under the standing raw-output
retention directive; none was required for an active restart. Retained:
effective inputs, logs, binaries, build logs, input SHA256s, per-HDF5 SHA256s
and compact `refined-results.json`/`mixed-results.json`, plus the one-off
evaluator. Deleted targets (relative to evidence root):

- `refined/output_00001`:115,613,694bytes.
- `refined/output_00002`:115,666,782bytes.
- `mixed-owner/output_00001`:64,169,380bytes.

The failed `mixed` attempt produced no raw dump. Deletions are not recoverable
from this workspace; regeneration requires rerunning the retained inputs.
No unrelated scratch, shared raw inputs or production outputs were removed.

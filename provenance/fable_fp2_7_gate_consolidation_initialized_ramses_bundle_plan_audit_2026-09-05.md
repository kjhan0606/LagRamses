# Fable F-P2.7 plan audit — gate consolidation and initialized-RAMSES

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)
Auditor: Fable
Mode: read-only plan audit; no files, tests, builds, or jobs modified/launched

## Verdict

`CONDITIONAL APPROVE`

Approve the native SNRT-enabled build gate, the initialized-RAMSES smoke, and
the checkpoint commit as mandatory. Approve the provenance index and
production-source transition in reduced form. Reject artifact-moving churn in
this bundle. The operator's instruction was to audit and then implement; the
conditions below are applied to the plan.

## Required corrections applied

- Corrected the inventory to 19 native runners, 19 smoke programs with 11
  callers and 8 orphans, six phase0 references, and eight divergent mirror
  modules.
- Corrected the premise that existing runners perform four duplicated full
  SNRT CUDA builds. They do not; the SNRT-enabled build is new required
  coverage.
- Reordered the work: checkpoint first, then consolidated gate, then
  initialized-RAMSES smoke, then phase0 runner transition.
- Reduced D1 to an index/taxonomy. Prompt consolidation and historical
  relabeling are deferred.
- Removed superseded JSON/provenance moves from D5. The index classifies them
  without path churn; large storage decisions remain separate.
- Added hard-fail requirements for missing `mpirun`, missing `nvcc`, or fewer
  than two MPI ranks; named conservation assertions, injection counts, and
  elapsed times are required in the gate summary.
- Added the initialized-smoke constraints: `&STELLAR_ENRICHMENT_PARAMS`,
  `feedback_mode='legacy'`, hydro-only level 3, future `tout/aout`,
  `nstepmax=2`, new run directory, and Slurm GPU allocation.
- D3 now ports four runners to production sources and labels phase0-only
  drivers as fixtures; it does not create another source mirror.

## Mandatory / deferred / dropped

Mandatory: checkpoint commit, one `SNRT=1 USE_CUDA=1` bundle gate, the
initialized-RAMSES smoke and diagnostic failure twin, and production-source
ownership for production-relevant stellar controls.

Deferred: prompt folding, historical relabeling, physical deletion of the
phase0 mirror after replacement tests, `.quarantine_hdf5`, the virtual
environment, and large JAX validation trees.

Dropped as redundant: moving superseded JSON/provenance files, a separate
narrative gate report, and any audit for individual runners or tiny repairs.

## Audit governance

The current bundle has one plan audit. It receives one end audit after the
single consolidated implementation/evidence pass. No intermediate audits are
planned for file moves, runners, or bounded repairs.

## Explicit non-approvals

This is engineering approval only. It does not approve physical stellar/AGN
SEDs, yields, SNIa/DTD/PISN, dust physics, live production hydro evolution,
HDF5 restart, distributed AMR scaling, or publication claims.

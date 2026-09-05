# Claude Opus 5 F-P1.2 final bundle audit — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5 (primary; read-only)
Prompt: `provenance/opus5_fp1_2_stellar_feedback_transaction_bundle_final_audit_prompt_2026-09-04.md`
Execution result: completed successfully; no files modified, jobs run, or runtime activated
Model cost reported: approximately USD 1.60

## Verdict

**PASS** — the F-P1.2 stellar/SNIa transaction bundle satisfies its own
engineering gate and is honestly scoped toward the project's production-ready
and publication-ready RT/stellar-AGN-feedback/dust objective.

## Closed items

- Generic and SNIa contributions are independent row-major scratch deltas;
  the mutating SNIa scatter adapter is absent from the mixed production prepare
  path.  The non-mutating builder signatures accept neither `unew` nor
  particle arrays, making the non-mutation boundary structural.
- Target bounds and uniqueness are enforced while legal virtual/reception
  rows remain eligible for RAMSES reverse exchange.
- The 4096-striped `omp_lock_t` is initialized before the OpenMP feedback
  region; the current row is reread under lock, progress is prepared before
  mutation, and the row/mass/progress commit has no fallible operation after
  the first shared write.
- Runtime uses `ndim==3`, energy slot `ndim+2`, all three momentum fields,
  conditional `idelay`, total metal, and active-element mappings with complete
  non-overlap validation.
- Generic and SNIa kinetic terms are computed independently, including bulk,
  cross, and source-momentum terms.  Zero returned mass with nonzero momentum
  fails closed.  The native opposed-momentum test quantitatively retains
  positive independent energy while merged net momentum is zero.
- SNIa is excluded from the generic SSP integrator, avoiding energy/mass
  double counting.  Delayed cooling receives only the SNII returned-mass
  tracer, not an invented energy reservoir.
- The direct bridge→field-map Makefile prerequisite is present and the
  recorded `-j4` object check and SNRT/CUDA dry-run evidence are consistent.
- Native evidence pins `OMP_NUM_THREADS=4`, rejects a single-thread setup,
  exercises real builder failures, checks expected full fields, and tests
  full-row same-cell accumulation.  Evidence wording distinguishes live
  production-array non-mutation from the independent transaction model.

## Deferred limitations — not blockers for F-P1.2

- MPI cross-rank atomicity for virtual/reception rows and distributed-neighbor
  transactions; RAMSES's existing reverse virtual-cell exchange remains the
  reconciliation mechanism.
- A hard process-crash exactly-once journal between hydro state and persisted
  particle/progress state.
- Live RAMSES integration execution and runtime activation.
- Physical yield/fate, SED, dust, radiation-pressure, delayed-cooling
  calibration, and publication gates.
- Dimensional generalization beyond the explicit `ndim==3` contract.

## Non-blocking recommendations

1. When the bridge helper is next touched, replace any reliance on Fortran
   logical short-circuit behavior around optional `total_metal_var` and
   `element_var` arguments with nested `if` blocks.
2. Align the native test's failure labels with its architectural zero-delta
   guarantee (the evidence document already states this accurately).
3. Make the Python commit-tail static guard whitespace-insensitive.
4. Add executable `deposit_one_star` integration coverage in the deferred
   live-runtime activation bundle.

These recommendations do not reopen F-P1.2. The hashed lock is the sole
implemented synchronization mechanism; the unimplemented named critical
fallback, process-crash journal, and MPI transaction are correctly deferred.

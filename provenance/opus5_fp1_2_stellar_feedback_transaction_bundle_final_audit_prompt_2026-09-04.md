# Claude Opus 5 final bundle audit request — F-P1.2

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only primary-auditor request. Do not edit files, run jobs, activate a
runtime, commit, or push. Two earlier Opus bundle-end attempts completed: the
first returned `CONDITIONAL PASS` for a missing direct Makefile prerequisite
and weak native evidence; after those fixes the second returned
`CONDITIONAL PASS` only for evidence integrity. The latest bounded conditions
are now addressed. This is the final primary re-audit, not a parallel audit.

Read:

1. `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`
2. `provenance/fp1_2_stellar_feedback_transaction_bundle_implementation_evidence_2026-09-04.md`
3. `provenance/opus5_fp1_2_stellar_feedback_transaction_bundle_end_audit_retry_2026-09-04.md`
4. `patch/lagRamses/stellar_ramses_runtime.f90`
5. `patch/lagRamses/stellar_ramses_bridge.f90`
6. `patch/lagRamses/stellar_ramses_field_map.f90`
7. `patch/lagRamses/stellar_snia_cell_deposition.f90`
8. `patch/lagRamses/feedback.kjhan3.f90`
9. `patch/lagRamses/stellar_progress_contract.f90`
10. `bin/Makefile`
11. `simulation/snrt/native/phase0/fp12_stellar_feedback_transaction_test.f90`
12. `simulation/snrt/tests/run_fp12_stellar_feedback_transaction.sh`
13. `simulation/snrt/tests/stellar_feedback_transaction.py`

The final objective is a production-ready and publication-ready high-level
RAMSES hydro feedback package covering RT, stellar/AGN feedback, and dust.
Assess only whether this F-P1.2 stellar/SNIa transaction bundle satisfies its
own engineering gate and is honestly scoped toward that objective.

Verify the complete production boundary: independent generic/SNIa row-major
scratch deltas; no mutating SNIa scatter during prepare; target bounds,
uniqueness, and legal virtual/reception rows; 4096-striped OpenMP locking,
under-lock reread, and no fallible post-write call; full field-map validation
including energy `ndim+2`, all momentum, SNII-only delayed tracer, metal, and
active elements; independent per-component kinetic energy and zero-mass
rejection; progress ordering; and explicit MPI/process-crash/live-runtime/
physical/dust non-goals.

Confirm the evidence-integrity corrections:

- `bin/Makefile:320` directly depends on `stellar_ramses_field_map.o` for the
  bridge, and the `-j4` object check plus SNRT/CUDA dry-run are recorded.
- The native runner exports `OMP_NUM_THREADS=4` and `OMP_DYNAMIC=FALSE`, and
  the native test checks `omp_get_max_threads()>1`.
- The native test performs actual invalid-volume, invalid-map, and `ndim=2`
  builder calls; checks zero delta; checks expected generic/SNIa fields;
  checks opposed generic/SNIa momentum with positive independent energy and
  near-zero merged energy; and checks complete-row same-cell accumulation.
- The evidence and Python output distinguish the architectural builder
  non-mutation/zero-delta test from a live production-array identity test.

Do not treat the lazy initialization hardening note or MPI/process-crash
journal as blockers for this bundle; they are explicitly deferred follow-up.
Do not require an unimplemented OpenMP critical fallback: the hashed lock is
the sole implemented synchronization mechanism and the plan says so.

Return a concise verdict within 1400 words using exactly one of:

- `PASS` — the F-P1.2 engineering bundle is accepted;
- `CONDITIONAL PASS` — list bounded mandatory follow-up conditions;
- `FAIL` — explain the blocking defect and a bounded repair.

Cite paths and approximate line numbers. End with closed items, deferred
limitations, and any required follow-up. Do not modify files or run commands.

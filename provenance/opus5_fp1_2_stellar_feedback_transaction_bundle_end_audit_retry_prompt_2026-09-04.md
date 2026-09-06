# Claude Opus 5 bundle-end re-audit request — F-P1.2 condition closure

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only primary-auditor request. Do not edit files, run jobs, activate a
runtime, commit, or push. The first bundle-end Opus call timed out after 300 s
without a verdict. The governed GPT-5.6-Sol backup then returned
`CONDITIONAL PASS` with two bounded conditions. Those conditions have now been
implemented and re-tested; this call is the primary Opus re-audit, not a
parallel audit.

Read these files:

1. `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`
2. `provenance/fp1_2_stellar_feedback_transaction_bundle_implementation_evidence_2026-09-04.md`
3. `patch/lagRamses/stellar_ramses_runtime.f90`
4. `patch/lagRamses/stellar_ramses_bridge.f90`
5. `patch/lagRamses/stellar_ramses_field_map.f90`
6. `patch/lagRamses/stellar_snia_cell_deposition.f90`
7. `patch/lagRamses/feedback.kjhan3.f90`
8. `patch/lagRamses/stellar_progress_contract.f90`
9. `bin/Makefile`
10. `simulation/snrt/native/phase0/fp12_stellar_feedback_transaction_test.f90`
11. `simulation/snrt/tests/stellar_feedback_transaction.py`
12. `simulation/snrt/tests/run_fp12_stellar_feedback_transaction.sh`

The project's final objective is a production-ready and publication-ready
high-level RAMSES hydro feedback package covering RT, stellar/AGN feedback,
and dust. Assess only whether this F-P1.2 stellar/SNIa source-to-cell
transaction bundle satisfies its stated engineering acceptance gate and is
honestly scoped toward that objective.

Recheck the earlier closed items: generic and SNIa row-major scratch deltas;
no mutating SNIa production scatter during prepare; bounds/uniqueness and
legal virtual/reception rows; 4096-striped OpenMP lock, current-row reread,
and no fallible call after the first shared mutation; full field map including
energy `ndim+2`, momentum, delayed-cooling SNII tracer, total metal, and active
elements; independent momentum/kinetic-energy accounting; progress ordering;
and explicit MPI/process-crash/live-runtime/physical non-goals.

Specifically verify the two backup conditions:

1. `stellar_ramses_bridge.o` now has a direct
   `stellar_ramses_field_map.o` Makefile prerequisite, and the dependency is
   valid for a parallel build.
2. The native test now uses actual builder failures for invalid volume,
   invalid field map, and `ndim/=3`, checks unchanged row/mass/progress,
   checks expected generic/SNIa mass, all momentum components, energy, metal,
   delayed tracer, and elements, checks opposed generic/SNIa momenta retain
   independent kinetic energy, and checks full-row same-cell accumulation
   under the striped lock.

Return a concise final verdict within 1400 words using exactly one of:

- `PASS` — the F-P1.2 engineering bundle is accepted;
- `CONDITIONAL PASS` — list bounded mandatory follow-up conditions;
- `FAIL` — explain the blocking defect and a bounded repair.

Cite file paths and approximate line numbers for material findings. End with a
short checklist separating closed items, deferred limitations, and required
follow-up. Do not modify files or run commands.

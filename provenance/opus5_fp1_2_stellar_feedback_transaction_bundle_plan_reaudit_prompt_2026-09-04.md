# Claude Opus 5 plan re-audit request — F-P1.2

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only primary-auditor request. Do not edit files, run jobs, activate a
runtime, commit, or push. The first Opus review inspected this bundle but
stopped without a final verdict, so GPT-5.6-Sol was invoked only as the
governed backup and returned `APPROVE WITH CHANGES`. The amended plan now
contains those mandatory changes. This is the primary Opus re-audit, not a
parallel audit.

Read these files first:

1. `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`
2. `provenance/gpt56sol_fp1_2_stellar_feedback_transaction_bundle_plan_backup_audit_2026-09-04.md`
3. `patch/lagRamses/stellar_ramses_runtime.f90`
4. `patch/lagRamses/stellar_ramses_bridge.f90`
5. `patch/lagRamses/stellar_snia_cell_deposition.f90`
6. `patch/lagRamses/feedback.kjhan3.f90`
7. `patch/lagRamses/stellar_progress_contract.f90`
8. `patch/lagRamses/stellar_ramses_field_map.f90`
9. `patch/lagRamses/stellar_enrichment_contract.f90`
10. `patch/lagRamses/stellar_enrichment_driver.f90`

Audit the amended plan against the project's final objective: a
production-ready and publication-ready high-level RAMSES hydro feedback
package covering RT, stellar/AGN feedback, and dust. For this bundle, assess
only whether the proposed stellar/SNIa source-to-cell transactional boundary
is physically and computationally justified, implementable in the current
Fortran/OpenMP wiring, and honestly scoped. Check especially:

- complete scratch delta for `unew(cell,variable)`, `mp_after`, and
  `indtab_after`, with no mutating SNIa call during preparation;
- valid local target and one MPI owner, named OpenMP synchronization, current
  row revalidation inside the lock, and no fallible operation after the first
  shared mutation;
- full field-map validation including density, total energy at `ndim+2`,
  momentum, delayed-cooling SNII returned-mass tracer, total metal, and active
  elements;
- independent generic/SNIa momentum and kinetic-energy accounting, including
  opposed momenta and zero-mass rejection, without merged-net-momentum loss;
- progress ordering, failure-injection evidence, same-cell concurrency tests,
  Makefile/syntax evidence, and the explicitly deferred process-crash/MPI
  transaction limitations;
- whether the bundle is the correct next high-level feedback step and does not
  accidentally claim yield-source, high-mass-fate, SED, dust, or live-runtime
  closure.

Return a concise final verdict within 1200 words using exactly one of:

- `APPROVE` — implementation may start;
- `APPROVE WITH CHANGES` — list mandatory changes before implementation;
- `REJECT` — explain the blocking defect and a bounded replacement.

If the listed evidence is sufficient, stop reading and issue the verdict; do
not spend time on unrelated historical work. Cite file paths and approximate
line numbers for material findings. End with a short acceptance checklist.

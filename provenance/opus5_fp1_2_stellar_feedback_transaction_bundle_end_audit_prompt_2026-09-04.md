# Claude Opus 5 bundle-end audit request — F-P1.2

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only primary-auditor request. Do not edit files, run jobs, activate a
runtime, commit, or push. This is the bundle-end audit after the amended
F-P1.2 plan was implemented and focused evidence was captured. GPT-5.6-Sol
was used only as the prior governed backup when the first Opus plan review did
not issue a verdict; do not invoke or emulate a parallel auditor here.

Read these files first:

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

Audit the completed bundle against the project's final objective: a
production-ready and publication-ready high-level RAMSES hydro feedback
package covering RT, stellar/AGN feedback, and dust. Assess only this bundle's
stellar/SNIa source-to-cell transaction and its evidence; defer unrelated
historical work.

Check in particular:

- whether generic and SNIa contributions are fully staged in row-major
  `unew(cell,variable)` scratch deltas, with proposed `mp_after` and
  `indtab_after`, and with no mutating SNIa production scatter during prepare;
- target resolution, index bounds, target-list uniqueness, and correct
  treatment of legal virtual/reception rows. Do not demand exactly one MPI
  owner: cross-rank atomicity is explicitly deferred and virtual rows are
  reconciled by RAMSES's existing reverse exchange;
- the named/striped OpenMP synchronization, current-row reread under lock,
  and absence of fallible calls after the first shared mutation;
- full runtime field-map validation: density, total energy at `ndim+2`, all
  three momentum components, delayed-cooling SNII returned-mass tracer, total
  metal, and active elements, including non-overlap;
- independent generic/SNIa momentum and kinetic-energy accounting, including
  opposed momenta and zero-mass rejection, with no merged-net-momentum loss;
- progress ordering, failure-injection and same-cell evidence, native compile
  and Makefile wiring evidence, and whether limitations are honestly stated;
- whether the implementation accidentally claims physical yield/fate closure,
  SED, dust, live runtime, MPI transaction, or process-crash journal closure;
- whether any mandatory correction is needed before this engineering bundle can
  be marked conditionally complete and the next high-level feedback bundle can
  be planned.

Return a concise final verdict within 1400 words using exactly one of:

- `PASS` — the bundle satisfies its stated engineering acceptance gate;
- `CONDITIONAL PASS` — list bounded mandatory follow-up conditions;
- `FAIL` — explain the blocking defect and a bounded repair.

Cite file paths and approximate line numbers for material findings. End with a
short checklist separating closed items, deferred limitations, and any required
follow-up. Do not modify files or run any command.

# Claude Opus 5 pre-implementation plan audit — F-P1.5-R

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Read-only plan audit. Do not edit files, run jobs, activate a runtime,
commit, or push. The preceding F-P1.2 stellar transaction bundle is closed
with Opus `PASS`. The next bundle has not been approved or started; audit only
the plan's validity.

Read:

1. `provenance/fp1_5_agn_effective_efficiency_convention_bundle_plan_2026-09-04.md`
2. `provenance/opus5_fp1_agn_ledger_transaction_bundle_end_audit_2026-09-04.md`
3. `patch/lagRamses/sink_particle.kjhan.f90`
4. `patch/lagRamses/snrt_ramses_driver.f90`
5. `patch/lagRamses/snrt_agn_source.f90`
6. `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
7. `simulation/snrt/snrt_core/sink_diagnostic.py`
8. `bin/Makefile`
9. `provenance/production_publication_readiness_plan.md`

Evaluate whether F-P1.5-R is the right next bounded high-level feedback task
for the final objective: a production-ready and publication-ready RAMSES
hydro package covering RT, stellar/AGN feedback, and dust. Check:

- whether the raw/effective efficiency mismatch is real and materially affects
  Lbol/photon-budget consistency;
- whether a shared pure helper can reproduce the existing writer semantics
  without silently changing the approved/legacy physics contract;
- whether the proposed inputs, invalid-value/fallback policy, MAD/eddington
  boundary, and `365.25 d`/unit semantics are complete and feasible in the
  current Fortran wiring;
- whether using the effective value in both coarse-ledger luminosity and SNRT
  photon production is the correct ownership boundary, while preserving raw
  provenance;
- whether the proposed Makefile, source-order, arithmetic, negative, and
  transaction-regression evidence is sufficient and proportionate;
- whether the plan accidentally opens runtime activation or makes a physical
  AGN SED/dust/yield claim, and whether its deferred limitations are honest;
- whether any prerequisite should be added, removed, or reordered before
  implementation.

Return a concise verdict within 1200 words using exactly one of:

- `APPROVE` — the plan is justified and implementation may begin only after
  explicit user approval;
- `APPROVE WITH CHANGES` — list mandatory bounded changes before approval;
- `REJECT` — identify the blocking scope/physics/feasibility defect and a
  bounded replacement plan.

Cite paths and approximate line numbers. End with an ordered acceptance
checklist. Do not modify files or run commands.

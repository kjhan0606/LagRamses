# Claude Opus 5 plan audit — F-P1.2 stellar feedback transaction bundle

Audit the proposed next implementation bundle in read-only mode. Do not edit
files, commit, push, launch RAMSES or other jobs, enable runtime flags, or
create approval artifacts.

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`), branch `main`.
Plan: `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`.

The final purpose is a production-ready and publication-ready high-level
hydrodynamics stack focused on radiative transfer, stellar/AGN feedback, dust,
and coupled source terms. The proposed bundle is a bounded feedback
correctness step. It must not select a physical yield source, resolve the
40–120 M☉ fate gap, activate runtime production, or claim physical closure.

Read the plan, the current production-readiness and feedback roadmaps, and the
relevant source paths, especially:

- `patch/lagRamses/stellar_ramses_runtime.f90`
- `patch/lagRamses/stellar_ramses_bridge.f90`
- `patch/lagRamses/stellar_snia_cell_deposition.f90`
- `patch/lagRamses/stellar_progress_contract.f90`
- `patch/lagRamses/stellar_enrichment_contract.f90`
- `patch/lagRamses/stellar_cell_deposition.f90`
- `patch/lagRamses/feedback.kjhan3.f90`
- `patch/lagRamses/stellar_enrichment_sources.mk` and the relevant
  `bin/Makefile` rules

Evaluate the plan, not merely whether the existing code compiles. Return
exactly one decision: `APPROVE`, `APPROVE WITH CHANGES`, or `REJECT`.

Check in particular:

1. Whether combining generic stellar and SNIa source increments in a temporary
   row-major target state is a real improvement over the current direct
   field-by-field `unew` writes, and whether all hydro fields, particle mass,
   and progress state remain unchanged on any pre-commit failure.
2. Whether the proposed energy accounting cleanly separates source/internal
   event energy from returned-mass bulk kinetic energy, source-momentum cross
   terms, and source-momentum kinetic energy, without double counting.
3. Whether the explicit SNII-only delayed-cooling reservoir convention is
   scientifically honest and sufficiently scoped, rather than silently
   presenting returned mass as a delayed-cooling energy model.
4. Whether field-map validation, unit conversion, finite/post-update checks,
   target ownership, OpenMP critical/serial semantics, and MPI boundaries are
   concrete enough for one feasible bundle. Do not demand generic AMR/HDF5
   hardening that is not a direct dependency.
5. Whether the proposed tests can demonstrate transaction rollback, mixed
   generic+SNIa closure, row-major orientation, exact-once progress ordering,
   and channel ownership without synthetic data opening production.
6. Whether the order and acceptance criteria are aligned with the final
   project purpose while preserving the current fail-closed physical-source,
   dust, live-coupling, and publication gates.

List mandatory plan changes, risks, deferred issues, and the exact acceptance
conditions. Distinguish engineering/transactional closure from physical
stellar-source or hydro-runtime approval. GPT-5.6-Sol is backup-only and must
not be invoked as an automatic parallel auditor for this plan.

# Claude Opus 5 plan audit request — F-P1.5 AGN ledger transaction bundle

Audit in read-only mode. Do not edit files, commit, push, launch RAMSES, or
enable any runtime flag. The project root is `/gpfs/kjhan/LRD_JWST` and the
repository is `kjhan0606/LagRamses` on branch `main`.

Read:

- `provenance/audit_governance_amendment_2026-09-04.md`
- `provenance/fp1_agn_ledger_transaction_bundle_plan_2026-09-04.md`
- `provenance/production_publication_readiness_plan.md`
- `provenance/feedback_implementation_plan.md`
- `patch/lagRamses/AGN_COARSE_STATE.md`
- `patch/lagRamses/sink_particle.kjhan.f90` around `AGN_feedback` and
  `dump_agn_coarse_state`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_agn_source.f90`
- `simulation/snrt/snrt_core/sink_diagnostic.py`
- `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
- `simulation/snrt/tools/reproduce_fable_sn_agn_findings.py`
- relevant AGN/SNRT tests and Makefile/source-set declarations

Assess the proposed bundle against the final project purpose: a
production/publication-ready lagRamses high-level hydrodynamics stack for
radiative transfer, stellar/AGN feedback, dust, and coupled source terms.

Answer with one verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`, and explain:

1. whether the bundle is scientifically and technically in scope while G2
   physical stellar yields remain blocked;
2. whether the proposed canonical key, duplicate policy, raw/effective
   efficiency distinction, pre-reset boundary, and algebraic checks are
   sufficient and non-circular;
3. whether an all-or-nothing group transaction actually prevents partial
   photon deposition and retry double counting in the current Fortran/JAX
   wiring, including failure paths and OpenMP/MPI implications;
4. missing acceptance tests, schema fields, or provenance/hash bindings that
   would make the plan unsafe or impossible to audit;
5. exact changes required before implementation, if any.

Keep the distinction between arithmetic/transactional evidence and physical
AGN hydro closure. Do not demand generic AMR/HDF5/checkpoint work unless it
directly affects this source transaction. Do not treat the static nine-group
AGN ledger as proof of a live physical AGN model.

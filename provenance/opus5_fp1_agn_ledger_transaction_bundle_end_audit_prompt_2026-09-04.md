# Claude Opus 5 bundle-end audit prompt — F-P1.5 AGN ledger transaction

You are the sole active auditor for this bundle.  Perform a read-only audit;
do not edit files, launch RAMSES, submit Slurm jobs, or use external web
search.  Inspect the current worktree at `/gpfs/kjhan/LRD_JWST` and the exact
files below, not a stale `/home` mirror:

- `provenance/fp1_agn_ledger_transaction_bundle_plan_2026-09-04.md`
- `patch/lagRamses/snrt_agn_source.f90`
- `patch/lagRamses/snrt_agn_source_smoke.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_agn_locator.f90`
- `patch/lagRamses/sink_particle.kjhan.f90` around
  `dump_agn_coarse_state` and `AGN_feedback`
- `bin/Makefile`
- `simulation/snrt/snrt_core/sink_diagnostic.py`
- `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
- `simulation/snrt/tools/p7_convert_sinkprops.py`
- `simulation/snrt/tools/audit_agn_coarse_ledger.py`
- `simulation/snrt/tests/agn_ledger_transaction.py`
- `simulation/snrt/data/agn_coarse_state_transaction_fixture.jsonl`
- `simulation/snrt/data/agn_coarse_ledger_transaction_audit.json`
- `patch/lagRamses/AGN_COARSE_STATE.md`
- `simulation/snrt/P4_AGN_RATE_LEDGER.md`
- `simulation/snrt/P7_NATIVE_SINKPROPS_LEDGER.md`

Project purpose: a production-ready and publication-ready lagRamses
high-level hydro stack for radiative transfer, stellar/AGN feedback, and
dust, with auditable source physics and safe coupling to the RAMSES runtime.
This bundle is deliberately narrower: it is an AGN source-ledger and
multi-group transaction prerequisite.  It does not select an AGN SED, approve
obscuration/escape fraction, close hydro feedback, promote the blocked G2
stellar yields, or claim a persistent cross-restart journal.

Assess independently whether:

1. the canonical reader validates required finite fields, pre-reset
   instantaneous semantics, 365.25-day conversion, Bondi/Eddington/Lbol
   algebra, raw versus effective efficiency, idle effective efficiency zero,
   semantic same-key duplicate collapse, and fail-closed conflicts;
2. the production Fortran writer is called exactly once before AGN deposition
   and every accumulator reset, and its emitted fields match the reader;
3. `snrt_agn_deposit_transaction` truly has a prepare/validate/commit barrier
   for all groups, with no partial intensity mutation on failure;
4. `accounted_mass` advances only after a successful transaction, is keyed by
   stable `idsink` across sink-array reorder, is serial with respect to shared
   state, and respects the local-leaf/single-MPI-owner contract;
5. `SNRT_RT_ENABLE` is latched once per process and the production Makefile
   actually links the complete SNRT module graph under the existing CUDA
   runtime gate;
6. the tests and audit JSON are sufficient arithmetic/transactional evidence
   and do not overclaim physical AGN or cross-step closure.

Pay special attention to the known limits from the plan audit: no run UUID or
dump counter (therefore rewind conflicts fail closed), live driver raw
`eps_sink` versus effective ledger mismatch remains open, deferred `Esave`
cross-step re-emission remains open, and `accounted_mass` is not a durable
crash journal.  Check that these are preserved as limitations rather than
silently declared solved.

Return a concise but evidence-backed report with:

- one verdict exactly `PASS`, `CONDITIONAL PASS`, or `FAIL`;
- closed criteria and file/line references;
- any blocking correctness or wiring defects;
- explicitly deferred/non-goal issues;
- whether the bundle is safe to close as an engineering prerequisite while
  keeping runtime activation and publication claims gated;
- a short recommended next-bundle scope if the verdict is conditional.

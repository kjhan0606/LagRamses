# Claude Opus 5 bundle-end audit request — F-P1.5-R AGN effective efficiency

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Date: 2026-09-04 (KST)

This is the final re-audit after the first Opus verdict `FAIL` and the
subsequent Opus verdict `CONDITIONAL PASS`.  The bounded repairs already made
include correcting the native MAD-high expected inflow from `1.0` to the
defined `min(0.2,1.0)=0.2`; refreshing the fixture and evidence/hash; making
the audit regex understand the block-form commit; consuming the retained
cursor only after a successful transaction; recording the supplied-ledger
upper-bound and rollover caveats; re-scoping the Python algebra banner; and
documenting the zero-Eddington policy.  The latest repair additionally carries
`efficiency_status` and `efficiency_contract_ok` through `AgnCoarseState`,
emits both in the P4 coarse CSV, and refuses non-promotable coarse rows before
writing an artifact.  A subsequent C1 repair keeps status-flagged,
non-promotable rows readable in the append-only ledger while excluding them
from `AgnCoarseState` promotion and refusing the P4 artifact.  Judge the
current files and current evidence, not any superseded audit state.

You are the sole primary bundle-end auditor.  This is a read-only audit:
do not edit files, run commands, launch jobs, set runtime flags, commit, or
push.  Do not invoke another auditor.  Use the repository files and the
recorded evidence below; if a claim cannot be verified, identify it rather
than assuming it.

The project's final purpose is a production-ready and publication-ready
high-level RAMSES hydro package covering RT, stellar/AGN feedback, and dust.
This F-P1.5-R bundle is deliberately narrower: it must close the AGN
effective-efficiency and supplied-inflow convention between the coarse-state
writer and the SNRT source driver, while preserving fail-closed physical
boundaries.  It must not be treated as AGN SED, hydro, dust, or live-runtime
approval.

Read:

1. `provenance/fp1_5_agn_effective_efficiency_convention_bundle_plan_2026-09-04.md`
2. `provenance/fp1_5_agn_effective_efficiency_convention_bundle_implementation_evidence_2026-09-04.md`
3. `provenance/opus5_fp1_5_agn_effective_efficiency_convention_plan_reaudit_2026-09-04.md`
4. `patch/lagRamses/snrt_agn_efficiency.f90`
5. `patch/lagRamses/snrt_agn_efficiency_smoke.f90`
6. `patch/lagRamses/sink_particle.kjhan.f90`
7. `patch/lagRamses/snrt_ramses_driver.f90`
8. `patch/lagRamses/snrt_agn_source.f90`
9. `patch/lagRamses/snrt_agn_source_smoke.f90`
10. `patch/lagRamses/AGN_COARSE_STATE.md`
11. `bin/Makefile`
12. `simulation/snrt/snrt_core/sink_diagnostic.py`
13. `simulation/snrt/tools/audit_agn_coarse_ledger.py`
14. `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
15. `simulation/snrt/tests/agn_effective_efficiency.py`
16. `simulation/snrt/tests/agn_ledger_transaction.py`
17. `simulation/snrt/tests/run_fp15_agn_efficiency.sh`
18. `simulation/snrt/data/agn_coarse_ledger_transaction_audit.json`
19. `simulation/snrt/data/agn_coarse_state_transaction_fixture.jsonl`
20. `provenance/opus5_fp1_5_agn_effective_efficiency_convention_bundle_conditional_audit_2026-09-05.md`

Audit these specific questions:

- Is there exactly one pure, RAMSES-independent efficiency resolver, called by
  both writer and driver, with no duplicated MAD formula?  Are spin-disabled
  default, spin-enabled uninitialized `eps_sink`, invalid raw efficiency,
  invalid rates, zero Eddington, and invalid `X_floor` explicit and
  fail-closed/non-promotable where required?  Are raw/resolved `(0,1)` and
  effective `[0,1)` conventions internally consistent, including the explicit
  spin-disabled raw-zero diagnostic exception?
- Does the writer use the helper's effective coefficient for `Lbol` and emit
  raw/resolved/effective/status/mode fields at the pre-feedback/pre-reset
  boundary?  Does the driver use the same coefficient for every photon group
  and avoid the former hidden raw-efficiency clamp?
- Does source accounting use only the increment of cumulative
  `min(dMBH_coarse,dMEd_coarse)` keyed by `idsink`, while using `dMsmbh` only
  for a one-sided retained-mass check?  Is gas-floor/Esave slack represented
  honestly?  Does the accounting marker advance only after the complete
  multi-group transaction commits, and is same-step/reordered-sink behavior
  preserved?
- Does the photon-budget API and its smoke compile unambiguously define the
  first positional argument as supplied inflow rather than retained BH mass?
  Does its direct efficiency boundary remain fail-closed for unity/super-unity
  input?
- Are the reader, converter, fixture/audit, helper SHA256 provenance, and
  direct Makefile prerequisites coherent?  Do the recorded Python/native
  tests, default full build, SNRT/CUDA driver build, dry-runs, JAX venv, and
  `git diff --check` support the engineering gate?  Distinguish actual builds
  from dry-runs and do not infer live CUDA/hydro correctness from compilation.
  Specifically verify that `AgnCoarseState` carries both efficiency status and
  contract arrays, that the P4 coarse CSV emits them, and that a false contract
  is rejected before an artifact is written.  Also verify that a status-flagged
  non-promotable row remains readable for diagnostics but cannot be promoted
  into `AgnCoarseState` or a P4 artifact.  The legacy sinkprops path must remain
  explicitly mode-unresolved rather than being mistaken for helper parity.
- Are the stated limitations honest: no runtime activation, no AGN SED or
  hydro/dust physical closure, no durable crash journal, and no physical
  publication claim from the arithmetic fixture?

Return at most 1400 words and use exactly one final verdict:

- `PASS` — the F-P1.5-R engineering bundle is accepted;
- `CONDITIONAL PASS` — list bounded mandatory follow-up conditions;
- `FAIL` — explain the blocking defect and a bounded repair.

Cite paths and approximate line numbers.  End with closed items, deferred
limitations, and required follow-up.  Do not modify files or run commands.

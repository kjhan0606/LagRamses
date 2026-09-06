# Claude Opus 5 audit prompt: independent Fable reproduction and priority plan

You are the final independent gate auditor for `/gpfs/kjhan/LRD_JWST`.
Perform a read-only audit. Do not edit files, launch RAMSES, submit/cancel
jobs, or change the implementation.

Read completely:

1. `provenance/fable_sn_agn_feedback_audit_2026-09-02.md`
2. `provenance/fable_sn_agn_independent_reproduction_2026-09-02.md`
3. `simulation/snrt/tools/reproduce_fable_sn_agn_findings.py`
4. `simulation/snrt/tests/test_fable_sn_agn_reproduction.py`
5. `provenance/feedback_implementation_plan.md`
6. `provenance/production_publication_readiness_plan.md`

Spot-check the source evidence named by the reproduction tool, especially:
`bin/Makefile`, `simulation/snrt/tests/run_g1_native_contract.sh`, the
compiled and native stellar runtime/interpolation/increment files, the HDF5
backup/restore files, the legacy Sedov expressions, the SNRT AGN driver and
sink files, the AGN diagnostic writer, and the referenced JSON/CSV metadata.

Audit questions:

- Are F1--F17 independently supported by current repository evidence rather
  than merely copied from Fable?
- Are the three `partially_reproduced` dispositions conservative and justified?
- Are any finding, line/evidence claim, numerical counterexample, or status
  materially wrong or overclaimed?
- Does the P0--P3 order correctly put production blockers before physical-table
  promotion and production reruns?
- Does the plan distinguish the native mirror PASS from the compiled runtime,
  and does it preserve the current B3 jobs as diagnostic-only work?
- What exact changes, if any, are required before this audit can be accepted?

Return a concise report with:

1. `PASS`, `CONDITIONAL`, or `BLOCK`;
2. an evidence-based summary;
3. any required corrections, grouped by severity;
4. whether the revised priority order is accepted.

Do not infer that a static source audit proves a dynamic runtime magnitude.
Treat the production/publication bar as fail-closed.

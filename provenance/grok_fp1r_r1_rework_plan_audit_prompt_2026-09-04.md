# Grok plan re-audit: F-P1R R1 evidence rework before R2

Act as the active bundle-start plan reviewer for the amended F-P1R plan in
`/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work read-only. Do not edit
files, run jobs, select physical sources, contact authors, or redistribute
data. Audit the amended plan and its governance, not implementation execution.

Current evidence HEAD is `9e422e5`; R1 implementation is `a514fd5`. Read:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/opus5_fp1r_r1_converter_fixture_audit_2026-09-04.md`
- `provenance/gpt56sol_fp1r_r1_converter_fixture_adjudication_2026-09-04.md`
- `simulation/snrt/tests/yield_converter.py`
- `simulation/snrt/tools/audit_fp1_physical_package_admission.py`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. F-P1R remains pre-admission evidence
hardening. The real state must remain fail-closed: zero physical nodes,
unresolved `[0.8,1.0]` and `[40,120] M_sun` seams, no canonical conversion or
runtime deposition, and a blocked LC18 publication gate.

Return exactly one plan decision: `APPROVE`, `APPROVE WITH CHANGES`, or
`REJECT`. Evaluate whether the amended R1 rework is a bounded, necessary,
non-destructive prerequisite for R2 and whether the plan correctly classifies
the remaining Opus findings as later hardening.

## Required R1 rework

GPT-5.6-Sol returned `REWORK R1` because the current fixture hashes only
`simulation/snrt/config` and `simulation/snrt/data`, compares before the final
post-restore check, and never directly calls the restored real physical-package
audit. The amended plan now requires:

1. staged-source per-file/composite hashes in the fixture snapshot;
2. after all synthetic seams are restored, a direct real
   `audit_physical_package_admission()` assertion for
   `blocked_no_qualified_physical_package`, `canonical_conversion_allowed is
   False`, zero physical nodes, and no selection, plus the existing real
   converter rejection and three absent output paths;
3. the final config/data/staged-source comparison only after those post-restore
   checks, covering the complete fixture window.

R1 rework must not change real contracts/data, select a source, or activate
conversion/runtime feedback. R2 must remain blocked until this evidence
rework is independently tested and Opus audits the completed rework. AGY is
retired and must not be called.

Assess whether these requirements are technically sufficient and feasible,
whether “staged source” is scoped to the repository's actual source manifest
without mutating inputs, whether the direct audit assertion is compatible with
the current fail-closed contract, and whether the order/stop rules are clear.
List mandatory amendments, risks, and valid deferrals. End by stating whether
the R1 rework may begin and whether R2 may begin after its Opus audit.

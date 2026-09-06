# GPT-5.6-Sol adjudication: F-P1R R1

Act as the independent adjudicator for Claude Opus 5's implementation-stage
audit of F-P1R step R1 in `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`).
Work read-only: do not edit files, launch jobs, select physical sources,
contact authors, or redistribute data. You may run bounded non-writing checks
if useful. Inspect the actual checkout, not just the audit summaries.

Implementation commit: `a514fd5`.
Opus audit record:
`provenance/opus5_fp1r_r1_converter_fixture_audit_2026-09-04.md`.
The accepted bundle plan and Grok amendments are:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/grok_fp1r_bundle_plan_audit_2026-09-04.md`

The final project goal is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack for radiative transfer, stellar/AGN
feedback, and dust. R1 is only evidence hardening before physical yield/fate
activation and must preserve zero physical nodes, unresolved `[0.8,1.0]` and
`[40,120] M_sun` seams, no canonical conversion/runtime deposition, and a
blocked LC18 publication gate.

Return exactly one adjudication status: `ACCEPT R1`, `REWORK R1`, or `BLOCK`.
`BLOCK` is reserved for an unsound converter admission boundary or a live
fail-closed bypass. `REWORK R1` is for a material R1 evidence gap that must be
fixed before R2. Do not convert harmless test-quality findings into a block.

## R1 implementation contract

Review `simulation/snrt/tests/yield_converter.py` and
`simulation/snrt/tools/convert_yield_rows_to_canonical.py`. The fixture must
patch only the four converter-module seams named in the plan; use temporary
contracts/output paths; derive the admitted mapping from the non-writing
proposal; reach the real positive converter write path; verify sidecar and
mapping hashes; exercise mapping-content mutation with recomputed and
unrecomputed declarations plus hash-only mutation; fail before outputs; use
`finally` restoration; preserve existing blocked tests; and prove config/data
hash invariance and post-restore fail-closed behavior.

## Opus conditions to adjudicate

1. Distinguish the equality-guard exception from the declared-hash exception
   in the mutation assertions.
2. After seam restoration, directly assert the real physical-package audit is
   `blocked_no_qualified_physical_package` with
   `canonical_conversion_allowed is False`, and include staged-source hashes.
3. Widen the invariance window to include the entire fixture and the
   post-restore blocked check.
4. Wrap the pre-existing `main()` source-contract test seam in `try/finally`.

Classify each condition as required before R2, suitable for the later R1/R2
hardening record, or unnecessary. Check for hidden false positives such as
assertions removed under `-O`, a mutation case that never reaches the intended
guard, or a temporary fixture that accidentally changes repository files.
Independently assess the driver claim that `YIELD_CONVERTER_TEST_OK` passed.
End with the exact adjudication status, a concise rationale, required fixes (if
any), and whether R2 may begin. AGY is retired and must not be called.

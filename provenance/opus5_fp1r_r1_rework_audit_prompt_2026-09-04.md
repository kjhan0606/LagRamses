# Claude Opus 5 R1 rework completion audit request

Act as the independent scientific and implementation auditor for the completed
F-P1R R1 evidence rework in `/gpfs/kjhan/LRD_JWST`
(`kjhan0606/LagRamses`). Work read-only. Do not edit files, run shell commands
or jobs, select physical sources, contact authors, or redistribute data. Use
Claude Opus 5's own judgment.

Audited rework commit: `848f328`.
Original R1 implementation: `a514fd5`.
The rework was required by GPT-5.6-Sol's `REWORK R1` adjudication and Grok's
amended plan review:

- `provenance/gpt56sol_fp1r_r1_converter_fixture_adjudication_2026-09-04.md`
- `provenance/grok_fp1r_r1_rework_plan_audit_2026-09-04.md`
- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. F-P1R is pre-admission evidence hardening. The
real state must remain fail-closed: zero physical nodes, unresolved `[0.8,1.0]`
and `[40,120] M_sun` seams, no canonical conversion/runtime deposition, and a
blocked LC18 publication gate.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` only for an unsound converter admission boundary or a live
fail-closed bypass. Use `CONDITIONAL PASS` for an R1 evidence contract gap
that must be addressed before R2. Findings explicitly deferred by the plan
are not a conditional failure. R2 may begin only on an unconditional `PASS`.

## Rework acceptance contract

Inspect `simulation/snrt/tests/yield_converter.py`,
`simulation/snrt/tools/convert_yield_rows_to_canonical.py`,
`simulation/snrt/tools/audit_fp1_physical_package_admission.py`,
`simulation/snrt/tools/audit_g2_source_package_fingerprints.py`, and
`simulation/snrt/tools/validate_fp1_source_identity_rights.py`.

The fixture must now:

1. snapshot before patching all `simulation/snrt/config` and
   `simulation/snrt/data` files plus only the 65 files listed in
   `external/g2_candidates/acquisition_manifest_v1.json`, with the existing
   code-owned LC18 per-file hashes and composite
   `3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`;
2. use the existing manifest-scoped fingerprint implementation and lock,
   without recursing over unlisted external files or writing any source,
   manifest, fingerprint, or tracked audit file;
3. restore all four synthetic converter seams in `finally`, prove their
   restored identities, and only then call the real Python
   `audit_physical_package_admission()` with its default real contract—not its
   writing CLI and not a synthetic path;
4. assert the real audit is
   `blocked_no_qualified_physical_package`, with false canonical conversion,
   runtime deposition, production, and publication flags, zero physical nodes,
   and null selection;
5. retain the real converter rejection and prove all three temporary output
   paths remain absent; and
6. compare all snapshots only after these post-restore checks, covering the
   whole synthetic fixture window. The positive path and the three mapping
   mutation cases must still be present and fail before writes.

Assess actual code and test order, hidden false positives, whether the direct
audit is genuinely the restored real function, whether manifest scope is
correct, and whether R1 remains evidence-only. The driver observed
`YIELD_CONVERTER_TEST_OK` under the normal interpreter; treat it as a claim to
evaluate, not as a substitute for source review. Classify optional findings
(distinguishing exception messages, old main seam `try/finally`, assert under
`-O`, extra sidecar fields) as deferred unless they violate the above
contract. End with the verdict, any mandatory fixes, and a direct statement on
whether R2 may begin.

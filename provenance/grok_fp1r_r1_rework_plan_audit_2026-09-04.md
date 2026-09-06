# Grok plan re-audit: F-P1R R1 evidence rework before R2

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Reviewer: Grok CLI, read-only bundle-start plan re-audit
Evidence HEAD at audit: `9e422e5`
R1 implementation: `a514fd5`
Prompt: `grok_fp1r_r1_rework_plan_audit_prompt_2026-09-04.md`

## Decision

**APPROVE WITH CHANGES**

Grok confirmed that GPT-5.6-Sol's `REWORK R1` is correct: the converter
admission boundary is sound, but the original R1 evidence contract did not
cover the complete staged-source and post-restore fail-closed state. The
proposed rework is bounded, necessary, feasible, and non-destructive. It does
not select a physical source, create nodes, resolve the fate seams, or activate
conversion/runtime feedback.

## Mandatory amendments applied

- **M8 — staged-source inventory:** hash only files listed in
  `external/g2_candidates/acquisition_manifest_v1.json`, confined to that tree,
  and include the existing code-owned LC18 per-file hashes and composite
  `3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`. Do not
  recurse over unlisted files or write manifest/source/fingerprint artifacts.
- **M9 — real audit entrypoint:** after restoring all four converter seams,
  call the Python `audit_physical_package_admission()` function with the
  default real contract. Do not invoke its writing CLI or a synthetic path.
  Assert restored binding identity, blocked status, false conversion/runtime/
  production/publication flags, zero nodes, null selection, the existing real
  converter rejection, and three absent output paths.
- **M10 — snapshot order:** take the snapshot before patching and compare
  config/data/staged-source hashes only after all post-restore audit and
  converter checks.
- **M11 — acceptance/lineage:** make the R1 acceptance gate name these checks
  and record the rework against evidence `9e422e5` and implementation
  `a514fd5`, retaining `db1bb66`/`25bd05f` as historical lineage.

## Stop rule

R1 rework may begin after these amendments are recorded. R2 remains blocked
until the rework is independently tested and Claude Opus 5 returns an
unconditional `PASS`; a conditional or negative Opus result requires GPT
adjudication and another Grok plan review. The direct admission check must use
the function, not `main()`, because the latter writes the tracked audit JSON.
The unbounded external-tree hash, per-run nonce, and synthetic contract audit
are explicitly disallowed. AGY is retired and is not part of the chain.

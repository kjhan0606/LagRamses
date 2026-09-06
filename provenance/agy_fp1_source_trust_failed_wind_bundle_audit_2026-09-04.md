# AGY bundle-end audit: F-P1 source trust and LC18 failed-wind cross-check

Date: 2026-09-04
Model: AGY `gemini-3.8-flash-high`
Prompt: `agy_fp1_source_trust_failed_wind_bundle_audit_prompt_2026-09-04.md`
Scope: uncommitted delta from base `296dd0c` in `/gpfs/kjhan/LRD_JWST`
Mode: read-only; no files/jobs/downloads/author contact

## Verdict

**PASS.** AGY found no critical or high-severity finding and judged the bundle
complete enough for the driver to design the next implementation bundle. It
confirmed that the expected blocked state is intentional: zero physical nodes,
the four physical hard blockers retained, and runtime deposition, canonical
conversion, production, and publication all disabled.

## Findings

### Low / operational: Python version dependency

AGY identified `zip(..., strict=True)` in
`simulation/snrt/tools/adapt_g2_candidate_sources.py:223` as requiring Python
>=3.10. The project runner uses the `/gpfs` `.venv` (Python 3.13.11), so the
focused runner passes; host `/usr/bin/python3` 3.9 fails. Suggested follow-up:
document the requirement or add an explicit version assertion. This is an
operational portability item, not a bundle admission bypass.

### Informational / process: dirty worktree binding

The admission audit hashes active filesystem bytes rather than Git objects or a
clean worktree. AGY verified that current unresolved intervals, missing gates,
zero physical nodes, and compiled empty production identity still prevent
promotion. It recommends committing the verified bundle coherently before the
next implementation bundle. This is retained as a process requirement, not a
claim that the current gate is open.

## Independent checks reported by AGY

AGY confirmed the code-owned lock profile, exact candidate/release/DOI/Zenodo
identity, five-file inventory, hashes and composite fingerprint
`3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`; rejection
of substitutions, coherent rewrites, symlinks, non-regular files, malformed
types/dates/paths, and runner exceptions; exact sidecar artifact pinning and
blocked-publication invariant; parser reuse; 108/108 join; 52/56 split; 19
duplicate phase rows collapsed to 845 nodes; 3--8 unique phases per model;
zero age/mass monotonicity violations; and 96 table7 records plus 12 explicit
nulls.

It also independently reproduced the source distinction: the 52 successful
rows have nonzero release Wind tables, all 56 failed rows have positive summary
wind but zero release Wind tables, CDS endpoint wind is positive for 53 failed
rows and zero/quantized for 3, and summary-to-CDS discrepancies are retained
without silent reconciliation. The four physical blockers remain unchanged.

## Disposition

The two observations are inputs to the next bundle plan. No immediate change
is required before the current bundle is recorded and committed. Opus 5 is the
second independent bundle-end auditor; comparison and reproduction occur only
after its result is available.

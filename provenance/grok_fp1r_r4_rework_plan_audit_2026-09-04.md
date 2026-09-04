# Grok plan audit: F-P1R R4 evidence rework

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Auditor: Grok Build, model `grok-4.6-build`, read-only  
Prompt: `grok_fp1r_r4_rework_plan_audit_prompt_2026-09-04.md`

## Decision

**APPROVE.** The proposed M-R4-1 rework is the minimum scientifically and
technically justified completion of the already approved F-P1R R4 boundary.
The driver may implement it now. R3 remains blocked until the rework is
tested and Claude Opus 5 returns an unconditional R4 `PASS`.

## Rationale

The defect is real but not a publication bypass: on a path mismatch no terms
bytes are opened, while the payload currently claims that a candidate record
was parsed from locked bytes. An explicit not-read/path-mismatch sentinel
fixes the overloaded `None` state at its root. The locked-path reader must
replace that sentinel with its specific success, I/O, JSON, or missing-
candidate result.

The proposed regression additions are sufficient:

- `requirements["publication_terms_record_parsed"] is False`;
- `source_terms_lock["record_source"] ==
  "not_available_due_to_terms_error"`;
- the existing denial, path, hash, and named source-rights blockers remain
  asserted.

Terms bytes remain the sole production rights authority. The rework does not
edit the terms catalog, select a physical source, create nodes, enable
conversion/runtime deposition, or broaden into RAMSES infrastructure. LC18
remains empty-approval review-only. The extra OSError branch, attestation-
mutation message branch, detached non-authoritative display labels, and
absolute `/gpfs` paths are non-blocking. AGY has no role.

After implementation, run the focused tests and deterministic G2 preflight,
then obtain a clean Opus 5 R4 implementation audit before R3. Record commit
lineage for the R4 implementation.

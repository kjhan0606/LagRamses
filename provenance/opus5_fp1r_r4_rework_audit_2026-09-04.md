# Claude Opus 5 audit: F-P1R R4 rework (M-R4-1)

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Auditor: Claude Opus 5 CLI, read-only  
Audited implementation: `ac73d52`  
Prompt: `opus5_fp1r_r4_rework_audit_prompt_2026-09-04.md`  
Prior audit: `opus5_fp1r_r4_publication_terms_audit_2026-09-04.md`  
Adjudication: `gpt56sol_fp1r_r4_publication_terms_adjudication_2026-09-04.md`  
Plan review: `grok_fp1r_r4_rework_plan_audit_2026-09-04.md`

## Verdict

**PASS.** M-R4-1 is closed with no mandatory fixes. F-P1R R3 may begin.
AGY was not called and has no active role.

## Checks established

1. `terms_error` is initialized to
   `publication_terms_not_read_path_not_code_locked` before the path check.
   A locked-path call unconditionally replaces it with the reader's success
   or specific I/O/JSON/catalog/candidate error; it cannot leak into a valid
   locked read.
2. The same `read_bytes()` buffer supplies both SHA-256 and JSON parsing.
   Malformed JSON and missing candidates retain their specific errors, while
   the valid locked catalog reports parsed=true and the locked-byte record
   source.
3. The wrong-path regression now asserts parsed=false,
   `record_source=not_available_due_to_terms_error`, and the named
   `publication_terms_record_unavailable` blocker. Mutated bytes,
   review-only, empty approval, forged dictionary, gate mutation, and
   post-evaluation terms mutation remain fail-closed; temporary fixtures do
   not touch repository terms.
4. LC18 passes no digest or detached rights record, and its live evidence
   remains `allowed=false`, `publication_ready=false`, and review-only. No
   export caller, physical node, canonical row, or runtime activation was
   added.
5. Active governance contains Opus 5, Grok, and GPT-5.6-Sol; AGY is explicitly
   retired and historical reports are provenance only.

The driver's focused tests and expected `G2_PREFLIGHT_BLOCKED` claim are
consistent with the checked-in runner/evidence. Opus did not re-execute tests
or modify files because this was a read-only audit.

## Non-blocking findings

- The LC18 test could assert the locked-path hash and record-source fields more
  directly.
- Locked-file `OSError` coverage, the attestation-specific mutation message,
  detached non-authoritative display labels, and absolute `/gpfs` evidence
  paths remain optional/deferred maintenance items.
- The bundle plan should be updated to record R4 lineage and replace its
  interim “rework required” header now that this PASS is recorded.

## R3 decision

The required R4 rework is complete and the audit is unconditional. F-P1R R3
same-run high-mass evidence freshness may begin. No AGY audit is required.

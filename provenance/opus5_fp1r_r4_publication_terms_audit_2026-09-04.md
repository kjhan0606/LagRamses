# Claude Opus 5 audit: F-P1R R4 publication-terms boundary

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Auditor: Claude Opus 5 CLI, read-only  
Implementation: `3e4a4f5`  
Prompt: `opus5_fp1r_r4_publication_terms_audit_prompt_2026-09-04.md`

## Verdict

**CONDITIONAL PASS.** The code-owned byte-authority boundary, LC18 wiring,
forged-gate refusal, terms catalog, active governance, and fail-closed
physics/admission state are correct. One material but non-bypassing evidence
accuracy defect remains, so the F-P1R stop rule requires GPT-5.6-Sol
adjudication before R3 begins.

## Confirmed controls

- The evaluator reads the terms bytes, hashes them internally, and parses
  `sources[candidate_id]` from those same bytes. No caller-supplied digest or
  detached `source_record` remains in the evaluator API or its callers.
- The private `_test_locked_terms_profile` seam is isolated to synthetic
  in-memory tests; synthetic rights fields are stored in the hashed fixture
  bytes. The test covers positive admission, review-only and approval
  failures, mutated bytes, wrong paths, malformed JSON, and a missing
  candidate record without touching repository terms.
- LC18 uses the new API with an empty approval and remains review-only with
  `allowed=false`, `publication_ready=false`, and `review_use_only=true`.
  The live terms catalog remains unedited and reports the unresolved CDS
  rights state.
- `require_publication_allowed()` rejects a plain/shallow-copied dictionary,
  detects evaluated-gate mutation, and re-evaluates current terms bytes. No
  production export caller was introduced.
- The active governance names Opus 5, Grok, and GPT-5.6-Sol and explicitly
  retires AGY; historical AGY records are provenance only.
- The fail-closed state remains unchanged: zero physical/canonical rows,
  blocked publication/runtime deposition, unresolved fate seams and failed-
  wind anomaly, and preserved 52/56, 48/4, 53/3, and 101/7 accounting.

## Mandatory fix

**M-R4-1 — false parse attestation on the unlocked-path branch.** When
`terms_path` does not match the code-owned path, the evaluator skips the file
read but currently leaves `terms_error=None`. The resulting blocked payload
therefore says `publication_terms_record_parsed=true` and
`record_source=candidate_record_parsed_from_locked_terms_bytes`, although no
file was opened and the parsed record is empty. This is not a publication
bypass because path, hash, and record-derived requirements are false, but it
is an inaccurate rights attestation at the surface R4 is intended to secure.

Use an explicit not-read/path-mismatch sentinel (cleared only after a
successful locked read), making the parsed-record requirement false and the
record source unavailable on the wrong-path branch. Extend the wrong-path
test to assert both fields.

## Non-blocking observations

1. The test does not reach the code branch that reports an already-evaluated
   gate was mutated; add a blocked-gate upgrade/mutation case.
2. The locked-file `OSError` path is not tested.
3. LC18 display-only source identity labels still come from the adapter's
   detached record, although they are correctly marked non-authoritative;
   sourcing them from the gate would improve single-sourcing.
4. Absolute `/gpfs` paths in the checked-in evidence are a portability
   consideration already deferred by the plan.

## R3 decision

R4 authority logic needs no redesign, but this is a conditional result. Run
the GPT-5.6-Sol adjudication, incorporate any confirmed evidence fix, and
only then begin F-P1R R3. AGY was not called.

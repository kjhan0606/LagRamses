# Next bundle plan: F-P1 identity and publication closure

Date drafted: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Status: Fable APPROVE WITH CHANGES; mandatory changes incorporated; implementation complete; bundle-end audit pending
Fable record: `fable_fp1_identity_publication_closure_plan_audit_2026-09-04.md`

## Authorization and parent evidence

This is the next implementation bundle after the completed F-P1 admission-
closure bundle. The user has authorized proceeding after the previous pause,
but implementation begins only after this plan is reviewed by Fable. The
parent implementation is `033799a2d2ea8618877596122f02a2007d8d64bb`; its
bundle-end audit/provenance record is commit `1db1db4`.

The parent bundle remains accepted only as evidence/admission hardening. No
physical source has been selected, no physical source nodes exist, and all
production, publication, canonical-conversion, and runtime-deposition flags
remain false.

## Purpose fit

The final purpose is a production-ready and publication-ready lagRamses
high-level hydrodynamics stack focused on radiative transfer, stellar/AGN
feedback, and dust. This bundle closes integrity boundaries that must hold
before any stellar-yield/fate package can be selected or any derived review
artifact can be published. It does not choose the unresolved 40--120 M_sun
physics, fabricate missing wind/energy/momentum data, or activate runtime
feedback.

## Bundle work packages

### B1 (P0) — bind selected package identity to executable evidence

- In the pure physical-package selection predicate, require the selected
  package SHA256 to equal the `package_fingerprint_sha256` returned by the
  selected candidate's passed `source_identity_and_rights` executable
  validator.
- Make the validator registry reject a passed identity report with a missing or
  malformed package fingerprint. Require all nine required gate reports to be
  present, identity-matched, and passed on the positive path; a list of gate
  IDs alone is not executable evidence.
- Keep the code-owned selection state and real contract review-unselected.
- Add positive and negative synthetic-registry fixtures proving that a
  self-consistent but source-unrelated package hash is rejected, while the
  synthetic registry is restored in `finally` and no synthetic approval file
  is written.

### B2 (P0) — define and bind the source-node mapping fingerprint

- Add a code-owned canonical JSON serialization helper for the exact
  `snrt-fp1-source-node-row-mapping` document. Its sorted-key, fixed-indent,
  terminal-newline bytes are the only bytes hashed. Normalize numeric types
  before serialization, use `allow_nan=False` and `ensure_ascii=True`, and
  reject or canonicalize negative zero so semantically identical mappings do
  not acquire different fingerprints.
- Extend the physical-package selection record with an explicit mapping
  document and require its schema, row count, source-node coverage, approval
  identity, source-node-contract hash, selected package hash, and canonical
  asset hash to be well-typed and internally consistent.
- Compute `source_node_mapping_sha256` from those exact mapping bytes in both
  admission and `convert_yield_rows_to_canonical.py`; require the converter's
  generated mapping to equal the admitted mapping document before writing any
  asset.
- Add a non-writing proposal mode that computes the would-be asset/mapping
  evidence without creating repository outputs. The normal admitted path must
  still embed the asset hash and refuse to write until the exact selection
  mapping matches.
- Add tests for mapping mutation, row reordering/duplicate coordinates,
  unknown node IDs, package-hash mismatch, and hash mismatch. The blocked
  checked-in selection remains `null` for both mapping fields.

### B3 (P1) — executable publication gate for derived CDS review artifacts

- Add a code-owned publication-rights predicate that requires an exact,
  hash-locked source-terms record, explicit derived-artifact publication
  approval, explicit redistribution permission for derived artifacts,
  verified production-license status, attribution, and
  `review_use_only == false`.
- Call it from the LC18 failed-wind cross-check and emit a structured
  `publication_gate` result with named blockers. Current Limongi CDS evidence
  must remain blocked because its catalogue redistribution terms are not
  explicitly identified and no derived-artifact approval exists.
- Prove that changing only the review label or report's `publication_ready`
  field cannot open the gate. Do not modify the current terms record to claim
  permission and do not publish/copy CDS-derived artifacts.
- Restrict gate inputs to the code-locked terms hash, the terms record, and the
  explicit approval record. Provide synthetic terms and lock values as test
  parameters only in memory; restore all mutable test registries in
  `finally` and assert config/data hash invariance.

### B4 (P1) — complete LC18 successful-control diagnostics

- Add symmetric successful-control counts for zero and positive CDS terminal
  wind, including the current four successful zero endpoints.
- Add the all-model zero/positive partition and assert the current total of
  seven zero endpoints (three failed and four successful).
- Keep the existing signed/relative residual definitions and all fail-closed
  anomaly/blocker semantics unchanged.

## Acceptance and stopping conditions

- F1 synthetic mismatch is rejected before any positive selection can pass;
  the real code-owned state remains unselected and the eight outstanding
  validators remain outstanding.
- The mapping hash is computed from one shared canonical serialization and the
  converter rejects any admitted-selection/mapping disagreement before file
  creation.
- The G2 manifest-scoped review fingerprint remains explicitly distinct from
  the executable source package identity; the selected package field is bound
  to the latter and never promoted by the former.
- The current LC18 report has a structured publication gate with a false
  decision and explicit rights blockers; no mutable label can override it.
- Successful-control and all-model CDS terminal-wind zero counts are
  symmetric and regression-tested: 4 successful zero, 3 failed zero, 7 total.
- Existing unresolved intervals, four admission blockers, zero physical nodes,
  all false production/publication/conversion/deposition flags, and unsent
  inquiry status remain unchanged.
- Focused tests, population/fate contract, G2 preflight with its expected
  `G2_PREFLIGHT_BLOCKED` result, JSON regeneration, fixture hash invariance,
  compilation, and `git diff --check` pass.
- No source selection, author contact, CDS redistribution, runtime feedback,
  or unrelated RAMSES infrastructure work occurs.
- After the whole bundle is implemented, obtain one AGY
  (`gemini-3.8-flash-high`) audit and one Claude Opus 5 audit, independently
  reproduce/triage both, and record the comparison before drafting the next
  plan. A FAIL result activates the previously agreed `gpt-5.6-sol`
  re-audit path; Grok remains excluded while unavailable.

## Explicitly deferred

The lower-priority F4/F5/F6/F7/F9/F10 hygiene observations, Python 3.9
compatibility, debug-compiler warnings, the unresolved 40--120 M_sun source
decision, the remaining eight physical validators, author inquiry, physical
node construction, energy/momentum realization, and runtime deposition remain
later bundles unless a new approved plan promotes them.

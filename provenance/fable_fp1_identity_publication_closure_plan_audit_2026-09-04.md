# Fable plan audit: F-P1 identity and publication closure

Date: 2026-09-04
Model: Claude Fable (`--model fable`)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Plan audited: `fp1_bundle_plan_identity_publication_closure_2026-09-04.md`

## Decision

**APPROVE WITH CHANGES.** Fable confirmed that B1--B4 target real findings in
the committed F-P1 admission closure bundle, are correctly ordered, fit in one
coherent bundle, and remain within the review-only boundary. This is not
physical-source, publication, or runtime approval.

## Mandatory changes incorporated

- The executable validator registry must reject a passed identity report with
  a missing or malformed package fingerprint. The selection predicate must
  use all nine identity-matched, passed reports, not a self-declared list of
  gate IDs; the synthetic identity runner must return the fixture fingerprint.
- The shared mapping serializer must normalize numeric types, reject NaN, use
  `allow_nan=False` and `ensure_ascii=True`, and reject or normalize negative
  zero. Admission and the converter must use the same exact bytes.
- The admitted mapping must contain the asset hash, and the converter must
  compare the exact mapping before writing. A non-writing proposal mode is
  required to compute would-be evidence without modifying the repository.
- The publication gate must be an executable, hash-locked predicate over the
  terms record and explicit approval. It must compute both `publication_ready`
  and `review_use_only`, classify the JSON and inquiry packet as internal
  review-only, and remain closed for the current Limongi CDS record. A label
  mutation must not open it; synthetic rights/lock values are test parameters,
  not mutations of module constants.
- B4 must measure and assert the 4 successful / 3 failed / 7 total zero-CDS
  endpoint partition. The existing relative-difference null denominator is a
  separate quantity and remains unchanged.
- Reconcile the G2 manifest-scoped review fingerprint wording: it is not a
  package identity claim, while the selected package field is bound to the
  executable source-identity fingerprint.

## Verified scope and deferred risks

The plan preserves the current unresolved intervals, four blockers, zero
physical nodes, false production/publication/deposition flags, and no author
contact. F4--F10 hygiene, Python 3.9, compiler warnings, exact VizieR terms
resolution, the eight remaining validators, 40--120 M_sun physics, node
construction, energy/momentum realization, runtime deposition, and decoupling
package publication from production remain deferred.

Package-level publication readiness remains coupled to production readiness in
the existing admission predicates; this bundle's B3 gate is specifically for
CDS-derived review artifacts. Generated JSON absolute paths remain a known
deferred reproducibility issue.

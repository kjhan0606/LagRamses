# Claude Opus 5 bundle-end audit — F-P1 admission closure and LC18 cross-check

Date: 2026-09-04  
Model: Claude Opus 5 (`--model opus`)  
Repository: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Audited commit: `033799a` (`Close FP1 admission coupling and cross-check gates`)  
Verdict: **CONDITIONAL PASS**

## Method limitation

This Opus session had read-only file tools but no shell tool. It inspected the
committed source, generated JSON, Fortran mirrors, tests, and provenance, but
did not execute `git show`, tests, or the claimed regeneration/M8 commands.
Those execution claims were independently reproduced by the driver and are
recorded in the bundle completion document.

## Findings

### F1 — MEDIUM, blocking before selection: package identity is not bound

`evaluate_physical_package_selection()` validates
`selected_package_sha256` and `source_node_mapping_sha256` as 64-hex strings
and checks node-to-selection agreement, but does not compare the selected
package hash to
`candidate_report[selected_id]["verified_gate_evidence"]`
`["source_identity_and_rights"]["package_fingerprint_sha256"]`. The mapping
hash has no computed definition or source binding. A future selected contract
can therefore make an internally consistent but source-unrelated package hash
pass the positive branch. This is latent today because code-owned selection is
false and eight required gates have no registered validator.

Required next-bundle action: bind the selected package hash to the executable
validator's verified composite fingerprint and define/compute the mapping
fingerprint from the exact node mapping. Update the positive-path fixture so a
fake `"a" * 64` hash cannot pass.

### F2 — MEDIUM, blocking for publication: CDS-derived review artifacts lack a
technical publication gate

`g2_source_use_terms_evidence_v1.json` marks the Limongi--Chieffi CDS catalogue
as public-with-citation for research but with no identified catalogue licence
and `production_license_status: not_approved`. The LC18 cross-check correctly
labels the material `review_use_only` and `authoritative_for_verdict: false`,
but the committed JSON embeds 108 CDS phase histories/96 structures and the
inquiry packet embeds a 56-row derived table without an executable publication
barrier. The same terms file has a more resolved VizieR precedent for a
different catalogue.

Required next-bundle action: either resolve the exact LC18 VizieR redistribution
terms with evidence, or add a technical publication/derived-artifact gate that
cannot be bypassed by the review-only label. Runtime remains safely blocked.

### F3 — LOW-to-MEDIUM, diagnostic completeness

The committed LC18 report contains 7 zero CDS terminal-wind rows: 3 failed and
4 successful control models. `failed_wind_anomaly` reports the failed 3, while
`successful_release_control` omits its symmetric zero count. Add the count and
assert it, so the control narrative does not imply all zero endpoints are in
the failed group.

### F4–F10 — low/deferred or observation

- F4: two rights checks compare code-owned constants to duplicated literals;
  the published-file path is byte-lock dependent and its negative test is
  intercepted by the source-byte lock. Keep as hygiene/coverage work.
- F5: exact path pinning occurs before the generic escape helper in the fate
  and physical admission layers; current tests exercise pinning, not those
  helper escape branches. The source-rights confinement is genuinely tested.
- F6: `selected_candidate_hard_blockers` is necessarily empty after the
  selection predicate rejects any selected blocker; this is defensive
  redundancy, not a live defect.
- F7: generated JSON contains absolute `/gpfs/kjhan/LRD_JWST` paths, so
  byte-identical regeneration is checkout-path dependent.
- F8: the audited commit hash `033799a` was not recorded in the plan,
  bundle, or roadmap at audit time.
- F9: the all-required-nodes policy is fail-closed but true domain coverage is
  deferred to the unimplemented coordinate-hull validator.
- F10: standalone F-P1H-E re-audits source-node evidence but consumes other
  upstream readiness values from locked JSON; the outer F-P1 audit re-runs the
  terminal-deposition audit.

## Confirmed sound

Opus confirmed the single coupling predicate, overclaim/stale/partial/approval
identity failures, code-owned evidence locks, nine-required-gate reachability,
LC18 counts and invariants, no CDS substitution or reconciliation tolerance,
real lstat symlink confinement, M8 isolation, and the absence of a live
production/publication/runtime bypass. It corrected the wording: there are
nine required gates, one registered executable validator, and eight outstanding
validators—not nine implemented validators.

## Disposition

The current review-only state is scientifically and operationally justified.
F1 and F2 must be carried as gating items before source selection/publication;
F3 and F5 are next-bundle work; F4/F6/F7/F9/F10 are deferred hygiene or scope
items. No production flag is opened by this audit.

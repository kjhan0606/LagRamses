# Grok audit: F-P1R bundle-start plan

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Reviewer: Grok CLI, final model reported as `grok-4.6-build`
Audit type: read-only bundle-start plan audit
Parent evidence HEAD: `db1bb66`
Implementation parent: `25bd05f`

## Decision

**APPROVE WITH CHANGES**

Grok found F-P1R to be a justified, bounded pre-admission hardening bundle
aligned with the final production-ready/publication-ready lagRamses purpose.
The bundle may proceed after the mandatory plan amendments below are recorded.
No F-P1R code or evidence was changed by this audit, and no tests or jobs were
run by the reviewer.

## Scope judgment

The plan correctly preserves the fail-closed state: zero physical source
nodes, unresolved `[0.8,1.0]` and `[40,120] M_sun` seams, no canonical
conversion or runtime deposition, and a blocked LC18 publication gate. Grok
did not require physical source selection, 40–120 M_sun fate resolution, CDS
licensing/redistribution, or runtime activation in this bundle.

R1 is safe only as a temporary/in-memory fixture around the existing converter
predicates. It must not add a production bypass or relax the blocked physical
package contracts. R2's distinction between parsed Table 5 exact-zero values
and physical-zero inference is scientifically appropriate; the 3 failed CDS
zeros must remain separate from the 56-model BR26 failed-release anomaly. R3
must cover both physical-package and fate consumers in the same invocation.
R4 must derive the candidate record from the bytes that it hashes, rather than
accepting a detached caller record or caller digest as authority.

## Mandatory amendments applied to the driver plan

- **M1 — R1 isolation:** patch only converter-module seams, use temporary
  outputs/in-memory admitted mapping and selection, restore all seams in
  `finally`, and prove the real repository path still blocks.
- **M2 — R1 mapping cases:** derive the admitted mapping from the non-writing
  proposal/generated document; test matching, recomputed-hash mutation,
  unrecomputed mutation, and hash-only mutation, with all mismatches failing
  before writes.
- **M3 — R2 semantics:** use parsed exact-zero terminology at Table 5
  precision; expose `0.01 M_sun`, the `0.005 M_sun` half-bin, and
  `physical_zero_inferred: false`; retain 52/56, 48/4, 53/3, and 101/7;
  do not name the 3 failed CDS zeros as the BR26 anomaly; leave historical
  audit reports unchanged.
- **M4 — R3 freshness:** cover both admission consumers and assert the
  same-invocation post-regeneration SHA-256 in both admission JSONs; use no
  per-run nonce; retain blocked/review-only status.
- **M5 — R4 byte authority:** read/hash locked terms internally and parse the
  candidate from those same bytes; reject caller digest/record authority and
  malformed, missing, mutated, or wrong-path inputs; retain LC18 review-only.
- **M6 — acceptance hashes:** keep staged-source bytes unchanged; require a
  second byte-identical high-mass regeneration; allow only the intended live
  LC18 semantic/name changes; do not require all 248 config/data hashes to be
  unchanged across R2.
- **M7 — lineage/order:** retain `db1bb66` as parent evidence and `25bd05f` as
  implementation parent; execute R1 first, sequence R2/R4 because both touch
  the LC18 tool, and integrate R3 last.

The amended plan is
`provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`.
Claude Opus 5 remains the implementation-stage auditor for each completed
R1–R4 step; AGY is retired and is not part of this audit chain.

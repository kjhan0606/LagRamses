# Claude Opus 5 G2 source-package staging audit — 2026-09-02

Auditor: `claude-opus-5` CLI  
Scope: read-only audit of the G2 candidate-source review and Sukhbold
lossless projection. No canonical rows, commit, push, or RAMSES integration
was created by the auditor.

## Verdict

`G2 SOURCE-PACKAGE ENGINEERING VERDICT`: **CONDITIONAL PASS**  
`G2 SOURCE-PACKAGE SCIENTIFIC VERDICT`: **BLOCK**

The scientific block is intentional and remains correct: no candidate has
project physics approval, the Sukhbold branch is solar-only and starts at
9 M☉, the 40--120 M☉ fate/yield/remnant/lifetime/energy/momentum values are
not promoted, and the review projection emits zero canonical rows and allows
zero runtime deposition.

## Confirmed

- The Sukhbold source contract pins the staged files and SHA256 values;
  candidate audit and projection policy are fail-closed.
- The complete candidate audit reports `candidate_review_only` and
  `production_ready: false`; LC18 remains quarantined for its failed-model
  wind/table inconsistency and missing age/energy/momentum fields.
- The 26 Sukhbold review records preserve 13 wind and 13 terminal components,
  stable tracked plus untracked mass, and null lifetime, release-age,
  decay-complete mass, and canonical momentum fields.
- No runtime reachability to the review projection was found, and all bounded
  candidate/projection tests passed.

## Actionable findings and disposition

1. **E1 — top-level candidate-audit fail-closed propagation.**
   `simulation/snrt/tools/audit_g2_candidate_sources.py` records nested
   manifest/inline-parser failures but does not propagate them to the
   top-level status/blockers or return a non-zero exit code. A corrupted hash
   therefore exits 0. Fix by propagating nested failures and adding a
   negative test.
2. **E2 — projection component/source-column coupling.**
   `build_g2_sukhbold_channel_projection.py` chooses the data column from the
   dict key but reports ownership from the contract field. A swapped contract
   label could relabel records without changing the source data. Fix with an
   explicit component-to-source-column identity guard and a mutation test.
3. **E3/E4/E5 — projection evidence/closure/coverage.**
   The reduced-vector closure check currently uses an `untracked` value already
   derived from the same total and is therefore not independent. The persisted
   projection lacks the source package file fingerprints, and branch counts do
   not reconcile model mass sets. Fix by retaining source fingerprints,
   adding independent source-side closure evidence, and reporting/rejecting
   branch mass-set mismatches without fabricating rows.
4. **S1 — high-mass review coverage.**
   The current projection parses only the Z9.6 yield directory; W18/N20
   high-mass tables are counted but not represented as component vectors.
   Extend the review-only parser to expose those source branches and their
   coverage/missing fields. This must not authorize interpolation or
   production deposition.

## Post-audit disposition

The actionable engineering items in this audit were implemented and rerun:

- **E1 closed:** aggregate manifest/inline parser integrity failures now set
  `candidate_review_blocked_input_integrity` and the CLI returns exit code 2;
  `g2_candidate_sources.py` mutates a manifest hash and verifies this path.
- **E2 closed:** the projection rejects a component/source-column mismatch;
  `g2_sukhbold_channel_projection.py` includes the mutation regression.
- **E3/E4 closed:** projection records retain the four source-package
  fingerprints and carry independent source mass-budget evidence without
  claiming exact closure.
- **E5/S1 closed for review coverage:** the audit parses the available W18/N20
  high-mass yield tables and the seven W18 high-mass implosion-wind vectors.
  The projection now exposes 19 high-mass review records. Missing branch
  masses remain explicit in `missing_high_mass_yield_masses`; no values are
  interpolated or invented.

The remaining documentation/package findings (candidate matrix hashes for
unstaged packages, one lower-seam convention, and LC18 rationale wording) are
not physical source approval and remain secondary review work. The scientific
verdict therefore remains **BLOCK**.

Additional documentation findings were noted: split the LC18 ranking rationale
from approval language, use one explicit lower-seam convention, and fill
candidate package hashes where available. These are secondary to E1/E2/S1 and
remain review-only.

## Next step

Implement E1, E2, and the high-mass review-only S1 parser in that order. Re-run
the candidate and projection tests, then attach the next bundled Opus audit.
The scientific gate remains blocked until source ownership, decay horizon,
age-resolved wind semantics, transition seams, canonical momentum, licensing,
and project physics approval are resolved.

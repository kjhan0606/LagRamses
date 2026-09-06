# Fable plan audit: F-P1 admission-closure bundle

Date: 2026-09-04
Model: Claude Fable (`--model fable`)
Prompt: `fable_fp1_next_bundle_plan_admission_closure_audit_prompt_2026-09-04.md`
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Decision

**APPROVE WITH CHANGES.** Fable judged the bundle aligned with the final
RT/stellar-AGN-feedback/dust high-level-hydrodynamics purpose, correctly
ordered, feasible without source selection or runtime activation, and honest
about unresolved physical blockers. It found no need to split the bundle, but
required M1--M8 design specifications before implementation.

## Fable's required changes

- **M1 P0 — code-owned admission state:** replace the contradictory metallicity
  guard with a code constant for the current review-unselected state. The
  positive path must require code-selected state, contract agreement, and all
  nine registered validators; a JSON edit alone must not open it.
- **M2 P0 — pure coupling predicate:** require fate/sidecar production and
  publication to agree with physical-package production/publication,
  conversion, deposition, node count, source-node readiness, deposition
  readiness, approval identities, and selected-candidate blockers. Test both
  stale-sidecar and stale-package directions and report the computed coupling.
- **M3 P0 — pure selection evaluation:** factor selection guards into a pure
  function. Use an in-process synthetic registry restored in `finally` to test
  each guard, the nine-way conjunction, and default-contract unreachability;
  never write synthetic approval artifacts.
- **M4 P1 — evidence lock semantics:** move all five physical evidence digests
  to code. Test code-lock/on-disk/descriptive-contract triple agreement and
  document the relock chain through the physical contract, sidecar, and audits.
- **M5 P1 — LC18 acceptance:** raise controlled errors for age, mass, and PSN
  violations; add negative fixtures; choose a differential test against the
  existing G2 phase audit; define signed/relative residual denominators as BR26
  summary wind with null for zero denominator.
- **M6 P1 — commit boundary:** record the coherent commit hash in provenance
  and roadmap; after commit, regenerated tracked audits must be byte-identical;
  external source mirrors remain untracked.
- **M7 P2 — reachability:** document/test license, published-file, and hash-read
  error paths. Published-file mismatch may co-occur with byte mismatch; registry
  hash failures must become controlled registry errors.
- **M8 P2 — fixture isolation:** hash all config/data inputs before and after
  the test matrix; preserve zero nodes, all false flags, two unresolved
  intervals, four blockers, and unsent inquiry status.

## Priority adjustment and feasibility

The order is correct through commit, coupling, evidence lock, selection-path
tests, LC18 fail-closed behavior, and adversarial tests. Fable specifically
requires the M4 lock profile before M3 positive-path tests, because otherwise
the tests exercise self-certified evidence. Fable allows the source-node rights
binding to be deferred or implemented minimally through the approved-node
fixture. It recommends a differential test rather than a G2 refactor for the
duplicated phase aggregation.

## Missing acceptance conditions added by Fable

The two unresolved fate intervals and the four blocker strings must be
byte-identical before and after the fixture matrix. Relative residuals must
declare their denominator and emit null when the denominator is zero. The
default contract must remain unreachable for positive selection, and removing
any one of the nine required passing validators must block selection.

## Deferred and out of scope

Physical source selection, the 40--120 M_sun decision, author contact, CDS
redistribution, remaining physical gate validators, coordinate-hull work,
runtime source-node deposition, momentum/energy realization, and G2 tool
refactoring remain deferred. Optional hygiene includes Python >=3.10
documentation/assertion, MD5 FIPS spelling, hex symmetry, archive reopen CRC,
and parent-directory symlink checks.

## Authorization boundary

Implementation may begin after M1--M8 are written into the next plan. This
approval is for the admission/evidence hardening bundle only; it is not
approval of a physical source or runtime feedback activation.

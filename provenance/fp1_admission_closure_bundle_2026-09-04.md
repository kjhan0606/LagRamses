# F-P1 admission-closure and cross-check hardening bundle

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Status: implementation complete; bundle-end AGY/Claude Opus 5 audit pending

## Scope

This bundle implements the Fable-approved M1--M8 plan in
[`fp1_next_bundle_plan_admission_closure_2026-09-04.md`](fp1_next_bundle_plan_admission_closure_2026-09-04.md).
It hardens the stellar-feedback source admission boundary and the review-only
LC18 failed-wind cross-check. It does not select a yield source, create
physical source nodes, enable runtime deposition, contact an author, or reopen
completed RAMSES infrastructure.

The previous trust-root bundle is fixed at commit `e000295`:
`Harden FP1 source trust and failed-wind evidence`.

## Implemented controls

- F-P1H-E now has a code-owned birth-metallicity selection state. The checked-in
  contract must agree with that state; editing JSON alone cannot select a
  package.
- Physical-package selection is evaluated by a pure, side-effect-free
  predicate. It requires all nine executable gates, non-empty unique physical
  nodes, matching package fingerprints, all upstream gates, complete approval
  flags, and the admitted status. The test-only synthetic validator registry
  is restored in `finally` and writes no project evidence.
- The five physical-package evidence digests are code-owned and are checked
  against both the contract and on-disk bytes. The audit reports code-lock and
  contract-declared digests separately.
- F-P1 fate admission now uses a single coupling predicate across the fate
  map, sidecar, physical package, source nodes, and terminal deposition. It
  rejects sidecar overclaim and stale physical-package states and reports the
  component readiness vector and shared approval identities.
- LC18 phase-history violations (duration/age, non-increasing mass, and PSN
  terminal phase) now raise a controlled error with structured diagnostics.
  The audit computes signed and BR26-summary-wind-relative residuals, derives
  one-to-one/blocker invariants from inputs, and labels CDS rights evidence and
  non-source-attested phase ordering as non-authoritative/contractual.
- Source-rights tests cover candidate-root/source-contract symlinks, duplicate
  candidate identities, license/published-file failure paths, runner and hash
  failures, and malformed inputs.

## Verification evidence

The following completed successfully:

- focused F-P1 tests for fate admission, physical-package admission,
  source identity/rights, and LC18 cross-check;
- Python compilation of changed tools and tests;
- `bash simulation/snrt/tests/run_fp1_population_fate_contract.sh`, ending in
  `FP1_POPULATION_FATE_CONTRACT_OK`;
- `bash simulation/snrt/tests/run_g2_preflight.sh`, ending in the expected
  `G2_PREFLIGHT_BLOCKED` state;
- regeneration of tracked F-P1 audit JSON, followed by a second regeneration
  with byte-identical hashes;
- an M8 matrix over 248 files under `simulation/snrt/config` and
  `simulation/snrt/data`, with path set and SHA256 manifest unchanged before
  and after the complete fixture matrix;
- `git diff --check`, JSON loading, and invariant assertions.

## Preserved fail-closed state

The current staged evidence remains review-only:

- fate unresolved intervals remain exactly
  `low_mass_lifetime_seam: [0.8, 1.0]` and
  `massive_terminal_fate_seam: [40.0, 120.0]`;
- the physical-node inventory is empty;
- physical-package and fate production/publication/canonical-conversion/
  runtime-deposition flags are all false;
- the LC18 candidate retains the four unchanged blockers: failed-model Wind
  anomaly requiring author/corrected release, missing age-resolved wind,
  missing per-node injected-energy mapping, and missing canonical
  momentum/deposition;
- no canonical rows or physical nodes are emitted and the inquiry remains
  unsent;
- the genuine staged source fingerprint remains
  `3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`.

## Audit boundary

Per-step audits are intentionally not used. This complete bundle is now ready
for independent AGY (`gemini-3.8-flash-high`) and Claude Opus 5 audits. Their
findings will be independently reproduced and triaged before the driver drafts
the next bundle; Fable will then review that next plan for final-purpose fit,
scientific/technical justification, and feasibility.

# F-P1 admission-closure and cross-check hardening bundle

Date: 2026-09-04  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Status: historical implementation complete; AGY PASS; Claude Opus 5
CONDITIONAL PASS; AGY has no continuing role; next-bundle implementation was
paused pending explicit approval
Implementation commit: `033799a2d2ea8618877596122f02a2007d8d64bb`
Audit records: `agy_fp1_admission_closure_bundle_audit_2026-09-04.md`,
`opus5_fp1_admission_closure_bundle_audit_2026-09-04.md`, and
`fp1_admission_closure_bundle_audit_comparison_2026-09-04.md`

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
  package. The contract has nine required gates; only one executable validator
  is currently registered and the other eight remain intentionally outstanding.
- Physical-package selection is evaluated by a pure, side-effect-free
  predicate. It requires all nine required gates, non-empty unique physical
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

Per-step audits were intentionally not used for this historical bundle. The
complete bundle received an independent pre-retirement AGY
(`gemini-3.8-flash-high`) **PASS** and a Claude Opus 5 **CONDITIONAL PASS**.
AGY is no longer an auditor and must not be scheduled again. The driver
independently reproduced the execution claims and triaged the findings; the
detailed records are listed at the top of this document. F1 package-hash
binding and F2 CDS-derived publication gating are blocking candidates for the
next bundle, while F3 control statistics and lower-priority hygiene items are
separately recorded.

No next implementation bundle has started. The driver will wait for explicit
user approval before beginning it; after approval, its plan remains subject to
the agreed Fable review for final-purpose fit, scientific/technical
justification, and feasibility.

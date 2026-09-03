# Next implementation bundle: F-P1 admission closure and cross-check hardening

Date drafted: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Status: Fable-approved implementation complete; AGY PASS; Claude Opus 5
CONDITIONAL PASS; next-bundle implementation paused pending explicit approval
Implementation commit: `033799a2d2ea8618877596122f02a2007d8d64bb`
Audit records: `agy_fp1_admission_closure_bundle_audit_2026-09-04.md`,
`opus5_fp1_admission_closure_bundle_audit_2026-09-04.md`, and
`fp1_admission_closure_bundle_audit_comparison_2026-09-04.md`

## Purpose fit

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. This bundle strengthens the stellar-feedback
source admission boundary so that future yield/fate data cannot be promoted by
an editable sidecar, a disconnected gate, or a successful-but-invalid review
tool. It is directly relevant to physical feedback correctness and provenance;
it does not reopen completed RAMSES topology, header, HDF5, restart, or generic
distributed-runtime work.

The bundle does not select a new stellar source, create physical nodes, send the
author inquiry, or enable runtime deposition. It is a prerequisite hardening
bundle for later physical yield implementation.

## Precondition

Before new implementation begins, finalize the currently audited source-trust /
LC18 bundle as one coherent Git commit. The commit is a provenance boundary for
the code-owned trust root; the external source mirror remains gitignored and
must continue to be verified by the committed lock profile. No production flag
is changed by this commit.

This precondition was satisfied by commit `e000295` ("Harden FP1 source trust
and failed-wind evidence").

## Part A — couple the production admission gates

### M1 (P0) — code-owned admission state

Replace the contradictory metallicity-domain guard with a code constant
declaring the current state as review-unselected. The default contract must
declare the same state. The positive path is reachable only when the code
constant declares a selected state and the contract agrees and all nine
required gates pass through approved executable validators. A JSON edit alone
must not open the path. The current tree has one registered validator and eight
outstanding validator implementations; it does not claim nine implemented
validators.

### M2 (P0) — pure coupling predicate

1. In `audit_fp1_fate_admission.py`, require the evaluated physical-package
   report to agree before F-P1 can report either `production_ready` or
   `publication_ready`: physical-package readiness, publication readiness,
   canonical conversion, runtime deposition, physical-node count, source-node
   production readiness, deposition allowed, equal approval ids, and empty
   selected-candidate blockers must be mutually consistent.
2. Add a regression fixture in which the fate map/sidecar is made ready while
   F-P1H-E remains blocked; the result must be a named controlled failure. Keep
   the real contract unresolved and all runtime consumers disabled.
3. Test the converse stale-package/sidecar direction and report the computed
   coupling booleans in the audit JSON.

The source-node approval/rights binding is optional in this bundle and may be
deferred to the approved-node fixture, as Fable permits.

## Part B — make physical-package selection testable and non-circular

### M4 (P1) — lock evidence before selection tests

1. Move the five physical-package evidence SHA256 declarations into a
   code-owned lock profile, or equivalently prove an immutable signed/committed
   trust root. Editable contract values may remain descriptive but cannot be
   the sole evidence identity.
2. Add a triple-agreement test: code lock equals on-disk bytes equals the
   descriptive contract value. Document the relock chain through the physical
   contract, sidecar digest, and regenerated data audits.

### M3 (P0) — pure selection evaluation

1. Factor selection guards into a pure function over candidate report, node
   inventory, source nodes, selection, approval, and upstream readiness flags.
   Use an in-process synthetic registry patched only for the test and restored
   in `finally`; never write synthetic approval artifacts.
2. Keep the real contract unselected and add a proof that the default contract
   is unreachable for positive selection. Removing any one of nine passing
   validators must block selection.
3. Add named tests for each selection guard: candidate qualification, non-empty
   nodes, package/mapping hashes, node fingerprint agreement, upstream gate
   conjunction, approval flags, and status.
4. Scope or rename the all-candidate `unique_hard_blockers` report so a future
   selected package cannot appear clean while an unselected candidate's
   blockers are ambiguously presented as the selected result.

## Part C — make the LC18 cross-check a true fail-closed evidence tool

### M5 (P1) — fail-closed LC18 evidence

1. Raise a controlled `Lc18FailedWindCrosscheckError` when phase age/duration,
   total-mass monotonicity, or PSN-terminal invariants fail; retain the detailed
   violation list in the error/report path where possible.
2. Add negative synthetic-record fixtures proving each new raise and print the
   violation details to stderr. Use a differential test against the existing
   G2 phase-history audit for the duplicated aggregation; do not refactor G2 in
   this bundle.
3. Compute `one_to_one` and `hard_blockers_unchanged` from the verified inputs,
   rather than emitting literal `true` values.
4. Add signed and relative residual statistics (including failed-row sign
   counts) to distinguish definitional endpoint differences from random source
   disagreement. Use BR26 summary wind as the denominator and emit null for a
   zero denominator. Do not introduce a reconciliation tolerance or substitute
   CDS values for BR26 values.
5. Mark the CDS rights fields as `authoritative_for_verdict: false` and record
   the phase-order mapping as a project contract assumption where it is not
   source-attested.

## Part D — adversarial and controlled-error coverage

### M7 (P2) — reachability and controlled errors

Extend the tests for the source-rights and admission registries to cover:

- candidate-root and source-contract root symlinks;
- duplicate candidate records in the acquisition manifest;
- a runner returning all-true requirements while reporting a blocker;
- controlled registry/validator hash-read failures;
- controlled malformed physical-package report assembly;
- license and published-file identity failures on fixtures that reach those
  checks without being intercepted earlier by a byte mismatch.

Registry hash-read failures must become controlled registry errors. The
published-file mismatch may be tested as a co-occurring byte-mismatch set when
the byte lock necessarily intercepts it first.

### M8 (P2) — fixture isolation and invariants

Hash every config/data input read by the tools before and after the whole test
matrix. Assert the existing genuine staged fingerprint, two unresolved fate
intervals, four blocker strings, zero-node inventory, all false production /
publication/conversion/deposition flags, and unsent inquiry status remain
unchanged. Do not write to the external mirror, download anything, or contact
an author.

## M6 (P1) — commit boundary

Before the next implementation commit, finalize this current audited bundle as
one coherent Git commit. Record its commit hash in this document, the bundle
provenance document, and the roadmap. After committing, regenerate the tracked
data audits and require byte-identical output. Keep the external source mirror
gitignored and untracked.

## Acceptance and stopping conditions

- current bundle is committed coherently before the next implementation commit;
- all focused tests, `run_fp1_population_fate_contract.sh`, and
  `run_g2_preflight.sh` pass, with the expected final `G2_PREFLIGHT_BLOCKED`;
- every new malformed, stale, or contradictory state fails with a controlled
  named error;
- the two unresolved fate intervals and four hard-blocker strings are
  byte-identical before and after the complete fixture matrix;
- the default contract remains unable to select a package, and removal of any
  one of the nine required passing validators blocks selection;
- relative residual statistics declare BR26 summary wind as denominator and
  emit `null` for a zero denominator;
- the real physical-node inventory remains empty and all production,
  publication, canonical-conversion, and runtime-deposition flags remain false;
- no candidate is promoted and no author inquiry is sent;
- `git diff --check`, JSON validation, and Python compilation pass;
- after the entire bundle is complete, obtain independent AGY
  (`gemini-3.8-flash-high`) and Claude Opus 5 audits. Reproduce and triage their
  findings before drafting the subsequent bundle plan.

The current plan was reviewed by Fable as **APPROVE WITH CHANGES**. Its M1--M8
changes are implemented and verified locally. AGY returned **PASS** and Claude
Opus 5 returned **CONDITIONAL PASS**; the findings and driver reproductions are
recorded in the three audit records listed above. This approval is not
physical-source or runtime feedback approval. F1 package-hash binding and F2
CDS-derived publication gating are carried forward as blocking items, with F3
control statistics and lower-priority hygiene items explicitly triaged.

No subsequent implementation bundle has started. Per the current operating
instruction, the driver is waiting for explicit user approval before starting
the next bundle.

## Out of scope

New physical yield/fate package selection, the unresolved 40--120 M_sun source
decision, author contact, CDS redistribution, runtime source-node deposition,
momentum/energy realization, and unrelated historical infrastructure remain
later bundles. They cannot be closed by these schema and admission tests.

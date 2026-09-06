# Next bundle plan: F-P1H-E executable validator and admission closure

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent bundle: F-P1R, closed at `343619e`
Status: **H5 implementation complete; final bundle-end Claude Opus 5 audit pending**

## Authorization and governance

The user explicitly approved starting the next bundle after F-P1R. The user
then replaced Grok with Claude Opus 5 as the active auditor. Opus 5 reviews
this bundle's plan and each completed work package. AGY is retired and will
not be called. A conditional or negative Opus result may use the existing
GPT-5.6-Sol adjudication path when specifically invoked before a later plan
is written.

The plan audit was attempted on 2026-09-04 with the read-only Grok CLI from
this project. Both `grok-4.6` and `grok-4.5` returned HTTP 402
`Grok Build usage balance exhausted`, including a no-tool smoke test. No Grok
verdict exists yet; this is an external service/quota blocker, not a plan
approval.

On 2026-09-04 the user explicitly instructed the driver to start
implementation despite the unavailable Grok verdict, and subsequently
replaced Grok with Opus 5. This is a direct governance record for this
bundle; the missing Grok review remains historical and no Grok approval is
claimed. The active H1 implementation audit is being performed by Opus 5.

## Final-purpose fit

The final purpose is a production-ready and publication-ready lagRamses
high-level hydrodynamics stack for radiative transfer, stellar/AGN feedback,
dust, and their coupled source terms. This bundle is a prerequisite for that
purpose: it makes the stellar physical-package gate executable and auditable
before any source values can enter the feedback runtime.

The bundle is deliberately not a physical-source approval. The current tree
must remain fail-closed with zero physical nodes and canonical rows, unresolved
`[0.8,1.0)` and `[40,120] M_sun` seams, false production/publication/
conversion/deposition flags, and no runtime feedback activation. Existing
candidate files are review inputs only; no missing energy, momentum, lifetime,
decay, remnant, or fate value may be inferred.

## Why this is the next bundle

F-P1R closed converter-path, parsed-zero, same-run freshness, and
publication-terms evidence. The F-P1H-E contract currently names nine
scientific/runtime gates but only the source-identity gate is code-registered.
The other eight are therefore reported as absent rather than as explicit
executed blockers. Before an approved source package can be admitted, this
leaves a dead production-opening surface and allows future declarative evidence
to be mistaken for executed qualification. The high-mass review also still
needs the direct cross-source/phase-history evidence to retain its fail-closed
meaning.

## Work packages

### H1 (P0) — register all nine executable gate identities

1. Extend the code-owned validator registry to cover the nine exact contract
   gates: source identity/rights; coordinate hull/population; fate structure/
   remnant; lifetime/wind history; terminal mass/species closure; decay epoch/
   projection; energy/momentum/deposition; pair instability; and runtime
   invariance/reproduction.
2. Give every validator an exact requirement set, candidate identity, code-file
   SHA256, controlled `pass`/`blocked` result, blocker list, and package
   fingerprint semantics. A validator may return `pass` only when it has
   inspected the authoritative candidate bytes and all of its requirements;
   review metadata, a non-empty status label, or a hash-only JSON artifact may
   never pass.
3. Use repository-relative, confined inputs and existing candidate adapters or
   audits where they already provide the source semantics. Do not duplicate a
   parser merely to create a gate report. Missing physical assets must produce
   named `blocked_missing_authoritative_*` findings, not synthetic zeros.
4. Make registry and contract lists exact in both directions. A missing,
   extra, mis-bound, unregistered, stale, or code-hash-mismatched validator
   must fail closed before selection.

### H2 (P0) — bind source-node rights to executed identity evidence

1. Remove the future promotion dependency on a node-local string such as
   `research_use_status` or `redistribution_status`. A physical node can be
   eligible only if its source candidate maps to the executed, code-owned
   source-identity/rights validator and that validator passes for the same
   source bytes and package fingerprint.
2. Preserve the current LC18 rights boundary and non-authoritative terms
   labeling. Do not convert review-use or scientific-use permission into
   redistribution permission, and do not change the terms catalog to claim
   rights that are not present.
3. Add an adversarial test for a node/contract that says `approved` while the
   executed rights validator is absent, blocked, or bound to another package.
   It must fail before any canonical or runtime output can be written.

### H3 (P0) — make candidate admission evidence complete and explicit

1. Require every candidate's admission record to name all nine validator IDs,
   including validators whose result is blocked. The report must distinguish
   `missing_evidence`, `validator_blocked`, `validator_error`, and `pass`.
2. Connect each validator to the existing candidate-specific source audits and
   preserve their actual blockers. The current staged candidates must remain
   blocked for their documented source gaps; no candidate may acquire a
   production qualification merely because a gate report exists.
3. Keep the pure selection predicate and its cryptographic coupling. A
   selected package still requires all nine passed reports, non-empty physical
   nodes, exact package/source-node/mapping fingerprints, upstream fate and
   deposition readiness, and matching approval identities.
4. Keep the generated admission artifacts deterministic and repository-relative
   where applicable. The normal run must continue to report zero physical nodes,
   null selection, and all approval flags false.

### H4 (P1) — close high-mass evidence tool semantics

1. Make the LC18 failed-wind cross-check exit non-zero or return a controlled
   blocked result when age, mass, or terminal-phase violations are present;
   reporting a non-zero violation count with exit 0 is not an audit pass.
2. Reuse the existing phase-history aggregation, or explicitly factor it into a
   shared implementation, so the 845-row/19-duplicate accounting cannot drift
   between two independent copies.
3. Preserve the outcome-specific parsed-zero terminology and the 52/56,
   48/4, 53/3, and 101/7 accounting. Keep the three failed CDS parsed zeros
   separate from the unresolved BR26 failed-release anomaly.
4. Keep source precision, physical-zero inference, branch wind consistency,
   and cross-source disagreement explicit. No failed-wind correction or
   terminal energy/momentum inference is permitted.

### H5 (P0) — regression and evidence window

The bundle must add or update tests for:

- all nine registry entries, exact requirement sets, code hashes, and controlled
  missing/blocked/error results;
- rights self-certification, package mismatch, source-node mismatch, and
  validator substitution;
- complete candidate evidence with current review-only blockers;
- every selection guard through an isolated temporary fixture, while proving
  the real repository remains blocked and no project artifact is written;
- LC18 violation exit behavior, shared phase aggregation, parsed-zero wording,
  and unchanged accounting totals;
- deterministic regeneration, focused F-P1 population/fate tests, the full
  `/gpfs` F-P1 runner, expected `G2_PREFLIGHT_BLOCKED`, Python/native syntax or
  compilation checks, and `git diff --check`.

The fixture window must snapshot and compare the manifest-listed staged source
bytes and tracked `simulation/snrt/config`/`data` paths after all temporary
seams are restored. It must not recurse over unlisted external files, rewrite
the external manifest, create physical nodes, or modify the candidate archives.

H5 implementation is now complete. It adds direct registry tests for
unregistered validators, gate/validator mis-binding, malformed reports, and
stale validator-code hashes. The G2 phase-history artifact is compared to the
live report in full, the LC18 cross-check artifact has the same stored-vs-live
guard, and the G2 report carries the existing caveat that intermediate-burning
phase ordering is a project-contract assumption rather than source-attested
data. Synthetic converter/yield fixtures remain test-only and cannot alter
the production admission boundary.

## Ordered implementation and audit boundary

1. H1 registry and contract wiring.
2. Claude Opus 5 read-only audit of H1; the first result was conditional and
   its two hardening conditions were repaired. The re-audit returned `PASS` in
   `provenance/opus5_fp1h_e_h1_reaudit_2026-09-04.md`.
3. H2 rights/source-node binding and H3 candidate evidence wiring.
4. Claude Opus 5 audit of H2/H3 as one implementation stage.
5. H4 high-mass cross-check semantics.
6. Claude Opus 5 audit of H4.
7. H5 registry adversarial coverage and stored-vs-live evidence freshness.
8. Full G2/F-P1 regression and bundle closure record.
9. Claude Opus 5 final bundle-end audit.

## Acceptance and stopping conditions

The bundle is accepted only if the code and contract name the same nine
validators, every current candidate has explicit executed blocked evidence,
rights cannot be self-certified, and all relevant failure paths are controlled.
The full runner and deterministic evidence checks pass while preserving the
review-only state. No physical source package is selected, no canonical row is
created, no runtime deposition is enabled, and no RAMSES evolution is
launched.

If any authoritative asset is missing, the result is a stronger named blocker
and the next physical-data requirement is recorded; it is not a reason to
invent a value or to promote a candidate. After this bundle, implementation
pauses for the next explicit user approval before another bundle begins.

## Explicit deferrals

Corrected/redistributable multi-Z and multi-rotation source acquisition,
author resolution of the LC18 anomaly, physical node population, age-resolved
wind data, fate/remnant/decay/energy/momentum/deposition approval, the
`[40,120] M_sun` physics decision, SNIa hard-crash journal and distributed
MPI/neighbour deposition, dust-mixture selection, and generic AMR/HDF5/restart
hardening remain outside this bundle.

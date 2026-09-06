# Claude Opus 5 F-P1H-E bundle-end audit

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5 CLI, read-only
Prompt: `provenance/opus5_fp1h_e_bundle_end_audit_prompt_2026-09-04.md`

## Verdict

**PASS.** Opus found no blocking defect in H1–H5. No files were edited, jobs
were run, or RAMSES build/evolution was launched.

## Findings by work package

### H1 — registry and nine-gate symmetry

- `fp1_gate_validator_blocks.py` is the registry source for nine gate IDs and
  requirement sets; `fp1_gate_validator_registry.py` registers the source
  validator plus eight explicit blocked adapters.
- The admission audit checks exact code/contract gate identity and requirement
  symmetry. The registry enforces exact report keys, identity, requirements,
  outcome/status consistency, fingerprints, and recomputed validator-code
  hashes.
- The eight unavailable adapters cannot pass, emit no package fingerprint,
  and report `authoritative_validation_available: false`.

### H2 — rights binding

- Approved nodes resolve through the code-owned locked candidate profiles,
  execute the registered source-rights validator, require verified source
  files, and require an equal package fingerprint.
- Unknown sources, blocked/substituted validators, and package mismatch fail
  closed. Empty/review contracts cannot promote.
- Node-local rights strings are not an admission trust root.

### H3 — candidate evidence and selection

- Every candidate names all nine gate validators and receives an explicit
  `missing_evidence`, `validator_blocked`, `validator_error`, or `pass`
  status.
- Selection requires all nine verified pass reports, complete physical nodes,
  exact package/source-node/mapping fingerprints, and upstream readiness.
- Current state remains four candidates, one LC18 source-rights pass plus
  eight blockers, three candidates with nine blockers, zero physical nodes,
  and null selection.

### H4 — phase-history semantics

- Both consumers use `fp1_limongi_phase_history.py`, including its fail-closed
  source-order, duration/age, total-mass, cumulative-wind, PSN, duplicate, and
  phase-label invariants.
- Accounting remains 108 models, 845 unique phase rows, 19 exact collapsed
  duplicates, and 52/56, 48/4, 53/3, 101/7 outcome/parsed-zero counts.
- Parsed zero is not interpreted as physical zero. Failed-wind anomalies,
  cross-source residuals, source precision, and no-inference semantics remain
  explicit; no wind, energy, momentum, or reconciliation value was invented.

### H5 — adversarial and freshness evidence

- Direct registry tests cover unregistered, gate/validator mis-bound,
  malformed-report, and stale validator-code-hash cases.
- G2 and LC18 checked-in reports are compared to live reports in full.
- The G2 report carries the caveat that intermediate-burning phase ordering is
  a project-contract assumption, not source-attested data.
- Synthetic source fixtures are confined to tests; the production H2 boundary
  remains intact.

## Driver evidence accepted by Opus

The driver reported successful focused tests, `compileall`, `git diff --check`,
the full F-P1 runner, and the full G2 preflight. The G2 preflight terminal
state was `G2_PREFLIGHT_BLOCKED`, which is the expected admission result while
no authoritative physical source package is qualified.

## Non-blocking findings

1. The audit module repeats the nine gate names as a literal set instead of
   deriving them directly from `GATE_REQUIREMENTS`; traced drift directions
   still fail closed.
2. The test-only converter monkeypatch is not wrapped in `try/finally`.
3. The pure selection helper trusts provenance of reports supplied by its
   production caller; that caller already validates registry provenance.
4. This closure record was not present at the time of the audit; it is now
   supplied as `fp1h_e_validator_admission_bundle_closure_2026-09-04.md`.

These are low/informational maintainability items, not conditions for bundle
closure and not unsafe physical promotion paths.

## Limitations

Opus did not rerun the Fortran/G2 or F-P1 runners in read-only mode. The driver
executed them independently. Opus also did not rehash every external archive;
the repository's lock-pinned source identities and generated evidence were
checked for consistency.

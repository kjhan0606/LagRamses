# F-P1H-E executable validator and admission closure

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent bundle: F-P1R, closed at `343619e`
Bundle plan: `provenance/fp1h_e_validator_admission_bundle_plan_2026-09-04.md`
Final auditor: Claude Opus 5 CLI, read-only

## Final verdict

**PASS.** F-P1H-E is closed as an executable, fail-closed admission and
evidence bundle. Claude Opus 5 replaced Grok as the active auditor and found
no blocking defect in the H1–H5 implementation.

This is not physical-source approval. The repository remains deliberately
blocked: four review candidates, zero physical source nodes and canonical
rows, null package selection, and false production/publication/
canonical-conversion/runtime-deposition flags. The `[0.8,1.0)` lifetime seam
and `[40,120] M_sun` high-mass fate seam remain unresolved. No missing
physical values were inferred.

## Closed work packages

- H1 registers exactly the nine F-P1 gates, with exact requirements and
  strict code-hash/report validation. The eight unavailable physical gates
  are explicit never-passing adapters.
- H2 binds an approved node to the code-owned source-rights validator,
  verified source bytes, and an equal package fingerprint. Node-local rights
  declarations cannot promote data.
- H3 requires complete nine-gate evidence for every candidate and preserves
  the current LC18 one-pass/eight-blocked and three-candidate nine-blocked
  state.
- H4 factors the Limongi phase-history aggregation into one shared,
  source-order fail-closed implementation. It preserves 108 models, 845
  unique phase rows, 19 exact collapsed duplicates, and the 52/56, 48/4,
  53/3, and 101/7 accounting without wind, energy, momentum, or
  cross-source inference.
- H5 adds direct registry adversarial coverage for unregistered,
  mis-bound, malformed-report, and stale-code-hash validators; full
  stored-vs-live checks for the G2 and LC18 reports; and an explicit
  intermediate-burning-order provenance caveat in the G2 report.

## Verification evidence

On GPFS, the driver passed the focused H1–H5 tests, Python `compileall`, and
`git diff --check`. It also passed:

```text
FP1_POPULATION_FATE_CONTRACT_OK
G2_PREFLIGHT_BLOCKED
```

The latter is the expected successful fail-closed terminal state: the G2
preflight reaches all validation stages and refuses to open production while
no authoritative physical package is qualified.

The final Opus audit is recorded at:

`provenance/opus5_fp1h_e_bundle_end_audit_2026-09-04.md`, with prompt at
`provenance/opus5_fp1h_e_bundle_end_audit_prompt_2026-09-04.md`.

## Non-blocking follow-up observations

Opus noted four maintainability items that do not invalidate this closure:

1. derive the audit module's literal `REQUIRED_GATES` set from the registry
   source of truth;
2. make the test-only converter monkeypatch exception-safe with `try/finally`;
3. optionally make the pure selection helper independently validate the
   code-registered provenance of supplied gate reports;
4. retain the final Opus audit record itself as this closure artifact.

Items 1–3 are deferred to a future explicitly approved bundle because they
cannot enable promotion through the current production caller and changing
them would require another complete audit cycle. Item 4 is closed by this
record.

No commit, push, RAMSES build, or evolution run was performed as part of this
closure.

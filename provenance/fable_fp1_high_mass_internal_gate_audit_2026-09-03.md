# Fable audit — F-P1 40–120 M☉ internal gates

Date: 2026-09-03
Auditor: actual `fable` through `claude -p --model fable`
Mode: read-only inspection of `/gpfs/kjhan/LRD_JWST`; regenerated artifacts
were written only under `/tmp`
Prompt: `fp1_high_mass_internal_gate_dual_audit_prompt_2026-09-03.md`

## Verdict

- Top level: **CONDITIONAL PASS for engineering / BLOCK for physics**.
- F-P1H-A--E engineering controls: **CONDITIONAL PASS**.
- Physical 40--120 M☉ gap: **BLOCK**.
- Production readiness: **BLOCK**.
- Publication readiness: **BLOCK**.

All eight F-P1 Python suites, candidate-grid coverage, the ifx policy unit test,
and `git diff --check` passed. Eight regenerated audit artifacts were
byte-identical and all checked sidecar hashes matched. Source parity remained
blocked on `production_linked_build_evidence`.

## Findings

1. **High:** no Fortran channel-3 fate consumer exists. The widened window is
   integrated without reading a source-node outcome; only the top-level gate
   prevents deposition. A linked-binary negative test and a 121 M☉ coverage
   negative are needed.
2. **High:** the compiled identity is not in the recorded production binary;
   F-P1H-A currently has source-regex and standalone-unit evidence only.
3. **Medium:** converter/asset “binding” accepts arbitrary 64-hex node-contract
   and mapping hashes instead of hashing the actual files. Energy/momentum
   metadata are only checked as non-empty.
4. **Medium:** package `passed_gate_ids` is self-declared and no future-capable
   per-gate evidence evaluator exists. No physical-node validator enforces the
   84-field null/zero/direct-collapse contract.
5. **Medium:** failed-remnant count, tolerance result, radioactive admissibility,
   and package erratum requirement contain literal conclusions; corresponding
   tests are partly tautological.
6. **Medium:** the F-P1 runner is not called by an aggregate runner.
7. **Medium:** the F-P1 provenance says current source parity passes although it
   is now blocked.
8. **Low:** `currently_resolved_terminal_intervals_msun` overstates an 8--40 M☉
   interval that is only a configured, not physics-approved, window.
9. **Low:** the future runtime fate-policy input/promotion route is not wired.
10. **Low:** generated-audit hash pinning is reproducible but fragile.

Fable independently reproduced the same headline counts as AGY. It additionally
found that all five failed N20 nodes lack wind records and that seven isotopes
outside the current 6-dex warning threshold differ by 1--3 dex at 80/100 M☉.

## Required remediation

Before F-P1H-F, Fable requires computed evidence rather than literals, actual
contract/mapping hash verification, a source-node validator, a per-gate evidence
evaluator, aggregate-runner wiring, and wording corrections. Before runtime
activation it also requires the channel-3 fate consumer, coverage negative,
clean linked rebuild/parity, and a defined fate-policy input route.

## Candor assessment

The main documents correctly say the physical gap is unresolved. Fable flagged
only the phrase implying an implemented runtime fate filter and the stale
source-parity sentence as overclaims.

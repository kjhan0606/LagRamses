# Audit governance amendment

Effective: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Current auditor roles

- **Claude Opus 5:** implementation-stage auditor. Each completed
  implementation step receives an independent algorithm, wiring, physical
  legitimacy, and evidence audit before the next step is accepted.
- **Grok:** bundle-start plan auditor. Before a new implementation bundle
  begins, Grok reviews the driver's plan for final-purpose alignment,
  scientific/technical justification, feasibility, scope discipline, and
  acceptance criteria. The driver waits for that plan decision.
- **AGY:** retired from the active auditor roster effective on this date. No
  new AGY audit is commissioned or used as an approval gate.

Historical AGY reports remain immutable provenance for work performed before
retirement. In particular, the F-P1 identity/publication closure AGY `PASS`
is a historical bundle-end result, not a future approval authority.

This amendment supersedes any earlier plan wording that says an AGY audit is
pending, should be commissioned, or is required to open a gate. Such wording
is historical workflow text only; it is not an executable instruction. The
historical files and their findings are retained for provenance and are not
rewritten or deleted.

If a future Opus audit is negative or conditional under the existing project
rule, GPT-5.6 Sol may perform the adjudication re-audit. Any resulting
remediation becomes part of the driver's next plan, which Grok must review
before implementation. Grok is currently the replacement for the previously
assigned Fable plan-review role.

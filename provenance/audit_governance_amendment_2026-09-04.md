# Audit governance amendment

Effective: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Current auditor roles

- **Fable:** primary plan auditor. Before a new implementation bundle begins,
  Fable reviews the driver's plan for final-purpose alignment,
  scientific/technical justification, feasibility, scope discipline, and
  acceptance criteria. Fable is also the end-of-bundle backup if Opus 5
  cannot perform or decide the implementation audit.
- **Claude Opus 5:** primary end-of-bundle auditor and backup plan auditor.
  At the end of each completed implementation bundle, Opus 5 performs one
  bundled read-only audit of the algorithm, wiring, physical legitimacy, and
  evidence. If Fable cannot perform or decide the plan audit, Opus 5 reviews
  the plan instead.
- **Grok:** removed from the active auditor roster by the latest user
  directive. The exhausted-quota attempt remains historical provenance only;
  no Grok approval is required or claimed for future bundles.
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

GPT-5.6 Sol is not an automatic independent auditor. It is a backup auditor
and may be called only in either of these cases: (1) Opus 5 does not perform
the audit or cannot issue a verdict, or (2) the operator does not regard the
Opus result as sufficiently trustworthy and explicitly requests a separate
check. It is not invoked merely because a bundle is under review, and it is
not run in parallel by default. Any remediation resulting from an invoked
backup audit becomes part of the driver's next plan, which Opus 5 reviews
before implementation.

## Latest directive

On 2026-09-05 the user explicitly exchanged the primary and backup auditors at
bundle start: Fable is now primary for plan authorization and Claude Opus 5 is
the backup. At bundle end, Claude Opus 5 remains primary and Fable is the
backup. Grok, AGY, and automatic parallel audits remain disabled.

# Project instructions: lagRamses high-level RT / feedback / dust

## Completed test output retention

Operator directive (2026-09-10): after a simulation test AND its stage
evaluation finish, remove its raw simulation outputs. Retain the effective
inputs, logs, build identity, compact measured results and a cleanup manifest.
Keep checkpoints still needed for an active restart or CPU/GPU comparison
until that evaluation ends. Resolve exact file targets before deleting;
do not sweep unrelated workers' runs or production/scientific archives.

This workspace tracks `kjhan0606/LagRamses`. Confirm the actual working
directory and origin before mutations; do not infer project identity from
the directory name alone.

## Planning audits

Operator-approved reduction (2026-09-10): Fable is not a routine per-bundle
gate. Request a plan review for a new physical model, a substantive change
to conservation or module coupling, or a scientific uncertainty the driver
cannot resolve confidently. One review may cover multiple implementation
bundles within the same design. Continuing approved work, extending data
coverage within that design, bounded fixes, regression tests and documentation
do not require another call. The driver performs end evaluations; optional
review suggestions are not new mandatory completion conditions. The already
reviewed low-Z PARSEC plan proceeds through implementation and verification
without another review unless a substantive design issue meets these triggers.
Existing preapprovals remain valid.

Per operator directive (2026-09-09), every future plan-audit prompt must ask
**Q-GOAL** (alignment with the final simulation-ready RT/feedback/dust goal)
FIRST, then **Q-LEAN** (excessive instrumentation or gates), and request answers
in that order before detailed technical findings. Include the approved scope
and existing completion criteria. See
[audit cadence](provenance/audit_cadence_amendment_2026-09-05.md).

These are two questions within the existing review, not new gates or audits.
Do not reopen an approved bundle for this instruction. Preserve the latest
operator-defined auditor roles, cadence and preapprovals.

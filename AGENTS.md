# Project instructions: lagRamses high-level RT / feedback / dust

## General physical-model completion scope

Operator decision (2026-09-11): completion means a galaxy-formation effective
model using justified, declared literature-based approximations, with
consistent mass/energy/element accounting and working RT/feedback/dust
coupling inside its stated domain. Direct microscopic calculations are
long-term research, NOT prerequisites or new gates for this closeout.
See `provenance/general_physics_completion_plan_2026-09-11.md` and
`provenance/microscopic_physics_long_term_backlog.md`.
Do not relabel comparison limits as calibrated physics, remove physical
admission checks, or claim implementation complete merely from this scope
approval. Existing scientific limitations remain documented.

COLIBRE-type galaxy-observable calibration and resolution convergence form
a separate science plan (`provenance/galaxy_calibration_science_plan.md`),
not additional prerequisites for the approved three-step implementation.

## Short-term M5 transport and coupling objective

Operator decision (2026-09-14): the short-term objective is to select and
implement conservative transport of M5 moments, then connect it to the
existing gas and dust physics. M5 is the angular representation/closure;
spatial transport and radiation-matter coupling are separate design choices.
Retain the current conservative kinetic M5 transport and S_N path as comparison
options. A direct moment Riemann solver is a candidate, not an established cure
for the observed angular/spatial artifacts.

Judge transport by propagation and flux errors, shadows, conservation,
positive/realizable states, and runtime/memory cost, not spherical appearance
alone. Preserve physical anisotropy in inhomogeneous density fields. Reuse the
existing controls documented in
`provenance/snrt_m5_origin_controls_2026-09-14.md`; do not create an expanding
benchmark or gate program as a substitute for implementation.

Reuse the existing gas and dust source-term interfaces where applicable,
checking photon, energy and momentum exchange consistently with the selected
reduced-light-speed convention and enabled physics. Standalone fixed-gas M5
tests do not establish coupled simulation readiness: a bounded integrated
gas/dust run must verify the connection. Unrelated physics extensions and
galaxy calibration remain outside this short-term objective.

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

Current operator-role update (2026-09-16): Fable performs plan reviews and
bundle-end audits; Claude Opus 5 performs code implementation/structure
audits. The default Fable plan-audit timeout is 600 seconds (10 minutes). A
timeout is recorded as no verdict, not as approval or rejection.

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

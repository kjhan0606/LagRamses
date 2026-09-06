# Audit cadence amendment

Effective: 2026-09-05
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Major-bundle rule

Implementation work is organized into larger coherent bundles.  A bundle
should contain the related algorithm, production wiring, boundary semantics,
tests, evidence, and documentation needed to make one meaningful engineering
or physics decision.  A single helper edit, test correction, wording change,
or bounded audit repair is not a new bundle.

The plan auditor is called once before each substantial bundle and the
implementation auditor once at its end, not after each micro-step or
intermediate test. Fable is primary for the plan review and Claude Opus 5 is
its backup; Claude Opus 5 is primary for the end review and Fable is its
backup. During implementation, local tests and static checks are ordinary
engineering evidence, not separate audit events.

If a bundle-end audit returns bounded conditions, the driver records them and
folds all related remediation into one larger repair/closure bundle.  A single
re-audit is requested only after that repair bundle's implementation, focused
evidence, and documentation are complete.  GPT-5.6-Sol remains a reserve only
when the active auditor chain cannot issue a verdict or the operator
explicitly requests an additional check.

An earlier audit is justified only for a user-requested review or an urgent
safety/physics boundary that must be resolved before continuing.  Otherwise,
the next bundle starts only after the operator's approval, and its scope is
set at the level of a complete high-level RT/feedback/dust engineering task.

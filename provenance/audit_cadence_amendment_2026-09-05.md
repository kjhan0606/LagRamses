# Audit cadence amendment

Effective: 2026-09-05  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

## Major-bundle rule

Implementation work is organized into larger coherent bundles.  A bundle
should contain the related algorithm, production wiring, boundary semantics,
tests, evidence, and documentation needed to make one meaningful engineering
or physics decision.  A single helper edit, test correction, wording change,
or bounded audit repair is not a new bundle.

Claude Opus 5 remains the sole active auditor, but it is called once at the
end of each substantial bundle, not after each micro-step or intermediate
test.  During implementation, local tests and static checks are ordinary
engineering evidence, not separate audit events.

If a bundle-end audit returns bounded conditions, the driver records them and
folds all related remediation into one larger repair/closure bundle.  A single
re-audit is requested only after that repair bundle's implementation, focused
evidence, and documentation are complete.  GPT-5.6-Sol remains a backup only
when Opus cannot issue a verdict or the operator explicitly requests it.

An earlier audit is justified only for a user-requested review or an urgent
safety/physics boundary that must be resolved before continuing.  Otherwise,
the next bundle starts only after the operator's approval, and its scope is
set at the level of a complete high-level RT/feedback/dust engineering task.

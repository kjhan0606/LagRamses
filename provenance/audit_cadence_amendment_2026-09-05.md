# Audit cadence amendment

Effective: 2026-09-05
Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001` (`kjhan0606/LagRamses`)

## Current policy: selective Astra review (operator update 2026-11-05)

This section supersedes the historical per-bundle review/re-audit cadence,
external end-auditor assignment and approval-wait wording below. Existing
operator preapprovals remain effective; do not add an approval wait.

The 2026-11-05 operator update replaces Fable with the internal Astra Codex
sub-agent (`gpt-6-astra`) for every pending or future Fable audit assignment,
including plan/design review, scientific or technical consultation,
bundle-end/closeout audit, and fallback review. This changes the auditor only;
it does not increase review cadence or reopen prior approvals. Claude Opus 5
retains its separately assigned code implementation/structure audits. Use the
internal Astra agent, not an external CLI. Give a plan review up to 600 seconds
(10 minutes); a timeout is recorded as no verdict, never as approval or
rejection. Completed Fable reports remain historical records; do not rewrite
or relabel their verdicts.

- Call Astra for a new physical model, a substantive conservation-law or
  module-coupling design change, or a scientific question the driver cannot
  confidently resolve. Explicit operator requests also take precedence.
- Do not call routinely for continued implementation of an approved design,
  data-coverage extensions within that design, bounded bug fixes, regression
  tests, documentation or cleanup. A substantive scientific change discovered
  during such work is evaluated against the triggers above, not its task label.
- One plan review can cover several bundles within the same design. Bounded
  corrections are checked by the driver without an automatic re-audit.
- When a review is needed, ask Q-GOAL first and Q-LEAN second. Do not promote
  optional recommendations to mandatory completion conditions or extra gates.
- The driver performs implementation/end evaluation. No routine external
  end review or automatic backup audit is required.
- The low-Z PARSEC precision plan has already been reviewed. Implement and
  verify that scope without a further call unless a substantive new design
  issue meets the triggers above.

Keep proportionate physical-correctness and regression checks. This reduces
external review frequency, not conservation requirements or test-output
retention rules.

## First questions in every plan audit (operator directive 2026-09-09)

Every future plan-audit request must ask these questions FIRST, in this order,
and request answers in the same order before detailed technical findings:

1. **Q-GOAL:** Does this bundle directly advance the final objective: physically
   justified, simulation-ready RT, stellar/AGN feedback and dust in lagRamses?
   Identify the actual runtime capability it delivers and any work unrelated
   to the approved objective or necessary completion criteria.
2. **Q-LEAN:** Is the plan overinstrumented or over-gated? Identify unnecessary
   new instrumentation, duplicate tests/reviews, excessive gates and artificial
   subdivision. Recommend removing, merging or deferring them while retaining
   the minimum evidence needed for physical correctness and safe execution.

Supply the final objective, approved scope and existing completion criteria
as context. For each question require a short conclusion and concrete reasons;
then proceed to the existing scientific/technical and feasibility review.
Do not turn optional improvements into new mandatory completion conditions.
These are priority questions WITHIN the existing single plan audit, not two
new gates, audit calls or approval waits. The Astra role replacement above
does not create extra audits; do not reopen an already approved bundle solely
for this update. Any still-pending request that names Fable is reassigned to
Astra under the same existing trigger, scope, and approval state.

## Historical major-bundle rule (cadence superseded above; auditor assignment superseded 2026-11-05)

Implementation work is organized into larger coherent bundles.  A bundle
should contain the related algorithm, production wiring, boundary semantics,
tests, evidence, and documentation needed to make one meaningful engineering
or physics decision.  A single helper edit, test correction, wording change,
or bounded audit repair is not a new bundle.

The old assignment called the plan auditor once before each substantial
bundle and the implementation auditor once at its end, not after each
micro-step or intermediate test. Fable/Claude Opus 5 were the historical
plan-review pair, and Claude Opus 5/Fable were the historical end-review
pair. Effective 2026-11-05, Astra replaces Fable at those assigned call points
under the current policy above. During implementation, local tests and
static checks are ordinary engineering evidence, not separate audit events.

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

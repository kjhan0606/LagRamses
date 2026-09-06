# Opus 5 DUST-5 end audit

Actual model claude-opus-5; session 30893bf5-d0d8-4e52-8265-6d8402b5d391;
duration 381 s. Read-only review, no edits/jobs/subagents/web. This is a
summary of the report with separate driver dispositions.

## Verdict: CONDITIONAL PASS

Two documentation-only conditions; no code repair requested.

Verified: Planck units/Bose series, total-power-by-node-sum construction,
bath-relative source, log-T-in-power inversion/range guard, upwind face signs
and reciprocal-flux conservation, old-field escape, thin-cell source response,
normalized angular weights, separate native FP64 energy state, photon
conversion at fixed frequency, and all admission/error paths. Native thick-
branch `1-exp(-tau)` can differ from expm1 near the branch threshold; this
does not violate transmission+loss conservation and is non-blocking here.

Only the success branch writes caller-owned state/diagnostics; local trials
always start at the same old state. Three bitwise rollback checks exercise
late nonconvergence, CFL and graph failures. Zero/weak-source and two-compiler
physical-opacity differential evidence are adequate for the bounded plan.
Build graph additions are minimal and consistent; module is dependency-free
and not invoked by the live driver. Evidence correctly distinguishes native
operator compilation from full-link/live/science approval.

## Conditions, both applied

1. State that exact energy closure is at the nonlinear fixed point; committed
   steps close only to the finite stop tolerance (1e-9, measured 7.407e-10).
   Added to the native guide, including the zero-primary old-inventory scale.
   Future gas/force coupling must budget this residual explicitly.
2. Intel defines are flag-path exercise only, because this source contains
   no preprocessor conditionals. Added to the evidence; no extra physics or
   production branch coverage is inferred.

Driver precision notes: agreement of two spectral moments is useful evidence,
not a mathematical proof that arbitrary spectral permutations are excluded.
The direct algebra/code review also matters. The zero-source smoke has zero
old inventory; it is not a source-off/nonzero-stock decay experiment. Neither
point requires expanding this bounded audit into a new test framework.

No numerical rerun is needed for these documentation-only changes: source and
all successfully tested numerical paths are unchanged. No further audit round.

Deferred: dust abundance/depletion/mixture science, persistent live state and
spectral layout/resolution, primary/gas/force transaction and RSLA derivation,
AMR/MPI/restart, and production-scale cost of the native source interpolation.
Reciprocity guarantees conservative graph flux, not Cartesian geometric
validity; an actual live mapper must supply correct geometry and boundaries.

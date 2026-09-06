# Codex end audit: G4 thermal-atlas domain guard bundle

Date: 2026-09-06
Project: `/gpfs/kjhan/LRD_JWST`
Base: `be70375`

## Verdict

**PASS within the bounded engineering scope; conditional for
production/publication.**

P5 now rejects out-of-domain epoch and initial gas-state inputs before JAX
execution/output creation, while retaining the low-level clamped interpolator
only for explicitly controlled offline callers. Runtime thermal states that
leave the admitted temperature interval are already represented as bound hits
and fail P5 validation. The full G4 mapping/source/thermal closure gate passes.

## Scope check

The narrow thermal atlas was not extrapolated or promoted to broad redshift
coverage. No physical dust mixture, DTM/depletion prescription, source-cell
model, live dust state, gas/force feedback, or cosmological qualification is
claimed.

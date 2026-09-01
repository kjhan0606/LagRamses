# B2 audit — independent Claude Opus 5 review

Date: 2026-09-01  
Root: `/gpfs/kjhan/LRD_JWST`  
Mode: read-only; no calculations executed by the auditor  
Scope: B2 only, not an overall production/publication verdict

## Verdict: CONDITIONAL PASS

The algorithm repair is real and correct; the recorded numbers reproduce from
the source as read. The remaining defects are in the gate, not the result.

### Confirmed

- The atom-inventory attenuation cap and dead direct/time-averaged branch are
  genuinely removed. `gas_absorption_scale` is a returned compatibility
  diagnostic only and is never applied to attenuation.
- The absorbed-count-to-rate conversion divides by the same time-averaged H I
  density used for optical depth. The H relaxation and its time average are a
  correct C2-Ray-style closure.
- The H ledger is nonzero away from the fixed point and closes there; it is not
  a vacuous arithmetic identity, although it is related to the fixed-point
  residual in photon units.
- The He backward-Euler, collisional, recombination, dust-energy, and secondary
  terms are internally consistent by inspection.
- The preserved 12/16/20-iteration artifacts show geometric contraction:
  maximum residual `8.618e-4`, `1.436e-4`, `2.396e-5`. Twenty iterations is the
  first tested multiple of four below `1e-4`.
- P4/P5 reject fewer than 20 iterations and hard-gate species ledgers,
  primary-absorption closure, residual, finiteness, positivity, and fraction
  bounds.
- Solver A and Solver B are not aliases, but they share transport, H
  relaxation, cross sections, and recombination coefficients. Their
  differential is a wiring cross-check; the analytic radius is the independent
  physics reference.
- Dust and secondary controls activate real channels. The shadow control is a
  valid full-wrapper wiring check but intentionally reproduces the P1
  transport-only setup with inert gas.
- The zero-He divide-order repair in `multiphysics.py` is safe.

### Mandatory for B2 closure

1. Gate maximum fixed-point residual and H-ledger L1 error across all four
   Solver-A runs and Solver B, not only the baseline Solver-A run.
2. Narrow or evidence the helium claim. Every B2 run has zero helium, so B2
   does not validate He transport/chemistry despite the implementation being
   internally consistent by inspection.
3. Relabel cap-active fraction and minimum scale as structural regression
   invariants rather than measured physical diagnostics; the retired cap makes
   them constant by construction.

### Later-gate improvements, not B2 blockers

- Apply strict positive density floors to the sibling conservative solvers.
- Document the shared A/B core and consider tightening the differential.
- Add a reactive-gas shadow at 20 iterations.
- Add magnitude bands to dust and secondary fixture deltas.
- Make the B2 CLI itself reject fewer than 20 iterations.
- Add a relative fixed-point measure for nearly neutral cells and remove the
  hard-coded Solver-B 18 eV energy.

B2 PASS must not be interpreted as overall production or publication PASS;
later RT findings, physical assets, source models, yields, licensing, and live
RAMSES coupling remain independently gated.

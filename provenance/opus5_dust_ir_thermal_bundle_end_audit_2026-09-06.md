# Opus 5 DUST-2 bundle-end audit — 2026-09-06

## Scope and verdict

Read-only audit of the DUST-2 grain thermal balance and one-pass IR ledger in
`/gpfs/kjhan/LRD_JWST`.  The project objective was supplied explicitly:
production/publication-ready RT, stellar/AGN feedback, and dust physics, with
SNRT/P5 serving as an auditable static validation layer rather than the live
RAMSES implementation.

**Verdict: CONDITIONAL PASS**

The implementation is sound as a labelled static JAX/P5 candidate.  No
physics, ledger-leak, or double-counting blocker was found.  It must not yet
be promoted to transported IR, native Fortran, live RAMSES, or production
dust physics.

## Validated strengths

- F1–F9 from the Fable plan audit are implemented; F10 was correctly treated
  as optional.  Explicit out-of-band energy prevents false closure, including
  the cold-dust case.
- Kirchhoff power uses `4*pi*C_abs(E)*B_E(T)` with dimensionally consistent
  integration and exact group-edge endpoints.  Thermal mean photon energies
  are used and are validated inside their groups.
- The CMB bath uses the same emission table and `T_bg=2.7255/a`; the solver
  returns a JIT-safe range mask and never silently reports clipped solutions
  as valid.
- DUST-1 heating and the IR ledgers use the same subcycle timestep weights;
  the thermal path does not enter gas energy, primary radiation, or the DUST-1
  absorption/scattering ledgers.  Thermal-on/off P5 HII and gas energy were
  bitwise identical in the supplied evidence.
- Source-table SHA-256, dust mass/H, group edges, code manifest, and payload
  hashes are bound fail-closed.  Instrumentation is limited to scalar closure,
  range, and optical-depth diagnostics.

## Findings and required actions

1. **Moderate — hashed algorithm description is inaccurate.**  The sidecar
   says inverse/log-power interpolation, while runtime uses linear-in-power
   interpolation on a log-temperature grid and 32-step bisection.  Correct
   the description and regenerate the sidecar.
2. **Moderate — test gates are incomplete.**  Add a genuine multi-IR-group
   operator and photon-number test, compare DUST-1 absorption/scattering and
   heating ledgers with `atol=0`, and run a zero-dust thermal-flag control.
3. **Moderate — thermal/opacity binding has no negative test.**  Assert
   rejection for wrong source hash, wrong dust mass/H, and non-v3 opacity.
4. **Moderate — document the CMB excess spectral approximation.**  The
   implementation distributes `P(T)-P(T_bg)` using fractions at `T`, rather
   than the exact differential group split.  It conserves total energy but
   conservatively under-reports tracked IR near the CMB; this is acceptable
   only for the labelled candidate and must be resolved before transport.
5. **Moderate process — register DUST-2 in `provenance/INDEX.md`.**
6. **Low — provenance clarity.**  Record that the DUST-1 sidecar was
   regenerated only because code/payload hashes changed and that its numeric
   physics is unchanged.
7. **Low — avoid presenting energy closure and power residual as independent
   checks; they are algebraically the same residual in this one-pass path.**
8. **Low — document the 5 K floor (thermal path valid for z approximately
   above 0.83) and label the current optical-depth scalar as a conservative
   source-weighted upper bound.**
9. **Low — consider gating `n_H<=0` invalidity on dust-bearing cells and
   recording the declared closure tolerance attribute.**

## Disposition

Findings 1–5 are in-scope corrections before this bundle is committed.
Findings 6–9 are documentation/hygiene or later refinement and are recorded
without widening the bundle.  Deferred promotion gates remain: physical
absorption/depletion normalization, stochastic/PAH heating, dust–gas exchange,
source obscuration, IR self-absorption and recursive transport, native/live
coupling, MPI/restart qualification, and matched source/geometry convergence.

Recommended next bundle: **DUST-3 IR source promotion study**, only after the
CMB differential split is corrected or explicitly accepted in its plan, with
an IR group-resolution decision and a non-convergence exit.

## Operator disposition after audit

Findings 1–5 were addressed before commit: the builder and regenerated sidecar
now describe the actual interpolation/bisection; tests cover two IR groups,
per-group photon-number/energy closure, zero-dust operation, DUST-1 ledger
immutability, source-hash and mass mismatches, and non-v3 rejection; the CMB
excess approximation is documented; and DUST-2 is registered in
`provenance/INDEX.md`.  The focused suite was rerun successfully.  The
candidate status and all later physical/native/live promotion gates remain
unchanged.

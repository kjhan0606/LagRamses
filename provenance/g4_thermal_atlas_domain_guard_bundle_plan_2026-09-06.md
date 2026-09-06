# G4 thermal-atlas domain guard bundle plan

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Work location: `/gpfs`

## Objective

Remove a production-facing ambiguity in the G4 thermochemistry boundary. The
JAX interpolation primitive may retain edge clamping for deliberately bounded
offline controls, but a P5 run must not silently turn an epoch or initial gas
state outside the admitted thermal-atlas domain into an in-domain value.

## Implementation scope

1. Add a host-side `ThermalAtlas.validate_runtime_domain` check for scale
   factor, positive finite `n_H`, and positive finite initial temperature.
2. Call it immediately after the validated atlas and static input are loaded,
   before JAX transport or output creation.
3. Count an internal-energy temperature outside the admitted table range as a
   thermal bound hit during subcycling, preserving the existing fail-closed
   P5 validation path.
4. Extend the existing G4 dust runner with scale-factor and initial-temperature
   rejection paths. No new matrix, run directory, or validation framework is
   added.

## Explicit non-scope

This does not widen the atlas, extrapolate cooling data, change the metal-only
physics, or claim redshift coverage beyond the admitted table. The current
atlas remains the narrow `a=0.20849--0.20851` output-00017 bracket. Physical
DTM/depletion, source-cell convergence, live RAMSES coupling, and cosmological
production remain open gates.

## Acceptance

The normal derived-dust P5 control remains green; an out-of-range scale factor
and initial temperature are rejected before output creation; and the complete
G4 mapping/source/thermal gate remains green on CPU/JAX.

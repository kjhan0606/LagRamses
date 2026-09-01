# P0 Radiation Initialization Policy

Status: Draft v0.1  
Date: 2026-08-28

## Decision

Existing RAMSES-RT output is not required to start S_N development or to run
the first S_N calculation on a hydro snapshot. S_N creates and evolves its
own angular radiation field from the gas state and source catalogue.

The minimum snapshot input is AMR geometry, gas conservative variables, and
the stellar and/or sink source catalogues. RT checkpoint output is optional.

## Initial state

1. Initialize every S_N angular intensity to zero, unless a prescribed
   metagalactic UV background is enabled.
2. Deposit stellar and AGN luminosities directly into the S_N emissivity
   field at the first radiation step.
3. Start the initial chemistry test with a documented neutral H/He state, or
   with a controlled equilibrium initializer derived from density, temperature,
   and the selected UV background.
4. Evolve radiation, chemistry, heating, and radiation momentum forward from
   that state. No M1 radiation field is imported.

## Role of existing RT output

Existing M1 RT fields are useful only for:

- direct M1--S_N comparisons at matched gas and source states;
- optional chemistry initialization; and
- regression tests of the RAMSES source-export path.

Their absence must never block the standalone S_N solver, the rasterizer, or
the first static AGN radiation calculation.

## Consequence for candidate outputs

The SIDM AGN output may be used after its binary payload is independently
validated. Its lack of an RT checkpoint is not a defect for S_N. The missing
completion marker remains a separate integrity concern.

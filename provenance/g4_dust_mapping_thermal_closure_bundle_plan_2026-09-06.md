# G4 dust mapping and thermal-ledger closure bundle plan

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Work location: `/gpfs`

## Objective

Close the next bounded G4 input/ledger boundary for the project's final
objective: a production/publication-ready high-level radiation-transfer,
stellar/AGN-feedback, and dust path in `lagRamses`. This bundle verifies that
static dust abundance reaches SNRT from an explicit source and that the
one-pass dust thermal ledger is observable without being confused with gas
heating or H/He absorption.

## Implementation scope

1. Use one resolver for the yt staging boundary. A direct dust abundance is
   authoritative; a metallicity/DTM mapping is accepted only as a complete
   pair; a partial pair fails closed. No fields remain an explicit zero-dust
   pilot control, not an inferred physical abundance.
2. Preserve `metallicity_solar` and `dust_to_metal` in P5 output when present,
   and record the abundance contract/origin. The derived path must use the
   exact cell product and must not introduce a redshift law, depletion law,
   metallicity floor, or hidden rescaling.
3. Run one compact CPU/JAX gate covering the P4 mapping/read-back contract,
   source-to-dust binding regression, and the existing P5 physical/reference
   dust plus one-pass thermal ledger controls.

## Explicit non-scope

This is not approval of the Draine/WD01 mixture, an astrophysical DTM law,
depletion, grain-size distribution, stochastic heating, source obscuration,
IR recursive transport, live RAMSES dust state, dust force, or AMR/MPI/restart
coupling. The external opacity and thermal assets remain candidate/reference
inputs under their existing binding and status contracts.

## Acceptance

The gate must show:

- derived abundance round-trips as `metallicity_solar*dust_to_metal`;
- partial mappings are rejected at both the field-map and array resolver
  boundaries;
- P5 preserves the composition fields and origin in HDF5;
- the existing source/opacity/thermal ledgers close, the zero-dust control is
  unchanged, and IR remains explicitly `recorded_not_transport_reemitted`.

The result is a conditional engineering pass, not a physical G4 promotion.

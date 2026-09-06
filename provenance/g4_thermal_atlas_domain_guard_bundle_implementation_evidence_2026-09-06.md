# G4 thermal-atlas domain guard bundle implementation evidence

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Date: 2026-09-06

## Delivered changes

- `ThermalAtlas.validate_runtime_domain` checks the admitted scale-factor,
  density, and temperature axes without clamping.
- P5 invokes the check before constructing the JAX atlas or writing an output.
- Subcycle states whose internal energy implies a temperature below or above
  the admitted thermal range are counted as thermal bound hits and therefore
  cannot pass P5 validation.
- The existing `p5_dust_runner.py` now asserts both scale-factor and initial-
  temperature rejection, with no rejected output file left behind.

## Verification

```text
simulation/snrt/tests/run_g4_dust_closure.sh
```

The gate passed with:

```text
P4_INGESTION_OK format=v3 shape=4x3x2 sources=2 groups=2
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
G4_DUST_CLOSURE_PASS tests=3 mapping=explicit thermal=one_pass backend=cpu
```

The P5 test also passed the negative checks for `a=0.20848` and `T=10^10 K`.
The in-range atlas remains `a=[0.20849,0.20851]`,
`log10(n_H/cm^-3)=[-10,4]`, and `T=[10,10^9] K`.

## Limits

This is a fail-closed domain guard, not a thermal-atlas expansion or a
production redshift qualification. Physical dust and feedback promotion
still require broader audited source data, convergence, and live coupling.

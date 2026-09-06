# P0 completion record: static S_N RT foundation

## Delivered scope

- A static Cartesian, multigroup discrete-ordinates transport core with vacuum boundaries, explicit transport, and exact local absorption.
- Carlson level-symmetric S4, S6, and S8 angular rules, normalized for photon-number moments.
- Generic point-source deposition, with host float64 conversion of physical photon luminosities before transfer to the float32 JAX field.
- H I, He I, and He II Verner cross sections; H/He non-equilibrium bookkeeping and photon-budgeted coupling.
- CPU-executable B01 primordial Stromgren, B03 opaque-clump shadow, and B04 symmetric crossing-beam benchmarks.
- P0 input contract, benchmark matrix, radiation initialization policy, and TPU memory model already recorded in the P0 documents.

## Reference validation

Run from this directory:

```bash
JAX_PLATFORMS=cpu .venv/bin/python tests/p0_smoke.py
```

Baseline measured on 2026-08-28 with JAX 0.11.1:

- `Rion / RS = 0.4369` after 128 B01 steps.
- B03 shadow transmission relative to the transparent counterpart: `0.0000` at 160 steps.
- B04 symmetry-plane flux factor: `5.433e-09` after 64 steps.

## P0 exit criteria met

- No snapshot or RAMSES RT output is required to initialize the solver.
- S4/S6/S8 have the expected 24/48/80 directions, unit total weight, and zero first angular moment.
- The three P0 benchmark paths compile and execute on CPU with finite fields.
- cgs source and diagnostic quantities avoid float32 intermediate overflow.
- B01 applies the local absorption operator exactly, so its time step is limited by the directional transport CFL rather than cell optical depth.

## Deliberate deferrals to P1+

- Angular and spatial convergence studies, including S4/S6/S8 comparison at production resolution.
- Dust, scattering, IR re-emission, X-ray secondary ionization, and thermal/hydrodynamic energy coupling.
- Local implicit chemistry, TPU multi-device sharding, and production benchmark throughput.
- RAMSES/zoom-snapshot ingestion and all dual-AGN science applications.

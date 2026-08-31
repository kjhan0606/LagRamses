# P3 implicit chemistry and TPU execution

## Local implicit H/He recombination

After photon-driven primary and secondary ionizations, the P3 closure solves case-B recombination with a backward-Euler update. At fixed electron density the H II, He II, and He III fractions are analytic. The remaining one-dimensional electron-density constraint is solved by 24 fixed bisection iterations per cell.

This avoids a data-dependent Newton loop, preserves non-negative fractions, and maps directly to XLA. In the stiff pure-H validation (`alpha_B n_e dt=2.58`), the result is `x_HII=0.457505`, matching the analytic backward-Euler root.

## Spatial sharding

`XShardings` uses a `NamedSharding` mesh over the Cartesian x dimension:

- Intensity `[group, direction, x, y, z]` is split only in `x`.
- Group scalar fields `[group, x, y, z]` use the corresponding partition.
- Groups and directions remain replicated, so the sweep remains static and branch-free.
- The reference global stencil is retained; XLA SPMD is responsible for x-boundary halo communication.

The validation runs with two virtual CPU devices and obtains the same field as the single-device reference. It validates sharding semantics, not TPU throughput or inter-host networking. A production TPU run must use a real device mesh with an x dimension that divides the grid exactly.

## Reproduction

```bash
XLA_FLAGS=--xla_force_host_platform_device_count=2 JAX_PLATFORMS=cpu .venv/bin/python tests/p2_p3_validation.py
```

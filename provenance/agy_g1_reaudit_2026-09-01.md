# AGY G1 gate re-audit — final

- Date: 2026-09-01 (Asia/Seoul)
- Auditor: AGY / Gemini Antigravity CLI (`gemini-3.1-pro-high`)
- Scope: G1 native units, cumulative source contract, interpolation, field
  mapping, restart/repeat progress, and the native/JAX differential matrix
- Project root: `/gpfs/kjhan/LRD_JWST`
- Mode: read-only plan audit
- CLI session: `4842`

## Verdict

**PASS.** AGY judged the G1 exit criteria satisfied and cleared G2 table
generation. The previous high-severity coverage finding was closed by the
six-query differential matrix.

## Final findings

- `age_yr -> age_gyr`, the RAMSES code-time convention including explicit
  `aexp**2`, and Msun/erg/g cm/s conversions are explicit and consistent.
- `C(age_gyr + dt_gyr) - C(age_gyr)` is used for source increments. Negative,
  nonfinite, and out-of-domain cases fail closed; endpoint clamping is absent.
- Full Cartesian corner presence, duplicate tuples, nonfinite values, and
  cumulative monotonicity are audited before runtime use.
- The versioned NVAR=18 map covers density, momentum, energy, total metal,
  delayed cooling, and H–Fe fields without overlap.
- `indtab` progress is committed after deposition; repeat, restart, and abort
  behavior is covered by the native test. Atomic MPI/checkpoint coupling is
  explicitly deferred to G5.
- The Fortran and CPU-JAX differential tests now cover exact low/high corners,
  mass/metallicity/age boundary faces, and an interior point. The recorded
  maximum relative difference is approximately `2.13e-16`.

## Evidence

- `simulation/snrt/tests/run_g1_native_contract.sh` completed with
  `G1_NATIVE_CONTRACT_RUN_OK`.
- `simulation/snrt/data/g1_native_jax_differential.json` records six passing
  queries on JAX `0.11.1` CPU.
- `stellar_ramses_runtime.f90` syntax-compiled against the existing RAMSES
  module files. It is still an integration candidate, not a clean production
  executable.

## Limitations and out of scope

HDF5 integrity, physical yield-paper approval, RT/dust/AGN physics, live
MPI/AMR restart qualification, and G2–G7 were not audited. The external
backup mirror comparison was not used by AGY; the local source provenance and
syntax check remain recorded separately.

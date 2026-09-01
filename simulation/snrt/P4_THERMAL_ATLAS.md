# Offline metal thermal atlas

The runtime atlas is a provenance-enforced, UVB-free metal-only product. It
does not contain equilibrium primordial cooling or photoheating. SNRT evaluates
atomic H/He cooling from the live non-equilibrium ion fractions and adds local
S_N photoheating exactly once. The complete B1 rationale and validation are in
[`B1_THERMAL_COUPLING.md`](B1_THERMAL_COUPLING.md).

At runtime, `snrt_core.jax_thermal_atlas` brackets scale factor and performs
linear interpolation in `log n_H` and `log T` on a solar-metallicity table;
the result is multiplied analytically by cell `Z/Zsun`, including exactly zero
metallicity. The accepted v3 provenance contract fixes the component,
input and generator checksums, Grackle/data revisions, UVB exclusion, sign
convention, analytic metallicity application, and CMB-floor choice. Legacy v1
full-equilibrium atlases and the defective four-dimensional v2 atlas are
rejected by `read_thermal_atlas`.

The pinned source is
[`CloudyData_noUVB.h5`](../../external/grackle/CloudyData_noUVB.h5), revision
`928696482fbe15d9bac4382de6134d95568f099c` of the Grackle data repository.
Build the atlas with:

```bash
./.venv/bin/python tools/build_metal_thermal_atlas.py \
  --source-data ../../external/grackle/CloudyData_noUVB.h5 \
  --expected-source-sha256 0abe25cceeb5c0825381c5f17059982a9a2cdd27ce369a475c559fba6a8fa106 \
  --scale-factors config/p6_thermal_atlas_scale_factors.txt \
  --output data/production_metal_thermal_atlas_v2.h5
```

The atlas mean molecular weight is only a neutral-primordial staging aid.
Runtime heat capacity is computed from the evolved H/He state. Production
snapshots must provide usable pressure/temperature and cell-wise metallicity;
the atlas equilibrium-temperature fallback is not a physical substitute.

The metal CMB term continuously subtracts the source coefficient evaluated at
`T_CMB`. This deliberately removes Grackle revision `f93091f`'s optimization
that stops the subtraction above `100 T_CMB`, avoiding a finite cutoff step;
the deviation is recorded in the atlas provenance.

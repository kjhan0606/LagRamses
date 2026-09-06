# AGY audit request: B1 thermal coupling gate only

Date: 2026-09-01
Project root: `/gpfs/kjhan/LRD_JWST`
Auditor: AGY / Gemini Antigravity CLI, requested model `gemini-3.1-pro-high`

Perform an independent, read-only audit of **B1 thermal coupling only**. Do not
audit or summarize the complete G0--G7 program. Do not edit files, launch batch
jobs, copy/delete data, commit, fetch, pull, or merge. The GitHub `main` branch
is known to be one commit ahead of this checkout; report that as an integration
condition, not as a reason to broaden this gate audit.

The gate objective is to determine whether SNRT now has a physically defensible
and correctly wired thermal/ionization model in which:

1. live non-equilibrium H/He species control primordial atomic cooling and
   collisional ionization;
2. RT photoabsorption supplies photoionization and photoheating exactly once;
3. the tabulated contribution is metal-only, linearly scaled with metallicity,
   contains no UV background or photoheating, and applies its CMB floor with
   the intended sign;
4. species transitions, absorbed photons, and thermal energy have mutually
   consistent units/signs/ledgers;
5. provenance makes legacy equilibrium/UVB atlases fail closed.

Read at least these files and relevant imports/callers:

- `simulation/snrt/B1_THERMAL_COUPLING.md`
- `simulation/snrt/snrt_core/primordial_cooling.py`
- `simulation/snrt/snrt_core/implicit.py`
- `simulation/snrt/snrt_core/thermochemistry.py`
- `simulation/snrt/snrt_core/multiphysics.py`
- `simulation/snrt/snrt_core/conservative_primordial.py`
- `simulation/snrt/snrt_core/thermal_atlas.py`
- `simulation/snrt/snrt_core/jax_thermal_atlas.py`
- `simulation/snrt/tools/build_metal_thermal_atlas.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tools/p6_run_conservative_thermochemical_pilot.py`
- `simulation/snrt/tests/b1_thermal_coupling.py`
- `simulation/snrt/tests/p2_p3_validation.py`
- `simulation/snrt/tests/p4_thermal_atlas.py`
- `simulation/snrt/tests/p8_sharded_conservative.py`
- `simulation/snrt/data/b1_validate_conservative_primordial_thermal_atomic.json`
- `manifests/lrd_jwst_external_assets.json` records
  `grackle_cloudy_data_no_uvb` and `production_thermal_atlas_v1`
- `manifests/production_readiness_manifest_v1.json` thermal-atlas metadata

Independently inspect the actual HDF5 attributes/datasets and hashes when useful.
The pinned source is `external/grackle/CloudyData_noUVB.h5`; the generated atlas
is `simulation/snrt/data/production_metal_thermal_atlas_v1.h5`. Compare formulas
and data semantics against the locally recorded exact Grackle revision
`f93091ff8456962d7017a5bff7472945a30e3dad` rather than trusting documentation
claims. Tests and the stored convergence report are evidence, not proof.

Focus on scientific/algorithmic legitimacy and wiring. Check especially:

- cooling/heating signs and cgs volume-rate units;
- density powers and species factors in H/He processes;
- consistency between collisional ionization transitions and ionization-energy
  losses;
- case-B choice and coefficient consistency;
- electron-density closure and boundedness of the implicit solve;
- splitting/order errors between RT, chemistry, and thermal update;
- metal-table interpolation, redshift/scale-factor handling, CMB subtraction,
  and metallicity scaling;
- accidental double counting of primordial cooling, UVB, or photoheating;
- behavior at table/domain boundaries and unsupported states;
- CPU/JAX/sharded-path parity and diagnostics;
- whether the tests materially exercise the production P5/P6 paths.

Return a compact but detailed report with:

1. one verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED` for B1 only;
2. confirmed-correct design elements;
3. findings ranked critical/high/medium/low, with exact file and line references;
4. validation gaps versus actual algorithm defects clearly separated;
5. mandatory fixes for B1 closure;
6. explicitly deferred physics that belongs to later dust/RT production gates;
7. whether the atlas may be promoted from `independent_gate_audit_pending`, while
   separately respecting the unresolved data-license approval.

Do not equate B1 passage with overall production/publication readiness.

# AGY re-audit request: B1 thermal coupling closure only

Date: 2026-09-01
Project root: `/gpfs/kjhan/LRD_JWST`
Auditor: AGY / Gemini Antigravity CLI, requested model `gemini-3.1-pro-high`

Perform a read-only **B1-only re-audit**. Do not audit the complete G0--G7
program. Do not edit files, run batch jobs, copy/delete artifacts, commit,
fetch, pull, or merge. GitHub `main` is known to be one commit ahead of this
checkout; keep that as a later integration condition.

Read the original request and finding record:

- `provenance/agy_b1_thermal_coupling_prompt_2026-09-01.md`
- `provenance/agy_b1_thermal_coupling_audit_2026-09-01.md`

Then inspect the revised implementation and actual artifacts, especially:

- `simulation/snrt/snrt_core/thermal_atlas.py`
- `simulation/snrt/snrt_core/jax_thermal_atlas.py`
- `simulation/snrt/tools/build_metal_thermal_atlas.py`
- `simulation/snrt/tests/b1_thermal_coupling.py`
- `simulation/snrt/tests/p4_thermal_atlas.py`
- `simulation/snrt/tests/p4_hdf5_staging.py`
- `simulation/snrt/B1_THERMAL_COUPLING.md`
- `simulation/snrt/P4_THERMAL_ATLAS.md`
- `simulation/snrt/P5_THERMOCHEMISTRY.md`
- `simulation/snrt/data/production_metal_thermal_atlas_v2.h5`
- `simulation/snrt/data/b1_validate_conservative_primordial_thermal_v2.json`
- thermal records in `manifests/lrd_jwst_external_assets.json` and
  `manifests/production_readiness_manifest_v1.json`

Verify independently, rather than trusting the docs, that:

1. format v3 has no metallicity dataset/dimension and its physical table shape
   is exactly `(a, n_H, T)`;
2. host and JAX runtimes apply one analytic `Z/Zsun` multiplier, including
   exact `Z=0`, and negative/non-finite metallicity cannot become physical
   cooling/heating;
3. off-grid scalar/vector tests would have failed the old logarithmic
   interpolation and now exercise both host and JAX paths;
4. CMB coefficient subtraction is continuous at all tabulated temperatures;
5. provenance explicitly records that this is a deliberate physical-continuity
   deviation from Grackle revision
   `f93091ff8456962d7017a5bff7472945a30e3dad`, whose reference implementation
   skips subtraction above `100 T_CMB`;
6. old format-v2 atlas v1 is quarantined and cannot be loaded by production
   readers or selected by current batch/runtime paths;
7. no dimension/broadcast/sharding regression was introduced.

Recorded local evidence to verify includes B1/P4/P4-HDF5/P2-P3/P5-dust/P8
tests passing, deterministic 32-cubed coarse/fine convergence, atlas SHA256
`b1290d930b22ed049d6d3c5ed47ce56ecf3e0d2e693b39792740768e80fdf6ac`,
generator SHA256
`6fe00f80795ac948ed2386512011014932f559c48e7a1455116031852c0f5280`,
and convergence-report SHA256
`0f1bd1f115a88484bc9233dffa52939b30739bba2377050dfa1f0e5a938ef0cf`.

Return:

1. one B1-only verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCKED`;
2. disposition of each original critical/high/medium finding;
3. any new algorithm defect separately from validation gaps, with exact file
   and line references;
4. mandatory fixes, if any;
5. whether the code/provenance portion may clear
   `independent_gate_audit_pending` while the separate upstream data-license
   status remains unresolved.

Do not equate B1 closure with overall production/publication readiness.

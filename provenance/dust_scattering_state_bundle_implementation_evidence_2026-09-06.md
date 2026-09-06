# Dust scattering/state bundle implementation evidence (2026-09-06)

## Scope and disposition

This record covers D1--D3 of
[`dust_scattering_state_bundle_plan_2026-09-06.md`](dust_scattering_state_bundle_plan_2026-09-06.md).
The plan received Fable's `CONDITIONAL APPROVE`; required amendments are
recorded in [`fable_dust_scattering_state_bundle_plan_audit_2026-09-06.md`](fable_dust_scattering_state_bundle_plan_audit_2026-09-06.md)
and are applied here.

Disposition: `conditional_candidate`.  This is a static JAX SNRT closure and
not an astrophysical production approval.  Native Fortran dust channels, live
RAMSES force injection, IR re-emission, grain-temperature evolution,
aggregate STAR+AGN dust closure, and AMR/MPI/restart qualification remain
deferred.

## Implemented contract

- `snrt_dust_opacity_v1` remains an absorption-only `reference_control`.
- `snrt_dust_opacity_v2` remains a source-bound absorption-only candidate.
- `snrt_dust_opacity_v3` carries absorption and scattering cross sections,
  absorption/scattering-weighted energies, the scattering-weighted first two
  angular moments, and the diagnostic transport-corrected scattering cross
  section.  The phase status is explicit:
  `phase_isotropic_candidate`; source-bound status is
  `candidate_scattering_isotropic`; unbound control status is
  `reference_scattering_control`.
- The v3 sidecar records the isotropic-candidate momentum overestimate bound
  `1/(1-g)` per group and marks the `g=1` group explicitly unbounded.  The
  extinction/albedo residual is enforced at `1e-2`, while the measured raw
  moment residual and its `1e-4` rounding envelope are recorded.
- The Draine source parser now retains `albedo`, `⟨cosθ⟩`, and `⟨cos²θ⟩`.
  `K_abs × M_dust/H` remains authoritative because the published columns are
  rounded; the independent extinction/albedo relation is recorded with its
  measured residual.  The raw moment inequality allows the measured `9e-5`
  rounding envelope and rejects larger violations.
- The static snapshot stores `dust_relative_abundance_origin` as either
  `direct` or `metallicity_solar_times_dust_to_metal`; the latter requires both
  fields and an exact declared product within the input tolerance.  No
  dust-to-metal prescription is invented by the kernel.
- `advance_with_isotropic_scattering` performs the exact local constant-source
  isotropic solve after the existing explicit upwind transport.  Scattering
  is within-group and frequency-conserving.  `absorbed_intensity` remains
  absorption-only for the H/He and primary-absorption ledgers.
- P4/P5 expose separate cumulative scattered-photon and
  absorption/scattering/total momentum datasets.  The historical total-force
  dataset is retained as an explicitly labeled alias, preventing a consumer
  from adding a total to its components. Momentum uses physical `c` while
  transport continues to use its configured reduced light speed. Old sidecars
  cannot enable scattering and a v3 sidecar cannot run with scattering
  disabled.

## Files changed for this bundle

- `simulation/snrt/snrt_core/dust.py`
- `simulation/snrt/snrt_core/transport.py`
- `simulation/snrt/snrt_core/multiphysics.py`
- `simulation/snrt/snrt_core/thermochemistry.py`
- `simulation/snrt/snrt_core/snapshot.py`
- `simulation/snrt/tools/build_draine_dust_opacity.py`
- `simulation/snrt/tools/p4_run_transport_pilot.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tools/p4_stage_hdf5_level15.py`
- `simulation/snrt/tools/refine_static_rt_input.py`
- `simulation/snrt/tools/p4_attach_pilot_sources.py`
- `simulation/snrt/tests/dust_opacity.py`
- `simulation/snrt/tests/draine_dust_opacity.py`
- `simulation/snrt/tests/source_sed_dust_closure.py`
- `simulation/snrt/tests/p4_dust_runner.py`
- `simulation/snrt/tests/p5_dust_runner.py`
- `simulation/snrt/tests/p4_ingestion.py`
- `simulation/snrt/P2_MULTIPHYSICS.md`
- `simulation/snrt/P0_OUTPUT_CONTRACT.md`
- `simulation/snrt/P4_INGESTION.md`
- `simulation/snrt/P4_DUST_OPACITY.md`
- `simulation/snrt/P5_THERMOCHEMISTRY.md`

## Verification

The following commands passed on CPU JAX in `/gpfs/kjhan/LRD_JWST`:

```text
tests/dust_opacity.py                 DUST_OPACITY_TEST_OK
tests/draine_dust_opacity.py         DRAINE_DUST_OPACITY_TEST_OK rows=812 groups=9 scattering_groups=9
tests/source_sed_dust_closure.py     SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
tests/p4_ingestion.py                 P4_INGESTION_OK
tests/p4_dust_runner.py               P4_DUST_RUNNER_TEST_OK ... scattering=isotropic
tests/p5_dust_runner.py               P5_DUST_RUNNER_TEST_OK
tests/p4_hdf5_staging.py              P4_HDF5_STAGING_OK
tests/refine_static_input.py          REFINE_STATIC_INPUT_TEST_OK
```

The tests cover zero-scattering reduction, pure-scattering photon
conservation, beam redistribution, analytic mixed coefficients, source/
absorption steady-state invariance, S4/S8 scattering conservation, scattering
momentum, mixed P4/P5 execution, source-bound stellar and AGN synthetic SED
identities, HDF5 abundance-origin round-trip plus a negative consistency case,
staged-sidecar hash pinning, and both sides of the v1/v2 versus v3 activation
fence.

`git diff --check` and Python bytecode compilation of the changed modules also
pass.  An unrelated `tests/p2_p3_validation.py` invocation was not counted as
part of this bundle because its pre-existing device-count assertion requires
two visible JAX devices on the host.

## Staged candidate hash

The unbound physical Draine v3 reference control is:

`external/draine_wd01_rv31/p0_dust_opacity_rv31_photon_index1_scattering.json`

Its source table SHA-256 is
`b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c` and the
sidecar SHA-256 at evidence generation is
`61350545eea164c8db94ff830abd7f9e57cd7efc6f1e36f389d61627d364b9da`.

The sidecar is a reproducible reference control, not a promoted stellar/AGN
mixture.  A source-bound v3 sidecar must be rebuilt from the selected source
SED and will carry the source identity, raw input hash, edge hash, builder
hash, code manifest, and canonical payload hash.

## Opus 5 closure disposition

Opus 5 performed the 2026-09-06 bundle-end audit and returned
`CONDITIONAL PASS`; the full record is
[`opus5_dust_scattering_state_bundle_end_audit_2026-09-06.md`](opus5_dust_scattering_state_bundle_end_audit_2026-09-06.md).
The C1–C9 required repairs were applied: force-channel names and aliases are
explicit, the anisotropy bound is serialized and documented, all real dust
abundance producers preserve their origin, the independent opacity residual is
gated, the measured moment residual is recorded, the steady-state/S8/mixed
tests were added, the staged sidecar is hash-pinned, malformed v3 metadata
fails with `ValueError`, and a source-bound AGN v3 fixture was added. C10
documentation cleanup is complete for the implemented contract. The result
remains a `conditional_candidate`; it is not an astrophysical production
approval.

The repaired record was re-audited by Opus 5; see
[`opus5_dust_scattering_state_bundle_reaudit_2026-09-06.md`](opus5_dust_scattering_state_bundle_reaudit_2026-09-06.md).
That audit returned `CONDITIONAL PASS` and identified only the closure repairs
now applied in the working tree: source-bound moment provenance, legacy
non-zero-dust fail-closed reading, decisive malformed-v3 coverage, correct AGN
SED attribution, complete API/documentation semantics, producer-list accuracy,
and reproducible unbound-v3 provenance.

The final Opus 5 audit is recorded at
[`opus5_dust_scattering_state_bundle_final_audit_2026-09-06.md`](opus5_dust_scattering_state_bundle_final_audit_2026-09-06.md).
It returned `CONDITIONAL PASS`: the static DUST-1 candidate is physically and
technically closed for its declared scope, while the remaining IR/grain/live
hydro work stays deferred. Its only required R2/R5 documentation and test
items were applied afterward: the canonical staging documents now describe
format v3 and legacy fail-closed behavior, the format-version and legacy
negative tests are present, and the null-v3 test is decisive.

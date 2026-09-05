# F-P2.1 source-closure verification and honest metadata — implementation evidence — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Bundle: F-P2.1 repair of the F-P2 source-SED → photon-ledger → dust-closure
boundary. The driver approved this bundle on 2026-09-05.

Status: implementation and local evidence complete; Claude Opus 5 bundled
read-only end audit returned `CONDITIONAL PASS`. No commit or push is part of
this evidence capture.

## Delivered scope

### R1 — truthful explicit-SED metadata

`tools/p4_build_agn_photon_ledger.py` now constructs the explicit tabulated
SED and parameterized pilot metadata separately. The explicit path records:

- the validated source-SED identity and raw input hash;
- the actual configured group intervals as its support intervals;
- explicit normalization/escape-fraction units and interpretation; and
- intrinsic per-`L_bol` energy/photon moments separately from escaped moments
  and source totals.

It no longer emits the pilot Sazonov citation, the pilot Lyman normalization,
or the pilot 10-eV support claims in the explicit path. The built-in Sazonov
shape remains a null-identity `reference_control_parameterized_pilot`.

### R2 — independent Draine closure verification

`tests/source_sed_dust_closure.py` now independently recomputes

`∫ q_E κ_abs(E) dE / ∫ q_E dE`

and

`∫ E q_E κ_abs(E) dE / ∫ q_E κ_abs(E) dE`

from CSV/table rows, without calling the production integration routine. The
test deliberately uses SED samples, Draine samples, and group edges that do
not coincide. It also compares the result with a closed-form
`q_E = constant`, `κ_abs = 4e-21 E^(-1/2)` reference. The test covers the
source-bound v2 loader and runner, escape-fraction bookkeeping, and rejection
of tampered source, edge, Draine-table, and builder hashes.

### R3 — enforced provenance

The v2 dust sidecar now records and the loader re-hashes:

- the source SED input named by `source_sed_contract.input_path`;
- the exact group-edge file;
- the Draine source table; and
- the builder source file.

P4 and P5 validate the photon-ledger group-edge path/hash before output
creation and pass the expected source and edge contract into the dust loader.
Validated closure status/schema, source identity/hash, and group-edge hash are
written to output HDF5 attributes. Null-identity v1 reference controls remain
usable without a source-bound match and are rejected when attached to a bound
source ledger.

## Evidence executed on `/gpfs`

All commands below were run from `/gpfs/kjhan/LRD_JWST/simulation/snrt` with
the repository `.venv` and CPU JAX:

```text
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
MERGE_PHOTON_SOURCE_LEDGERS_TEST_OK components=2 sources=3 mixed_dust_gate=1
DUST_OPACITY_TEST_OK groups=1 weighted_energy_ev=7
DRAINE_DUST_OPACITY_TEST_OK rows=812 groups=9 max_consistency=8.078e-04
AGN_PHOTON_LEDGER_TEST_OK p0_groups=9 legacy_groups=5 hard_xray_group=positive subthreshold_opacity=zero edges=exact
P4_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
P0_SED_CLOSURE_OK synthetic_groups=4 metadata_groups=9 species=3 subthreshold_opacity=zero
STELLAR_PHOTON_LEDGER_TEST_OK sources=3 groups=9
AGN_NINE_GROUP_ARTIFACT_OK hard_q=3.52996e+51 hard_to_soft_q=0.219425 hard_supported_sed_fraction=0.0437846 hard_bolometric_fraction=0.00794605
```

The changed Python modules and focused test compile with `py_compile`, and
`git diff --check` is clean. The canonical AGN validator reports 27/27
criteria true (`AGN_NINE_GROUP_PASS`), including the top-level and nested
group-edge hash checks.

The canonical artifacts were regenerated with the current builder/runner and
validated together. Current SHA-256 values are:

```text
simulation/snrt/config/p0_photon_group_edges_ev.txt
  d28f78f1703730c6c0b9a7d183edfe0c5e6337979e737ce002a572b66fc53ff1
simulation/snrt/data/p4_pilot_agn_photon_ledger.csv
  bfdab21cffc3fb9a7d02e7d6a3f2e892fd6e7d87076022ecb6b2fda902e19a4e
simulation/snrt/data/p4_pilot_agn_photon_ledger.json
  86cbb8cd193137c6e7bbd66f8f48a36176bddade0f26e666e429c455f79c7fdb
simulation/snrt/data/p4_coeval_static_rt_input_agn9.h5
  c25720164947df6a7f64b01d898401cf00e1420e74281f233323d92d2316a8b2
simulation/snrt/data/p4_coeval_static_rt_input_agn9.json
  52aa155ce61e736088974ed71582690961b5406cd3eeb698c87b1f9e74dae6b7
simulation/snrt/data/p4_validation/p4_agn9_stage4_0p001myr.h5
  e591dc15daac0c59aecc3252c5254f217ab0ba678edf7663c630795834294171
simulation/snrt/data/agn_nine_group_external_assets.json
  38168ae0fd45425adf3b0a49d4094d5b31019dd2dc452d471175991691c819b3
simulation/snrt/data/agn_nine_group_validation.json
  7689641efb8136854a41985fba989fea62cc0abcb36131663645b7f174a9555e
```

The canonical pilot remains a reference-control artifact: its Sazonov-style
shape has null source identity and is not an approved AGN SED. No production
source activation, live RAMSES run, mixed STAR+AGN dust admission, or
publication claim was introduced.

## Deliberately deferred

Mixed STAR+AGN aggregate dust admission remains a later bundle, including
aggregate metadata provenance and `aexp` preservation. Physical stellar/AGN
SED and escape/obscuration approval, the `[40,120] M_sun` yield seam, dust
scattering, grain temperature, IR re-emission, destruction/growth, full
radiation pressure, live RT–RAMSES coupling, and production convergence are
outside F-P2.1.

## End audit disposition

Opus accepted the three F-P2.1 repair objectives and found no blocker. The
bundle remains conditional because the next closure-integrity bundle should:

1. hash and enforce the complete closure-code dependency set, including
   `snrt_core/sed.py` and the AGN-ledger closure dependencies;
2. bind or recompute-check the v2 sidecar payload, with an array-tampering
   negative test;
3. replace free-text dust status with a fixed vocabulary and propagate the
   validated dust metadata/table/builder hashes into P4/P5 outputs; and
4. add canonical explicit-SED nine-group coverage and a clean-tree attestation
   to the canonical validator.

These recommendations are recorded, not started. The full report is
[`claude_opus5_fp2_1_source_closure_verification_bundle_end_audit_2026-09-05.md`](claude_opus5_fp2_1_source_closure_verification_bundle_end_audit_2026-09-05.md).

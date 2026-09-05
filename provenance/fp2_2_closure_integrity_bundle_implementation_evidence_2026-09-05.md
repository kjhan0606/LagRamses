# F-P2.2 closure-record integrity and explicit canonical coverage — implementation evidence — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Bundle: F-P2.2, pre-approved after the F-P2.1 Opus 5 conditional-pass audit.
Status: implementation/evidence complete; bundled Claude Opus 5 end audit
recorded as `CONDITIONAL PASS`. No commit or push is included in this bundle.

## Delivered scope

### I1 — complete closure dependency manifests

`snrt_core/provenance.py` provides deterministic raw-file SHA-256,
role/path/hash manifest, and canonical JSON payload-hash primitives.

- Explicit AGN photon metadata binds the ledger builder, `sed.py`,
  `primordial.py`, and the integrity helper.
- Source-bound Draine v2 metadata binds the dust builder, `sed.py`, `dust.py`,
  and the integrity helper.
- The source boundary requires exactly those roles, exact resolved paths, and
  current file hashes; missing, substituted, extra, or stale entries fail
  closed.
- Explicit-source photon metadata re-hashes its source SED input before the
  closure is accepted.

### I2 — payload/status/output integrity

Explicit photon metadata and v2 dust metadata carry
`sha256_canonical_json_without_payload_sha256_v1` plus `payload_sha256`.
The v2 loader and explicit photon closure validator reject modified closure
arrays even when the input-file hashes are unchanged. Dust status is limited
to schema-consistent `reference_control` (v1) and
`candidate_source_sed_matched` (v2). P4/P5 outputs now record the metadata,
payload, source-table, and builder hashes; v1 reference controls leave the
source-bound fields empty by design.

### I3 — explicit canonical nine-group control

The synthetic, non-physical `data/p4_explicit_agn_sed_control.csv` covers the
same pinned P0 nine-group edges as the pilot. It is converted to an explicit
source ledger, rebound into a static input, and run through the same bounded
P4 transport path. `validate_agn_nine_group_ledger.py --source-mode explicit`
and its artifact test now validate this path separately from the historical
Sazonov parameterized pilot. The validator records `git_head`, a boolean
working-tree cleanliness flag, and a SHA-256 of the current porcelain status;
the canonical artifact therefore does not pretend that the dirty development
tree is a committed release.

## Evidence executed on `/gpfs`

Focused commands were run from `/gpfs/kjhan/LRD_JWST/simulation/snrt` with
the repository `.venv`, CPU JAX, and no live RAMSES job:

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
AGN_NINE_GROUP_EXPLICIT_ARTIFACT_OK mode=explicit criteria=all_true
```

The changed Python modules and focused tests compile with `py_compile`, and
`git diff --check` is clean. The final pilot and explicit validator runs are
required to be `AGN_NINE_GROUP_PASS`; each has 29 criteria after the explicit
mode/working-tree additions. The pilot artifact test performs the current-file
hash sweep; the explicit artifact test currently checks the canonical criteria
and attestation shape only, a gap recorded by the Opus audit and assigned to
F-P2.3.

The canonical asset manifest contains both the historical pilot control and
the explicit engineering control. The explicit source SED is a wiring fixture,
not an adopted AGN spectrum. The static inputs remain zero-dust controls; the
source-bound v2 dust path is tested in temporary fixtures and is not promoted
to production.

## Current boundary and exclusions

The bundled Claude Opus 5 audit is recorded at
`provenance/claude_opus5_fp2_2_closure_integrity_bundle_end_audit_2026-09-05.md`.
Its verdict is `CONDITIONAL PASS`: I1/I2 and the core explicit routing are
closed, while I3 remains conditional on repairing the stale explicit asset-
manifest digest, adding a non-vacuous explicit provenance/freshness test, and
making working-tree attestation enforcement symmetric. These repairs are the
proposed F-P2.3 bundle and await driver approval; no F-P2.3 work was performed
here.

This bundle closes metadata/code/payload integrity and canonical coverage; it
does not approve physical stellar or AGN SEDs, escape/obscuration, dust-to-
metal normalization, mixed STAR+AGN aggregate dust closure, `[40,120] M_sun`
yields, scattering, grain temperature, IR re-emission, destruction/growth,
full radiation pressure, live RT–RAMSES coupling, production convergence, or
publication claims.

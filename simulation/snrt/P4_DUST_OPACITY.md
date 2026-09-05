# P4 dust-opacity and heating contract

The static gas input already carries a dimensionless
`dust_relative_abundance` field, but that field is not an opacity. A run may
activate dust only when it also supplies a versioned JSON sidecar. The
historical `snrt_dust_opacity_v1` sidecar is a labeled reference control;
source-bound candidate runs use `snrt_dust_opacity_v2` and must match the
photon-ledger source identity.

## Sidecar schema

The required arrays are:

```json
{
  "schema": "snrt_dust_opacity_v1",
  "group_edges_ev": [1.0, 10.0],
  "absorption_cross_section_per_h_cm2": [1.0e-21],
  "absorption_weighted_energy_ev": [4.0],
  "reference_mixture": "declared dust mixture",
  "opacity_source": "versioned opacity table or calculation",
  "spectral_weighting": "declared photon-number SED or group weighting"
}
```

A source-bound v2 sidecar retains those fields and adds the identity of the
validated source SED:

```json
{
  "schema": "snrt_dust_opacity_v2",
  "status": "candidate_source_sed_matched",
  "source_sed_identity": "sha256 of the path-free SED contract",
  "source_sed_sha256": "sha256 of the raw SED input bytes",
  "source_sed_contract": {"schema": "snrt_source_sed_v1", "...": "..."},
  "group_edges_path": "path to the exact edge file",
  "group_edges_sha256": "sha256 of the exact group-edge file",
  "closure_code_manifest": [{"role": "...", "path": "...", "sha256": "..."}],
  "payload_hash_scheme": "sha256_canonical_json_without_payload_sha256_v1",
  "payload_sha256": "sha256 of the canonical sidecar payload",
  "source_sed_group_energy_fraction_of_lbol": ["per-group values"],
  "source_table": {"path": "Draine input", "sha256": "..."},
  "builder": {"path": "builder source", "sha256": "..."}
}
```

The loader re-hashes the source SED, edge file, Draine table, builder source,
and every expected closure-code manifest role. It also validates the canonical
sidecar payload self-hash, so changing the opacity or heating arrays without
rebuilding the sidecar is rejected. The v2 status is restricted to
`candidate_source_sed_matched`; v1 is restricted to `reference_control`. When
a photon ledger is source-bound, the P4/P5 runners also require the v2
identity, raw SED hash, and edge hash to match the ledger and require exactly
the same group edges. A v1 sidecar is intentionally rejected in that case; it
remains available for the null-identity pilot controls.

`absorption_cross_section_per_h_cm2` is the absorption cross section per H
nucleus for the declared reference dust mixture. The static input abundance
scales that reference mixture cell by cell. The weighted energy is the mean
energy of photons absorbed by dust in that group, using the same opacity
weighting used to produce the group cross section; it is not silently replaced
by the gas/stellar photon mean energy. The loader checks units, dimensions,
monotonic positive group edges, non-negative opacity, in-band weighted energy,
exact agreement with the photon-ledger group boundaries, and non-empty
provenance strings.

The implementation is in
[`snrt_core/dust.py`](snrt_core/dust.py). The P4 runner accepts
`--dust-opacity-metadata PATH`; without that option, a non-zero static dust
abundance is rejected instead of being silently ignored. With the option, the
validated opacity and the static abundance are passed to the JAX multiphysics
kernel. Dust receives the opacity-weighted absorbed photon energy and does not
alter H/He ion fractions.

The source-bound builder is invoked with the same explicit SED contract used
by the AGN ledger:

```bash
python tools/build_draine_dust_opacity.py \
  --source external/draine_wd01_rv31/kext_albedo_WD_MW_3.1_60_D03.all \
  --group-edges config/p0_photon_group_edges_ev.txt \
  --sed-table PATH/to/source_sed.csv \
  --sed-bolometric-fraction 0.95 \
  --output data/source_bound_dust_opacity.json
```

For every group it evaluates
`integral(q_E kappa_abs(E) dE) / integral(q_E dE)` and stores the
corresponding absorption-weighted photon energy. The SED, Draine table, and
group boundaries are joined at every sample/boundary and the base/refined
quadrature is compared at the declared `5e-6` tolerance; no missing support or
silent extrapolation is accepted.

## Scope gate

This stage implements dust absorption, local dust heating, and an
absorption-only momentum diagnostic. Dust scattering, temperature-dependent
re-emission into the IR group, and full radiation-pressure coupling remain
separate physics gates. The existing P4 transport artifacts remain zero-dust
controls. A physical Draine/WD01 Milky-Way \(R_V=3.1\) candidate is now staged
under `/gpfs/kjhan/LRD_JWST/external/draine_wd01_rv31`; its P0 sidecar was
generated with the auditable \(dN/dE\propto E^{-1}\) reference weighting. This
sidecar is not yet a source-specific production closure: the stellar/AGN
mixture and the cell dust-to-metal prescription still have to be selected and
recorded.

The synthetic loader and heating check is
[`tests/dust_opacity.py`](tests/dust_opacity.py).
The staged physical-table parser and nine-group Draine closure are checked by
[`tests/draine_dust_opacity.py`](tests/draine_dust_opacity.py).
The source identity, source-weighted v2 closure, runner binding, and negative
paths are checked by
[`tests/source_sed_dust_closure.py`](tests/source_sed_dust_closure.py).

The P4 output records the dust photon ledger separately at
`diagnostics/cumulative_dust_absorbed_photons_cm3`, the absorption-only force
rate at `rates/dust_momentum_rate_dyn_cm3`, and its time integral at
`diagnostics/cumulative_dust_momentum_g_cm2_s`. H/He photon closure excludes
the separately recorded dust absorption.

P4 and P5 also record `dust_opacity_metadata_sha256`,
`dust_payload_sha256`, `dust_source_table_sha256`, and
`dust_builder_sha256`. The last three are populated for v2 source-bound
sidecars and remain empty for the intentionally unbound v1 reference control.

The v2 path is an engineering/source-provenance closure, not an astrophysical
approval: the SED, dust-to-metal normalization, escape/obscuration model, and
all non-absorption dust physics remain later promotion gates.

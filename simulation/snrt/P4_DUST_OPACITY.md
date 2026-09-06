# P4 dust-opacity and heating contract

The static gas input already carries a dimensionless
`dust_relative_abundance` field, but that field is not an opacity. A run may
activate dust only when it also supplies a versioned JSON sidecar. The
historical `snrt_dust_opacity_v1` sidecar is a labeled absorption-only
reference control; source-bound absorption candidates use
`snrt_dust_opacity_v2`. The scattering candidate is
`snrt_dust_opacity_v3` with `phase_function=phase_isotropic_candidate` and
status `candidate_scattering_isotropic`; source-bound v2/v3 runs must match
the photon-ledger source identity.

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

A v3 sidecar retains the v2 provenance and adds the Draine scattering channel:

```json
{
  "schema": "snrt_dust_opacity_v3",
  "status": "candidate_scattering_isotropic",
  "phase_function": "phase_isotropic_candidate",
  "scattering_cross_section_per_h_cm2": ["per group"],
  "scattering_weighted_energy_ev": ["per group"],
  "scattering_angle_cosine": ["scattering-weighted <cos(theta)>"],
  "scattering_angle_cosine_squared": ["scattering-weighted <cos^2(theta)>"],
  "transport_corrected_scattering_cross_section_per_h_cm2": ["per group"],
  "isotropic_candidate_momentum_overestimate_factor": ["1/(1-g), or null when unbounded"],
  "isotropic_candidate_momentum_bound_unbounded": ["per group boolean"]
}
```

The loader re-hashes the source SED, edge file, Draine table, builder source,
and every expected closure-code manifest role. It also validates the canonical
sidecar payload self-hash, so changing the opacity or heating arrays without
rebuilding the sidecar is rejected. The v2 status is restricted to
`candidate_source_sed_matched`; v1 is restricted to `reference_control`; and
v3 is either the unbound `reference_scattering_control` or the source-bound
`candidate_scattering_isotropic` status. When a photon ledger is source-bound,
the P4/P5 runners require the sidecar identity, raw SED hash, and edge hash to
match the ledger and require exactly the same group edges. Old v1/v2 sidecars
cannot enable scattering, and a v3 sidecar cannot run with scattering disabled.

The Draine extinction/albedo residual is required to be at most `1e-2`; the
measured raw moment inequality residual is recorded with its rounding envelope.
The isotropic candidate is deliberately not an HG or delta-Eddington closure.
For a measured scattering anisotropy `g`, applying the full scattering opacity
can overstate the transport/momentum effect by at most `1/(1-g)` relative to
the transport-corrected coefficient. A group with `g=1` is marked unbounded;
the candidate therefore remains a labeled control, not a physical phase
function approval.

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
`--dust-opacity-metadata PATH` and `--dust-scattering {off,isotropic}`;
scattering defaults to `off`. Without the metadata option, a non-zero static dust
abundance is rejected instead of being silently ignored. With the option, the
validated opacity and the static abundance are passed to the JAX multiphysics
kernel. Dust receives the opacity-weighted absorbed photon energy; scattering
is within-group and does not alter H/He ion fractions.

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

This stage implements dust absorption, local dust heating, candidate isotropic
within-group scattering, and separate absorption/scattering momentum
diagnostics. Temperature-dependent re-emission into the IR group and full
radiation-pressure coupling remain separate physics gates. The existing P4
transport artifacts remain zero-dust controls. A physical Draine/WD01
Milky-Way \(R_V=3.1\) candidate is now staged
under `/gpfs/kjhan/LRD_JWST/external/draine_wd01_rv31`; its P0 sidecar was
generated with the auditable \(dN/dE\propto E^{-1}\) reference weighting and
the scattering-enabled v3 sidecar records the raw angular moments. These
sidecars are not yet astrophysical production approvals: the stellar/AGN
mixture, cell dust-to-metal prescription, source obscuration, and phase-function
promotion still have to be selected and recorded.

The synthetic loader and heating check is
[`tests/dust_opacity.py`](tests/dust_opacity.py).
The staged physical-table parser and nine-group Draine closure are checked by
[`tests/draine_dust_opacity.py`](tests/draine_dust_opacity.py).
The source identity, source-weighted v2/v3 closure, runner binding, and
negative paths are checked by
[`tests/source_sed_dust_closure.py`](tests/source_sed_dust_closure.py).

The P4 output records the dust photon ledgers separately at
`diagnostics/cumulative_dust_absorbed_photons_cm3` and
`diagnostics/cumulative_dust_scattered_photons_cm3`. The total dust force rate
is exposed at `rates/dust_total_momentum_rate_dyn_cm3` (the
legacy-compatible `rates/dust_momentum_rate_dyn_cm3` is an alias with the same
total semantics), while the components are exposed at
`rates/dust_absorption_momentum_rate_dyn_cm3` and
`rates/dust_scattering_momentum_rate_dyn_cm3`, with separate time integrals.
Consumers must use the total or sum the two components, never both.
H/He photon closure excludes the separately recorded dust absorption and
scattering channels.

P4 and P5 also record `dust_opacity_metadata_sha256`,
`dust_payload_sha256`, `dust_source_table_sha256`, and
`dust_builder_sha256`. The last three are populated for v2/source-bound-v3
sidecars and for the pinned unbound v3 reference control; they remain empty for
the historical unbound v1 reference control.

The local exponential-integral helpers use a float32 thin-cell series branch;
the sidecar contract is validated in float64 before runtime casting. The
source-table provenance retains the extinction/albedo consistency residual;
no separate group-averaged extinction array is used by transport, which
consumes the explicit absorption and scattering arrays.

The v2/v3 paths are engineering/source-provenance closures, not astrophysical
approvals: the SED, dust-to-metal normalization, escape/obscuration model,
anisotropic phase function, and live hydro coupling remain later promotion
gates. The separate DUST-2 thermal/IR candidate is documented in
[`P6_DUST_THERMAL.md`](P6_DUST_THERMAL.md); it does not change this opacity
contract or recursively transport its recorded IR source.

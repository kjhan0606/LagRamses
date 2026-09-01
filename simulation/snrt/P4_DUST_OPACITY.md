# P4 dust-opacity and heating contract

The static gas input already carries a dimensionless
`dust_relative_abundance` field, but that field is not an opacity. A run may
activate dust only when it also supplies a versioned
`snrt_dust_opacity_v1` JSON sidecar.

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

The P4 output records the dust photon ledger separately at
`diagnostics/cumulative_dust_absorbed_photons_cm3`, the absorption-only force
rate at `rates/dust_momentum_rate_dyn_cm3`, and its time integral at
`diagnostics/cumulative_dust_momentum_g_cm2_s`. H/He photon closure excludes
the separately recorded dust absorption.

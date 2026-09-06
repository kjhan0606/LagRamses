# F-P2 source SED and dust spectral-closure implementation evidence — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Bundle: F-P2 source SED → photon ledger → dust closure, parent `bd0411b`.
Status: implementation/evidence complete; candidate engineering path only;
Opus 5 end-of-bundle audit returned `CONDITIONAL PASS`.

## Delivered boundary

The bundle now has one validated source-spectrum contract in
[`simulation/snrt/snrt_core/sed.py`](../simulation/snrt/snrt_core/sed.py):

- CSV columns are `energy_ev` and `energy_fraction_per_ev`.
- The energy grid is finite, positive, strictly increasing, and the energy
  fraction is finite and non-negative.
- The integrated fraction is checked against the declared `L_bol` fraction.
- `q_E = f_E / (E_eV * 1.602176634e-12 erg/eV)` is integrated on a union grid
  containing every source sample and every configured group boundary.
- A path-free identity includes the raw SED hash and contract fields; the
  path is retained only as provenance.

The explicit AGN converter in
[`simulation/snrt/tools/p4_build_agn_photon_ledger.py`](../simulation/snrt/tools/p4_build_agn_photon_ledger.py)
now accepts `--sed-table` and `--sed-bolometric-fraction`. It serializes the
source identity/hash, support, per-group energy fraction, photon moments, and
the H I/He I/He II absorber-weighted closure. The built-in Sazonov-style
converter remains a null-identity `reference_control_parameterized_pilot`.

The Draine builder in
[`simulation/snrt/tools/build_draine_dust_opacity.py`](../simulation/snrt/tools/build_draine_dust_opacity.py)
now emits source-bound `snrt_dust_opacity_v2` metadata when given the same SED.
For each group it computes
`∫q_E κ_abs(E)dE / ∫q_E dE` and the corresponding absorbed-photon mean
energy, while retaining the Draine input and exact group-edge hashes.
The existing v1 `E^-1` sidecar remains a reference control.

P4/P5 runners pass the photon identity/hash into the dust loader. A v2 sidecar
must match it and the exact group edges; a v1 sidecar is rejected for a
source-bound photon ledger. Output HDF5 provenance records the binding status
and schema. The mixed-ledger tool records every component identity, creates a
new aggregate identity and gas closure, and sets
`component_only_sidecars_allowed=false` so a STAR-only or AGN-only dust
sidecar cannot be attached to a mixture.

## Evidence executed on `/gpfs`

The following focused tests passed after the implementation:

```text
SOURCE_SED_DUST_CLOSURE_TEST_OK source_bound_v2=1 agn_explicit_sed=1
MERGE_PHOTON_SOURCE_LEDGERS_TEST_OK components=2 sources=3 mixed_dust_gate=1
DUST_OPACITY_TEST_OK
AGN_NINE_GROUP_ARTIFACT_TEST_OK
```

The focused suite also retained passing P4/P5 dust-runner, Draine-table,
AGN-ledger, P0 SED-closure, and stellar-ledger checks. Negative coverage
includes incomplete SED support, source-identity mismatch, v1-for-bound-source
rejection, group-edge mismatch, source-ID collision, and refusal to create
outputs on rejected mixed-ledger inputs.

The canonical nine-group transport artifact was regenerated only because the
runner now records the source/dust binding attributes. The static HDF5 input
remained byte-identical; the transport output and its manifest hash were
updated together. The current relevant SHA-256 values are:

```text
data/p4_pilot_agn_photon_ledger.json       66e4a2968be9db8f770a3e40e3724150ef8efe4ed4218c7003296719d4530088
data/p4_coeval_static_rt_input_agn9.json   e6f60326138a2e8833ac100341940b18f6777369680cd69b2cc0b0d75f3a1812
data/p4_validation/p4_agn9_stage4_0p001myr.h5
                                             47a5958e24414e3bc56ca4ee90c6ae7477a43ea4d36d992b548184f544552b0e
data/agn_nine_group_external_assets.json    d08edc0abc63759a4b30f0ca4b04681ec9d2bed19516ce3557b1a14b0c206889
data/agn_nine_group_validation.json         26c5acd8008e1921ba413a037fea14830535b066e0eb8529747f16c823155f03
```

No live RAMSES job, production source activation, or publication artifact was
authorized by this bundle.

## Limits and deferred work

No BPASS/CCSN/AGB physical source was promoted, no AGN SED or obscuration model
was selected, and the `[40,120] M_sun` yield seam remains in the medium-term
review plan. Dust scattering, grain temperature, IR re-emission,
destruction/growth, full radiation pressure, live RT–RAMSES coupling, and
production-scale convergence remain later gates. The source-bound sidecar is
therefore `candidate_source_sed_matched`, not `approved`.

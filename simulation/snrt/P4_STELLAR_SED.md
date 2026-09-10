# P4 stellar SED and photon-ledger contract

The native transitional checkpoint now has a decoded stellar metadata
catalogue, but it does not contain a stellar spectrum. The converter therefore
keeps the SED choice as an explicit input and refuses to invent a luminosity
from age, mass, or metallicity alone.

## Input table

`tools/p4_build_stellar_photon_ledger.py` accepts a rectangular CSV table with
these required columns:

```text
age_myr,metallicity_solar,energy_ev,photon_rate_per_msun_per_ev_s
```

The last column means photons s^-1 eV^-1 Msun^-1 and is normalized to initial
stellar mass. An energy or luminosity spectrum is intentionally not accepted;
the energy-to-photon conversion must happen in the auditable SED preparation
step before this converter is run. The table must contain every combination
of its age, metallicity, and energy samples, with positive age, metallicity,
and energy and a non-negative photon rate.

The historical RAMSES-RT source tree was also checked for a reusable SED
directory. Its reader expects `RAMSES_SED_DIR` (or `sed_dir`) containing
`metallicity_bins.dat`, `age_bins.dat`, and the Fortran-unformatted
`all_seds.dat`; the latter stores wavelength spectra in solar-luminosity per
Angstrom per solar mass, not this converter's photon-number-per-eV contract.
No such data files were found in the migrated or registered external paths.
The presence of the historical reader source therefore does not certify a
usable SED asset.

Age and metallicity are interpolated linearly in log10 coordinates. The
photon-number spectrum is interpolated linearly, preserving tabulated zeroes.
Out-of-range sources are rejected by default; `--clamp-table-range` is an
explicit, recorded alternative. A zero native metallicity requires the
explicit `--metallicity-floor-solar` option. The solar mass fraction, selected
mass field, source scale factor, and source-side escape factor are all written
to the metadata sidecar.

## P0 groups and closure

The default boundaries are recorded in
[`config/p0_photon_group_edges_ev.txt`](config/p0_photon_group_edges_ev.txt):

```text
0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500, 2000, 10000 eV
```

These are nine groups. Each source spectrum is integrated with the trapezoid
rule on the SED grid augmented by every group boundary. The aggregate spectrum
weighted by the selected stellar mass is passed through the existing
`sed_weighted_group_closure`, so the sidecar contains the same photon-number
weighted H I, He I, and He II cross sections and photoelectron excess energies
used by the RT microphysics.

A stellar-only table may legitimately emit zero photons in a hard group. Such
groups are retained in the CSV with zero `q_group_N_s`; the sidecar records
`empty_source_group_zero_photons` and uses a geometric-mean energy with zero
opacity only as an inactive transport placeholder. That placeholder must not
be reused when a later AGN or combined source ledger injects photons into the
group; the combined source population needs its own aggregate closure.

The implementation is
[`tools/p4_build_stellar_photon_ledger.py`](tools/p4_build_stellar_photon_ledger.py),
and its synthetic contract test is
[`tests/stellar_photon_ledger.py`](tests/stellar_photon_ledger.py). A successful
run emits the standard contiguous `q_group_0_s` through `q_group_8_s` columns,
which are readable by `snrt_core.source_ledger` without further SED inference.

## Current status

The official BPASS v2.2.1 binary `imf135_300` HDF5 product is now staged at
`/gpfs/kjhan/LRD_JWST/external/bpass_v2.2.1` with a verified checksum. The
direct adapter
[`tools/p4_build_bpass_stellar_photon_ledger.py`](tools/p4_build_bpass_stellar_photon_ledger.py)
reduces its 100,000-wavelength node spectra to group moments before applying
log-age/log-metallicity interpolation. It generated the candidate ledger
[`data/feedback_transition_phase0_output_00011_bpass_stellar_photon_ledger.csv`](data/feedback_transition_phase0_output_00011_bpass_stellar_photon_ledger.csv)
for all 42,342 native stars.

That ledger is deliberately marked `candidate_bpass_stellar_photon_ledger`,
not a science input. The run explicitly used BPASS's lowest grid metallicity
as a floor for 338 zero-metallicity particles, clamped 42,004 sources below
the BPASS metallicity range, clamped 178 sources younger than 1 Myr, padded
the 0.01--0.123984 eV interval with zero photons, and assumed the HDF5 spectra
are normalized per initial stellar Msun. Escape fraction was set to the
baseline 1.0. These assumptions are all recorded in the metadata sidecar and
must be revisited before a STAR+AGN production merge.
The staged HDF5 adapter and candidate artifact are checked by
[`tests/bpass_stellar_ledger.py`](tests/bpass_stellar_ledger.py).

## Native evolving stellar source (2026-09-10)

The live Fortran source is separate from the historical P4 snapshot CSV
adapter above. [`tools/build_bpass_native_sed.py`](tools/build_bpass_native_sed.py)
exports the pinned Galacticus BPASS spectrum, whose upstream converter
already removed the original 1e6-Msun burst normalization. It integrates
photons per second per **initial** Msun. The native model uses linear age
in Myr and linear metallicity **mass fraction**, not log coordinates.
Intervals crossing age knots are integrated exactly for that declared model;
0--1 Myr explicitly holds the first spectrum. There is no additional age/Z
extrapolation, IMF reweighting or invented flux outside 1--100000 Angstrom.

The existing v2 asset and default exporter output are unchanged. Optional
`--energy-moments` produces v3 with both Q and energy rate (eV/s/initial Msun)
from the same positive piecewise-linear photon integrand. Q uses trapezoids;
energy analytically integrates that same integrand times hc/lambda. Escape
fraction applies to both once. Native age/Z interpolation and time integration
operate separately on Q and E; the injected mean is the final interval E/Q.

Select the retained native v3 comparison with:

```sh
export SNRT_ALLOW_REFERENCE_CONTROL=1
export SNRT_GROUP_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
export SNRT_STELLAR_SED=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_stellar_sed_bpass_independent_v3.nml
export SNRT_SPECTRAL_MODEL=hhe_maxent64_v1
```

This requires the existing SNRT-enabled **DUST_LIVE=0, CHIMES=0** band build
and the other normal RT/run settings. It is not an independent launch recipe.
The paired source deposits actual BPASS energy while preserving the fixed
state encoding `E = Eref*N + shift`. FP32 packet rounding is bounded using
actual source energy, not Eref. The v3 loader rejects nonphysical moments,
fixed spectral mode and energy-discarding callers. MPI/restart identity binds
the full Q/E table and semantics; changing even energy alone rejects restart.

This remains an **independent radiation-population comparison**: BPASS's
binary IMF/support do not match the selected feedback population. It does
not establish common-population eleven-element returns, a complete stellar
SED shape within each group, dust/CHIMES spectral coupling, or production
approval for combined stellar+AGN physics. See
[implementation and bounded MPI2 restart evidence](../../provenance/snrt_stellar_energy_implementation_2026-09-10.md).

# Stellar and dust asset baseline

Status: candidate inputs staged on `/gpfs`; neither candidate is promoted to a
science run.

## Stellar population

The selected raw asset is
`external/bpass_v2.2.1/SSP_Spectra_BPASSv2.2.1_bin-imf135_300.hdf5`, from the
[BPASS v2.2.1 Zenodo record](https://zenodo.org/records/6338460). The file is
the binary-star product with the broken IMF slopes -1.30 over 0.1--0.5 Msun
and -2.35 over 0.5--300 Msun. Its Zenodo MD5 is
`eac7f8ef432ff86b3d7ab8dc32c1b8b1`; the staged SHA256 is
`b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3`.

The HDF5 axes are 13 log-relative-solar metallicities, 51 ages in Gyr, and
100,000 wavelengths in Angstroms. `spectra` declares `L_sun/Hz`. The adapter
converts a spectrum with
`q_E = L_nu * L_sun / (h * E_eV)`, reduces node moments, and interpolates the
moments in log age and log metallicity. It requires explicit acknowledgement
that the SSP is per initial stellar Msun. The native `output_00011` candidate
also requires explicit young-age clamping and metallicity flooring because
the native gas is mostly below the BPASS grid; these decisions are retained in
the candidate metadata.

## Dust opacity

The selected raw table is Draine's WD01 carbonaceous-silicate Milky-Way
`R_V=3.1` D03-renormalized table
[`kext_albedo_WD_MW_3.1_60_D03.all`](https://www.astro.princeton.edu/~draine/dust/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all).
Its staged SHA256 is
`b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c`.
The opacity builder uses the table's declared `K_abs` and `M_dust/H`, checks
the rounded `C_ext/H * (1-albedo)` consistency, and produces the P0 sidecar
`external/draine_wd01_rv31/p0_dust_opacity_rv31_photon_index1.json` with
photon-number weighting proportional to `E^-1`. Its staged SHA256 is
`bde87ee9a4785bd47903c64cfff4a8609a48d8f8efe1fc42272944d67988bb07`.

The sidecar is absorption-only. Dust scattering, IR re-emission, grain
temperature evolution, and a source-specific stellar/AGN spectral weighting
remain separate promotion gates.

## Coevality gate

The BPASS stellar candidate is tied to the stopped native checkpoint at
`a=0.148540709098256`. The available AGN candidate is tied to output 00017 at
`a=0.20849776467675274`; the strict merger rejects that combination. A merged
ledger exists only as a labeled `noncoeval_integration_control`, with AGN IDs
shifted by `+1000000`. It is not a physical snapshot.

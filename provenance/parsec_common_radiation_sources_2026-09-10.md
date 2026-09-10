# Matched PARSEC photon input inspection

Workdir `/gpfs/kjhan/LRD_JWST`, LagRamses. This is input inspection,
not a completed native common-population source.

Primary paper https://arxiv.org/html/2501.12917v2 (section3.5); author
data https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/ .
The paper describes ATLAS9/PHOENIX/WM-basic/PoWR atmospheres; the public
QH files supply five integrated photon-number constraints, not the full
spectra or their energy moments. OII means O+, not a138-eV group.

Matching public files acquired without authentication, using a per-command
`curl --insecure` because the author's TLS chain still fails verification.
No credentials or global TLS changes. Hashes bind bytes, not independent
author authentication. Original physical inputs are retained.

| Archive in `.parsec-sources.9ZIUlJ/` | SHA256 |
| --- | --- |
| Z0.008_Y0.263_photons.zip | f6c78dc93e9aa9118dca8733983251e47040dd97cb3afa55230a6f981029a183 |
| Z0.014_Y0.273_photons.zip | 3f0caaa8f29fd7707f1538ca9fc8694bb0b5ad475db6350a789db59d4f0ec8d6 |

All90 feedback mass/Z nodes have photon tables. Across every inspected row,
QWerner >= QHI >= QHeI >= QOII >= QHeII. A64-point positive Gauss-Legendre
Planck photon-shape integration, normalized separately to their differences
on nominal11.2/13.6/24.6/35.12/54.4-eV edges, gives high-band energy divided
by tabulated bolometric luminosity in [2.76e-23,.98383073] atZ.008 and
[3.16e-24,.98098478] atZ.014. Thus this initial feasibility calculation
does not require a negative sub-LW residual. This is NOT evidence that the
reconstructed energy spectrum equals the original atmosphere spectrum.

First photon ages range12.36--235.94yr atZ.008 and13.10--221.39yr atZ.014.
The photon endpoints precede feedback track endpoints. Most differences
are only a few1e-5 of lifetime; the maximum atZ.014 is432.557yr (14Msun).
However twoZ.008 nodes are materially different:

| M/Msun | Missing end interval/yr | Lifetime fraction | Missing track bolometric energy fraction | final logT(photon)-logT(track) |
| --- | --- | --- | --- | --- |
| 60 | 5137.136 | .00131306 | .00206450 | -.32786254 |
| 70 | 11876.895 | .00327191 | .00497347 | -.30132679 |

Their endpoint effective temperature changes by about a factor2. A small
bolometric fraction does NOT bound HeII photon error. Therefore these two
tails must not be silently filled by the planned small-sampling last-state
hold or called measured radiation. Check the combined author archive for
more complete source rows before choosing a declared missing-tail closure.
No source data have been modified to remove this discrepancy.

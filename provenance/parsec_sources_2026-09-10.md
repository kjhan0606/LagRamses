# PARSEC input inspection (not implementation completion)

Primary paper: https://arxiv.org/abs/2501.12917 (Costa et al. 2025).
Author database: https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/ .
Files were downloaded under `.parsec-sources.9ZIUlJ/` in this project.
The author's TLS chain failed verification; public unauthenticated downloads
used a per-command `curl --insecure` fallback. No credentials were sent and
no global TLS configuration was changed. Hashes pin acquired bytes; they are
not independent author authentication. Paper and public index were checked.

| Asset | SHA256 |
| --- | --- |
| all_ejecta.zip | 49c0f0dac42ffb9643afc82eefe28793b44556cf6bbf0ed7d6df651a5aaa9b63 |
| Z0.006_Y0.259_tracks.zip | 2705d921e016780a892b526cf33c833175cce53361ea4e13fe826278b99df772 |
| Z0.006_Y0.259_photons.zip | bdaf060fb592e6161f5aa07fbd036125666b9c59cfd3d8e4e2b586f358bdceac |
| Z0.008_Y0.263_tracks.zip | aa93c7c44c05fe30411c89b6078329e90716b555859163226c415cff3c631c8f |
| Z0.014_Y0.273_tracks.zip | 83c31adb07b1d2280745c64b4dd601ae865ed73b90e34c0b0facd505410855e7 |
| readme.html | d4dcbc8db36e6d71780a16d9eb21a260d9ec5db3fbc23e8bcab67be247525e31 |

The full track archive is 1,282,140,114 bytes and full photon archive
115,670,480 bytes; only the relevant single-Z archives were downloaded.
These are original physical input archives, not disposable simulation dumps.

## Findings

- Separate wind and total tables: terminal isotope mass = total minus wind.
  Every inspected terminal isotope difference is nonnegative (13 Z branches).
- Full-chemical closure residual `(Min-Mbar-sum(ejecta))/Min` ranges from
  about -0.0582% to +0.1978%. This is not automatically a roundoff estimate.
  The eleven tracked elements must fit the gas return separately; do not
  normalize away a missing physical inventory.
- At Z=0.008 and 0.014, all eleven-element wind and terminal sums fit the
  baryonic mass return. Their remaining mass can remain untracked. Some
  other Z branches have small negative residuals from subtracting printed
  masses or terminal components and require explicit treatment.
- Mrem is gravitational, Mbar is baryonic. Neutrino losses do not supply gas.
- Terminal fate/yields cover >=14 Msun; the 2--9 Msun tracks terminate at
  early AGB, not complete TP-AGB envelope return. Nonrotation is a model
  coordinate, not a general rotation average.
- Fe is decay-inclusive; this is not a fresh radioactive inventory. Wind
  initial-mixture entries missing from some nuclear columns cannot be called
  physically absent isotopes. Eleven-element projection keeps residual mass.
- Photon files contain only five cumulative photon rates, not full spectra
  or the nine native Q/E pairs. They cannot directly close group 1. Their
  sampling and endpoint ages also differ slightly from the distributed
  evolutionary tracks (e.g. Z=.006 120 Msun end ages differ by 33.4 yr).
- Track surface columns cover H/He/C/N/O/Ne/Mg; initial heavy abundances
  supply unchanged Si/S/Ca/Fe under the author's wind prescription. Nuclear
  minor-isotope omissions are not a fully time-resolved isotope network.

## Phase-history numerical inspection

For Z=.008/.014 the source tracks contain 121,513 / 141,867 rows over their
45 terminal nodes each. Cumulative trapezoidal surface-abundance integration
against lost stellar mass differs from the published integrated tracked
wind yield by at most 0.1493% / 0.1148% of initial stellar mass. This is a
source sampling/aggregation difference, not a floating-point claim.
Calibrating each monotone cumulative elemental history to its published
wind endpoint leaves the sum <= cumulative lost mass at every inspected
track point. Such calibration, if used, must be named and quantified.
Rounded AGE values repeat near terminal burning, with changes of order
1e-9 Msun in the example 120 Msun tracks. Combine identical-age cumulative
knots deterministically; do not create zero-duration divisions.

## Energetic boundaries requiring an explicit choice

HW02 https://2sn.org/DATA/HW01/bulk_yields.txt supplies PISN E_expl for
He-core masses 64--133.3 Msun. Woosley2017
https://arxiv.org/abs/1608.08939 Table1 supplies PPISN pulse kinetic energy
at He masses 34--62 Msun. **Its 64 Msun row is a complete PISN**, not a
measured PPISN point. A PPISN energy closure on 62--64 cannot silently be
called same-fate table interpolation. The planned feedback comparison must
specify its endpoint treatment; no tabulated wind velocity was found.

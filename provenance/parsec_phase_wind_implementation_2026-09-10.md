# PARSEC phase-dependent wind energy — bounded PASS

Driver end evaluation: implemented, physical inputs converted, Intel/GNU
native tests and MPI2/OMP2 continuous/exact-restart run PASS. One Fable
[plan review](parsec_phase_wind_fable_2026-09-10.md), no extra end audit.
Group4 phase-speed connection is delivered as a selectable comparison;
radioactive inventories remain unfinished. No new approval wait or default
change. Workdir `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.

## Physical prescription and runtime connection

Existing `build_parsec_pair_feedback.py` gains `--wind-model phase_escape_f22_v1`.
It rejects `--wind-km-s` for that named model. `fixed` remains the default
and still requires an explicit speed. Native v4 cumulative wind energy is
consumed through the existing stellar source/thermal feedback bridge; no
Python in RAMSES, new carrier, second radial kick, CR-wind acceleration or
wind dust-shock claim. No RAMSES namelist semantics changed, so its two
generators need no new key. Makefile VPATH was untouched.

For hot tracks use the electron-scattering-corrected escape-speed law of
[Fichtner2022 equations1--3](https://arxiv.org/html/2201.07244v2), with
nonrotating factor1 and (Z/.02)^.13. Our simplified branch classifier is
explicit: Teff<=10000K has10km/s; hotter X_H<.4 has d1.6; other hot states
d2.6. This is NOT the paper's HRD/optical-depth classifier. Cool A/yellow
supergiants are lumped with slow red outflows; bi-stability, eruptive LBV
and dense atmospheres are not resolved. Extending this scaling above the
author's158Msun grid to600Msun is a comparison assumption. No clipping of
Gamma_e or arbitrary velocity floor is used.

The actual90PARSEC nonrotating14--600Msun/Z.008/.014 tracks determine M,L,
Teff and surface H. R is derived from L/Teff and checked against RSTAR.
Gamma_e is NOT the track's EOS-Gamma column. Integrate endpoint-average
specific kinetic energies against the SAME rescaled mass increments as the
positive elemental-margin balancing. Include cumulative energy in the
existing1e-4 compression criterion. The constant-speed branch retains its
exact analytic integral and old knots; both original fixed output hashes
are byte-identical. No duplicated mass/element conversion pipeline.

Mass/element endpoints, five fates, baryonic remnants and all explosion
energies remain EXACTLY unchanged. Energy-sensitive compression adds knots;
intermediate interpolated material histories can differ within the same
declared1e-4 per-final-component approximation. This is not an exact
pointwise equality claim between the two compressed material curves.

Native SED loader accepts the two exact fixed/phase model IDs and still
requires equality with the actual prepared feedback table. All consumed
feedback energies and SED moments remain MPI/restart-bound. No implicit
restart migration. The matched high-mass source still is NOT a full SSP.

## Luminosity-unit correction discovered by the radius check

The initial candidate was rejected before writing any input. Using nominal
IAU2015 Lsun3.828e33 produced a systematic .2340588% radius mismatch. Actual
RSTAR/Teff/logL imply3.84598265e33 with~1.3e-9 relative row spread, consistent
with [PARSEC Bressan2012 Table3](https://arxiv.org/pdf/1208.4498) Lsun3.846e33
and its older Stefan constant. BOTH wind and PARSEC SED converters now use
the source's unit. Maximum remaining R discrepancy2.2559997e-6 is below the
tighter1e-5 check; the bound was not loosened to admit the wrong unit.

The earlier SED evidence used3.828e33 and remains historical. Measured five
Q constraints are unchanged; bolometric normalization, sub-Werner residual
and full-track Planck tails are corrected. Corrected fixed and phase SED
files are byte-identical EXCEPT their model line. Do not assert equality
to the older nominal-unit payload. BPASS and other stellar converters are
unaffected. The adopted luminosity unit is now in both manifests.

## Reproduction and compact evidence

Root `.parsec-wind.Soq8Rf/`, source archives remain in `.parsec-sources.9ZIUlJ/`.
Generate with existing converters:

```
python3 simulation/snrt/tools/build_parsec_pair_feedback.py \
  --source-dir .parsec-sources.9ZIUlJ --output NEW_FEEDBACK \
  --wind-model phase_escape_f22_v1
python3 simulation/snrt/tools/build_parsec_native_sed.py \
  --source-dir .parsec-sources.9ZIUlJ --feedback NEW_FEEDBACK \
  --output NEW_SED --escape-fraction 1
```

Set `PHASE0_YIELD_TABLE=NEW_FEEDBACK/yields.dat`, existing
`high_mass_history_path=NEW_FEEDBACK/history.nml`, and
`SNRT_STELLAR_SED=NEW_SED/source.nml`, with absolute paths and existing
PARSEC v4 channel/IMF admission. Outputs refuse overwrite.

- Feedback90nodes,26641rows including26281wind knots; maximum compression
  error9.999989633225664e-5, energy included. Speed10--3346.520km/s,
  maximumGamma_e.9895600555. All three closure branches occur.
- Existing `parsec_pair_history_test.f90` passes Intel and GNU bounds/FPE
  checks, now additionally comparing every node's material/remnant and
  terminal energy with the fixed baseline. All five fates, time/Z
  telescoping and IMF64/128 remain tested. Existing radiation smoke passes
  all four IMFs and six actual feedback mismatch rejections.
- Retained `check_sources.py/check-sources.log`: elementary law/admission,
  constant-speed analytic limit, original fixed bytes, corrected SED
  equality, independent90-node energy integration. Maximum relative energy
  discrepancy7.327471962526033e-15. No new test framework.
- For1e4Msun Chabrier SSP, Z.011,20Myr: wind returned994.116338113309Msun,
  total returned1537.98680082280, remnants670.823145723004Msun unchanged.
  Wind energy2.156563200597274e52erg versus fixed9.883852574240858e51erg.
  CCSN2.861600964241552e52 and pair3.258050887085736e52erg unchanged.
- Corrected SED90nodes30440knots, native identity670619doubles. Relative
  bolometric closure1.83449e-14, maximumQ-constrainedE/L.97922622. Chabrier
  1e4Msun lifetime Q1.155645565286312e66 and E7.071216795408662e66eV.

## Live run and driver verdict

Binary `ramses_parsec_wind3d`, SHA256
`3f62d81e5fb3257343c56d0943efdf193bf38161791ab7bf078ab3daa49bc127`.
IntelMPI/ifx checked CPU/OpenMP, NVAR19/NENER1/NVECTOR32,
SNRT1/DUST_LIVE0/HDF51,USE_FFTW0,FDMDEBUG1; incremental build reused the
same-option `.parsec-sed.ohf9MX` objects, executable copied into this root.
Old executable retained. New native-fixture edits do not alter binary physics.

Effective `live/physical.nml` and `restart/physical.nml` retained, with
`environment.sh`. Noncosmological periodic4^3,4steps,MPI2/OMP2, gravity/SF,
PARSEC phase winds/CCSN/pair thermal feedback, advectiveCR SNfraction.1,
matched correctedQ/E and HHe node-FS2010. Dust/CHIMES/AGN/sinks/ordinary
cooling and CR-SF support intentionally off. NENERCR included in total gas
energy. Initial density.001,H.74,He.25,Fe.01,thermal pressure1e-8,CR1e-8;
physical units and all hydro/SF options are in the retained namelist.

noutput1,aout2,tout1e30,foutput2,fbackup1000000: two fresh snapshots and
one restart final. No unintended scheduled outputs. Free space58TB atlaunch.
Wall49.640s fresh,40.971s restart2->4. All82physical datasets EXACT,
including38hydro,17particle,3SNRT datasets and AMR/gravity/etc.; all numeric
fields finite, all SNRT attributes and physical clocks exact. Total
gas+star mass.001exact; minrho.0009995350152522352, minthermal6.074272338414694e-7,
meanCR4.557530053623415e-8. PrimaryN/E code-cell sums3.5465533128269 /
20.70626361499916 (not volume-integrated physical photons/erg).

`evaluation.txt` records the measurements and checks. No MG nonconvergence,
Fortran/MPI abort, or source failure. Pre-existing noncosmological SFRD NaN
is a diagnostic, not a physical field. This is NOT a closed isolated-energy
box or galaxy calibration. Bounded implementation PASS; no universal
atmospheric/phase accuracy or radioactive-decay completion claim.
Evaluated raw outputs are removed under the operator's retention rule;
see [cleanup manifest](parsec_phase_wind_raw_cleanup_2026-09-10.md).

## Physical input identities

Feedback yields SHA256 `4135dbcfb73f19f71f552203eb688c5cca98694dac7de57eb1dd3a43e1770b5c`;
history `49cf5470abe2ea80b3b2a727757fbbd8d363a3b0c30f49691dd3848de71fa7f9`.
SED nodes `8a477e9de2d3360f34fc213ff485735e29eb9e47ad75f7b69ce15e851c29dcfa`;
SED nml `e46c72f3d373343ab2ead0d6e8d49b336f484266dda2ddf8c8392100d7f45bf1`.

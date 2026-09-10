# PARSEC common high-mass metallicity extension

Status: bounded source/native/MPI2 live+exact-restart PASS; driver evaluation complete.
Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
This extends an existing comparison population, not all medium groups or a
full SSP / low-Z population / publication calibration.

## One plan review and disposition

[Fable](parsec_metallicity_extension_fable_2026-09-10.txt) returned PROCEED
WITH CORRECTIONS, Q-GOAL yes and Q-LEAN reduce duplicate admission. It
identified the current Z>.014 rejection as a real enriched-gas limitation.
Adopted: source-grid authority in the offline converter; structural native
validation and exact binding to prepared feedback, not a second enumerated
Fortran grid allowlist; shared material/radiation bracket tolerance; old90
source blocks exact subset of new225; one phase-wind MPI2 comparison, not
separate fixed/phase galaxies. Seven changed-source admissions and one
actual five-Z/two-Z mismatch use the existing radiation fixture.

Prechecks requested by the review are complete: all45 photon/track nodes
exist at each new Z and photon end ages precede track ends. PPISN He cores
are34.4338--56.6574Msun atZ=.017 and34.7683--34.9722 atZ=.02, within the
existing W17 domain; Z=.03 has no pair events. No new energy bridge or cap.
The phase-wind model passes all135 new tracks with maximum electron Gamma
.9827024/.9340330/.8012443 at the respective Z. No Eddington clipping.

## Actual source/grid and implementation

[Costa2025](https://arxiv.org/html/2501.12917v2) and the
[author database](https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/)
supply the physical tracks, ejecta and five photon constraints.
New `build_parsec_pair_feedback.py --metallicity-grid metal_rich_five`
selects .008/.014/.017/.02/.03, with original Y=.263/.273/.279/.284/.302.
`solar_pair` is the unchanged default. Source mass is14--600Msun and the
native IMF denominator remains .08--600Msun. Each selected Z has45 actual
mass nodes, yielding225 nodes and66216 native rows (65316 wind knots).

All new wind/terminal eleven-element budgets fit the published baryonic
return without renormalization or remnant edits. Other Z branches with
negative residuals were not inserted as zero, clamped or declared solved.
Z=.01 is still not a source node; .008--.014 is the existing same-age
cumulative mixture. The explicit list in the manifest is not just a range.

`build_parsec_native_sed.py` infers exactly this grid from the actual
feedback manifest and binds source/feedback hashes. Existing Q5/Planck
within-band closure, PARSEC Lsun and full-track Planck missing tails remain.
New photon tail fractions are6.97e-6--1.41e-4 of source lifetime, not an
assumption that the measured photon file reaches death. Final SED84883knots,
max bolometric closure1.83449e-14 and constrained-energy/Lbol .9792262.

Native `snrt_parsec_source` reads structurally valid2--512nodes, sorted Z/M,
at least two Z values, finite nonnegative moments and in-band N/E. It then
requires exact M/Z/fate/age/model equality to the prepared feedback source
before publishing photons. No arbitrary table becomes physically approved
because a structural load succeeds. The existing named physics models and
reference-control status checks remain. Source identity holds1869175doubles
(about15MB); it is startup/dump metadata, not extra per-cell material.

`source_metallicity_bracket` is now shared by high-mass feedback and stellar
RT: same32epsilon relative-to-Z exact-node snap, otherwise linear-Z weights,
no physical extrapolation. This repairs their prior slight roundoff-rule
difference; it does not change interior two-Z physical data. The direct
Fortran module dependency was added to Makefile; VPATH order is unchanged.
No RAMSES NML control/carrier/layout was added, so neither namelist generator
needs a new field. Select existing yield/history/SED file paths together.

## Source acquisition

Retained originals in `.parsec-sources.9ZIUlJ/`, not disposable raw outputs.
Public unauthenticated track downloads again needed per-command
`curl --insecure` due the author's TLS chain. No credentials or global TLS
change. The installed curl lacks --no-clobber; that first command failed
before a download. Successful commands first checked destination absence.
Hashes pin acquired bytes, not independent author authentication.

| Added file | SHA256 |
| --- | --- |
| Z0.017_Y0.279_tracks.zip | a8617cb04eeff29eec752d2fbc3cd2dd993e2b29bf19a87160f034b6d7094d0b |
| Z0.02_Y0.284_tracks.zip | e6a53fdd54cb7690210affcbf3127d0def54d540608c9d75c5a8e87cc812d1c5 |
| Z0.03_Y0.302_tracks.zip | d974343f929d30839cb454b1d8a50ec1d15ab45d0f0b0a18284a51f107be5c28 |
| Z0.017_Y0.279_photons.zip | e38befb640ab2f78e75419b8acd580ee226b91be1740bab04e43388110bc51cc |
| Z0.02_Y0.284_photons.zip | c0e7171fba4de5b6a3b7bccdfe2b31ced926021716e9c59d2088f47e5e6db7e8 |
| Z0.03_Y0.302_photons.zip | 4e99bb71e0043e3b7cb5698c570ed8d26af238ea889cc631f4a18085fe93676c |

Only the three exact photon ZIP members were extracted, without overwrite,
from retained all_photons.zip SHA256
`df31fd6b9fe49abb776ba1701fcc1aef19a407c3ab95768c58537f11b4175f53`.
Old two-Z per-file source bytes stay pinned and are still used directly.

## Completed native evidence

Retained work root `.parsec-metal.TNrmgs/`:

- `superset.log` / `check_superset.py`: every old90 material row and photon
  node block is an exact byte subset of new225. Separately rebuilt default
  `two-regression/` yields/history compare byte-identically to prior two-Z
  phase outputs. No redundant fixed-five conversion was needed.
- `native_tests.sh`, `gnu/`: current-source GNU bounds/FPE builds;225 real
  material nodes exercise66CCSN/80FSN/48PPISN/24PISN/7DBH, correct five-channel
  ownership, same-age elemental/energy mixtures, SSP/time splits, old two-Z
  source and LC18/AGB/effective-Ia shared-bracket regressions all pass.
- Four IMFs, new Z endpoints and all three added interior brackets, upper
  bounds, near-node tolerance, changed actual source and two-Z/five-Z mismatch
  pass the existing native radiation fixture. Intel build also passes,
  `intel-radiation.log`. Harmless underflow/denormal flags occur in an extreme
  Planck tail; invalid/zero/overflow are trapped by the GNU tests.
- Initial `bin` build failed at link because a pre-existing CUDA object
  cache mixed with requested CPU options. No simulation used that build.
  Rebuilt in the known same-option `.parsec-sed.ohf9MX` cache and copied the
  resulting binary into this work root; old binaries remain unchanged.
  Build log `build-matched.log`: Intel MPI/ifx/OpenMP, NVAR19/NENER1,
  NVECTOR32, SNRT/HDF5, DUST0/CHIMES0/USE_CUDA0/USE_FFTW0, bounds checked.
  Binary SHA256 `febc8136908cc08b1242402230d4001b3c394176dbde17267313decc0304ac58`.
  A later source comment and explicit Makefile dependency do not change its
  compiled behavior. No new generic source/binary audit campaign was added.

## Integration and driver end evaluation: complete

Effective absolute inputs:
`/gpfs/kjhan/LRD_JWST/.parsec-metal.TNrmgs/live/physical.nml` and
`/gpfs/kjhan/LRD_JWST/.parsec-metal.TNrmgs/restart/physical.nml`.
Noncosmological4³, MPI2/OMP2, four steps and restart from step2; HHe node-FS
RT, actual stellar phase winds/CCSN/pair feedback and thermal+advective-CR
source fraction .1. Gravity and star formation on; ordinary cooling,
dust, CHIMES, AGN/sinks and CR-SF pressure support deliberately off.
Density .001, gas/CR pressure1e-8, birthZ=.025, H=.682/He=.293 (source linear
Y mixture); initially metals are placed in Fe as a coupling-test IC, not
a realistic full abundance pattern or a calibrated galaxy.
Noutput1, aout2, tout1e30, foutput2, fbackup1e6; two fresh dumps and one
restart dump expected, estimated<=40MB each/120MB total, free102TB checked.
Keep the step2 checkpoint until the restart and driver comparison complete.

Completed fresh four steps in199.006s and restart2->4 in189.154s. Most wall
time is startup/source validation; reported evolution timers are19.806s and
10.070s. This is not a GPU/CPU performance comparison or a new tuning gate.
All82 physical datasets are exact across uninterrupted/restarted runs;
all numeric fields, full SNRT identity attributes and physical clocks pass.
Actual stellar birth Z is .025 throughout this fixture. Gas+star total mass
.0010000000000000002, minimum density .0009995530895406085, minimum thermal
energy density5.780221289567662e-7 and mean CR3.874730479903429e-8.
Code-cell summed photon N3.0661351779180204 and E16.456604738547853 are finite
and group-supported. No MG nonconvergence, source-staging or RT failure.
The old noncosmological SFRD diagnostic NaN is not a nonfinite evolved field.
This is not a closed isolated-energy box or a calibrated galaxy/IMF model.

Artifacts: input `yields.dat` SHA256
`be3c0347cc4932f13e64ad88881d1e0ba0b0f16e890b4fd42c53c520450a30b6`,
`history.nml` `2b470bac05517c5a0c0ba18b761abaf498947b9f923139882ab94c2e599b84c1`;
SED `nodes.dat` `5ce078d6d45bd8ec9c04bd42e32bf339d6aaa5b9f0ffaa68fddfee5901a90a52`.
Retain `evaluate.py/evaluation.txt`, effective NML/environment, source inputs,
all logs and binary. The three32,825,104-byte raw files were removed by
[exact cleanup](parsec_metallicity_raw_cleanup_2026-09-10.md); no checkpoint
is still needed for this comparison. No extra external end audit is pending.

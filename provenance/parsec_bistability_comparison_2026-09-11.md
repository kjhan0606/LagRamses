# PARSEC velocity-only bi-stability comparison

Operator authorized named approximation models on 2026-09-11. This implements
only the bundle5 wind-velocity comparison, not full phase-wind/LBV completion.
Workdir `/gpfs/kjhan/LRD_JWST`; verified origin
`git@github.com:kjhan0606/LagRamses.git`. Read shared context and project AGENTS.
Main owns aggregate design review, the native loader and integration. No new
review gate, MPI simulation, commit, raw-output deletion or default change.

## Named law and scientific limits

Offline selector: `--wind-model phase_escape_f22_bistability_v1`.
Exact native model ID (51 characters):
`parsec2025_w17_hw02_phase_escape_f22_bistability_v1`.

Retain the existing nonrotating F22-inspired electron-scattering escape law:

```
Gamma_e = 0.2*(1+X_H)*L/(4*pi*c*G*M)
v_inf = d * sqrt(2*G*M/R*(1-Gamma_e)) * (Z/0.02)^0.13
```

For `Teff <= 10000 K`, keep the existing 10 km/s outflow. For hotter
`X_H < 0.4` states keep `d=1.6`. For hotter `X_H >= 0.4` states use:

| Temperature | d |
|---|---|
| 10000 < T <= 22500 K | 1.3 |
| 22500 < T < 27500 K | 1.3 + 1.3*(T-22500)/5000 |
| T >= 27500 K | 2.6 |

[Fichtner et al. 2022, section2.1, equations1--4](https://arxiv.org/html/2201.07244v2)
anchors the escape-speed expression and the existing OB/helium factors, but
does not supply this cool-H-rich factor or transition rule.
[Vink, de Koter & Lamers 2001, section8, equations24--25 and adjacent text](https://arxiv.org/abs/astro-ph/0101509)
specifies 2.6 on the hot side and 1.3 on the cool side, with a critical
22500--27500 K bracket. Their recipe chooses a side using the jump position,
not our ramp. **Linear interpolation in temperature is the explicitly
authorized comparison closure, not an equation from either paper.** No
Vink mass-loss formula or mass-loss multiplier is applied.

The 1.3 continuation below the paper's 12500 K grid edge to the inherited
10000 K boundary is another explicit comparison extrapolation. The existing
10000 K discontinuity is retained; there is no second bi-stability jump.
The H-rich proxy is not a resolved OB-supergiant, optical-depth or LBV
classification. F22's grid ends near158 Msun; the inherited extension to
600 Msun is not validated by that grid. Neither paper validates this closure
over every mass/luminosity state in the supplied PARSEC population.

Admitted table grids remain `solar_pair` (Z=.008,.014) and
`metal_rich_five` (also .017,.02,.03), nonrotating initial masses14--600 Msun.
The helper rejects Z outside [.008,.03]. The builder rejects
`precision_eleven` with this selector: no new low/extreme-Z wind extension.
The source-coordinate domain is a bounded named comparison, not a claim of
universal wind accuracy. Rotation, altered evolution, eruptive LBV mass-loss
histories, dense-wind atmospheres and resolved wind subphases remain absent.

## Implementation and identity

Only `simulation/snrt/tools/build_parsec_pair_feedback.py`, its new focused
test `simulation/snrt/tests/parsec_bistability_comparison.py`, and this note
were edited. Existing `phase_wind_speed` is unchanged. The new helper reuses
its validation/escape speed and adjusts only the H-rich velocity. Branch
diagnostics have five explicitly ordered entries for the new selector;
old selectors retain their three-entry manifests and existing model IDs.

The existing positive, endpoint-calibrated PARSEC material history and the
existing endpoint-average specific kinetic-energy quadrature are retained:
`DeltaE = MSUN*DeltaM*(v_left^2+v_right^2)/4`. Energy remains part of the
1e-4 per-final-component compression criterion. No new mass/element yield,
remnant, terminal age, fate, explosion energy, pulse timing or mass-loss rate
is created. Changing energy can select different retained material knots,
so equality of compressed material histories at arbitrary times is not
claimed; exact material endpoints and terminal histories are tested.

The manifest records the primary references, chosen ramp, branch order and
all major extrapolations. Both history model_id and model_coordinates change.
The existing SED builder copies the model ID from the feedback manifest;
no SED-builder edit is needed. Radiation Q/E payload remains unchanged apart
from the identifying model line. Native feedback/SED matching must continue
requiring exact identity, with no restart migration or relaxed matching.

## Main-owned native connection

Add only the exact model ID above to the accepted model list in
`patch/lagRamses/snrt_parsec_source.f90`. It fits the existing character(64)
loader field and character(128) feedback identity fields. No new Fortran API,
hydro field, NVAR, Makefile option, runtime namelist key or generator key is
required by this converter change. Existing native wind energy is thermalized
through the source bridge; do not add another radial kick or photon-energy
subtraction for this change.

Use the existing source-interface settings, with matched absolute paths:

```
PHASE0_YIELD_TABLE=COMPARISON/yields.dat
SNRT_STELLAR_SED=MATCHED_SED/source.nml
feedback_mode='channel_resolved'
fate_policy='user_selected_model_v1'
high_mass_preset='source_consistent'
high_mass_history_path='COMPARISON/history.nml'
yield_source_basis='per_star_cumulative'
population_model='single_star_ssp'
imf_id=2, imf_mass_min_msun=0.08, imf_mass_max_msun=600
binary_fraction=0
channel_mass_min_msun=14,1,14,3,14
channel_mass_max_msun=600,8,600,8,600
use_wind=T, use_agb=F, use_snii=T, use_snia=F, use_pisn=T
```

These are the existing bounded population settings, not a complete run
namelist. Existing build support for stellar enrichment and, when consuming
the matched SED, SNRT is required. Native history remains v4,
`wind_linear_terminal_step`, with channels1/3/5. Matched SED remains v4,
`match_feedback_high_mass_only`, `photon_number_and_energy_v1`. Main still
owns loader rebuild and integrated-run preflight; this note authorizes no run.

## Verification and retained evidence

Private directory: `/gpfs/kjhan/LRD_JWST/.parsec-bistability.rzxSDk`.
Inputs are the existing checksum-pinned `.parsec-sources.9ZIUlJ` archives.
All generated inputs and logs are retained; no raw simulation outputs created.

- 90-node actual conversion:26687 rows; all90 wind-energy endpoints change.
  Independent integration maximum relative error8.215650382226158e-15.
- 225-node actual conversion:66468 rows; all225 energy endpoints change.
  Independent integration maximum relative error9.547918011776346e-15.
- Both legacy selectors regenerated: yields.dat, history.nml AND manifest.json
  byte-identical to retained `.parsec-wind.Soq8Rf/input-fixed` and `input-phase`.
- Exact material endpoints and terminal rows; all five diagnostic branches;
  branch mass/energy fraction closure; constant-speed analytic limit;
  real source-time knots; independent cumulative-energy compression check.
- Radius, Gamma_e, H abundance, nonfinite and Z-domain admission tests.
  CLI conflict, unknown selector and extreme-Z rejections create no package.
- Actual matched SED generation:90 nodes,30440 knots, bolometric relative
  error1.8344902269312498e-14. Q/E payload equals retained phase SED except
  the model line; required native interface flags checked.
- Existing retained GNU `source_test` passed on the new90-node physical
  package and regenerated fixed baseline. This is a native reader/source
  calculation, not a newly rebuilt runtime or MPI simulation. Binary SHA256:
  `5a91f57f9860d555dc0b48cd537b2487a5af44e6810b65fe6e6227da31415046`.
  Five fate counts12,23,24,24,7. Reference SSP wind energy
  2.0897716214505891e52 erg; returned1537.9868008227991 Msun,
  remnant670.82314572300379 Msun; CCSN/pair energies unchanged.

Logs: `converter-final-test.log`, `converter-five-final-test.log`,
`native-source-test.log`. Native SED loader admission of the new model remains
main-owned and is not claimed tested by the old retained source executable.

| Payload | SHA256 |
|---|---|
| comparison/yields.dat | 7ab6405a1ba42c788443c327c5c548ec872d2312f194892b1b27c87bdd74d68b |
| comparison/history.nml | 264229aad5390114d4a560ef9e7afb92bf2a67d6cae2f5bd1d6f2ad022d8d5b9 |
| comparison-five/yields.dat | ec52df49f0910374299795b58bf3cebe6539e404c02beb351e0afa77173b0ddb |
| comparison-five/history.nml | 92b272c2d3fcc2df382f293fcec6bf354f054b9d952e28347fc39448396e6859 |
| sed-comparison/nodes.dat | 360a4301a7bfcaebe2f3d5e654a6d2775fd94308d863bc5d6d3eb14506f935b1 |

Reproduce into fresh output directories (converter refuses overwrite):

```sh
python3 -B simulation/snrt/tools/build_parsec_pair_feedback.py \
  --source-dir .parsec-sources.9ZIUlJ --output NEW_COMPARISON \
  --wind-model phase_escape_f22_bistability_v1
python3 -B simulation/snrt/tools/build_parsec_native_sed.py \
  --source-dir .parsec-sources.9ZIUlJ --feedback NEW_COMPARISON \
  --output NEW_SED --escape-fraction 1
python3 -B simulation/snrt/tests/parsec_bistability_comparison.py \
  --source .parsec-sources.9ZIUlJ --package NEW_COMPARISON \
  --baseline .parsec-wind.Soq8Rf/input-phase \
  --sed NEW_SED .parsec-wind.Soq8Rf/sed-phase
```

For the225-node test use `--metallicity-grid metal_rich_five` and a phase
baseline built with that same grid. Additional test options `--legacy-fixed`
and `--legacy-phase` each take regenerated and retained directory pairs.

## Updated native SED/feedback binding smoke — PASS

Main added the exact bi-stability model case; Gibbs subsequently owns the
shared loader for mixed low-mass work. Tested a private snapshot without any
edits to Gibbs's files or to `.pah-catalytic.Av4xCW` (including its executable).
The model ID length was checked from the actual string: **51 characters**.

Evidence directory: `.parsec-bistability-native.ljpuEY/`.
Copied eight existing Intel stellar/dust support objects and their module
files plus amr_parameters.mod from `.pah-catalytic.Av4xCW`; byte comparison
after testing confirms these originals still match the copies. Recompiled
private copies of the current spectral contract, updated PARSEC loader,
stellar SED frontend and the unchanged
`patch/lagRamses/snrt_parsec_source_smoke.f90` using Intel ifx2025.3.0,
`-qopenmp -fpp -DWITHOUTMPI -O0 -g -check bounds -traceback -fpe0`.
The reused support objects retain their original compilation flags; the new
debug flags do not retroactively apply to them. No MPI runtime or RAMSES
simulation was launched, and no shared build directory was used as output.

The snapshot contains private physical feedback/SED inputs and the existing
reference-control spectral contract. Only the private SED node-file pathname
was changed to its private copy; photon payload and SHA identity are unchanged.
`build.sh`, `test.sh`, `snapshot-inputs.sha256`, `workspace-sources.sha256`,
`effective-inputs.sha256`, `binary.sha256`, `build.log`, `matched.log` and
`mismatch.log` retain the reproduction commands, configuration and identities.

- Actual bi-stability SED + actual bi-stability feedback: exit0,
  `PARSEC_COMMON_NATIVE_PASS IDENTITY_VALUES=670619`.
- Native IMFs0/1/2/4, same-age metallicity mixtures, Q/E interval telescoping,
  spectral support, pre-birth/dead-source limits, no metallicity extrapolation,
  mandatory energy output and deterministic identity all pass.
- All seven existing failed-rebind cases pass: fate, death age, source mass,
  history version, channel mass domain, model identity and source metallicity.
  Every failed rebind leaves Q/E publication disabled; valid rebinding works.
- Separate process with new SED and legacy phase-wind feedback: exit0 in
  rejection mode, no Q/E publication. The fixture's pre-existing marker is
  `PARSEC_ACTUAL_GRID_MISMATCH_PASS`; in this test the two grids are identical,
  so the exercised mismatch is the wind model identity, not the grid.

Snapshot loader SHA256:
`3730b18d9b9396f833950128f5d7ee1582609b99b98c4edc45b86e3e9d1c6a1f`.
Unmodified smoke source SHA256:
`98c3bb9f05ef212c622b06fbfa4fadc6a3867660b9b407257e1c96a54de4dbba`.
Private native smoke executable SHA256:
`c2a4822c7baefe8e5ea225988187f4e2c5e0f4f4de938ea548ac18865276d227`.

This closes the requested velocity-only native binding check. It does not
claim verification of subsequent Gibbs mixed-population changes. Main owns
the aggregate design review; no additional audit or simulation gate added.

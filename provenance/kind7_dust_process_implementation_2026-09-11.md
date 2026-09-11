# Kind7 dust processes: implementation and driver evaluation

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
**PASS for approved bundle1's bounded native connection.** This is not
completion of the remaining PAH/Fe, full-SSP, microscopic-Ia or isotope tasks.
[Plan and selective Fable disposition](kind7_dust_process_bundle_plan_2026-09-11.md).

## Delivered coupling

- Kind7 admits existing C/silicate source condensation, unresolved SN
  destruction, sublimation and relative dust dynamics. Kind5/6 retain their
  old restrictions. Fe/PAH remain excluded. The narrow existing-sink profile
  admits condensation only, not shocks, sublimation or relative dynamics.
- Competitive photo absorption accumulates actual node-weighted four-phase
  energies and signed directional moments on the accepted CVODE trajectory.
  Positive/negative moment counters are separate inside the positive solve.
  The old API and checkpoint carrier widths are unchanged.
- Absorption gives each phase momentum E*n/c using physical c. Midpoint
  mechanical work is debited once from its material heat. Chemistry uses the
  staged row's phase-aware kinetic energy, including existing nonthermal
  terms, and preserves staged mechanical energy.
- Stationary primary scattering is disabled for kind7 relative grains.
  Current per-ray node spectra supply photon-weighted group-grey transport
  opacities to the existing nine-group conservative moving-scatter solver.
  This is not node-resolved moving scattering or cross-group Doppler transfer.
- Existing mass/source, gas-element reconciliation and material/IR operators
  retain their ownership. Opacity is frozen during photo and updated after
  material evolution for the next call. Split sublimation is permitted with
  relative motion; coupled-IR sublimation remains coadvected-only.
- Both setup frontends were updated together. Make VPATH is unchanged;
  the CHIMES runtime has one explicit moving-scatter module dependency.

## Integration defects fixed without relaxing budgets

Absolute moving IR must not inherit the net-bath lower-energy rejection.
Its C/silicate callback now uses the already implemented analytic DL01 cold
material and per-band Planck continuation. Phase enthalpy reconstruction,
temperature extraction and mass exchange use that same mixture energy.
Below the sublimation table, unchanged mass is accepted only if evaluating
erosion at the warmer first knot gives exactly zero represented mass loss.
No bath heating, reaction-table extrapolation or finite sublimation-rate
extension is invented. CHIMES's separate grain-temperature bounds remain.

An unchanged cgs mass cannot be converted back through ordinary arithmetic
and then treated as a strictly one-way loss. Under production ifx-O3, the
isolated40000-case reproduction gave1250 apparent mass increases with the
survival-ratio expression, zero with an explicit equality/identity branch.
That branch now precedes rescaling nonzero erosion. Strict phase no-growth
and energy acceptance remain; failure logs identify the rejected operation.

## Measured evidence

Private work/effective inputs/logs: `.kind7-dust.PakmEZ/`.

| Evidence | Result |
| --- | --- |
| Native thermochemistry |264 PASS assertions, including node N/E reconstruction, exact per-phase capture/impulse, signed anisotropy, zero-step and rollback |
| Native dust mass/backend | Existing suites PASS, including source/shock/latent/phase dynamics and cold absolute-IR two-step energy/rollback regression |
| Frontend CLI/GUI |49 pass,1 unavailable-display skip |
| Coadvected actual-source run | MPI2/OMP2,2 steps; max element8.471e-16, charge4.270e-16, IR1.662e-12; gas+star mass residual zero in FP64 calculation |
| Relative actual-source run | MPI2/OMP2,2 steps; max element6.353e-16, charge4.863e-16, mass2.169e-16, IR6.359e-10; nonzero phase drift and stellar mass return |
| Restart |442 datasets at each of steps1/2 and all compared physics/clock attributes bitwise equal |
| Formatting | `git diff --check` PASS |

Final relative executable `ramses_kind7_identity3d` SHA256:
`b3395e1d5585f2bf1de887902d30e2c72295144b874104754ee0d695ec4398c3`.
NVAR199, NENER0, DUST_DYNAMICS, CHIMES/HDF5, CPU/OpenMP; no VPATH change.
The restart seed was made by the preceding corrected cold-material build.
Its442 datasets equal the final build's fresh step1 exactly; step2 also
matches, so the comparison does not conceal a different initial state.

The actual population is Chabrier .08--600Msun with matched PARSEC14--600
source coverage, not a full SSP. Native source re-evaluation at Z=.02 and
age0--5.148513Myr confirms wind condensation and nonzero **PPISN** energy
entering the SNII+pair shock input. Ordinary SNII is zero in that interval;
no ordinary-SNII dust signal is claimed. Source availability is not a
separate measured destruction-rate calibration. The relative profile has
different split physics/timesteps from the coadvected profile; their dust
losses are not an accuracy/parity comparison.

Reports retained: `coadvected-evaluation.json`, `relative-evaluation.json`,
`restart-comparison.json`, `source_check.log`, native/build/frontend logs,
effective namelists/environment, executables and text snapshot metadata.
After evaluation/restart completed,9 raw dump directories (including
superseded failed-run prefixes), about566MiB, were deleted. `cleanup.json`
and raw hashes record exact targets; no active checkpoint remains required.

No GPU, cosmological, general AMR-subcycling convergence, calibrated galaxy
or implicit sublimation/photo-opacity claim follows from these bounded runs.
No routine external end audit was added; the driver performed this evaluation.

# DUST-3 implementation evidence

Workspace / project: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses.
Base commit: 80c680f. Plan approved conditionally by Fable; operator approved
execution before plan review. Status: implemented static gray candidate;
Opus bundle-end verdict CONDITIONAL PASS. No native/live/production promotion.

## Implementation

- `snrt_core/dust_ir.py`: stable CMB differential energy/photon emission;
  shared-interpolant conservation; explicit six-face vacuum escape;
  energy-density transport via the existing S_N kernel; bounded local
  reabsorption iteration; local and global convergence rejection.
- `tools/p6_run_dust_ir.py`: validated P5 frozen-primary source, exact static
  input and opacity/thermal hash binding, same-table Planck opacity,
  exclusive output creation and compact balance/iteration diagnostics.
- P5 adds `static_input_sha256`; its numeric physics and the DUST-2 core
  and sidecars are unchanged. No new opacity sidecar schema was introduced.
- `P7_DUST_IR_TRANSPORT.md` describes units, approximations and invocation.

## Tests performed

CPU JAX 0.11.1 float64 in `simulation/snrt/.venv`:

`JAX_PLATFORMS=cpu .venv/bin/python tests/dust_ir_transport.py`

Passed: analytic two-IR-band emission, nonnegative energy/photon differences,
13.1 K background with 1e-42 erg/cm3/s primary heating, exact zero sentinel,
tracked+outside closure (1e-10 relative), transparent and anisotropic six-face
escape balance, reabsorption, zero dust/source, and explicit CFL,
nonconvergence and temperature-range rejection.

Uniform synthetic cube, unit box/light speed, primary heating 1e-24 per
volume/time, duration 2, two bands each with absorption coefficient 0.8:

| Run | Stored energy | Boundary escape | Outside escape | Balance relative | Stationarity |
|---|---:|---:|---:|---:|---:|
| 4^3, S4, CFL .4 | 2.93512426e-25 | 8.51312578e-25 | 8.55174995e-25 | 6.52e-10 | 5.52e-4 |
| 4^3, S4, CFL .2 | 2.91841925e-25 | 8.50706915e-25 | 8.57451158e-25 | 6.74e-10 | 3.47e-4 |
| 8^3, S4, CFL .4 | 2.82283546e-25 | 8.62669374e-25 | 8.55047079e-25 | 6.26e-10 | 1.58e-4 |
| 4^3, S8, CFL .4 | 2.94067510e-25 | 8.49759320e-25 | 8.56173169e-25 | 6.52e-10 | 5.28e-4 |

Primary injection = 2e-24; all local solves required at most 29 iterations.
Stored-energy differences vs baseline: dt half -0.57%, mesh double -3.83%,
S8 +0.19%. Mesh doubling also halves the CFL timestep; this is combined
space/time refinement, not a measured spatial order. Stationarity measures
field change per step relative to max(injected step energy, old inventory).
These are bounded comparisons, not universal accuracy tolerances.

`JAX_PLATFORMS=cpu .venv/bin/python tests/p5_dust_runner.py`

Passed existing DUST-2 controls plus actual P5-to-P6 invocation at both opacity
temperatures, output assertions, and rejection of a valid but differently
hashed static snapshot before creating an output.

| Planck reference | Band sigma_abs/H cm2 | Integrated reprocessed energy erg | Outside fraction of emitted | Balance relative |
|---|---:|---:|---:|---:|
| 20 K | 6.34602197e-25 | 4.481458e32 | 0.664065 | 6.590e-11 |
| 50 K | 1.83550573e-24 | 1.296204e33 | 0.664065 | 4.767e-11 |

Physical experiment: same 2^3 static input and frozen final P5 heating;
duration 1e10 s, c_hat/c=.01. The approximately 2.89x opacity/reprocessing
spread demonstrates linear optically thin propagation of the opacity change,
not nonlinear trapping sensitivity. The synthetic cube exercises the latter.
About 66.4% of
emission leaves through the spectral complement, principally below .01 eV.
This is an optically thin wiring study, not a science simulation.

## Input assets

Raw Draine SHA256: b56680cc38b85f051f20c4405303e8c480cc9bec714fd5ba722a257a40ae840c.
Opacity SHA256: 7521ef988a47b590f375f49cdedf375109f5ee306968e54749b38f5e43a1faa8.
Thermal sidecar: e0065f2b6de47b43f1b2739ff7fedb4ba90f074d17f775955672e00f42da5259.
Generated test inputs/outputs are temporary; commands reproduce them. Study
outputs record actual static/P5/sidecar/raw-table and implementation hashes.

## Post-audit bounded repair verification

Opus verdict: CONDITIONAL PASS; report and driver dispositions are in
`opus5_dust_ir_transport_end_audit_2026-09-06.md`. Added distinct thermal versus
iteration failure reasons, strict total-power admission, four guard tests,
dataset units and limitation labels. The two original suites were rerun:

```text
DUST_IR_TRANSPORT_TEST_OK two_IR_groups=1 weak_CMB=1 failures=7
P6_PHYSICAL_IR Tref=20 stationarity=0.0102481 tau_cell=2.88348e-07 balance=6.590e-11
P6_PHYSICAL_IR Tref=50 stationarity=0.0102481 tau_cell=8.34011e-07 balance=4.767e-11
P5_DUST_RUNNER_TEST_OK groups=9 dust_model=metadata
```

The tabulated energy results above remained unchanged. At the finite end
time the physical field still changes (stationarity 0.01025); no strict
steady-state result is claimed. The initial fixed source is stationary, but
the radiation field evolves toward a steady state.

## Remaining promotion scope

Fixed band opacity, finite spectral coverage/free outside escape, omitted IR
scattering, frozen primary heating, single-temperature dust. No dust force,
gas exchange, native/live/MPI/restart qualification or spectral convergence.
Spectral-complement free escape underestimates trapping in this fixed-opacity
model. Long-wavelength coverage and temperature/frequency-consistent opacity
are prerequisites for science claims or native IR promotion.

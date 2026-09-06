# DUST-4 spectral IR implementation evidence

Project `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, base `5ab016b`.
Status: implemented static spectral candidate; Opus CONDITIONAL PASS,
all bounded conditions implemented and both complete suites rerun successfully.
Operator has preapproved continued work absent material blockers.

## Implementation and scope

Extended `snrt_core/dust_ir.py` with the spectral constructor, and existing
P6 runner/two dust test files. No new builder, schema, validator or test file.
Full pinned Draine-domain quadrature (4-point Gauss in log-energy bins),
node-identical emission/absorption opacity, exact discrete total-power curve,
inserted CMB bath node, stable weak-excess source, per-node photon energies.
Reuses the DUST-3 conservative transport/reprocessing step unchanged.

P6 spectral mode is explicit, rejects a gray reference-temperature flag,
retains primary nine-group input binding and exclusive output creation, and
records spectral nodes/weights/opacity/T/domain, code and input hashes.
Frequency-summed cubes and one-dimensional spectral inventory avoid storing
directional fields. Gray results remain controls. Plan review conditions and
driver dispositions: `fable_dust_spectral_ir_plan_audit_2026-09-06.md`.

## Executed compact tests

CPU JAX 0.11.1 float64, existing `.venv`, no simulation jobs submitted.

`JAX_PLATFORMS=cpu .venv/bin/python tests/dust_ir_transport.py`

Passed the existing gray controls and spectral extensions: wider-domain
constant-opacity Stefan-Boltzmann normalization (1e-6 tolerance), node-wise
full-c Kirchhoff identity (1e-12), bath insertion/deduplication, positive
1e-42 weak-excess energy/photon closure, invalid opacity/CMB rejection.
Full suite was executed once; after a finite-power guard was added, only
`spectral_checks()` was rerun, with identical numbers. No new dt/mesh/angular
study was commissioned; historical gray checks remain in the existing suite.

Manufactured dusty cube: 2^3 cells, width 1e18 cm, effective reference-mixture
density nH*dust_relative = 1e6 cm^-3, T_CMB=13.1 K, S4, CFL .4, full c,
duration 4 dx/c. Identical primary heating in both runs is independently
integrated on the dense raw opacity grid from B_E(20 K)-B_E(13.1 K).
Peak box optical depth at 20 K is .759; this is a moderate-depth numerical
experiment with physical opacity, not an approved astrophysical simulation.

| Quantity | 136 nodes (4 bins/decade) | 264 nodes (8 bins/decade) |
|---|---:|---:|
| Stored energy (erg) | 2.0883963775e45 | 2.0883925634e45 |
| Reprocessed energy (erg) | 2.7542085843e45 | 2.7541658118e45 |
| Mean diagnostic T (K) | 21.17696746 | 21.17742326 |
| Stored below .01 eV (erg) | 1.0930388332e45 | 1.0928452780e45 |
| Final below-.01-eV absorption rate (erg/s) | 6.71625924e36 | 6.71000640e36 |
| Global balance relative | 7.40710e-10 | 7.40762e-10 |
| Maximum source iterations | 33 | 33 |

Spectral-complement escape is exactly zero in this finite-domain model.
Stored change = 1.8264e-6 relative; reprocessed change = 1.5530e-5;
max temperature change = 2.1523e-5. These satisfy the predeclared 2% energy
and 1% temperature comparisons; no third frequency resolution was needed.
These are one-case sensitivities, not a universal production error bound.
Reprocessed energy is internal recycling, not another global sink/source.

`JAX_PLATFORMS=cpu .venv/bin/python tests/p5_dust_runner.py`

Passed actual P5->P6 spectral invocation alongside unchanged gray 20/50 K
results, inactive gray-flag rejection, output inventory closure and input
hash-mismatch rejection. The spectral run has 136 ordinates, stores 66.3898%
of its field energy below .01 eV, and closes energy to 5.339e-11. This P5
cube is optically thin: moving its former complement into the radiation
field mainly produces boundary escape, not recovered strong trapping.
After adding the all-temperature-node sidecar comparison, this existing
runner suite was rerun successfully. Maximum relative power-curve difference
is .0120929 (1.20929%), including old log-T interpolation error at the newly
inserted bath node. It does not enter the new total-power inversion.

## Final post-audit verification

Opus report/dispositions: `opus5_dust_spectral_ir_end_audit_2026-09-06.md`.
After all bounded repairs, both complete commands above passed.
Analytic blackbody relative errors: 1.26807e-7 at default 4 bins/decade,
3.24927e-11 at 8. Kirchhoff is a cross-module check of the same Planck law;
the 1e-270 absolute tolerance covers only the old builder's negligible
x=700 Wien floor. This does not supply independent physical normalization;
the analytic blackbody assertion does.

One 136-node half-Courant comparison gives .00575276 relative stored-energy
change (0.5753%) and .00214215 mean-temperature change (0.2142%); below the
declared 2%/1% checks. No spectral mesh/angular matrix was added. Report this
as a bounded timestep sensitivity, not a universal accuracy bound.

P5 spectral output now reports both sidecar differences: .0120929 including
the inserted bath versus .00387798 on original nodes only. The latter
isolates discrete quadrature against the trapezoid thermal sidecar. The
emission-weighted cell tau is 1.56917e-7. Maxima over all frequencies are
explicitly named as such, because nearly unpopulated UV nodes dominate them.
No physical state changed in these reporting repairs; all original spectral
energy/refinement numbers above are unchanged.

## Remaining physical boundaries and handoff

The raw opacity domain is finite: 1.23984198e-4--1.23984198e4 eV.
Conditional low-tail/total-power upper estimate is at most 1.29032e-5 on
the admitted grid, assuming opacity below the data never exceeds the lower
endpoint. It is not a measured bound, an extrapolated input, or automatically
a relative bound on arbitrarily weak CMB-excess emission. The analogous
upper-endpoint Wien-tail estimate is reported in natural-log space.

Source temperature interpolation is only exactly Planckian at thermal
nodes. Mixture remains single-temperature, with fixed dust, no IR scattering,
frozen primary heating, and no force or gas-energy exchange. Reduced-c energy
inventory is not physical full-c LTE density. Native/live AMR, force/gas
coupling and dust evolution remain open. DUST-3's missing *in-table* far-IR
channel and arbitrary reference-temperature absorption are resolved in the
new spectral mode, not retroactively in gray mode and not as a claim of
complete physical dust modeling. Next useful work is native dust operator
integration, with explicit conservative coupling semantics, not more Python
gate infrastructure.

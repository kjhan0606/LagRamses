# DUST-3: bounded IR transport promotion study

Authorization: operator approved the next bundle after DUST-2 (80c680f).
Workspace / project: /gpfs/kjhan/LRD_JWST, kjhan0606/LagRamses.
Final objective: production/publication-ready RT, stellar/AGN feedback and
dust physics in lagRamses. This bundle supplies the missing reabsorption and
transport experiment needed before native IR coupling can be justified.

## Deliverable

A callable static JAX IR energy transport step and a small executable study
runner consuming a validated DUST-2 thermal sidecar, its matching v3 opacity,
the static gas input, and a P5 dust-heating snapshot. The primary heating is
explicitly frozen at the supplied rate; this is a stationary-source numerical
study. It does not qualify concurrent primary RT or live hydro evolution.

1. Correct the CMB differential emission in the new transport closure:
   per-group power and photon number are differences of the tabulated
   emission at T and T_CMB. Interpolate tabulated group power and photon rate
   linearly against total power, and integrate segment slopes across the
   excess-power interval to avoid cancellation for weak heating. Total tracked
   plus out-of-band power closes. Keep DUST-2 historical controls intact.
2. Transport directional **energy density**, not a photon count multiplied
   by a temperature-dependent mean energy. This avoids energy creation when
   photons move between cells with different dust temperatures. The existing
   S_N upwind/analytic-local-absorption kernel is linear and reusable with
   energy units. Export emitted photon counts as a separate source diagnostic.
3. Retain the current nine-group contract and transport only its complete
   IR group(s) below 1 eV. Out-of-band excess escapes in a separate ledger.
   Compute fixed, explicitly temperature-labelled Planck-weighted absorption
   from the same raw Draine table (default reference temperature 20 K).
   Omit IR scattering in this bounded absorption/re-emission study and label
   that assumption. No reuse of AGN-weighted opacity as IR opacity.
4. In each step, transport the old energy once. Iterate the local constant
   re-emission source from primary heating plus the same step's absorbed IR
   energy / dt until convergence (damping 0.5, bounded iteration count).
   Re-evaluate from the old state for every iterate; never repeatedly advance
   time during iteration. Reject nonconvergence, negative/nonfinite state,
   CFL violation, invalid mixture/edge bindings and thermal table overflow.
5. Record energy in the field, primary injected energy, physical outward
   boundary flux, escaped out-of-band energy, reprocessed energy (diagnostic
   only), and the global balance residual. A nonconverged result cannot be
   published as successful. Output creation refuses an existing target.

## Compact checks and evidence

One focused test covers an analytic synthetic two-IR-band thermal model
(including CMB-dominated inputs), a transparent transport energy balance,
a dusty reabsorption case, zero source/dust, and deliberate CFL/iteration/
thermal-range failures. The same homogeneous cube and heating are compared
at dt/2, twice linear resolution, and S4/S8; report differences in stored and
escaped energy without inferring universal production tolerances. Exercise
the runner against a small physical-sidecar snapshot. Reuse the existing
transport kernel; avoid a new gate framework or broad audit matrix.

Plan audit: Fable, necessity/feasibility/overinstrumentation and physical
contract. End audit: Opus 5; Fable is backup only if Opus cannot judge.

## Limits and sources

Fable conditions adopted: explicit six-face outflow and normalized angular
source; local plus global convergence; in-step self-absorption diagnostic;
13.1 K weak-heating test; explicit below-band escape interpretation;
20/50 K opacity comparison; stationarity metric. Add static-input hash to
P5 output so the study can enforce geometry/source snapshot association.
Transparent balance includes stored energy; no finite-crossing-time escape
claim is imposed on a continuously emitting upwind system.

Single-temperature equilibrium; fixed reference-temperature band opacity;
finite photon bands with explicit out-of-band escape; static Cartesian grid;
float64 CPU reference. Grain-size/PAH physics, source obscuration, spectral
convergence, dust-gas exchange, momentum coupling, and native/live/MPI/restart
promotion remain later work. Fixed band opacity is a study approximation,
not a frequency-resolved radiation model.

- Draine raw optical data: https://www.astro.princeton.edu/~draine/dust/dat/mix/
- CMB heating/background distinction: da Cunha et al. 2013,
  https://arxiv.org/abs/1302.0844

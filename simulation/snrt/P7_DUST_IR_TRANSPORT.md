# DUST-3 excess IR energy transport study

Status: static gray candidate. The executable study reads a validated P5
snapshot and freezes its final dust heating rate for the requested duration.
This resolves transport and reabsorption under that fixed source; native
Fortran, evolving primary RT, momentum deposition and live hydro remain open.

`snrt_core/dust_ir.py` transports directional energy density with the existing
S_N upwind/analytic-absorption kernel. The angular source is normalized by the
quadrature weights. Each local source iteration starts from the same old
radiation state, solving equilibrium from primary heating plus reabsorbed IR.
The bounded iteration rejects nonconvergence, invalid states and table
overflow. It checks local fixed-point error and global energy balance.

The exact balance used here (within numerical tolerances) is

`primary injected = stored IR + boundary escape + spectral-complement escape`.

Reprocessed energy counts absorption events and is recycled; it is not another
injected or escaped energy term. Boundary escape is computed independently
from outgoing directions on the six faces using the transport light speed.

The CMB is an external bath. Only excess emission is transported, using
differences of tabulated group power and group photon rate at T and T_CMB.
Segment slopes parameterized by total emitted power retain weak positive
heating even below floating-point precision relative to the bath power.
This defines a new interpolation for the DUST-3 path; DUST-2 remains the
historical one-pass control. Emitted photon counts are source diagnostics;
the transported quantity is energy, so changing grain temperature cannot
reassign energy to photons already in flight.

Reported grain temperature is diagnostic: zero marks inactive cells, and
very weak heating rounds to the CMB temperature although the differential
emission remains positive. Reconstructing a source SED from this temperature
can lose that excess. Use the differential emission rates instead.

The current nine-group contract is unchanged. Only complete groups ending
at or below 1 eV are admitted. Its single [0.01, 1] eV IR band misses much of
the cold far-IR spectrum below 0.01 eV. That complement (plus the negligible
high-energy tail above 1 eV at these temperatures) is recorded as freely
escaping energy. This underestimates IR reabsorption and trapping under the
fixed-opacity model; it must be resolved before science or native promotion.
No spectral convergence or universal dust-mixture approval
is inferred from this study.

Absorption uses `integral(sigma_abs B_E)/integral(B_E)` in each band, from the
same pinned raw Draine table. The reference temperature is fixed for a run,
and 20/50 K results expose its sensitivity. IR scattering is omitted in this
bounded absorption/re-emission experiment. The output records band opacity,
reference temperature, hashes, energy balance, stationarity and the maximum
in-step fraction of a constant source absorbed locally. Large optical depths
may need a different nonlinear solver; exhausting iterations is an error.

The fixed absorption mean and temperature-dependent emissivity do not ensure
frequency-resolved detailed balance away from the reference temperature.
The 20/50 K comparison measures this experiment's sensitivity, not a universal
error bound. Consistent frequency/temperature treatment and long-wavelength
coverage are prerequisites before native IR or science promotion.

```sh
.venv/bin/python tools/p6_run_dust_ir.py \
  --input STATIC.h5 --p5-heating P5.h5 \
  --dust-opacity-metadata OPACITY_V3.json \
  --dust-thermal-metadata THERMAL.json --output NEW_IR.h5 \
  --duration-s 1e10 --opacity-temperature-k 20
```

P5 must record the exact static-input, opacity and thermal-sidecar hashes;
older outputs missing static-input binding are refused. Regenerate the small
P5 control instead of asserting an unverified geometry match. Outputs use
exclusive creation, CPU float64 and fail before successful publication on
convergence errors. The output code hashes identify the implementation.

Sources: [Draine optical data](https://www.astro.princeton.edu/~draine/dust/dat/mix/)
and [da Cunha et al. 2013](https://arxiv.org/abs/1302.0844) for the distinction
between CMB heating and emission relative to the background.

## DUST-4 opt-in spectral extension

Add `--ir-mode spectral` to the invocation and omit `--opacity-temperature-k`
(it is rejected in this mode). `--spectral-bins-per-decade 4` is the default;
8 is the comparison resolution. Gray remains the default historical control.
No primary photon-group edges or AGN photon ledger are changed.

The secondary field has four Gauss-Legendre ordinates per log-energy bin,
with explicit .01 and 1 eV breakpoints. It spans the full pinned raw Draine
table, about 1.23984e-4--1.23984e4 eV. At each ordinate the same log-log
interpolated absorption cross section enters emission and absorption.
Quadrature weights multiply emission/energy, NEVER the absorption
coefficient. Photons are diagnostic energy divided by that ordinate's E.
The [Draine source description](https://www.astro.princeton.edu/~draine/dust/dustmix.html)
documents the finite 1 cm--1 Angstrom coverage and the dust-mass normalization.

The thermal power curve is the sum of these exact discrete emission
channels, not a rescaled old sidecar curve. Its admitted temperature grid
comes from the validated sidecar, with the CMB inserted and deduplicated.
At temperature nodes, including the bath, the closure obeys node-wise
Kirchhoff balance. Between nodes it uses the common-power interpolation:
total energy closes exactly, but the spectrum need not equal a Planck
spectrum at the diagnostic interpolated temperature. Weak positive excess
still uses the stable slope integration, not subtraction of near-equal baths.

All in-domain thermal energy, including below .01 eV, is now transported and
available for reabsorption. The finite-domain model has **zero complement
escape**, not proof of zero emission outside the raw domain. A low-energy
Rayleigh-Jeans tail estimate assumes sigma below the data is bounded by its
lower endpoint; a high-energy integrated Wien estimate assumes the analogous
upper-endpoint bound, recorded in natural-log relative units to avoid numeric
underflow. Both are explicitly conditional estimates, not measured bounds
or extrapolated opacities used by the solver. The low-end estimate is about
1.3e-5 of finite-domain total power at the worst 5 K node of this mixture.
It is a bound on total emitted power, not automatically on arbitrarily weak
CMB-excess power. No full-frequency physical completeness is claimed.

Spectral output keeps the existing format, with `ir_mode=spectral` and
`status=static_spectral_candidate`. `spectral/` contains energies, weights,
cross sections, thermal nodes and domain; the energy-density and emitted
photon cubes are summed over frequency. `energy_per_frequency_erg` retains
the spectrum of total stored energy. `sidecar_power_relative_difference_max`
compares against the old curve (including its log-T interpolation at the
inserted bath). All outputs retain static/P5/sidecar/raw/code hash binding.
`sidecar_original_nodes_power_relative_difference_max` excludes that inserted
node, separating quadrature changes from the old temperature interpolation.
Spectral all-frequency tau/self-absorption maxima are explicitly suffixed
`_all_frequencies`: high-energy nodes can dominate even when carrying no
thermal radiation. `emission_weighted_cell_tau` and
`emission_weighted_in_step_self_absorption_fraction` instead weight by the
cumulative emitted energy per node/cell (zero sentinel for no emission).
Neither is a universal trapping factor. The Kirchhoff test is cross-module
consistency; Stefan-Boltzmann supplies the independent normalization at both
4 and 8 bins/decade. The LTE comparison's tiny absolute tolerance covers the
historical Planck builder's x=700 floor in negligible Wien channels.

For BOTH modes, reduced-c inventory is a conservative transport quantity:
with unchanged physical source power and absorption rate c_hat*kappa*U,
stationary U scales as c/c_hat. Do not interpret it as physical full-c LTE
energy density or use it for force/gas coupling without an explicit RSLA
coupling derivation. The discrete LTE test and moderately opaque spectral
comparison use full c. Spectral refinement here is a bounded case comparison,
not a publication/production error guarantee. The physical P5 cube remains
optically thin; newly covered radiation mainly becomes boundary escape,
not evidence of strong trapping. Single equilibrium-temperature mixture,
frozen primary heating, omitted scattering, fixed dust and no native/live
coupling remain in force. The next handoff is native transport and coupled
energy/momentum semantics, not another infrastructure validation framework.

## DUST-5 native operator

`patch/lagRamses/snrt_dust_ir.f90` implements the secondary spectral thermal
source and conservative transport/reprocessing in FP64 Fortran. It is listed
in the SNRT module build graph, but is **not called by the live RAMSES driver**.
It does not reinterpret the existing nine-group FP32 photon state as energy.
The compiled array API is a native candidate, not live dust activation.

`snrt_dust_ir_initialize` accepts positive frequency nodes (eV), integration
weights (eV), absorption cross sections (cm2/H), strictly increasing thermal
nodes (K) including the exact CMB. It computes the Planck powers itself;
frequency quadrature generation/data admission remain caller responsibilities.
No sidecar loader or implicit dust-abundance model is hidden in this module.

`snrt_dust_ir_advance` takes directions `(3,ndir)`, weights normalized to ONE
(the primary native angular convention is different), reciprocal neighbors
`(6,ncell)` ordered x-/x+/y-/y+/z-/z+, a common cell width (cm), dt (s), c_hat
(cm/s), reference-mixture density nH*relative_dust (cm^-3), and primary heating
(erg/cm3/s). Zero neighbor means vacuum. It is an equal-width cell-set API:
do not pass incomplete MPI ghosts or coarse/fine faces as vacuum boundaries.

State arrays are energy `(nfreq,ndir,ncell)` in erg/cm3 per normalized angular
direction, diagnostic temperature `(ncell)` K and cumulative emitted photons
`(nfreq,ncell)` photons/cm3. Frequencies are energy channels, not primary
ionizing photon groups. Each call returns per-step escaped, absorbed and
primary energy (erg), local/global residuals and iteration count. Reabsorption
is recycled energy, not an additional source. **Only success commits all
state and diagnostics.** Nonzero status leaves all inout arrays and
diagnostics bit-identical; the caller must inspect it before proceeding.
Table/config/shape/state/CFL/thermal-range/nonconvergence statuses are distinct.
The committed step conserves energy to the nonlinear stop tolerance, not
machine precision: 1e-9 relative to injected energy here (observed 7.407e-10);
with zero injection the scale is old radiation inventory, with a tiny floor
for an empty state. Exact closure holds at the converged fixed point.
Future live gas/force coupling must explicitly budget this finite residual.

The algorithm matches the reference bath-relative slope source, log-T
interpolation in total power, thin-cell response, old-state outgoing flux,
zero-primary old-inventory scale, and half-relaxed local iteration. It has no
gas, momentum or primary-photon side effects. The same finite-domain, thermal
interpolation and RSLA interpretation limits apply as in DUST-4.

```sh
JAX_PLATFORMS=cpu .venv/bin/python tests/dust_ir_transport.py --native
```

This opt-in check compiles the actual module and one native smoke driver with
available GNU/Intel compilers in private temporary directories. A plain
numeric fixture supplies the same 136-node manufactured cube as DUST-4;
Fortran computes its own source table and evolves the field. It checks
energy/temperature/photon differences to 1e-8, conservation to 1e-9,
nonconvergence/CFL/neighbor-reciprocity rollback, and zero/weak source.
Default invocation stays compiler-free. Compilation/build-graph inclusion is
not a full production link, MPI/AMR/restart or live-runtime qualification.
Persistent state layout/spectral resolution, dust-density mapping and the
coupled primary/gas/force transaction remain explicit subsequent design work.

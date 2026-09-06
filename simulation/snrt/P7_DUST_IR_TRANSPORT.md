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

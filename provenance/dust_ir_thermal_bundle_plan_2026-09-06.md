# DUST-2 grain thermal balance and IR re-emission bundle plan — 2026-09-06

## Purpose and final-project alignment

This bundle advances the project's final objective: a physically defensible,
production/publication-ready high-level radiation-transfer path coupled to
stellar/AGN feedback and dust in `lagRamses`.  It closes the first missing dust
energy channel after DUST-1, while keeping the current static JAX/SNRT scope
honest.  It does not claim that the chosen dust mixture, depletion law, source
obscuration, or live RAMSES coupling is scientifically approved.

DUST-1 now measures dust absorption and candidate within-group scattering.  It
deliberately leaves the absorbed energy in a separate dust-heating ledger.
DUST-2 turns that ledger into a named thermal/emission contract.  Absorbed
energy must not be relabeled as an IR photon source without a grain thermal
prescription and a frequency distribution.

## Scope

### D2.1 External thermal/emissivity admission contract

Add a versioned `snrt_dust_thermal_v1` sidecar contract with:

1. exact `group_edges_ev` and its hash;
2. a declared list of configured IR output groups, each a complete interval of
   the photon-group contract;
3. a strictly increasing grain-temperature grid in K;
4. reference-mixture emitted power per H nucleus as a function of grain
   temperature, in `erg s^-1 H^-1`;
5. temperature-dependent emitted-energy fractions for every configured IR
   group and temperature-dependent emission-weighted mean photon energy for
   every configured IR group, in eV;
6. a declared reference dust mixture, opacity/emissivity source, grain-size
   and charging assumptions, and the treatment of emission outside the
   configured groups;
7. source/table/builder/code hashes and a payload hash, using the existing
   provenance conventions.  The thermal sidecar also carries the active
   absorption source-table hash and reference dust mass per H nucleus so the
   runner can prove that thermal and opacity data describe the same mixture.

The sidecar is accepted only when the emission fractions and photon energies
are finite and non-negative, each photon energy lies inside its group interval,
and the configured-group fractions plus a tabulated
`untracked_energy_fraction` sum to one within a declared tolerance at every
temperature.  The untracked fraction is an explicit far-IR/out-of-band
ledger, not an unreported loss.  The power curve must be finite, positive, and
strictly increasing; the zero-temperature sentinel is handled separately.
The loader must reject missing group edges, duplicate IR indices, mismatched
edge hashes, malformed temperature grids, an out-of-range photon energy, or an
unproven sidecar provenance chain.  The active opacity sidecar's source-table
hash and reference dust mass per H must match, or the thermal attachment fails
closed.

The sidecar is a physical-input candidate, not an approval.  For the first
equilibrium candidate, the power curve, group fractions, and photon energies
may be derived from the admitted Draine absorption table by Kirchhoff's law
(absorption cross section times the Planck function), under an explicitly
declared single-temperature equilibrium assumption.  The derivation builder,
quadrature tolerance, temperature grid, and convention must be hashed.  A
synthetic table is allowed only for unit and wiring controls and must be
labelled as such; it may not be promoted as a physical substitute.

### D2.2 Fixed-shape local grain thermal balance

Implement a pure static-JAX operator that receives the DUST-1 absorbed dust
power density, the cell dust abundance, and the epoch's CMB background
temperature, then solves

``P_emit(T) × n_H × dust_relative_abundance = P_absorbed + P_CMB_absorbed``

using the admitted sidecar power curve.  The first implementation uses a fixed
count (32-iteration) bracketed bisection in log-temperature space (inverse
interpolation is an allowed simpler implementation).  It must return a finite
temperature, the emitted IR energy density/rate, a per-group IR photon source
rate computed from the temperature-dependent emission-weighted energies, and
an out-of-range mask.

The zero-input branch is exact: zero absorbed power gives zero newly emitted IR
energy and photons, with a temperature sentinel of exactly 0 K rather than an
invented equilibrium.  The CMB term is evaluated from the same admitted
absorption data, so the thermal table must extend to or below the run's
background temperature.  Cold coverage is a load-time error; hot cells return
an out-of-range mask that the runner raises after the step.  No silent edge
clamping is permitted.  The operator must never modify H/He ionization, gas
thermal energy, or the DUST-1 absorption and scattering ledgers.

The operator runs once per existing thermochemical subcycle, using the same
dt-weighting as the DUST-1 heating ledger, and accumulates tracked IR and
untracked energy plus per-group photon ledgers.  The exported temperature is
the final subcycle value.  The first closure is local and one-pass: IR output
is recorded as a re-emission ledger and is not recursively re-transported in
the same step.  This avoids double counting and avoids an unannounced
nonlinear transport iteration.  A later bundle may feed the IR source back
into transport after a matched source/geometry convergence study.

### D2.3 P5 output and runner wiring

Add opt-in `--dust-thermal-metadata PATH` to the static P5 path.  The flag is
required when non-zero dust is present and IR re-emission is requested; absent
the flag, the existing DUST-1 heating-only behavior remains the explicit
control.  A thermal sidecar with an incompatible group contract, opacity
mixture hash/mass, or malformed outside-band fraction fails closed.  The JAX
operator returns masks; only the runner converts out-of-range hits into a
hard error.

Record only the compact outputs needed by the physics closure:

- `thermal/dust_grain_temperature_k`;
- `thermal/dust_ir_reemitted_energy_erg_cm3`;
- `thermal/dust_ir_untracked_energy_erg_cm3`;
- `rates/dust_ir_reemission_erg_cm3_s`;
- `sources/dust_ir_photon_rate_cm3_s` by configured group;
- scalar attributes for thermal schema, IR group indices, sidecar/payload/code
  hashes, energy-closure tolerance, the background temperature used, maximum
  cell IR optical depth, the fraction of dust-bearing cells with optical depth
  above one, and the explicit statement `recorded_not_transport_reemitted`.

The re-emitted tracked plus untracked energy must equal the admitted absorbed
dust energy plus the CMB absorbed term for the one-pass closure, modulo the
declared numerical tolerance.  The closure test must also compare the
table-evaluated emitted power at the solved temperature against the integrated
absorbed power; it must not pass merely by copying the input ledger.  The
existing dust heating ledger remains the pre-emission absorbed-energy
diagnostic; it is not added again to the gas equation and is not counted as an
IR source.

## Compact evidence gates

One bundle-level runner and focused tests are sufficient:

1. sidecar schema, payload/hash, edge binding, monotonic power and fraction
   normalization checks;
2. one-zone equilibrium recovery for an analytic synthetic thermal table;
3. zero absorption, zero dust, out-of-range temperature, malformed fractions,
   and mismatched-edge negative paths;
4. emitted-energy and photon-number closure for one and multiple IR groups;
5. temperature monotonicity with absorbed power and dust abundance;
6. P5 output/attribute checks proving that gas heating, H/He chemistry, and
   DUST-1 absorption/scattering ledgers are unchanged, including the zero-dust
   control and the recorded-not-transported semantics.  Extend the existing
   `tests/p5_dust_runner.py`; do not add a second source-bound subprocess.

Evidence must include tolerances, dtype, group-edge hash, input/sidecar/code
hashes, and an explicit statement that the IR source was recorded but not
recursively transported.

## Explicit non-goals and later work

This bundle does not select a universal dust-to-metal/depletion law, derive a
grain-size distribution from metallicity, add grain charging or stochastic
heating, include dust-gas collisional exchange, model source obscuration, add
anisotropic IR transport, modify native Fortran, couple dust force to live
RAMSES, or qualify MPI/restart/production cosmological runs.  It also does not
invent missing IR data.  If the admitted absorption table cannot support a
traceable Kirchhoff-derived equilibrium table, the deliverable is the
fail-closed contract and a blocked candidate, not a synthetic production
substitute.  Small-grain/PAH stochastic emission remains outside this bundle.

The next later promotion gate is a source/mixture-specific physical review of
the admitted thermal/emissivity asset and a matched transport convergence
study.  Native/live coupling follows only after that review.

## Over-instrumentation fence

Do not add another audit matrix, full cosmological run, AMR/MPI harness,
independent dust database, or duplicate source ledger here.  The minimum
implementation is one loader, one fixed-shape thermal operator, one P5 opt-in
path, and the compact tests above.  Existing DUST-1 evidence and provenance
helpers are reused.

## Entry and exit criteria

**Entry:** DUST-1 scattering/state bundle is committed and pushed, with its
conditional status and deferred physical assumptions recorded.

**Exit:** a traceable thermal sidecar contract is either admitted for the
candidate controls or explicitly rejected; the static one-pass grain thermal
operator and P5 outputs pass the compact gates; no DUST-1 or H/He ledger is
double-counted; and the result remains labelled `conditional_candidate` until
the physical dust mixture/emissivity review is complete.

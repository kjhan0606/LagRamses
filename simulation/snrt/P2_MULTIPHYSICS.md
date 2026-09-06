# P2 multiphysics closure

## Dust

`DustModel` accepts group absorption cross sections in `cm^2/H` and a cell-wise relative dust abundance. No Milky-Way opacity is hard-coded: production runs must provide a redshift- and metallicity-appropriate table, for example sampled from Draine model data. Dust receives the absorption-weighted energy of the photons absorbed by the dust fraction and does not change H/He ion fractions; conservative scattering adds no dust heating.

The audited sidecar loader is documented in
[`P4_DUST_OPACITY.md`](P4_DUST_OPACITY.md). It additionally accepts a
dust-absorption-weighted photon energy per group, which is used for dust
heating and absorption momentum. A v3 candidate also stores the
scattering-weighted energy and uses a named isotropic within-group closure for
the separate scattering ledger and momentum diagnostic. The P4 and P5 runners
require `--dust-opacity-metadata` whenever the static input has non-zero dust
abundance; otherwise they use an explicit zero-dust control. The momentum
diagnostic is not yet coupled to hydro or a full radiation-pressure closure.

The current scattering candidate is deliberately isotropic and records the
measured Draine angular moments without claiming to use an anisotropic phase
function. The DUST-2 candidate in
[`P6_DUST_THERMAL.md`](P6_DUST_THERMAL.md) admits a Kirchhoff-derived
single-temperature table, records tracked and out-of-band IR energy and
thermal photon rates inside the P5 subcycle, and includes the CMB floor. It is
explicitly one-pass (`recorded_not_transport_reemitted`) and cannot alter
gas/H/He state; recursive IR transport, stochastic heating, and physical
mixture approval remain deferred.

## X-ray secondaries

The production closure now bilinearly interpolates the Furlanetto--Stoever
(2010) fast-electron tables in electron energy and H II fraction. It tracks
heating, H I/He I/He II secondary ionization, and excitation as separate energy
channels and renormalizes their small Monte-Carlo/table-rounding residual to
exact local energy closure. Below the 10 eV table floor all energy becomes
heat; above 9937.21 eV the asymptotic tabulated fractions are held fixed. There
is no special branch at 100 eV.

Secondary ionization is opt-in at the solver API. If a tabulated target
species is numerically absent in the actual cell, its assigned energy is
returned to heat rather than creating an ionization of a nonexistent species.
The availability mask is fixed from the timestep's incoming species inventory
so the global opacity iteration does not acquire a channel-switching
discontinuity as He II is first created.
Excitation energy is treated as escaping line radiation and is deliberately
not returned to the gas thermal budget; P4/P5 record this policy and enforce a
separate photoelectron-energy ledger.

The tables assume primordial H/He, equal H II and He II fractions, and
negligible He III. Applying them to a state that departs strongly from that
composition is an explicit approximation, not a new arbitrary-composition
fit. The exact source, license, hashes, continuity result, and P5 effect
measurement are documented in
[`SECONDARY_IONIZATION_VALIDATION.md`](SECONDARY_IONIZATION_VALIDATION.md).

## Thermal/hydro source term

`ThermalState` stores gas internal-energy density. `advance_radiative_energy()` applies net radiative heating/cooling as a conservative local source term and returns an ideal-gas temperature. Hydrodynamic fluxes, shocks, and metal cooling remain the host hydro solver's responsibility.

## Validation

- A unit-optical-depth dust cell transmits `exp(-1)` and receives the exact missing photon energy.
- A primordial neutral H/He cell exposed to 1 keV photons produces
  `x_HII=0.090666` from the converged FS2010 primary-plus-secondary rate
  update, with both He ledgers and the photoelectron-energy ledger checked.
- The P2/P3 validation script also checks finite gas-heating output.

The production multiphysics path now uses the B2 C2-Ray-style
time-averaged-opacity closure documented in
[`B2_PRODUCTION_SOLVER_VALIDATION.md`](B2_PRODUCTION_SOLVER_VALIDATION.md).
The former per-cell atom-inventory attenuation cap is retired.

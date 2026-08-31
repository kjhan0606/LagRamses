# P2 multiphysics closure

## Dust

`DustModel` accepts group absorption cross sections in `cm^2/H` and a cell-wise relative dust abundance. No Milky-Way opacity is hard-coded: production runs must provide a redshift- and metallicity-appropriate table, for example sampled from Draine model data. Dust receives the full absorbed photon energy and does not change H/He ion fractions.

Scattering and IR re-emission are intentionally not in this live transport step. They require a separate angle- and frequency-coupled closure rather than treating scattering as absorption.

## X-ray secondaries

For photoelectron energy at or above 100 eV, the baseline uses the high-energy Shull--van Steenberg fractions for heating, H I ionization, He I ionization, and excitation. The code keeps excitation as a separate energy channel. Below 100 eV all excess energy becomes heat and no secondary ionization is added.

The more accurate Furlanetto--Stoever treatment depends on both energy and ionized fraction in a non-separable way. It remains the planned table-interpolation replacement for this baseline; the present analytic closure must not be used to claim precision X-ray spectra below its stated range.

## Thermal/hydro source term

`ThermalState` stores gas internal-energy density. `advance_radiative_energy()` applies net radiative heating/cooling as a conservative local source term and returns an ideal-gas temperature. Hydrodynamic fluxes, shocks, and metal cooling remain the host hydro solver's responsibility.

## Validation

- A unit-optical-depth dust cell transmits `exp(-1)` and receives the exact missing photon energy.
- A neutral H cell absorbing 1 keV photons produces `x_HII=0.185493` from primary plus secondary ionizations.
- The P2/P3 validation script also checks finite gas-heating output.

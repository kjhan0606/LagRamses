# P0 Requirements: TPU S_N Radiation Transport

Status: Draft v0.1  
Date: 2026-08-28  
Scope owner: Paper-III dual-AGN project

## 1. Scientific objective

Quantify whether the M1 closure changes the physical and observational
interpretation of high-redshift massive SMBH and dual-AGN systems. The solver
must resolve anisotropic AGN radiation, mutual illumination, dense-cloud
shadows, and energy-dependent absorption more faithfully than M1.

The primary comparison is between:

1. no AGN radiation;
2. Kimm-style multigroup M1 radiation; and
3. multigroup discrete-ordinates (S_N) radiation.

The target observables are the SMBH accretion history, dual-AGN duty cycle,
gas column-density distribution, ionized-bubble topology, and inputs for
JWST and X-ray synthetic observations.

## 2. Scope and non-goals

P0 defines a standalone JAX/XLA solver for fixed, block-structured Cartesian
domains. P1 begins with static radiative transfer on converted lagRamses zoom
snapshots. A live RAMSES-to-TPU call at every AMR step is explicitly out of
scope: dynamic AMR, host-device transfers, and asynchronous load balancing
are incompatible with the initial static TPU execution model.

P0 does not replace lagRamses hydrodynamics and does not replace SKIRT for
Monte Carlo observational post-processing. A coupled local RHD re-simulation
is a later milestone, contingent on a material M1--S_N difference.

## 3. Numerical contract

### 3.1 Mesh and parallel model

- Static, nested Cartesian blocks; no runtime mesh refinement.
- Default block shape: 32 x 32 x 32 active cells with fixed halo layers.
- Each production domain is padded to a compile-time block layout. Inactive
  cells use masks rather than shape changes.
- The spatial block axis is sharded across TPU cores. Photon-group and
  angular axes are static compile-time dimensions.
- Input and output checkpoints use HDF5 in CGS units. A lagRamses converter
  performs conservative AMR-to-nested-grid rasterization before a TPU run.

### 3.2 Angular transport

- Use level-symmetric discrete ordinates quadratures.
- Development quadrature: S4, 24 directions.
- Reference quadrature: S6, 48 directions.
- Angular-convergence check: S8, 80 directions, on selected subdomains.
- Start with a positivity-preserving first-order upwind finite-volume scheme;
  add a monotonic second-order reconstruction only after the reference tests
  pass.
- Directions are processed in fixed batches of eight and reduced into photon
  number, energy, and flux moments. The quasi-static solver therefore does
  not retain all G x N_omega intensity fields in memory simultaneously.
- The initial transport mode is explicit with a reduced-speed-of-light study.
  The production candidate is backward-Euler transport solved by a matrix-free
  Jacobi iteration plus diffusion-synthetic acceleration for optically thick
  groups.

### 3.3 Photon groups

The baseline contains nine groups. Group-averaged cross-sections and excess
photoheating energies are regenerated from the adopted source SEDs.

| Group | Energy range | Principal interaction |
| --- | --- | --- |
| IR | 0.01-1.0 eV | Dust absorption, re-emission, IR pressure |
| OPT | 1.0-5.6 eV | Dust absorption and direct pressure |
| FUV | 5.6-11.2 eV | Dust and photoelectric heating |
| LW | 11.2-13.6 eV | H2 dissociation |
| HI | 13.6-24.59 eV | HI ionization |
| HeI | 24.59-54.42 eV | HI and HeI ionization |
| EUV | 54.42-500 eV | HeII ionization and photoheating |
| SX | 0.5-2.0 keV | X-ray ionization and heating |
| HX | 2.0-10.0 keV | Deep X-ray penetration and heating |

Photons above 10 keV are excluded from P1 and their luminosity is recorded as
an escaping-energy diagnostic. They may become a separate group only if the
AGN SED and Compton-heating tests show that they affect the target gas.

### 3.4 Sources and gas coupling

- Stellar sources: tabulated age- and metallicity-dependent SEDs.
- AGN sources: luminosity from instantaneous BH accretion, with a separately
  recorded intrinsic SED and radiative efficiency.
- Initial chemistry state: H, H+, He, He+, He++, H2, and electrons.
- Required gas processes: photoionization, photoheating, H2 dissociation,
  dust absorption, radiation pressure, X-ray secondary ionization, and
  Compton heating.
- Every update records absorbed photon number, deposited thermal energy, and
  radiation momentum. Conservation is checked globally and per group.

## 4. Data contract

The rasterized input must contain density, temperature, velocity, metallicity,
dust-to-metal ratio, ion fractions, and static block geometry. A source table
must contain source ID, type, position, velocity, luminosity per group, and
the source time stamp. The solver outputs group-integrated photon density,
flux, photoionization/heating rates, radiation force, chemistry increments,
and the conservation ledger.

The initial science configuration is a quasi-static solve on each selected
snapshot. Time-dependent transport is enabled only after the AGN variability
time, radiation crossing time, and local recombination time have been compared
for the selected domain.

## 5. Verification plan

Each milestone must retain a deterministic reference result and a conservation
ledger. The required tests are:

1. analytic Stromgren sphere and expanding ionization-front tests;
2. single opaque clump shadow;
3. two-source crossing-beam and mutual-shadow test;
4. dusty slab with radiation pressure;
5. X-ray illuminated neutral cloud with a fixed microphysical reference; and
6. clumpy dual-AGN nuclear domain compared with a long-characteristic solve.

For selected scientific domains, the S6-to-S8 change in ionized volume,
escape fraction, radiation force, and mass-weighted temperature must be below
10 percent. The M1-to-S6 difference is reported without imposing a direction.

## 6. P0 deliverables and gates

P0 is complete when the following are fixed:

1. this numerical and physical contract;
2. canonical benchmark initial conditions and reference outputs;
3. the exact lagRamses snapshot fields needed by the rasterizer;
4. a TPU memory model for S4 and S6 angular batching; and
5. the criteria for promoting snapshot post-processing to local RHD.

Gate A occurs after the first real zoom snapshot comparison. If S6 changes
the key dual-AGN or LRD selection metrics by less than 10 percent relative to
M1, the project stops at validated post-processing. Otherwise, P4 local RHD
is authorized.

## 7. Risks

- S_N ray effects: control with S4/S6/S8 convergence, angular filtering, and
  rotated quadrature experiments.
- TPU memory: use direction batching and nested domains; do not allocate the
  full angular intensity field when a quasi-static solve is sufficient.
- Stiff chemistry: use cell-local subcycling and a positivity-preserving
  semi-implicit update before attempting a fully coupled implicit solve.
- Subgrid uncertainty: vary dust-to-metal ratio, AGN SED, and radiative
  efficiency independently of the transport-method comparison.

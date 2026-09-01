# B1 thermal-coupling disposition

Status: **B1 PASS** on 2026-09-01. AGY's first audit found an off-grid
metallicity defect; format v3 was implemented, locally revalidated, and passed
the required B1-only AGY re-audit. Upstream data-license approval remains
independently pending.

## Defect closed

The migrated v1 atlases were generated with Grackle flags
`primordial_chemistry=1`, `metal_cooling=1`, and `UVbackground=1`. They stored
the complete equilibrium thermal rate after solving the H/He state. Adding
that rate to SNRT photoheating while separately evolving non-equilibrium H/He
double-counted primordial and UVB physics and made the rate inconsistent with
the live ion fractions.

`read_thermal_atlas` now accepts only format v3 with a complete provenance
group and the exact component contract:

- `thermal_component=metal_only`
- `primordial_rates_included=false`
- `uv_background_included=false`
- `photoheating_included=false`
- `metallicity_application=analytic_runtime_multiplier`
- heating-positive/cooling-negative sign convention

The migrated `p4_validation_thermal_atlas.h5` and
`p6_thermal_atlas_a0p20849776.h5` are intentionally rejected. The first
metal-only `production_metal_thermal_atlas_v1.h5` is also quarantined: it
tabulated a log-metallicity axis and linearly interpolated rates in that
coordinate, violating the declared linear-Z model between grid points.

## Runtime model

The non-photo thermal source is now

`R_nonphoto = R_HHe(x_HII, x_HeII, x_HeIII, n_e, T, a) + R_metal(n_H, Z, T, a)`.

`R_HHe` is evaluated directly from the live SNRT state. It includes
collisional excitation, collisional ionization cooling, case-B recombination
cooling, free-free cooling, and CMB Compton heating/cooling. The coefficients
are ported from Grackle 3.4.2-dev revision
`f93091ff8456962d7017a5bff7472945a30e3dad`. The same Abel collisional
ionization coefficients are coupled to both the multiphysics and conservative
H/He chemistry solvers; their discrete transition counts are included in the
species ledgers.

`R_metal` is generated only from `CoolingRates/Metals/Cooling` in the pinned
no-UVB Grackle data file. The atlas stores one solar-metallicity `(a, n_H, T)`
table and the host and JAX runtimes apply `Z/Zsun` analytically, including
exactly zero metallicity. It continuously subtracts the cooling coefficient at
`T_CMB`; this deliberately removes Grackle revision `f93091f`'s two-dex cutoff
optimization and its finite step, and is declared in provenance. The atlas
contains no primordial table, UVB photoionization, or photoheating.

## Pinned assets

- Grackle data repository revision:
  `928696482fbe15d9bac4382de6134d95568f099c`
- `CloudyData_noUVB.h5` SHA-256:
  `0abe25cceeb5c0825381c5f17059982a9a2cdd27ce369a475c559fba6a8fa106`
- generated atlas:
  `data/production_metal_thermal_atlas_v2.h5`
- generated atlas SHA-256:
  `b1290d930b22ed049d6d3c5ed47ce56ecf3e0d2e693b39792740768e80fdf6ac`
- generator SHA-256 embedded in the atlas:
  `6fe00f80795ac948ed2386512011014932f559c48e7a1455116031852c0f5280`

## Acceptance evidence

`tests/b1_thermal_coupling.py` verifies the Grackle ionization coefficients,
an `x_HII` cooling-rate sweep, CMB sign reversal, discrete H/He chemistry
ledgers, host/JAX analytic metallicity scaling at off-grid scalar and vector
values, continuous high-temperature CMB subtraction, legacy-atlas rejection,
provenance tamper rejection, and rejection of the HM2012 UVB source as a
metal-only input.

The 32-cubed conservative coarse/fine run in
`data/b1_validate_conservative_primordial_thermal_v2.json` passed with:

- mean coarse/fine `|Delta log10 T| = 2.37063e-4`
- `x_HII` L1 difference `= 6.45379e-5`
- maximum fine fixed-point residual `= 4.47035e-7`
- fine H/He I/He II ledger relative errors
  `= 5.77972e-6 / 3.76322e-7 / 1.12988e-8`
- no temperature-ceiling hits

The first AGY B1 audit and its blocked disposition are recorded in
`../../provenance/agy_b1_thermal_coupling_audit_2026-09-01.md`; the successful
closure audit is in
`../../provenance/agy_b1_thermal_coupling_reaudit_2026-09-01.md`. The
code/test/provenance audit condition is closed. Production promotion remains
blocked by the independent upstream data-license approval.

## Remaining scope limits

The metal table assumes the no-UVB Cloudy ionization model and linear solar
abundance scaling. It does not evolve individual metal ions in the local SNRT
radiation field. That approximation must be disclosed and tested as a model
sensitivity before a paper claim that depends strongly on metal-line cooling.
The atlas mean molecular weight is a neutral-primordial initialization aid;
runtime heat capacity comes from the live H/He state. Production hydro inputs
should carry a usable temperature/pressure and metallicity field rather than
relying on the atlas equilibrium-temperature fallback.

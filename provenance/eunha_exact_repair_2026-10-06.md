# Equilibrium exact cooling repair — 2026-10-06

Repository: `/gpfs/kjhan/chimes-cvode-gpu-20261001`, origin
`kjhan0606/LagRamses`. Scope: repair the existing equilibrium Grackle-table
thermal integrator, not introduce a chemical network or change rate data.

## Findings and correction

- Removed interpolation of independently normalized corner Y functions.
  Bilinear interpolation in log density/metallicity coordinates now precedes
  construction of the specific-energy/net-rate curve. Each linear energy
  segment evolves analytically under `du/dt=-Lambda(u)`.
- Cooling and heating use the same curve. Its zero-rate roots are approached
  asymptotically without crossing; energy-knot traversal handles long steps.
- Evolve RAMSES `T/mu` as specific energy directly. Do not use the old
  three-iteration physical-temperature reconstruction in the update.
- Removed stored Y tables and their redshift-dependent reconstruction.
- The bundled multi-z binary retains a placeholder redshift followed by six
  log-domain limits. The reader previously treated that placeholder as the
  density minimum, shifting all limits. This was reproduced with floating
  point traps in the mass-density conversion and fixed in the reader.
- Explicitly reject exact/compare with `-Dgrackle`, a non-5/3 gamma, or
  `neq_chem=true`, instead of silently bypassing or mixing models.

## Verification

`module load gnu13`

`bash simulation/snrt/tools/test_eunha_exact.sh /gpfs/kjhan/chimes-cvode-gpu-20261001/patch/cuRamses/aux/grackle_multi_z.bin`

Passed with array bounds checks and floating-point invalid/zero/overflow
traps. These are bounded single-thread unit tests, not RAMSES evolution.

- Six synthetic curves: constant cooling, constant heating, stable
  equilibrium, exponential cooling, varying table molecular weight, and
  corner-dependent equilibrium positions. Analytic agreement and split-step
  agreement required relative error <=1e-10.
- Long-step equilibrium approach from both directions, boundary stops, and
  zero-timestep identity passed.
- Existing multi-z table: 48 states at z=0,3,6, four densities and four
  energies; finite positive endpoints and split-step agreement passed.
  Confirmed table ranges: nH=1e-8..1e4 cm^-3, Z/Zsolar=1e-6..10,
  T=10..1e9 K, ten redshift snapshots spanning z=0..10.
- Initial Slurm unit-test submissions 414463/414464 were pending and cancelled
  after verifying their exact WorkDir; no science jobs were affected.
- `git diff --check` passed.

## Limits

Exactness is for the declared piecewise-linear specific-energy curve, not
the continuous physical rates. Density, abundance and background are frozen
within a source update. This path remains equilibrium chemistry with total
metallicity, not nine-species NLTE or element-resolved cooling. Boundary
clamping remains: zero metallicity uses the minimum-metallicity plane and
z>10 uses the highest-z plane. UV, shielding and dust assumptions of the
rate-generation data still require confirmation before a dust/RT-off science
launch. No complete RAMSES rebuild/evolution or performance claim is made
by these unit results. No namelist keys/defaults were added or changed.

# B3 timestep and spatial-resolution validation

Date: 2026-09-02
Stage status: canonical 2-by-2 matrix running; Claude Opus 5 audit required

## Predeclared contract

Before the canonical matrix was launched, the acceptance contract was fixed in
`config/b3_timestep_spatial_gate.json`. The primary comparison is the 32 cubed,
Courant 0.1 baseline against a factor-two 64 cubed prolongation at the same
Courant number. The finer cell width makes the physical transport timestep
exactly one half. The gate requires the relative change in volume-mean xHII to
be below 2%, with the denominator defined as the larger absolute mean.

The complete matrix is:

| Case | Grid | Courant | Physical timestep relative to baseline | Role |
| --- | ---: | ---: | ---: | --- |
| `n32_courant0p1` | 32 cubed | 0.10 | 1 | baseline |
| `n32_courant0p05` | 32 cubed | 0.05 | 1/2 | coarse time-only control |
| `n64_courant0p1` | 64 cubed | 0.10 | 1/2 | primary simultaneous dx/2 and dt/2 result; same-dt spatial control against the preceding row |
| `n64_courant0p05` | 64 cubed | 0.05 | 1/4 | fine time-only control |

All runs use the canonical nine photon groups, 0.1 Myr, S4, reduced light
speed 0.01c, float64, FS2010 secondary ionization, point deposition, source
target 0.25, one thermal subcycle, 24 thermal iterations, and 32
time-averaged-absorption iterations. Every run must independently pass the
fixed-point, H/He ledger, thermal closure, photoelectron-energy, root, bound,
fraction, and finite-value gates in the contract.

## Spatial-control interpretation

The 64 cubed input is a conservative, piecewise-constant factor-two
prolongation of the same staged 32 cubed gas state. Each coarse source is split
equally among its eight children, preserving total group luminosity and the
centre of luminosity. This isolates static solver resolution sensitivity but
is not an independently reconstructed higher-resolution hydrodynamic
snapshot. It therefore cannot establish hydro-resolution convergence or live
radiation-hydrodynamic feedback convergence.

## Execution provenance

The canonical array is Slurm job `330190`, launched by
`b3_timestep_spatial_matrix.sbatch` from repository commit
`7e6dab63d87707dc4ee1749f242d3a809191cc00`. Each task checks the commit and
SHA256 of the coarse input, refined input, photon metadata, and thermal atlas
before execution. The matrix results, quantitative comparison, artifact
hashes, and final Opus 5 verdict will be added only after all four tasks
finish.

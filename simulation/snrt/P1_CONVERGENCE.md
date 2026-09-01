# P1 convergence and chemistry validation

## Numerical choices fixed by this phase

- Use exact local attenuation after each explicit transport/source update. The time step is limited by the directional transport CFL, not by `c dt kappa`.
- Use S8 (80 directions) as the minimum production quadrature for source-shadow geometry. S4 and S6 remain useful low-cost diagnostics only.
- Use at least 64 cells across the 128 pc B01 reference cube at the current first-order upwind spatial discretization.

## B03 angular convergence

Fixed 48-cells-per-side geometry; clump optical depth is four; the observable is downstream transmission relative to the transparent calculation.

| Quadrature | Directions | Transmission |
| --- | ---: | ---: |
| S8 | 80 | 0.1333715 |
| Gauss-Legendre x azimuth | 128 | 0.1337046 |
| Gauss-Legendre x azimuth | 192 | 0.1339602 |

S8 differs from the 192-direction reference by 0.44%. The 128-to-192 direction change is 0.19%.

## B01 spatial and analytic convergence

The box side is fixed at 128 pc and the integration time at 22,265.6 yr. The
corrected control uses `ĉ=0.03c` and 12 steps per linear cell count, so every
grid reaches the same physical duration. At that time the infinite-light-speed
hydrogen-only solution is `R(t)/R_S=0.5501305`. Helium fields remain present in
the state layout, but the 20 eV fixture lies below the He I ionization threshold;
helium therefore remains neutral and does not explain the residual. This is an
RSLA and numerical-discretization comparison against the H-only analytic result.

| Grid | Rion / RS | Rion / Ranalytic(t) |
| --- | ---: | ---: |
| 32 cubed | 0.5144022 | 0.9351 |
| 48 cubed | 0.5114129 | 0.9296 |
| 64 cubed | 0.5094979 | 0.9261 |

Every grid is within the declared 10% analytic-radius threshold. The maximum
spatial spread is 0.96% relative to the 64-cells-per-side result. The older
`0.01c`/four-steps-per-cell values were self-consistent but about 20% below the
time-dependent analytic radius; they are superseded rather than promoted as a
convergence result.

## Photon and chemistry checks

- A single opaque cell with unit optical depth absorbs `1-exp(-1)=0.63212055` photons, and H II rises by precisely the same amount.
- A one-zone H/He photoionization step with `Gamma dt=1` gives `x_HII=x_HeII=0.63212055` and leaves He III at zero.
- A zero-photo-rate recombination step keeps H and He fractions finite and the helium fractions closed.

## Remaining P1 boundary

This P1 check covers a historical primordial UV path. The production
multiphysics RSLA matrix and its explicit `0.01c` error bound are documented in
[`RSLA_REFINEMENT_VALIDATION.md`](RSLA_REFINEMENT_VALIDATION.md). Dust,
scattering, metal chemistry, temperature evolution, and hydrodynamic coupling
are outside P1.

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

## B01 spatial convergence

The box side is fixed at 128 pc, the integration time is 22,265.6 yr, and the number of steps is scaled with the linear grid size.

| Grid | Rion / RS |
| --- | ---: |
| 32 cubed | 0.4368803 |
| 48 cubed | 0.4399448 |
| 64 cubed | 0.4393452 |

The maximum spread is 0.70% relative to the 64-cells-per-side result.

## Photon and chemistry checks

- A single opaque cell with unit optical depth absorbs `1-exp(-1)=0.63212055` photons, and H II rises by precisely the same amount.
- A one-zone H/He photoionization step with `Gamma dt=1` gives `x_HII=x_HeII=0.63212055` and leaves He III at zero.
- A zero-photo-rate recombination step keeps H and He fractions finite and the helium fractions closed.

## Remaining P1 boundary

This validation establishes the static, primordial UV core. Dust, scattering, X-ray secondary ionization, metal chemistry, temperature evolution, and hydrodynamic coupling remain P2 work.

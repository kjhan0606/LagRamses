# NGC 2403 rotation curves at fixed baryonic geometry

We compare the QUMOND circular speed with 64 SPARC measurements between 1.06 and 20.87 kpc. We use the same observed radii and define the residual as model speed minus observed speed. The three radial intervals contain 32, 10, and 22 measurements. Each mean gives equal weight to each radius within the interval. The observational errors and the radial correlations enter neither the means nor a likelihood. [SPARC](https://doi.org/10.3847/0004-6256/152/6/157)

## Periodic boundary

We keep the physical cell width at 0.125 kpc and enlarge the box along all three axes. The stars follow the SPARC light profile. The gas follows 44 H I measurements digitized from Fig. 9 of Fraternali et al. (2002). We multiply the H I density by 1.33 for helium and taper the profile to zero between 21.86 and 23.86 kpc. The FR02 gas profile approximates the observed gas distribution. The stellar and gas exponential scale heights are 0.2414 and 0.2000 kpc. The sampled gas mass is 4.2431 billion solar masses in every run. [Gas observations](https://doi.org/10.1086/340358) · [Source parameters](../rar_ngc2403_fr02_recovery_20261108/gate512/source_contract.json)

We solve the discrete Newtonian and QUMOND field equations for the prescribed density, with zero imposed external field and the simple interpolating function at an acceleration scale of 1.2 × 10⁻¹⁰ m s⁻². Reflection symmetry of the centered disc permits inversion of the periodic seven-point Laplacian in one octant. We evaluate the midplane inward acceleration at 64 azimuths per observed radius and convert the median to a circular speed. At 64 kpc both fields reproduce the compiled lagRamses rotation curves for both stellar mass-to-light ratios (M/L) to within 7 × 10⁻¹¹ km/s at every radius. The larger boxes use the same uniform-grid operators.

| Box side (kpc) | QUMOND speed at 20.87 kpc (km/s) | Mean residual at 1.06–5 kpc (km/s) | Mean residual at 10–20.87 kpc (km/s) |
| ---: | ---: | ---: | ---: |
| 64 | 118.94 | +7.34 | −8.85 |
| 96 | 124.13 | +7.39 | −6.35 |
| 128 | 125.27 | +7.40 | −5.79 |
| 192 | 125.82 | +7.41 | −5.52 |
| 256 | 125.94 | +7.41 | −5.46 |

The table uses M/L = 0.7 in solar units at 3.6 μm. The outermost observed speed is 134.0 ± 3.12 km/s. Increasing the box from 64 to 256 kpc raises the modeled speed by 7.01 km/s at the last radius. Increasing the box from 192 to 256 kpc changes the speed by at most 0.125 km/s across the 64 radii for M/L = 0.5 and 0.7. The last change measures sensitivity to the box size at fixed cell width. Vertical resolution and the limiting isolated solution require separate tests. [Residuals by radius](box_size/residuals.pdf) · [Numerical values](box_size/zones.csv)

At 256 kpc the SPARC distance of 3.16 ± 0.16 Mpc and inclination of 63° ± 3° give a mean residual range of +1.79 to +12.55 km/s in the inner interval for M/L = 0.7. The outer interval ranges from −12.57 to +1.01 km/s. We vary the two quantities coherently across all radii, hold the angular matter profile and angular box size fixed, and retain the adopted deprojection for the mass map. Distance then multiplies the model speed by the square root of the distance ratio. Inclination rescales the observed speed by the ratio of sines. The extrema describe the specified rectangle in distance and inclination. The residuals retain additional uncertainty from the adopted gas profile, spatial variation in M/L, and radial covariance.

The NFW curve in the residual figure holds the original halo contribution fixed in each of four radial folds. Each fit used the 48 radii outside the 16-radius test block in NGC 2403. The stellar M/L is 0.5. Only the Newtonian baryonic term changes with the box. Cosmological predictions for the halo mass and concentration require independent constraints. [Original NFW fits](../rar_ngc2403_fr02_recovery_20261108/comparison/results/comparison.csv)

## Stellar mass from photometry

De Blok et al. report single-disc mean M/L values of 0.29 for a Kroupa initial mass function and 0.41 for a diet-Salpeter initial mass function in NGC 2403. We apply each value uniformly to the same SPARC light profile and solve both field equations again in the 256 kpc box. The two values represent alternative stellar population assumptions. The original photometry also shows a radial colour gradient. The constant M/L calculations omit the radial variation in the original mass model. [De Blok et al. (2008), Sect. 6.2 and Tables 3–4](https://doi.org/10.1088/0004-6256/136/6/2648)

| M/L | Mean residual at 1.06–5 kpc | Mean residual at 5–10 kpc | Mean residual at 10–20.87 kpc | Speed at 20.87 kpc |
| ---: | ---: | ---: | ---: | ---: |
| 0.29 | −20.95 | −24.87 | −19.44 | 113.39 |
| 0.41 | −11.07 | −17.04 | −14.93 | 117.42 |
| 0.50 | −4.72 | −11.86 | −11.80 | 120.23 |
| 0.70 | +7.41 | −1.82 | −5.46 | 125.94 |

All entries are in km/s, apart from M/L. At M/L = 0.41 the mean residual remains negative in every interval across the distance and inclination ranges. Even the most favorable common geometry gives −6.39 km/s in the inner interval and −8.69 km/s in the outer interval. [Rotation curves and residuals](population/rotation_curves.pdf) · [Regional statistics](population/zones.csv) · [Accelerations and individual radii](population/rings.csv)

We test a common velocity scale by measuring the factor needed to set the mean residual to zero separately in each interval. The factor is the sum of observed speeds divided by the sum of QUMOND speeds over the same radii. At M/L = 0.41 the inner and outer intervals require factors of 1.1249 and 1.1264. The two averages admit almost the same normalization. The catalogued distance and inclination permit a maximum relative increase of about 5.1% under the adopted geometry. At M/L = 0.7 the inner interval requires a factor of 0.9308 and the outer interval requires 1.0428. A common rescaling preserves the ratio of the two factors. Comparing the inner and outer intervals tests the curve shape independently of a common speed normalization. Individual radii retain residuals after matching an interval mean. [Velocity scale at each M/L](velocity_scale.csv)

## Reproduction

| Grammar CPU job | Source commit | Calculation |
| --- | --- | --- |
| 1324354 | `52ef919` | 64, 96, 128, and 192 kpc at M/L = 0.5 and 0.7 |
| 1324357 | `2893519` | 64 and 256 kpc at both M/L values |
| 1324358 | `2893519` | Radial residuals and coherent geometry changes for all five boxes |
| 1324370 | `668ce59` | Independent M/L = 0.29 and 0.41 field solves, 64 kpc reference check, and four-curve comparison |
| 1324373 | `a717f9e` | Inner and outer velocity scales and logarithmic acceleration residuals |

All five jobs completed with exit code 0. The two grid records preserve the masses and the comparison with the compiled lagRamses reference. The residual and population directories each include a checksum manifest for the original result tables and figures. The calculations use a prescribed disc at a single epoch. The conclusions apply to the adopted gas profile, stellar distribution, and force law for NGC 2403. A comparison of cosmological galaxy formation requires live galaxy simulations.

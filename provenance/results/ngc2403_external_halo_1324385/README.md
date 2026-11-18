# A photometry-based halo for NGC 2403

We compare the NGC 2403 rotation curve with two prescribed mass models. The stars and gas have the same density distribution in both models at each stellar mass-to-light ratio (M/L). The Newtonian model adds a spherical NFW dark halo. QUMOND uses the existing three-dimensional field solution with zero imposed external field. The comparison covers the same 64 SPARC radii at a fixed cell width of 0.125 kpc in a 256 kpc periodic box. [Stellar, gas, and force inputs](../ngc2403_fr02_boxes_20261118/README.md)

## Halo calibration

We fix the halo using the NGC 2403 single-disc Kroupa normalization of 0.29 in the 3.6 μm band. The resulting stellar mass is 2.9183 billion solar masses. Dividing by 1.06 places the stellar mass on the Chabrier scale. We invert the central-galaxy relation of Moster et al. at redshift zero using Table 1 and Eq. (2). At the SPARC distance of 3.16 Mpc the resulting total-matter-equivalent mass in a collisionless simulation is 1.5809 × 10¹¹ solar masses within 111.12 kpc. The mean enclosed density is 200 times the critical density for a Hubble constant of 70.4 km/s/Mpc. [Photometric calibration](https://doi.org/10.1088/0004-6256/136/6/2648) · [IMF conversion, Eq. 2](https://doi.org/10.1088/0067-0049/214/2/15) · [Stellar-to-halo mass calibration](https://arxiv.org/pdf/1205.5807)

The WMAP7 concentration relation gives a concentration of 8.4135 for the mass-selected population at redshift zero. We assign a fraction 1 − 0.046/0.272 of the enclosed NFW mass to dark matter and add the computed Newtonian force of the observed disc. The calculation assumes an uncontracted halo and omits a separate diffuse baryonic component. One stellar photometry anchor sets the halo for all four disc M/L scenarios. Changing the disc M/L at a given distance tests the force response with the halo held fixed. [Concentration calibration, Eq. 24 and Table A5](https://arxiv.org/pdf/1411.4001#page=19)

We vary distance and inclination over the catalogued central value and both error endpoints. Distance rescales all physical disc lengths, the physical ring radii, and the anchor stellar mass at fixed angular profile and surface brightness. Inclination changes the deprojection of the observed velocity. We retain the adopted deprojection of the mass map. We also vary stellar mass at fixed halo mass by ±0.15 dex in the forward population relation and concentration at fixed mass by ±0.11 dex. The concentration offset uses the scatter reported by Dutton and Macciò as a sensitivity scale. Every choice applies coherently to all radii. The extrema summarize this parameter grid. Probabilistic inference would also require a halo mass prior, the selection function, and observational covariance. [Concentration scatter](https://arxiv.org/pdf/1402.7073) · [All parameter combinations](scenarios.csv)

## Velocity and acceleration residuals

The mean velocity residual is the unweighted mean of model minus observed speed, in km/s. The radial intervals contain 32, 10, and 22 measurements. The table gives the NFW and QUMOND means in each interval at the nominal distance, inclination, and halo parameters.

| M/L | NFW 1.06–5 kpc | QUMOND 1.06–5 kpc | NFW 5–10 kpc | QUMOND 5–10 kpc | NFW 10–20.87 kpc | QUMOND 10–20.87 kpc |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.29 | −26.46 | −20.95 | −37.02 | −24.87 | −39.79 | −19.44 |
| 0.41 | −18.86 | −11.07 | −32.64 | −17.04 | −37.85 | −14.93 |
| 0.50 | −13.61 | −4.72 | −29.49 | −11.86 | −36.43 | −11.80 |
| 0.70 | −2.95 | +7.41 | −22.85 | −1.82 | −33.34 | −5.46 |

At M/L = 0.7 the full NFW grid gives outer mean residuals between −46.41 and −19.26 km/s. The QUMOND geometry grid gives −12.57 to +1.01 km/s. In the same outer interval the nominal mean logarithmic acceleration residuals are −0.251 dex for NFW and −0.036 dex for QUMOND. The plot shows measured velocity errors in grey and the NFW population-relation sensitivity at the nominal geometry in blue. The lower panels show the nominal acceleration residuals against the computed baryonic acceleration. [Comparison figure](comparison.pdf) · [Regional results](zones.csv) · [Individual radii](rings.csv)

The last ring has a speed of 134.0 ± 3.12 km/s at 20.87 kpc. For a spherical dark halo the required mass inside the ring equals the radius times the difference of the squared observed and Newtonian baryonic speeds divided by the gravitational constant. At M/L = 0.7 the required dark mass is 7.44 × 10¹⁰ solar masses. The halo calibrated from photometry supplies 3.27 × 10¹⁰ solar masses. Imposing the adopted median concentration relation and fitting this single ring gives a total-matter-equivalent halo mass of 7.90 × 10¹¹ solar masses. Only the enclosed mass follows directly from the circular-speed decomposition. Extrapolation to the halo radius depends on the assumed density profile. The kinematic estimate enters none of the predictions. [Enclosed masses and speed-error sensitivity](endpoint_mass.csv)

Grammar jobs 1324384 and 1324385 completed the comparison. The first job used revision `df4ef0f`. The second used revision `7f3c61c` and added the outer-ring mass calculation. The halo parameters, radial predictions, regional statistics, and scenario tables agree byte for byte across both runs. [Result checksums](results.sha256)

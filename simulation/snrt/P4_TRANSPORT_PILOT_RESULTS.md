# P4 coeval static transport result

## Input contract

The coeval input is `data/p4_coeval_static_rt_input.h5`. It uses output 00017
gas and the audited output 00017 instantaneous AGN photon ledger at
`aexp = 0.208497764676753`. The 32-cubed mesh is a volume-conservative
resampling of every AMR leaf intersecting the target cube. Its coverage weights
are `0.9999999999999998` to `1.0000000000000002`.

The mesh has mean `n_H = 0.145638 cm^-3` and maximum `n_H = 44.4322 cm^-3`.
The photon ledger has ten in-cube sinks, but the ID 10 AGN dominates its
luminosity. The spectrum is the documented unobscured Sazonov-style baseline:
`nu L_nu(13.6 eV) = 0.1 L_bol` and escape fraction one.

## Angular comparison

Both runs use the same physical duration, `6.3697581 Myr`, reduced light speed
`0.01 c`, directional CFL `0.4`, fixed-temperature H/He chemistry, P3
24-iteration implicit recombination, zero dust, and P2 secondary ionization.

| Quadrature | Directions | Mean x_HII | Max x_HII | Cells with x_HII >= 0.1 | Max absolute delta x_HII versus P8x16 |
| --- | ---: | ---: | ---: | ---: | ---: |
| S8 | 80 | 0.001168 | 0.4000 | 0.0001526 | 0.0310 |
| product 8x16 | 128 | 0.001171 | 0.4310 | 0.0000916 | reference |

S8 differs from the 128-direction reference by mean absolute `x_HII` of
`1.68e-4` and a local maximum of `0.0310`. S8 is therefore the minimum
retained angular order for the next 32-cubed static P4 run. This is not a
general convergence claim beyond this grid, source spectrum, duration, and
gas state.

## Escape-fraction sensitivity

Holding gas, intrinsic SED, duration, and S8 fixed while changing only the
unresolved nuclear escape fraction from one to `0.1` changes the mean `x_HII`
from `0.001168` to `0.0001786`. The fraction of cells above `x_HII = 0.01`
falls from `0.02457` to `0.001831`.

The global maxima, both near `0.4`, occur in the dominant AGN injection cell
at index `(14, 17, 19)` with `n_H = 11.23 cm^-3`; that 1.2-kpc source cell is
unresolved and is not an escape diagnostic. Excluding all ten source cells,
the maximum falls from `0.1071` at escape fraction one to `0.04124` at `0.1`,
and no non-source cell remains above `x_HII = 0.1` in the attenuated case.
This is a sensitivity result, not an obscuration-model calibration.

## Limits

- The AGN SED normalization and escape fraction require an obscuration and
  sensitivity study before interpreting an LRD or JWST observable.
- Temperature is held fixed; heating is diagnostic and does not change
  hydrodynamics.
- No stellar photon sources, dust model, or metal chemistry are available in
  this snapshot.

## Artifacts

- `data/p4_coeval_static_rt_input.h5`: coeval P4 gas and photon-source input.
- `data/p4_coeval_static_rt_input.json`: AMR coverage and provenance.
- `data/p4_coeval_transport_s8_6p37myr.h5`: S8 result.
- `data/p4_coeval_transport_product_8x16_6p37myr.h5`: angular reference.

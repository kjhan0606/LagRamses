# LC18 failed-wind source inquiry packet

Date: 2026-09-04
Status: internal review packet; **not sent**
Scope: Boccioli & Roberti (2026) LC18 branch versus Limongi & Chieffi (2018)

Publication classification: **internal review only**. The executable
`fp1_lc18_failed_wind_crosscheck` publication gate reports
`allowed=false`, `publication_ready=false`, and `review_use_only=true` because
the CDS redistribution terms and a derived-artifact approval record are not
present. Editing this packet's label or the generated JSON cannot authorize
publication or redistribution.

## Reproducible finding

The joined review artifact is
`simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`. It joins all 108
LC18 coordinates one-to-one by initial mass, [Fe/H], and rotation. The
Boccioli--Roberti summary identifies 52 successful and 56 failed models.

- All 52 successful models have positive summary winds and nonzero release
  Wind tables. Their maximum summary-minus-Wind-table element-sum residual is
  0.007183005956 M_sun.
- All 56 failed models have a positive summary wind but an exactly zero release
  Wind table. They therefore remain unresolved; no missing wind was rebuilt.
- The independent LC18 CDS table5 endpoint gives a positive cumulative wind
  for 53 failed models and a value rounded to zero for three. Its summary-wind
  differences exceed the nominal 0.005 M_sun table5 half-bin for all 108
  models, including all 52 successful controls. The largest absolute
  difference is 1.5476 M_sun. Thus CDS is corroborating evidence that the
  failed stars underwent mass loss, but it is not a numerically identical
  replacement for the release summary.
- CDS table5 yields 845 unique phase rows after 19 exactly duplicated rows are
  collapsed. Each model has 3--8 unique phases, a PSN endpoint, strictly
  increasing cumulative age, and non-increasing total mass.
- CDS table7 supplies pre-SN structure for 96 models. The other 12 remain
  explicit nulls. Its binding energy is not interpreted as injected explosion
  energy.

## Byte identity

- Boccioli--Roberti `LC18.zip`: SHA256
  `249aea46713ab41cad7e8d7406835c205e4f02d36958e113a8f2231f81ebef5e`,
  Zenodo MD5 `c1d21fcbdf7ed200344881f8a47a211b`.
- Boccioli--Roberti `README`: SHA256
  `66dae04f90bf7b96460199a7ebbedf0126c1e70ab840e1693a8a326cd7ae2316`.
  Lines 30--38 state that `*_Wind.txt` contains wind yields reported in the
  original papers/pre-SN files and that only Post/PreSN yields are set to zero
  for failed explosions.
- Limongi--Chieffi CDS `table5.dat`: SHA256
  `c02278a38671e66c42efeb88fc6b2119ffefc72712cf7443128c09ffd5f5122d`.
- Limongi--Chieffi CDS `table7.dat`: SHA256
  `165201557b9093cab56b978e75b4d5989bc6e12491da54cc2bb26d20d99e6b15`.
- Acquisition manifest: SHA256
  `ebbbfa42357e6a729d19d839f0e1c826ef33cb4705d895b82c5e96e88493f764`.

The Boccioli--Roberti release is CC-BY-4.0. No explicit CDS catalogue
redistribution license has been identified for the LC18 comparison tables, so
the derived cross-check remains internal review evidence and is not a
production asset.

## Questions for the source authors

1. Are the all-zero LC18 `*_Wind.txt` columns at the 56 failed-explosion
   coordinates intentional? The README explicitly limits the failed-model zero
   policy to Post and PreSN yields, while the summary reports positive wind
   mass at every affected coordinate.
2. Does `M_wind` in `Summary_table_LC18.txt` come directly from the original
   LC18 pre-SN model files, and should it be used when a failed-model Wind table
   is zero? If so, is a corrected per-element/per-isotope Wind release planned?
3. Why does `M_wind` differ from `M_initial - M_total(PSN)` in the published
   LC18 CDS table5 beyond its nominal printed-mass half-bin at every coordinate?
   Are the underlying pre-SN model revisions, mass-loss definitions, or
   endpoint conventions different?
4. Can the per-model explosion energies and their semantics (diagnostic,
   asymptotic, or injected), terminal mass cuts, fallback/remnant masses, and
   the BHSN exception list be released in machine-readable form?
5. For failed models, is a remnant mass or fallback/direct-collapse convention
   available that preserves the pre-terminal wind separately from terminal
   ejecta?

## Affected failed-model rows

Columns are metallicity label ([Fe/H]), rotation in km/s, initial mass,
Boccioli--Roberti summary wind, release Wind-table element sum, CDS cumulative
wind `M_initial-M_total(PSN)`, cumulative terminal age in years, and CDS iron
core mass. `null` means table7 has no row and is not an inferred zero.

| Z | v_rot | M_initial | BR26 summary wind | BR26 Wind sum | CDS wind | CDS age | CDS Fe core |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| D (-3) | 0 | 80 | 0.8979 | 0.0 | 0.10 | 3.62455e+06 | 1.81 |
| D (-3) | 0 | 120 | 1.5476 | 0.0 | 0.00 | 3.10352e+06 | null |
| C (-2) | 0 | 13 | 0.1644 | 0.0 | 0.00 | 1.75461e+07 | 1.40 |
| C (-2) | 0 | 15 | 0.3577 | 0.0 | 0.20 | 1.43173e+07 | 1.08 |
| C (-2) | 0 | 80 | 2.1923 | 0.0 | 1.40 | 3.66965e+06 | 1.78 |
| C (-2) | 0 | 120 | 18.2290 | 0.0 | 17.00 | 3.14942e+06 | null |
| B (-1) | 0 | 80 | 40.5478 | 0.0 | 40.10 | 3.64244e+06 | 2.18 |
| B (-1) | 0 | 120 | 69.8521 | 0.0 | 69.30 | 3.11952e+06 | null |
| A (+0) | 0 | 15 | 1.8754 | 0.0 | 1.70 | 1.41598e+07 | 1.43 |
| A (+0) | 0 | 20 | 12.5349 | 0.0 | 12.46 | 9.74299e+06 | 1.10 |
| D (-3) | 150 | 13 | 1.3680 | 0.0 | 1.20 | 2.04861e+07 | 1.61 |
| D (-3) | 150 | 20 | 2.0197 | 0.0 | 1.80 | 1.13388e+07 | 1.91 |
| D (-3) | 150 | 25 | 2.3607 | 0.0 | 2.10 | 9.87237e+06 | 1.84 |
| D (-3) | 150 | 60 | 18.4853 | 0.0 | 18.10 | 4.73062e+06 | 1.73 |
| D (-3) | 150 | 80 | 26.8263 | 0.0 | 26.30 | 3.93721e+06 | null |
| D (-3) | 150 | 120 | 45.5513 | 0.0 | 45.20 | 3.33333e+06 | null |
| C (-2) | 150 | 15 | 2.1552 | 0.0 | 2.00 | 1.67899e+07 | 1.85 |
| C (-2) | 150 | 20 | 2.6047 | 0.0 | 2.40 | 1.15167e+07 | 1.91 |
| C (-2) | 150 | 25 | 0.3957 | 0.0 | 0.10 | 8.9221e+06 | 1.92 |
| C (-2) | 150 | 60 | 24.3529 | 0.0 | 24.00 | 4.63782e+06 | 2.22 |
| C (-2) | 150 | 80 | 32.0305 | 0.0 | 31.50 | 3.95572e+06 | null |
| C (-2) | 150 | 120 | 48.4396 | 0.0 | 47.70 | 3.39916e+06 | null |
| B (-1) | 150 | 20 | 6.2672 | 0.0 | 6.10 | 1.16914e+07 | 1.87 |
| B (-1) | 150 | 25 | 8.9226 | 0.0 | 8.80 | 9.07726e+06 | 1.88 |
| B (-1) | 150 | 30 | 4.3665 | 0.0 | 4.10 | 7.57182e+06 | 2.00 |
| B (-1) | 150 | 80 | 42.7682 | 0.0 | 42.40 | 3.94293e+06 | 2.32 |
| B (-1) | 150 | 120 | 77.7270 | 0.0 | 77.30 | 3.35242e+06 | null |
| A (+0) | 150 | 15 | 8.7571 | 0.0 | 8.69 | 1.69896e+07 | 2.80 |
| A (+0) | 150 | 20 | 12.3652 | 0.0 | 12.29 | 1.18501e+07 | 1.85 |
| A (+0) | 150 | 30 | 18.8768 | 0.0 | 18.80 | 7.80758e+06 | 1.61 |
| A (+0) | 150 | 40 | 26.4011 | 0.0 | 26.30 | 6.13369e+06 | 1.86 |
| D (-3) | 300 | 15 | 1.3479 | 0.0 | 1.20 | 1.8571e+07 | 2.07 |
| D (-3) | 300 | 20 | 0.2352 | 0.0 | 0.00 | 1.20667e+07 | 1.93 |
| D (-3) | 300 | 25 | 11.8416 | 0.0 | 11.70 | 9.57991e+06 | 1.89 |
| D (-3) | 300 | 60 | 22.2333 | 0.0 | 21.90 | 4.65962e+06 | 2.53 |
| D (-3) | 300 | 80 | 27.7302 | 0.0 | 27.20 | 4.05092e+06 | null |
| D (-3) | 300 | 120 | 38.0575 | 0.0 | 37.20 | 3.3792e+06 | null |
| C (-2) | 300 | 13 | 1.6611 | 0.0 | 1.50 | 2.21254e+07 | 1.92 |
| C (-2) | 300 | 15 | 1.4539 | 0.0 | 1.30 | 1.79243e+07 | 1.97 |
| C (-2) | 300 | 20 | 3.3588 | 0.0 | 3.20 | 1.24977e+07 | 2.03 |
| C (-2) | 300 | 25 | 11.9659 | 0.0 | 11.80 | 9.65479e+06 | 1.92 |
| C (-2) | 300 | 30 | 15.0899 | 0.0 | 14.90 | 7.97357e+06 | 1.84 |
| C (-2) | 300 | 60 | 22.9757 | 0.0 | 22.60 | 4.79623e+06 | 1.77 |
| C (-2) | 300 | 80 | 33.4911 | 0.0 | 33.00 | 4.10502e+06 | null |
| C (-2) | 300 | 120 | 50.3346 | 0.0 | 49.60 | 3.52754e+06 | null |
| B (-1) | 300 | 13 | 2.3894 | 0.0 | 2.30 | 2.19474e+07 | 1.90 |
| B (-1) | 300 | 15 | 3.9165 | 0.0 | 3.80 | 1.78794e+07 | 1.99 |
| B (-1) | 300 | 25 | 6.6392 | 0.0 | 6.50 | 9.62792e+06 | 1.92 |
| B (-1) | 300 | 30 | 14.2110 | 0.0 | 14.00 | 8.20453e+06 | 2.05 |
| B (-1) | 300 | 60 | 32.7377 | 0.0 | 32.50 | 4.77463e+06 | 2.23 |
| B (-1) | 300 | 80 | 48.2494 | 0.0 | 47.90 | 4.00094e+06 | 2.20 |
| B (-1) | 300 | 120 | 79.9214 | 0.0 | 79.50 | 3.38673e+06 | 1.71 |
| A (+0) | 300 | 15 | 8.8345 | 0.0 | 8.77 | 1.83069e+07 | 2.77 |
| A (+0) | 300 | 20 | 11.8982 | 0.0 | 11.82 | 1.25991e+07 | 1.94 |
| A (+0) | 300 | 30 | 18.9072 | 0.0 | 18.80 | 8.48107e+06 | 1.65 |
| A (+0) | 300 | 40 | 26.3268 | 0.0 | 26.20 | 6.61381e+06 | 1.91 |

## Project disposition

This packet establishes the anomaly and the exact source questions; it does
not establish a correction. The four Boccioli--Roberti admission blockers are
unchanged, physical-node inventory remains empty, and canonical conversion,
runtime deposition, production, and publication remain disabled.

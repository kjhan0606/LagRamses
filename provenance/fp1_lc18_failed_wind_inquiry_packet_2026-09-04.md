# LC18 failed-wind source inquiry packet

Date: 2026-09-04; updated 2026-09-05 after operator-authorized source review
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
- The independent LC18 CDS table5 endpoint gives a parsed positive cumulative
  wind for 53 failed models and a parsed exact-zero endpoint for three; the
  successful controls contain 48 parsed-positive and 4 parsed exact-zero
  endpoints. These are exact differences of the parsed masses; no pipeline
  rounding is performed and no physical zero or upper limit is inferred.
  F6.2 specifies a 0.01 M_sun printed decimal step, NOT an attested physical
  precision or rounding rule. The three failed parsed exact-zero endpoints remain inside the
  unresolved BR26 zero-Wind release anomaly and do not define or resolve it.
  The summary-wind differences exceed both the printed-format half-bin and
  the three-significant-digit sensitivity half-bin for all 108 models,
  including all 52 successful controls. The largest absolute difference is
  1.5476 M_sun. Thus CDS is corroborating evidence that the failed stars
  underwent mass loss, but it is not a numerically identical replacement for
  the release summary.
- CDS table5 yields 845 unique phase rows after 19 exactly duplicated rows are
  collapsed. Each model has 3--8 unique phases, a PSN endpoint, strictly
  increasing cumulative age, and non-increasing total mass.
- CDS table7 supplies pre-SN structure for 96 models. The other 12 remain
  explicit nulls. Its binding energy is not interpreted as injected explosion
  energy.
- Original LC18 integrated winds are available and positive for all 108
  coordinates, including all 56 BR26 failed models: sum table9 for <=25 M_sun
  and table8 set-R wind-only for >25 M_sun (LC18 section 5). They have NOT
  replaced BR26 winds. For (M, [Fe/H], v) = (20, 0, 0), (120, -1, 300),
  (120, -3, 0), the original sums are respectively 12.459052587158395,
  79.50545602731917, 0.34122894710262286 M_sun; BR26 Wind is zero in each.

## Precision correction and unresolved mass-definition question

On 2026-09-05 the publisher's [Table 5 MRT](https://content.cld.iop.org/journals/0067-0049/237/1/13/revision1/apjsaacb24t5_mrt.txt)
was retrieved read-only. SHA256 is
`2d1787477d6c8160967c8b4289e4ce476cde88e540e9469a80553f67237e19a3`.
All 864 data rows have identical whitespace-separated tokens to the staged
CDS table. The publisher also specifies F6.2: no higher-precision replacement
or explicit nearest-rounding rule was found in this source. The CDS ReadMe
retrieved on the same date matches the staged SHA256
`674df04c58b87111c8276a3e1ad827bf0310194640143a4371740e8f97fe8682`.

All 864 total masses lie on a three-significant-digit grid: 60 rows at step
0.01, 776 at 0.1, and 28 at 1 M_sun. This is an observation, not proof of the
source rounding algorithm. Retain two descriptive comparisons:

| Comparison, 108 distinct coordinates | Above 0.005 M_sun printed-format half-bin | Above mass-dependent three-digit sensitivity half-bin |
|---|---:|---:|
| Original LC18 integrated wind vs initial-minus-PSN | 89 | 2 |
| BR26 summary wind vs initial-minus-PSN | 108 | 108 |

The sensitivity assumes nearest rounding with half-bin
`0.5 * 10**(floor(log10(M_PSN)) - 2)`; it is NOT an approved tolerance and
does not include an attested isotope-sum uncertainty budget. Remaining
original-wind-minus-endpoint outliers are:

- (M=30, [Fe/H]=-1, v=300): +0.052420896102670866 M_sun (half-bin 0.05).
- (M=120, [Fe/H]=0, v=300): -0.08899943941717936 M_sun (half-bin 0.05).

In 106/108 models, `(BR26 summary wind - original LC18 wind) /
(M_initial - original LC18 wind)` lies in 0.009--0.011; the median is
0.010060988970999758. [LC18 section 2](https://openaccess.inaf.it/bitstream/20.500.12386/28562/2/Limongi_2018_ApJS_237_13.pdf)
describes a late-stage inner-zone boundary at 99% of the total mass. An
atmosphere/profile-boundary convention is therefore a hypothesis to ask
about, not a confirmed cause or correction. The two fraction outliers are
(40, -1, 0): 0.014255817194057633 and
(120, -3, 150): 0.0046494621372653476.

The 2026-09-04 Opus parsed-zero PASS remains a historical implementation
audit, not approval of physical source closure. Its 0.01/0.005 M_sun physical
precision wording is superseded by this review and the live reports. Do not
rewrite the historical auditor's verdict or turn the sensitivity into PASS.

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
- Limongi--Chieffi CDS `table8.dat`: SHA256
  `bdfc317ca12a377f545964424dba4a666eab964292d0ebfa8b8c9641f073f218`.
- Limongi--Chieffi CDS `table9.dat`: SHA256
  `64b4bbe536514dfd164bb0e971d2d9163043a767a2ec2de552c76bf131b5088d`.
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
3. Is the approximately 1% excess over original LC18 integrated winds related
   to excluding an outer atmosphere from the BR26 pre-SN profile mass? Which
   mass boundary and stellar-model revision define `M_presn` and `M_wind`?
   Is atmosphere mass intentionally assigned to Wind, and, if so, at what
   age and composition? Please also explain the two fraction outliers above.
4. Can the per-model explosion energies and their semantics (diagnostic,
   asymptotic, or injected), terminal mass cuts, fallback/remnant masses, and
   the BHSN exception list be released in machine-readable form?
5. For failed models, is a remnant mass or fallback/direct-collapse convention
   available that preserves the pre-terminal wind separately from terminal
   ejecta?
6. For the LC18 Table 5 authors: what is the actual mass precision and rounding
   rule before F6.2 formatting? Can higher-precision PSN masses be supplied?
   All entries lie on a three-digit grid, but this alone does not establish
   nearest rounding. Can the two original-wind closure outliers above be
   explained by source precision, endpoint definitions, or model revisions?
7. Can original LC18 winds be combined with BR26 terminal ejecta without an
   atmosphere or remnant correction? For every one of the 52 successful
   models, `M_initial - original LC18 wind - BR26 Post element sum - BR26
   mass cut` is positive: 0.0501418214872682--0.6072574340308559 M_sun.
   A naive wind replacement therefore leaves missing mass rather than a
   self-consistent combined package.

The two author teams have different ownership: questions 1--5 and 7 concern
the BR26 release; question 6 concerns the original LC18 tables. This is a
prepared inquiry, not an email or an issue submission. No external message
has been sent and no reply is implied.

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

Implementation verification, 2026-09-05 (existing focused checks only):

- `python simulation/snrt/tests/g2_limongi_phase_mass_history.py`:
  `G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK`.
- `python simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`:
  `FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK`.
- `python simulation/snrt/tests/fp1_physical_package_admission.py`:
  `FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK`.
- Existing phase-history and cross-check JSON reports regenerated. The
  cross-check now binds the current metallicity-domain contract hash. Source
  checksums and raw phase histories are unchanged; the cross-check test
  verifies no mutation of source tables or config/data during evaluation.
  These are software/evidence checks, NOT a physical-package PASS.

Original LC18 winds are a viable **integrated-wind candidate**, not an
approved drop-in correction or a time-resolved feedback source. Use the same
wind supplier for successful and failed nodes if this candidate is selected;
do not switch suppliers only where BR26 Wind is zero. Preserve the original
mass/[Fe/H]/rotation coordinates, absolute metallicities and isotopic epoch.
The grid covers 13--120 M_sun, [Fe/H]=0,-1,-2,-3 and v=0,150,300 km/s;
it does not cover the project's entire Z=0--0.139 target or select a rotation
distribution. Table8's >25 M_sun wind-only convention does not impose LC18's
terminal fate prescription on a different explosion model.

Next physical work, within the existing source-admission bundle: obtain or
approve age-dependent composition and wind-speed/energy prescriptions;
resolve the combined wind+terminal+remnant mass budget; select the rotation
population and low/high-Z supplements. A mass-only phase history cannot
supply isotope timing or kinetic energy, and dust/RT sources require their
own composition/spectral mapping. No new converter, gate or runtime path is
justified until these physical inputs are supplied. Do not renormalize wind,
allocate the unexplained 1% by hand, or invent a failed-model remnant mass.

This packet establishes the anomaly and the exact source questions; it does
not establish a correction. The four Boccioli--Roberti admission blockers are
unchanged, physical-node inventory remains empty, and canonical conversion,
runtime deposition, production, and publication remain disabled.

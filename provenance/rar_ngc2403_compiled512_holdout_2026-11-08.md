# Compiled NGC 2403 force fields and radial holdouts

Grammar CPU analysis `1279415` completed in 1m50s using analysis source
`f00f5f2a6ca43ea1222953565e9b448b29adaeb7`. It reads the existing two
512^3 snapshots from `451793`, verified by Syntax postprocessor `456010`.
Snapshot digests are checked again during reading. No simulation was rerun.

The 64 selected SPARC radii span 1.06--20.87 kpc. The source uses fixed
stellar M/L=0.5 and the previously prescribed THINGS-like gas distribution.
Each of four contiguous blocks excludes 16 radii from the NFW fit; its two
parameters are fitted on the remaining 48. QUMOND uses the fixed input
parameters. All NFW optima lie inside the declared parameter grid.

| Diagnostic over excluded radii | Fixed QUMOND | Analytic NFW |
|---|---:|---:|
| Velocity RMS residual (km/s) | 13.333 | 5.160 |
| Mean velocity residual (km/s) | -11.230 | +0.942 |
| Median absolute velocity residual (km/s) | 12.614 | 3.282 |
| RMS log10(g_model/g_observed) (dex) | 0.09912 | 0.04135 |
| Median log10(g_model/g_observed) (dex) | -0.11039 | -0.00530 |
| Diagonal chi-squared per point | 326.66 | 33.58 |

The NFW prediction has a smaller RMS residual in the aggregate and in
three of four excluded blocks. QUMOND has the smaller RMS residual in the
2.68--4.97 kpc block (6.329 versus 7.495 km/s). The fitted NFW values across
folds are V200=97.75, 98.75, 113.75, 98.5 km/s and c200=12.4, 12.15, 9.0,
12.1. Exact per-radius predictions, folds and source digests are retained
in the CSV and JSON results.

These residuals describe one prescribed galaxy source. Population intrinsic
RAR scatter is unmeasured. The compiled baryonic acceleration differs from
the SPARC mass reconstruction by a median -0.02989 dex. NFW has two fitted
parameters per fold while QUMOND remains fixed. Most training radii are
shared between folds, and diagonal velocity errors omit radial covariance,
distance/inclination, stellar-population and gas uncertainties. The large
diagonal chi-squared values remain explicit; these scores do not supply a
calibrated likelihood or a cosmological model-selection significance.

The NFW contribution is an isolated spherical analytic acceleration added
to the periodic baryonic force. A live CDM halo and a self-consistent
cosmological Lambda-CDM/QUMOND comparison remain subsequent work. Periodic
box effects, source reconstruction and existing AMR self-force limitations
also remain separate from the exact discrete force checks.

Original output and checksums are preserved in
`/gpfs/kjhan/RAR_QUMOND/analysis/ngc2403_compiled512_cdm_nfw_holdout_1279415`.
Project copies are in `reports/ngc2403_compiled512_cdm_nfw_holdout_1279415/`.
Compact results are tracked in
`provenance/results/rar_ngc2403_holdout_1279415/`.

Visual inspection found the original right-panel title and legend touched
the plotted curve. The renderer moves them into the empty upper-right
region, retaining the original figure. The replacement figure is generated
directly from the saved CSV, without reading the force files or fitting again.
Grammar render `1279417` completed and the replacement PDF/PNG passed visual
inspection. The CSV is byte-identical to the analysis input. The figures are
tracked as `provenance/figures/rar_ngc2403_compiled512_holdout.pdf` and `.png`;
their digests and render source revision are recorded in the result bundle.

## Why agreement with the reference does not reproduce the observed curve

Grammar postprocessing job `1279459` completed in 9 s using source
`bec88445555928415db829d5f88ad942f5adbf11`. It reads the checksum-pinned
64-radius CSV from `1279415`; no force calculation or fit was repeated.
Let `f(g)` be the plotted empirical reference, `g_b` the computed baryonic
field, and `g_s` the SPARC baryonic acceleration. At each radius,

`log10(g_Q/g_obs) = log10(g_Q/f(g_b)) + log10(f(g_b)/f(g_s)) + log10(f(g_s)/g_obs)`.

| Term | Mean over all 64 radii (dex) | Mean over 13.42--20.87 kpc (dex) |
|---|---:|---:|
| QUMOND minus reference at computed baryonic field | -0.019146 | -0.001266 |
| Reference shift between the two baryonic inputs | -0.033875 | -0.072037 |
| Reference at SPARC baryonic acceleration minus observed | -0.032320 | -0.058308 |
| QUMOND minus observed | -0.085341 | -0.131611 |

The identity closes to `7.64e-17 dex`. These equal-weight arithmetic means
add; the individual medians do not generally add. The decomposition explains
the apparent coexistence of agreement with the reference and an offset from
NGC 2403. It is algebraic accounting, not a causal attribution or an error
budget. In particular, the second term is not the response of a new QUMOND
solution with a SPARC-matched source, and the first term does not isolate
boundary or discretization errors. The observed points belong to this one
galaxy. The result does not measure population intrinsic scatter.

All per-radius terms, statistics, immutable input, analysis source and
checksums are retained in `results/rar_ngc2403_residual_split_1279459/`.

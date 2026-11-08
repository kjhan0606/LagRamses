# Recovery of existing NGC 2403 force analysis

Syntax job `451793` allocated one A100 GPU. Both compiled, uniform
512-cubed Newtonian/QUMOND force cases completed, including FFT and QUMOND
diagnostic markers, with no solver-error markers. Each stored snapshot is
28,684,370,544 bytes. The wrapper failed later in `analyze.py` because
`torch.quantile(errors, .9)` rejected the large input tensor.

Analysis commit `c9c7938f5f6471c9257ee44cd3b1089031a69001` replaces this
statistic with an exact NumPy quantile using linear interpolation. The
original lower-median statistic, force/source equations, mask, and `1e-6`
acceptance thresholds are unchanged. Snapshot hashes are streamed instead
of allocating a whole-file byte string. Unused arrays are released before
the next stage. An optional exclusive output path preserves existing data.

Recovery script `provenance/tools/recover_ngc2403_512_analysis.sbatch` reads
the original run directory without modifying it and extracts analysis code
from the explicit Git revision. Syntax job `455930` requests one A100,
2 CPUs and 128 GiB for GPU postprocessing only. New output is
`/gpfs/kjhan/RAR_QUMOND/analysis/ngc2403_512_recovery_455930`. At recording,
the job is pending resources and has no GPU allocation. Actual recovered
force/source results are not yet established.

The force comparison uses the same prescribed baryon source in Newtonian
and QUMOND gravity. It does not include a CDM halo, a native AMR particle
calculation, or cosmological evolution. The full original run and failed
postprocessing logs are retained.

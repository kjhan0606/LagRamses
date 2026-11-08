# NGC 2403 graph-derived gas source at matched force resolution

The preceding goal turn explained existing RAR offsets but did not advance
state. This turn recovers the already prepared FR02 Figure 9 calculation and
then compares it with the existing THINGS-like field at the same 512-cubed
resolution. It retains the first Lambda-CDM/QUMOND science direction; the
analytic NFW comparator remains a prescribed halo, not cosmological evolution.
The full hydro, live-CDM, BIMOND/lagCAMB and FDM objectives remain open.

Syntax accounting reports old job 454073 FAILED, exit 1:0, elapsed zero,
with empty output/error files and no run directory. A metadata-only shell
probe reproduces failure to reach the end of strict oneAPI startup; staging
the vendor setup before enabling errexit/nounset returns zero and locates
mpirun. The replacement follows the startup sequence of the completed force
calculations in 451793. That job's later analysis failure was separately
recovered by 456010; the original wrapper remains FAILED, not relabelled.

This is existing approved source/force work. No Fortran physics is changed.
The pinned binary is source 65eeb24ddb43e22e6b6acc48b6758c386d506413,
SHA256 7baf42cc884f33e9babb09b21aca5c46c380b3c32186cd562d20971c6414061c.
FR02 preparation 1274676 has a passing original manifest. Its gas is the
44 digitized markers, distance rescaled to 3.16 Mpc and helium multiplied
by 1.33, with a declared outer taper. Stars have M/L=0.5 and the corrected
0.24143-kpc exponential height. The source is a candidate, not identified
as the exact SPARC mass-model input. Heights, geometry, taper and periodic
boundary remain assumptions.

The new runner freezes namelist templates and the independent oracle from
an explicit analysis commit. It admits only 128 or 512 cubed uniform grids,
first the short 128 source-ingestion check and then 512 for a same-resolution
comparison to 451793. The oracle is the validated 512 operator reconstruction
with FR02 contract and level coverage checks. It uses historical source-only
telemetry explicitly, never fabricating a second-solve diagnostic. Poisson
source RMS and force median/p90 thresholds remain 1e-6; exact discrete
operator agreement is not continuum or isolated-boundary accuracy.

Both cases are static baryon-source Newtonian/QUMOND solves on Syntax Slurm;
the Newtonian case supplies gbar and is not Lambda-CDM. GPU FFT is enabled;
hydro state is a uniform static tracer with feedback/cooling inactive for this
force experiment, not the eventual evolving hydro programme. There are no
particles, sinks, dust or ray tracing. Per case nstepmax=1, noutput=1,
aout=1.1, tout=1e100, foutput=1 and fbackup=1e9 intentionally yield one dump.
At 128, measured comparable dumps are 448291224 bytes; reserve 2 GiB total
and require 4 GiB free. At 512, 28684370544 bytes per dump, reserve 100 GiB
and require 150 GiB free. Current shared free space exceeds 140900 GiB.
Effective absolute nml paths and runtime storage are recorded by the runner.
Original failed jobs and all outputs are retained. Grammar Slurm handles
CPU observational comparisons. No computation was run on a login host.

Status at preparation: no replacement job submitted yet. Updated job IDs,
results and limits will be appended after execution.


## Successful recovery and high-resolution execution

The frozen analysis/runner revision is 371307f73afcabe0a9f1ce44764599dd5ef70f38.
Syntax job 456155 completed the 128-cubed paired run in 49 s, exit 0:0.
Both independent fields pass: Newtonian median/p90 relative force error
2.78744e-15 / 6.20767e-15, QUMOND 4.15684e-15 / 7.72007e-15. Newtonian
source Laplacian RMS relative error is 3.76977e-14. All 128^3 cells were
updated in both phantom-source calls, with no skipped/missing neighbors.
The force error distribution uses 376872 cells in the stated radial shell;
these metrics are discrete operator checks, not continuum estimates.

The pending short job's partition list was expanded to a100,a100_pcie,
which share the admitted A100 architecture and permit this account. Its ID
was preserved. Dependent 512-cubed job 456162 started on syn102 after 456155
succeeded. It retains the same physical input tables and compiled binary.
No failed prerequisite was bypassed. Scheduler forecasts changed as resources
became available; observed execution takes precedence over earlier estimates.

The Grammar comparison program is committed in 48ec3e6 and has not yet run.
It fixes the 64 observed radii, source-free parameters and sampling method,
reuses the immutable 1279415 baseline CSV, and repeats the same four radial
NFW holdouts for FR02. The same stellar table is byte-identical in both
sources. Both radial profile and taper change in the gas intervention; the
result must not be presented as a pure mass-normalization effect. The
THINGS branch already applies the 1.33/1.36 helium rescale, so there is no
separate helium-factor change to attribute. New plots compare at the same
radii to avoid confusing different baryonic x coordinates.

Compact checkpoint: results/rar_ngc2403_fr02_recovery_20261108/.

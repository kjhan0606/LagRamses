# lagRamses operational tools

These scripts configure production and benchmark calculations for the
sources selected first by `bin/Makefile` through
`VPATH = ../patch/lagRamses:...`.

## `dmo_benchmark_setup.py`

Builds a matched-phase DMO campaign for the implemented DE and modified
gravity models. The default `model` IC mode imports transfer functions from
`/home/kjhan/BACKUP/CAMB/CAMB` for every parameter-matched lagCAMB model.
The generator uses `/home/kjhan/BACKUP/LagMUSIC/music/build/MUSIC` and fixes
the absolute amplitude through `force_pnorm`. The normalization is derived
from the lagCAMB linear spectrum at the starting redshift and retains the
common primordial amplitude rather than imposing a common `sigma_8`.
For DMO 2LPT ICs, the generator also derives `vfact_scale` from the
model-to-LCDM ratio of lagCAMB's total velocity-to-density transfer. This
corrects the initial velocity without changing the density or displacement
field. Generation stops if that ratio varies by more than `1e-4` over
`0.001 <= k <= 1 h/Mpc`; such a model needs a scale-dependent velocity
kernel rather than one scalar correction.

```bash
python3 patch/lagRamses/aux/dmo_benchmark_setup.py \
  --outdir /gpfs/kjhan/Hydro/DE_nonstd/DMO_bench_v1 \
  --make-ics
```

- Default: `500 Mpc/h`, `256^3`, `levelmax=14`, 2LPT at `z=49`
- Default validation suite: LCDM, F5, F6, N1, N5, and Symmetron A
- Uses one random seed and identical phases for all model-specific ICs
- Supports the old shared LCDM protocol through `--ic-mode shared`
- Enables `match_aout` so every model is dumped at the same expansion factor
- Exposes the coarse cosmological step through `--aexp-step-limit`
  (default `0.1`; reduce it for temporal-convergence tests)
- Records every exact transfer match and LCDM fallback in `campaign.json`
- Uses `mpirun -np 1` for IC generation by default

## `plot_dmo_pk_ratios.py`

Plots matched-phase power-spectrum ratios together with the parameter-matched
lagCAMB linear prediction, and writes CSV and JSON validation products. By
default it interpolates in `log(P)` and `log(a)` when legacy outputs do not
share an exact epoch; use `--nearest` for campaigns generated with
`match_aout`.

```bash
python3 patch/lagRamses/aux/plot_dmo_pk_ratios.py CAMPAIGN \
  --models f5 f6 n1 n5 sym_a --nearest
```

## Uniform-resolution convergence

Generate consistent `64^3`, `128^3`, and `256^3` campaigns:

```bash
python3 patch/lagRamses/aux/dmo_resolution_ladder.py \
  --root /gpfs/kjhan/Hydro/DE_nonstd/DMO_resolution_v1
```

The ladder also supports targeted extensions. For example, a common-anchor
L8--L9 F5 convergence pair with only the two required analysis dumps is:

```bash
python3 patch/lagRamses/aux/dmo_resolution_ladder.py \
  --root /gpfs/kjhan/Hydro/DE_nonstd/DMO_f5_anchor9 \
  --levels 8 9 --models lcdm f5 --phase-anchor-level 9 \
  --scalar-eps 1e-5 --output-redshifts 5 0
```

Add `--make-ics` for resolutions that are safe to generate on the current
host. By default the ladder places the random seed at the finest requested
level and obtains every coarser white-noise field by LagMUSIC restriction.
This common phase anchor is required for a nonlinear convergence test:
independently putting the same integer seed at L6, L7, and L8 generates
different Fourier realizations. Override the default only with
`--phase-anchor-level LEVEL`, where `LEVEL` must not be below the finest
requested resolution. LagMUSIC's distributed restriction RNG is not yet
decomposition independent, so the ladder uses one MUSIC rank at every level
below the phase anchor; the anchor-level IC may still use its parallel
resource preset. `check_grafic_phase_matching.py` directly verifies the
generated `ic_deltab` low-k Fourier phases after correcting the deterministic
cell-centre offset between grid resolutions. Each campaign contains
`submit_all.sh` for concurrent Slurm submission.
After runs complete, compare every resolution with lagCAMB at `z=0`:

```bash
python3 patch/lagRamses/aux/measure_dmo_pk.py \
  CAMPAIGN/L8_256/lcdm/output_00008 \
  --nmesh 256 --kmax 0.5

python3 patch/lagRamses/aux/plot_dmo_resolution_convergence.py \
  /gpfs/kjhan/Hydro/DE_nonstd/DMO_resolution_v1
```

`measure_dmo_pk.py` is the precision DMO estimator: it deposits particles on
a common CIC mesh, interlaces a half-cell-shifted mesh, deconvolves the CIC
window, and writes both raw and mass-weighted Poisson-shot-subtracted spectra
on identical integer-|k| shells at every simulation resolution. Validation
uses the raw spectrum by default because these ICs are perturbed particle
lattices, not Poisson samples; `--shot-noise poisson` is an optional
diagnostic only. The cheap on-the-fly RAMSES diagnostic is NGP on the base
grid and should only be selected with `--pk-estimator runtime-ngp`; it is not
suitable for a 0.1% convergence claim.

Only resolutions for which every model has a measured spectrum at the
requested scale factor are included. The JSON report gives both the full
`k <= 0.5 h Mpc^-1` residual and a separate large-scale residual (default
`k <= 0.2 h Mpc^-1`), together with the worst shell's `k` and signed
residual for each theory and adjacent-resolution comparison. Resolution-pair
records also report the number of failing shells and the largest contiguous
`kmax` for which every shell from the fundamental mode satisfies the target;
change the large-scale limit with `--large-scale-kmax`. The default acceptance
limit is 0.1% (`--residual-target 0.001`). Use
`--require-resolution-pass` and/or `--require-theory-pass` when a failed
P(k) check must also produce a nonzero process exit status.

## Hu-Sawicki f(R) solver performance

The f(R) implementation is resolved by `bin/Makefile` through VPATH from
`patch/cuRamses/force_fine.kjhan.f90`.  Its uniform-domain fast path:

- rescales the previous scalar field by the change in the homogeneous
  `fR_bar` before using it as the next Newton initial guess;
- performs eight spectral Newton/Helmholtz corrections before local
  red-black Newton-GS relaxation;
- reads existing same-level neighbors directly and invokes the Morton/CIC
  coarse-fine lookup only when a same-level cell is absent;
- uses static OpenMP scheduling for equal-size uniform-grid chunks; and
- evaluates the common Hu-Sawicki `n=0` and `n=1` inversions with division
  and `sqrt`, retaining the general real-power expression for other `n`.

Solver completion and performance are not validation criteria by themselves.
Any optimization must be checked against an unmodified matched-phase run
with `measure_dmo_pk.py`; the maximum shell residual, not an RMS summary,
must remain below the campaign's 0.1% acceptance threshold.

## nDGP solver performance

The uniform-domain nDGP operator split needs twelve diagonal scalar samples
per cell to construct the mixed Hessian. The FFT path caches the twelve
edge-neighbour octs once per vector chunk and then reads the corresponding
cells directly. The general Morton/parent-CIC sampler remains the fallback
when an AMR neighbour is genuinely absent. The fifth-force gradient likewise
reuses the six face neighbours already gathered by the routine.

On the `64^3` N1 z=0 regression (4 MPI ranks x 4 OpenMP threads), these
changes reduced the RAMSES total timer from 19.253 s to 8.052 s. The common
CIC P(k) file was byte-identical. At `128^3`, the maximum CIC P(k) change was
`1.31e-5%`, well below the 0.1% acceptance threshold.

The remaining repeated work is cached at the scalar-solve scope. Six face
and twelve edge grid neighbours are now constructed once at the beginning
of each nDGP/Galileon solve and reused by all nonlinear corrections. The
cache is conservatively rebuilt at the next coarse step so a load balance
cannot leave stale active-grid indices.

The shared spectral Helmholtz correction retains its FFTW plans and
real/complex work arrays and precomputes the one-dimensional discrete
Laplacian eigenvalues. Multi-rank uniform grids now use a separate FFTW-MPI
slab path through `512^3`: owned RAMSES cells are exchanged with the FFT
slabs, transformed with the same seven-point eigenvalues, and returned
without replicating the full field on every rank.

Controlled same-node baseline/candidate tests gave:

- N1 `64^3` z=0: total 9.055 -> 6.671 s and nDGP 5.294 -> 2.960 s;
- N1 `128^3` z=0: total 165.75 -> 125.04 s and nDGP
  110.69 -> 69.96 s; and
- N1 `256^3`, first 20 steps (8 MPI x 4 OpenMP): total
  127.46 -> 101.06 s and nDGP 103.85 -> 76.51 s.

For the full `64^3` and `128^3` runs, both the native spectra and common
interlaced `256^3` CIC spectra were identical at z=5 and z=0. The L8 short
test reproduced every timestep, nonlinear iteration count, and reported
residual exactly. Full F5 and Symmetron-A `64^3` regressions also produced
identical native and common-CIC spectra at both redshifts; their runtimes
were unchanged within run-to-run noise.

The distributed-path F5 `128^3` same-node A/B retained exactly the same
210,464 Newton-GS sweeps and byte-identical z=5/z=0 common-CIC spectra while
reducing total time from 762.35 to 519.99 seconds and `fR-solve` from
727.22 to 487.04 seconds. The largest ID-matched z=0 particle-position
difference was `3.21e-14` box units. Separate N1 and Symmetron-A `64^3`
regressions were also byte-identical in common-CIC P(k), covering both the
zero-mass Poisson and massive Helmholtz forms.

## Symmetron solver performance

After symmetry breaking, the strict `1e-6` field residual can require
hundreds of Newton-GS sweeps. Each sweep already gathers the six same-level
neighbour grids, so the uniform-domain path now reads those cells directly
and calls the Morton/parent-CIC sampler only at a genuine AMR boundary. The
fifth-force gradient uses the same direct-neighbour rule.

On a controlled same-node Symmetron-A `64^3` z=0 A/B run (4 MPI ranks x
4 OpenMP threads), the total timer fell from 68.89 to 17.03 seconds and the
Symmetron timer from 65.01 to 13.12 seconds. The z=5 and z=0 common
interlaced `256^3` CIC spectra were byte-identical.

The face-grid topology is now also built once per scalar solve and reused by
all f(R)/Symmetron spectral corrections and Newton-GS sweeps. A second
same-node `64^3` A/B test reduced F5 from 101.49 to 22.54 seconds
(`fR-solve`: 97.59 to 18.66 seconds) and the already optimized Symmetron run
from 17.15 to 6.70 seconds (`symmetron`: 13.24 to 2.84 seconds). For both
models, the z=5 and z=0 common CIC spectra remained byte-identical.

At `128^3`, the fully optimized Symmetron-A run completed in 97.15 seconds
instead of 1749.50 seconds; its Symmetron timer fell from 1718.34 to
65.19 seconds. The exact z=5 and z=0 common CIC files were again
byte-identical.

The corresponding F5 `128^3` total timer fell from 2850.22 to 727.31
seconds and `fR-solve` from 2815.58 to 694.45 seconds. Its exact z=5 and
z=0 common CIC files were also byte-identical.

### F5 residual-tolerance audit

The production validation runs retain the strict relative field residual
`fR_eps=1e-6`.  Controlled same-node tests quantify the possible speed and
power-spectrum cost of relaxing that stopping criterion; they do not redefine
the 0.1% power-spectrum validation target.

At `64^3`, changing `fR_eps` from `1e-6` to `1e-5` reduced the total timer
from 43.57 to 28.18 seconds and the F5 timer from 39.29 to 23.68 seconds.
The z=0 common-`256^3` CIC spectrum changed by at most `0.00398%`.  A more
aggressive `1e-4` tolerance reduced the total timer to 15.12 seconds but
changed z=0 P(k) by as much as `0.03614%`, consuming too much of the 0.1%
error budget to adopt without finer-resolution tests.

At `128^3` on one node with 8 MPI ranks and 4 OpenMP threads per rank,
`1e-6 -> 1e-5` reduced the total timer from 744.06 to 386.98 seconds and
the F5 timer from 710.66 to 353.71 seconds.  Newton-GS work fell from
210,464 to 95,382 sweeps.  The maximum raw-spectrum change was negligible
at z=5 and `0.01108%` at z=0 for `k <= 0.5 h Mpc^-1`.  This makes `1e-5`
a promising performance setting.

The full `256^3` A/B test subsequently passed: `1e-5` reduced total time
from 13,831.17 to 5,425.03 seconds (2.55x) and Newton-GS work from
1,043,465 to 316,122 sweeps. Relative to the strict `1e-6` run, the z=0
raw P(k) changed by at most `0.03056%` over both `k <= 0.2` and
`k <= 0.5 h Mpc^-1`; z=5 changed by less than `6e-8%`. Thus `1e-5` is
approved for the targeted L8--L9 extension, while the strict run remains
the numerical reference and the spatial-resolution residual is reported
separately.

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

Add `--make-ics` for resolutions that are safe to generate on the current
host. By default the ladder places the random seed at the finest requested
level and obtains every coarser white-noise field by LagMUSIC restriction.
This common phase anchor is required for a nonlinear convergence test:
independently putting the same integer seed at L6, L7, and L8 generates
different Fourier realizations. Override the default only with
`--phase-anchor-level LEVEL`, where `LEVEL` must not be below the finest
requested resolution. Each campaign contains `submit_all.sh` for concurrent
Slurm submission.
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
`k <= 0.2 h Mpc^-1`); change the latter with `--large-scale-kmax`. The
default acceptance limit is 0.1% (`--residual-target 0.001`). Use
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

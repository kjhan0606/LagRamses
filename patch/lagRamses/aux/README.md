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
host. Each campaign contains `submit_all.sh` for concurrent Slurm submission.
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

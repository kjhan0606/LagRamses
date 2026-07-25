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

```bash
python3 patch/lagRamses/aux/dmo_benchmark_setup.py \
  --outdir /gpfs/kjhan/Hydro/DE_nonstd/DMO_bench_v1 \
  --make-ics
```

- Default: `500 Mpc/h`, `256^3`, `levelmax=14`, 2LPT at `z=49`
- Uses one random seed and identical phases for all model-specific ICs
- Supports the old shared LCDM protocol through `--ic-mode shared`
- Enables `match_aout` so every model is dumped at the same expansion factor
- Records every exact transfer match and LCDM fallback in `campaign.json`
- Uses `mpirun -np 1` for IC generation by default

## `plot_dmo_pk_ratios.py`

Plots matched-phase power-spectrum ratios and writes CSV and JSON validation
products. By default it interpolates in `log(P)` and `log(a)` when legacy
outputs do not share an exact epoch; use `--nearest` for campaigns generated
with `match_aout`.

```bash
python3 patch/lagRamses/aux/plot_dmo_pk_ratios.py CAMPAIGN \
  --models f5 f6 n1 n5 sym_a --nearest
```

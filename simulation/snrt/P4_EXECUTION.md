# P4 execution contract

## Thermal model

P4 accepts only the v3 UVB-free metal atlas described in
[`P4_THERMAL_ATLAS.md`](P4_THERMAL_ATLAS.md). The old equilibrium
`build_grackle_atlas.py` workflow is retired from the execution path because it
mixed equilibrium primordial cooling, HM2012 UVB terms, metal cooling, and
separately applied SNRT photoheating.

AMR covering-grid padding cells with non-positive density or pressure are not
treated as resolved gas. Provided the region also contains valid hydro cells,
P4 replaces only those padding cells with the atlas minimum `n_H` and minimum
`T`; this is a declared near-vacuum array floor. A target whose physical cells
lack usable pressure or temperature is an ingestion validation artifact, not a
production thermal state.

## Target selection

1. Read documented sink ID, mass, and position from `sink_00016.info`.
2. Use the most massive non-boundary sink only to center a 64-cubed local hydro probe.
3. Select the 32-cubed subvolume with maximal mean hydrogen density.
4. Restage that subvolume from the AMR/hydro files into canonical HDF5.

The density criterion, not sink mass, determines the target. The sink only
avoids an impractical whole-box AMR scan.

## Source ledger

`sink_00016.csv` has two columns not defined by `sink_00016.info`. P4 does not
identify either as an accretion rate and emits no AGN photons from it. The
output roster lists every sink inside the selected volume for a later audited
mapping. No star-particle source data are present in this P4 input.

## Run

```bash
sbatch build_p6_thermal_atlas.sbatch

./.venv/bin/python tools/p4_stage_high_density.py \
  --info /gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/info_00016.txt \
  --sink-info /gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/sink_00016.info \
  --thermal-atlas data/production_metal_thermal_atlas_v2.h5 \
  --scale-factor 0.208497764676753 \
  --output data/p4_high_density_rt_input.h5 \
  --scratch .p4_scratch \
  --manifest data/p4_high_density_manifest.json \
  --sink-roster data/p4_high_density_sinks.csv
```

The canonical HDF5 is an ingestion artifact, not a dual-AGN science result. An
AGN RT run still requires an audited instantaneous-accretion-rate column and a
declared AGN SED/group conversion.

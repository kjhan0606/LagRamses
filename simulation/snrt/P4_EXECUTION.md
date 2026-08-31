# P4 execution contract

## Thermal model

The local `Cooling_Grackle` generator evolves a Grackle ionization-equilibrium chemistry network and stores net rate plus mean molecular weight on `(n_H, Z, T)`. It must run **offline**, once per requested scale factor. `build_grackle_atlas.py` combines those subtables into a runtime atlas on `(a, n_H, Z, T)`. The simulation/RT runtime only interpolates this atlas; it never invokes Grackle.

The atlas is used only to invert hydro pressure into temperature and provide background cooling/heating. It does not export H/He species fractions, so SNRT starts those fractions neutral and subsequently evolves them with its own photon-conserving H/He chemistry. Local S_N source photo-heating and photoionization are not encoded in the atlas and must never be double-counted.

AMR covering-grid padding cells with non-positive density or pressure are not treated as resolved gas. Provided the region also contains valid hydro cells, P4 replaces only those padding cells with the atlas minimum `n_H` and minimum `T`; this is a declared near-vacuum floor needed by the static-array RT contract. If the output pressure is zero throughout a selected sample, P4 can explicitly use the atlas net-rate-zero temperature fallback, but such a target is an ingestion validation artifact rather than a physical thermal-state result.

## Target selection

1. Read documented sink ID, mass, and position from `sink_00016.info`.
2. Use the most massive non-boundary sink only to center a 64-cubed local hydro probe.
3. Select the 32-cubed subvolume with maximal mean hydrogen density.
4. Restage that exact subvolume by uniform `arbitrary_grid` resampling from AMR/hydro files into canonical HDF5.

The final density criterion, not sink mass, determines the target. The sink merely avoids an impractical whole-box AMR scan.

## Source ledger

`sink_00016.csv` has two columns not defined by `sink_00016.info`. P4 does not identify either as an accretion rate and emits no AGN photons from it. The output roster lists every sink inside the selected volume for a later audited mapping. No star-particle source data are present in this P4 input.

## Run

```bash
./.venv/bin/python tools/build_grackle_atlas.py \
  --generator /home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/grackle_cooling_grid \
  --grackle-data /home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5 \
  --scale-factors config/p4_thermal_atlas_scale_factors.txt \
  --work-directory data/grackle_subtables \
  --output data/p4_thermal_atlas.h5

./.venv/bin/python tools/p4_stage_high_density.py \
  --info /gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/info_00016.txt \
  --sink-info /gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/sink_00016.info \
  --thermal-atlas data/p4_thermal_atlas.h5 \
  --scale-factor 0.208365 \
  --output data/p4_high_density_rt_input.h5 \
  --scratch .p4_scratch \
  --manifest data/p4_high_density_manifest.json \
  --sink-roster data/p4_high_density_sinks.csv
```

The canonical HDF5 and manifest are an ingestion artifact, not a dual-AGN science result. An AGN RT run requires an audited instantaneous-accretion-rate column and a declared AGN SED/group conversion.

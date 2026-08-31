# P4-A: RAMSES hydro snapshot ingestion

## Scope

P4-A makes a RAMSES hydro state usable by the static S_N radiation solver. It does not start the dual-AGN science analysis and it does not depend on a RAMSES `rt` output. Radiation sources are a separate, explicit catalogue.

## Canonical staging format

`snrt_core.snapshot` defines a versioned HDF5 input file. The file has proper-cgs fields:

- `/gas/hydrogen_number_density_cm3` and `/gas/helium_number_density_cm3`: nuclei number density.
- `/gas/temperature_k`: gas temperature.
- `/gas/dust_relative_abundance`: dimensionless multiplier for the supplied dust cross section.
- `/ionization/x_hii`, `/ionization/x_heii`, `/ionization/x_heiii`: initial primordial number fractions.
- `/sources/cell_index` and `/sources/photon_luminosity_s`: optional source positions and photon-number luminosities per RT group.
- `cell_width_cm` and `grid/left_edge_cm`: proper grid geometry.

The source catalogue is deliberately not inferred from gas cells, particle masses, BH accretion rates, or a RAMSES `rt` output. The science workflow must supply an audited stellar/BH SED-to-group luminosity conversion.

## Current snapshot audit

The metadata at `/gpfs/kjhan/Hydro/Sidm/Agn/Run0/run_cdm/output_00016/info_00016.txt` loads as a `yt` `RAMSESDataset` with a 1024-cubed root domain and redshift 3.7993. Its field index cannot be built when particle files are present: the particle header gives a two-value `nstar_tot` record where the stock `yt` frontend expects one value. A project-local hydro-only view, containing copied text metadata and links to only `amr`/`hydro` rank files, builds successfully. The available hydro fields are density, pressure, three velocities, and magnetic left/right fields; this output has no cooling fields.

This is a particle/frontend compatibility boundary, not evidence of bad hydro data. `stage_ramses_hydro_only` creates the safe view below a caller-specified project scratch subdirectory and uses the standard AMR/hydro reader. No source snapshot file is changed. An audited particle frontend patch remains a separate future task if particle catalogues must be read through `yt`.

## Adapter contract

`stage_ramses_with_yt` and `stage_ramses_hydro_only` require an explicit `RamsesFieldMap`. The map supplies density and either a temperature field or thermal pressure plus a declared mean molecular weight. It never guesses field names, units, or thermal composition. For the current output, use `('gas', 'density')`, `('gas', 'pressure')`, and an explicitly justified `mu`; do not use yt's derived temperature without recording its composition assumption. The adapter converts density to neutral primordial H/He with the declared hydrogen mass fraction (default 0.76), then writes the canonical HDF5 file only after the full grid has been assembled.

An initial staging run should use a modest uniform region and record these items alongside the HDF5 file:

- snapshot path and output number
- bounding box, selected AMR level, and proper cell width
- field map and each source-field unit
- hydrogen fraction, dust normalization, and initial-ionization recipe
- source-catalogue construction and SED group definitions

## Local dependencies and check

Install optional staging dependencies only in `snrt/.venv`:

```bash
./.venv/bin/pip install -r requirements-ramses.txt
JAX_PLATFORMS=cpu .venv/bin/python tests/p4_ingestion.py
```

The check validates cgs H/He conversion, geometry, source-group catalogue, and HDF5 round-trip. It does not touch the production RAMSES snapshot.

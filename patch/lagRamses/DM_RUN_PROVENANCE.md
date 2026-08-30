# Dark-matter run provenance

Every normal output now receives a master-rank
`dm_run_provenance_<output>.txt` sidecar next to its copied `namelist.txt` and
`compilation.txt`.  It is written from `output_amr.kjhan.f90`; it never calls
the sink merger, changes a force, or changes the CDM/SIDM/FDM dynamics.

The file has one `key = value` record per line.  It identifies the active
realization as follows:

- `cdm`: particle-mesh (`pic=.true.`) with neither SIDM nor FDM enabled;
- `sidm`: the Monte-Carlo SIDM branch, with its cross-section and scattering
  controls;
- `fdm`: the Schrodinger--Poisson branch, with axion and wave/HJM controls;
- `none`: no particle DM branch was active.  This is deliberately not relabelled
  as CDM.

FDM and SIDM are already mutually exclusive in `read_params`; encountering
both flags at output is fatal.  Failure to create, write, flush, or close the
sidecar is also fatal before the output can receive its `COMPLETE` marker.

## Relation to the SMBH capture ledger

`sink_particle.kjhan.f90` remains a model-agnostic recorder of the numerical
sink-capture state.  Its event UID should be joined to the surrounding output
run by the output time, copied namelist, compilation record, and this DM
sidecar.  Do not copy SIDM scattering physics or FDM wave quantities into the
sink merger routine.

For a CDM/SIDM/FDM comparison, use separate runs with common initial galaxy
and SMBH conditions.  Their capture events need not occur at identical steps;
compare their provenance-bound ensembles rather than forcing a one-to-one
event match.  FDM still requires the separate
`fdm_outer_wave_provenance_<output>.txt` and full wave snapshots before an
outer response can be evaluated.

The sidecar establishes only run identity and active numerical-model settings.
It is not a physical capture classification, a local density profile, or a
coalescence-time estimate.

Validate a sidecar before joining it to a capture ensemble:

```bash
python3 patch/lagRamses/aux/validate_dm_run_provenance.py \
  output_00042/dm_run_provenance_00042.txt
```

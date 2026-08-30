# Raw resolved-physics inventory

`output_amr.kjhan.f90` writes
`resolved_physics_inventory_<output>.txt` after every raw normal-output
component has been flushed and before the root `COMPLETE` marker is created.
The inventory changes no force, sink state, or CDM/SIDM/FDM evolution.  It is
only a fail-closed index for the later model-specific postprocessor.

The companion `dm_run_provenance_<output>.txt` sidecar records the configured
SMBH merge radius.  A non-compacting zoom is represented only by
`smbh_merge_radius_cells = 0` together with
`smbh_compaction_mode = no_finite_radius_rmerge_zero`.  This configuration provenance
does not by itself establish resolved binary evolution or a physical delay.

The record identifies the selected output, active DM model, raw snapshot
directory, particle/hydro/potential/sink prefixes, and the availability of the
stars, gas, and active dark-matter channels.  A particle dump does not by
itself establish a stellar force channel, so its status is
`requires_particle_classification` until postprocessing applies the explicit
particle-type convention.  Gas is `available` only when the hydro branch is
active.  `none` is never relabelled as CDM.

The potential prefix is accompanied by `potential_checkpoint_status`:
`validated`, `unvalidated`, or `absent`.  It is the normal-output Poisson
restart-validity state, not a source-decomposed SMBH force measurement.

The V1 inventory intentionally records these as `unavailable`:

- source-decomposed SMBH force ledger;
- time-resolved conservation ledger;
- cumulative SIDM scattering ledger.

Those data are not produced by a normal output today.  A downstream tool must
not replace them with zero or call the model-specific result complete.  The
existing SIDM parameter/probability sidecar is configuration provenance, not a
cumulative scattering history.

For FDM, the inventory points to the full `fdm_<output>.out*` snapshot and,
when enabled, the compact outer-wave provenance.  It records only
`resolved_wave_only` accounting; analytic FDM drag cannot be registered beside
a resolved wave wake.

Run the schema check only after the root output has its `COMPLETE` marker:

```bash
python3 patch/lagRamses/aux/validate_resolved_physics_inventory.py \
  output_00042/resolved_physics_inventory_00042.txt
```

Schema validity is not model-specific physics acceptance, a calibrated delay,
or SMBH coalescence evidence.  It only tells the postprocessor exactly what it
must derive from the raw output and what is still absent.

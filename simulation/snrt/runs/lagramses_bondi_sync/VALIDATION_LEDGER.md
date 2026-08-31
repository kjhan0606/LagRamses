# lagRamses Bondi sink-coordinate validation

## Active source and change

- Active VPATH precedence selects `patch/lagRamses/sink_particle.kjhan.f90` before the `cuRamses` fallback.
- `kjhan_sync_sink_particle_coordinates` now copies the cloud-averaged `xsink` and `vsink` into the canonical `PTYPE_SINK` particle.
- The call is placed after `kjhan_create_cloud(1)` and before `make_tree_fine`, so the particle tree is rebuilt after the coordinate update.
- `bondi_hoyle` now selects the unique canonical `PTYPE_SINK` by its validated sink index rather than requiring exact `xp == xsink` equality.
- The binary used for this run is `/home/kjhan/BACKUP/lagRamses/bin/ramses_final3d`, built with `make HDF5=1 ramses` from that active source tree.

## Runtime gate

- Slurm job: `459034` on Grammar `debug`.
- Input: one-sink, Bondi, `nrestart=1` regression seeded from the read-only `target2_baseline_o1/output_00001` checkpoint.
- AGN was enabled only to exercise the ledger and sinkprops writers; both feedback coupling fractions were zero.
- The simulation completed in 19 s and wrote `output_00002/COMPLETE`.
- Every emitted `NaN_CHK sink_post_bondi` count was zero.
- Native outputs: `sink_00000.dat`, `sink_00001.dat`, and `agn_coarse_state_v1.jsonl`.

## Ledger checks

- `p7_convert_sinkprops.py` converted both native sinkprops files successfully after the final source change.
- The two converted records had finite values for every numeric field, including sink mass, Bondi rate, Eddington rate, inflow rate, saved feedback energy, position, and velocity.
- Direct AGN coarse-state JSON contained exactly two finite records: coarse steps `0` and `1`, both for sink ID `1`.

## Scope

This is an active-binary integration and no-NaN regression.  The one-sink symmetric setup does not impose a deliberately nonzero `xp - xsink` mismatch; a dedicated adversarial coordinate-mismatch regression remains the next stronger test.

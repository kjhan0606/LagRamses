# lagRamses moving-sink Bondi validation

## Active implementation

- Source used: `/home/kjhan/BACKUP/lagRamses/patch/lagRamses/sink_particle.kjhan.f90`.
- The active VPATH places this source before the `patch/cuRamses` fallback.
- Canonical `PTYPE_SINK` coordinates are synchronized from cloud-averaged `xsink` and `vsink` before rebuilding the particle tree.
- `bondi_hoyle` selects canonical sinks by their validated sink index, not by exact floating-point equality of `xp` and `xsink`.

## Isolated runtime gate

- Slurm job: `459035` on Grammar `debug`, completed in 17 s.
- Binary: `/home/kjhan/BACKUP/lagRamses/bin/ramses_final3d`, rebuilt after both source changes with `make HDF5=1 ramses`.
- Input: 512-sink Bondi test with `drag=.true.`, based on the read-only compatible `target_drag2_512_baseline_o1/output_00001` checkpoint.
- AGN was enabled only to exercise native sinkprops and AGN coarse-state writers.  Both feedback coupling fractions remained zero.
- Runtime gates passed: `output_00002/COMPLETE`, zero `sink_post_bondi` NaN counts, `sink_00000.dat`, `sink_00001.dat`, and `agn_coarse_state_v1.jsonl`.

## Results

- `p7_convert_sinkprops.py` converted both native files: 512 sinks each, 1,024 source rows in total.
- Direct AGN coarse-state output contains 1,024 rows.
- Every numeric value in both 1,024-row ledgers is finite.
- Sink IDs are identical across the two native records.
- All 512 sinks moved between the two records.
- Position-displacement statistics: maximum `5.59543800e-11`, mean `5.59543575e-11` code units.

This explicitly exercises a moving-sink case, rather than only the stationary one-sink integration gate.

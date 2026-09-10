# Evaluated cold PAH/Fe integration: raw cleanup

Project: `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Operator authorizes cleanup of test raw outputs after evaluation. Driver
evaluation and restart comparison are complete; see
[implementation record](medium_physics_implementation_2026-09-10.md).

Exact directories under `.medium-cold.LVkaiC/`:

| Directory | HDF5 file SHA256 |
| --- | --- |
| `live.S2wuGA/output_00001` | `913f950223476b5277e07663da4f017e4e13cd405981b1e17cd67680b1d2abe5` |
| `live.S2wuGA/output_00002` | `8dd537955f119f34fe91c92d785c5a811e7209b24a959a6ef7f014309e75e891` |
| `restart.QmAkTN/output_00001` | `913f950223476b5277e07663da4f017e4e13cd405981b1e17cd67680b1d2abe5` |
| `restart.QmAkTN/output_00002` | `2e99d441604c4f7c002275610ee21e1e223b1ab65380077de4a720b2efe56efe` |

Each directory's corresponding `data_0000N.h5` (60,045,368 bytes) and
`cooling_0000N.out` (2,213,720 bytes) are removed: 8 files,
249,036,352 bytes (237.50 MiB). No active calculation or checkpoint consumer
remains. No other raw directory is included.

Inputs, environment, logs, binaries, build logs, headers and compact results
are retained. Raw data have no retained backup and require rerunning the
preserved input to reproduce. Retained completion markers are historical;
these output directories are no longer valid restart inputs after cleanup.

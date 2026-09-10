# Evaluated Fe test raw-output cleanup

Driver evaluation: see `medium_physics_implementation_2026-09-10.md`.
Status: completed. Post-removal scan found no matching HDF5/cooling raw
products under the exact root; retained directory size is 1.1 MiB.
Operator authorized removal after the test and stage evaluation. These are
only this bundle's completed/failed test products, not other production runs.

Root: `/gpfs/kjhan/LRD_JWST/.medium-fe-live.S0OD4h`.
For each directory below remove only its named `data_NNNNN.h5` and
`cooling_NNNNN.out`. Each HDF5 file is 60,045,368 bytes and each cooling
table 2,213,720 bytes. Total: 18 files, 560,331,792 bytes (about 534 MiB).

| Directory | HDF5 SHA256 |
|---|---|
| `short/output_00001` | `d20a5f671475643124f5f9126baac4cf7a098b049d72e3e8a814fe737437e834` |
| `verified/output_00001` | `255fb669b2c71b36d33856cc05fb79d08779b47c848e0ab91a5a7e3309295282` |
| `verified/output_00002` | `261267210c378b2566604fc4363d743142b1a39fb6d93dffec1355160fb336e7` |
| `verified/restart/output_00001` | `255fb669b2c71b36d33856cc05fb79d08779b47c848e0ab91a5a7e3309295282` |
| `verified/restart/output_00002` | `7d6fdc0b9c9dd6c41d5dfe31bd49faa79279e709ea3a69eb7bc5ec89f97e7ded` |
| `final/output_00001` | `5aae63d92c4f55586ba90047ba9b85d21a46ee5e855e4ddb8629ef3bb6ab1af0` |
| `final/output_00002` | `e8d0b2630b422a58c9b5c4a57b1f886e08b7d4227a4468cfb60e41dcf3b77ae1` |
| `final/restart/output_00001` | `5aae63d92c4f55586ba90047ba9b85d21a46ee5e855e4ddb8629ef3bb6ab1af0` |
| `final/restart/output_00002` | `288238bfd115b63f764812c3825dc66c27aa2ac0ef70d9c11a37c07f3bb33d30` |

Namelists, environment files, logs, headers/build records, final executable,
source captures and compact results are retained. No active comparison or
restart needs these dumps. Removed raw files have no retained backup here;
reproduction requires rerunning the retained inputs. Output directory marker
files are historical metadata only, not an available restart after cleanup.

# PAH H2 evaluated raw-output cleanup

Performed after native/live/restart tests and driver end evaluation passed,
under the operator's standing retention instruction. Only three resolved
regular HDF5 files were permanently unlinked; no recursive sweep, production
archive or other worker's outputs were touched.

Base: `/gpfs/kjhan/LRD_JWST/.pah-h2.syBrDR/`.

| Removed relative path | Bytes | SHA256 before removal |
| --- | ---: | --- |
| live/output_00001/data_00001.h5 | 12709280 | 7a4f3a0f34e041572f06b30df936ded3546430f9ebc5be636fcd5bcfb642a9f1 |
| live/output_00002/data_00002.h5 | 12709280 | 78ebbfa9a9670b7982d5f85b4c38d23b925de8c06eba6e851121f5780dd2f21a |
| restart/output_00002/data_00002.h5 | 12709280 | 07e70efe5d27b07253c892a0b8d2e2f032127cb7403c7ba8a7cfde128b55c9d0 |

Total **38,127,840 bytes**. Retain effective namelists/env/physical source
tables, output metadata, all logs, checked binary/build identity, native
smoke evidence and `evaluation.txt`/`evaluate.py`. Different whole-file hashes
reflect run/header metadata; the checked physical datasets match exactly.

No trash recovery. These output directories are no longer usable checkpoints
even though metadata/COMPLETE markers and directory symlinks remain. Rerun
the retained small fixture if raw arrays are needed again. An unused control
namelist was removed before execution per the lean plan review; no control
simulation was launched.

# D03 spectral evaluated raw-output cleanup

Native and live/restart driver evaluation is complete. Under standing user
instructions, remove exactly the following three raw HDF5 files, each
7,700,608 bytes (total23,101,824bytes), after confirming these hashes.

| Absolute target | SHA256 |
| --- | --- |
| `/gpfs/kjhan/LRD_JWST/.d03-spectral.LTodPp/live/output_00001/data_00001.h5` | `8b96f2f81a9bbf27ca3184998d68408f40b2ff0c0c6d3002fa299ec70e253a51` |
| `/gpfs/kjhan/LRD_JWST/.d03-spectral.LTodPp/live/output_00002/data_00002.h5` | `50c4cca2a01cb6c66320cd41f0bbffd49ace140c6e953208f006270772174b42` |
| `/gpfs/kjhan/LRD_JWST/.d03-spectral.LTodPp/restart/output_00002/data_00002.h5` | `9fce359c885e48c80fa3ac7499f7a5bfb46e8c65cb860ae6168e5d7f0ecfd0db` |

Retain binary, effective inputs/environment, original source tables,
manifests, logs, compact evaluation and evaluator. Restart output1 is a
symlink to the first target, not a fourth raw dump. Removal is permanent;
future restart requires regeneration. No other output/source/archive is
authorized by this manifest.

Execution DONE: pre-removal sizes/hashes matched; all three raw targets are
confirmed absent. Only23,101,824bytes of evaluated simulation raw data removed.

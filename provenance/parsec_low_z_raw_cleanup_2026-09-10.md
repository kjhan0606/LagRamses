# PARSEC low-Z evaluated raw-output cleanup

State: COMPLETE. All three exact targets were permanently deleted and
their absence verified;165313728bytes removed.
Operator requests deletion after tests AND driver evaluation. MPI2/OMP2
fresh4 steps and restart2->4 completed;82 physical datasets are exactly
equal, all numeric datasets finite, SNRT attributes and physical clocks
equal, mass conserved and density/thermal/CR positive. Retained evidence:
`.parsec-lowz.zbQdfU/evaluation.txt`, effective NML/environment, source
packages/manifests, binary identities, logs and
[driver disposition](parsec_low_z_precision_implementation_2026-09-10.md).
No further active comparison requires these checkpoints.

Exact regular, nonsymlink HDF5 targets beneath
`/gpfs/kjhan/LRD_JWST/.parsec-lowz.zbQdfU/`:

| Relative path | Bytes | SHA256 before deletion |
| --- | ---: | --- |
| live-gradual/output_00001/data_00001.h5 | 55104576 | d28e4a361cb43cd58b44caa64962aed9da9ea5274f6c9674054437f3baba5645 |
| live-gradual/output_00002/data_00002.h5 | 55104576 | c691c4b530bb137af0b313107e490cf72919df9204b8282757c8bca37ff5d273 |
| restart-gradual/output_00002/data_00002.h5 | 55104576 | 42694c04075b09fc62d52a25b9509154b608a5997d732527d5fd369c8befc614 |

Total165313728bytes. `restart-gradual/output_00001` is a link to the fresh
checkpoint, not a fourth raw file. Output-directory metadata is retained.
The failed FTZ `live/` run produced no raw dump. Physical source ZIPs,
generated input tables, logs and other workers' runs are NOT deletion targets.
Deletion is permanent; regenerating data requires rerunning retained inputs.

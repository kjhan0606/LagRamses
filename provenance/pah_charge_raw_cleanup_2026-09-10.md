# PAH charge comparison raw-output cleanup

Root: `/gpfs/kjhan/LRD_JWST/.pah-charge-live.9ReUG2`.
Authority: operator completed-test retention directive, project AGENTS.md.
Evaluation: bounded charged comparison and exact continuous/restart PASS;
see `medium_physics_implementation_2026-09-10.md`. The earlier `evolution`
experiment used a superseded carrier-mass conversion and is excluded from
acceptance evidence. Its diagnosis is complete; no checkpoint is needed.

Exact raw targets below, each `data_NNNNN.h5` (60,657,048 bytes) and matching
`cooling_NNNNN.out` (2,213,720 bytes), total 377,224,608 bytes (~359.75 MiB):

| Directory relative to root | HDF5 filename | HDF5 SHA256 |
| --- | --- | --- |
| evolution/output_00001 | data_00001.h5 | ceb910395bf3342626f4a78879b5807a44d320995e93f46796ba0ee9afb690e6 |
| evolution/output_00002 | data_00002.h5 | fc32d6186149ce6e0b79f10c6dececb44a69bde3dc5f40e56b76d299ccacb4ca |
| fixed.9jxVyS/output_00001 | data_00001.h5 | b624dde9b014dcb6c11e25cb27e84987cbd6c9c463335a5b3a1ab2cf25721b41 |
| fixed.9jxVyS/output_00002 | data_00002.h5 | 8d4aaea682840cca6c48077fa9c6fd047534185e0f894fbbab5c01d0d4abb05b |
| restart.TALRDn/output_00001 | data_00001.h5 | b624dde9b014dcb6c11e25cb27e84987cbd6c9c463335a5b3a1ab2cf25721b41 |
| restart.TALRDn/output_00002 | data_00002.h5 | ae30078cf51d5bb2c9309e9032ba15c8bb21d8b31d9e339c6d17b9c0578beb6b |

Retain all namelists, environment files, logs, executables/build logs, output
metadata, original optical/chemistry inputs and compact `physical-check.txt`
and `restart-check.txt`. Whole-file restart hashes differ because metadata
differs; the compared physical datasets are exactly identical. No unrelated
run or scientific/production archive is targeted. Raw files are not backed
up here; regeneration requires rerunning the retained inputs and binary.

Status: all twelve listed raw files removed; scoped post-deletion search
found no remaining data HDF5/cooling dumps in these three runs. Metadata
directories (including historical COMPLETE markers) are retained, but are
no longer usable restart checkpoints. Logs/inputs/builds/results are intact.

# Evaluated raw output cleanup

Standing operator authorization in AGENTS.md: remove simulation raw output
after the test and driver evaluation complete. Driver evaluation is recorded
in `snrt_secondary_spectrum_implementation_2026-09-10.md`.

Removed permanently (not trash), three files,4,161,000bytes. No running
simulation or pending restart needed these files. Inputs, logs, binary,
evaluator and compact result remain in `.rt-node-secondary.WxSPuz/`.

Paths below are relative to `/gpfs/kjhan/LRD_JWST/.rt-node-secondary.WxSPuz`;
each was1,387,000bytes. Exact filenames were resolved and hashed before
deletion; no recursive directory or wildcard removal was used.

| File | SHA256 before removal |
| --- | --- |
| live/output_00001/data_00001.h5 | acad41c6b5ac66efd780cf751f8156e2c9bfb4f502e3af050424ded12eeaf85b |
| live/output_00002/data_00002.h5 | 681666060c88596de2dea6565fa632f425ff656c60f4db35dd548b12574ad705 |
| restart/output_00002/data_00002.h5 | cfe49772959b1e068384eac8729c8d0bcc412e36ebd49d26ff8aec8872f62b0e |

Retained output-directory metadata/COMPLETE markers and the restart symlink
do not constitute usable checkpoints after this removal. Reproduction must
run afresh from the retained input, not attempt to reuse those directories.

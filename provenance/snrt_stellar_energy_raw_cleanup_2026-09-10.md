# Evaluated stellar Q/E raw-output cleanup

Performed after the native/live/restart tests and driver end evaluation
passed, under the operator's standing raw-output retention instruction.
Only these exact regular HDF5 files were unlinked, permanently. No recursive
cleanup, production archive or other worker's directory was touched.

Base: `/gpfs/kjhan/LRD_JWST/.rt-stellar-energy.kzyggZ/`.

| Removed relative path | Bytes | SHA256 before removal |
| --- | ---: | --- |
| completed/output_00001/data_00001.h5 | 1387000 | 56dc2ceebfdb77f93e6869f44e0c3d95bc13ce91c02fe3ab72e69ab644622eda |
| completed/output_00002/data_00002.h5 | 1387000 | 194d9797784b7fdd47a5174a08af1bf557fe8ca5aa159a2cbb7a59dc209f0514 |
| completed-restart/output_00002/data_00002.h5 | 1387000 | a1f955314d68421b3861ccc64699cb6487878f754aa7c4121fd260845ce41533 |

Total **4,161,000 bytes**. Retained: effective namelists/env/source tables,
all logs/build identity and executables, output metadata, `evaluation.txt`,
`evaluate.py`, native smoke evidence and the implementation report. Whole
HDF5 hashes differ across continuous/restart due to run/header metadata;
the 82 physical datasets and physical clocks match exactly, as recorded.

These are no longer usable checkpoints, regardless of retained COMPLETE
markers or output-directory symlinks. There is no trash recovery; repeat the
small simulation with retained inputs if raw arrays are needed again.

# PARSEC five-Z evaluated raw cleanup

Status: COMPLETE; all three exact files removed and absence verified.
Operator requested raw deletion after tests AND driver evaluation. The
[completed evaluation](parsec_metallicity_extension_implementation_2026-09-10.md)
finds all82 physical datasets and actual SNRT identity/physical clocks exact
between uninterrupted and restarted outputs, with finite/positive fields.
Both runs have finished; their step2 checkpoint is no longer needed.
Only the following exact raw HDF5 files under
`/gpfs/kjhan/LRD_JWST/.parsec-metal.TNrmgs/` are in scope.

| Relative file | Bytes | SHA256 |
| --- | ---: | --- |
| live/output_00001/data_00001.h5 | 32825104 | a4328f34d37416b2161dbb6982faeaff0d3ab9cb2bbd50484e5c601c53ba6756 |
| live/output_00002/data_00002.h5 | 32825104 | d9a8628eeba2ccd116e23cd14a8a884e914fae4a7d86bdb6b314d48917400365 |
| restart/output_00002/data_00002.h5 | 32825104 | 4d185b45b923da5015d4dfd67fdf5ed32580944dddd43d8f1cd2917f6d9e15bf |

Total98,475,312bytes. Restart/output_00001 is a link to the original live
checkpoint, not a fourth raw copy. Retain all physical source archives,
derived input packages, effective NML/environment, logs, source/native tests,
binary/build identity, evaluation.txt and small output metadata. No other
worker/run or production archive is swept. Removal is permanent; retained
inputs permit rerunning, not byte-for-byte recovery of deleted HDF5 files.

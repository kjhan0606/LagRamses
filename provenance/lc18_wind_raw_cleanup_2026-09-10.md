# Evaluated LC18 wind raw cleanup

Status: COMPLETE; all eight exact files removed and absence verified.
Operator requested deletion of evaluated simulation raw outputs. Only the
eight files below, rooted at `/gpfs/kjhan/LRD_JWST/.lc18-native.X4wP2A/`,
are in scope. Prior evaluation is in
[real-source progress](real_source_integration_progress_2026-09-07.md#lc18-actual-wind---native-feedback-explicit-approximation-2026-09-07):
fresh/restart complete, material identity and source budgets agree, velocity
roundoff separately quantified; changed-speed restart rejected. These old
checkpoints are not needed for the current LC18 prompt input or PARSEC run.
Input tables, histories, logs, binaries and existing compact evidence remain.
Total4,335,488bytes. Removal is permanent; rerunning retained inputs is the
reproduction path, not a promise of recovering the same HDF5 bytes.

| Relative file | Bytes | SHA256 |
| --- | ---: | --- |
| fresh/output_00001/data_00001.h5 | 448752 | 03cfe2aea8a6e04eae3c1ae6416611d85a7cb1629dbb1fedb6d19df5d6ad864a |
| fresh/output_00002/data_00002.h5 | 506096 | 70ac273a380a96271ad38ac8bd8fa9aff2b3df3a0df008845f1d330f86e0f016 |
| fresh/output_00003/data_00003.h5 | 563440 | b0a881ac968538211b50e7f9787f4d00c2d603e911c66d41ead719b0ef191ff3 |
| fresh/output_00004/data_00004.h5 | 620784 | 4546958e8b43a4f32244e549df9d5259bfe20a9bf7535a325b9dd2b0af250cc8 |
| restart/output_00002/data_00002.h5 | 506096 | 70ac273a380a96271ad38ac8bd8fa9aff2b3df3a0df008845f1d330f86e0f016 |
| restart/output_00003/data_00003.h5 | 563440 | 8efa8bd8748926b2e12f71f22e853b2c3fd921ca305482316b8116658e2dbd04 |
| restart/output_00004/data_00004.h5 | 620784 | e92a628cbc1660001242ec4ed0069a8d73a8d782782e83a6fbeeb71c3f0c7fd2 |
| reject-speed/output_00002/data_00002.h5 | 506096 | 70ac273a380a96271ad38ac8bd8fa9aff2b3df3a0df008845f1d330f86e0f016 |

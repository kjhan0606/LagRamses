# PARSEC pair-feedback evaluated raw-output cleanup

Driver end evaluation is complete; native source tests, continuous/restart
79-dataset equality and conservation are recorded in
[implementation](parsec_pair_feedback_implementation_2026-09-10.md).
Under the operator's standing retention instruction, remove only these
three evaluated raw files, each 6,456,376 bytes (total 19,369,128 bytes).

| Absolute target | SHA256 before removal |
| --- | --- |
| `/gpfs/kjhan/LRD_JWST/.parsec-pair.dnTIS0/live/output_00001/data_00001.h5` | `1bea026d72399e5dc9674cd69bf9a34b9ebfea68bc93fdc5a0ece345896dec2f` |
| `/gpfs/kjhan/LRD_JWST/.parsec-pair.dnTIS0/live/output_00002/data_00002.h5` | `37182908a1f14547cc595ab8dbb87b13e36a4ba33afbcff53bc124da4de3ce5f` |
| `/gpfs/kjhan/LRD_JWST/.parsec-pair.dnTIS0/restart/output_00002/data_00002.h5` | `053abc85695b48228182f1d9380ab3854a5f72c0403d463f5afe4d15d98915f9` |

Retain effective inputs, source packages (including rejected initial
candidate), manifests, logs, binaries and compact evaluation. Restart's
output1 symlink refers to the first target and is not a fourth dump.
Removal is permanent; future restart requires regeneration. No source
archive, active checkpoint, other run or production output is a target.

Execution: all three hashes/sizes matched immediately before exact-path
removal; all three raw paths are confirmed absent afterward. Cleanup DONE.

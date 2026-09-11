# Evaluated test snapshot cleanup

Scope: `/gpfs/kjhan/LRD_JWST/.chimes-sources.r1TaD8` only. Operator's standing
instruction permits removing raw test outputs after evaluation. Inputs,
logs, binaries, build logs, measured JSON results and this manifest remain.
Snapshots below are either a completed stellar regression or superseded
failed/contaminated integrations. No active restart depends on these copies.
Deletion is not recoverable locally; rerun retained inputs to regenerate.

## First cleanup set

Delete each listed `output_NNNNN` directory, including its raw auxiliary
outputs. SHA256 below identifies the HDF5 data file within that directory.

| Directory relative to run root | HDF5 SHA256 |
|---|---|
| stellar/output_00001 | 3f2f44bf7c4ff10882c2e86a79f847b43470e54365001639820ff25b2f7f4f07 |
| stellar-diagnostic/output_00001 | 3f2f44bf7c4ff10882c2e86a79f847b43470e54365001639820ff25b2f7f4f07 |
| stellar-positive/output_00001 | 3f2f44bf7c4ff10882c2e86a79f847b43470e54365001639820ff25b2f7f4f07 |
| stellar-positive/output_00002 | eb9b779776587884ba3713d893d27f454013d8e2599c084f07db76d0d568ffc7 |
| coupled-final/output_00001 | 2766b66ee756b2f16e2676b893e621578779ea45b9da307be6964d8cb38ffeed |
| coupled-live/output_00001 | e800814743298dc225049cf8af662dab52ba53be5e4f0f8d14ad4ffc485e12cc |
| coupled-live/output_00002 | 59e8a53e4d0d9b3dd5725a0ed85fa88011066bd8002e3efb531821695bf9f0d7 |
| restart/output_00001 | e800814743298dc225049cf8af662dab52ba53be5e4f0f8d14ad4ffc485e12cc |

First set: removal completed. No unrelated run was touched.

## Final cleanup set

Final driver evaluation completed before removal. Both successful final
runs and diagnosed/superseded attempts are included; the latter are not
relabeled as passes. Retained compact results include `final-results.json`,
`ready-results.json`, `resolved-halfbox-results.json`, `integration-results.json`
and `stellar-positive-results.json` in the private run root.

| Directory relative to run root | Bytes removed (all files) | HDF5 SHA256 |
|---|---:|---|
| coupled-diag/output_00001 | 64165138 | 7ba2e4e2f7cf5f481c7288e5c4c29d8edd459f93ba0772b6434fe5ac433d185f |
| coupled-bounded/output_00001 | 64167714 | 7dd03c430d0c4c358aa9040b93d4bb3412c8d61f03f58b49ce80376909f2641e |
| coupled-verified/output_00001 | 64165138 | 31efc3b19b244c8ba6aeed7cf8bfb118b77795837477f938448e0a7a5c596a34 |
| restart-resolved/output_00002 | 64178690 | 78ce4a80fc5281bb608c2c2cbfa86b3314fe88d410fdaa51c3d2462c7ccbf5a5 |
| restart-verified/output_00001 | 64165138 | 31efc3b19b244c8ba6aeed7cf8bfb118b77795837477f938448e0a7a5c596a34 |
| restart-resolved/output_00001 | 64165138 | 6f2d3458a2426c8b0de401446fe59b5b1f2520dc7ff5ab1a022b819153fdccdf |
| restart-diag/output_00001 | 64165138 | 31efc3b19b244c8ba6aeed7cf8bfb118b77795837477f938448e0a7a5c596a34 |
| coupled-resolved/output_00002 | 64178690 | 25821dc2a162c65549c9bf7f8a02abefb705ffffe40cf8ed8bc636613ba9bb64 |
| restart-budget/output_00001 | 64165138 | 31efc3b19b244c8ba6aeed7cf8bfb118b77795837477f938448e0a7a5c596a34 |
| coupled-resolved/output_00001 | 64165138 | 6f2d3458a2426c8b0de401446fe59b5b1f2520dc7ff5ab1a022b819153fdccdf |
| coupled-ready2/output_00001 | 64167714 | 5c1270586bc2f5a5af986489672d9919be866dae1e49bef2320197a86db955be |
| coupled-ready2/output_00002 | 64180930 | 634302321fafcb7c921058dcd17551ec4df51d502768233ed4dbdabb10e8a368 |
| restart-ready2/output_00001 | 64167714 | 5c1270586bc2f5a5af986489672d9919be866dae1e49bef2320197a86db955be |
| restart-ready2/output_00002 | 64180930 | 2abdaf159512f7eddd7d1b9b7cd11f0af59a6b343f8b6ae133c7d08a20640f9b |

Final set total: 898378348 bytes in14 directories. Input copies needed by the
completed restarts are now included; no active calculation depends on them.

Final set: removal completed. All22 listed snapshot directories across both
sets are removed. Source data banks, other workers' directories and all
retained inputs/logs/results are untouched.

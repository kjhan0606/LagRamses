# Evaluated PARSEC phase-wind raw output cleanup

Operator: delete test raw outputs after tests AND stage evaluation complete.
Driver evaluation PASS is recorded in
[implementation](parsec_phase_wind_implementation_2026-09-10.md).
MPI2 continuous/restart comparison completed; no checkpoint is still needed.
Only the three exact files below are in scope. No production archive, source
asset, other worker run, input, binary or compact evidence is removed.

Root: `/gpfs/kjhan/LRD_JWST/.parsec-wind.Soq8Rf/`.
Each file13098976bytes; total39296928bytes (39.30MB).

| Relative file | SHA256 before deletion |
| --- | --- |
| live/output_00001/data_00001.h5 | cbd98b73b0248173df2d6c6ecadafd68f960e3bff630f5173ce03f2fdf57862c |
| live/output_00002/data_00002.h5 | ea83bcfa9131648c4c898d3b4267b716f2946fad288dcbb8586857a728bc7b33 |
| restart/output_00002/data_00002.h5 | e4b339da52825cd23bf0dca4039ae41b25f67b496a2dc4c9fb4af3ea35baec66 |

Status: COMPLETE. Exact-file deletion succeeded; absence of all three files
was verified. Recovered39296928bytes of logical raw-file storage.
Retained: both effective namelists/environment, physical source packages and
archives, source/build/binary identities, runtime/native logs,
`evaluate.py/evaluation.txt`, `check_sources.py/check-sources.log`, and driver
review. The restart/output_00001 symlink is not a separate snapshot.
Deletion is permanent; reconstruct outputs by rerunning retained inputs.

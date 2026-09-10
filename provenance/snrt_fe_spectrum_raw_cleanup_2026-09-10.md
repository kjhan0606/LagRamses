# Static Fe spectral raw output cleanup

Status: COMPLETE. All three exact files permanently deleted; absence verified.
Deletion authorized by the operator's retention policy. Driver evaluation
completed before selecting these targets. No built-in recovery copy remains.
Root `/gpfs/kjhan/LRD_JWST/.fe-spectral.MDeWoz/`.

| Raw path relative to root | Bytes | SHA256 before deletion |
|---|---:|---|
| live-fixed/output_00001/data_00001.h5 | 7623736 | 7471b35ef331ca7a7c01afa7627680ceaeaa8174605730ac0cc15134bfaae6cd |
| live-fixed/output_00002/data_00002.h5 | 7638584 | 60115d371f1ba66ca2de156dc7e67a04ad19d1d10ac40477d170a6949b638e4b |
| restart/output_00002/data_00002.h5 | 7638584 | 4a5cec49e62edde3dd2f5495501d58fe1b0781ee79bf1920ba1b21dcb4b18c7f |

Total22900904 bytes, three files. Restart output1 is a symlink to the first
directory, not a fourth raw file. Failed `live/` produced no raw snapshot.
Keep all effective namelists/source controls, environment, physical source
archives, generated optics, logs, binaries, compact `evaluation.txt`,
`evaluate.py`, source-resolution results and this manifest. No input/physical
archive or another worker's output is part of the cleanup.

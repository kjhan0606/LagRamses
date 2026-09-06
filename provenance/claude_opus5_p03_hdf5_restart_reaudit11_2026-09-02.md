# P0.3 HDF5 Stellar Restart Gate — Independent Audit 11 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict at audit time: BLOCK

This audit was run after the AMR/sink checked-read implementation and before
the follow-up correction recorded below.  It confirmed that the AMR and
ksection changes were real, but found a new correctness defect in the
`nindsink` consistency condition and identified missing source-to-binary
provenance in the durable contract.

## Findings

| Finding | Audit-11 result |
|---|---|
| Sink payload extent/status | Closed in source; the audit noted that the then-current fixture did not exercise a nonzero sink payload. |
| `nindsink` consistency | **HIGH defect:** equality with `maxval(idsink)` rejects legitimate post-merger checkpoints because `nindsink` is monotonic while a merged high-ID sink can disappear. The binary restart convention accepts `nindsink >= maxval(idsink)`. |
| AMR coarse/distributed payload | Closed: coarse `son/cpu_map/flag1`, distributed `xg/son_flag`, same-ncpu payload, and ksection CPU boxes use checked readers. |
| ksection safety | Closed: range checks, runtime allocation equality, checked tree arrays, and two-rank positive/oversized-dimension negative coverage. |
| ADM payload | Checked in source but not runtime exercised by the fixture. |
| Runtime negative evidence | Genuine for the AMR negative and existing count/ksection negatives; sink payload and `nindsink` consistency were not yet exercised. |
| S1 uninterrupted star continuation | Remains open; nonzero `PTYPE_STAR` state was still h5py-injected, with the evidence flag correctly false. |
| Source/build provenance | **HIGH gap:** the durable contract recorded binary hash/mtime but not the HDF5 source hashes, git HEAD, or dirty-worktree state. |
| Hydro/gravity payload | Medium-high residual: `uold`, `phi`, and related payload readers remain outside this P0.3 checked-read surface. |
| Per-rank checked readers | Medium residual: early validation returns in the 1-D checked readers lack MPI consensus and can produce a root-only diagnostic. |

## Required correction

1. Change the sink guard to reject only `nindsink_file < maxval(idsink)` and
   retain the file value when it is larger.
2. Add a real one-row sink payload to the temporary checkpoint so the normal
   linked restart executes every sink dataset reader, plus an
   `nindsink=0` negative.
3. Record source hashes for `ramses_hdf5_io.f90`, `restore_hdf5.f90`,
   `backup_hdf5.f90`, `output_part.f90`, and the build Makefile together with
   git HEAD and dirty status.

The audit did not close P0.3. A follow-up audit is required after these
corrections; the S1 continuation and hydro/gravity scope remain separate
residuals even if the corrected restart surface passes.

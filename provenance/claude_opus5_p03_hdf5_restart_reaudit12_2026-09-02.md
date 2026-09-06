# P0.3 HDF5 Stellar Restart Gate — Independent Audit 12 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

The audit-11 corrections are real, present in the linked binary, and bound by
durable source/build provenance. The corrected `nindsink` rule now preserves a
larger monotonic file value and rejects only a value below the active sink-ID
maximum. A complete temporary `nsink=1` row is read by the normal linked
restart, and the `nindsink=0` negative is diagnostic-tied and self-guarded.

## Closed findings

- AMR coarse, distributed, and same-ncpu payloads use checked readers with
  MPI-consensus extent/status handling; the truncated `son_flag` 4-rank
  negative is genuine runtime evidence.
- ksection dimensions, allocation equality, tree arrays, and CPU boxes are
  checked; the real two-rank positive and oversized-dimension negative remain
  passing.
- All sink payload datasets use checked rank/extent/status readers, and the
  normal linked restart reaches the post-sink `HDF5 particle restore done.`
  marker.
- Binary evidence now records SHA256 for the relevant HDF5/particle sources,
  the Makefile, git HEAD, dirty state, and binary identity. Independent hashes
  match the durable JSON artifact.
- The audit-11 `nindsink` regression is closed.

## Remaining condition for unconditional P0.3 PASS

The three stellar ledger attributes `nstar_tot`, `mstar_tot`, and `mstar_lost`
are still read through unchecked helpers in `restore_part_hdf5`. They must be
converted to checked scalar attribute reads with MPI consensus and fail-closed
diagnostics. ADM has checked source reads but no runtime fixture; S1
uninterrupted `PTYPE_STAR` continuation remains honestly open; hydro/gravity
payload checking and per-rank diagnostic consensus remain carried-forward
scope. These do not erase the conditional closure of the AMR/sink/restart
surface, but the stellar ledger gap must be resolved before final P0.3 PASS.

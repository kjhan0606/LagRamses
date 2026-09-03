# P0.3 HDF5 Stellar Restart Gate — Independent Audit 14 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

Audit 13's H1/H2 condition is closed. The complete header clock/cosmology
surface now uses checked scalar/string/array helpers with MPI consensus, and
`dtold/dtnew` are read at the checkpoint extent before bounded assignment. The
linked binary and durable evidence include the current source hashes, and the
missing-`aexp` negative is runtime-backed.

## Closed findings

- `nstar_tot`, `mstar_tot`, and `mstar_lost` are checked scalar attributes with
  consensus and a missing-ledger negative.
- Sink reader and writer payloads are exercised with `nsink=1`; the monotonic
  `nindsink` rule and below-maximum negative are correct.
- AMR coarse/distributed payloads and ksection tree/CPU-box data are checked;
  the 4-rank truncated-`son_flag` negative is genuine.
- Source-to-binary provenance is bound by five source hashes, Makefile hash,
  binary hash/size/mtime, git HEAD, and dirty-worktree state.

## Remaining condition

The `stellar_state_schema_version` attribute itself is still read by the
unchecked `hdf5_read_attr_int` after an existence-only probe. A malformed
rank/extent or failed `H5Aread` can therefore overflow or compare an
undefined value. Convert this schema-marker read to the checked scalar reader,
use MPI consensus/fail-closed abort, and add a malformed-attribute negative.

Explicitly carried forward: S1 uninterrupted `PTYPE_STAR` continuation,
ADM runtime coverage, hydro/gravity payload checked reads, and non-root
diagnostic consensus. The plan narrative also needs to point to the latest
audit rather than old reaudits.

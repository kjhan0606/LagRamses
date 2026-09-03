# P0.3 HDF5 Stellar Restart Gate — Independent Audit 13 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

Audit 12's stellar-ledger condition is closed. The linked binary and durable
JSON now cover checked `nstar_tot`, `mstar_tot`, and `mstar_lost` attributes,
including a diagnostic-tied missing-`nstar_tot` 2-rank negative. Sink reader
and writer coverage, `nindsink` monotonicity, AMR/ksection payload checks, and
source-to-binary provenance remain closed.

## Closed in this audit

- A complete temporary `nsink=1` row is read by the normal linked restart;
  sink writer output is checked for all required fields and exact lengths.
- `nindsink_file < maxval(idsink)` is the only rejection condition, preserving
  valid post-merger checkpoints where the file counter is larger.
- The distributed AMR checked collective readers use MPI-consensus extent and
  status validation, with a real 4-rank truncated-`son_flag` negative.
- The binary evidence records the relevant source/Makefile hashes, git HEAD,
  dirty state, and binary hash/size/mtime.

## Condition carried to the next audit

The `/header` time, cosmology, and schedule attributes are still read through
unchecked helpers, including `time`, `aexp`, `nstep_coarse`, `dtold/dtnew`,
and `ordering`. The audit classified this as HIGH because a malformed
checkpoint can silently restart at the wrong epoch; it also identified the
unchecked `dtold/dtnew` extent as a possible overflow vector. Hydro/gravity
payload checked reads, ADM runtime coverage, S1 uninterrupted star
continuation, and non-root diagnostic consensus remain explicitly carried
forward rather than claimed closed.

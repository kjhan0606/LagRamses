# P0.3 HDF5 Stellar Restart Gate — Independent Audit 15 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

Audit 14's schema-marker condition is closed: the marker now uses the checked
scalar attribute reader with MPI consensus, and a non-scalar schema attribute
is rejected by the linked 2-rank negative. Header clock/cosmology reads,
stellar ledger attributes, sink/AMR/ksection/particle payloads, sink writer
coverage, and source-to-binary provenance remain closed and runtime-backed as
previously recorded.

## Remaining gate conditions

- The restore path still uses unchecked chunked/1-D reads for hydro `uold_*`
  and gravity `phi`, `f_*`, `scalar_gr`, and optional `psi_re/im` payloads.
  A malformed checkpoint can therefore silently inject incomplete state.
- `hdf5_read_dataset_1d_*_checked` does not MPI-reduce rank-local failures
  before its collective HDF5 transfer, so a non-root failure may be
  fail-closed by `MPI_ABORT` without a diagnostic; this is a diagnosability
  residual, not a newly observed hang.
- ADM has checked source reads but no runtime coverage, and S1 uninterrupted
  `PTYPE_STAR` continuation remains intentionally false.
- Several restore branches still call `clean_stop`; rank-local branches can
  deadlock through collective `MPI_FINALIZE` if reached asymmetrically.

## Next priority

Convert the hydro/gravity chunk and 1-D payload reads to checked readers with
exact expected cell extents and add a malformed `uold` or `phi` linked
negative. Then replace rank-local restart `clean_stop` calls with the
fail-closed abort path and address optional ADM/error diagnostics.

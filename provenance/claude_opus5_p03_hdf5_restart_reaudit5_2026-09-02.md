# P0.3 HDF5 Stellar Restart Gate — Independent Audit 5 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

The `npart_per_cpu` implementation is substantially correct, but P0.3 remains
open. S1, the real-star continuation gate, still blocks a final production or
publication PASS.

## S2/S2r assessment

`hdf5_read_dataset_all_int_checked` is MPI-collective-safe: dataset open and
every HDF5 stage are reduced with `MPI_Allreduce(MAX)` before the next stage or
the broadcast, and cleanup is guarded. `restore_hdf5.f90` now rejects negative
entries and a count sum different from `npart_total_file` before using the
partition. Missing-dataset and inconsistent-sum negatives run at `nproc=2`
and assert field-specific diagnostics. The same-`ncpu=2` fixture checks
per-rank values.

## S2/S2r residuals

1. The negative-entry branch has no runtime fixture; a negative array can have
   the correct sum, so the sum test does not cover it.
2. There is no wrong-length `npart_per_cpu` fixture that exercises checked
   extent status 14; the current inconsistent case has length one and reaches
   the sum check.
3. The same-`ncpu=2` view is Python-injected metadata from a one-rank output,
   not a real two-rank writer, and this injection is not yet labeled in JSON.
4. `ncpu_file` and `npart_total_file` are obtained through unchecked HDF5
   attribute helpers, leaving S8 on the critical path.

## Observed versus injected

The linked writer and post-restart HDF5/native-binary payloads are observed;
the nonzero release values and PTYPE_STAR marking are injected into a
temporary checkpoint. The ncpu=2 metadata view and all corruptions are also
fixture-injected. The evidence correctly flags the true PTYPE_STAR
continuation equivalence as false.

## Remaining findings

- **S1 HIGH:** no real-PTYPE_STAR uninterrupted-vs-restart continuation
  equivalence; `indtab` no-double-counting remains unmeasured.
- **S2r HIGH:** add negative-entry, wrong-length, and real two-rank writer
  coverage; label injected metadata.
- **S8 MEDIUM/HIGH:** checked HDF5 attribute reads are required for header
  values controlling allocation and restart branch selection.
- **S9 MEDIUM:** `ngrid_per_cpu` remains unchecked and lacks partition
  validation.
- **S3 MEDIUM:** ADM particle reads remain relaxed.
- **S5/S6 MEDIUM:** source hashes/HEAD/dirty state are absent and a failed
  runtime may leave a stale passing artifact.
- **S4 LOW/MEDIUM:** binary parser still assumes eight header records.
- **S7 LOW:** synthetic and literal evidence fields need clearer labels; the
  scope string omits the nproc=2 cases.

## Next closure order

1. Checked header attributes (`ncpu`, `npart_total`) and explicit count
   negative/wrong-length fixtures.
2. Checked AMR partition and ADM reads; source/evidence failure provenance.
3. Real star formation plus uninterrupted/restart equivalence for S1.

Until S1 and the residual S2/S8/S9 items are closed, P0.3 remains
`CONDITIONAL PASS`.

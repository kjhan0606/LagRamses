# P0.3 HDF5 Stellar Restart Gate — Independent Audit 6 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

The checked header attribute path is now closed for the values that control
particle restart allocation and branch selection. P0.3 remains open because
S1 real-star continuation and residual S2r/S3/S9/provenance work remain.

## Closed in this iteration

- `hdf5_read_attr_int_checked` and `_int8_checked` validate open, dataspace,
  scalar shape, and read status; scalar and length-one attributes are accepted.
- `ncpu`, `nvar`, and `npart_total` call sites reduce status across MPI ranks
  before shared abort, and reject invalid `ncpu_file`/particle totals.
- `npart_per_cpu` checked read, nonnegative-entry and sum consistency checks,
  nproc=2 missing/sum negative cases, and same-ncpu per-rank values are
  present and pass.

## Remaining S2r gaps

1. No runtime fixture directly exercises a negative `npart_per_cpu` entry.
2. No wrong-length `npart_per_cpu` fixture reaches checked extent status 14.
3. The same-ncpu=2 metadata view is rewritten in Python from a one-rank
   checkpoint, not emitted by a real two-rank writer, and is not explicitly
   labeled injected in the evidence artifact.

## Remaining findings

- **S1 HIGH:** no real `PTYPE_STAR` uninterrupted-vs-restart continuation
  equivalence; no-double-counting via `indtab` remains unmeasured.
- **S3/S9 MEDIUM:** ADM datasets, `ngrid_per_cpu`, and other AMR partition
  reads remain relaxed and lack partition validation.
- **S8 residual MEDIUM:** other allocation-controlling attributes remain on
  relaxed readers.
- **S5/S6 MEDIUM:** source hashes/HEAD/dirty state are absent and a failed
  runtime can leave a stale passing artifact.
- **S4 LOW/MEDIUM:** binary comparison still hard-codes eight header records.
- **S7 LOW:** runtime scope omits nproc=2 fixtures and evidence labels need
  further cleanup.

## Required next actions

Add the negative-entry and wrong-length runtime fixtures, label the injected
same-ncpu metadata view, then close the remaining checked AMR/ADM and
provenance items before attempting S1 with a real star-forming continuation.

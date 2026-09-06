# P0.3 HDF5 Stellar Restart Gate — Independent Audit 4 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

The `npart_per_cpu` checked-read implementation is substantially correct, but
P0.3 remains open. S1, the real-star continuation test, still blocks
publication readiness.

## S2 status

`hdf5_read_dataset_all_int_checked` returns distinct status codes for its
HDF5 stages, reduces failures with `MPI_Allreduce(MAX)` before any rank can
enter the broadcast, and cleans up handles with guards. `restore_hdf5.f90`
wires it to `npart_per_cpu` with a field-specific diagnostic and abort. The
missing-dataset negative and the same-`ncpu=2` per-rank fixture pass.

## Remaining S2 gaps

- All negative fixtures still run at one process, so collective fail-closed
  behavior is not directly exercised under MPI failure.
- The count array is only shape-checked. It does not reject negative entries
  or verify `sum(npart_per_cpu) == npart_total_file`; an in-range but wrong
  partition can silently misroute particles.
- The same-`ncpu=2` view is Python-injected by rewriting header and partition
  metadata on a one-rank output rather than emitted by a real two-rank writer.

## Observed versus injected

Observed from linked runtime: the initial zero payload, post-restart HDF5 and
native-binary streams mapped through the emitted descriptor, nproc=4 per-rank
diagnostics, and negative exit/diagnostic behavior. Injected by the fixture:
all nonzero release values, `ptypep=1`, the ncpu=2 metadata view, and all
corruptions.

## Residual findings

1. **S1 HIGH:** no real-`PTYPE_STAR` uninterrupted-versus-restart equivalence;
   `indtab` no-double-counting remains unmeasured.
2. **S2r HIGH:** multi-rank fail-closed negatives and count-array
   non-negativity/sum validation are required.
3. **S8 MEDIUM:** HDF5 header attribute reads, including `ncpu_file`, swallow
   status even though they control allocation and restart branch selection.
4. **S3 MEDIUM:** ADM particle reads remain on the relaxed reader and are
   untested.
5. **S9 MEDIUM:** `ngrid_per_cpu` partition reads remain unchecked.
6. **S5/S6 MEDIUM:** source hashes/HEAD/dirty state are absent and a failed
   runtime can leave a stale passing artifact.
7. **S4 LOW/MEDIUM:** the binary parser still assumes eight header records
   outside the descriptor.
8. **S7 LOW:** synthetic/literal evidence fields need clearer labeling.

## Closure criteria

For S2, run missing and wrong-length `npart_per_cpu` negatives at `nproc>=2`,
reject negative counts and wrong sums with field-specific diagnostics, and
generate the ncpu=2 checkpoint from a real two-rank writer. Overall PASS also
requires S1, checked header/AMR/ADM reads, and source/evidence provenance
closure.

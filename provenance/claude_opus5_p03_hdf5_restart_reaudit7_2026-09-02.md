# P0.3 HDF5 Stellar Restart Gate — Independent Audit 7 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

P0.3 remains open. S2/S2r is closed for fail-closed behavior; S1 real-star
continuation still blocks a final PASS.

## Closed

- `hdf5_read_dataset_all_int_checked` is MPI-safe: every HDF5 stage is
  `MPI_Allreduce(MAX)`-reduced before the next stage and before the broadcast,
  with guarded cleanup.
- `restore_hdf5.f90` checks negative `npart_per_cpu` entries before the sum,
  and checks the sum against `npart_total_file`.
- At `nproc=2`, the missing, negative-entry, wrong-length, and sum-mismatch
  fixtures reach their own field-specific diagnostics.
- Critical `ncpu`, `nvar`, and `npart_total` attributes use checked readers and
  rank-consensus before abort.
- Evidence labels writer-observed payloads, injected restart payloads, and
  injected same-ncpu metadata honestly.

## Remaining findings

1. **S1 HIGH:** no real `PTYPE_STAR` uninterrupted-vs-restart continuation
   equivalence; `indtab` no-double-counting remains unmeasured.
2. **S2r residual MEDIUM:** no real two-rank writer and no test for a
   sum-consistent but wrongly partitioned count array.
3. **S8 MEDIUM:** `noutput`, `nlevelmax_file`, and other allocation-controlling
   attributes remain relaxed; `noutput` sizes `tout/aout` reads.
4. **S9 MEDIUM:** `ngrid_per_cpu` and other AMR partition arrays remain on
   relaxed readers with no extent or sum validation.
5. **S3 MEDIUM:** ADM particle arrays and sink metadata remain relaxed.
6. **S5/S6 MEDIUM:** source hashes/HEAD/dirty state are absent and failed
   runtime can leave stale passing JSON.
7. **S4 LOW/MEDIUM:** binary parser retains the hard-coded eight-record header
   offset.
8. **S7 LOW:** literal/synthetic fields and the effective MPI coverage need
   clearer evidence labels.

## Next action

Close S9 with checked AMR partition reads and validation, then address ADM,
remaining metadata/provenance, and finally build the real-star continuation
fixture required for S1.

# P0.3 HDF5 Stellar Restart Gate — Independent Audit 10 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

The ksection OOB defect that caused audit 9's FAIL is closed. P0.3 remains
open because sink/AMR payload status handling and the real-star continuation
gate are not complete.

## Verified closures

- **N2 ksection:** file dimensions are read into locals with checked,
  MPI-consensus attributes; range and strict equality against allocated
  runtime arrays are checked before reshape; `ksec_wall/ksec_next` are checked
  reads; real two-rank writer/restart and oversized-dimension negative run.
- **N1 nsink:** checked/bounded `nsink`, checked `nindsink`, and two-rank
  oversized-sink negative.
- **N3 nlevelmax_file:** checked and bounded before level-indexed arrays.
- **S9:** `ngrid_per_cpu` checked with negative/wrong-length nproc=2 fixtures.
- **S2/S2r:** checked/count-validated `npart_per_cpu`, collective negatives,
  and same-ncpu rank checks.

## Residual blockers

1. **S1 HIGH:** no real `PTYPE_STAR` uninterrupted-versus-restart
   continuation; all nonzero state remains h5py-forged and the evidence flag
   is correctly false.
2. **S3a HIGH:** approximately twenty sink datasets after the `nsink` guard
   still use relaxed readers. An in-range `nsink` with a short `msink` can be
   silently accepted. `nindsink` is later overwritten by `maxval(idsink(...))`.
3. **AMR HIGH:** `xg_*`, `son_flag`, and `cpu_map` use collective readers that
   do not validate extent/status; `bisec_cpubox_min/max` is an unchecked
   reshape.
4. ADM `dark_energy_int`/`dark_h2_frac` remain relaxed; some allocation and
   schedule metadata remains unchecked.
5. Binary header offset, source hashes/dirty state, stale artifact behavior,
   and evidence labels remain incomplete.

## Next order

Route sink datasets and AMR payload through checked readers with extent/range
validation and add in-range truncation negatives. Then finish ADM/metadata and
source provenance before attempting the real-star continuation equivalence.

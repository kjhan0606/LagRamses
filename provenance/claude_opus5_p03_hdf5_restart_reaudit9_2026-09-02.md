# P0.3 HDF5 Stellar Restart Gate — Independent Audit 9 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: FAIL

The linked contract test passes, but a newly claimed ksection closure is
unsafe and the ksection runtime branch has no coverage.

## Findings

### N1 — closed with residual sink payload risk

`nsink` now uses checked attribute read, MPI consensus, and an `nsinkmax`
bound; the two-rank oversized-sink negative reaches its diagnostic. However,
the subsequent `idsink` and approximately twenty sink datasets still use
relaxed readers, and `nindsink` is overwritten by `maxval(idsink(...))`.

### N2 — ksection dimensions can cause heap OOB (HIGH, introduced)

The reachable same-ordering ksection branch reads file `ksec_kmax` and
`ksec_nbinodes` into the module globals after `init_amr` has already allocated
`ksec_wall`, `ksec_next`, and `ksec_indx` using the runtime dimensions. It then
writes file-sized reshapes into those existing arrays without reallocation or
an equality check. The new bounds against `ncpu` are insufficient; for
non-power-of-two `ncpu`, `2*ncpu-1` is not a general k-ary tree bound. The
flattened wall/next datasets are also relaxed reads.

The ordering-enabled ksection branch has zero runtime coverage. The variable
restart ksection build path is internally unreachable because the outer
condition guarantees the inner mismatch test.

### N3 — unchecked `nlevelmax_file`

The file level count is read without status and is used to index arrays sized
by the requested runtime `nlevelmax`.

### Other open items

- S1: no real `PTYPE_STAR` uninterrupted-vs-restart continuation;
- S3: ADM reads relaxed;
- S4: binary parser assumes eight header records;
- S5/S6: source provenance and stale artifact handling incomplete;
- S7: literal pass markers and effective MPI coverage need cleanup;
- `noutput_file` controls discarded `tout/aout` reads through a relaxed helper;
- AMR payload (`son`, `cpu_map`, `flag1`, `ksec_*`) lacks complete checked
  status and range validation.

## Required closure

1. Read ksection dimensions into locals, validate exact equality with existing
   runtime allocations (or safely reallocate), and checked-read wall/next;
   reject unsupported metadata before any reshape.
2. Add an ordering-enabled ksection runtime fixture and a negative dimension
   case; remove/dead-code-proof the impossible variable-ksection branch.
3. Check `nlevelmax_file`, remove or check discarded output schedule metadata,
   then address sink/ADM/AMR payload readers and provenance.
4. Build the real-star continuation equivalence required for S1.

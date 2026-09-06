# P0.3 HDF5 Stellar Restart Gate — Independent Audit 8 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

S9 is closed as scoped: all three `ngrid_per_cpu` reads use the checked
collective helper, reject negative entries, and the nproc=2 negative fixtures
reach their own diagnostics. Header `ncpu/nvar`, particle `npart_total`, and
sum-validated `npart_per_cpu` are likewise protected. P0.3 remains open.

## Newly found high-severity gaps

### N1 — unchecked `nsink` can overrun fixed sink arrays

`nsink` is read unchecked in the sink restore path and then used in slices and
loops without a bound against `nsinkmax`. A corrupt checkpoint can request
more sinks than the allocated arrays hold.

### N2 — unchecked k-section dimensions control allocation/reshape

`ksec_kmax` and `ksec_nbinodes` are read without checked status and drive
allocation and reshape into fixed k-section structures. The current fixture
does not set `ordering`, so this path is not exercised.

### Additional metadata gap

`noutput_file` sizes the temporary `tout/aout` attribute reads through a
relaxed helper, although those values are discarded. Removing the unnecessary
read or checking the metadata is required.

## Scientific blocker

S1 remains HIGH: the linked writer creates no real `PTYPE_STAR` and no
uninterrupted-versus-restart comparison measures released mass, metals,
energy, or `indtab` no-double-counting. The evidence labels the nonzero
checkpoint payload as injected and the continuation flag as false.

## Required next actions

1. Add checked/bounded `nsink` restore and a corrupt-`nsink` negative fixture.
2. Add checked k-section metadata and reject impossible dimensions before
   allocation/reshape, with an ordering-enabled negative fixture.
3. Check or remove `noutput_file`/`tout/aout` reads and extend provenance.
4. Continue to S1 real-star continuation after all restart safety guards are
   closed.

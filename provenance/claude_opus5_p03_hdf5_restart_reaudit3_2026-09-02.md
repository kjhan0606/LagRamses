# P0.3 HDF5 Stellar Restart Gate — Independent Audit 3 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Verdict: CONDITIONAL PASS

P0.3 remains open for production/publication purposes. Five of the six prior
FAIL findings are genuinely closed. The remaining scientific evidence gap is
explicitly declared, and two residual code/provenance gaps must be addressed
before a full PASS.

## Prior finding status

1. Default runs overwriting the evidence JSON: **closed**. Only
   `--fortran-runtime` writes durable evidence.
2. Requiring the initial writer payload to be zero: **closed**. Zero-ness is
   recorded, while post-restart nonzero payload is asserted.
3. Tautological binary order: **substantially closed**. Record mapping is
   derived from the emitted `part_file_descriptor.txt` and checked by bitwise
   equality.
4. Status-swallowing non-release particle reads: **partially closed**. The
   ordinary particle fields, birth epoch, and metallicity now use checked
   readers; `npart_per_cpu` and optional ADM fields remain open below.
5. Nonspecific negative checks: **closed**. Schema, `indtab`, and truncated
   `ptypep` cases require both nonzero exit and field-specific diagnostics.
6. Fixture-injected nonzero/PTYPE_STAR state: **not closed**. The evidence
   labels this honestly and sets
   `true_ptype_star_continuation_equivalence: false`.

## Directly observed versus injected

The production linked writer emits the initial zero payload. The linked HDF5
reader restores the temporary nonzero state, and the linked HDF5 and binary
writers independently re-emit the values bit-for-bit. This proves the
serialization chain. The origin of the nonzero values and `PTYPE_STAR` is
still a Python mutation of the temporary checkpoint; no linked RAMSES run in
this gate has yet created a star and advanced the release cursor.

## Residual findings

### S1 — No true star-formation continuation equivalence (HIGH, open)

The gate has not compared an uninterrupted linked run against a
checkpoint-and-restart branch for released mass, metals, energy, and `indtab`.
Moreover, the restore guard is controlled by `star .or. sink`, not by the
restored per-particle type, so the injected `PTYPE_STAR` is not itself proof
that the star feedback path was exercised.

Closure requires a linked fixture that produces a real `PTYPE_STAR`, advances
the stellar release driver, and asserts uninterrupted versus restart
equivalence at a common final time.

### S2 — `npart_per_cpu` read remains unchecked (HIGH, code)

`hdf5_read_dataset_all_int` discards HDF5 open, dataspace, read, and broadcast
status. Its values determine `npart_loc` and `offset_part`, including the
same-`ncpu` branch. A malformed count array can therefore misroute subsequent
particle reads.

Closure requires a status-returning checked variant, fail-closed restore,
same-`ncpu` multi-rank coverage, and a corrupted/missing `npart_per_cpu`
negative fixture with a field-specific diagnostic.

### S3 — Optional ADM particle reads remain relaxed (MEDIUM)

`dark_energy_int` and `dark_h2_frac` still use the relaxed reader when ADM is
enabled. Sink attributes also need an explicit status policy.

### S4 — Binary header count is hard-coded (MEDIUM)

The descriptor-derived field mapping still uses `record_index = 8 + ivar - 1`.
The header-record count should be emitted and read as part of the descriptor.

### S5 — Source provenance is incomplete (MEDIUM)

The artifact records the binary hash and mtime but not the audited source
hashes, `HEAD`, or dirty-tree state. The binary was built from uncommitted
working-tree edits, so source-to-binary linkage must be explicit.

### S6 — Failed runtime can leave a stale prior artifact (MEDIUM)

If the runtime test raises, the prior passing JSON remains on disk. Write an
explicit failed evidence record or remove the artifact on runtime failure.

### S7 — Synthetic-only keys need clearer labeling (LOW)

Synthetic round-trip fields and literal guards should be labeled as such so
they cannot be mistaken for production evidence.

## Required next order

1. Close S2: checked `npart_per_cpu`, same-`ncpu` multi-rank test, and a
   field-specific negative case.
2. Close S3/S4/S5/S6/S7 as evidence and robustness cleanup.
3. Close S1 with real star formation and uninterrupted/restart equivalence.

Until S1 and S2 are closed, P0.3 remains `CONDITIONAL PASS`, not a final PASS.

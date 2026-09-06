# Claude Opus 5 P0.3 HDF5 restart-state gate audit

Date: 2026-09-02  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: Claude Opus 5 (read-only independent gate audit)  
Scope: P0.3 HDF5 stellar restart state, not a production simulation

## Verdict

**CONDITIONAL PASS — P0.3 cannot be closed yet.**

## Confirmed implementation strengths

- The three continuation fields `tpp`, `mp0`, and `indtab` are the complete
  stellar release state used by the current runtime. The transient progress
  object is rebuilt on each call and `indtab` is exported only after
  deposition.
- The HDF5 packing path and the binary restart path use the same active-particle
  predicate and ordering. The existing packing-count guard protects the
  `npart`/active-particle correspondence.
- MPI offsets, rank-count changes, zero-particle ranks, unrestored slots, and
  binary-restart compatibility reuse the existing particle dataset path and
  were judged structurally consistent.
- The committed-progress cursor preserves the no-double-counting invariant
  across a restart when the serialized state is valid.

## Findings

1. **HIGH — runtime evidence is empty.** The linked Fortran fixture produced
   `release_dataset_lengths: {"tpp": 0, "mp0": 0, "indtab": 0}`. Therefore no
   packing loop executed and no non-empty value fidelity was tested. The
   `np.all(empty)`-style zero-fixture assertion is tautological.
2. **HIGH — Fortran fail-closed branches were not executed.** Missing-schema or
   missing-dataset rejection was checked by a Python reimplementation and
   source guards, but not by an actual Fortran restart with a deliberately
   incomplete HDF5 file.
3. **MEDIUM — HDF5 I/O status is swallowed.** The shared
   `ramses_hdf5_io.f90` helpers do not propagate every HDF5 status. A
   present-but-truncated `tpp`, `mp0`, or `indtab` dataset can therefore leave
   an output buffer invalid while the restart continues. Presence checks alone
   do not provide fail-closed I/O semantics.
4. **LOW — evidence metadata needs tightening.** The fixture strips
   `stellar=.true.` and the synthetic `bitwise_synthetic_roundtrip` field must
   represent an actual non-empty byte comparison, not only a fixed byte-size
   check. The scope wording must distinguish “no production output opened”
   from “no runtime output opened”.

## Required closure actions

- Run a linked Fortran writer/reader roundtrip with at least one active stellar
  particle and value-level assertions for all release fields.
- Add a multi-rank and `ncpu_file != ncpu` case, including zero-particle ranks.
- Execute a Fortran negative restart test by removing the schema marker or a
  required release dataset and require a nonzero exit.
- Add binary-to-HDF5 bitwise cross-check coverage for the same active-particle
  fixture.
- Propagate and validate HDF5 read status and dataset extent for all three
  release fields before accepting the restored state.

The current bounded synthetic and linked tests are useful prerequisites, but
they do not establish P0.3 closure or production/publication readiness.


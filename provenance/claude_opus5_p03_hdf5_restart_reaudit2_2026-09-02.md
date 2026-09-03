# P0.3 HDF5 Stellar Restart Gate — Independent Audit 2 (Claude Opus 5)

Date: 2026-09-02 · Root: `/gpfs/kjhan/LRD_JWST` · Read-only audit.

## Scope

The audit inspected the production HDF5 writer/reader, the checked HDF5
helper, the P0.3 contract test, and its evidence artifact. The reported
command was:

```text
python3 simulation/snrt/tests/p0_hdf5_stellar_restart.py --fortran-runtime
```

## Verdict: FAIL

P0.3 remains open and is not production/publication ready. The Fortran read
path is materially stronger than at the previous audit, but the evidence
claim is not reliable and the central nonzero-writer assertion is inverted.

## Directly proven versus fixture-injected

If the linked-runtime command is run successfully, the implementation and
fixture directly exercise:

- schema attribute and `tpp`, `mp0`, `indtab` writes in active-particle order;
- checked rank, extent, hyperslab, read, and close status for release-critical
  fields and fail-closed aborts;
- reader-to-memory-to-writer fidelity for values placed in the checkpoint;
- negative missing-schema and bad-extent runs;
- a 1-to-4 rank restore with a zero-particle rank.

However, the nonzero release values and `PTYPE_STAR` are injected by Python
into a temporary checkpoint. The linked writer initially emits zeros, the
fixture allocates through the sink path, and no actual star forms. The
binary/HDF5 comparison also compares all-zero payloads.

## Findings

### H1 — Evidence artifact contradicts the claimed run

`simulation/snrt/data/p0_hdf5_stellar_restart_contract.json` currently records
`production_fortran_runtime_roundtrip: false` and the synthetic-only scope.
The test rewrites that file on every run, so a default test run can overwrite
the evidence from a linked-runtime run. The artifact is untracked, and no
build hash/mtime is recorded.

### H2 — The writer is asserted to be zero

The test raises unless the initial production writer release payload is all
zero. A legitimate nonzero writer payload would fail the gate, so the prior
HIGH finding about writer nonzero fidelity has been inverted rather than
closed.

### H3 — Nonzero-state helper is dead code

`_write_phase0_runtime_table` exists but is not called; `_run_binary` removes
`PHASE0_YIELD_TABLE` when no explicit table is passed. The test therefore does
not drive a production writer to create a nonzero state.

### H4 — No uninterrupted-versus-restart continuation equivalence

The test does not compare an uninterrupted branch with a checkpoint/restart
branch for released mass, metals, energy, and `indtab`. The no-double-counting
cursor invariant remains source reasoning rather than a measured result.

### M1 — Non-release particle fields remain status-swallowing

`x_*`, `v_*`, `mass`, `identity`, `levelp`, and `ptypep` still use relaxed
HDF5 helpers that discard `ierr`. A malformed `ptypep` or `levelp` can therefore
silently misalign the particle type and release state that were checked.

### M2 — Negative tests only check nonzero exit

The missing-schema and bad-extent tests do not assert their specific failure
diagnostics; unrelated launch or namelist failures could satisfy them.

### M3 — Evidence fields are tautological or hard-coded

Record order is compared to a Python constant, while several booleans are
literal `True` or initialized without being derived from observed data. The
artifact must distinguish `writer_observed` from `fixture_injected` values.

### M4 — Multi-rank coverage is shallow

Only rank-local counts are checked. Per-rank values, an `ncpu_file > 1` input,
and non-contiguous active slots are not covered. The compact writer should
explicitly document that holes are not representable in this checkpoint
schema.

### L1 — Legacy-checkpoint migration path is undocumented

Fail-closed rejection of a legacy DM-only checkpoint when stellar restart is
enabled is scientifically correct, but a converter or explicit migration
procedure must be documented.

### L2 — Scope wording is ambiguous

The runtime fixture launches the production binary, while the scope says no
production output was opened. It should say no production *data* was used and
identify the runtime output as a fixture.

## Required closure criteria

1. Run the linked test and persist matching evidence with binary hash/mtime.
2. Drive a nonzero writer payload and assert distinct nonzero `tpp/mp0/indtab`.
3. Repeat binary/HDF5 comparison on that nonzero payload and derive binary
   record order from file metadata rather than a Python constant.
4. Add uninterrupted-versus-checkpoint/restart continuation equivalence with
   at least one real `PTYPE_STAR`.
5. Route all restart-critical particle fields through checked readers.
6. Assert diagnostic-specific negative failures and add truncated `ptypep`.
7. Add per-rank values, `ncpu_file=4 -> ncpu=1`, and compact-slot wording.
8. Remove tautological booleans and label observed versus injected evidence.
9. Document legacy-checkpoint migration.

Items 1–5 are required before P0.3 can close or P0.4 can be promoted.

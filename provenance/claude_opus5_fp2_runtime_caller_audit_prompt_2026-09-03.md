# Claude Opus 5 audit request — F-P2 SNIa runtime caller and automatic evidence

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, commit, push, launch
RAMSES, or change runtime activation. This is an independent physics and code
audit after the automatic evidence bundle was strengthened on 2026-09-03.

The review target is the F-P2 SNIa runtime caller and its population/physical/
thermal contracts, not unrelated historical gates. Inspect the source, native
mirror, Makefile wiring, source-parity and production-linked-build evidence,
the F-P2 audit JSON, and the invoked tests. Treat the current checkout as the
only authority and cite exact file paths and line numbers.

Required questions:

1. Is the three-group runtime namelist handoff correctly bound to one immutable
   source commit and one approval id, and does `phase0_initialize` remain
   fail-closed while `runtime_activation_allowed` is false?
2. Is interval DTD evaluation restart-deterministic, and is the WD-reservoir
   debit/returned-mass/energy/tracked-ejecta/momentum ledger transactional and
   free of duplicate ownership with the generic SSP driver?
3. Does the persisted-mass reconstruction correctly account for a normal
   retry/restart after the generic population ledger has advanced? Distinguish
   this normal-restart evidence from hard-crash exactly-once semantics; do not
   treat the latter as proven unless an atomic journal exists.
4. Is AMR leaf-cell location, MPI ownership, row-major `unew(n_local_cells,nvar)`
   indexing, OpenMP accumulation, and scratch-to-destination transaction
   correct? Distinguish the actual one-cell NGP runtime path from the separate
   weighted multi-cell bridge test.
5. Are physical closure claims publication-safe: returned mass, terminal
   remnant, tracked and residual ejecta, total energy including bulk kinetic
   energy, signed/vector momentum convention, and element deposition?
6. Do the tests actually run the claimed positive and negative cases, and does
   production source parity plus the linked `ramses_final3d` build prove the
   new module is wired into the production binary? Identify any test that is
   only static, native-mirror-only, or a build smoke test.

Return a structured report:

- top-level verdict: PASS, CONDITIONAL PASS, or BLOCK;
- a table of findings with severity (blocker/major/minor), status
  (implemented/verified/partial/missing), and file:line evidence;
- explicit separation of engineering wiring, physical validity, and
  publication-readiness;
- acceptance gates required before enabling SNIa in production;
- any claims in the local evidence/provenance that are too strong.

Do not recommend activating SNIa merely because the scaffolding compiles. The
safe expected state may be a coherent, production-linked, runtime-gated
implementation with hard-crash journal and multi-cell/MPI qualification still
open.

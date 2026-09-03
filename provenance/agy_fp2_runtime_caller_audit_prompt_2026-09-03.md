# AGY audit request — F-P2 SNIa runtime caller and automatic evidence

Audit `/gpfs/kjhan/LRD_JWST` read-only with the requested AGY/Gemini model.
Do not edit files, commit, push, launch RAMSES, or change runtime activation.
This is an independent physics/code audit after the automatic evidence bundle
was strengthened on 2026-09-03. Use generous reasoning/time and inspect the
actual checkout, not only this prompt.

Review the F-P2 SNIa runtime caller, its population/physical/thermal contracts,
native mirror, Makefile/source-parity wiring, linked production-build evidence,
F-P2 audit JSON, provenance, and the tests that the runner actually invokes.
Cite exact file paths and line numbers.

Assess all of the following:

1. Three-group runtime namelist handoff, immutable source-commit and approval
   binding, ordered loading, and fail-closed `phase0_initialize` behavior while
   `runtime_activation_allowed` is false.
2. Interval DTD mathematics and restart determinism; WD-reservoir ownership and
   debit; returned mass, energy, tracked/residual ejecta, terminal remnant, and
   signed/vector momentum ledger closure; no duplicate generic-driver return.
3. The persisted-mass reconstruction used after a normal restart/retry. State
   precisely whether hard-crash exactly-once is proven; an atomic journal must
   not be assumed from normal-restart tests.
4. AMR leaf-cell lookup, MPI ownership, row-major `unew(n_local_cells,nvar)`
   indexing, OpenMP accumulation, and scratch transaction. Separate the actual
   one-cell NGP runtime path from weighted multi-cell bridge evidence.
5. Whether the positive/negative tests are really invoked by the shell runner,
   and whether source parity plus a linked `ramses_final3d` build establishes
   production wiring rather than merely native-mirror coverage.
6. Whether the current claims meet production-ready/publication-ready physics
   standards, including the remaining multi-cell/MPI, hard-crash, full-yield,
   and momentum-convention qualifications.

Return:

- top-level verdict: PASS, CONDITIONAL PASS, or BLOCK;
- findings table with severity, status (implemented/verified/partial/missing),
  and file:line evidence;
- separate engineering, physical, and publication-readiness verdicts;
- exact gates needed before enabling SNIa;
- overclaims or ambiguities in the local audit/provenance wording.

Do not recommend production activation simply because the implementation builds.
The correct result may be a production-linked but runtime-gated implementation
with hard-crash journal and weighted multi-cell/MPI qualification still open.

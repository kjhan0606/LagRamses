# AGY F-P2 runtime-caller audit — 2026-09-03

Auditor: AGY / Gemini Antigravity CLI, requested model `gemini-3.8-flash-high`.
The audit was run in plan mode with a 30-minute print timeout against
`/gpfs/kjhan/LRD_JWST`. No repository edit, commit, push, or RAMSES run was
requested. Full external artifact:
`/home/kjhan/.gemini/antigravity-cli/brain/e2bb44b5-cb6f-4939-af8e-485027668d75/fp2_snia_runtime_caller_audit.md`.

## Verdict

- Evidence bundle and implementation: **CONDITIONAL PASS**.
- Engineering readiness: **CONDITIONAL PASS**.
- Physical baseline: **PASS, baseline only**.
- Production SNIa activation: **BLOCK**; keep `enable_snia=.false.` and
  `runtime_activation_allowed=false`.
- Publication readiness: **BLOCK**.

## Findings

AGY independently verified the ordered three-group namelist handoff, common
40-hex source binding and approval id, interval DTD additivity, AGB-only WD
reservoir debit, no duplicate generic-driver SNIa return, normal restart mass
reconstruction, AMR leaf lookup, row-major `unew` bridge transaction, OpenMP
updates, source parity, and linkage of the production `ramses_final3d` build.

The following remain open:

1. **Hard-crash exactly-once — major/partial.** Normal restart/retry arithmetic
   is tested, but there is no atomic journal joining `unew`, particle mass, and
   `indtab` checkpoint state.
2. **Runtime spatial deposition — major/partial.** The bridge supports weighted
   multi-cell inputs, but the connected caller passes one local target with
   weight 1.0 (NGP). No genuine MPI boundary exchange or multi-cell runtime
   qualification exists.
3. **Production negative execution — major/missing.** The runner does not start
   the linked production binary with SNIa enabled and assert the fail-closed
   initialization result.
4. **Test environment naming — minor/partial.** The isolated contract test uses
   `FP2_SNIA_RUNTIME_CONTRACT`, while the production caller reads
   `PHASE0_SNIA_RUNTIME_CONTRACT`.
5. **Publication physics — blocker/missing.** Net yields are zero, metallicity
   dependence is unity, and the adopted baseline is a single HESMA n100 model.

## Overclaims identified by AGY

- The old blocker-reconciliation sentence saying thermal coupling is “not
  called by the runtime” is obsolete after the caller connection.
- `production_ready: true` alongside `runtime_activation_allowed: false` is
  ambiguous and should be renamed or qualified as a production-baseline
  approval, not activation readiness.
- Weighted two-cell bridge evidence must not be described as if the active
  runtime caller already uses a multi-cell kernel.

AGY's final recommendation is to keep activation blocked until the hard-crash,
multi-cell/MPI, production-negative, and scientific/publication gates close.

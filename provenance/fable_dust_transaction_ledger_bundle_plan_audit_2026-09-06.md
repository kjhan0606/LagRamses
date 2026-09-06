# Fable bundle-start audit: DUST-8 CUDA ledger to FP64 RAMSES handoff

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Auditor/model: Fable alias through Claude CLI (`--model fable`)
- Mode: read-only plan audit; no edits, compilation, jobs, web search, or delegation
- Date: 2026-09-06
- Verdict: **CONDITIONAL APPROVE**

## Audit result

Fable judged DUST-8 necessary: DUST-7 has a real A10-tested four-species
CUDA ABI but no Fortran caller, so its ledgers are not load-bearing in the
native path. The proposed prepared-cell adapter and direct H/He handoff are
technically realizable. The scope was judged only mildly instrumented; the
separate ledgers are justified by the two closure identities.

## Conditions adopted

1. The native zero-dust routing must be explicit and auditable, with a named
   `ZERO_SCAFFOLD` mode and a once-per-run log/evidence record. It must not be
   described as active dust transport or a dust-physics result.
2. The zero-dust path must have an old-three-species versus new-four-species
   equivalence gate at FP32 round-off, including H/He ledgers and the actual
   thermochemistry input.
3. The FP32-to-FP64 tolerance must be documented as a formula involving
   FP32 epsilon, group count, and substep count rather than an unexplained
   hard-coded number.
4. The direct H/He CUDA ledger must bypass host reconstruction/repartition of
   the assigned total; source grep and a negative test must support this.
5. Dust ledgers are validated and quarantined as trial diagnostics. No dust
   thermal, momentum, abundance, or other persistent state is committed, and
   a negative test must show the dust state is unchanged.
6. The old ABI must remain source-compatible and be present in the native
   link. P4 sidecar loading, dust-to-metal mapping, AMR/legacy changes remain
   outside this bundle.

## Minimal evidence requested

Fable requested a focused A10 smoke for finite/non-negative/shape and both
closure identities, an old/new equivalence table with the tolerance formula,
a dust-state-unchanged negative test, a native link listing both ABI symbols,
and an evidence record containing the scaffold mode, commit/GPU facts, and
tolerance parameters.

The audit was obtained after removing the accidental `--bare` invocation and
the inherited `CLAUDE_OPTS` from the audit process. The initial long
file-reading invocation produced no response; the successful audit used the
same plan facts in a tool-disabled read-only prompt. This provenance note
does not treat the failed invocation as an audit result.

# Claude Opus 5 F-P1.2 bundle-end re-audit — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5 (primary; read-only)
Prompt: `provenance/opus5_fp1_2_stellar_feedback_transaction_bundle_end_audit_retry_prompt_2026-09-04.md`
Execution: `claude --model opus --permission-mode plan --allowed-tools Read --output-format json --no-session-persistence`
Result: completed successfully after the prior 300-second primary timeout
Model cost reported: approximately USD 2.19

## Verdict

**CONDITIONAL PASS**

Opus judged the F-P1.2 production-source transaction correctly implemented and
the two GPT-5.6-Sol conditions materially met, but identified two bounded
evidence-integrity conditions before final engineering closure.

## Conditions and disposition

### 1. Direct parallel-build dependency — closed

`bin/Makefile:320` now gives `stellar_ramses_bridge.o` a direct
`stellar_ramses_field_map.o` prerequisite, matching the bridge's module import
at `patch/lagRamses/stellar_ramses_bridge.f90:19`.  The field-map target is in
the source graph and the runtime edge remains explicit.  The actual focused
object target was checked with `make -C bin -j4 ...` and the SNRT/CUDA dry-run
contains the field-map, bridge, runtime, and feedback objects in the link
graph.

### 2. Native evidence integrity — addressed

The native runner now exports `OMP_NUM_THREADS=4` and `OMP_DYNAMIC=FALSE`, and
the test rejects a single-thread configuration before claiming synchronization
evidence.  Its same-cell case updates a complete `unew(target,1:nvar)` row
under the striped lock.

The native test now exercises actual builder failures for invalid volume,
invalid field map, and `ndim=2`, asserting zero returned delta.  It also checks
independent expected generic/SNIa density, all momentum components, total
energy, total metal, delayed tracer, and element values, plus opposed
generic/SNIa momenta whose net momentum is zero while per-component kinetic
energy remains positive.

The evidence wording was corrected: builder failure cannot mutate production
row/mass/progress state because neither non-mutating builder accepts those
arrays; the native test checks zero-delta failure behavior and a separate
state model, rather than claiming a live production-array identity test.  The
Python banner now says `failure_model_identity` and `same_cell_model`.

## Confirmed closed implementation items

- Generic and SNIa contributions are separate row-major scratch deltas; the
  mutating SNIa scatter helper is not called by the mixed production prepare
  path.
- Target bounds and uniqueness are checked while legal virtual/reception rows
  are retained for RAMSES reverse exchange.
- The 4096-striped OpenMP lock re-reads the current row, progress is prepared
  before mutation, and the complete row/mass/progress assignments have no
  fallible post-write operation.
- Runtime uses `ndim==3`, energy slot `ndim+2`, and validates density,
  momentum, energy, delayed-cooling, total-metal, and active-element
  non-overlap.
- Generic and SNIa kinetic terms are independent, including bulk, cross, and
  source-momentum terms; zero-mass/nonzero-momentum input fails closed.
- Energy ownership is not double-counted: the generic SSP path excludes the
  SNIa channel, while delayed cooling receives only the SNII returned-mass
  tracer.
- MPI cross-rank atomicity, process-crash journaling, live integration,
  physical yield/fate calibration, SED/dust closure, and dimensional
  generalization remain explicitly deferred.

## Notes retained for later work

- The named OpenMP critical fallback described in the original plan is not
  implemented or tested; the plan has been amended to describe the hashed lock
  as the sole implemented synchronization mechanism.
- The lazy initialization branch inside `phase0_feedback` should be hardened
  in a later bundle, although serial initialization currently precedes the
  OpenMP caller and no current runtime defect was assigned.
- A process-crash journal and distributed/MPI transaction remain outside
  F-P1.2.

No files were modified and no jobs were run by this audit.  The two bounded
conditions above were implemented afterward and passed the focused tests; a
fresh primary audit is still required before marking the bundle fully closed.

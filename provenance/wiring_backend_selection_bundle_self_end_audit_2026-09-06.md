# Codex end audit: wiring and backend-selection bundle

Date: 2026-09-06
Base: `85cfa45`

## Verdict

**PASS within bounded bundle scope; conditional for production/publication.**

The requested SNRT, stellar/AGN feedback, and dust wiring was rechecked using
positive/negative source assertions and the consolidated native/production
gate. The gate passed production link, MPI transaction coverage, A10 CUDA
multigroup transport, dust-ledger receiver checks, and fail-closed negatives.

The only code change in this follow-on bundle removes a dead AGN `gpu_sink`
autotune/switch block. Its namelist spelling is retained but explicitly does
not select a backend. This prevents an unconnected flag from being treated as
feedback acceleration evidence.

The no-auto-selector disposition is correct for the inspected code: SNRT has
no CPU transport equivalent, stellar/AGN feedback has no CUDA equivalent, and
dust has no OpenMP transport equivalent. OpenMP in the active paths is host
parallelism or an FP64 receiver, not a drop-in replacement for CUDA transport.

## Required future conditions

An actual OpenMP/CUDA selector may be reconsidered only after each candidate
backend has an identical state-transition and AMR/MPI ownership contract,
deterministic cross-backend tolerances, conservation/rollback differential
tests, and runtime device-affinity handling. DUST-8's `ZERO_SCAFFOLD` and
candidate IR APIs remain non-live.

Fable did not return a verdict after two corrected/compact five-minute
read-only attempts; that is recorded as unavailable rather than approval.

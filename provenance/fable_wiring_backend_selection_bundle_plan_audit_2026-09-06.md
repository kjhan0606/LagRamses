# Fable plan-audit record: wiring/backend-selection bundle

Date: 2026-09-06
Project: `/gpfs/kjhan/LRD_JWST`
Requested auditor: Fable via Claude CLI

## Status

**UNAVAILABLE — no Fable verdict obtained.** This file is a provenance record
of the failed audit attempts, not an approval and not a substitute for an
external review.

## Attempts

1. The first invocation used the previously documented `--no-subagents` and
   `--tools` options. The installed Claude CLI rejected `--no-subagents` as an
   unknown option before contacting the model. The current CLI help was then
   inspected and the invocation was corrected to `--allowed-tools`.
2. A corrected read-only invocation with `--permission-mode plan`,
   `--permission-prompts none`, and `Read,Grep,Glob` was allowed 300 seconds;
   it exited with no response.
3. A compact read-only invocation with the complete plan facts in the prompt
   and no file-reading dependency was allowed another 300 seconds; it also
   exited with no response.

No invocation was allowed to edit files or run jobs. No Fable judgment is
inferred from the timeouts.

## Operator disposition

The bundle is not treated as Fable-approved. Codex proceeded only with the
pre-authorized bounded work: source-level wiring assertions, the existing
consolidated native/production gate, and removal of the demonstrably dead
`gpu_sink` feedback autotune block. The no-auto-selector decisions remain
explicitly conditional on future equivalent backend implementations and are
subject to the Codex end audit.

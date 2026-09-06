# Claude Opus 5 P0.4 focused re-audit — 2026-09-02

Auditor: Claude Opus 5 CLI, high effort, read-only.

Scope remained strictly P0.4 stellar-feedback admission control. Generic
HDF5/restart, AMR, ksection/CPU boxes, gravity, RT, dust, AGN, and unrelated
working-tree changes were excluded.

## Verdict

**PASS (P0.4 admission control only).**

The auditor independently recomputed the production binary SHA-256 as
`46e408e60ab289eccffde6b00c3e7f93c4669e4ada25f455448bd2354deedb05`
and found it matched the machine-readable build evidence. It also verified
that all P0.4-relevant source hashes in that evidence match the current files.

## Finding closure

- C1 CLOSED: feedback worker errors cross the OpenMP boundary through per-thread
  slots and only the parent thread calls `MPI_ABORT`.
- C2 CLOSED: production main performs enabled-channel mass-window coverage
  preflight before adaptive integration.
- C3 CLOSED: missing stellar namelist maps to error 1005, invalidates the
  parameter set, and is unit-tested.
- C4 CLOSED: output records the actual loaded path/row count and effective IMF,
  population, channel-enable, element-enable, and mass-window settings.
- C5 CLOSED: the real production binary executes two negative admission paths;
  missing table exits 1 and incomplete coverage exits 121.
- C6 CLOSED: NaN, negative, transactional rollback, missing-group, SNIa, and
  PISN unit cases are present.
- C7 CLOSED: all parsed physical table columns are checked finite before row
  commit.
- C8 CLOSED: the production and native increment drivers both reject nonfinite
  channel windows.

No new P0.4 blocker was found. The remaining generic RAMSES `clean_stop` zero
exit behavior is out of scope and recorded in the supporting-infrastructure
long-term backlog; it cannot admit stellar feedback physics.

## Next mandatory feedback work

The auditor confirmed that the active roadmap prevents the binary population
and fate model, SNIa DTD, and PISN/PPISN population decision from being skipped.
They remain F-P1 through F-P3, followed by integrated channel realization
F-P4. Each receives its own Claude Opus 5 physics/code audit.

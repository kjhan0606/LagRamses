# Claude Opus 5 P0.4 fail-closed runtime audit — 2026-09-02

Auditor: Claude Opus 5 CLI, high effort, read-only.

Scope was restricted to P0.4 stellar-feedback admission control.  HDF5,
restart internals, AMR, ksection/CPU boxes, gravity, RT, dust, AGN, and unrelated
dirty-worktree changes were explicitly excluded.

## Initial verdict

**CONDITIONAL PASS.**  The auditor found no path that silently selected
synthetic yields, clamped an out-of-domain query, or enabled unimplemented
binary/SNIa/PISN physics.  It judged all six P0.4 criteria substantively met,
with hardening and evidence conditions.

## Findings and owner disposition

| ID | Auditor finding | Disposition |
|---|---|---|
| C1 | feedback error called collective `clean_stop` inside OpenMP | accepted; collect worker error and abort from parent thread |
| C2 | enabled-channel/window coverage was deferred until first star | accepted; production startup coverage preflight |
| C3 | absent/misspelled stellar namelist silently used defaults | accepted; group is mandatory and tested |
| C4 | info file re-read the environment and omitted enable flags | accepted; record actual loaded identity, rows, flags, and elements |
| C5 | build smoke did not execute an admission branch | accepted; add production-binary negative runner and narrow evidence wording |
| C6 | NaN/negative/transaction/PISN unit cases missing | accepted; extend production-source unit tests |
| C7 | loader admitted nonfinite non-age columns | accepted; finite check all parsed physical columns before commit |
| C8 | production driver lacked mirror finiteness guard | accepted; add the same guard in both trees |

The auditor also noted that generic RAMSES namelist failures return zero through
legacy `clean_stop`.  The stellar namelist is now rejected, so this cannot admit
feedback physics.  Changing global termination semantics is outside P0.4 and is
recorded in `supporting_infrastructure_long_term_backlog.md`.

The audit explicitly agreed that population normalization/fate ownership, the
SNIa DTD, and PISN/PPISN eligibility are later **active mandatory feedback
physics**, not P0.4 conditions and not omitted work.  It found the roadmap
substantively adequate and requested durable links from the main plan; those
links are now present.

Re-audit after the focused remediation is required before final P0.4 closure.

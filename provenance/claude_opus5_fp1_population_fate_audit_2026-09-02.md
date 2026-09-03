# Claude Opus 5 F-P1 population/fate audit — 2026-09-02

Auditor: Claude Opus 5 CLI, high effort, read-only.  Scope was restricted to
F-P1 population/fate semantics and excluded RT, dust, AGN, SNIa DTD/PISN
implementation, and generic infrastructure.

## Verdict

**BLOCK.**  The scientific self-assessment was judged honest, but the auditor
did not accept the initial `ENGINEERING CONTRACT COMPLETE` headline.

All reported commands, production linkage, source-parity result, and binary
SHA-256 `188ddd9fc698730b58fe4b3736c395c124c42a234e466bd5dc436786a59c6665`
were independently reproduced.

## In-scope findings

| ID | Finding | Auditor severity/disposition |
|---|---|---|
| D1 | The implemented Chabrier low and high branches are discontinuous at 1 Msun; the high branch lacks the continuity amplitude.  This materially overweights massive stars. | Blocking code/physics defect; fix now and replace the self-copying JAX check with an actually independent shape/integral check. |
| D2 | Production `stellar_yield_audit` requires `sum(tracked)==returned`, contradicting the declared `sum(tracked)<=returned` untracked-ejecta residual contract.  Production deposition also omits that residual from generic metals. | Blocking production-contract defect. |
| D3 | `stellar_population_ledger` is linked, but no production runtime caller executes it.  Ownership and over-return guards are therefore unit-test-only. | Blocking implementation defect; do not describe linked symbols as operational enforcement. |
| D4 | Production table audit lacks native duplicate-coordinate, complete-grid, energy-monotonicity, remnant-ownership, and age-zero checks. | Blocking production/native semantic mismatch. |
| D5 | The `feedback_mode='legacy'` namelist branch commits only the mode and silently discards other values parsed in that group. | Minor transactional defect. |

The auditor accepted the explicit source-basis rejection, mandatory IMF
support configuration, SNIa/PISN double fail-closed layers, build evidence,
and the honesty of the physical-source blocker.

## Scientific blockers retained

- No population/fate source is approved with immutable version/license/hash,
  published IMF convention, lifetime/fallback/remnant map, binary assumptions,
  and complete source-basis sidecar.
- The provisional defaults leave 40--120 Msun stars with winds but no terminal
  owner.  For the current Kroupa shape the auditor estimated this as about
  6.8% of initial SSP mass, so it cannot be hidden in the living-mass residual.
- PISN is currently labelled a terminal remnant owner even though a genuine
  PISN leaves no remnant.  It is harmless while fail-closed and must be decided
  explicitly at F-P3 together with PPISN semantics.

## Required report corrections

Until D1--D4 are independently closed, describe the ledger as linked but not
production-exercised, the old JAX result as transcription consistency rather
than an independent shape validation, and mass closure as a definition plus
nonnegativity/ownership checks.  F-P1 remains open.

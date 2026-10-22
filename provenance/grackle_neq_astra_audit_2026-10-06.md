# Astra read-only audit — Grackle NEQ validation

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001`
Audit date: 2026-10-06
Auditor: Astra (`gpt-6-astra`), read-only review; no files changed and no
calculations launched.

## Verdict

- **Declared-scope Grackle integration verification: CONDITIONAL PASS.** The
  bounded verification bundle may be reported as complete after the corrections
  and exclusions below are stated. Do not describe it as unqualified completion
  of every Grackle physics validation.
- **Production qualification: FAIL / not yet qualified.**

## Audit opinion for the completion email

> CONDITIONAL PASS for completion of the declared Grackle NEQ
> integration-verification scope. Evidence supports the tested chemistry
> interface budgets, restart, AMR, MPI/OpenMP consistency and UV response. The
> six-step real-IC run ends at `a=0.01851` and does not reach reionization. This
> is not production qualification: stiff-regime accuracy, actual feedback
> coupling, later cosmological evolution, hydro positivity and throughput
> remain unresolved.

## Evidence and limits

The retained reports, logs and summaries support the adapter comparisons,
element/charge checks, tested MPI/OMP consistency, restart continuation, AMR
behavior, UV response and early cosmological integration within their stated
fixtures. These results do not establish general physical accuracy or parallel
scaling. The six-step cosmological comparison ended at `a=0.01851` (approximately
`z=53`), well before reionization. Its maximum field difference was computed in
HDF leaf-array order; coordinates were not cross-matched. Astra reviewed the
retained evidence after raw HDF5 snapshots had been removed and did not
independently recompute those snapshot statistics.

The following remain separate qualification blockers or explicitly excluded
physics claims:

- Actual stellar/SN ejecta deposition with production-approved yields was not
  exercised; a synthetic inventory restart is not a physical event test.
- Later UV/reionization evolution was not tested.
- Stiff-regime timestep convergence is a substantive chemistry-accuracy gap.
- The initial pre-cooling negative-energy state is reproduced with Grackle
  disabled, but the no-Grackle control is only one step; later deterioration
  cannot be causally separated from prior-step cooling history. Positive final
  dumps do not prove accurate recovery through the temperature-floor path.
- Cooling throughput was measured in one configuration, not production-scale
  scaling. The high cooling fraction is an operational warning, not by itself
  a correctness failure.

The added verification was judged proportional to the distinct interfaces and
the observed hydro diagnostic; broad additional testing is not required to
close this explicitly bounded report, but the limitations above must remain
visible.

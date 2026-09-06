# Claude Opus 5 retry: bounded plan audit

Read-only audit. Do not edit, build, commit, push, launch jobs, or enable
runtime flags. Return a verdict and concise reasons within 1200 words.

Repository: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`, branch `main`.

Read only these relevant materials:

1. `provenance/fp1_agn_ledger_transaction_bundle_plan_2026-09-04.md`
2. `provenance/audit_governance_amendment_2026-09-04.md`
3. `patch/lagRamses/AGN_COARSE_STATE.md`
4. `patch/lagRamses/snrt_agn_source.f90`
5. `patch/lagRamses/snrt_ramses_driver.f90`
6. `simulation/snrt/snrt_core/sink_diagnostic.py`
7. `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
8. `simulation/snrt/tools/reproduce_fable_sn_agn_findings.py` around F10--F14
9. `patch/lagRamses/sink_particle.kjhan.f90`, lines 6130--6330 and
   2170--2430 only

Judge the plan against the final purpose of a production/publication-ready
lagRamses high-level RT, stellar/AGN feedback, dust, and coupled-source
stack. Return exactly one of `PASS`, `CONDITIONAL PASS`, or `BLOCKED`.

Check: (a) whether this is an in-scope, bounded pre-G2 engineering bundle;
(b) whether `(nstep_coarse,sink_id)` duplicate handling and raw versus
effective efficiency semantics are adequate; (c) whether the proposed
all-or-nothing multi-group source transaction actually prevents partial
mutation and retry double counting in the current Fortran wiring; (d) whether
OpenMP/MPI, source-order, and fail-closed implications are covered; and (e)
which acceptance tests or scope limits are missing. Distinguish arithmetic
and transaction evidence from physical AGN-hydro closure. Do not demand
generic AMR/HDF5 work and do not treat the static nine-group ledger as a live
physical AGN approval.

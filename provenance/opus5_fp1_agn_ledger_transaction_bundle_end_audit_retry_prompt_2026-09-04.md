Read-only, no edits/jobs/web.  You are Claude Opus 5, sole auditor.  Give a
verdict exactly PASS, CONDITIONAL PASS, or FAIL in <=1200 words.  Audit the
current `/gpfs/kjhan/LRD_JWST` worktree for the bounded F-P1.5 AGN ledger
transaction bundle.  Project goal is production/publication-ready
lagRamses high-level hydro for RT, stellar/AGN feedback, and dust; this bundle
only closes an AGN source-ledger/transaction prerequisite and must not claim
AGN SED, obscuration, hydro closure, G2 stellar yields, or durable restart
journaling.

Read these exact files:

- `provenance/fp1_agn_ledger_transaction_bundle_plan_2026-09-04.md`
- `patch/lagRamses/snrt_agn_source.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/sink_particle.kjhan.f90` (dump_agn_coarse_state and
  AGN_feedback)
- `bin/Makefile`
- `simulation/snrt/snrt_core/sink_diagnostic.py`
- `simulation/snrt/tools/audit_agn_coarse_ledger.py`
- `simulation/snrt/tools/p4_build_agn_rate_ledger.py`
- `simulation/snrt/tools/p7_convert_sinkprops.py`
- `simulation/snrt/data/agn_coarse_ledger_transaction_audit.json`

Check only these acceptance questions:

1. Does the reader enforce finite required fields, explicit pre-reset
instantaneous markers, 365.25-day conversion, min(Bondi,Eddington),
raw/effective efficiency distinction (effective=0 idle allowed), and
Lbol algebra; does it collapse semantic same-key duplicates and fail closed on
conflicts/rewind ambiguity?
2. Is the Fortran writer called exactly once before AGN_blast and both mass
accumulator resets, with matching markers/fields?
3. Does the transaction prepare every group and commit no intensity state if a
later group fails; does the driver advance accounted_mass only after success?
4. Is accounting keyed by stable idsink and remapped safely after sink-array
reorder, with a serial shared-state source loop, local-leaf/single-MPI-owner
contract, and one-time SNRT_RT_ENABLE latch?
5. Does bin/Makefile link the required SNRT module graph under the existing
SNRT=1 requires USE_CUDA=1 gate?
6. Do the audit JSON limitations correctly preserve the known open issues:
no run UUID/dump counter, live raw eps_sink versus effective ledger mismatch,
cross-step deferred Esave re-emission, and no durable crash journal?

Report closed criteria with short file/line references, blockers, deferred
issues, and whether this bounded engineering bundle is safe to close while
runtime activation and publication claims remain gated.  Do not treat the
synthetic fixture as a physical run.

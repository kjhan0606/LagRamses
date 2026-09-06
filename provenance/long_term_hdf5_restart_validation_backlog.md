# Long-term HDF5 restart validation backlog

Date: 2026-09-02  
Project: `/gpfs/kjhan/LRD_JWST`

This backlog separates general lagRamses checkpoint hardening from the
stellar/AGN feedback implementation gates.  The items below are useful for a
future production restart qualification, but they do not block the focused
P0.3 stellar-feedback state contract.

## Existing implementation baseline

- `cd3d64b`: variable-CPU HDF5 decomposition and ksection gate.
- `4cc2e2b`: HDF5 particle restart free-slot handling.
- `255d2cb`: modified-gravity scalar checkpoint state.
- `77402d3`: output-clock resynchronization after restart.

These capabilities predate the current feedback work and must not be counted
again as new feedback deliverables.

## Deferred validation

1. Validate exact rank and extent for generic hydro `uold_*` and gravity
   `phi`, `f_*`, `scalar_gr`, and optional `psi_re/im` datasets, including
   malformed-checkpoint negatives.
2. Replace rank-local restart `clean_stop` paths with an MPI-safe abort policy
   and test asymmetric failures.
3. Check every optional HDF5 existence-query status before using its result.
4. Add runtime coverage for ADM checkpoint payloads.
5. Improve per-rank failure diagnostics for collective HDF5 reads.
6. Requalify ksection/CPU-box and same-/different-CPU AMR redistribution at
   production-representative MPI counts.

## Feedback-gate boundary

P0.3 is limited to persistence and continuation of the state required to
apply stellar feedback exactly once: particle type, birth/release time,
initial stellar mass, release cursor, and stellar mass ledgers.  Sink/BH
state belongs to the later AGN live-coupling gate.  Generic hydro, gravity,
AMR decomposition, and checkpoint-corruption qualification remain in this
backlog unless a focused feedback test demonstrates a direct dependency.

The broad findings in P0.3 audits 9--15 remain useful historical engineering
notes, but their generic restart findings are not feedback gate blockers.

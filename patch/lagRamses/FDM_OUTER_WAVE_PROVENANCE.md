# FDM outer-wave provenance output

Set `fdm_outer_ledger=.true.` in `&fdm_params` to write one
`fdm_outer_wave_provenance_<output>.txt` file in each normal FDM output
directory.  The option is off by default.

The current V3 file is a compact raw diagnostic.  It records the output epoch,
axion mass, effective code `hbar`, HJM/wave seam settings, dual-soliton seed
switch and parameter values, total leaf-cell wave mass, and the global FDM
mass-current integral.  On wave levels the current uses
`hbar Im(psi* grad psi)`; on HJM levels it uses `rho grad(S)`.  Only cells with
a complete same-level central stencil enter the current integral.  The file
therefore reports the stencil coverage fraction explicitly.

V3 additionally records `mpi_ncpu`, the exact number of rank-suffixed
`fdm_<output>.outNNNNN` and matching AMR shards expected for the output.  A
consumer of a time-resolved wave diagnostic must require every suffix from
`00001` through `mpi_ncpu`; discovering one matching pair is not sufficient.
It also records `restart_parent_output` as a restart hint.  This field does
not provide a run UUID or prove an uninterrupted restart branch, so outputs
from separate executions must remain conditional unless a stronger solver-side
lineage attestation is added.

The dual-soliton fields are runtime configuration provenance, not a claim that
the supplied coherent state has relaxed.  They allow the outer analysis to
bind a V2/V3 output to its materialized two-core seed after the separate input
preflight has checked the complete namelist and `ic_sink` rows.  The Python
reader continues to accept legacy V1 records, but those lack this runtime
attestation and cannot alone verify a controlled dual-soliton initialization.

This output never changes `psi`, gravity, sink evolution, or any drag force.
It records `analytic_fdm_drag_enabled = .false.` and
`force_accounting = resolved_wave_only`; it must not be used to justify adding
an analytic FDM drag term to a calculation with a resolved wake.

It is not a calibrated kpc-to-pc delay, radial FDM closure, or complete wave
ledger.  The per-rank `fdm_<output>.out*` snapshots remain mandatory for the
time-resolved radial profile, local current, core centre, dipole/quadrupole,
granule, seam, boundary, and force-ledger postprocessing required by the
pure-FDM handoff contract.  If that provenance cannot be reconstructed, the
event is censored rather than assigned zero delay.

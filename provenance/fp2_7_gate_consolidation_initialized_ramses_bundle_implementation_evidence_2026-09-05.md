# F-P2.7 implementation evidence — incomplete runtime qualification

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Checkpoint: `55b1d772173b70b62ba0d215d73f46988db91f96`; subsequent edits remain
in the worktree. Status: **partial evidence; D4 and bundle-end audit pending**.

## Recorded native results

The completed bundle gate summary is
`simulation/snrt/build/snrt_bundle_gate_summary.txt` (ignored runtime output).
It records `SNRT_BUNDLE_GATE_PASS`, a 170.308-second SNRT/CUDA build,
four linked symbols, thermochemistry, spectral/checkpoint, two-rank MPI
transaction and CUDA multigroup checks, and production-negative execution.
Negative families report 4 thermochemistry, 10 spectral, 11 transaction/config
cases and the production/SNIa rejection checks. Conservation thresholds are
`2e-6` for the CUDA photon budget and `1e-13` for spectral group sum versus
Lbol. Driver source-route checks are labelled `STATIC_SUPPORTING_CHECK`.

The current binary hash was independently reread during entry reconciliation:

```text
37bc576d44f540dc5b07d495458a702614f2aa9e4efd0f0ca4ec1c45320e75eb  bin/ramses_final3d
```

D3 changed the four mirror-building runners to compile stellar modules from
`patch/lagRamses` and moved eleven test programs into
`simulation/snrt/tests/fixtures/phase0`. The two other callers use the moved
fixture paths. The remaining mirror modules were retained. Prior execution
in this work session recorded G1 native/JAX agreement, G2 population-ledger,
F-P2 SNIa native contract checks and F-P1.2 transaction success. G2 preflight
returned its expected `G2_PREFLIGHT_BLOCKED` physical-admission result.

The CPU P0 linked harness passed and G1 passed normal mode at that intermediate
CPU binary. The subsequent SNRT/CUDA build replaced `bin/ramses_final3d`.
Consequently, the current P0 validator correctly reports
`STELLAR_SOURCE_PARITY_BLOCKED blocked=production_linked_build_evidence`.
Do not claim simultaneous current P0 CPU and SNRT/CUDA binary certification.
Repeatedly rebuilding the shared binary to alternate these two checks does
not close this evidence-identity conflict; qualification needs distinct
retained build identities when that workflow is next changed.

## D4 runtime evidence remains outstanding

Slurm job `332448` failed before runtime validation because the submitted
script was executed from a Slurm spool copy and derived `BINARY=//bin/ramses_final3d`
from `BASH_SOURCE`. This was an infrastructure failure, not a RAMSES physics
result; no baseline or injected case ran and no output directory was created.
The runner was corrected to use the immutable `/gpfs/kjhan/LRD_JWST` project
root and run directory explicitly, then shell-checked and resubmitted as job
`333139`. It is currently `PENDING` for scheduler `Priority`, with one node,
two MPI tasks, two GPUs on partition `a10`, and a ten-minute limit. No D4
runtime pass or injected rollback is evidenced yet.

The intended template and runner are under
`simulation/snrt/runs/fp2_7_initialized_ramses_smoke/`. The template specifies
hydro, legacy feedback mode, a fixed level 3 grid, two coarse steps, disabled
PIC/Poisson and a future output schedule. These are intended settings only;
the effective per-case namelist cannot be attested before job execution.

Actionable scheduler disposition: retain job `333139` and recheck its runtime
logs after start. Completion requires direct baseline and injected runtime
logs, executable identity, and confirmation that no unexpected dump was
produced. A pending job or static marker cannot satisfy D4.

## Outstanding closure

F-P2.7 is not marked complete and has no bundle-end verdict. The provenance
index and this single evidence record carry the pending state forward.
The next bundle's operator approval is recorded separately; it does not turn
this pending runtime evidence into a pass. The path-repair patch and
resubmission are recorded in the current worktree; D4 remains open until the
job produces both required runtime cases.

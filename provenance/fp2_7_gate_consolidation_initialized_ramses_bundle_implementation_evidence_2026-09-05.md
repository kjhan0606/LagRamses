# F-P2.7 implementation evidence — initialized RAMSES runtime qualification

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Checkpoint: the D4 runner and evidence update are being recorded on top of
`d6e677c`. Status: **D4 passed; physical source admission and live production
qualification remain outside this gate**.

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

## D4 initialized-RAMSES runtime evidence

Slurm job `332448` failed before runtime validation because the submitted
script was executed from a Slurm spool copy and derived `BINARY=//bin/ramses_final3d`
from `BASH_SOURCE`. This was an infrastructure failure, not a RAMSES physics
result; no baseline or injected case ran and no output directory was created.
The runner was corrected to use the immutable `/gpfs/kjhan/LRD_JWST` project
root and run directory explicitly, then shell-checked and resubmitted as job
`333139`. It reached the binary but failed during PMI2 initialization; this was
replaced by the Intel `mpirun` launcher in the runner. Job `333201` then
demonstrated the runtime behavior but exposed a harness-only expectation bug:
diagnostic fail-closed returns normally with status 0. The runner now treats
the exact runtime markers as authoritative for that diagnostic case.

The intended template and runner are under
`simulation/snrt/runs/fp2_7_initialized_ramses_smoke/`. The template specifies
hydro, legacy feedback mode, a fixed level 3 grid, two coarse steps, disabled
PIC/Poisson and a future output schedule. The corrected launcher records the
Intel MPI launcher and runs from each case directory, so RAMSES receives the
case-local effective namelist.

Job `333211` is the completed D4 qualification. Direct case evidence is under
`simulation/snrt/runs/fp2_7_initialized_ramses_smoke/job_333211/`:

* baseline: `D4_CASE baseline status=PASS return_code=0`, with
  `SNRT_RT_TRANSACTION_COMMIT_PASS`, `SNRT_RT_CLOSURE_PASS`, and `Run completed`
  in `baseline/ramses.log`;
* injected receiver failure: `D4_CASE injected status=PASS return_code=0`,
  with both `SNRT RT transaction rollback: class=receiver` and
  `SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver` in `injected/ramses.log`;
* launcher: Intel `mpirun`, two MPI ranks, two A10 GPUs, and recorded
  binary/effective-namelist SHA-256 files in each case;
* output safety: no `output_*` directory exists in either case; the runner
  completed with `D4_PASS job_id=333211`.

This closes the initialized-runtime D4 gate only. It does not qualify live
stellar, AGN, dust, or physical high-mass yield coupling, and it does not
resolve the retained CPU-versus-CUDA binary identity caveat above.

## Bundle disposition

F-P2.7 D4 is **PASS for initialized SNRT runtime behavior**. The bundle is not
publication-ready: physical source admission, live multi-source feedback, and
the retained build-identity separation still require their explicitly scoped
follow-up gates. The next implementation bundle is G3/G4 physical source/SED
and dust runtime qualification, subject to the available canonical data.

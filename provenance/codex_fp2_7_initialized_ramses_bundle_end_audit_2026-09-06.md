# Codex end audit — F-P2.7 initialized RAMSES qualification

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Scope: D4 initialized RAMSES SNRT runtime qualification only.

## Verdict

**PASS within the declared D4 scope.** The gate now has direct two-rank Slurm
runtime evidence for both the normal and diagnostic-failure cases. This is an
engineering runtime qualification, not approval of a physical stellar/AGN/
dust production model.

## Evidence checked

Slurm job `333211` completed with exit code `0` on the `a10` partition using two
MPI ranks and two GPUs. The effective namelist was copied into each new case
directory and its SHA-256, the template hash, executable hash, compiler, and
launcher were recorded.

* The baseline case recorded `D4_CASE baseline status=PASS return_code=0`,
  `SNRT_RT_TRANSACTION_COMMIT_PASS`, `SNRT_RT_CLOSURE_PASS`, and `Run
  completed`.
* The injected receiver case recorded
  `SNRT RT transaction rollback: class=receiver` and
  `SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver`. Its return code was zero by
  design: diagnostic mode records the fail-closed result and exits normally.
* Neither case created an `output_*` directory.

The prior runner failures were correctly classified as infrastructure/harness
issues: the first used an invalid spool-derived path, the second used an
incompatible PMI2 launch, and job `333201` exposed the incorrect nonzero-exit
expectation. These were repaired without weakening the production transaction
path.

## Residual conditions

This verdict does not close physical yield/SED admission, live stellar/AGN/dust
coupling, restart equivalence, distributed-AMR qualification, or publication
convergence. The shared executable still carries the documented CPU-versus-
CUDA build-identity caveat; a future release gate must retain distinct build
identities. D4 also does not authorize a large RAMSES production run.

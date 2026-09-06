# G5 SNRT hydro-state restriction bundle — implementation evidence

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)

## Scope

This bounded G5 wiring repair closes one live-coupling ordering hole. The
normal hydro `upload_fine(ilevel)` runs before the SNRT receiver. SNRT can then
commit a new thermal/chemistry state into `uold`; without a second restriction,
the parent AMR level could retain the pre-SNRT state. The repair adds a
requested-SNRT-only restriction immediately after `snrt_ramses_advance_level`
and before the diagnostic hook.

This is an engineering repair. It does not promote physical stellar yields,
SEDs, dust, SNIa/PISN, or a large RAMSES run.

## Changed path

- `patch/lagRamses/amr_step.jaehyun.f90`: import the existing rank-uniform
  `snrt_agn_rt_requested()` latch and call `upload_fine(ilevel)` after a live
  SNRT step when hydro is enabled.
- `simulation/snrt/runs/fp2_7_initialized_ramses_smoke/run_fp2_7_initialized_ramses.sbatch`:
  accept an explicit `D4_BINARY` override so a newly qualified binary can be
  tested without overwriting the shared historical binary.
- `tests/run_snrt_hydro_restriction_unit.sh`: compact source-order regression;
  this is supporting evidence, not a live-runtime gate.

## Evidence

```text
SNRT_HYDRO_RESTRICTION_WIRING_OK receiver_line=946 post_restriction_line=953 diagnostic_line=956
```

The production SNRT+CUDA link was rebuilt with `SNRT=1 USE_CUDA=1` from the
current `/gpfs` tree. The qualification binary is:

```text
/gpfs/kjhan/LRD_JWST/simulation/snrt/build/g5_wiring_qualification/ramses3d
sha256 eebc442fef604cf628c0d320df26afd17353b49f228749c319b3fc1a07061a38
```

The initialized-RAMSES two-rank GPU smoke was submitted as Slurm job `333253`
using that binary in a new job directory. At record time it was pending queue
allocation; baseline and injected-case results are therefore intentionally not
claimed here.

## Limits and disposition

The repair is eligible for the next initialized-RAMSES runtime qualification
after job `333253` completes. The current dust driver remains
`ZERO_SCAFFOLD`, physical G2 source admission remains blocked, and this bundle
does not authorize production evolution or publication claims.

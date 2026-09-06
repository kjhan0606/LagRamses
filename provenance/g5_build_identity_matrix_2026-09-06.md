# G5 CPU/SNRT build identity matrix

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)
Purpose: preserve independent CPU production-link and SNRT/CUDA qualification
artifacts so one build cannot silently invalidate the other.

## Matrix

| Role | Binary | Build contract | SHA-256 | Status |
|---|---|---|---|---|
| CPU production-linked stellar/P0 | `/gpfs/kjhan/LRD_JWST/bin/ramses_final3d` | `P0_BUILD=1 bash simulation/snrt/tests/run_p0_production_linked_contract.sh` | `5eaef4bdc321abbe2494a481d229b7170221c0be623d572f722b3f326253bff2` | `STELLAR_SOURCE_PARITY_PASS`; P0 smoke pass |
| SNRT + CUDA reference/control | `/gpfs/kjhan/LRD_JWST/simulation/snrt/build/g5_startup_contract_admission/ramses3d` | `make -C bin EXEC=... SNRT=1 USE_CUDA=1 ramses` | `e30a938397781ee02eb060222937007f659bed4a0e5d268c1cee125901e019f0` | native contracts pass; D4 Slurm rerun pending |

The two artifacts are intentionally retained at different paths. The CPU
build was regenerated after the D2/D3 consolidation commit; the SNRT/CUDA
artifact was built and recorded under the G5 startup-contract qualification
directory. No claim is made that the two binaries are interchangeable.

## Verification

The CPU validator and its regression test returned:

```text
STELLAR_SOURCE_PARITY_PASS blocked=none
STELLAR_SOURCE_PARITY_GATE_OK status=pass differing_shared=11
```

The CPU no-argument startup smoke returned the expected RAMSES namelist
diagnostic. The SNRT/CUDA binary is the one submitted to D4 Slurm job
`333272`; that job remains a pending infrastructure run and is not counted as
passed here.

The current checkout contains unrelated in-progress worker changes, so this
matrix is a qualification record rather than a clean-release claim. Before a
physical production run, freeze the complete source tree, rebuild the selected
binary from that clean commit, and recalculate this matrix.

## Boundary

This closes the CPU-versus-CUDA evidence identity ambiguity. It does not
promote physical stellar yields, SEDs, SNIa runtime activation, AGN SED
selection, or live dust. The bounded reference/control path remains the only
simulation path admitted by the current contracts; physical dusty feedback
still requires the recorded G2--G5 assets and live receivers.

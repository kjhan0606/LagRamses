# G5 simulation-ready preflight

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)

The read-only launcher preflight is
`simulation/snrt/runs/simulation_ready_preflight.sh`. It submits no job and
creates no run directory. The caller must select a simulation class
explicitly.

## Results

Using the initialized-RAMSES effective namelist
`simulation/snrt/runs/fp2_7_initialized_ramses_smoke/effective.nml`:

```text
reference-control -> SIMULATION_READY_PASS
legacy-comparison -> SIMULATION_READY_PASS
physical-production -> SIMULATION_READY_BLOCKED reason=physical_production_manifest_blocked (exit 2)
```

The reference-control path verifies the retained SNRT/CUDA binary, the
reference nine-group contract, the Furlanetto--Stoever contract, hydro/legacy
namelist settings, and the intentional future output schedule
(`noutput=1`, `aout=2.0d0`, `tout=1.0d30`, `foutput=fbackup=1000000`).
The legacy path verifies a CPU binary without SNRT symbols. Both paths print
binary/namelist hashes and the repository identity for the launch record.

The physical path is intentionally not bypassable: the current production
manifest remains `blocked_until_all_required_assets_are_approved`, and the
live dust driver remains `ZERO_SCAFFOLD`. This preflight therefore establishes
simulation readiness for bounded reference/control and legacy-comparison
calculations, not physical dusty production.

The post-build identity repair also makes the P04/SNIa production-negative
runner accept an explicit `P04_BINARY`/`SNRT_PRODUCTION_BINARY` path. The
SNRT/CUDA bundle gate now passes its freshly linked binary directly, so the
negative test cannot accidentally inspect the separately retained CPU
`bin/ramses_final3d`. The override path was exercised with the retained G5
SNRT/CUDA binary and returned `P04_PRODUCTION_NEGATIVE_OK`.

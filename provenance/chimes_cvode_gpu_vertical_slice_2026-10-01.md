# CHIMES CVODE GPU vertical slice (experimental, not production qualified)

Workspace: `/gpfs/kjhan/chimes-cvode-gpu-20261001`, branch
`agent/chimes-cvode-gpu-20261001`, parent lagRamses `9bec08c`.
Original repository and `/gpfs/kjhan/chimes-cuda-realrhs.1PGaf0` are untouched.
No commit or push. No simulation snapshots produced.

## Dependency finding and scope

Actual prefix: `/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release`.
Headers identify SUNDIALS **5.8.0**, DOUBLE, 64-bit indices. It provides CPU
CVODE, serial vectors, dense matrices and dense/SPGMR solvers, but no
`libsundials_nveccuda`. No system module or loader-cache installation was found.
The job supplies the prefix explicitly and prints `ldd` and library hashes.

Source snapshot is CHIMES `a58e5c0` plus the six existing local modifications
from `.chimes-transition.YdVRvD/chimes`. `YdVRv0` in the request did not exist.
The saved `simulation/snrt/data/chimes_cvode_gpu_baseline.patch` is **inherited
work**, not new GPU implementation. All edits are in the dedicated clone
`.cvode-gpu/chimes`. The additional factory patch does not change CHIMES structs.

The official SUNDIALS GPU model leaves solver control on the CPU and puts
vector/RHS/linear operations on device. A CUDA vector cannot simply replace
CHIMES' serial vector: its RHS, roots and state extraction use `NV_Ith_S`,
host arrays and host callbacks; CPU `SUNLinSol_Dense` rejects CUDA vectors.
See the installed SUNDIALS 5.8 source `src/sunlinsol/dense/sunlinsol_dense.c`,
and [SUNDIALS GPU documentation](https://sundials.readthedocs.io/en/v7.6.0/sundials/GPU_link.html).
The [official batched example](https://github.com/LLNL/sundials/blob/main/examples/cvode/cuda/cvRoberts_block_cusolversp_batchqr.cu)
shows concatenated systems under a single CVODE. Adopting that here would
share step/order/error control across cells, unlike today's independent solves.

For a full device path, separately build the pinned SUNDIALS with CUDA enabled
(`ENABLE_CUDA=ON`, matching double/index settings, supported CUDA compiler),
and implement device RHS, callbacks and a compatible block linear solver.
Do not assume this dependency exists or switch integrators to bypass it.
In particular, the installed 64-bit index ABI is incompatible with SUNDIALS
5.8's cuSPARSE/cuSOLVERSP batch-QR backend (requires 32-bit indices per
`cmake/SundialsBuildOptionsPost.cmake`). That backend would need a separate
consistent rebuild of dependencies, or a different supported linear solver;
silently changing index types in the current libraries is not valid.

## Implemented slice

Keep per-cell CVODE BDF and its original numerical Jacobian. Coalesce pending
dense LU setups through a custom **SUNLinearSolver backend**, using cuBLAS
`DgetrfBatched`. This is not a custom ODE integrator. Species, thermal and
photon RHS and triangular solves remain on the CPU. No chemistry tables are
uploaded by this slice. This tests dependency/production dispatch and whether
the available OMP callback batch can profitably offload any CVODE work.

Build `CHIMES=1 USE_CUDA=1 CHIMES_CUDA=1`, using the factory-patched CHIMES
library. Config `SNRT_CHIMES_CVODE_BACKEND=cpu` (or unset) retains the original
CPU path. `cuda_batched_lu` explicitly enables this experimental backend.
An unknown setting, missing build option, or unpatched CHIMES library rejects
initialization. No namelist key has changed; mkrun/nml generators need no edit.

At most the active OMP team can arrive concurrently, not all `nleaf` cells.
The first arrival waits up to 50 microseconds for peers; groups are split by
matrix dimension and capped at 64. Fewer than 4 compatible calls use CPU LU.
There is no barrier requiring all cells to arrive, so unequal work, retries,
network changes or the final partial wave do not deadlock.

Each eligible batch tries the existing process-local stream pool once.
An unavailable lease immediately uses CPU. The pool selects the device using
MPI local rank; this backend neither changes that policy nor gathers across
MPI ranks. Exactly one batch coordinator per rank is active. The pinned and
device workspace is reusable per leased slot and bounded by the largest
observed group (maximum 64 matrices of order 256 per slot). Stream scratch is
private; DMA completes before the lease is released. No claim of cross-rank
GPU load balancing or stream-count optimization is made.

The original host matrix is unchanged until a successful GPU factorization
is copied back. CUDA failure, singular LU or invalid pivot data causes CPU
factorization from that unchanged matrix and the original SUNLS status is
returned. cuBLAS one-based pivots are explicitly converted to SUNDIALS
zero-based pivots. Each caller owns its CVODE state and matrix throughout.

Existing bridge retries, positivity/element/electron reconciliation, thermal
roots and return codes remain in control. The Fortran OMP loop, MPI collective
rejection and later all-or-nothing trial commit have no edits. This preserves
their ownership structurally; an actual MPI RAMSES integration has not yet
been tested and must not be inferred from a native bridge test.

## Reproduce and evaluate

Submit from this worktree:

```
sbatch simulation/snrt/tools/chimes_cvode_gpu_test.sbatch
```

The script uses H100/H200, one GPU, eight CPUs, 16 GiB, at most 15 minutes,
and a 240-second runtime cap. It compiles the private source and real stream
pool and bridge, prints runtime dependencies and hashes, then tests:

- FP64 pivoted dense LU, singular status and OMP1 CPU fallback;
- 24 real CHIMES dark cells with mixed molecular fraction, T, nH and CR rate;
- CPU/GPU abundance, temperature and thermal-energy parity, positivity,
  nucleus and charge conservation, and native bridge return status;
- occupied-stream fallback and rejected-input preservation.

Species comparison allows `1e-17 + 1e-6*abs(CPU abundance)`; temperature and
thermal energy relative error <1e-6; nucleus relative error <1e-8 and absolute
charge <1e-10. These comparisons do not loosen CHIMES integration tolerances
or the bridge's acceptance criteria.

The dense benchmark explicitly synchronizes arrivals and is synthetic. The
real CHIMES test has **no per-callback barrier**, so actual callback grouping
is measured. Table load, pool initialization and cold LU timing are printed
separately; steady timings include gathering, packing, transfers and copies.

Photochemistry hooks deliberately abort if called in this dark-only test.
Photon variables, illuminated bands, event-root transitions, MPI commits,
multiple ranks/cards and full production speedup are **not verified** here.
This slice alone does not complete full CHIMES GPU acceleration. In particular,
the measured CPU cost was dominated by RHS; a fast LU cannot remove that cost.

## Measured result / disposition

Job **409496** completed 0:0 on **syn09 / NVIDIA H100 NVL**, 12 s allocation.
Log: `.cvode-gpu/job-409496.log`; executable and private library:
`.cvode-gpu/build-409496/{test,libchimes.so}`. Actual loader paths in that log
confirm the release SUNDIALS above and CUDA/cuBLAS 12.8.1. Test SHA256:
`bf853fe3eeac550553639abb2ab577b2c1d84a3f96d984ed194a805a7d9318f5`.

| Measurement | CPU OMP8 | Experimental GPU | CPU/GPU speedup |
|---|---:|---:|---:|
| Synthetic LU, 100 rounds, eight matrices/round | 0.0896365 s | 0.194321 s | 0.4613 |
| Real CHIMES bridge, 24 manufactured dark cells | 0.0120458 s | 0.0366036 s | 0.3291 |

The real chemistry comparison actually performed **814 GPU LU setups in
128 GPU batches**, with 60 CPU LU fallbacks; maximum batch 8. No device
errors were recorded. Pool saturation forced 874 CPU setups, of which 831
belonged to otherwise eligible batches; the remainder were too small. This
fallback run took 0.017636 s and reproduced the reference results.

There were **22 successful cells**. Their 157 abundances, temperature and
derived thermal energy matched the CPU results exactly in this test.
Maximum nucleus relative error was 2.2204460492503131e-16 and absolute net
charge 2.3417289710648721e-20. Synthetic pivoted LU maximum solution error
was 5.773159728050814e-15; singular LU and OMP1 fallback checks passed.

**Two cells (indices 5 and 17) failed in both CPU and GPU modes with native
bridge status 41**, and both modes left their trial outputs at the original
inputs. These share nH=0.11, T=1e4 K, initial H2 abundance 0.05, dt=1e10 s.
This is a retained CPU-baseline RHS/integration limitation, not a successful
evolution test for those states. No tolerance or physics was changed to hide
it. The pass marker refers to parity of 22 accepted states plus identical
rejection for two states, not success of all 24 physical integrations.

Startup measured table loading **0.0975583 s**, CUDA pool initialization
**0.336941 s**, and the first GPU LU round **0.0229404 s**. Steady GPU timings
follow those initializations and include queue wait, packing, H2D/D2H and
factorization. No chemistry-table GPU upload exists in this slice. This is
a short one-node bounded comparison, not a full simulation performance claim.

**Do not promote this backend for production performance.** It is correct
on this bounded accepted-state sample but approximately 3.04x slower for the
actual callback workload. Default CPU is unchanged. An independent-cell
host batch is safe for transaction ownership, but opportunistically gathering
only OMP8 LU calls does not expose enough work to amortize device transfers.

Full acceleration remains blocked on a resident **RHS + cooling + batched
solver** data path. A larger host batch would require an explicit gather/
solve/scatter interface at the tile/leaf level; preserving per-cell adaptive
error control versus a common CVODE step must be decided explicitly. The
large synthetic thermal interpolation speedup from the separate agent cannot
be used as evidence that this production callback path is fast.

Failure history is retained: 409492 compilation (`RTLD_DEFAULT`, corrected
to `dlopen(NULL)`); 409493 native test linkage (stream-pool cleanup hooks for
unused MG/scalar/particle state; isolated test now supplies no-op finalizers);
409494 rejected the two failing CPU baseline cells. 409496 explicitly checks
matching failure statuses/unchanged outputs without labeling them successful.
The test-only finalizers are not used by the production Makefile.

## Changed files and integration boundary

- `bin/Makefile`: optional build and symbol export, test linkage, rebuild
  bridge on config changes; VPATH order is unchanged.
- `patch/lagRamses/snrt_chimes_bridge.c`: explicit config and factory-ABI
  check before chemistry initialization; default CPU remains available.
- `patch/lagRamses/snrt_chimes_cvode_cuda.{cu,h}`: batch LU adapter and lease
  lifetime, CPU fallback, bounded workspaces.
- `simulation/snrt/data/chimes_cvode_gpu_{baseline,factory}.patch`: inherited
  source state and minimal new external-CHIMES factory hook respectively.
- `simulation/snrt/tests/fixtures/phase0/chimes_cvode_gpu_test.cpp` and
  `simulation/snrt/tools/chimes_cvode_gpu_test.sbatch`: bounded evidence.
- This report.

Native bridge/config/backend compiled and ran; a full lagRamses binary, MPI
all-or-nothing transaction, illuminated 9-group photons, band/dust-relative
paths and thermal-root events were **not** run. The production wiring is
implemented as an opt-in experiment, not qualified for merging/deployment.
There are no dumps or raw simulation outputs to delete; compact logs remain.

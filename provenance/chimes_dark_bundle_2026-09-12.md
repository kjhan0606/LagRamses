# CHIMES dark integration: audit synthesis and bounded implementation

Scope: operator-approved P0/P1 bundle, not new physical completion gates.
Project: `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.

## Independent audit synthesis

Reports are retained in `.chimes-parallel-audit.Yugsr1/`: Fable, Grok 4.6,
Sol (gpt-5.6-sol), AGY (gemini-3.8-flash-high), Astra (gpt-6-astra).
All completed. AGY only reviewed supplied excerpts because headless file
permissions failed; its coverage is narrower. Astra identified a concrete
accumulator initialization defect and unnecessary explicit candidate work.
Fable favored solver counters/Jacobian before workspace reuse; Grok/Sol/AGY
favored reuse. Neither negligible allocation cost nor large reuse speedup
has been measured. Worker elapsed/wall ~3.9 is NOT OMP1/4 scaling or CPU
utilization. CUDA whole-solver rejection is not established by IR timings.

Priority: (1) deterministic accumulation, (2) exact redundant dark work,
(3) setup/RHS/linear-cost evidence before workspace/Jacobian investment,
(4) scheduling only if imbalance warrants it, (5) bounded GPU-resident FP64
batch pilot if warranted. No stream-count tuning, tolerance relaxation,
hot-mode spectrum reduction, or new routine audit campaign.

## Implemented

- `update_rate_vector` resets ALL species creation/destruction accumulators
  before reaction tables add to them. Inactive species also receive additions.
  Removed compact-index explicit initialization and redundant active-only RHS
  clearing. Sparse active mapping is still used for state/ODE indexing.
- If `explicitTolerance <= 0`, skip the impossible explicit candidate, not
  the implicit solve. Initial constraint checks and coefficient construction
  remain. Positive-tolerance explicit behavior is preserved apart from the
  corrected deterministic accumulation.
- Dark secondary budget returns after setting its outputs to exact zero.
  Molecular photon coefficients retain zero initialization but skip photon
  shielding interpolation at `N_spectra=0`. Molecular cooling columns,
  collisional dissociation, case B, roots, retries and admission are unchanged.
- Actual external source and `simulation/snrt/data/chimes_native_receiver.patch`
  agree (`git apply --reverse --check` passes). Pinned upstream is
  `a58e5c0311993b51abc63d84fff2958a0104f6d6`; receiver ABI remains 6.
- Existing mixed probe promoted to a tracked fixture with direct OMP4/serial
  state, temperature, elapsed-time and status comparison. Added a tiny native
  poisoned-accumulator regression; these are not a new test framework.
- Namelist semantics and Makefile VPATH order unchanged.

## Native evidence

Logs and binaries: `.snrt-performance.jDx9Iz/bundle-*`.

- Neutral 9-state and mixed 6-state before/after printed double-precision
  states: maximum temperature relative and species/H absolute difference 0.
  These comparisons do not establish universal trajectory equivalence.
- Mixed OMP4 vs serial directly compares all outputs and passes.
- Existing thermochemistry smoke passes, including accepted upward root,
  no event after requested endpoint, and atomization/atomic remainder.
- Existing long-interval photo/CMB regression passes. Initial invocation
  omitted its mandatory hex fixture argument and failed before physics;
  corrected invocation uses `chimes_photo_long_interval.hex`.
- Poisoned sparse mapping (active indices 0,2, old compact initialization
  clears 0,1): old retained library FAILS at species 2; corrected library
  PASSES twice. This demonstrates the initialization contract defect, not
  evidence that previous cosmological outputs were corrupted. The initial
  test harness hit Intel libimf/glibc dlopen ordering; linking libm explicitly
  with `--no-as-needed` resolves it before the native comparison.

Old library retained at `.snrt-performance.jDx9Iz/libchimes.bundle-before.so`.
The earlier binaries dynamically load the external library, so historical
logs are the baseline unless the retained old library is explicitly selected.

## One-step integration

Job 542031, grammar[109-112], MPI32 x OMP4, CPU only, 12-minute cap.
Inputs: `.snrt-performance.jDx9Iz/bundle-profile/run.nml` and `run.sbatch`.
Same effective physics as 542016; only stale two-step comment corrected.
This node allocation differs from grammar[031-034]; wall-time comparison is
not a controlled same-node scaling experiment.

RAMSES SHA256: `092087a09b19b24b7c3229b99bd8ebe04c4e38f48b587c32fdc0a38c39df9c6d`.
libchimes SHA256: `f491790e1dd2705bb1248f56c6bd288b16e06788bd35035aabd4faa1e13ecb9a`.
NML SHA256: `ad63a14662da74be9a89c73dad73d2b2bd9c80ffa5bc5f2c5d62816be1494bc3`.

No scheduled dumps: nstepmax=1, noutput=1, aout=1.1, tout=1e100,
foutput=fbackup=1000000. Expected raw output 0; free space at launch 91 TiB.
Job COMPLETED, exit 0:0. RAMSES wall 325.976843 s, allocation 337 s,
TotalCPU 08:20:25 (8.3403 core-hours). All 32 CMB commits and one IR commit
pass; reported mcons=econs=0; no chemistry/IR rejection. No `output_*`
directory exists: raw simulation output count 0, deleted 0. Inputs, logs,
binary and old/new library identities retained.

Compared with job542016: cold wall 197.777 -> 193.346 s (-2.24%); summed
dark worker 776.141 -> 760.636 s; TotalCPU 08:23:17 -> 08:20:25 (-0.57%).
Overall wall 312.5697 -> 325.976843 s (+4.29%). Primary decision collective
7.674 -> 31.863 s; IR wall 65.828 -> 59.945 s. Different allocated nodes
and variable collective timing preclude causal overall speedup claims.
This bundle closes correctness/redundant-work fixes, NOT the remaining
large chemistry cost. Native timing improvements are small and unreplicated.

Next meaningful performance decision remains a bounded split of solver
setup, RHS/Jacobian and linear-solve costs using existing probes/counters.
Do not infer allocation dominance, impose new gates, or implement full GPU
CHIMES solely from reviewer votes. No commit/push performed in this bundle.

# F-P2.8 driver-faithful SNRT harness and mkrun GUI evidence

Date: 2026-09-07 (Asia/Seoul)

## Scope

This bundle adds two bounded items:

1. `mkrun.py --mode gui` (with the existing `--gui` spelling retained as a
   compatibility alias). The AMR base and maximum levels are collected in one
   side-by-side form and submitted as one validated pair; Back removes the
   pair as one answer batch.
2. A native RAMSES harness that reaches the production
   `amr_step -> snrt_ramses_advance_level` call. Baseline and forced receiver
   failure cases exercise the coupled source/RT transaction.

The harness uses `SNRT_DRIVER_TEST_SEED_SOURCE=1`, an explicit nonproduction
driver-only seed, because the current fresh sink-creation path is not a valid
SNRT source fixture. The driver accepts this seed only when the existing
nonproduction diagnostic mode is explicitly enabled; it is never enabled by
a normal namelist or by a production invocation. The test does not claim a
physical sink-formation or AGN SED validation.

## Code and build boundary

- Project: `/gpfs/kjhan/LRD_JWST`
- Remote: `git@github.com:kjhan0606/LagRamses.git`
- Active source order: `VPATH = $(PATCH):../patch/cuda:../patch/oct_tree:../patch/cuRamses:../$(SOLVER):../aton:../hydro:../pm:../poisson:../amr`, with `PATCH=../patch/lagRamses` and `SOLVER=hydro`
- Build: `make -C bin -j1 SNRT=1 USE_CUDA=1 USE_FFTW=0 ramses`
- Binary: `/gpfs/kjhan/LRD_JWST/bin/ramses_final3d`
- Binary SHA256 for the recorded run: `4c2b9597f785901f3956a7c1039644d67f965bc10bf8a64e0d1e110559506d0d`

The binary was rebuilt from the current working tree with the active
`patch/lagRamses` driver and the full SNRT module graph. `gpu_hydro` is false
in this small harness to avoid coupling this source transaction test to the
separate CUDA hydro inter-buffer path; SNRT still requires and detects the
CUDA device.

## Native execution

Runner:

`simulation/snrt/tests/run_snrt_agn_driver_faithful_smoke.sh`

Configuration:

`simulation/snrt/config/snrt_agn_driver_faithful_smoke.nml`

Reference contracts:

- `simulation/snrt/config/snrt_group_contract_reference_control_v1.nml`
- `simulation/snrt/config/snrt_secondary_table_contract_v1.nml`

Recorded run root:

`simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_20260907T080854_2554271`

Baseline log markers:

```text
SNRT_DRIVER_TEST_SEED_SOURCE applied: NONPRODUCTION
SNRT_RT_TRANSACTION_COMMIT_PASS level=3 iteration=1 residual=  0.0000E+00
SNRT_RT_CLOSURE_PASS level=3 leaves=512 photon_nonnegative=1 species_simplex=1 thermal_finite=1 unassigned_code=  0.0000E+00
active sources: 1
```

Injected log markers:

```text
SNRT_DRIVER_TEST_SEED_SOURCE applied: NONPRODUCTION
SNRT RT transaction rollback: class=receiver level=3 iteration=1 residual=  1.7977+308
SNRT_RT_DIAGNOSTIC_FAIL_CLOSED class=receiver level=3
```

The runner returned:

```text
F-P2.8_CASE baseline PASS return_code=0
F-P2.8_CASE injected PASS return_code=0 expected_mpi_abort=1
F-P2.8_PASS run_root=/gpfs/kjhan/LRD_JWST/simulation/snrt/runs/fp2_8_agn_driver_faithful_smoke/job_20260907T080854_2554271
```

No RAMSES output directory was created. The runner rejects fatal/error markers
even when `clean_stop` returns status zero. The injected case may return
through the diagnostic fail-closed path without a nonzero launcher status; the
rollback markers and the absence of a commit marker are authoritative.

## GUI verification

Command:

`python3 -B -m unittest patch/cuRamses/aux/test_ramses_run_gui.py -v`

Result: 17 tests passed; 1 real-widget test skipped because this node had no
usable graphical display. The test suite covers `--mode gui` dispatch,
side-by-side level-pair replay and widget semantics when a display is
available, CLI/GUI byte equality, validation, Back behavior, preview-only
operation, and safe confirmed saves.

## Boundary and follow-up

An earlier attempt with `create_sinks=.true.` reached
`patch/lagRamses/sink_particle.kjhan.f90:kjhan_quenching` and terminated with
SIGSEGV before a sink was created. That is a separate fresh sink-formation
blocker and is intentionally not hidden by this harness. Fresh sink creation,
restart/migration persistence, MPI source ownership, physical AGN SED/MAD,
and live dust are outside this evidence and remain unapproved for production
science.

Follow-up: `sink_formation_prerequisites_2026-09-07.md` corrects the crash
location to the caller's post-quenching rho_star access and records the
unsupported-input rejection plus a successful zero-source formation scan
with Poisson enabled. Actual physical sink formation remains unvalidated.

The Opus5 bundle-end audit and its post-audit closure are recorded in
`provenance/opus5_f_p2_8_bundle_end_audit_2026-09-07.md`.

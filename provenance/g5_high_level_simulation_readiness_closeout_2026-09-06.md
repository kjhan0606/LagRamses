# High-level RT/feedback/dust simulation-readiness closeout

Date: 2026-09-06 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)

## Decision

The repository now has two explicitly admitted executable smoke profiles:

| Profile | Decision | Evidence |
|---|---|---|
| SNRT/CUDA reference-control | `READY_FOR_BOUNDED_QUALIFICATION` | D4 Slurm job `333272`, two MPI ranks/two A10 GPUs, startup-contract and receiver-rollback PASS |
| CPU legacy comparison | `READY_FOR_BOUNDED_COMPARISON` | two-rank `nstepmax=2` smoke, `feedback_mode='legacy'`, finite-state diagnostics and `Run completed` |
| Physical dusty production | `BLOCKED` | The active checkout's 2026-10-21 production preflight found 17 release/admission blockers; the 2026-10-30 source-terms recheck confirms physical-data licenses remain unresolved. The A40 current-source comparison `442563` is still pending. Selected dust live/restart paths have bounded evidence; cosmological qualification remains open. |

The first two profiles are qualification/comparison runs, not a scientific
production release. The physical profile cannot be enabled by selecting a
different binary or namelist.

## Verification completed

The final bounded checks completed with PASS:

```text
G4_DUST_CLOSURE_PASS tests=3 mapping=explicit thermal=one_pass backend=cpu
SNRT_NATIVE_DUST_CONTRACT_RUN_PASS
SNRT_NATIVE_DUST_MAPPING_RECEIVER_RUN_PASS
SNRT_NATIVE_DUST_TRANSACTION_ALL_OK
FP1_PHYSICAL_PACKAGE_ADMISSION_TEST_OK
G2_SOURCE_SELECTION_GATE_TEST_OK
PARKED_AGB_SOURCE_CHECKS_RAN
G2_LIMONGI_PHASE_MASS_HISTORY_TEST_OK
FP1_LC18_FAILED_WIND_CROSSCHECK_TEST_OK
```

The reference-control preflight used the recorded SNRT/CUDA binary
`simulation/snrt/build/g5_startup_contract_admission/ramses3d` with SHA-256
`e30a938397781ee02eb060222937007f659bed4a0e5d268c1cee125901e019f0` and the
effective initialized-RAMSES namelist. The CPU comparison used
`bin/ramses_final3d` with SHA-256
`5eaef4bdc321abbe2494a481d229b7170221c0be623d572f722b3f326253bff2`.

## Physical blockers that remain intentional

- The production manifest still lacks approved, checksummed registrations for
  the physical stellar yield table, stellar SED, AGN SED, and dust model.
- `p4_snapshot_hdf5` remains a not-migrated, hash-pending reference payload;
  the existing small P4-derived inputs are validation/control artifacts, not a
  substitute for that physical input contract.
- The thermal atlas has an approved physics provenance record but its license
  status is still pending explicit confirmation.
- The September baseline's statement that the live RAMSES driver lacked a
  persistent dust field and dust receivers is superseded for the selected
  source-active profile by the 2026-10-20 integrated job 430823: it exercised
  live two-size dust, M5 source/absorption, grain heating, CPU CHIMES, and
  mixed-AMR coupling with photon/energy accounting. That bounded smoke does
  not qualify every dust thermal/momentum/scattering/IR path, cosmological
  evolution, production-scale AMR, or a release binary; those remain separate
  qualification requirements.
- The current checkout is dirty because unrelated migration/high-mass work is
  in progress. A physical release requires a clean named commit and a rebuild
  of the selected binary from that commit.
- The 40--120 M☉ source seam and AGB source choice remain review-only; parked
  source readers do not constitute physical approval.

These are physical-release conditions, not reasons to repeat the already
passing native contract smokes. The next implementation work is therefore
physical asset admission and live receiver activation, followed by the G5
hydro/restart/AMR qualification matrix.

## Reproduction entry points

```text
simulation/snrt/runs/simulation_ready_preflight.sh
simulation/snrt/tests/run_g4_dust_closure.sh
sbatch simulation/snrt/tests/run_snrt_bundle_gate.sbatch
```

The bundle gate and its child test runners resolve scratch paths from the
checkout containing the scripts; the gate refuses direct login-node execution.
The A100 Slurm wrapper is the supported invocation for the current checkout.

The latest record-keeping repair was committed and pushed as `7330087`.

After that record update, the integrated gate was rerun from the current
workspace and returned `SNRT_BUNDLE_GATE_PASS` at commit `7330087`.  It
reported production-link PASS in 178.207 s, AGN partition PASS, five native
symbols, dust-ledger receiver PASS, four thermochemistry negative cases, ten
spectral negative cases, eleven transaction rollback/configuration negative
cases with two-rank MPI coverage, A10 CUDA PASS, production-negative PASS,
and diff/conservation checks PASS.  The temporary linked test binary had
SHA-256
`684dec40632e94058f5fd35028aeda20985418eb37a365bd33303f31f7f6ada0` and
was removed by the gate after the run.

## Current G0 audit scope correction — 2026-10-10

The current checkout is `/gpfs/kjhan/chimes-cvode-gpu-20261001` with origin
`kjhan0606/LagRamses`. The external-assets registry is a project-wide
migration/comparison inventory, not a list of assets all selected into one
production run. The G0 auditor previously required production license and
provenance metadata, and rejected blocked/comparison statuses, for every
registry entry. It now applies those admission checks only to IDs in
`required_production_assets`; an absent selected ID reports
`required_asset_not_registered` without a redundant metadata-missing error.
Registry-wide duplicate-ID and absolute-path integrity checks remain.

Verification: `simulation/snrt/tests/production_manifest.py` passes, including
an unselected `not_migrated` legacy comparison asset (must not block), and
selected legacy, missing metadata, missing registry entry, blocked-status,
dirty-worktree and commit-mismatch cases (must fail closed). The auditor now
checks the declared repository path, live Git HEAD and actual worktree state
against the manifest's existing source-identity fields. `git diff --check` is
clean. A read-only audit of the current manifest returns
`production_gate_pass=false` with 13 concrete blockers: pending thermal-atlas
license; missing metadata for the selected executable, CUDA smoke binary and
P4 input; the manifest's recorded dirty-source state and a source-commit mismatch
(the auditor also enforces a clean live worktree);
four unregistered physical stellar/AGN/dust assets; and the P4 input's missing
path/checksum plus `not_migrated` status. These are not waived by narrowing
the audit, and no external registry, physical asset, license record or source
table was changed.

The 2026-09-06 statement that the live driver had no dust state/receivers is
historical and superseded by the integrated relative-motion and effective
population evidence in `dust_remaining_implementation_2026-09-09.md` and
`effective_population_execution_2026-09-11.md`. Those bounded
noncosmological CPU/OpenMP results do not convert this G0 manifest to PASS or
qualify the cosmological production profile.

## Current checkout binding correction — 2026-10-11

The active checkout for this continuation is
`/gpfs/kjhan/chimes-cvode-gpu-20261001`, origin
`git@github.com:kjhan0606/LagRamses.git`, HEAD `2a88e02`. It is a distinct
worktree from `/gpfs/kjhan/LRD_JWST` (same origin, HEAD `9bec08c`). The
`manifests/production_readiness_manifest_v1.json` present in this checkout
still binds `repository.path`, `project_root`, source registry, and environment
to `/gpfs/kjhan/LRD_JWST` and recorded commit `45cee683...`. Consequently,
running `audit_production_manifest.py` with that manifest causes its live Git
identity probe to inspect the old checkout, not this active source tree. Its
reported 13 blockers are not a verified G0 result for HEAD `2a88e02` and must
not be used as such.

This does not establish G0 PASS for the active checkout: its source worktree
is dirty, and it has no manifest/environment pair bound to its current root and
identity. The old manifest's `p4_snapshot_hdf5` requirement also conflicts
with `simulation/snrt/P4_EXECUTION.md`, which calls P4 an ingestion artifact,
not a production scientific input. Do not edit the external-asset registry or
carry that requirement into a current-checkout manifest without re-evaluating
the actual production profile. Existing physical-source/license concerns stay
unresolved until their provenance is independently bound to the active
checkout. This is a provenance correction, not a new qualification gate.

The active bounded Grackle reproduction is also not a high-level-physics
qualification: job 418474's pinned binary is
`c16b06efca712df30a44fa8c312bab707dc77669db6c29824edb5783e8a58ebe`, while
its effective namelist selects `radiation_transport='none'`,
`dust_mass_enabled=.false.`, `agn=.false.`, and `delayed_cooling=.false.`. Its
yield-table and IC paths reference `/gpfs/kjhan/LRD_JWST`; the run is correctly
labelled exploratory and is testing cosmological Grackle/restart failure
reproduction, not RT/feedback/dust production. The canonical
`bin/ramses_final3d`, `simulation/snrt/build`, and local SNRT `.venv` are absent
in the active checkout; the job-specific executable is not a substitute for
those production identities. Any future production candidate must bind its
effective inputs to explicit external paths, checksums, provenance and
approved status rather than inherit this exploratory profile.

## Current coupled M5 evidence and separate Grackle repair — 2026-10-11

The active checkout now also has successful bounded **co-moving dust + M5 +
actual AGN source + CPU CHIMES** evidence. Slurm job 419549 is `COMPLETED 0:0`
in 41m27s, and its native log contains `Run completed`. Its effective namelist
selects `snrt_mn`, moment order 5, live two-size dust, `chimes_neq_v1`, and
AGN/sink coupling with CPU hydro. The acceptance result explicitly says
`reference_controls=NONPRODUCTION`, `momentum_assert_cells=16896` and
`raw_outputs=0`. The stored momentum increment, photon/energy accounting and
positive actual-source/absorption checks are detailed in
[the existing M5 execution record](snrt_m5_radiation_momentum_2026-10-11.md).
This supersedes any inference that no coupled live M5 receiver has run in
this checkout, but does not turn a reference-control fixture into an approved
physical production source package.

The scopes remain distinct: this M5 fixture is noncosmological, uses a short
three-step trajectory and does not qualify the rejected relative-dust/AGN
layout. Grackle job 419596 instead evolved the saved cosmological RT/dust-off
profile from coarse step 380 and reproduced the composition closure
rejection in interval 391, terminating FAILED at 21:29 KST. The native trace
captures an admitted pre-hydro state becoming inadmissible across the
hydro-update window, before chemistry. That run does not test the M5
coupling, and the successful M5 run cannot waive its chemistry admission
failure. The corrected Grackle runtime candidate has passed actual-routine
AMR, face transport and thermal regressions. Its initial bounded corrected
replay 419694 completed intervals 380 and 381 without rejection, but the
measured cost predicted exceeding its unextendable three-hour limit. The
driver deliberately cancelled it, preserving logs and checkpoint, and
submitted the already verified used-face optimization as 419714. The latter
is RUNNING on H100 at 22:08 KST with the identical effective inputs and
stopping criterion, and four-hour headroom. The optimization passes all
three-axis native face regressions; passage of the actual interval-391
failure remains pending and is not inferred from those native passes. See
[the existing Grackle execution record](galaxy128_grackle_calibration_baseline_2026-10-08.md).

No new gate or completion bundle is introduced by this reconciliation.
The active physical-source identity/provenance concerns remain unresolved
where previously documented, and the old-root production manifest is still
not a release identity for this checkout. Scientific calibration and
resolution convergence remain the separate operator-approved science plan.

### Selected input roles checked against the actual files

Read-only inspection confirms that 419549's BPASS v3 wrapper explicitly uses
`population_binding='independent_radiation_reference'`, `status='reference_control'`
and an empty approval ID. Its group and dust contracts likewise explicitly
retain reference-control status. These are actual input roles, not inferred
from the word "smoke" in a filename.

A feedback-matched alternative **does exist** in the already used external
package: `/gpfs/kjhan/LRD_JWST/.effective-population.gQodLz/inputs/source.nml`,
selected by that package's recorded environment. It is v5,
`population_binding='common_imf_mixed_evolution'`, IMF1/population2, and
binds the existing mixed-lowmass node file. It also explicitly says
`status='reference_control'`, empty approval ID, and no late radiation after
the track/terminal. Thus the source connection is not universally absent;
neither this declared truncated effective-population source nor the
independent BPASS control is silently relabelled a full approved physical
SSP. The previously approved effective-model scope and bounded matched
source evidence remain valid. No same-author library, resolved microscopic
physics, old full-hull admission policy, or new test program is imposed by
this input-role check. No external source, wrapper, registry or approval
record was changed.

### Grackle repair result reconciled — 23:43 KST

The previously pending replay 419714 has now terminated `COMPLETED 0:0`,
including its real compute step, after 1h36m38s. Its native log reaches 392
and `Run completed` without the interval-391 Grackle rejection, MG failure,
or nonzero NaN counters; the pinned executable/effective-input hashes match.
The conservative AMR/face-composition repair is therefore closed for that
selected cosmological RT/dust-off model. It is not a coupled M5 or a global
physical-production qualification, and does not waive the input roles or
active-checkout release-identity issues above. No new raw dump was generated;
the saved step-380 restart remains necessary and is preserved.

### Release-admission entry point completed — 2026-10-12

`simulation_ready_preflight.sh --mode physical-production` no longer ends
in an unconditional historical refusal. It calls the existing G0 auditor
with the actual launch checkout, executable and namelist. Admission requires
the existing clean source/selected-asset license/provenance rules, a locked
environment recipe, `production_status=ready`, and an approved `release_record`
binding the selected executable, actual file hashes and physics-selection
environment. A manifest naming another checkout is rejected **before** any
Git or asset audit there; its historical asset failures are not propagated
as current-checkout requirements. Native RAMSES environment records explicitly
use `backend.runtime=native_ramses` with compiler/MPI identities and the
existing hashed lock mechanism; P0/JAX records retain their pinned-version
check. No new gate, monitoring framework or physical model was introduced.

The existing production-manifest tests PASS, including native admission,
changed binary/namelist/environment, absent approval, unselected executable,
missing native environment fields and disabled clean/metadata policies. A
tiny clean synthetic checkout also exercises the actual shell entry point:
the synthetic approved release admits, and a pending release rejects.
These are launch-identity tests, **not** scientific validation or approval
of any real source. The current historical manifest still rejects with
`production_launch_root_mismatch`; no physical asset, external registry,
license or approval record was edited. README documents the interface.

Re-reading `effective_population_execution_2026-09-11.md` confirms that the
Kroupa mixed-v5/DTD-N100 profile and matched PARSEC effective-population source
were already selected and exercised within their declared approximation
domain. The driver's repeated model-selection question is withdrawn. This
does not relabel the source wrapper's reference-control role as an approved
full SSP or approve a physical release automatically.

### Selected effective-population M5 connection — 2026-10-12

The previously selected Kroupa/effective SSP2 matched mixed-v5 source, empirical
Ia population, existing-BH AGN, co-moving two-size dust and CPU CHIMES have now
completed a bounded **M5** coupled trajectory in this active checkout. Job419903
and its actual simulation step both terminate `COMPLETED 0:0`, coarse step4,
with the profile-specific acceptance marker. This is distinct from the earlier
AGN/BPASS reference fixture and from the RT/dust-off Grackle science run.

Existing source/receiver accounting, finite-state and stored-momentum checks
pass, and raw dump count is zero. Actual stellar births reach256. Total
injected photons4.471853742514118e67 exceed an intentionally excessive
AGN-only rest-mass/lowest-edge upper bound by571604.51, establishing nonzero
stellar-to-M5 injection without attributing the aggregate ledger to stars
by assumption. Actual AGN luminosity is also positive. No due Ia event is
claimed in this ~0.206Myr trajectory; the prior nonzero event/debit test remains
its evidence. See the existing
[M5 execution record](snrt_m5_radiation_momentum_2026-10-11.md) for input/binary
identities, source attribution, ledger residuals and declared limits.

A launcher affinity defect was found and corrected: two MPI ranks shared the
same OpenMP core places under unbound srun. The reusable runner now starts
with explicit per-local-rank CPU slices and OMP places (generic Slurm core
binding alone was later disproved by419934). The active job completed after its second
rank's main/worker affinities were separated inside the existing allocation,
without restarting or modifying physical inputs. Its mixed-affinity12m56s
runtime is not a controlled scaling measurement. The failed auxiliary probe
for an absent system gdb is not a failed simulation step.

This closes that bounded selected-source connection, not the whole physical
release. The inputs retain explicit reference-control/truncated-source roles;
the active-checkout release identity is not approved, and this no-dump run
does not independently prove M5 restart behavior. No parked microscopic
research, old full-library admission policy or new qualification ladder is
made mandatory by these limits.

### M5 restart ownership repair — 2026-10-12, initial strict FAIL

Actual same-MPI coupled HDF5 restart419937 exposed a missing SNRT primary/IR
ownership-transfer hook in forced restart load balancing. The existing
`snrt_regrid` typed payload is now copied through the old expanded peer maps
before cpu ownership changes.419946 verifies that the previously missing
photon inventory is restored; its continuous and resumed simulations both
complete0:0. The corrected launcher's distinct actual MPI/OMP affinities are
also observed. After exact cell-coordinate alignment, all gravity datasets
are bitwise equal. Checkpoint star identities are preserved; new stars require
the declared exact birth-event/label bijection, not file-index pairing.

The full strict restart-state comparison nevertheless remains **FAIL**:
fine HII differs by2.154720087250983e-13 absolute (maximum relative6.1542e-10),
and59 hydro datasets exceed rtol1e-12/atol1e-25. CPU-only control419969 and
source-order control419975 do not remove the residual; the ineffective source
sorting patch is reverted. No criterion is loosened and no broad readiness
is inferred from conservation receipts alone. See the existing
[M5 execution record](snrt_m5_radiation_momentum_2026-10-11.md) for exact jobs,
native repair, comparison semantics, retained active checkpoint and evaluated
raw cleanup. This bounded repair advances the real restart connection; it
does not resolve physical-input roles or active-checkout release approval.

### M5 strict restart residual resolved in traced execution — 2026-10-12

First-interval position-matched native traces420002 show exact restored gas,
chemical, dust/IR and primary radiation inputs and exact source injection.
The first difference occurs in transport, from combining face increments in
MPI-owner-dependent groups. The retained fix integrates each face in fixed
quadrature order, reverses six separate physical face channels, and combines
them in a fixed cell-face order. Shared-face equal/opposite flux accounting
and coarse-receiver SUM are retained; reverse traffic and active-face scratch
increase explicitly, without adding persistent angular storage.

Default-auto same-profile420008 now completes0:0 and passes the unchanged
rtol1e-12/atol1e-25 comparator:449 datasets,427 bitwise,39 exact contract
attributes and exact clocks. After transport all traced radiation entries
match exactly; after material all ion fractions match exactly. This replaces
the earlier strict FAIL for that bounded static64-cell MPI2 case, not a
universal coarse/fine or rank-count reproducibility claim. The temporary
trace was removed; current trace-free build420011 is verified, and final
coupled rerun420012 is under evaluation. Physical-input/reference roles and
active-checkout release approval remain separate, explicit unresolved items.

Trace-free final420012 subsequently completed0:0 in3m05s with the same449-
dataset PASS,427 bitwise fields and39 exact contract attributes. Actual
stellar return is identical in both branches; native photon/energy accounting
and finite-state guards pass, and default-auto CPU/GPU closure counters are
both nonzero. The declared same-MPI/static64-cell numerical issue is closed.
The temporary recorder is absent from the retained source. Evaluated raw
cleanup420015 is submitted only after this driver evaluation, excluding the
ongoing calibration/science archives. See the existing M5 execution record
for exact binary identity, bounded claim, resource costs and cleanup manifest.
Cleanup420015 has now completed0:0:42 exact evaluated files/1.63GiB removed,
with checksums and removed=true retained in the current build's deletion
manifest. Compact PASS evidence and all inputs/logs/builds remain; science
outputs and the active Grackle calibration checkpoint remain untouched.

### Current mixed-AMR evidence and release-candidate binding — 2026-10-12

Current trace-free binary420011 also completes source-active mixed-level
job420024 (MPI2/OMP2, A100)0:0 in2m57s. The inspected mesh has56 level2
coarse leaves and64 level3 fine leaves, with successful M5/material commits
on both levels, positive source/absorption/dust heat and417 final stars.
Native finite-state, momentum and nonconvergence guards pass; no raw dump
was produced. This closes the bounded mixed-level regression for the face-
accumulation correction, not cosmological qualification, AMR restart
equivalence or an analytic inter-level flux/convergence study.

`manifests/native_m5_release_candidate_2026-10-12.json` records the actual
checkout, binary/options, effective inputs and existing source-contract roles,
plus exact hashes of restart and mixed-AMR evidence. It deliberately retains
pending source/release approvals and is not a production launch manifest.
The old-root G0 manifest remains historical; its unselected P4 snapshot,
thermal atlas and JAX requirements do not transfer automatically to this
native profile. See the existing M5 execution record for identities and
limits. No new gate, physical approval, commit or push is inferred.

### Current-checkout M5+dust+CHIMES integrated A40 smoke — 2026-10-20

Job430823 completed0:0 on one A40 (MPI1/OMP8; 6m58s including build). The
source-active effective-population profile ran M5 order5 across mixed AMR
levels, live two-size dust and CPU CHIMES, with56 coarse and64 fine leaf cells.
All eight level commits completed; the photon-ledger residual/source ratio is
at most1.12e-15 and the largest printed M5 energy-ledger residual is
3.55e-22. The RAMSES trajectory took128.15s; no raw dump was written.

This is an integrated interface smoke, not a performance comparison, restart
qualification or physical production release: source/SED/material contracts
are reference controls. The build used checkout f459f6d plus then-dirty source;
binary and CHIMES-library hashes are retained in
`.m5-current-source-active-dust-chimes-build-430823/build.sha256`, but a full
source snapshot/clean release commit was not produced. The result therefore
does not update the native production-candidate binding above or close its
release/license conditions.

### Native CHIMES thermochemistry closure smoke — 2026-10-21

The production Makefile's `snrt_thermochemistry_smoke` target was linked with
`CHIMES=1`, `CHIMES_CUDA=0`, and the A40 integrated build's current CHIMES
objects/libraries, then run with the selected CPU backend and the same pinned
reaction, spectral and molecular tables. `SNRT_NATIVE_THERMOCHEMISTRY_OK`
passed, including the actual 157-species element/charge closures and the new
negative-trace test: a `-256*epsilon*Fe` roundoff residue is reconciled with
elemental closure, while `-1024*epsilon*Fe` rejects without publishing a
modified state. The default standalone H/He-only script was also run; its
compile mode omits `SNRT_CHIMES`, so it does not execute those conditional
CHIMES assertions. The explicit production-target run closes that coverage
gap for this bounded code path. This is a native CPU smoke, not a live
simulation or production release qualification; the linked objects came from
the same dirty source tree recorded for job430823, not a clean release commit.

### Current-checkout physical-production admission audit — 2026-10-21

The read-only `simulation_ready_preflight.sh --mode physical-production` was
run against the candidate executable and continuous-run namelist with the
manifest-declared environment bindings. Launch identity was checked, but the
gate correctly rejected production admission. After removing caller-vs-
manifest environment differences, 17 substantive blockers remain: physical
release/admission is not approved; there is no clean release commit, the
repository is dirty, and the candidate source snapshot differs at
`bin/Makefile`; explicit data terms remain unresolved for CHIMES data, Draine
D03, the KL16 mirror, and PARSEC VMS; release conditions 1--3 remain open;
and per-release license records are absent for the ordinary yields, stellar
SED, AGN groups, dust optics, and SNIa contract. The three manifest conditions
are to (1) record selected external-data license/provenance for this exact
release, (2) build and qualify a binary from the exact clean release commit,
and (3) pass launch admission with an actual-root native environment/selected-
asset manifest. Public availability is not treated as a license. This audit
does not approve or waive any physical input and creates no new gate.

### Current-source CPU M5+dust+CHIMES mixed-AMR smoke — 2026-10-30

The M5 source-active runner now checks its effective `noutput`, `aout`,
`tout`, `foutput`, and `fbackup` policy for both run profiles before launch.
The change is commit `24f62e7`.

CPU build job `1243592` completed in 4:03 on grammar `normal` with
`CHIMES_CUDA=0`, `USE_CUDA=0`, `SNRT=1`, `DUST_LIVE=1`, `CHIMES=1`,
`HDF5=1`, and `NVAR=187`. Its clean source worktree is bound to commit
`24f62e7f362e79b1c6b5ffc8fc304715e3c71a9b`; the binary SHA-256 is
`561e4ee63e2dc422990e2453260b05f047c3fe2a1107f2b8035abd2ae01395b2`.
The build's peak RSS was `26902024K`.

Run job `1243595` completed on grammar `normal` (MPI 1, OpenMP 4) in 154.44 s
of RAMSES runtime, with peak step RSS `786996K`. The effective-population
profile exercised M5 order 5, source-active mixed AMR, live two-size dust and
CPU CHIMES for four coarse steps. Levels 2 and 3 each committed four times;
the run had 56 coarse and 64 fine leaf cells. Across eight photon-ledger rows,
the six positive-source rows have maximum residual/source ratio `3.32e-15`;
zero-source rows have zero residual. No recoverable CVODE warning, CHIMES
cell rejection, nonfinite state, or nonconvergence marker was found. The runner reported
`M5_COUPLED_AMR_INTERFACE_PASS`.

The checked effective namelist has `noutput=1`, future `aout=2.0d0`,
`tout=1.0d30`, and `foutput=fbackup=1000000`; no raw dump was produced. Its
SHA-256 is `dbd00a880b477d0bdd4b905c68dae69582d81774cae317d0f100ecd5c1999dc9`.
The run's input inventory SHA-256 is
`376a7c404503cd62eace0bec674d8deea6956571b0ead174698fbf6679b55095`.

This closes a current-source CPU integration smoke and its output-schedule
preflight coverage. It remains a bounded nonproduction run: AGN/dust/source
contracts are reference controls, no CUDA execution was performed, and the
native release-candidate manifest, external-data terms, and physical-release
conditions remain open.

### Current-source CPU M5 effective-population restart equivalence — 2026-10-30

Job `1243610` ran CPU-only on grammar `normal` (MPI 1, OpenMP 4) using the
binary built by job `1243592` from clean source commit
`24f62e7f362e79b1c6b5ffc8fc304715e3c71a9b` (binary SHA-256
`561e4ee63e2dc422990e2453260b05f047c3fe2a1107f2b8035abd2ae01395b2`). Both
the continuous four-coarse-step trajectory and the trajectory restarted from
coarse step 2 completed. The first Slurm wrapper exited nonzero only in its
post-run comparator: the comparator incorrectly required `ncpu=2` for this
MPI-1 run. The comparator now derives rank count from the source checkpoint,
requires both final outputs to match it, and reports that rank count in scope.

The corrected comparison returned `M5_COUPLED_RESTART_PASS`: 449 datasets
compared, 439 bitwise-identical, and all 10 floating differences within the
predeclared `rtol=1e-12`, `atol=1e-25` tolerance (maximum tolerance ratio
`3.52e-4`; zero violations). All 39 contract attributes and six clock fields
matched exactly; physical grid ordering and particle identity also matched.
The resumed log contains both `SNRT_HDF5_RADIATION_RESTORE_PASS` and
`STELLAR_SOURCE_RESTART_IDENTITY_PASS`. The exact report is retained at
`.m5-effective-restart-1243610/restart-comparison.json`; launch, effective
namelists, execution logs, input inventory and raw HDF5 checksums remain beside
it. After evaluation, the three raw HDF5 dumps (303,845,800 bytes total) and
their output directories were removed; the exact targets are recorded in the
run's `cleanup-manifest.json`.

This qualifies same-rank CPU restart equivalence for this bounded
effective-population/M5/live-dust/CPU-CHIMES reference profile. It is not a
GPU comparison or a production release; source and dust/AGN contracts remain
reference controls and the release-candidate blockers above remain open.

### Exact-current-commit CPU M5+dust+CHIMES AMR smoke — 2026-10-30

Build job `1243662` completed on grammar `normal` with 4 CPUs, 32 GiB and
`M5_BUILD_CUDA=0` in 3:38 (peak RSS `26293816K`). Its isolated source tree was
clean at commit `3628b23e9aa5e905fa7ef8231c2f6e98298c342a`; the binary
SHA-256 is `ee0280258f2985710f88b2b4e78aeb2e31a5e5e75b51d52bcb0785b627ee5362`.

Simulation job `1243669` completed on grammar `normal` with MPI 1/OpenMP 4,
CPU-only, in 2:49 Slurm elapsed (RAMSES runtime 154.80 s; peak step RSS
`780592K`). It used that exact-commit binary for four coarse steps with
`AMR_INTERFACE_CHECK=1`, M5 order 5, live two-size dust, CPU CHIMES and the
selected effective stellar population. The runner returned
`M5_COUPLED_AMR_INTERFACE_PASS`: 56 coarse and 64 fine leaf cells, eight M5
level commits, eight photon-ledger rows, and maximum positive-source
residual/source `3.32e-15`. All 75 emitted NaN checks were zero. The effective
namelist used `nstepmax=4`, one future `aout=2.0d0`, `tout=1.0d30`, and
`foutput=fbackup=1000000`; no `output_NNNNN` dump was created.

The run and build identities, effective input, log, compact result, hashes and
no-output cleanup record are retained in `.m5-effective-amr-1243669/` and
`.m5-ramses-momentum-build-1243662/`. This confirms the bounded CPU
source-active AMR connection from the current commit. The available grammar
Slurm partitions exposed no GPU GRES, so this does not qualify the CUDA/hybrid
path. The selected AGN/dust/data contracts remain reference controls, and the
upstream data/license and production-release blockers remain explicit.

### Current CUDA qualification and source-term recheck — 2026-10-30

`sinfo --all` on the active grammar Slurm controller showed only `debug`,
`normal`, `normal1` and `normal2`, all without GPU GRES. A current-source GPU
job could not be submitted from this controller. The previous A40 integrated
run `430150` was built from `5ea1cb6`; comparison with current commit `3628b23`
shows intervening changes in `bin/Makefile`, the CHIMES bridge and CUDA RHS,
the CUDA RHS header, and `snrt_moment_ramses.f90`. Consequently, that older A40
run does not qualify today's hybrid/CUDA source. At that 2026-10-30
grammar-controller check, the latest exact-commit result was CPU-only and no
typed GPU partition was available there, so GPU qualification remained open.
The later syntax A40 submission is recorded below.

The official-source data-terms recheck is recorded in
`manifests/native_m5_release_candidate_2026-10-12.json` under
`primary_source_recheck_2026_10_30`. It found no new explicit reuse license
for the selected CHIMES, Draine D03, PARSEC VMS or KL16 data, so those terms
remain blocked. The production candidate was not promoted and no data-use
approval was inferred.

The preceding GPU-availability statement is a snapshot of the grammar
controller and predates the current syntax submission. Typed-GPU execution is
now queued as A40 job `442563`; its effective-population/HDF5 comparison and
current Slurm state are recorded in
[`m5_agn_cpu_reference_integration_2026-10-30.md`](m5_agn_cpu_reference_integration_2026-10-30.md).
The job remains pending resources, so current-source hybrid GPU qualification
is still open until the run and CPU/GPU state comparison complete. This does
not change the physical-data licensing blockers above.

### Repeated-event delayed-cooling timer regression — 2026-10-30

The focused Fortran unit now advances the active timer by 3 Myr before a
second SN event, changes the gas density through returned mass, and asserts
that the remaining delay stays at 17 Myr rather than restarting at 20 Myr.
It also checks subsequent non-SN mass return, time decrement/expiry, and
re-arming after expiry. `F90=ifx bash tests/run_stellar_feedback_policy_unit.sh`
passed on grammar with the Intel 2025.3 compiler (`delayed cooling runtime:
PASS`; the companion stellar-feedback policy unit also passed).

This directly tests the timer helper's repeated-event semantics. The separate
job1243549 restart smoke passed same-model HDF5 restore and rejected a changed
`t_diss`; that smoke does not exercise repeated feedback events. Neither check
is a source-active cosmological feedback qualification.

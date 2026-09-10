# Gas MHD completion: four approved bundles

Operator: execute 1 -> 2 -> 3 -> 4; bundle 4 is REQUIRED, not optional.
Final objective: physically justified, usable MHD in the existing lagRamses
RT/stellar+AGN feedback/dust simulation, with useful CPU/GPU execution.
Dust Lorentz/charging and field-aligned CR transport remain long-term.
No commits, pushes, long production jobs or destructive output operations
are implied by this approval.

1. Close gas numerics/cosmology: native Alfven/strong-shock and zero-field
   checks, homogeneous expansion, boundary semantics and different-MPI
   restart. Repair actual defects; reuse existing AMR/HDF5 evidence.
2. Couple the existing physical operators: gravity/cooling/RT/SF/feedback,
   sink/AGN energy and flux handling, coadvected dust/chemical carriers and
   nonthermal energies. Fix a supported profile, not every optional model.
   No new magnetized star-formation or jet field-injection model.
3. Freeze and exercise that CPU/MPI profile in a small integrated run with
   dynamic AMR, load balancing and restart; update setup tools and support
   boundaries. This is application-qualified, not universal production.
4. Thread-safe MHD batches and real CUDA numerical work using available
   stream slots, CPU execution when busy. Preserve CT/reflux correctness,
   check CPU/device parity and measure performance before enabling defaults.
   Do not call a dispatch stub or unrelated GPU workload MHD acceleration.

One coherent implementation sequence, not a new gate per helper/test.
Read-only planning review asks Q-GOAL FIRST, Q-LEAN SECOND, then concrete
physical/coding feasibility concerns. Existing approved scope is not reopened.
End verification is performed by the driver per the operator's later
instruction. Each new RAMSES calculation gets a fresh directory and an
explicit effective-input, physics and output-budget launch report.

Status: all four approved scoped bundles implemented and their native MHD
checks complete. Shared nonmagnetic build regression also passes. No blanket
production-physics approval is inferred.

## Planning review disposition

Fable read-only review: `fable_gas_mhd_completion_plan_2026-09-10.txt`.
Q-GOAL affirmative; Q-LEAN cuts adopted: production profile uses periodic
boundaries, so no physical-boundary test campaign. Changing MPI size is
tested once with bundle 3, not duplicated in 1. Setup changes only follow
actual new controls. GPU parity reuses numerical/integrated inputs.

Bundle 2 profile: channel-resolved stellar feedback, atomic cooling/RT,
coadvected bulk dust with existing material chemistry (not separate dust
momentum), existing accretion-powered AGN where its energy/flux contract is
verified; no SGS, no anisotropic CR, no non-ideal diffusion (eta_mag=0).
NENER=0 for the initial profile. Legacy feedback is not automatically admitted
under MHD merely because the gas-only reader needs its inactive namelist.
Sink accretion must retain face flux on the gas mesh; any flux removal or jet
field injection would require a new physical model, outside this bundle.
CUDA remains required, real numerical work, and off by default until measured.

## Bundle 1 results

Native test root: `.mhd-completion.MDLk0E/`. Periodic Alfven wave at t~0.09047:
L1 transverse-field error 1.1367701e-3 (16^3) -> 3.1695457e-4 (32^3),
ratio 3.59 (~1.84 order). This is a two-resolution smooth-wave check, not a
publication convergence study. Periodic Brio-Wu at 32^3 to t=0.02711524:
finite state, minimum rho 0.12465125, minimum thermal energy 0.09945547,
mass 0.5625, total energy 1.33125, discrete divB=0. It verifies shock
stability/conservation, not agreement with a high-resolution reference curve.
Zero-field MHD: rho=1, thermal energy >=0.8999999999999998, total E=0.905.

Homogeneous GRAFIC cosmology (`expansion-grafic`): a=0.1->0.10090210658,
code B follows sqrt(a/a_initial) to 5.55e-17 absolute, divB=0, positive heat.
Gas density equals omega_b/omega_m. The first `expansion` attempt had no IC
directory and was rejected before evolution; it is not a pass.
Real repairs: extended GRAFIC header record length now sets the gas plane
offset, and first MHD passive reads ic_pvar_00001 rather than ic_pvar_00004.
Physical-boundary neighbor lookup was also corrected to Morton and boundary
passives initialized, but that branch is outside the frozen periodic profile.

Benchmark definitions follow the primary Athena test descriptions:
https://www.astro.princeton.edu/~jstone/Athena/tests/cp-alfven-wave/cp-alfven.html
https://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
The wave here propagates along x, with the opposite velocity sign convention
to the cited left-going example; the analytic comparison uses x-t.

## Bundle 2 implementation notes

Independent review of Fable's proposed legacy-SN defect: the two cited
assignments retain E-KE and restore the new KE without a temperature operation.
They preserve B energy, so no speculative subtraction was made there.
The actual AGN receiver cap DID count B as heat. The generic receiver now
accepts an explicit unchanged energy reservoir; native MHD callers supply
cell magnetic energy, and primitive-gas validity and thermal capping exclude
it. Existing AGN tests plus added high-B cap/negative-heat cases pass.
Sink magnetic face flux remains on the gas mesh, unchanged by accretion.

Existing bulk-dust mass evolution itself forbids sinks/cosmology independently
of MHD. Do not silently remove this pre-existing closure restriction: exercise
two explicit existing profiles, (A) stellar+AGN+RT with the existing radiative
dust reservoir, and (B) isolated stellar+coadvected bulk dust evolution+RT.
This does not claim a new cosmological dust mass/accretion closure.

Bundle 2 live evidence: `stellar-dust-cpu` and `stellar-agn-cpu`, four steps,
MPI=2/OMP=2, NVAR=33/NENER=0, binary SHA256
93b3839969704d0e1653fecaf0e96a05a24d8d19998e9604dcff4f8c194cac3a.
Both completed, all 36 stored hydro/face fields finite. Minimum gas thermal
energy 9.67334315e-8 and 2.48949731e-7; max discrete divB 9.76e-19 and
1.41e-18 respectively. Actual particles: 2048 and 3788; AGN deferred energy
9.33028243e61 erg and sink mass 1.02814104e-5 code units. This is nonzero
source/ledger activity, NOT proof that the deferred reservoir was deposited.
The reference-control radiation and independent population caveats remain.
Earlier unsuffixed attempts used invalid backend `cpu` and were rejected;
successful runs use `openmp` and are not confused with those attempts.

Bundle 3 scope adjustment from an actual reader rejection: existing IR v4
forbids load balancing (`nremap/=0`) and requires RT on all AMR levels.
Keep that existing limit: qualify coupled AMR at fixed rank ownership, and
exercise MHD load balancing separately without persistent IR. No claim of
IR+dynamic-load-balance support is made by MHD qualification.

Bundle 3 completed: coupled AMR+restart final fields agree exactly (levels
1..4); changed-MPI gas restart maximum field difference 2.22e-16. Full small
summary and operator-requested raw-output cleanup inventory are in
`gas_mhd_raw_cleanup_2026-09-10.md`. Inputs, logs and binaries are retained.

## Bundle 4 results and driver evaluation

Real HLLD CUDA face batches share the existing pool; busy slots give work
immediately to the OpenMP CPU worker. No GPU-autotuning heuristic was added.
CPU stencil scratch is private; shared cell/face CT/reflux scatters are
protected together. Upstream LLF-fallback faces remain CPU computations.
Reconstruction and 2-D edge EMFs remain CPU: partial numerical offload, not
a full MHD-step GPU implementation. Flags are in both setup entry points.

CUDA native face test: 4096 faces, gamma=5/3 and 2, zero-field and varied
states, passives, plus deliberately occupied stream. Max normalized flux
errors 2.70264070e-15 and 1.70565514e-15; busy fallback passes. Test output:
`.mhd-cuda-build.wUUFN1/hlld-smoke.log` (`HLLD_DEVICE_SMOKE_OK`).

Frozen binaries (both NVAR=33, NENER=0, SNRT=1, DUST_LIVE=1):
- CPU `ramses_mhd_omp3d`: b45a104719523714ef5347bf3d34e356b90aacf6f893ce148574040e83b6799f.
- CUDA `ramses_mhd_hybrid3d`: decdb02a22337782930daa8e3382f027166316663843ee263af9d2b482c66787.

All live dispatch tests MPI=2, OMP=4, one CUDA stream/rank on RTX5000 Ada;
full actual inputs/logs under `.mhd-completion.MDLk0E/`:

| Check | Maximum field difference from previous serial CPU | Minimum heat | Maximum divB |
|---|---:|---:|---:|
| gas-omp4, mixed AMR + remap | 6.66e-16 | positive | roundoff |
| integrated-omp4, stars/RT/bulk dust AMR | 5.42e-20 | positive | roundoff |
| gas-hybrid, mixed AMR + remap | 6.66e-16 | 0.85204087 | 4.89e-15 |
| integrated-hybrid, stars/RT/bulk dust AMR | 2.17e-19 | 2.60915789e-8 | 2.39e-18 |
| alfven-hybrid, 32^3 / 32 steps | 1.78e-15 | 0.14999316 | 1.51e-14 |

Fields were coordinate-sorted before comparison. Gas reference NVAR=21
versus new NVAR=33: compare identical physical fields and map right faces
explicitly (22:24 -> 34:36). Coupled comparison includes all 36 fields.
All fields finite, all runs completed. Actual hybrid face counts are nonzero
on BOTH backends: gas rank1 23424 GPU / 48000 CPU / 125 busy batches;
coupled rank1 13392 GPU / 4464 CPU / 12 busy batches; Alfven rank1
1489536 GPU / 869760 CPU / 2265 busy batches. These are not dispatch stubs.

Performance decision: keep CUDA optional/off by default. The isolated
4096-face test measured CPU 7.572 ms versus CUDA 6.856 ms for 20 batches,
but gas-AMR Godunov timers were OMP4 0.101 s versus hybrid 0.110 s, with
rank variation and too little work for a universal performance conclusion.
Total elapsed times differ in I/O and concurrent-node load and are NOT a
GPU speedup measurement. No production acceleration claim follows from
the smaller total wall times. OpenMP is available; no autotuned default.

Driver end review: accepted within the declared NENER=0 periodic ideal-gas
and existing reference-profile scope. Conserved face layout, energy-source
separation, thread ownership, real device computation, busy fallback and
restart/AMR parity have direct native evidence. Remaining limitations are
the explicitly documented physical/profile boundaries, not new mandatory
micro-gates. No additional audit was launched. Setup tests: 47 discovered,
46 pass, one display skip. See `simulation/snrt/GAS_MHD.md` for operation.

Final housekeeping: raw HDF5 cleanup removed 23 files / 1,504,265,912 bytes
(1.401 GiB); retained test folder 46 MiB (logs, inputs, frozen executables).
After the live comparisons, dispatch logging only was reduced to base-level
`ncontrol` cadence, avoiding per-fine-substep logs. This does not alter the
numerical kernel/dispatch/scatter; final CPU/CUDA rebuilds check the change.

Final rebuilds exit 0: CUDA SHA256
190f10f55cb409011209e5f5d17e31daff5dc2136fcf15ec4f1bb11e8db7e0c7;
CPU MHD c6adfb8b518613947ff5e737c18680b206030545620e1ea003b03111953d2eda;
nonmagnetic hydro 6551b93abe8bbc82684a2268936170228e134b9b33df2b68d037325cd86a4aea.
The frozen earlier binaries above identify the live numerical comparisons;
the final MHD difference is report cadence only. AGN energy-reservoir smoke
rerun ends `AGN_NATIVE_CELL_COUPLING_SMOKE_OK` (log in `.gas-mhd.UDobi3/`).

# M5 stellar + AGN source attribution review — 2026-10-31

## Scope

Review whether the bounded effective-population M5 + live dust + CHIMES run
also exercised the existing accreting-AGN radiation path, and whether the
runner's result record proves that claim. This is a reference-control,
nonproduction integration; it does not approve the stellar or AGN spectra for
science use.

## Evidence from job 443809

- `.m5-effective-dust-chimes-443809/effective.nml` enables the effective
  stellar population SED, `sink`, `smbh`, `agn`, `sink_AGN`, and Bondi
  accretion. Its header describes the run as PARSEC radiation/feedback plus
  accreting AGN, live dust, and kind-7 CHIMES.
- `launch.txt` records the stellar SED enabled, hybrid CHIMES, and reference
  controls. `run.log` confirms the common-IMF stellar population binding,
  the nine-group reference-control contract, M5 order 5, and live dust.
- `agn_coarse_state_v1.jsonl` has records for coarse steps 0–16. Steps 2–16
  show positive accreted BH mass and positive bolometric luminosity; the
  feedback mode is `JET`. The luminosity declines from approximately
  `2.57e34` to `2.26e34 erg/s` over those records.
- The same run reaches the main step and emits nonzero aggregate M5 photon and
  energy source ledgers. Its final HDF5 dump was retained through the matched
  comparisons and removed afterward; the exact file inventory and hashes are
  in [the cleanup manifest](m5_a100_stream_pair_raw_cleanup_2026-10-31.json).

## Code wiring

The source path is connected in the current code:

1. `patch/lagRamses/sink_particle.kjhan.f90` partitions radiated Bondi
   accretion energy when `SNRT_AGN_MODEL=partition_reference_v1` is active and
   accumulates it in `agn_pending_erg`.
2. `patch/lagRamses/snrt_moment_ramses.f90` converts pending AGN energy into
   the configured photon groups and adds those photons and their energies to
   the M5 source arrays at the sink's owning leaf. It consumes the pending
   receipt only after the coupled transaction succeeds.
3. Stellar photons enter the same M5 source transaction through
   `stellar_photon_interval` when the effective-population SED is enabled.

This establishes a code-level route for both live stellar and AGN sources in
the tested profile. The BH sidecar confirms accretion and luminosity during
that run.

## Reporting gap and limits

Before the runner correction in this closeout, the effective-population result
string hardcoded `AGN_state=not_requested`, while the AGN-only path checked
only for a luminous accretion state. The result record therefore understated
the effective run and did not prove AGN photon deposition.

Commit `f260266` adds `SNRT_MN_AGN_SOURCE`, reporting the accepted AGN packet
count, photon count, and represented photon energy. These values are folded
into unused entries of the existing MPI-reduced M5 conservation receipt, so
the change adds no collective. The output appears only after the enclosing
M5/chemistry/dust transaction commits.

CPU build job `1249105` completed on grammar `normal` with
`M5_MOMENTUM_CLOSURE_BUILD_PASS` and `CHIMES_CUDA=0`; the executable and
recorded objects passed SHA256 verification and the executable has no CUDA
runtime dependency. CPU integration job `1249133` completed four coarse steps
on grammar `normal` with exit code 0. It formed 256 stars, loaded the effective
stellar SED, and reached successful M5/dust/CHIMES commits without raw dumps.
The BH sidecar records positive accretion and luminosity in `JET` mode. The
new receipt reports zero AGN photon packets throughout this run. That is
consistent with `agn_partition_release`: below `X_floor` the reference model
routes the released energy to the jet channel, leaving no radiative AGN
packet; `AGN_blast` is active. The nonzero aggregate M5 source is therefore
stellar in this execution.

The A100 stream comparison jobs `443809` and `443816` both completed with
exit code 0 using the same executable, effective profile, and three GPU
brokers; they differed in one versus three CUDA streams. Their application
times were 532.148 s and 436.964 s, respectively, a provisional 17.9% decrease
for three streams. This is not an accepted speedup: the final-dump comparator
reported `M5_COUPLED_FINAL_PAIR_FAIL_NUMERIC_TOLERANCE` at `rtol=1e-12`,
`atol=1e-25`. Of 449 datasets, 111 were bitwise equal and 115 exceeded the
strict tolerance (gravity 4, hydro 96, particles 3, sinks 10, SNRT 2); clock
checks had no tolerance violations. Consequently the timing remains
conditional until the divergence is understood.

After both pair evaluations completed, the three `output_00001` directories
for jobs `443660`, `443809`, and `443816` were removed. Their 42 exact files,
including the three 101,372,728-byte HDF5 dumps and SHA256 hashes, are recorded
in the cleanup manifest. Effective inputs, logs, result records, both
comparison JSONs, and build identity remain in place.

Those two GPU jobs used a pre-receipt binary, so their logs cannot establish
zero AGN packets. The current CPU binary from `1249105` contains the receipt
and job `1249133` exercised it. A first runner check, job `1249152`, stopped
before simulation because the AGN-only Chabrier history (`imf_id=2`) was
incorrectly checked against the effective-profile Kroupa id (`imf_id=1`). The
effective namelist was generated but no simulation started and no raw output
was written. The check is now profile-specific. Corrected CPU integration job
`1249156` completed four coarse steps on grammar `normal` in 1:36, exit code
0, with no raw outputs. Its result distinguishes the live state as
`JET_mode_zero_radiative_packets`, with zero accepted packets, photons, and
AGN photon energy; three sidecar records show positive BH accretion and
luminosity in JET mode. The run exercised 256 material-cell assertions, and
all recorded input hashes verify.

The older AGN-only CPU integration `1244059` contains a positive `THERMAL`
accretion state and a nonzero M5 source ledger, but predates the new receipt.
Current-source CPU job `1249188` closes this evidence gap: at coarse step 2,
the sidecar records positive accreted BH mass, `feedback_mode=THERMAL`, and
`Lbol=3.05619e40 erg/s`; the M5 receipt accepts one AGN packet representing
`1.90503e62` photons and `7.59140e51 erg`. The source-active photon and energy
ledger residuals are `3.15e-17` and `3.24e-16` in their printed conventions.
The three-step run exits successfully with the runner's AGN-only PASS sentinel
and no raw dumps. This confirms the thermal AGN receipt on the current CPU
binary under reference controls; it does not qualify the physical spectrum or
the stellar-plus-thermal-AGN combination. Do not treat the effective-profile
pass sentinel as an AGN-specific photon-deposition audit.

The first combined current-source attempt, job `1249312`, reached the
`Run completed` marker and formed stars, but its post-run assertion compared
physical AGN photon counts against the normalized code-unit M5 ledger source.
Those values are dimensionally different, so the runner failed without
indicating a RAMSES failure. The attempt produced no raw dumps and its logs are
retained. The check was corrected to use distinct source receipts in matching
units.

Build `1249326` compiled commit `e96ee64e49ae9131a1deaeb5d4634ff8efe44293`
with `CHIMES_CUDA=0`, and produced
`M5_MOMENTUM_CLOSURE_BUILD_PASS`. The executable contains both source receipt
labels and has no CUDA runtime dependency. Its dependent CPU run `1249327`
completed four coarse steps on grammar `normal` with exit code 0 in 1:45
Slurm elapsed time (90.589 s application time). It used one MPI rank, eight
OpenMP threads, M5 order 5, OpenMP transport/material, CPU CHIMES, live dust,
and the effective stellar population plus the thermal-seed AGN sink fixture.
The external Haardt-Madau field was disabled to isolate the two source paths.

At one M5 commit, the new receipts report 64 stellar sources producing
`7.71699705e66` photons and `9.29350206e55 erg`, alongside one accepted
thermal-AGN packet producing `1.56436650e70` photons and `6.23389269e59 erg`.
The runner reports `overlap_commits=1`; the same commit has nonzero gas and
dust heating, with printed photon-number and energy ledger residuals
`5.2640467126858e-20` and `1.0582739475971e-17`. The run exercised 256
momentum assertions, completed with the combined-source and integration PASS
sentinels, and retained zero raw outputs. Its effective namelist, launch
record, input hashes, sidecar and logs remain under
`.m5-effective-dust-chimes-1249327/`.

This reference fixture has `cosmo=.false.`. The legacy SFRD diagnostic divides
by a cosmological box volume and consequently prints NaN/Infinity in this
non-cosmological setup; SFRD is not an evaluated quantity here. The run's gas
state checks remain clean, and the SFRD line is not used as a success metric.
This source-attribution pass establishes simultaneous source routing and
accounting under reference controls only; it does not qualify either spectrum
or production physics.

The current evidence supports: (a) stellar SED execution, (b) live BH
accretion and luminosity, (c) a code-wired AGN-to-M5 path, and (d) a successful
aggregate M5/dust/CHIMES transaction with simultaneous stellar and thermal
AGN photons. It does not support a production claim: the spectra and dust
controls are reference-only, and physical package admission remains blocked by
the existing manifest audits.

The read-only `simulation_ready_preflight.sh --mode physical-production`
check was run with the CPU smoke binary from `1249326` and effective namelist
from `1249327`. It exited 2 with `production_gate_pass=false`, as required for
these reference artifacts. The release manifest still has
`status=operator_approved_effective_model_pending_release`, no release commit,
and `launch_admission_manifest=false`; it lists three unresolved release
conditions. Per-release license/provenance records are absent for ordinary
yields, the stellar SED, AGN groups, dust optics and SNIa, while data-use terms
for CHIMES, Draine D03, KL16 and PARSEC VMS remain unresolved. The probe also
rejected the smoke binary, namelist, environment and source snapshot as not
matching that release record, and found the repository dirty. These identity
mismatches are expected for the smoke run; the physical approval and source
terms remain genuine release blockers. See
`manifests/native_m5_release_candidate_2026-10-12.json` for their recorded
scope and exact conditions.

## Bounded verification results and current limits

The remaining source-routing checks are complete. Jobs `1249156` and
`1249188` separately verify the effective-population JET zero-packet and
AGN-only THERMAL receipt cases. Job `1249327` verifies simultaneous stellar
and THERMAL AGN source receipts on the CPU path.

The CUDA comparison also completed. Typed-A100 build job `443842` completed
successfully (one A100, four CPUs). Its dependent run `443906` completed on
one A100 with eight CPUs, `level_queue`, one GPU broker and one stream. The
AGN-only reference profile passed with no raw output; CHIMES recorded 781,416
GPU calls, 26,525 CPU fallbacks and zero device errors. The accepted AGN
receipt was one packet representing `1.90503e62` photons and `7.59140e51 erg`.
The printed SNRT energy-ledger residuals were `4.15e-17` and `-8.86e-17` for
the AGN-source step and its following transport step. This verifies execution
of the GPU CHIMES path and AGN-to-M5 receipt under reference controls, not
production physics.

Later CPU/A10 jobs `1249452` and `444137` used byte-identical effective
namelists (SHA256
`b63d24b1b9274f03151377b0e3145f2fbc7310e5ac0f05a0b3ba8554dfc1a349`) and
both passed the same three-step AGN-only source profile without raw output.
The CPU run used Grammar `normal` with eight CPUs and took 52:10; the A10
hybrid run used eight CPUs plus one A10 and took 1:08:57. The hybrid run was
about 32% slower in these allocations. Since the jobs ran on different
systems, this is an end-to-end configuration comparison rather than a
same-node hardware benchmark. It provides no evidence to make hybrid CHIMES
the default; CPU remains the default route.

The older A100 stream-count pair `443809`/`443816` remains a failed strict
final-state comparison (`M5_COUPLED_FINAL_PAIR_FAIL_NUMERIC_TOLERANCE`). Its
provisional timing difference is not an accepted speedup. Its old binaries
also predate the distinct source receipts, so those runs do not establish
AGN-specific deposition.

The current CUDA build and combined-source integration close that gap for the
reference profile. Typed-A100 build job `444209` built commit
`ccfe43b2b35224ff50b812feb0dd09782106278f` with `USE_CUDA=1` and completed
with `M5_MOMENTUM_CLOSURE_BUILD_PASS`; the executable SHA256 is
`ba2caceb023c86d933e6ce45199620d10a7bb21597d04490241cde8ebe7b38b0`.
Dependent job `444210` completed on one A100-SXM4-80GB, one MPI rank and eight
OpenMP threads. It used M5 order 5, live two-size dust, CHIMES non-equilibrium
cooling, one GPU broker/stream, and the asynchronous hybrid CHIMES backend.
The effective namelist exactly matches the CPU combined-source fixture
`1249327` (SHA256
`01e623647752e38b00ea6dde9fb2be95925c923a81ac5453e5dc4f93c67a4104`). The
four-step run completed in 128.118 s application time with 256 momentum
assertions, the combined-source attribution PASS, and the integration PASS.
At the joint source commit it recorded 64 stellar sources and one accepted
thermal-AGN packet: `7.71699705e66` stellar photons and `1.56436650e70` AGN
photons (`6.23389269e59 erg`). The M5 photon/energy ledger residuals printed
for that source step were `1.4361e-18` and `2.0204e-17`; for the following
transport step they were `3.5599e-19` and `-9.6428e-18`. These are subsystem
ledger values in the runtime's printed conventions, not a global gas+stars+
dust+radiation closure. The source step also records nonzero gas and dust
heating.

The hybrid path made 43,014 GPU RHS calls in 2,117 batches, with zero CPU
fallbacks, zero device errors, two peak in-flight batches, and one stream.
Total cooling took 126.969 s (99.1% of application time). The CPU fixture
`1249327` took 90.589 s, including 89.195 s in cooling; although the effective
namelists match byte-for-byte, the runs were on different systems, so this is
not a controlled same-node speed benchmark and does not establish a GPU
speedup. It does show that the hybrid CUDA path is active and numerically
completes this integrated reference transaction. Output clocks were audited
in the effective namelist (`noutput=1`, future `aout=2.0d0` and `tout=1.0d30`,
`nstepmax=4`, `foutput=fbackup=1000000`); no raw output directories were
created. The effective run is preserved at
`.m5-effective-dust-chimes-444210/`, with `inputs.sha256` SHA256
`1487b1d1b6ff9b5458a81a635ee1efbb691200c654422c9417be63c10388923b`.
Both profiles remain reference controls, not production qualification.

The RAMSES-wide `mcons` and `econs` printouts are not total feedback budgets:
`mcons` sums gas density over hydro leaves, while `econs` uses kinetic and
potential terms plus cosmological work and omits internal energy and explicit
radiation/feedback source terms (`../hydro/courant_fine.f90`,
`../amr/update_time.f90`). A nonzero `mcons` during star formation therefore
does not by itself imply loss of total baryonic mass; this diagnostic does not
close gas+stars+sinks. Likewise `econs` cannot certify total gas+dust+radiation
energy closure in these nonadiabatic runs. Use the subsystem ledgers for their
declared scopes and do not claim a global closure from these fields.

Production admission remains blocked as recorded above. The reference
spectra/dust controls and unresolved package provenance/data-use terms are
not promoted by these successful integration runs. The release manifest
must continue to reject physical-production admission until those conditions
are resolved.

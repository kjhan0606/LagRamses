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

## Remaining bounded verification

Jobs `1249156` and `1249188` separately verify the effective-population JET
zero-packet and AGN-only THERMAL receipt cases. Job `1249327` now verifies
simultaneous stellar and THERMAL AGN source receipts on the CPU path. CUDA
build job `443842` requests one typed A100 and remains pending with
`Priority` (last checked 2026-10-31). Its dependent runtime job `443906` is
queued with
`afterok:443842`; it holds M5 transport and dust to OpenMP while selecting
hybrid CHIMES (`level_queue`, one GPU broker, one CUDA stream). That comparison
must confirm actual GPU CHIMES work as well as the source/ledger checks; it is
not a production-physics approval. Keep physical admission blockers intact.

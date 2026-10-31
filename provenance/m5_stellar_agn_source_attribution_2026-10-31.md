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
  energy source ledgers. Its HDF5 output was retained for the matched-pair
  evaluation.

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

The runner's effective-population result string currently hardcodes
`AGN_state=not_requested`, while its AGN-only path checks for a luminous
accretion state. Thus the result record understates what the effective run
actually configured and observed.

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

The older AGN-only CPU integration `1244059` contains a positive `THERMAL`
accretion state and a nonzero M5 source ledger, but predates the new receipt.
Thus the new receipt has been compiled and exercised in the JET branch, while
a positive receipt from the THERMAL branch remains to be confirmed on the
current binary. Do not treat the existing effective-population pass sentinel
as an AGN-specific photon-deposition audit.

The current evidence supports: (a) stellar SED execution, (b) live BH
accretion and luminosity, (c) a code-wired AGN-to-M5 path, and (d) a successful
aggregate M5/dust/CHIMES transaction. It does not support a production claim:
the spectra and dust controls are reference-only, and physical package
admission remains blocked by the existing manifest audits.

## Next bounded closeout

After the queued three-broker/three-stream performance comparison no longer
reads the current runner, make its result fields distinguish stellar source,
BH accretion/luminosity, mechanical AGN feedback, and accepted AGN photon
packets. Verify that the current-build THERMAL AGN profile produces a
positive source receipt. Keep physical admission blockers intact. This is a
reporting and verification correction; it does not alter the selected
physical model.

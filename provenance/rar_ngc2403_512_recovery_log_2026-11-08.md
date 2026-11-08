# Historical NGC 2403 force-log compatibility

Syntax analysis job `455930` failed after 1m44s on 2026-11-08 at 19:53 KST.
It preserved the two completed snapshots from run `451793`; no simulation
was rerun. The analysis had passed its force/source numerical thresholds
before failing in the late diagnostic-log parser, but produced no final
result artifact, so no successful completed audit is claimed for that job.

The saved run used source `65eeb24ddb43e22e6b6acc48b6758c386d506413` and binary
SHA256 `7baf42cc884f33e9babb09b21aca5c46c380b3c32186cd562d20971c6414061c`.
Its log contains two source diagnostics, at steps 0 and 1 on level 9,
including nonzero phantom density and all 512^3 cells. It predates the
post-second-Poisson and packed-FFT diagnostics expected by the newer parser.
The analysis also used the old two-result parser return signature although
the current parser returns three collections.

The recovery explicitly allows the historical source-only format only when
both newer diagnostic families are entirely absent. Default parsing remains
strict; partially recorded newer diagnostics still fail. The analysis records
that post-Poisson/FFT telemetry is absent and does not claim source retention
was directly logged. It independently reconstructs both force fields and the
Newtonian Poisson source from the full stored snapshots, with unchanged
1e-6 source and force thresholds. Source log preflight now occurs before GPU
array allocation, and all three parser results are unpacked explicitly.

Runtime validation was pending at implementation commit `33b2a8f`. Syntax
job `456010` subsequently passed using a new postprocessing directory. The
failed analysis and original run remain unchanged; neither received a new
success marker.

The two 512^3 snapshots each supply 24,089,384 force comparison samples.
Newtonian median/p90 relative force differences are 1.03e-14 / 2.28e-14;
QUMOND differences are 1.68e-14 / 3.15e-14. The Newtonian source RMS relative
difference is 5.25e-13. These compare the compiled fields with an independent
reconstruction of the same periodic discrete equations and source. They do
not measure continuum or isolated-galaxy accuracy. The result explicitly
records two source-only diagnostics and zero post-Poisson/FFT diagnostics.
The final result, hash manifest and source identities are preserved under
`provenance/results/rar_ngc2403_recovery_456010/`.

The subsequent Grammar NFW holdout runner now accepts this separate completed
force-analysis directory. It requires both cases, all source hashes, matching
snapshot paths/digests and the unchanged 1e-6 force/source thresholds. It
verifies the snapshot digest again when reading forces. Analysis code is
exported from an explicit Git revision before execution. It never fabricates
completion files in the failed original wrapper directory.

The holdout question remains a fixed-source rotation-curve comparison: fit
the two analytic NFW parameters on 48 of 64 selected radii and evaluate the
other 16, in four contiguous folds, with fixed QUMOND parameters. This is
not a live Lambda-CDM cosmological calculation, and the two models have
unequal fitted freedom. The spherical isolated NFW acceleration is added
to the periodic baryonic field; no full periodic halo Poisson solve is claimed.

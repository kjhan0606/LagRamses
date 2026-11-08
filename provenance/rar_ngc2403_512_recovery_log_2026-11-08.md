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

Runtime validation of this compatibility correction is pending at this
commit. A new Syntax GPU postprocessing directory is required; preserve the
failed analysis and original run without adding a success marker to either.

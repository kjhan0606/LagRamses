# AGN uncapped-energy replay correction: compiled evidence

After the sink-cloud ordering correction, a ten-step MPI fixture still diverged
at the last thermal epoch by -7.117989499580268e-4. The old cell helper rebuilt
uncapped total energy from kinetic and internal components. Roundoff produced
a cumulative 4.9630836753181646e-24 code deferred-energy reservoir in two ranks.
Any positive saved reservoir selects replay, displacing a fresh thermal event.
One rank deposited 6.312606409329031e-7 code energy at the affected epoch, while
two ranks deposited approximately zero. Local budget tests alone did not expose
this branch change; matched epoch and spatial comparisons did.

Commit 4c13c58 retains trial_energy exactly, with deferred_energy=0, when the
temperature cap is inactive. Genuine capped deposition keeps the original
limited-energy/deferred-energy formula. There is no small-energy threshold or
changed feedback efficiency.

## Regression and integrated results

Isolated Grammar build 1275251 completed from 4c13c58 in 3m02s.
Compiled regression 1275256 reproduces the expected uncapped-energy failure
with the previous helper (exit 128) and passes with the new helper (exit 0):
925 uncapped cases with momentum and varied scale, plus an active cap requiring
positive deferred energy and total-energy closure. Both input-object and output
artifact hashes verify.

Ten-step one-/two-rank runs 1275252/1275253 completed in 2m30s/1m20s.
Cloud/position audit 1275257 finds 2109 moving clouds at all ten steps in both
runs. Spatial AGN comparison 1275259 finds 504 epoch/cell entries in each run,
all common, with L1 accepted-energy difference 2.1341259803868114e-21 code
energy. Thermal epoch comparison 1275254 finds a maximum absolute relative net
difference of 1.3565721501602264e-15, and last-epoch difference
-6.347914737542832e-16.

Combined exchange audits 1275258 (one rank) and 1275255 (two ranks) complete.
Both record 504 positive AGN receiver transactions and exactly zero deferred
energy. Requested/accepted totals are 2.6893390264983295e-6 and
2.68933902649833e-6 code energy respectively. Stellar mass and eleven-element
transaction checks also pass. All five integrated artifact manifests verify.

## Limits and follow-up

This establishes the measured correction within a ten-step synthetic uniform
CPU QUMOND hydro fixture. Stars/AGN/metals/cooling/heating are enabled; dust and
RT are disabled. The source bundle and artificial star-formation/BH settings
remain exploratory. These tests do not establish AMR-boundary accuracy, global
mass/element/energy closure, donor withdrawal and sink-generation budgets,
restart parity, a viable relativistic cosmology, or a galaxy RAR result.
Do not extend this diagnostic repeatedly as a substitute for those requirements.

The two-rank diagnostic directory previously reached about 300 MB, exceeding
the earlier 256 MiB estimate without exhausting the filesystem. The runner now
reserves 512 MiB plus a 1 GiB free-space margin, with no full snapshots requested.
Existing raw records, including the earlier divergent runs, remain preserved.

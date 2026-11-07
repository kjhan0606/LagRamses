# Whole cooling/heating call: global energy measurement

Source 7ae88e8 adds opt-in global before/after inventories at the main
cooling_fine call sites, including the later floor and velocity corrections
that occur after the accepted thermal-update records. Physics is unchanged.
Grammar build 1275327 completed in 2m51s; three-step one-/two-rank runs
1275328/1275329 completed in 21s/16s. Reconciliation 1275330 completed in 3s,
with its artifact/input SHA-256 manifests verified.

At each of three epochs, the analysis requires 32768 accepted thermal records
and matching owned leaf-cell coverage. Effective RHS and accepted increments
are summed in code units at the same epoch, then compared with the independent
whole-call gas-energy difference. It records discrepancies as unaccounted final
corrections rather than assuming local source agreement implies full closure.

For one rank the whole-call changes are 2.7952549956309228e-14,
-1.0827161610174939e-8 and -6.319912094376833e-8 code energy. The differences
from summed accepted updates are -2.158386545557106e-25,
5.955700410381799e-23 and -3.441071348220595e-22. For two ranks these
residuals are 3.78725809709366e-25, -9.926167350636332e-24 and
-4.632211430296955e-22. Both signs of effective heating/cooling are exercised.
All epochs pass the scale-aware comparison. No unexplained final-call energy
increment is resolved in this fixture; active nonzero floor corrections are
not established simply because the call was bracketed.

This complements the accepted stellar/element transport, global mass/BH
radiation and thermal AGN source-to-gas measurements. It is a synthetic
uniform periodic CPU fixture with exploratory SSP/BH/star-formation settings,
not a galaxy result, general AMR/MPI certification, pure emissivity measurement
or cosmic gravitational/expansion energy conservation. Dust and RT remain off.

The pending NGC 2403 fixed-source GPU calculation uses a separate immutable
build at 65eeb24ddb43e22e6b6acc48b6758c386d506413. Its clean checkout and binary
SHA-256 7baf42cc884f33e9babb09b21aca5c46c380b3c32186cd562d20971c6414061c
were rechecked. This hydro source branch must not be substituted for that
build when interpreting its eventual field results.

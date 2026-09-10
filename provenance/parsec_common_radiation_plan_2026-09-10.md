# PARSEC common high-mass radiation/feedback bundle

Final objective: physically justified, simulation-ready RT/stellar/AGN
feedback/dust. All remaining medium work is preapproved; continue without
per-helper approval waits. Ask Q-GOAL FIRST, Q-LEAN SECOND in the single
Fable plan review. Driver performs the end evaluation.

## Runtime capability and scope

Connect the existing PARSEC v4 14--600 Msun, Z=.008/.014, nonrotating
five-fate feedback population to radiation from the SAME initial-mass nodes,
metallicities and evolutionary tracks. Retain the .08--600 initial-mass
denominator and all four existing individual-star IMF choices; only the
14--600 contribution is supplied. Do NOT renormalize the covered subset,
fill the lower-mass radiation with zero-labelled physical data, or claim
closure of the full galaxy SSP. This is an explicit high-mass comparison.
Old BPASS independent and legacy SED options remain unchanged; this new
source must not be mixed silently with BPASS.

Public source: Costa et al. 2025, https://arxiv.org/html/2501.12917v2,
https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/ .
Downloaded matching Z photon archives accompany the already pinned track
and ejecta archives. Photon columns are instantaneous cumulative Q above
HI, HeI, HeII, OII and Werner thresholds, NOT five disjoint photon groups,
time-integrated Q, or a supplied nine-group spectrum.

## Declared spectral closure

Use differences of the five cumulative Q constraints for disjoint intervals
bounded by nominal 11.2,13.6,24.6,35.12,54.4 eV and infinity. Within each
interval use the positive photon Planck shape at the tabulated effective
temperature, normalized to that interval's actual photon rate. Integrate
over the true nine SNRT edges, not integer reassignment to a nearest group.
Energy moments follow this reconstructed spectrum, NOT nominal group means.
The remaining bolometric luminosity below 11.2 eV is distributed with the
same Planck prior. Reject negative residual bolometric energy, negative Q
differences, or unsupported source states; measure actual feasibility before
accepting this closure. Record that nominal source thresholds are rounded
atomic edges and this is NOT the original NLTE atmosphere spectrum. Published
five Q values and luminosity constrain the result; the within-band shape is
an explicit additional model, not inferred stellar feedback energy.
Retain/report radiation below .01 and above 10000 eV as escaped spectral
tails, not gas heating or reassignment into the boundary groups.

Use each track's own photon times; hold its first photon state back to age0
as a named pre-main-sequence approximation, clip at the exact feedback
terminal age, and only extend the last photon state to that age when the
source sampling mismatch is small and explicitly measured. No post-death
emission and no arbitrary SSP age extrapolation. Integrate positive Q/E
rates into cumulative histories, with a measured compression error if useful.

## Native wiring (not just an offline Python exercise)

A small per-source-node cumulative Q/E reader/integrator uses the existing
native IMF mass-fraction routine, nearest mass-cell boundaries INCLUDING
the existing 40-Msun seam, and same-age linear-Z mixtures. Thus its weights
are exactly feedback's initial-mass fraction divided by source-node mass.
Do not add per-cell spectral history. Check the loaded prepared feedback
table's model coordinates, M/Z/terminal-age/fate arrays and configured IMF,
channel windows before enabling this source; require actual v4 feedback,
not just a matching filename. Bind every consumed source value and its model
metadata in existing MPI/restart SED identity. Preserve legacy identities.
Integrate through existing stellar photon injection with actual Q/E and
existing spectral transport. No new RAMSES namelist key is planned; use an
explicit source-file selector and document admission. If a NML semantic
changes, update both generators. Preserve Makefile VPATH.

## Lean completion evidence

Reuse native stellar source tests to check physical Q constraints, positive
Q/E and bolometric closure with tails; native/source-independent IMF weights;
no photons after source death; interval telescoping and Z mixture;
rejection of population/feedback/source mismatches without partial publication.
One small MPI2/OMP2 continuous/restart evolution with actual PARSEC feedback
and radiation. Reuse existing builds and logs where possible; retain compact
evidence and remove only evaluated raw dumps. No new testing framework,
parameter campaign, certification gate or external end audit.

This closes the common HIGH-MASS population connection, not unmatched low-
mass tracks, microscopic Ia, pulse timing or universal atmospheric accuracy.

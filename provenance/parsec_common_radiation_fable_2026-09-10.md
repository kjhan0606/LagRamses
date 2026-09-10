# Fable plan review — driver disposition

Read-only `claude -p --model fable --permission-mode plan`, Read/Grep/Glob,
no `--bare`, no files/jobs/messages/subagents. Completed with exit0.
This is the driver's summary, not a verbatim transcript.

Verdict: **Proceed with corrections**. Q-GOAL yes for bounded high-mass
population; Q-LEAN trim. Keep this off by default and report HIGH-MASS ONLY.
Reuse offline spectral conversion, existing cumulative compression1e-4,
and extract/reuse the native source mass-cell edge builder. Preserve native
IMF weighting and lifetime-resolved integration. Bind actual feedback at
startup, not lazily at first radiation transaction. Reuse MPI/HDF5 identity.
Compare reconstructed group means to BPASS as a sanity comparison, not as
a matched physical truth (the two populations remain different).

Accepted the implementation/lean recommendations. The proposed conversion
was already offline. The source photon archives and formats were inspected
and pinned before use; see the source record. All actual Q differences are
nonnegative and bolometric residuals are positive. Therefore do NOT add the
suggested clipping or lower-edge-energy fallback: neither is needed and the
latter would silently change the selected spectral closure. Reject malformed
data. Fable's general 10--20% blanketing estimate was not independently
established for this package, and is not a measured accuracy claim. Its
statement that transport coefficients are fixed is outdated for the new
HHe/D03 spectral modes; actual Q/E must feed their evolving node absorption.

## Bounded source repair after inspection

Combined archive `all_photons.zip` SHA256
`df31fd6b9fe49abb776ba1701fcc1aef19a407c3ab95768c58537f11b4175f53`,
115670480bytes, contains byte-identical individual photon archives. Two
Z.008 photon tails miss material temperature changes. Do not treat them
as small last-state holds. Use the actual full PARSEC track's L and Teff
with an explicitly named `Planck_track_tail_v1` spectrum for ALL missing
end intervals. This supplies no invented measured Q; it is a separate
blackbody closure, with per-node missing duration and bolometric contribution
reported. At the join integrate each side separately, allowing a rate jump
without a fabricated smoothing interval. Only the measured photon-time
coverage uses the five Q constraints. Pre-first-photon hold remains named.
Both policies and every consumed cumulative Q/E value bind to restart.
This bounded source repair remains part of the reviewed comparison bundle;
it does not make a production/full-SSP or exact NLTE atmosphere claim.

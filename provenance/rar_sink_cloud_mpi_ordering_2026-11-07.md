# RAR sink-cloud MPI ordering correction

The no-dust/no-RT QUMOND hydro fixture exposed a decomposition-dependent sink
trajectory. Under the former create_sink rebuild order (make_tree_fine,
kill_tree_fine, virtual_tree_fine), migrated parents arrived after their level's
descent. In the two-rank run, 401 clouds remained at level 3 and 362 at level 4,
with only 1346 of the 2109 clouds at uniform base level 5. The one-rank run had
all 2109 at level 5. Particles were not destroyed; the coarse clouds were omitted
from the base-level move statistics used for the sink mean position.

Commit c01513c orders this loop as make_tree_fine, virtual_tree_fine,
kill_tree_fine. Received particles are therefore available before descent.
No force law, feedback efficiency, source table or timestep was changed.

## Compiled evidence

Grammar build 1275215 completed in 2m52s. Three-step runs 1275216 (one rank) and
1275217 (two ranks) completed in 23s and 16s. Their cloud/position audit 1275219
completed: both move-stat samples contain 2109 clouds at steps 1–3; raw mean
positions agree within approximately 3e-15 in unit-box coordinates. AGN spatial
comparison 1275226 finds exactly 112 receivers in each run, all common, with
L1 accepted-energy difference 3.2053248736429822e-22 code energy. Previously,
only 16 receivers were common and the difference was 4.1399567901722894e-7.

Epoch thermal comparison 1275227 finds relative net differences -9.9355e-16,
0, and -1.3566e-15, versus the former third-epoch 7.7313e-4 difference.
Combined local-exchange audits 1275228 (one rank) and 1275220 (two ranks) both
completed. The latter has 328 stellar transactions, 3608 element records and
112 positive AGN receivers; all requested 2.3502550743283348e-7 code AGN energy
was accepted, with zero deferred energy. All five audit artifact manifests pass
SHA-256 verification.

The initial census analysis 1275205 rejected stale finer-level list headers
in the before_rebuild phase. These are expected after gathering particles to
level 1 and are not a valid all-level tree. Reanalysis 1275214 retains these
observations while requiring valid unique entries after rebuilding and before
movement. Preserve the original failed analysis; this correction does not
certify arbitrary transient lists as valid.

## Scope

These are synthetic three-step uniform-grid CPU tests with seeded BHs and
artificially strong star formation. They establish the measured correction
within that fixture, not general AMR/MPI correctness, calibrated galaxy formation,
restart parity, donor/global energy accounting, or a cosmological RAR result.
Detailed source/input identities and measurements are preserved under
/gpfs/kjhan/RAR_QUMOND/audits/ and corresponding builds/runs.

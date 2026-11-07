# Native stage-history integration: first bounded implementation

Source 86d2b2c links separate current/prior Newtonian and QUMOND potential buffers
into native coarse/fine interpolation. The shared phi/phi_old solver storage and
ordinary cold-restart correction remain. An explicit selector is entered before
the Newtonian solve and second QUMOND solve; each solved stage is captured after
its force construction. Interpolation reads the selected stage's own pair.
Per-level physical times, startup validity and coarse epoch are recorded; first
capture sets that stage's prior equal to its current field. FFTW's periodic base
gauge is checked. Cell/grid topology and coordinates are frozen and checked.

Activation is opt-in RAR_QUMOND_STAGE_HISTORY=1. The implemented path rejects
MPI, cosmology, restart, load balancing, subcycling, CUDA builds, non-static or
non-two-level configurations, non-FFTW base solves and BIMOND. Changed allocation,
refinement/connectivity or grid positions abort rather than silently reusing
history. These restrictions bound this first implementation; they are not the
final lifecycle contract. Coupled source evolution remains unvalidated.

Grammar isolated build 1275333 completed in 2m55s. Initial execution 1275335
completed its old/new off cases but rejected the enabled input: nsubcycle=1
sets only the first supplied array entry, leaving the second level's default 2.
That failure and dependency-pending audit 1275336 are preserved; they are not
physics results or a successful validation. Corrected execution 1275343 uses
nsubcycle=1,1 explicitly in all three matched inputs, retaining all old files.
It completed old_off/new_off/new_on in 33s. Each case produced exactly two
snapshots, with six total; the 2 GiB reserve exceeded the roughly 547 MB estimated
snapshot total. Per-case effective inputs, source revisions and hashes remain.

Audit 1275344 completed in 5s. All 16 saved gravity arrays (two outputs, two
levels, phi and three force components) are bitwise equal between old_off and
new_off. Enabled fields are finite with maximum absolute change 0.2840484525899871
code units relative to new_off; this shows activation, not improved accuracy.
Twelve captured records cover both stages and levels, with stage-specific first
capture and advancing current/prior times. The result/artifact hashes verify.

This implements the guarded state/selection part of reviewed bundle 1. General
MPI ghosts, grid growth/regrid, relocation, persisted stage restart, subcycling
and arbitrary source epochs still require lifecycle integration. The phantom
source still skips missing fine neighbors. No shared-face composite source,
composite Poisson residual or synchronous hierarchy-force publication has been
implemented. Therefore this path is experimental and cannot certify AMR QUMOND
forces, evolving cosmological hydro, an RAR result or a relativistic MOND model.
The next reviewed spatial bundle remains necessary; do not remove these guards
merely to enable unsupported calculations.

# Native H/He + D03 grain spectral absorption and scattering

Final goal: simulation-ready, physically justified RT/stellar/AGN feedback
and dust in lagRamses. Operator preapproves all original medium-term tasks
and continuous execution. This bundle addresses group8's grain spectral
coupling. It does not claim to close CHIMES's different reaction-table
coupling, PAH high-energy physics, or all medium-term work.

## Q-GOAL first, Q-LEAN second

One plan review before implementation, driver end evaluation. Reuse existing
native transport/material tests and one tiny MPI/restart integration. No
new framework, job campaign, persistent spectral tensor, or approval wait.

## Runtime result

Add one explicit `SNRT_SPECTRAL_MODEL=hhe_d03_maxent64_fs2010_v1` comparison.
Reuse the existing positive64-node N/E reconstruction, Verner H/He and
node-resolved FS2010. Extend each node's competing opacity to four actual
graphite/silicate mass bins, with original D03 dielectric data, radii .01/.1
micron, densities2.2/3.8 and the existing1/3--2/3 graphite prescription.
Use the existing offline Mie builder to compile node absorption and transport
scattering Qsca*(1-g), not a rescaled representative-energy opacity. Source
range/units and all used coefficients must be bound in restart identity.

Per-cell inputs are only four mass columns plus immutable node coefficients,
not cell*direction*node persistent storage. Absorb once with total node
opacity; finite H/He caps return the rejected packets with their own E/N.
Dust is uncapped and only its actual captured energy heats grains. Rejected
H/He packets are not reassigned to dust a second time in this operator split.
Apply elastic delta-isotropic angular redistribution at those same nodes
after absorption, conserving N and E separately. Do not count scattering as
absorption or heating. Reconstructing the survivor spectrum for this split
is a declared two-moment closure, not exact spectral recovery.

Wire actual grain energy into the existing transactional material/IR solver;
fixed group means must not re-enter heating. Shared four-bin thermal material
temperature remains the existing comparison. Retain original fixed, HHe-only,
and FS-only modes and their defaults/checkpoint identities.

Initial admission: DUST_LIVE, composition/two-size D03 model, no CHIMES,
Fe, PAH, sublimation or relative dust dynamics. HHe/CIE-metal cooling remains
available. Relative-motion spectral force and the157-species reaction-table
closure are separate remaining work, not silently approximated by frozen
tables. No new namelist key; update both generators' guidance if necessary.

## Minimum evidence

Extend native transport smoke: zero grains equals HHe-only; zero gas gives
analytic grey extinction; actual D03 soft/hard packets differ; finite atomic
inventory; photon/energy closure; elastic scattering conservation and angular
relaxation; invalid coefficients fail atomically. Check node quadrature at
source edges against denser offline Mie integration, report limitations
rather than claim universal64-node convergence. Intel/GNU where available.
One small continuous/restart CPU MPI2/OMP2 evolution with actual D03 grains,
actual source photons and nonzero absorption. Check total radiation/material
energy and finite physical state, exact restart; retain compact evidence and
delete only evaluated raw outputs under the standing retention instruction.

Primary data/method: https://www.astro.princeton.edu/~draine/dust/dust.diel.html
and https://arxiv.org/abs/astro-ph/0308251 . D03 graphite approximation,
20K dielectric coefficients, finite nodes and delta-isotropic angular closure
remain declared limitations. No claim that these grains are PAHs or that
this implements photoelectric escape from grains.

## Pre-implementation source-resolution measurement

Direct Mie integration using the staged original dielectric files (miepython
3.3.0) found64-node uniform-dE absorption moments differ from1024-node by
4.62% for small graphite in54.42--500eV. At128 nodes, comparison with2048
nodes is0.541% there (large graphite0.172%); all tested500--10000eV bins
are below0.26%.256 nodes is not monotonically better at an edge. Therefore
use128 nodes for this NEW option, named
`hhe_d03_maxent128_fs2010_v1`, while preserving the old64-node models.
This is a measured resolution choice within the same bundle, not a new
gate or claim of universal edge convergence. Bind the immutable generated
node data by content SHA256, retaining the exact generated values in source;
avoid a new large per-cell or HDF5 attribute tensor.

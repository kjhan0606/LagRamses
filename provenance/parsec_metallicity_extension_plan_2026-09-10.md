# PARSEC common high-mass metal-rich extension: one implementation bundle

Final goal: physically justified simulation-ready RT/stellar/AGN feedback
and dust. All medium groups remain preapproved; no new approval wait.
Q-GOAL FIRST, Q-LEAN SECOND, one Fable plan review and driver end evaluation.
This advances the existing missing-M/Z and common-population groups, not a
new finish-line gate or a claim that all medium groups are finished.

## Bounded deliverable

Extend the existing nonrotating14--600Msun PARSEC feedback+matched Q/E source
from Z=.008/.014 to an explicit five-node grid .008/.014/.017/.02/.03.
Reuse native v4 history, cumulative mixture in linear Z, same IMF cells,
five-fate ownership, eleven-element channels and original Q5/Planck prior.
Both material and radiation must use precisely the same chosen grid; do
not extend one while silently extrapolating/omitting the other. Old two-Z
default artifacts and behavior remain available and byte-identical.

Actual all_ejecta.zip inventory check: all45 mass nodes at each added Z
have nonnegative eleven-element residuals in wind and terminal return.
No yield renormalization or change of the published baryonic remnant is
needed. Other Z grids have negative residuals at some nodes; they remain
outside this specific addition, NOT declared absent/zero or solved.
Z=.01 is not an interpolation node in either model. Existing interpolation
from .008 to .014 is an explicit cumulative mixture, not the missing Z=.01
source values. Adding .017/.02/.03 does not claim improved early-universe
low-Z coverage, full SSP radiation, binaries or resolved pulse chronology.

Physical inputs: Costa2025 https://arxiv.org/html/2501.12917v2 and author
database https://stev.oapd.inaf.it/PARSEC/Database/PARSECv2.0_VMS/.
Keep the pinned common ejecta/energy archives; acquire only three matching
track archives, and reuse matching photon archives from retained
all_photons.zip. Pin physical bytes. Never delete these original sources
as if they were disposable simulation raw output.

## Implementation and minimal evidence

One offline converter grid option (`solar_pair` default / `metal_rich_five`)
and matched SED inferred from the actual feedback manifest. Native reader
admits only these two exact grids and45 source masses each, then binds
actual coordinates/fates/ages/model to feedback as before. Remove only the
hardcoded90/two-Z limitation, not cross-input identity checks. Total225nodes
fit the existing512-node high-mass cap; no new carrier/index/restart format.
Existing fixed and named phase-velocity models remain explicit comparisons;
do not cap Eddington Gamma or fabricate missing phase velocities if a new
source fails their existing domain checks. Such a case is model-specific,
not permission to erase material or relax physical checks.

Extend existing native source fixtures to exercise new endpoints/interior
Z mixtures and reject unsupported/mismatched grids. Compare old two-Z bytes
and actual matched source budgets. One small MPI2/OMP2 radiation+feedback
integration and exact restart at a new interior metallicity is sufficient;
do not run a matrix of galaxies or add a new test framework. Before launch
report effective absolute NML/options/all output schedules and dump budget;
after driver evaluation delete exact raw HDF5, keeping inputs/logs/summary.
No RAMSES NML key or generator change expected; update both generators only
if implementation actually adds a native NML selection. Preserve VPATH.

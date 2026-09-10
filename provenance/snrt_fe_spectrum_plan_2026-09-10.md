# Static Fe comparison: spectral RT connection

Implemented with the [driver's review disposition](snrt_fe_spectrum_fable_disposition_2026-09-10.md):
static CHIMES-free layout is now explicit and both namelist generators are
updated. General Fe stellar/AGN irradiation remains unadmitted; the live
source is explicitly synthetic sub-eV. See the implementation evaluation.

Final goal: usable, scientifically explicit native RAMSES RT/feedback/dust,
not more instrumentation. Medium item8 already lists Fe spectral coupling.
Operator preapproved remaining medium implementations. Selective Fable review
is requested here because the primary RT/material coupling changes.

Q-GOAL first: is this existing missing connection useful toward that goal?
Q-LEAN second: can the implementation/evidence below be made smaller without
concealing a material scientific or conservation failure?

Scope: opt-in `hhe_d03_fe_maxent128_fs2010_v1`, using the EXISTING
`fe_electric_compare_v1` six-bin comparison (C small/large, silicate
small/large, Fe small/large). Do not admit complete ferromagnetic Fe optics.
Reuse Werner09/Henke causal dielectric + existing size Drude, mu=1 Mie:
electric and eddy response only, frozen dielectric, no spin absorption or
photoelectron escape. Generate 128 log nodes/group with the existing builder;
retain the old representative/IR tables and fixed-mode identities unchanged.
Measure actual Fe absorption/transport node resolution against a denser
quadrature before claiming its useful accuracy; sharp-edge convergence is
not universal. No new physical archive or simulation-scale runtime table.

Extend the existing positive maxent N/E dust operator from four to six
temporary mass columns. A single H/He/dust competing absorption event uses
the sum of actual per-bin mass opacity. Finite H/He rejection retains energy;
dust is uncapped. Existing elastic delta-isotropic node scattering follows,
without the old fixed-mean absorption/scattering paths. Actual absorbed
group energy is deposited once into the existing shared-temperature six-bin
material/IR transaction. This does NOT introduce six independent grain
temperatures, magnetic response, PAH physics, relative motion, sublimation,
or CHIMES spectral chemistry. Keep these combinations rejected and retain
all existing Fe model restrictions. Existing modes/defaults stay unchanged.

Inspection amendment: the existing Fe comparison is limited to 4 eV primary
absorption and 300 K material (dust_mass_physics), not a full-band heater.
The new nodal operator must reject nonzero absorption above 4 eV in an
Fe-bearing cell BEFORE publication; a group mean below 4 eV is not enough.
The 1--5.6 eV group's broad maxent reconstruction generally has such a tail
and must reject it, not clip/renormalize the spectrum. Existing fixed-mode
mean guard remains unchanged. A sub-1 eV live comparison can exercise the
wire; high-energy samples are rejection tests, not admitted Fe heat. This
does not close the missing full photoelectric/magnetic Fe model. If that
restriction makes this connection unhelpful, say so under Q-GOAL.

Bind both immutable nodal hashes and a new mode kind in native/HDF restart;
do not increase per-cell checkpoint carrier width. Keep the four-bin C ABI
valid, adding a bounded six-bin entry if needed. Native OpenMP fallback only;
no invented CUDA performance claim. No namelist key/meaning changes expected.

Evidence: extend existing backend tests for actual Fe-only and mixed opacity,
zero-Fe equivalence, photon/energy closure and invalid sixth-column rollback;
run existing old-mode regression and checkpoint checks. One small MPI2/OMP2
live/restart comparison with actual seeded Fe and shared IR closes the wire.
Driver evaluates; no routine external end audit. Delete evaluated exact raw
outputs, retain inputs/logs/compact results. No parameter sweep, new harness,
generic provenance gate or unrelated AMR review.

Review the plan against existing `snrt_band_spectrum.h`,
`snrt_openmp_transport.cpp`, `snrt_ramses_driver.f90`,
`dust_iron_radiation.f90`, `dust_iron_optics.f90` and
`simulation/snrt/tools/build_fe_grain_optics.py`. Identify only blocking
scientific defects versus optional refinements. Read-only; no file edits,
jobs, external messages, commits or subagents.

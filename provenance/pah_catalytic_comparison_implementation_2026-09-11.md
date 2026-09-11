# PAH catalytic H2 comparison: implementation and driver evaluation

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
User explicitly approved named approximation comparisons for bundles2--5.
**PASS for the bounded H13 abstraction connection**, not completion of
bundle2 carbon destruction, higher charges or Fe photoelectric physics.
[Plan and review disposition](pah_catalytic_comparison_plan_2026-09-11.md).

## Delivered model

`dust_pah_model='pah_h2_catalytic_v1'` is a strict superset of the existing
H2-capture comparison. H13 neutral/cation + H -> H12 + H2 uses one
barrierless site with sigma0.06 Angstrom^2 and Maxwell mean H speed.
No extra carriers: the existing3584 PAH state densities and CHIMES H/H2
receive all material changes. The closed-form finite-donor event count is
shared across both charges/excitation bins, after attachment/photo/IR and
before recombination. No equilibrium projection or unfunded clipping.

This is an explicit M13/Boschman hybrid: M13 H13 binding3.2eV, existing
H2 ground binding4.4781eV, unchanged PAH excitation and ground-state H2.
Net1.2781eV/event goes to gas heat under immediate accommodation. Neither
the outgoing H2 spectrum nor full H0--36 chemistry is resolved. The molecular
chemical-energy change closes the existing shared material/IR ledger.
See the primary references and limitations in the plan; no measured
fragment energy distribution is claimed.

Both setup frontends and their common tests include the selector, existing
build requirements and physical bounds. New reaction constants/thermal
speed and enable convention append to the restart identity only for this
selector. Old identity vectors/choices remain unchanged. Validity remains
noncosmo, coadvected, no Fe, gas10--10000K, primary<=13.6eV, fixed-group RT.

## Evidence

Private inputs/builds/logs/results: `.pah-catalytic.Av4xCW/`.

- Native finite shared donor: scarce/equal/abundant H, exact split-step
  identity for this reaction, number/charge/excitation/bond/H2 energy,
  zero dt/donor, invalid input rollback PASS.
- Native shared photo/IR/H step produces7.2036584320e-5 H2/cm3 in the
  isolated test, relative energy residual1.154e-16, H residual7.590e-19.
- One-off full split-cycle timestep sequence dt,dt/2,dt/4,dt/8 gives
  H2=1.3050940792e-5,1.3431699771e-5,1.3634508579e-5,1.3740668636e-5.
  Successive changes decrease; this is a split-error trend, not a general
  production timestep accuracy bound. Existing native dust suites pass.
- Frontends:49 PASS, one display-dependent skip (50 tests total).
- MPI2/OMP2 noncosmo64-cell run:4 steps; restart from step2 to4 PASS.
  Maximum element8.674e-16, charge2.533e-15, PAH-skeleton9.993e-16;
  IR energy residual8.696e-16. Gas H2 also evolves through CHIMES: its
  entire live increase is NOT attributed solely to this PAH reaction.
- Exact restart:7542 hydro +2 SNRT +8 gravity +12 AMR +1 header datasets
  match bitwise; compared step/time/NVAR attributes also equal.

Binary `ramses_pah_catalytic3d` SHA256
`8286a11a52589c45999a0937c16d85eaea7bdf41b6f5af0007ed1f1950cea30c`.
Intel ifx-O3/no-ftz, NVAR3771/NENER1, CHIMES/HDF5, CPU/OpenMP.
Initial clean parallel build exposed a missing pre-existing module edge:
`morton_keys.o` uses `amr_commons`. Added that exact dependency without
changing VPATH; subsequent clean-directory serial continuation succeeded.
`DEBUG=1` in the first attempt was not the Makefile's bounds-check selector;
no bounds-checked whole-runtime build is claimed.

Driver evaluation is complete. Raw cleanup is recorded in the private
cleanup manifest; retain inputs, logs, binaries, compact evaluations and
metadata. No routine external end audit was added.

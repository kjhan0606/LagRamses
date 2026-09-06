# P0 Benchmark Matrix: Multigroup S_N RT

Status: Draft v0.1  
Date: 2026-08-28

## 1. Rules

Every benchmark has one immutable YAML input, one versioned reference result,
and one conservation ledger. M1 is a comparison method, not the reference for
S_N correctness. The reference is an analytic solution where available, and a
converged long-characteristic calculation otherwise.

All tests run in both S4 and S6. S8 is required for tests B03, B04, B06, and
B07. A failure of photon-number conservation, positivity, or reproducibility
blocks subsequent physics additions.

## 2. Required test suite

| ID | Test | Physics isolated | Primary metrics | Acceptance criterion |
| --- | --- | --- | --- | --- |
| B01 | Static Stromgren sphere | HI photoionization and recombination | Stromgren radius, radial HI profile, photon ledger | radius within 2 percent of analytic value; global photon error below 1e-4 |
| B02 | Expanding I-front | time-dependent HI chemistry and reduced-c convergence | I-front radius versus time | within 5 percent of analytic/reference curve after reduced-c convergence |
| B03 | Opaque clump shadow | angular transport and absorption | shadow contrast, penumbra width, absorbed photons | S6-to-S8 contrast change below 10 percent; no negative intensity |
| B04 | Two-source crossing beams | mutual illumination and multi-source angular transport | flux vector field, central overlap intensity, shadow topology | S6-to-S8 change below 10 percent against long-characteristic reference |
| B05 | Dusty radiation-pressure slab | dust absorption, momentum coupling, IR group | absorbed luminosity and integrated force | energy and momentum ledger below 1e-4; force within 5 percent of reference |
| B06 | X-ray neutral cloud | EUV/soft-X/hard-X attenuation and secondary ionization | ionization, temperature, deposited energy versus column | S6-to-S8 and energy-group convergence below 10 percent |
| B07 | Clumpy dual-AGN nucleus | complete static solver, two AGN, dust, chemistry | escape fraction, column distribution, force, temperature PDF | all primary metrics change below 10 percent from S6 to S8 |

## 3. Canonical initial conditions

### B01 and B02: ionized hydrogen sphere

- Hydrogen number density: 1.0 cm^-3.
- Initial temperature: 100 K.
- Monochromatic HI ionizing source: 1.0e49 photons s^-1 at the box center.
- Case-B recombination coefficient: fixed at its 1.0e4 K value for the
  analytic comparison.
- Cubic domain: six analytic Stromgren radii per side.
- No dust, helium, hydrodynamics, or external background.

### B03: clump shadow

- Ambient hydrogen number density: 1.0 cm^-3.
- Spherical clump density: 1.0e3 cm^-3.
- Clump radius: 0.1 domain length.
- Plane-parallel HI ionizing source on one x boundary.
- Domain width: at least ten clump radii transverse to the illumination.
- Run both an aligned and a 22.5-degree rotated clump configuration to expose
  discrete-ordinate ray effects.

### B04: crossing beams

- Two equal isotropic HI sources, separated by 0.4 domain length.
- A pair of opaque clouds is placed so each source illuminates a distinct
  shadow boundary and the shadows overlap.
- The long-characteristic reference uses the same source spectrum and opacity
  tables. M1 is retained as an explicitly expected failure comparison.

### B05: dusty slab

- Uniform gas slab with a controlled UV optical depth of 1, 10, and 100.
- One-sided continuum source with fixed bolometric luminosity.
- Static gas for the first test; a separate prescribed-motion test checks the
  radiation force update.
- Absorbed UV energy, re-emitted IR energy, and integrated momentum are
  individually recorded.

### B06: X-ray cloud

- Neutral solar-metallicity cloud with hydrogen columns 1.0e20, 1.0e22, and
  1.0e24 cm^-2.
- Incident AGN power-law spectrum split over EUV, soft-X, and hard-X groups.
- Compare ionization and heating profiles with a high-resolution
  one-dimensional microphysics reference before enabling the 3D test.

### B07: clumpy dual-AGN nucleus

- Static 3D turbulent density cube with a controlled lognormal density PDF.
- Two AGN with adjustable separation, luminosity ratio, and obscuring column.
- Baseline separations: 100 pc, 500 pc, and 1 kpc.
- The suite reports source-resolved escape fractions and the force exerted by
  each AGN, not only their combined luminosity.

## 4. Standard reporting

Each run writes `summary.json` with git revision, TPU topology, precision,
mesh shape, quadrature, group table revision, source SED hash, wall time,
memory high-water mark, and all conservation errors. A benchmark result is
never accepted from a visualization alone.

## 5. Promotion sequence

1. B01 and B02 authorize H/He transport and chemistry development.
2. B03 and B04 authorize dual-source S_N development.
3. B05 and B06 authorize dust and X-ray source coupling.
4. B07 authorizes conversion of the first lagRamses zoom snapshot.

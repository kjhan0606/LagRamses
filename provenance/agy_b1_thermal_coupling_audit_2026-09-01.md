# AGY B1 thermal-coupling audit

Date: 2026-09-01  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: `/home/kjhan/local/bin/agy` (Gemini Antigravity CLI)  
Model: `gemini-3.1-pro-high`  
Scope: B1 thermal coupling only; read-only audit.

## Verdict

**BLOCKED**

## Confirmed-correct design elements

- Thermal and chemistry ledgers consistently use volumetric cgs rates and the
  heating-positive/cooling-negative sign convention.
- Collisional-ionization cooling and chemistry transitions use the same Abel
  coefficient implementation and the H I, He I, and He II ionization energies.
- The bounded chemistry and thermal solvers preserve positivity and photon
  accounting.
- Provenance checks reject the legacy equilibrium atlas and a UVB-bearing
  Cloudy source.

## Findings

1. **Critical:** `snrt_core/jax_thermal_atlas.py` interpolates rate values with
   weights from `log_metallicity_solar`. The source model is linear in `Z`, so
   this overestimates all off-grid metallicities; at the logarithmic midpoint
   between 1 and 10 solar it returns a 5.5 multiplier instead of 3.16.
2. **High:** `tests/b1_thermal_coupling.py` checks the linear-Z invariant only
   at exact metallicity grid points, which cannot expose the interpolation
   error.
3. **Medium:** the generator applies CMB subtraction only below the Grackle
   two-dex cutoff, leaving a finite step at the cutoff. AGY recommends
   continuous subtraction.

The first finding is an algorithm defect. The second is a validation gap. The
third reproduces the exact Grackle implementation at revision
`f93091ff8456962d7017a5bff7472945a30e3dad`, but is accepted here as a physical
continuity improvement that must be recorded as a deliberate deviation.

## Required closure

- Replace the tabulated metallicity dimension with one solar-metallicity table
  and apply metallicity as an analytic runtime multiplier.
- Add off-grid scalar and field-valued metallicity tests.
- Apply the CMB coefficient subtraction continuously and document the exact
  Grackle-reference deviation.
- Regenerate and rehash the atlas, rerun B1/P4/P5/P8 and 32-cubed convergence,
  then request a B1-only re-audit.

## Deferred scope

Non-equilibrium metal ions, dust IR re-emission, and live hydro fallback physics
remain later-gate work. This B1 result does not decide overall production or
publication readiness. Data-license approval remains independently pending.

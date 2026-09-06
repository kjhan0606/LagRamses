# Fable pre-bundle audit: dust scattering/state bundle (2026-09-06)

## Verdict: CONDITIONAL APPROVE

Fable judged this bundle necessary and on the critical path for G4 of the
production/publication readiness plan.  The current builder drops the Draine
table's scattering information and the current momentum diagnostic is
absorption-only.  The proposed separation of static scattering from grain
temperature, IR re-emission, live RAMSES force injection, and AMR/MPI/restart
qualification is physically and operationally appropriate.

The exact local isotropic-scattering solution proposed in the plan is
physically defensible for constant cell coefficients: the angular mean loses
only absorption, the directional field remains non-negative, and the operator
can conserve photon number.  It must be placed in the common transport path
used by both P4 and P5.

## Required amendments (applied to the plan)

1. Record scattering-weighted `⟨cosθ⟩`, `⟨cos²θ⟩`, and diagnostic
   `(1−⟨cosθ⟩) C_sca` per group.  Name the initial closure
   `phase_isotropic_candidate` and report the anisotropy bias; do not silently
   assume Henyey-Greenstein or delta-Eddington.
2. Treat `K_abs × M_dust/H` as the authoritative absorption value because the
   printed Draine columns are rounded.  Record the independent
   `C_ext/H × (1−albedo)` comparison with an explicit tolerance rather than
   claiming exact equality.
3. Remove `native` from this bundle's exit evidence.  Native Fortran dust
   channels and native-vs-JAX parity remain deferred and must be named as such.
4. Put the scattering operator in the common transport/multiphysics path and
   keep `absorbed_intensity` absorption-only so existing H/He and P5 primary
   absorption ledgers remain valid.
5. Add fail-closed runner semantics: scattering is off by default;
   scattering-enabled sidecars cannot be run with scattering disabled, and
   absorption-only sidecars cannot be used to enable scattering.
6. Add a scattering-weighted photon energy to the sidecar and use it for the
   scattering momentum ledger.
7. Represent dust-state origin explicitly (`direct` or
   `metallicity_solar_times_dust_to_metal`) and check consistency when both
   fields exist; do not invent an unapproved dust-to-metal prescription.
8. Add a redistribution benchmark, not only conservation: a beam must
   isotropize under optically thick pure scattering, and a uniform absorbing
   medium must reach the source/absorption steady state independently of
   scattering opacity.
9. Use two separately source-bound stellar/AGN synthetic fixtures with common
   edges; defer an aggregate STAR+AGN dust sidecar.

## Optional guidance

- Store group-averaged extinction as a sanity diagnostic.
- Use numerically stable `expm1`/small-optical-depth branches.
- Exercise at least S4 and S8.
- Extend existing dust tests rather than proliferating test files.
- Update P4/P5 and multiphysics documentation with the new scope.

## Audit basis

The audit read the proposed plan, `provenance/production_publication_readiness_plan.md`,
`simulation/snrt/snrt_core/dust.py`, `multiphysics.py`, `transport.py`, the
Draine builder and staged table, and the P4/P5 contracts.  No web search,
code edit, build, or job execution was used by the auditor.


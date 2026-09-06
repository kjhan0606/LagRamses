# Dust scattering/state bundle plan (2026-09-06)

## Purpose and project boundary

This bundle advances the project's high-level radiation-transfer, stellar/AGN
feedback, and dust objective toward a production/publication-ready closure.
It is a static SNRT physics bundle first; it does not claim that the live
RAMSES hydro coupling is already complete.  The implementation must remain
usable for both stellar and AGN source ledgers and must preserve the existing
H/He chemistry contract.

The current dust path is an audited absorption-only candidate.  The staged
Draine/Weingartner-Draine Milky-Way (R_V=3.1) table contains extinction,
albedo, and the first two scattering-angle moments, but the loader discards
the scattering information.  This bundle closes that omission without
silently promoting the candidate to an astrophysical production choice.

## Bundle scope

### D1. Explicit dust state and opacity contract

1. Extend the dust model and source-bound sidecar with non-negative
   scattering cross section per H and the declared angular closure metadata.
2. Treat the table-declared `K_abs × M_dust/H` as authoritative for the
   absorption cross section because the printed columns are rounded.  Record
   the independent check against `C_ext/H × (1-albedo)` with a declared
   tolerance, and independently record `C_sca/H=C_ext/H × albedo`.
   Recompute both from the raw table and reject inconsistent or out-of-range
   albedo/moment data; the known `9e-5` moment inequality envelope is treated
   as rounded source-column uncertainty and recorded.
3. Compute and store scattering-weighted `⟨cosθ⟩` and `⟨cos²θ⟩`, together
   with the diagnostic transport coefficient `(1−⟨cosθ⟩)C_sca`.  The first
   closure is explicitly `phase_isotropic_candidate`; it must report the
   bound implied by ignoring the measured anisotropy and may not silently
   become a Henyey-Greenstein or delta-Eddington model.
4. Make the cell scaling explicit.  The reference-mixture abundance is not
   inferred in the kernel: a named input contract records whether it is a
   direct abundance or a `metallicity_solar × dust_to_metal` product, together
   with its floor/normalization.  Missing or ambiguous non-zero dust state
   remains a hard error.
5. Keep v1 absorption-only sidecars as labeled reference controls.  A
   scattering-enabled sidecar gets a new version/status and cannot be accepted
   by an absorption-only caller by accident.

### D2. Conservative static RT scattering

1. Add a named, fixed angular phase-function closure.  The first candidate
   candidate is isotropic scattering unless the raw-data moments and a
   quadrature-level phase closure are explicitly validated; no unrecorded
   Henyey-Greenstein or delta-Eddington assumption is allowed.
2. Apply scattering within each photon group only: no frequency redistribution
   and no absorption masquerading as scattering.  Use the exact local
   constant-coefficient isotropic solve after the existing explicit spatial
   transport step, so the local angular mean loses only the absorption part.
3. Return separate directional absorption and scattering diagnostics.  The
   dust heating ledger counts absorption only; the dust momentum ledger
   includes incoming-minus-outgoing photon momentum for absorption plus the
   selected scattering closure, using physical (c) for momentum while
   transport may use reduced (c).
4. Record scattering-weighted photon energy separately from the
   absorption-weighted energy and use the matching value in the scattering
   momentum ledger.
5. Enforce positivity, finite values, angular-weight normalization, and
   photon-number closure in the local operator.  H/He ionization and
   photoelectron ledgers continue to see only gas absorption; dust absorption,
   scattering, and residual photons remain disjoint ledger channels.

### D3. Evidence and promotion fences

The bundle is accepted only with compact tests covering:

- raw Draine parsing, absorption/scattering reconstruction, bounds, and
  source-bound payload/hash validation;
- zero-scattering reduction to the existing absorption-only operator;
- pure-scattering photon conservation and directional momentum transfer;
- mixed absorption+scattering energy/photon closure and positivity;
- forward/reverse symmetric quadrature cases and at least two S_N orders;
- P4/P5 runner output fields and negative paths for an old/incompatible
  sidecar;
- two separately source-bound stellar and AGN synthetic fixtures with
  identical group edges; an aggregate STAR+AGN dust sidecar is deferred.
- explicit runner flag semantics: scattering is `off` by default; enabling
  isotropic scattering with an absorption-only sidecar, or disabling it for a
  scattering sidecar, is a hard error.

The evidence must report tolerances, dtype, group-edge hash, source/table
hashes, and the exact code hashes.  It must not claim live/restart/MPI hydro
qualification.  The staged physical table remains `candidate` until the
mixture, dust-to-metal/depletion prescription, source obscuration model, and
publication-facing citations are approved.

## Explicitly deferred to the next dust bundle

The following are not to be approximated in D1--D3:

1. temperature-dependent grain charging/size evolution and IR emissivity;
2. thermal balance with dust-gas collisional exchange;
3. absorbed UV/optical energy re-emission into configured IR groups;
4. full radiation-pressure force injection into live RAMSES hydro;
5. native Fortran dust-channel implementation and native-vs-JAX parity;
6. aggregate STAR+AGN dust closure, AMR/MPI/restart qualification, and a
   production cosmological rerun.

Those items form the next dust/feedback bundle after this static channel
closes.  In particular, a physical IR re-emission model requires an explicit
grain-temperature/emissivity prescription and cannot be created by merely
relabeling absorbed energy as an IR photon source.

## Entry and exit criteria

**Entry:** the AGN partition reference bundle has its Opus repair record and
consolidated native/production gate PASS, which is now recorded separately.

**Exit:** D1--D3 implementation, compact JAX/runner evidence, and
source/hash manifests pass; the result is marked `conditional_candidate` until
the deferred astrophysical selections are approved.  A failed phase-function
or source-state assumption blocks scattering promotion rather than silently
falling back to absorption-only behavior.

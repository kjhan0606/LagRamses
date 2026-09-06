# Independent stage-2 audit request: FS2010 secondary ionization

You are the independent scientific and algorithmic auditor for stage 2 of the
SNRT production/publication-readiness roadmap. Work read-only in
`/gpfs/kjhan/LRD_JWST`. Audit physical correctness, algorithmic justification,
wiring, validation design, provenance, and claim scope. This is not primarily a
style, security, or generic bug audit.

The required verdict is one of `PASS`, `CONDITIONAL PASS`, or `BLOCK`. Return
`PASS` only if there is no remaining stage-2 blocker. List findings in severity
order with exact file/line evidence, explain the physical or numerical impact,
and state a concrete remedy for every non-PASS finding. Explicitly distinguish
stage-2 blockers from work validly deferred to later roadmap gates.

## Stage-2 acceptance contract

1. Replace the old discontinuous 100 eV Shull--van Steenberg branch with the
   Furlanetto--Stoever (2010) electronic tables and continuous interpolation.
2. Demonstrate 99.9/100.1 eV continuity, table limits, independent interpolation
   agreement, finite/nonnegative fractions, and local photoelectron-energy
   closure.
3. Carry heating, excitation, and H I, He I, and He II secondary ionizations
   through both multiphysics and conservative-primordial chemistry, species
   ledgers, diagnostic names, runners, and sharding contracts.
4. Pin exact table provenance, license, upstream commit, and per-file hashes;
   make B2 and P5 canonical artifacts fail closed if code or tables change.
5. Remeasure the P5 temperature and volume-mean xHII effect in a matched OFF/ON
   control that passes fixed-point, H/He ledger, thermal, and finite-state gates.
   This is an effect/wiring control, not a spatial-convergence promotion.
6. Any production-reachable numerical issue uncovered while doing this work
   must be either repaired with a justified test or declared a blocker.

## Files and artifacts to inspect

- `simulation/snrt/snrt_core/secondary.py`
- `simulation/snrt/data/furlanetto_stoever_2010/README.md`
- `simulation/snrt/data/furlanetto_stoever_2010/TABLE_MANIFEST.json`
- all 14 `simulation/snrt/data/furlanetto_stoever_2010/*.dat` tables and the
  vendored MIT license
- `simulation/snrt/snrt_core/multiphysics.py`
- `simulation/snrt/snrt_core/conservative_primordial.py`
- `simulation/snrt/snrt_core/conservative_hydrogen.py`
- `simulation/snrt/snrt_core/implicit.py`
- `simulation/snrt/snrt_core/thermochemistry.py`
- `simulation/snrt/tools/p4_run_transport_pilot.py`
- `simulation/snrt/tools/p5_run_thermochemical_pilot.py`
- `simulation/snrt/tools/validate_multiphysics_b2.py`
- `simulation/snrt/tools/validate_p5_secondary_ionization.py`
- `simulation/snrt/tests/secondary_furlanetto_stoever.py`
- `simulation/snrt/tests/secondary_furlanetto_stoever_artifact.py`
- `simulation/snrt/tests/coupled_photo_collisional_hhe.py`
- `simulation/snrt/tests/b2_multiphysics_artifact.py`
- `simulation/snrt/tests/p5_secondary_ionization_artifact.py`
- `simulation/snrt/data/furlanetto_stoever_validation.json`
- `simulation/snrt/data/b2_multiphysics_transport_validation.json`
- `simulation/snrt/data/p5_secondary_ionization_validation.json`
- `simulation/snrt/SECONDARY_IONIZATION_VALIDATION.md`
- `simulation/snrt/P2_MULTIPHYSICS.md`
- `simulation/snrt/B2_PRODUCTION_SOLVER_VALIDATION.md`
- `simulation/snrt/P5_THERMOCHEMISTRY_VALIDATION.md`
- the current git diff from HEAD `438aff4`

The two P5 HDF5 files named in the canonical P5 JSON are present on GPFS and
their SHA256 values are checked by the artifact test.

## Questions requiring explicit judgment

- Do the table columns and the implementation's use of `f_ion`, `f_heat`,
  `f_exc`, and species ionization counts preserve the intended FS2010 physics?
  Is the small all-channel renormalization scientifically defensible?
- Does bilinear interpolation and low/high boundary handling make the former
  100 eV discontinuity genuinely absent without introducing a new relevant
  discontinuity?
- Is H II, rather than an electron-fraction surrogate, the correct coordinate,
  and is the primordial-composition approximation disclosed with an honest
  claim boundary?
- Are H I, He I, and He II secondary counts charged exactly once to chemistry
  and to the photoelectron-energy ledger in both solver families?
- Is the nearest-root electron-density bisection in `implicit.py` a justified
  way to avoid the spurious high-ionization backward-Euler branch for a hot,
  exactly neutral, unilluminated cell? Inspect bracketing, multiple roots,
  vectorization, precision, and ledger consistency. Does it create a new
  production blocker or merely a measured performance cost?
- Is changing the B2 200 eV fixture from invalid H-only composition to
  `n_He/n_H=0.079` scientifically correct, and are the new regression bands
  supported rather than used to hide a failure?
- Does the 0.1 Myr P5 matched pair actually support the reported very small
  mean deltas, and are its limitations stated clearly enough for publication?
- Do canonical reports bind every production dependency needed to reject stale
  results? Are the source/license/upstream claims adequate?

Known later-gate items are not by themselves stage-2 blockers unless this
implementation makes their current state unsafe: RSLA/refinement is stage 3,
AGN nine-group recovery is stage 4, P5 spatial/timestep science convergence is
stage 5, full coupled H+He front validation is stage 6, dust closure is stage 7,
and the broader thermal stress matrix/performance optimization are later gates.

Conclude with a checklist mapping each of the six stage-2 acceptance items to
`closed`, `conditional`, or `open`, then give the single overall verdict.

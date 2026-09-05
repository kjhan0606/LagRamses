# AGN event-level entrainment — implementation evidence

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Baseline: `e427029`; status: implemented, native checks pass, Opus end review PASS.
Authority: operator approved the amended plan with "진행하시오.", including
donor-frame resolved loading, independently normalized lobes, missing-lobe
fallback, and the explicitly retained cold-loading donor-heating convention.
Plan: `agn_entrainment_bundle_plan_2026-09-05.md`.

## What changed in the actual native path

- `agn_feedback_deposition.f90` now supplies one cylinder membership/weight
  helper for both passes, a half-open donor containment test, two-lobe deltas,
  cold donor withdrawal, staged hydro-plus-composition deposition, a small
  scalar-index validator and pure donor-payload pack/unpack routines.
  The existing five-field cap/deferred-energy arithmetic is unchanged.
- `AGN_feedback` selects the declared gas layout without initializing any
  stellar source: legacy `nelt` fields at `ichem`, or all eleven stored
  channel-resolved fields. A stellar source being disabled does not erase the
  gas field from advection. The total metal slot is additional, not counted
  twice as an element. NVAR/overlap checks include reserved `idelay`,
  `ivirial`, `ixion`, `isgs` where enabled. `NENER/=0` is rejected before any
  AGN mutation. Current builds have `NENER=0`; delayed cooling is not disabled.
- `average_AGN` sums the two volume-weighted lobe normalizations separately.
  It locates donors with half-open cells, records potential duplicate hits
  under a short OpenMP critical region, and globally checks exactly one
  owner per active directed event before mass removal. This is source
  ownership protection, not a new domain-decomposition algorithm.
- Geometry/owner consensus precedes sequential donor loading. Each load uses
  the then-current donor state, including the 25% cap and composition. It
  removes bulk kinetic energy, not donor internal-energy density. The actual
  loaded mass, three donor-velocity components and all transported fractions
  are packed by the owner and distributed using the existing nsink-indexed
  MPI_SUM pattern. Nonowners contribute zeros. All new error exits follow
  MPI_MAX consensus; ranks with nAGN=0 still participate.
- `AGN_blast` uses that globally available donor velocity, never BH bulk
  velocity, for the resolved return. Each supported lobe receives half the
  loaded mass. The second kinetic moment is retained where both streams
  overlap. The old nominal EAGN, kinetic fraction and volume-weighted thermal
  fraction remain in use, with no new efficiency or BH recoil prescription.
- If either global lobe sum is zero, all distributed deposition for that
  event is skipped and the owner executes the whole donor-cell fallback.
  The fallback no longer requires vol_gas=0. It returns mass, momentum,
  kinetic energy and composition there and thermalizes EAGN. Shared-donor
  fallbacks are serial. Hydro and composition commit together per receiver.
- `bin/Makefile` declares the conditional stellar configuration dependency
  introduced by the layout choice. No Python file, new runner or general
  payload/geometry framework was added. Existing unrelated dirty changes
  remain untouched; no commit/push is part of this implementation request.

## Proportional verification performed

The existing `simulation/snrt/tests/run_fp15_agn_efficiency.sh` compiles the
actual production AGN helpers and `agn_feedback_deposition_smoke.f90` at
O0 with Intel run-time checks. It passed:

```text
SNRT_AGN_EFFICIENCY_OK
SNRT_AGN_SOURCE_OK
AGN_NATIVE_CELL_COUPLING_SMOKE_OK
SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK
```

The same depositor and smoke also passed at O3 in
`build/agn_entrainment/agn_entrainment_O3_smoke`. Native cases include:

- Existing below/crossing/already-hot cap and nonunit-volume saved energy.
- Unequal-volume/unequal-weight lobes with a midplane receiver and mixed
  kinetic/thermal injection; summed donor+receivers conserve mass and all
  three momentum components, and gas-energy change plus deferred energy
  equals EAGN. A uniform velocity boost and a cap-crossing event also pass.
- Missing-lobe fallback leaves distributed receivers unchanged and restores
  donor composition; no half event is injected.
- Donor abundances remain unchanged on removal; receiver abundance moves
  toward the donor mixture. Global species sums alone would miss the old bug.
- Two sequential withdrawals from a shared donor use the updated density,
  loading cap and composition. The unchanged delayed-cooling density is
  explicitly checked. The 4/3 donor specific-heat rise is tested.
- Valid legacy, channel-resolved and empty layouts; invalid extent/overlap;
  invalid composition leaves the hydro and scalar row unchanged; exact donor
  face ownership is unique. The prior Opus source-error-code nit is included.
- One synthetic owner payload plus one nonowner zero payload exercise actual
  pack -> sum -> unpack, matching the owner-plus-zeros exchange convention.
  This is a serial payload arithmetic test, NOT execution of MPI exchange or
  proof of full distributed evolution. MPI collective placement is inspected
  in the caller and the caller is compiled.

Small synthetic case absolute tolerances are 1e-11 for mass/species, 1e-10 for
momentum and 1e-9 for gas-plus-deferred energy under the boost. These are not
large-run production tolerances or a new convergence standard.

Changed `sink_particle.kjhan.f90` compiled successfully in project scratch:

1. SNRT + PHASE0_STELLAR_ENRICHMENT.
2. CPU + PHASE0_STELLAR_ENRICHMENT.
3. CPU without PHASE0_STELLAR_ENRICHMENT (legacy layout branch).

Compiler: `mpiifx -qopenmp -fpp -O0 -g`, `NDIM=3`, `NVECTOR=500`, `NPRE=8`,
`NVAR=18`, `NENER=0`, `SOLVER=hydro`, `LONGINT`, `QUADHILBERT`,
`OUTPUT_PARTICLE_POTENTIAL`, `USE_FFTW`; freshly built AGN modules before
existing `bin` RAMSES interfaces in the include path. An initial scratch
compile encountered an old `bin/snrt_agn_efficiency.mod`; rebuilding that
module in scratch resolved it without modifying shared build products.

`git diff --check` passed. No full RAMSES link, native MPI run, snapshot,
live evolution, long calculation or shared binary overwrite was performed.

## Limits and retained model choices

- Donor bulk frame and separately normalized lobes change the model relative
  to historical runs. The geometry axis remains the BH/spin axis. Removing
  this event's numerical recoil is not a physical BH recoil model.
- Cold loading leaves donor internal-energy density in place and raises its
  specific heat by up to 4/3 per withdrawal, outside the receiver temperature
  cap. Several withdrawals can compound that factor. No enthalpy advection,
  upstream positivity remedy or new thermal-mass partition was approved.
- Delayed-cooling, virial, ionization and SGS reservoirs remain untouched.
  Their specific values can therefore change with gas density; their mass
  loading/thermodynamic semantics are not qualified. No dust hydro field is
  declared by the current field map, and no dust prescription is introduced.
- The tests establish the focused mechanical event identities on a fixed
  receiver stencil, not every possible interacting source event in live
  RAMSES. In particular, general overlapping mass-weighted thermal/replay
  sources while other jets change receiver density are not qualified here.
- MPI ownership/exchange wiring is compiled and statically checked, not a
  restart, load-balance, AMR/MPI scaling or real-evolution claim. The previous
  live-SNRT+AGN serial fresh-start restriction and legacy+SNRT exclusion stay
  unchanged. No SED, spin/radiative/mechanical partition or AGB source choice
  is altered. This is not full production/publication qualification.

One Opus 5 bundle-end review covered the complete implementation, evidence,
approved scope and honest limits. No per-helper reviews were added.

While the end reviewer was running, the driver independently noticed that
the boosted event's original cap=4 did not actually defer energy. Only the
existing boosted/fallback smoke cases were tightened to cap=1 with an
explicit deferred>0 assertion. No production source was changed during that
review. The final verification reruns cover this strengthened test; the
reviewer may have read the preceding test version, which is not claimed as
independent review of the added assertion.

## End review and closure disposition

Full result: `opus5_agn_entrainment_bundle_end_audit_2026-09-05.md`.
Opus returned **PASS**, with no required in-scope repair. It independently
confirmed the mechanical-event algebra, native wiring, distinction between
stellar source activation and gas-field transport, collective ordering and
the stated limits. The final strengthened O0/check-all and O3 smoke both pass,
including explicitly positive deferred energy. No follow-up audit is needed.
The review's phrase "exact by construction" describes the real-arithmetic
identity; actual floating-point evidence is tolerance-based, not bitwise exact.

Nonblocking follow-ups retained without starting a new gate:

1. The unchanged finite-radius spatial-bin search may miss a coarse donor's
   cell centre even when the sink is physically inside that cell. The new
   owner-count check then refuses loading instead of silently skipping it.
   Live stencil qualification and more informative sink-specific diagnostics
   remain future AGN work; donor ownership over all AMR configurations is
   not claimed by a half-open containment test or the new consensus check.
2. `read_hydro_params` currently places `isgs=ichem+1`; with SGS enabled and
   two or more chemical fields this overlaps the chemistry block and the
   new layout guard rejects it. This is an explicit newly rejected combination,
   not a claim that SGS+chemical entrainment is supported. Ordinary disjoint
   delayed-cooling plus chemistry remains allowed.
3. The active build requires OpenMP. The pre-existing dead non-OpenMP branch
   of `AGN_blast` is not qualified or repaired here.
4. The automatic donor send/receive buffers scale as
   `2*(4+nscalar)*nsink` reals (25.6 MB at nscalar=12, nsink=100000 for 8-byte
   reals). Large-sink memory/stack qualification is deferred, not a blocker
   for this bounded native implementation.
5. The owner-plus-zero wording has been clarified above. This remains a
   serial payload check, not two live MPI ranks.

No production source was changed in response to these nonblocking notes.
This authorized implementation bundle is complete within its scope.
No commit, push, runtime activation or next implementation bundle was started.

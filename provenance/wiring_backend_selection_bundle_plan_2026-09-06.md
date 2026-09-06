# Cross-subsystem wiring and backend-selection bundle plan

Date: 2026-09-06
Repository: `kjhan0606/LagRamses`
Working tree: `/gpfs/kjhan/LRD_JWST`
Base commit: `85cfa45`

## Objective

Revalidate the production wiring of SNRT, stellar/AGN feedback, and dust, then
decide whether an OpenMP/CUDA automatic selector is scientifically justified.
The project objective remains a production/publication-ready high-level
hydrodynamic implementation of radiative transfer, stellar/AGN feedback, and
dust. A backend switch is admissible only when both backends implement the
same state transition, source/receiver contract, precision policy, AMR/MPI
ownership, and conservation ledger.

This is one consolidated bundle. It does not add per-stage telemetry or a new
live RAMSES evolution run.

## Wiring review

### SNRT

The expected transaction is:

```text
AGN source budget
  -> snrt_agn_source / photon groups
  -> snrt_ramses_driver
  -> AMR topology + prepared transport
  -> CUDA multigroup species/dust ABI
  -> raw/H-He/dust/returned ledgers
  -> FP64 ledger receiver
  -> H/He thermochemistry
  -> transaction commit or rollback
```

Evidence must verify that the native caller reaches
`snrt_transport_absorb_multigroup_prepared_dust_trial`, that the direct H/He
ledger reaches thermochemistry, and that the old host-side
`snrt_partition_absorption` path is not used by this caller. The existing
three-species CUDA ABI must remain linked for compatibility. `ZERO_SCAFFOLD`
must remain explicit: dust opacity is zero and dust has no persistent
thermal/momentum/abundance receiver in this bundle.

### Stellar and AGN feedback

The expected stellar path is:

```text
thermal_feedback / sub_thermal_feedback
  -> phase0_feedback
  -> stellar_ramses_runtime
  -> stellar source/yield/population ledgers
  -> stellar_ramses_bridge
  -> transactional uold/unew deposition
```

The expected AGN path is:

```text
amr_step -> AGN_feedback
  -> average_AGN (withdrawal/geometry/ownership)
  -> AGN_blast (thermal or jet receiver)
  -> agn_feedback_deposition primitives
  -> RAMSES hydro/scalar state
```

SNRT AGN photon production is a separate, explicitly selected source path;
legacy mechanical AGN feedback must not be silently combined with live SNRT
AGN source ownership.

### Dust

The current dust wiring is intentionally partial:

```text
SNRT prepared transport -> CUDA H/He+dust ledger -> FP64 ledger validator
                                      -> (H/He only) thermochemistry
                                      -> no dust state commit

snrt_dust_coupling / snrt_dust_ir -> native candidate APIs and smoke tests only
                                  -> not called by the live RAMSES driver
```

The review must preserve this distinction and must not describe candidate IR
or nonzero dust opacity as live dust feedback.

## Backend decision and bounded implementation

The review will record one of `AUTO_EQUIVALENT`, `EXPLICIT_ONLY`, or
`UNAVAILABLE` for each subsystem, using source-level evidence and the existing
native/production bundle gate.

| Subsystem | CUDA path | OpenMP/CPU path | Decision |
|---|---|---|---|
| SNRT | CUDA multigroup transport and DUST-7 ABI | OpenMP host orchestration only; no equivalent transport kernel | `EXPLICIT_ONLY`: CUDA required; no silent CPU fallback |
| Stellar feedback | no equivalent CUDA stellar source/deposition path | FP64 source/ledger/bridge with OpenMP locks/parallel caller | `EXPLICIT_ONLY`: CPU/OpenMP path; do not map unrelated GPU hydro flags to it |
| AGN feedback | `gpu_sink` is a legacy control/autotune surface, but no equivalent consumer is present in the inspected AGN routines | native AGN feedback routines and SNRT reference path are CPU/FP64 | `EXPLICIT_ONLY`: retain named model selection; no fake GPU switch |
| Dust | CUDA primary transport boundary; FP64 host ledger receiver | no equivalent OpenMP transport implementation; native IR is a non-live candidate | `EXPLICIT_ONLY`: retain the hybrid contract; no cross-backend auto switch |

The bounded implementation is therefore a backend capability/wiring record and
negative checks against false equivalence. The nonfunctional `gpu_sink`
benchmark/switch block is removed from `amr_step.jaehyun.f90`; the namelist
name remains as a compatibility field but is explicitly nonselecting. A
runtime auto-selector is not added in this bundle: it would either select a
nonexistent implementation or silently change the physics/precision contract.
A future selector may be introduced only after an equivalent CPU SNRT kernel,
CUDA stellar/AGN source path, and OpenMP dust transport path each pass
differential conservation tests.

## Evidence and acceptance gates

1. Source-level wiring map with positive and negative edge checks.
2. Existing consolidated native/production gate passes, including GNU/Intel
   receiver and A10 CUDA checks.
3. Explicit proof that `gpu_sink` does not claim to select the new AGN/stellar
   feedback path; if the dead control is retained, it is recorded as a
   separate legacy issue rather than promoted as auto-selection.
4. No live RAMSES evolution, nonzero dust activation, IR receiver, or restart
   claim is made by this bundle.
5. Codex end audit records the final disposition and deferred prerequisites.

## Deferred work

An actual automatic selector remains a later implementation bundle. Required
inputs are: a reference CPU SNRT transport with identical AMR/MPI semantics; a
GPU implementation of the same stellar/AGN feedback state transition; a CPU
OpenMP dust transport equivalent; deterministic cross-backend tolerances; and
runtime device/affinity handling. Until then, CUDA absence for SNRT remains a
fail-closed condition and no GPU flag is allowed to imply new feedback or dust
physics.

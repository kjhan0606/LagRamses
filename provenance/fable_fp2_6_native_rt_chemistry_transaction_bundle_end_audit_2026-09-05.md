# Fable fallback bundle-end audit — F-P2.6 native RT/chemistry transaction

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Fable
Mode: read-only source, build, and evidence audit

## Verdict

`FAIL`

The bundle had the right transaction/fixed-point architecture, but the
implementation/evidence combination was not yet production-ready at the
time of this audit. The failures below were treated as remediation blockers,
not as reasons to broaden the bundle into physical SED, yield, dust, or live
hydro validation.

## Severity-ranked findings

### B1 — coarse trial-buffer capacity was underspecified

The driver allocated the coarse correction buffer with active-slot capacity,
while the transaction commit requires the complete persistent intensity
capacity. A correction targeting a valid persistent slot could therefore be
rejected or lost when those capacities differ.

Reference: `patch/lagRamses/snrt_ramses_driver.f90` and
`patch/lagRamses/snrt_rt_transaction.f90` coarse-buffer allocation/shape
checks.

Disposition: fixed in the remediation by allocating the driver buffer with
`size(snrt_intensity,3)` and retaining an exact full-capacity commit check.

### H1 — zero-leaf ranks needed allocated payloads

An empty rank could have only a slot map allocated while photon and chemistry
arrays were absent. That is unsafe when the rank must enter the same MPI
transaction and topology collectives as ranks with leaves.

Reference: `patch/lagRamses/snrt_state.f90` initialization and growth path.

Disposition: fixed by allocating zero-length photon and H/He payloads and
covering begin/restore with the MPI zero-leaf smoke.

### H2 — persistent chemistry arrays are slot-indexed

The first implementation mixed leaf-position and persistent-slot indexing in
transaction begin/restore. With non-contiguous slots, rollback could restore
the wrong H/He state or H I mirror.

Reference: `patch/lagRamses/snrt_rt_transaction.f90` begin/restore loops.

Disposition: fixed by resolving `islot=leaf_slot(i)` for every persistent
chemistry and mirror access; the named non-contiguous-slot rollback smoke
passes.

### H3 — rank-local return hazards existed around collective transport

Unmapped faces, optional topology arguments, interface exchange, coarse
correction, and CUDA failure could cause one rank to return before peers
reached the next collective. The CUDA-availability test also needed to be a
collective preflight.

Reference: `patch/lagRamses/snrt_transport_step.f90`,
`patch/lagRamses/snrt_amr_topology.f90`, and
`patch/lagRamses/snrt_ramses_driver.f90`.

Disposition: fixed for the production prepared path with collective error
reductions before each return, collective topology requirements, zero-leaf
CUDA bypass, and driver `uold`/unit/CUDA preflights. Legacy scalar adapters
remain compatibility/diagnostic surfaces and are not on the production
driver path.

### H4 — evidence did not exercise all named failure routes

The initial evidence proved core rollback and convergence primitives but did
not show all partition/chemistry/receiver selectors, mixed inventory cap
behavior, or a driver-route reachability check.

Disposition: fixed at the evidence level by running the named stage/leaf
selectors for all three stages, adding a mixed H/He inventory CUDA case, and
adding a static check for every driver selector. A live initialized RAMSES
failure-injection run remains explicitly record-only because no hydro
evolution was authorized in this bundle.

### H5 — diagnostic failure injection needed an explicit mode

The normal production configuration correctly rejects test injection, but a
driver-reachable diagnostic path was not explicit enough for audit use.

Disposition: fixed with `SNRT_RT_TX_DIAGNOSTIC_MODE=1`; the driver labels this
mode as non-production, permits the named injection only there, and retains
the production rejection test.

### M1 — pre-partition scale was required for host tolerance

Using the residual working inventory as the scale after earlier groups had
consumed atoms could make the tolerance depend on group order and become too
small for the original cell inventory.

Disposition: fixed by passing the cell's pre-partition inventory scale to
`snrt_partition_absorption` and covering it in the native thermochemistry
smoke.

### M2 — CUDA residual guard was fixed-absolute

The old fixed remainder threshold was not scale invariant for FP32 inventory
updates.

Disposition: fixed with the documented `256 * FLT_EPSILON` scale-aware guard
based on absorbed and available cell magnitudes, plus a mixed-inventory CUDA
smoke.

### M3 — no-absorption cells should not fail temperature reconstruction

A positive-density cell with no absorbed photons does not enter this bundle's
chemistry source. Rejecting it solely because its hydro temperature could not
be reconstructed would create an unrelated false failure.

Disposition: fixed by skipping chemistry/temperature validation for cells
with no absorbed group signal; any actually absorbed cell still requires a
finite positive pre-heating temperature.

### M4 — trial-state finite/bound checks were incomplete

The prepared CUDA output and coarse correction must be checked before any
chemistry or commit path can consume them.

Disposition: fixed by rejecting non-finite/negative trial photon and
absorption arrays and non-finite coarse corrections, followed by finite,
non-negative absorption/heating and H/He simplex checks before convergence.

### M5 — reproducibility evidence was stale

The evidence did not yet contain hashes for the changed CUDA surface, the
actual current smoke output, GNU compile details, or an explicit iteration
cost bound. The relaxed-iterate residual semantics were also insufficiently
stated.

Disposition: remediation requires a complete evidence refresh: source and
binary hashes, current CUDA output, GNU/mpiifx/full-link records, worst-case
fixed-point cost (`up to 32` transport evaluations per prepared level, each
with its configured subcycles), and an explicit relaxed-iterate residual
definition.

### L1–L4 — build/config/ledger hardening

The transaction and transport module prerequisites needed `mpi_mod.o`; the
transaction configuration report should be emitted once; the environment
buffer should not truncate diagnostics; and the unassigned ledger should be
globally reduced.

Disposition: fixed in `bin/Makefile` and
`patch/lagRamses/snrt_ramses_driver.f90`; the refreshed build and MPI smoke
are required evidence.

### L5–L7 — non-blocking semantics and future scope

The audit requested explicit relaxed-residual semantics, documentation that
legacy scalar adapters are compatibility-only, and accounting for coarse
corrections with no currently mapped persistent slot.

Disposition: relaxed residual and adapter semantics are recorded in the
bundle documentation; no-slot coarse accounting remains a G5 AMR topology
scaling/coverage item and is not a blocker for this bounded local bundle.

## Acceptance-gate disposition at initial audit

| Gate | Initial disposition | Remediation |
|---|---|---|
| C1 transaction/rollback | FAIL | Slot/zero-leaf/coarse/collective fixes |
| C2 unassigned receiver gate | CONDITIONAL | Scale-aware cap and refreshed ledger evidence |
| C3 bounded fixed point | CONDITIONAL | Relaxed-residual semantics and refreshed smoke evidence |
| C4 native build/evidence | FAIL | Current CUDA/GNU/MPI/build/hash transcript |

## Scope boundary

This failure did not evaluate or reject physical AGN/stellar SEDs, yield
tables, DTD/PISN physics, momentum/thermal subgrid feedback, dust
scattering/IR/grain evolution, HDF5 restart integration, live RAMSES
evolution, distributed AMR scaling, or publication convergence. Those remain
separate high-level physics or G5/G6 bundles.

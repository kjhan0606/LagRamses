# Astra review: native QUMOND AMR implementation

Read-only Astra review, 7 November 2026. This is not a Fable verdict, a solver
implementation, or production approval. The active RAR goal remains matched
Lambda-CDM with a CDM halo versus full-field baryon-only QUMOND, then BIMOND/FDM,
with stellar/AGN feedback, metals and cooling/heating; dust/RT are excluded.

Q-GOAL: aligned. Stage-specific N/Q histories are necessary prerequisites but
cannot repair omitted phantom-source cells or one-way elliptic coupling alone.
Q-LEAN: implement the existing composite design in three bounded bundles;
additional proxy variants or repeated first continuous-source tests are not
needed. Existing continuous-source/timestep trials must remain acknowledged.

## Approved next implementation sequence

1. Integrate explicit current/prior N and Q potential state, level/physical-time/
   hierarchy-generation validity and same-stage boundary selection. Preserve
   ordinary phi/phi_old and the current cold-restart correction. Seed each
   stage's first prior from that stage's first solved field. Exact history
   restart needs versioned arrays/timestamps; legacy restart is a declared cold
   reconstruction, not exact continuation. Apply periodic gauge consistency.
   Lifecycle includes allocation growth, ghost exchange, child initialization,
   relocation, load balance and restart; unsupported paths remain rejected.
2. Implement a native leaf/subface operator and synchronous two-stage solve.
   Freeze actual deposited gas-plus-particle sources at one hierarchy epoch.
   Enumerate each same-level face and coarse/fine subface once. Reconstruct
   coarse potential tangentially to the fine subface and divide the N-potential
   difference by (H+h)/2 for the normal gradient. Reconstruct tangential
   components at that face before evaluating nu(|g_N|/a0). Use identical
   area-weighted flux with opposite signs in adjacent leaf volumes. Covered
   coarse cells are restriction/workspace, not extra physical volumes.
   The same gradient/divergence pair defines BOTH composite Poisson residuals.
   Native make_fine_bc_rhs uses coarse interpolated Dirichlet data and returns
   no interface flux to coarse rows: unchanged sequential MG is insufficient.
   Use level MG/FFTW as a preconditioner with outer defect correction/FGMRES;
   do not choose CG without establishing symmetry/positive definiteness.
   Publish all levels' forces only after both composite solves converge.
3. Extend that same operation to deterministic MPI subface ownership, topology
   rebuild after regrid/load balance, history checkpointing and coupled
   evolution. Start evolving activation with a common hierarchy timestep;
   nsubcycle=1 alone does not establish synchronized physical clocks. Time
   subcycling requires explicit time-centered interface flux semantics.

First native activation may be periodic, noncosmological, serial, fixed two-level
2:1 hierarchy, with unsupported modes rejected. This bounds first implementation;
it does not replace the eventual MPI/evolving/restart objective. Do not replace
native AMR by a global auxiliary uniform grid as the final solution.

## Source findings and validation scope

Current 76c3976 poisson_commons has shared phi/phi_old, not separate N/Q histories.
interpol_phi uses both shared fields. Current amr_step seeds cold restart history
but subsequently overwrites phi in its second solve. rho_fine uses phi and flags
as scratch; calling it early is not read-only hierarchy source preparation.
Coarse kicks and accretion precede fine recursion, so a late fine correction
cannot repair already consumed forces. Existing AMR measurement skips 17304 of
187072 fine source cells. The separate history trials do not remove this layer.

One validation bundle per implementation bundle: state/startup/lifecycle and
unchanged default path; manufactured plus compact-interface source at two
resolutions with both-stage composite residuals and pre-mean-removal shared flux
cancellation; then bounded evolving coupled MPI/regrid/restart with current
exchange inventories. Check particle-force consistency before long evolution.
These are implementation acceptance criteria, not new proxy research projects.

Boundary/external-field prescription, force softening/coupling, and a specified
MOND background/ICs remain scientific decisions for the corresponding science
runs. They do not prevent the bounded native implementation. No new solver
accuracy or viable relativistic MOND cosmology is established by this review.

Source isolation: codex/rar-qumond-amr-history starts from 76c3976 in
/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107. Existing source,
queued jobs, snapshots and prior trial guards remain unchanged.

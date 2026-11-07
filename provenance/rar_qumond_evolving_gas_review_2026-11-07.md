# Astra review: evolving gas with native composite gravity

Read-only Astra follow-up at source13d3b7f, under the source-tree substantive
coupling review requirement. This is not a Fable verdict. The full RAR goal
remains Lambda-CDM/QUMOND then BIMOND/FDM, and hydro with stars/AGN/metals/
cooling/heating, without dust or ray tracing.

Q-GOAL: approve direct evolution of actual gas density. Q-LEAN: retain the
existing Godunov gravity predictor and old-density half-source; replace gravity
scheduling, rather than adding another explicit initial half-kick.

The compiled driver is patch/cuRamses/adaptive_loop.jaehyun.f90, not the generic
amr/adaptive_loop.f90. Configure the evolving flag before runtime refinement,
then suppress topology changes/remap/defrag while preserving initial mesh,
lower-level restriction and ghost bookkeeping. Keep exact hierarchy checks.
Use serial periodic CPU, noncosmological gas, two fixed levels and common dt
for first activation; particles/subgrid/MPI/regrid/restart remain guarded.

At the old common epoch publish native g_n from actual leaf density before
any CFL/flux consumer. Skip legacy rho_fine in this path; it clears phi and
performs unrelated deposition. Populate rho_max from physical, not phantom-added,
density. Call newdt_fine once per level, then select a common minimum including
any conservative acceleration bound. It updates dtold, so do not call again
inside fine recursion. Preserve parent set_unew before fine descent and all
native flux/reverse exchanges. Do not replace old covered uold before parent
Godunov commits. Keep ctoprim g_n*dt/2 predictor and set_uold hidden source
0.5*dt*rho_n*g_n. That source is HALF a timestep, correcting imprecise previous
wording describing a full gravity source routine. Remove legacy pre-solve
old/new sync and per-level final old-force half-kick for this mode.

After the full hierarchy transports and commits, finish covered restriction,
validate raw positive finite leaf density, solve from new rho at t_{n+1},
publish all new forces, and apply one leaf-only synchydrofine1(+dt/2). Restrict
all conserved fields again and refresh hydro/gravity ghosts. Momentum source
is dt/2*(rho_n*g_n+rho_{n+1}*g_{n+1}). Each source operator preserves its
intermediate internal energy by adding the corresponding kinetic-energy
change. This is not a proof of total gas-plus-QUMOND field energy conservation
or exact discrete vanishing self-force; ordinary econs cannot be used as one.

Keep the single deepest update_time advance. Record t_start, assert final
t=t_start+dt, and rebase completed endpoint cache across the outer coarse-step
counter increment only after validating physical time, gas density and hierarchy.
Do not rotate previous histories twice for the same physical source time.

One bounded validation bundle: smooth nonuniform moving gas across refinement
interfaces, matched Newtonian/QUMOND and timestep-halved QUMOND at the same
endpoint. Require density change, independently reconstructed new-source
residuals, synchronized clocks, fixed topology, positive finite states, mass
closure across interface fluxes and actual momentum/kinetic source increments
matching the specified half-step formulas. Recheck default-off/frozen behavior,
without rerunning the disc survey. Releasing additional source modules later
requires explicit placement of their mass/energy updates relative to epochs.

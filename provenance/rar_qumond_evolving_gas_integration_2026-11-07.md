# Native common-time evolving gas implementation

This is the reviewed direct gas-coupling step, not another static field proxy.
Source e61e7f8 adds explicit RAR_QUMOND_COMPOSITE_EVOLVE, requiring native
publication/history, nonstatic serial CPU periodic gas, fixed two levels,
common timestep, no cosmology/particles/stars/sinks/cooling/chemistry/RT/SGS/
pressure-fix/eEOS/restart/remapping. Dust/RT remain off in the isolated build.
The original default and frozen-source modes retain their existing operation.
These restrictions describe first activation; eventual full source modules,
MPI/regrid/restart and declared cosmological IC/model are still required.

The compiled outer cuRamses driver configures the flag after initial mesh
construction and before the main loop. Runtime refine/flag/remap/defrag paths
are bypassed, while lower-level restriction and ghost bookkeeping continue.
Every native history entry checks exact stored grid/cell/son identity.
Legacy rho_fine is skipped only in evolving composite mode. At coarse entry,
covered gas and ghosts are refreshed, actual physical gas density supplies the
native source, both gravity stages publish before CFL/flux consumers, and
rho_max uses physical source density rather than phantom density. newdt_fine
runs exactly once per active level. All dtnew entries receive the same minimum
of these CFL/free-fall candidates, an explicit conservative acceleration bound,
and optional positive RAR_QUMOND_EVOLVE_DT cap.

The old complete hierarchy force is held through recursive Godunov transport.
Parent unew is initialized before fine descent and never reset after fine
reverse fluxes. Existing gravity prediction and set_uold's old-density HALF
source are retained. No explicit initial kick or legacy pre-solve correction
is applied. After all levels transport/commit, covered gas is restricted, raw
leaf states are validated before any density floor, both native stages solve
from the new density, and only physical leaves receive the new-force half
source through synchydrofine1. Conserved fields are then restricted again and
all gas/gravity ghosts refreshed. The single deepest update_time remains;
end-phase assertions require one common clock advance and unchanged dt.

The complete endpoint field is reused next step only after checking its time,
exact gas source and hierarchy. The coarse counter may be explicitly rebased,
without rotating prior N/Q physical-time histories again. Optional exports
record each newly solved source epoch, not repeated cache reuse. Optional
per-leaf old/new half-source records retain density source/receiver, force,
halfdt and actual before/after momentum/energy. This checks the defined kinetic
work exchange, not a conserved gas-plus-QUMOND field energy or zero self-force.

Initial build1275956 failed because the hydro Courant parameter was omitted
from the module import. Correction b45ff50 compiled in build1275970 (3m14s).
Its first execution1275971 failed before initialization: the compiled cuRamses
reader does not accept the generic hydro condinit_kind selector. Correction
8297489 introduces hydro_initial_condition in the actual lagRamses parameter/
reader and cuRamses initializer, and removes the unused generic initializer
change. The explicit rar_periodic_wave IC is rho=1+0.1cos(2pi*x)+0.05sin(2pi*y),
with region-specified constant pressure and velocity. Default selector region
preserves existing behavior. Failed source/build/run records are retained.
Only our superseded dependency-waiting wrappers1275969/1275978 were cancelled;
no running science job or unrelated queue was cancelled.

Build1275982 uses clean numerical source8297489. Dependent bounded execution
1275983 has Newtonian/QUMOND six steps at dt cap2e-4, and QUMOND twelve steps at
cap1e-4, targeting common endpoint0.0012. Two fixed levels4/5 and spherical
refinement radius0.25 place moving smooth gas across the interface. Each case
expects two snapshots at midpoint/endpoint (foutput3/3/6), with noutput1,
aout1.1,tout1e100,fbackup1e9. Six dumps have a100663296-byte upper estimate;
512MiB includes epoch/kick exports and1GiB additional free-space margin is
required. Explicit effective namelist paths and source/binary/output hashes
are retained in /gpfs/kjhan/RAR_QUMOND/runs/rar_evolving_gas_JOB/CASE.

Dependent audit1275984 checks all completed physical source epochs, actual gas
mass, source clocks, topology, half-source densities/forces/work increments,
matched initial states and common-time half-dt replay. Six saved snapshots
receive independent native source/Poisson/phantom/force/publication and positive
state reconstruction checks. Execution1275983 completed in15s and audit1275984 in8s. Its full artifact manifest verifies, as do all three execution input/output manifests; each source identity is8297489. Earlier b45ff50 regression1275979 and audits1275985/1275986 completed. Updated-source regression1275999 completed in18s and audits1276001/1276002 in4s/6s; their artifact manifests are checked separately below.


## Measured evolving-gas acceptance

All three cases contain5384 physical leaves on an unchanged native hierarchy.
Newtonian and QUMOND6-step runs and the12-half-step QUMOND replay reach0.0012,
with source/field and saved times matching to floating-point roundoff. Density
relative changes are1.2632548187546246e-4 (N),1.2515126068576821e-4 (Q),
1.2514204616608022e-4 (Q half). Maximum relative gas-mass residual is
2.220446049250313e-16 across all exported source epochs. The maximum per-leaf
momentum-source residual is5.739156437416237e-17 and energy/kinetic-work
residual6.017677811132396e-16 code units. Old-half output and new-half input
states join exactly: no unrecorded intermediate gas change is detected.

Halving dt changes final density/momentum/energy by relative L2
7.825616039111481e-8 /1.2999454123135646e-7 /2.3793576668430788e-8. This is a
same-endpoint replay, not a measured temporal order or long-term accuracy proof.
Newtonian versus QUMOND final density/momentum/energy differences are
3.997569344021054e-6 /0.004394320052256419 /0.00012020791654201173. These are
physical-model responses within this short smooth fixture, not an RAR result.

All six saved midpoint/endpoint checks independently reconstruct the native
operator and source from actual moving-gas snapshots, require positive finite
gas/internal energy, match published leaf fields and covered restrictions, and
pass both true residual bounds1e-9. Example QUMOND endpoint N/Q residuals are
9.068301637482354e-10 /9.740708267650551e-10; phantom/force relative differences
6.308104057421056e-16 /6.749515193363363e-17; saved density/source difference
exactly zero. Snapshot gas momentum is nonzero as intended; unchanged gas energy
is not asserted. The source-state verifier distinguishes fixed and evolving
modes explicitly. Compact results are recorded in
rar_qumond_evolving_gas_result_1275984.json; full provenance stays under
/gpfs/kjhan/RAR_QUMOND/audits/rar_evolving_gas_1275984.

This closes the first actual evolving-gas coupling within the reviewed fixed
serial hierarchy. It is not implementation of particle force coupling, stellar/
AGN source placement, cooling/metals, MPI, regridding, restart or cosmological
background/ICs. Those remain required next work. The ordinary econs diagnostic
omits native QUMOND field energy and is not used as total-energy evidence.

Updated-source default/frozen audit manifests1276001/1276002 verify. All sixteen
old/new-default-off gravity arrays are bitwise equal. The frozen-source native
publication verifier also passes. No disc survey rerun was needed for this
coupling change; numerical implementation and evolving results are evaluated
within their stated activation, and existing default/frozen behavior remains.

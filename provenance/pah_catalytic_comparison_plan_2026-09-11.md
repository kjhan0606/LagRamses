# PAH H2 catalytic comparison within approved remaining bundle2

Final goal: simulation-ready RT/feedback/dust. User explicitly authorizes
named approximation comparisons, keeping existing defaults unchanged.
This is one implementation within bundle2, not completion of carbon loss,
higher charges or metallic-Fe hard-photon charging.

## Physical and native scope

Add selector `pah_h2_catalytic_v1`, extending the existing H0--13, neutral/
cation M13/DL01 comparison and its optional H2 vacancy capture. Add only
the barrierless, one-site C24H13 + H -> C24H12 + H2 abstraction channel for
both existing charges. Boschman et al.2015 section2.1.1, equation3 and
TableA.1 give cross section0.06 Angstrom^2, thermal H speed and no barrier
for H13, with one reactive site. Primary source:
https://doi.org/10.1051/0004-6361/201323165
https://pure.rug.nl/ws/portalfiles/portal/72850158/aa23165_13.pdf

Retain M13 attachment rates and H13 binding3.2eV, not claim to reproduce
Boschman's full H0--36 network. Use the existing H2 ground binding4.4781eV.
The reaction releases1.2781eV per event under this hybrid thermochemistry.
Explicit immediate-accommodation comparison: preserve PAH excitation,
return ground-state H2 and put net bond release into gas thermal energy.
Incoming H kinetic energy stays in the gas/H2 bath; no additional kinetic
energy is created or removed. This is not a resolved H2 rovibrational or
outgoing translational spectrum. Charge is unchanged.

Use a single finite H donor for both charges and all H13 excitation bins.
Their common coefficient permits an exact bimolecular event count, with
population transfers proportional to initial H13 populations. Preserve
both reactants, avoid clipping an unfunded event, and commit all outputs
together after H nuclei, molecule number, charge and total energy checks.
The existing split photo/IR/H step handles attachment and abstraction as
explicitly ordered local operators. Do not claim exact simultaneous-cycle
kinetics; check timestep refinement for the combined cycle.

## Wiring and minimal evidence

Reuse the existing PAH mass carriers, CHIMES atomic/molecular H exchange,
gas particle heat capacity and restart identity. No new passive slots,
datasets, table framework or dispatch. Update both setup frontends and
their common tests; old selectors and their identity stay unchanged.
Native checks: exact finite donor, zero donor/zero dt, per-charge molecule
and H balance, excitation preservation, bond/gas/H2 energy conservation,
rollback, and integrated old/new selector behavior. Follow with the
existing bounded PAH MPI2/OMP2 live/restart fixture, not a new simulation
campaign. Clean raw only after evaluation is complete.

Selective plan review asks Q-GOAL FIRST and Q-LEAN SECOND. Its scope is
the new reaction/energy coupling, not reopening completed bundle1.

## Selective review disposition

[Fable returned CONDITIONAL](pah_catalytic_fable_2026-09-11.txt).
All three required design changes are adopted: catalytic is a strict
superset of the existing H2-capture selector; bind the Maxwell mean speed
constants and coefficient with a new identity tag only for that selector;
run the common-H13 abstraction after both charge loops and before
recombination. Structural excitation/per-charge tests belong in the native
smoke, not extra runtime gates. Gas accommodation is the maximum immediate
gas-heating allocation of the net reaction release; other allocations
would be separate choices, not inferred measurements. Its suggestion that
operator-order bias is small is not itself a measured timestep result.

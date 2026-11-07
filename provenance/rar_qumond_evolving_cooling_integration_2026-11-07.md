# Original cooling/heating coupled to evolving native QUMOND gas

This implements the reviewed direct thermal coupling at source c5eaa66, followed
by numerical accuracy control82bf272. The reviewed placement retains full-dt
level-local Lie-split cooling after transport/old-force half source and before
endpoint native publication/new-force half source. The cooling operator does
not change density or momentum in the admitted path; no additional gravity
solve is inserted solely for thermal exchange. Old covered coarse states remain
unchanged until parent transport commits. This is not a Strang scheme or a
claim of second-order full physics evolution.

Cooling is admitted only with evolving native publication, original non-Grackle
thermal solver, isothermal/delayed/neq/RT/Jeans/eEOS off. Frozen publication still
rejects cooling; PIC/star/sink/SGS/MPI/cosmology/regrid/restart exclusions remain.
The disabled polytrope now sets T2min exactly zero, avoiding zero times an
undefined nISM normalization. No cooling/emissivity/UV table rates were altered.

Each actual level call is bracketed by whole-owned-leaf inventories V2 with
step, physical time and level. The inventory is now safe with pic=false: it
neither sizes nor traverses unallocated particle arrays. Existing PIC inventories
retain V1 and identical accumulation. Evolving per-cell cooling receipts V2 add
step/time identity, avoiding constant-aexp aliasing and the old32768-cell
assumption. Independent analysis checks leaf membership/volumes per call,
unchanged mass/momentum/passives, accepted versus complete-call energy, and
intervening thermal energy between the two gravitational half sources. Initial
temperature clipping is reconstructed separately from the effective table RHS.

## Physics and output contract

The moving periodic wave retains identical density/velocity across controls,
with initially advected metal fraction0.002 (0.1solar for the code's .02
normalization). This is metal-dependent cooling, not enrichment. Cooling uses
haardt_madau=true with the code's default Courty UV parameterization, fixed
z0/aexp1 in this noncosmological path, J21=0, self-shielding off. It is not an
HM2012 table claim. Dust/RT remain off. Initial pressure1 selects cooling;
pressure0.001 exercises heating. T2max1e10 and table floor T/mu1e-2K are explicit;
no nonzero clipping correction is activated in these measurements. Units are
rho1e-24g/cm3, length3.0856775814913673e22cm and time3.8707679655184575e15s.

Final suite1276113 has N/Q six common steps at2e-4, Q twelve at1e-4, Q cooling-off,
Q heating, and Q tightened thermal substeps, all to0.0012 on5384 fixed leaves
(levels4/5). It is one bounded integration/precision bundle, not galaxy
formation or an N-versus-Q observational model comparison. Each case requests
two midpoint/endpoint dumps: noutput1,aout1.1,tout1e100,foutput3 (half-dt6),
fbackup1e9. Twelve dumps have a201326592-byte upper estimate;1GiB reserves all
snapshots, source epochs, kick records and logs, with another1GiB free-space
margin required. Effective input paths and immutable source/binary hashes are
retained in /gpfs/kjhan/RAR_QUMOND/runs/rar_evolving_gas_1276113/CASE.

## Accuracy finding and response

Build1276056 and initial original-varmax4 execution1276088/audit1276091 passed
exchange/field checks, but half-dt final gas-energy difference was2.2246%.
This was not accepted as thermal accuracy. Source82bf272 exposes the existing
solver's inverse fractional-temperature-step limiter as cooling_varmax>=4,
default4 unchanged. Higher values refine the same original integration;
the500-iteration rejection remains. New build1276097 completed in2m54s.

Execution1276099 uses128 with a256 precision replay. Audit1276105 correctly
FAILED the declared0.001 relative-energy criterion: half-dt difference
4.2831733346522685e-5 but substep-refinement difference0.001324128070541183.
The criterion was not weakened, nor were rates fitted to an outcome. The failed
audit and original execution remain preserved. Final execution1276113 uses256
and512, the same immutable binary, and completed in38s. Audit1276133 completed
in11s and passed. Half-dt relative gas-state differences are density
6.022013242885519e-7, momentum1.6647170834604253e-4 and energy
1.3093953048330516e-5. Thermal substep256-to512 differences are density
5.410231905702462e-9, momentum2.5143487681900446e-6 and energy
6.573909783876405e-4. Both energy comparisons meet the declared0.001 criterion
for this hot cooling fixture only, not arbitrary heating or galaxy conditions.

## Measured source and field checks

The Q cooling run has32304 cooling receipts, no heating receipts, and accepted
energy-0.5349510550383955 code units. Q heating has32304 positive receipts,
accepted+0.003803439522427016. Both levels have complete owned-leaf call coverage.
All cooling-call mass/momentum/passive integrals are unchanged. Across cases,
maximum whole-call energy/accepted-receipt difference is1.1987806580737725e-14
code energy; RHS/clipping reconstruction residuals are at most7.77e-16 per cell.
Clipping energy is exactly zero. Gravitational half-source energy residuals
remain<=5.80e-16, with maximum gas-mass relative residual2.22e-16. The Q
cooling-on versus off final total gas-energy field differs by relative L2
0.3701810476593755, so thermal coupling is physically active, not a no-op.

All twelve saved midpoint/endpoint source/field checks pass independently from
actual gas snapshots. Example Q endpoint N/Q residuals9.105835072148621e-10 /
9.81557218496351e-10; independent phantom/force differences6.245547113523431e-16 /
6.73983089860072e-17; source mismatch exactly zero. All states are positive and
finite. Full audit and all six execution input/output manifests verify. Compact
result: rar_qumond_evolving_cooling_result_1276133.json; full immutable evidence:
/gpfs/kjhan/RAR_QUMOND/audits/rar_evolving_gas_1276133. Archived audit scripts are
authoritative; later tracked tool wording clarifies scope without changing its
numerical calculation. Original result JSON is preserved unchanged.

## Regression and remaining scope

Latest default/frozen run1276114 and audits1276134/1276135 pass, including all
sixteen bitwise-equal default-off gravity arrays; manifests verify. A separate
legacy PIC/star/AGN/metals/cooling execution1276100 completed in23s and audit
1276106 in7s: V1 AGN inventories, uniform whole-cooling budget and all-element
stellar transactions pass, with verified artifacts. This checks preservation of
that existing path, not native coupling of particles or feedback. Default
cooling_varmax4 is preserved in this regression. No disc survey was rerun.

The full programme still needs native particle deposition/force/time coupling,
stellar/AGN/enrichment placement, MPI/regrid/restart and declared cosmological
background/ICs, then BIMOND and FDM. These results establish original thermal
coupling within serial fixed-hierarchy gas evolution; they do not certify a
calibrated cooling curve, population RAR, feedback, total gas-plus-QUMOND field
energy conservation or full scientific completion.

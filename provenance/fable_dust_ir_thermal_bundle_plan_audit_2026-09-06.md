# Fable pre-bundle audit: DUST-2 grain thermal balance and IR re-emission

- Date: 2026-09-06
- Workspace: `/gpfs/kjhan/LRD_JWST`
- Requested model: Claude CLI `--model fable`
- Reported model: `claude-fable-5-1`
- Mode: read-only plan audit; no edits, jobs, builds, or web access
- Plan audited: `dust_ir_thermal_bundle_plan_2026-09-06.md`
- Final objective supplied to auditor: production/publication-ready
  high-level RT coupled to stellar/AGN feedback and dust in lagRamses

## Verdict

**CONDITIONAL APPROVE.**

The bundle is necessary, physically defensible in the existing static JAX/SNRT
architecture, and appropriately lean. It is the next missing dust energy
channel after DUST-1. Implementation may begin after the amendments below are
incorporated into the plan. The result must remain a staged one-pass ledger
candidate, not a claim of recursively transported IR closure or astrophysical
dust approval.

## Findings requiring amendment

### F1 — blocking: permit and record out-of-band emission

The lowest configured edge is 0.01 eV (124 micrometres). Real cold dust emits
substantially below that edge; for a modified blackbody with beta approximately
2, the fraction below 0.01 eV is roughly 0.8 at 15 K, 0.5 at 20 K, 0.2 at
30 K, and 0.04 at 50 K. Requiring configured IR fractions to sum to one at
every temperature would reject or silently renormalize realistic tables.

Amend the contract so configured-group fractions plus a tabulated
`untracked_energy_fraction` sum to one. Record the untracked energy as a named
ledger and close `absorbed = tracked IR + untracked`. Do not treat it as an
unreported loss.

### F2 — blocking: add temperature-dependent emission photon energies

Power and energy fractions do not define a photon rate. The thermal sidecar
must carry per-group emission-weighted mean photon energy versus temperature,
in eV, and each value must lie within its group interval. The photon rate is
`group fraction * emitted power / group mean photon energy`; the source-ledger
mean energy must not be reused. This is especially important for the very
broad 0.01--1 eV group.

### F3 — blocking: include a CMB floor rather than aborting cold cells

Near-zero local absorption occurs in ordinary cells, and a hard cold-table
error would abort realistic runs. The equilibrium input must include the CMB
background term evaluated from the same admitted absorption/emissivity data at
the run epoch. The table must extend to or below the background temperature.
The cold boundary then becomes a load-time coverage check. The first candidate
does not need stochastic heating, but must label the single-temperature
assumption.

### F4 — required: derive thermal data explicitly from admitted absorption data

The prohibition on reusing the Draine table as emissivity data is too strong.
Kirchhoff's law permits deriving the equilibrium power curve, group fractions,
and emission-weighted photon energies from the admitted absorption cross
section times the Planck function, under an explicitly declared
single-temperature equilibrium assumption. The builder, quadrature tolerance,
temperature grid, and derivation convention must be hashed. Small-grain/PAH
stochastic emission remains deferred and must be stated as a limitation.

### F5 — required: bind thermal and opacity mixtures

The thermal sidecar must carry the active absorption source-table hash and dust
mass per H nucleus. The runner must compare these with the active opacity
sidecar and fail closed on mismatch. The historical v1 sidecar may not be
silently accepted if it lacks the source-table identity; restrict DUST-2 to
pinned v3 controls or read the raw sidecar field explicitly.

### F6 — required: place the operator inside the thermochemical subcycle

The current instantaneous dust heating rate and the cumulative heating ledger
are different quantities. Solving temperature from the final rate and copying
the cumulative energy would make closure tautological. Apply the thermal
operator per subcycle in the existing thermochemical loop, accumulating
emitted, tracked, untracked, and photon ledgers with the same dt weights as
the heating ledger. Test table-evaluated emitted power at the solved
temperature against absorbed power.

### F7 — required: return masks from JIT, raise in the runner

JAX-jitted code cannot perform the requested Python hard error. The operator
must return an out-of-range/coverage mask and the runner must raise after the
step. Do not use silent edge clamping for a production candidate.

### F8 — required honesty diagnostic for one-pass IR

One-pass IR recording is a valid staged contract, but IR self-absorption can be
order unity in compact dusty cells. Record only two scalar diagnostics already
derivable from current fields: maximum cell IR optical depth over configured
IR groups and the fraction of dust-bearing cells with optical depth above one.
Do not add a per-cell optical-depth history.

### F9 — minor numerical amendments

Use a strictly increasing power curve, define the zero-input temperature
sentinel as exactly 0 K, exclude sentinel cells from finite/monotonicity tests,
and mask zero dust abundance before division. The hot branch remains a masked
out-of-range result handled by the runner.

### F10 — optional simplification

Because the admitted power curve is strictly monotone, inverse interpolation
in log-log space may replace 32 bisection iterations. Bisection remains
acceptable if fixed-shape control flow is retained, but it is not required for
this bundle.

## Scope and instrumentation assessment

The plan is not over-instrumented. One loader, one fixed-shape thermal
operator, one opt-in P5 path, and the compact tests are proportionate. The
second source-bound stellar/AGN subprocess should be removed: thermal metadata
binds to the group/mixture contract, not to a source SED, and DUST-1 already
tests source identity. Extend the existing `tests/p5_dust_runner.py` instead.

Do not add per-cell optical-depth arrays, per-subcycle histories, separate
silicate/carbon temperatures, a stochastic mid-IR ledger, a second synthetic
schema, a full cosmological run, or a new source ledger.

## Required plan amendments before implementation

1. Replace the zero-untracked rule with tracked plus tabulated untracked
   energy closure and an explicit untracked ledger.
2. Add temperature-dependent emission-weighted mean photon energies.
3. Derive the thermal table explicitly by Kirchhoff's law from the admitted
   absorption table, with single-temperature equilibrium and hashed builder
   provenance.
4. Add source-table hash and dust mass per H binding against the active
   opacity sidecar.
5. Add the CMB term and require table coverage down to the run's background
   temperature; retain a masked hot out-of-range result.
6. Run thermal balance per thermochemical subcycle with dt-weighted ledgers.
7. Return JIT masks; raise only in the runner.
8. Record maximum IR optical depth and above-unity-cell fraction as scalar
   honesty diagnostics.
9. Tighten monotonicity/sentinel handling and replace the separate stellar/AGN
   thermal subprocess with an extension of the existing P5 runner test.

The report found no need for a new architecture. It judged the one-pass
recording useful only as a clearly labelled staged closure because it does not
yet alter transported radiation.

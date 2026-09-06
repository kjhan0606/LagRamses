# Fable plan audit — F-P1H-F physical-source admission and high-mass seam

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Audited plan: `fp1h_f_physical_source_admission_bundle_plan_2026-09-05.md`

## Verdict

**CONDITIONAL APPROVE**

Approve the intent, reject the original shape. The bundle is not
overinstrumented in audit cadence, but it would build a new admission profile,
a new converter path, and a new closure gate for a source package that does
not exist yet. The reduced plan is a legitimate next step. F-P2.7 is not
closed: its D4 initialized-RAMSES smoke is still pending in Slurm, its
source/fixture renames and bundle gate are uncommitted, and its evidence/end
audit are pending.

## Findings

- No operator population decision is recorded for IMF, single/binary model,
  metallicity domain, rotation/engine marginalization, chemistry scope, or
  Pop-III/PISN treatment.
- No authorized production package exists. The physical-node inventory and
  production selection remain empty; staged candidates are review inputs.
- The 40--120 M☉ interval is already fail-closed in the configuration,
  resolver, and production runtime. W18/N20 mixed outcomes correctly reject a
  universal direct-collapse rule. Nothing in the next bundle should weaken
  this boundary.
- Existing machinery already supplies most proposed Work Packages A and B:
  locked rights profiles, admission contracts, candidate adapters, canonical
  converter guards, and the G2 preflight. A new profile/report/sidecar would
  be a fourth identity layer.
- The existing F-P1 population-fate runner, G2 preflight, and F-P2.7 bundle
  gate already provide the native/regression surface. A new native closure
  runner would duplicate it.
- PPISN/PISN eligibility and ownership are a separate unresolved physics
  decision. They must not be designed as part of this admission-only step.

## Required reduction

### Entry gate — no implementation

Before code changes, record all three inputs:

1. F-P2.7 evidence with D4 pass or an actionable scheduler blocker;
2. an operator population decision naming IMF, single/binary treatment,
   metallicity domain, rotation/engine marginalization, chemistry scope, and
   Pop-III/PISN exclusion;
3. an authorized package naming exact source files, checksums, and rights.

If the second or third input is absent, the only permissible artifact is the
corresponding decision/authorization record. Do not build infrastructure that
can only reproduce the already-known blocked result.

### A' — reuse the existing admission machinery

Fill the existing admission-contract fields for the selected domain and
population, and extend the existing locked candidate profile for the
authorized package. Do not add a new profile schema, report schema, or
sidecar. The existing admission audit JSON remains the single machine-readable
report.

### B' — replace only satisfiable blocked adapters

For an authorized package, replace blocked adapters one-for-one while keeping
the current report contract, and add only the source-specific row emitter
needed by the existing canonical converter. Keep 40--120 M☉ blocked unless the
package covers every required axis and field. Emit zero rows and null selection
for missing physical or rights evidence. Do not edit source archives in place
or infer lifetime, wind, energy, momentum, decay, fallback, or remnant values.

### C' — reuse current regressions

Add no native runner. Reuse the F-P1 population-fate runner, G2 preflight, and,
after F-P2.7 closes, the bundle gate. Capture one implementation evidence
record. Individual assertions remain evidence components, not new project
gates.

## Deferred or removed

- Defer PPISN/PISN eligibility and channel ownership.
- Defer a Fortran source-node consumer, native row-reproduction acceptance, and
  runtime restart/retry/MPI invariance for this bundle.
- Remove the new admission profile document, its sidecar, and the new native
  closure regression.
- Remove the acceptance clause requiring native and converter row agreement.

## Final-purpose and overinstrumentation assessment

The reduced A'/B' work directly targets the actual G2 blocker and therefore
contributes to the production/publication RT--stellar/AGN-feedback--dust goal.
The removed work mostly restated existing blocked machinery or anticipated a
Fortran consumer that is not yet present. This reduction is necessary to avoid
another audit-heavy, evidence-heavy cycle without a physical input.

## Authorization boundary

After the three entry artifacts exist and the operator approves the reduced
plan, repository-local edits may touch only the existing admission contract,
locked rights profile, satisfiable adapter/row-emitter path, and one evidence
record. No PPISN/PISN design, inferred 40--120 M☉ value, new native runner,
HDF5/JAX/environment movement, runtime activation, or hydro evolution is
authorized.

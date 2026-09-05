# F-P2.7 gate-consolidation and initialized-RAMSES integration bundle plan

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)
Parent: F-P2.6 native RT/chemistry transaction bundle
Work location: `/gpfs`

Status: Fable plan audit returned `CONDITIONAL APPROVE`; all seven conditions
are applied below. The operator requested implementation after this audit.
The bundle proceeds in the reduced mandatory order defined below. No push,
large RAMSES evolution, or destructive storage cleanup is included.

## Project goal and decision boundary

The project goal is production-ready and publication-ready high-level RAMSES
radiative transfer coupled to stellar/AGN feedback and dust. This bundle is an
engineering consolidation bundle, not a new physical model. It must reduce
verification friction while preserving the native fail-closed, conservation,
MPI/CUDA, physical-model approval, and live-RAMSES safety boundaries required
before feedback production.

It must not expand into new yield, SED, DTD/PISN, dust, momentum-feedback, or
publication-science work. It must not replace native evidence with Python-only
or source-grep evidence.

## Motivation from the Fable operational audit

Fable classified the current workspace `OVERINSTRUMENTED`. The audited tree
contains 19 native shell runners, 19 native smoke programs (11 called and 8
orphaned), six phase0 references (four mirror-building runners and two
phase0-only drivers), eight divergent mirror modules, and no initialized-
RAMSES SNRT smoke. F-P2.6 used five audit events although the intended policy
is one plan audit plus one end audit. The earlier claim that four runners each
performed a full SNRT CUDA build was false; the SNRT-enabled build is new
coverage, not duplicated coverage.

## Deliverables

### D1 — one provenance index and evidence taxonomy

- Add `provenance/INDEX.md` mapping each active high-level bundle to its plan,
  implementation evidence, plan audit, end audit, and explicit non-approval.
- Define three statuses for records: `active`, `superseded`, and
  `long_term_record`; do not rely on filename suffixes alone.
- Keep one durable audit report per audit event. Prompt consolidation and
  historical relabeling are deferred; scratch copies under `~/.claude/plans`
  are not project evidence.
- Do not rewrite scientific conclusions or remove physical approval records.

### D2 — one native bundle gate

- Add `simulation/snrt/tests/run_snrt_bundle_gate.sh` and a corresponding
  `snrt_bundle_gate` Makefile target.
- The gate must call the existing native thermochemistry, spectral/checkpoint,
  transaction/MPI, and CUDA controls, then perform the `mpiifx`/CUDA production
  build and native-symbol check once. This is net-new SNRT-enabled build
  coverage, not a claim that the current runners already duplicate such a
  build. It must produce one concise machine-readable/text evidence summary
  rather than four independent narrative gates.
- Existing focused runners remain callable for local debugging, but their
  individual PASS markers must not be counted as separate project gates.
- The gate must preserve conservation checks, zero-leaf MPI coverage,
  production failure-injection rejection, and `git diff --check`/hash capture.
- Missing `mpirun`, missing `nvcc`, or an MPI run with fewer than two ranks is a
  hard failure, never a skip/pass. The summary must name the spectral group-sum
  versus Lbol conservation assertion and the CUDA photon-budget assertion,
  list injection-rejection counts by family, and record elapsed seconds for
  every stage on the first run.
- Static driver-route grep checks must be labelled `STATIC_SUPPORTING_CHECK`,
  never as live execution evidence.

### D3 — production-source ownership for native stellar controls

- Establish `patch/lagRamses/` as the production source of truth.
- Audit the six phase0 references currently used by
  `simulation/snrt/native/phase0/` and choose one reproducible strategy:
  adapt the four mirror-building runners to compile production sources, and
  move the 11 phase0-only drivers to an explicitly named fixture namespace.
  Do not regenerate another mirror. Do not delete the mirror until every
  affected runner has a passing replacement and a source identity record.
- Add a source-identity/parity check that fails if a production-relevant test
  runs against a stale unrecorded copy. A driver that fails against production
  sources is a recorded finding to fix, not a reason to retain the mirror.

### D4 — small initialized-RAMSES SNRT smoke

- Build a bounded, new-directory smoke using the linked SNRT binary and an
  explicitly audited effective namelist: one initialized hydro level, a few
  cells, one or two short steps, two MPI ranks where available, outputs
  suppressed by an intentional future schedule, and no large production data.
- The effective namelist must include `&STELLAR_ENRICHMENT_PARAMS`, use
  `feedback_mode='legacy'`, set `hydro=.true.`, `poisson=.false.`,
  `pic=.false.`, `levelmin=levelmax=3`, `tout=1.0d30`, `aout=2.0`, and
  `nstepmax=2`. The run must use a Slurm GPU allocation on a new directory;
  never run the two-rank CUDA smoke on a login node.
- Verify startup spectral/thermochemistry/transaction contracts, finite
  hydro state, photon/species/thermal closure, same-final-trial commit, and
  clean-stop behavior for one diagnostic failure injection. The injected run
  must be clearly non-production and must not overwrite an existing output.
- Record the absolute namelist, compiler/binary identity, MPI size, output
  policy, wall time, and log markers. A successful process exit alone is not a
  pass; reject NaN/negative state, unexpected output, or missing fail-closed
  markers.
- First determine from the existing RAMSES namelist/module interfaces whether
  a minimal native harness can be run safely. If it cannot, record the exact
  missing initialization contract and do not substitute a Python/static test.

### D5 — safe artifact consolidation and checkpoint

- Add `*__genmod.f90` to `.gitignore`; leave source and scientific artifacts
  untouched unless explicitly classified.
- **Removed from this bundle:** moving superseded JSON/provenance variants.
  The index will classify them without path churn; any later archive move is a
  separate operator-approved storage task.
- Do not delete or move the 299 GB `.quarantine_hdf5` partial, the virtual
  environment, or large JAX validation trees in this bundle. Mark them for a
  later operator-approved storage decision.
- After an explicit scope manifest excludes generated products and unrelated
  files, create the local checkpoint commit **before** D2/D4 for the currently
  implemented F-P1/F-P2 source and evidence so hashes resolve. Do not push to
  GitHub in this bundle unless separately requested.

## Acceptance gates

- `provenance/INDEX.md` identifies every active bundle and no active evidence
  is made discoverable only by an auditor scratch path.
- One bundle gate runs the required native controls and full link exactly once
  and reports their results without weakening any individual conservation or
  fail-closed assertion.
- No production-relevant native test silently compiles a stale phase0 mirror;
  the chosen ownership strategy is reproducible and hashable.
- The small initialized-RAMSES smoke either passes with direct runtime markers
  or produces a documented, actionable initialization blocker. It is not
  replaced by static grep or Python output.
- The index-only classification is reference-safe and does not touch large
  external data. Generated products are ignored rather than committed.
- The bundle's implementation plus gate/audit time is demonstrably smaller
  than the current per-runner/per-audit process, while the project still has a
  clear pre-feedback safety boundary.

## Explicit exclusions

No new stellar/AGN physical source selection, 40–120 M_sun yield resolution,
SNIa/DTD/PISN decision, dust scattering/IR/grain evolution, momentum/thermal
subgrid design, HDF5 restart integration, production hydro evolution,
distributed-AMR performance campaign, or publication convergence study is
part of F-P2.7.

## Fable plan-audit conditions applied — 2026-09-05

Fable returned `CONDITIONAL APPROVE`. The factual premises were corrected:
19 runners; 19 smoke programs with 11 callers and 8 orphans; six phase0
references; eight divergent mirror modules; no existing SNRT-enabled automated
build; 222 auditor-prefixed provenance files; and five F-P2.6 audit events.
The mandatory order is now checkpoint commit, D2, D4, then D3. Artifact moves,
prompt consolidation, historical relabeling, and large storage cleanup were
removed or deferred. D2 now hard-fails missing MPI/CUDA/two-rank coverage and
records named conservation assertions, rejection counts, and timings. D4
includes the required stellar-enrichment namelist, legacy feedback mode,
GPU/Slurm, output-clock, and clean-stop constraints. D3 ports four runners to
production sources and gives phase0-only drivers an explicit fixture scope;
it does not create a new mirror.

The full plan audit is retained in the auditor scratch report
`/home/kjhan/.claude/plans/f-p2-7-gate-consolidation-and-glowing-sunrise.md`;
the durable project record will be added to `provenance/` before implementation
evidence is finalized.

## Audit request and workflow

Before implementation, Fable must issue one plan verdict with the following
specific questions:

1. Does each deliverable directly reduce a real barrier to production-ready
   RT/stellar/AGN feedback/dust, or is it merely bookkeeping?
2. Does D2 reduce duplicated instrumentation without removing independent
   native safety checks?
3. Is D3 technically feasible without creating another source mirror?
4. Is D4 a realistic small initialized-RAMSES gate, and is its output policy
   safe under the shared RAMSES run rules?
5. Is the scope appropriately sized for one bundle, or should deliverables be
   deferred/merged to avoid another over-instrumented cycle?

The plan auditor returned exactly one classification at the top:
`APPROVE`, `CONDITIONAL APPROVE`, or `REJECT`. It must identify mandatory,
deferrable, and redundant work. Implementation begins only after conditions
are applied and operator approval, which the operator supplied with the
instruction to audit and then implement. At bundle end, use the normal single
end audit; no audit is requested for individual file moves or tiny repairs.

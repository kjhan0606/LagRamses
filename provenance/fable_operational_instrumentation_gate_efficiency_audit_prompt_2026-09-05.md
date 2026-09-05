# Fable operational audit request — instrumentation and gate efficiency

You are performing a read-only operational audit of the current project
workspace `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, branch `main`). Do not
edit files, delete or move artifacts, launch jobs, run tests, invoke
Python/JAX, or use network tools. Inspect the actual worktree and distinguish
pre-existing dirty-tree work, generated build products, audit records, and
the currently active native RT/chemistry bundle.

The project goal is production-ready and publication-ready high-level RAMSES
radiative transfer with stellar/AGN feedback and dust physics. The operator's
concern is that progress has been slowed by overproduction of instrumentation
and verification artifacts and by too many gates/audits relative to
implementation. Assess that concern directly and quantitatively where the
worktree permits; do not reward documentation volume as evidence quality.

Read the project-local instructions and current planning/evidence material,
then inspect at least:

- `AGENTS.md` and relevant project handover/planning files;
- the current F-P2.6 plan, implementation evidence, and Fable/Opus audit
  records under `provenance/`;
- the preceding F-P2.x plans/evidence and their audit prompts/results;
- `simulation/snrt/tests/`, `simulation/snrt/tools/`, and native smoke
  runners, including whether checks are Python-only, source-static, native
  Fortran/CUDA, or live RAMSES;
- generated compiler modules, binaries, logs, and other untracked products;
- git status, file counts, sizes, timestamps, and duplicated or superseded
  evidence where available.

Answer the following operational questions:

1. Is the project currently over-instrumented, under-instrumented, or
   appropriately instrumented for its production/publication goal? Classify
   the state as exactly one of `OVERINSTRUMENTED`, `BALANCED`, or
   `UNDER-INSTRUMENTED` at the top of the report.
2. Quantify the principal sources of delay: implementation versus audit time
   if records support it; number and size of prompts/reports/evidence files;
   repeated tests or duplicated controls; generated artifacts; and any
   high-cost checks that do not materially increase confidence.
3. Separate evidence that is mandatory before a live feedback/hydro run from
   useful but deferrable evidence and from obsolete/redundant evidence. Pay
   special attention to the distinction between native code tests, static
   route checks, Python bookkeeping validators, full production builds, and
   live initialized-RAMSES tests.
4. Evaluate whether the current audit governance is excessive. The intended
   policy is one plan audit before a bundle, one end-of-bundle audit after the
   bundle, with a primary auditor and a fallback only when the primary cannot
   issue a verdict or the operator requests confirmation. There should be no
   audit for every tiny edit. State whether this policy is sufficient and
   whether the current records follow it.
5. Design a lean, production/publication-safe operating policy: a small set of
   fast pre-commit controls, one native bundle gate, one end-of-bundle audit,
   evidence retention/index rules, time budgets, and explicit escalation
   triggers. Preserve scientific safety and fail-closed behavior.
6. Identify what should be consolidated, archived, or moved to long-term G5/G6
   work. Do not perform those mutations; give exact paths/categories and a
   safe operator-approved migration sequence.

Do not propose skipping conservation, native compiler/link checks, MPI/CUDA
coverage where the production path depends on them, physical-model approval,
or the small initialized-RAMSES gate required before a real feedback run.
Conversely, do not treat every generated file, Python validator, repeated
prompt, or narrative provenance copy as an independent scientific gate.

Return a concise but evidence-backed report with:

- exactly one top classification (`OVERINSTRUMENTED`, `BALANCED`, or
  `UNDER-INSTRUMENTED`);
- severity-ranked findings with paths and counts/timestamps where possible;
- a three-tier evidence table: mandatory now / bundle-level / deferred or
  redundant;
- a recommended operating policy and concrete next action;
- an explicit statement of what must not be deleted or weakened without
  operator approval.

This is an operational review only. Do not approve or reject the physical
AGN/stellar SED, yield tables, DTD/PISN, dust physics, feedback production,
HDF5 restart, distributed AMR scaling, or publication science itself.

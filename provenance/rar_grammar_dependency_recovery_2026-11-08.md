# Grammar dependency audit, 8 November 2026

All nine stale pending RAR jobs were retired after verifying the successful
replacement jobs, their recorded scope, and their artifact checksums. Failed
jobs and scientific outputs remain preserved. No failed prerequisite was
declared successful and no scientific acceptance criterion was removed.

| Retired pending jobs | Verified replacement | Reason |
|---|---|---|
| 1267690 | 1267704, supported by 1267712 | Independent FDM analysis already handled the serial-only log-marker failure. |
| 1268670, 1268673, 1268677, 1268683, 1268688 | 1268706 after build 1268696 | The corrected build and bounded-cloud reader supersede the five failed build branches. This reader alone does not validate active AGN feedback. |
| 1275336 | 1275344 after run 1275343 | Corrected synchronous-subcycle history test already passed. |
| 1275358, 1275359 | 1275361, 1275362 after build 1275360 | Corrected module dependencies already passed the bounded composite-operator test. |

The exact controller records and replacement-file digests are in
`results/rar_grammar_dependencies_20261108/retired_jobs.json`. Cancellation
was restricted to these exact pending IDs after rechecking owner, project,
dependency and state. This removes obsolete scheduler entries; it is not a
new simulation or expanded physics validation.

The user-wide Grammar audit also identified four completed duplicate DE
analyses and retired their pending predecessors: 1258610 and 1260851 are
superseded by 1267635 and 1267636; 1260396 and 1267631 by 1260833. The latter
completed analysis retains its negative one-percent scientific gate. The DE
results remain in their separate DE project. Full operational records are
under `/gpfs/kjhan/scheduler_audits/grammar_dependency_recovery_20261108_v1`
and `_v2`, with local copies in the RAR workspace's
`reports/grammar_dependency_recovery_20261108/`.

After these actions, the remaining queue has 38 DependencyNeverSatisfied
entries, 23 Dependency entries, one JobHeldAdmin entry, and four running
jobs. These belong to other projects. Examples of unresolved prerequisites
include failed source/phase provenance, failed numerical or physical gates,
an infeasible F5/Takahashi optimizer result, and timeouts. Their dependencies
were preserved; this audit does not claim those chains are recovered. The
three active VoidSim chains, the FDM assembly chain, and the PPI admin hold
were preserved.

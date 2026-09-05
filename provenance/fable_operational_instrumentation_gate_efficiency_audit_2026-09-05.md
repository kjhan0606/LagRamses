# Fable operational audit — instrumentation and gate efficiency

Date: 2026-09-05 (Asia/Seoul)
Workspace: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`, `main`)
Auditor: Fable
Mode: read-only; no code, test, build, job, deletion, or archive mutation

## Classification

`OVERINSTRUMENTED`

The bookkeeping/audit layer has grown faster than the native implementation.
At the same time, one important control is missing: there is no initialized-
RAMSES SNRT live smoke. The project therefore has too much low-value
instrumentation and too little of one high-value native integration test.

## Evidence-backed headline numbers

| Measure | Observed value |
|---|---:|
| Provenance files / size | 298 / 1.4 MB; no index |
| Auditor-named files | 187 (63%); 97 prompt files |
| Unreferenced provenance files | 138 (46%) |
| Commits since Sep 1 | 37; 9 touched native code, 19 touched no code |
| Uncommitted work | 64 modified + 116 untracked files |
| Python validators/tests vs native production | ~58k vs ~10k lines |
| Hand-run Python scripts | 97 of 99; only 2 pytest-collectable |
| Native smoke programs / shell runners | 19 / 14 |
| Initialized-RAMSES SNRT smoke | 0 namelist, sbatch, or Makefile target |

For the F-P2 series, seven bundles closed in 13.7 hours with 42 provenance
files, 15 audit reports, and 13 prompts. F-P2.6 used five audit events despite
the intended two-event policy; its rough wall-clock split was 30% first-pass
implementation and 70% audit/remediation. The first failed audit found real
defects, but the later follow-up/closure audits mostly rechecked single low-
severity fixes. An Opus model rejection plus a 900-second timeout was also a
recurring fixed cost.

## Severity-ranked findings

1. **High — audit governance drift.** The one-plan/one-end-audit/fallback-only
   policy was not followed. Re-audits of individual fixes became the dominant
   delay. The policy itself is sufficient for safety.
2. **High — record volume exceeds implementation value.** Many commits contain
   only record/approval changes; reports are duplicated in
   `~/.claude/plans` and `provenance/`, while provenance has no index.
3. **High — missing high-value integration gate.** The F-P2.6 route markers are
   `grep`/`rg` source checks, not a live initialized-RAMSES execution. This is
   the principal confidence gap to fix.
4. **Medium — stale duplicated native source.**
   `simulation/snrt/native/phase0/` mirrors 21 stellar modules; nine diverge
   from `patch/lagRamses`, including a 625-line runtime difference, and five
   runners compile the mirror.
5. **Medium — Python validator sprawl.** There are many hand-run/orphaned
   bookkeeping tools and JAX-derived artifacts that are not native RT gates.
6. **Medium — superseded variants.** Numerous `_final`, `_remediated`, `_v2`,
   and similar JSON variants are retained without an index.
7. **Low — generated/disk products.** Fable observed a 299 GB partial HDF5,
   1.1 GB venv, large validation trees, root `*__genmod.f90` products, and
   duplicated candidate binaries.
8. **Low/medium — dirty-tree risk.** F-P2 through F-P2.6 are uncommitted, so
   evidence hashes do not yet resolve to durable commits.

## Evidence disposition

| Tier | Keep/operate as |
|---|---|
| Mandatory before live feedback | Native Fortran thermochemistry, spectral/checkpoint, transaction/MPI, and CUDA smokes; conservation checks; `mpiifx`/CUDA production link and symbol check; `git diff --check`; hashes; physical-model approvals; **new initialized-RAMSES smoke** |
| Once per bundle | One consolidated native bundle gate and one plan audit plus one end audit; asset validators only when their guarded asset changes |
| Defer/archive after operator approval | Single-fix re-audits; static route grep counted as a test; orphaned scripts; superseded JSON; duplicate prompts/reports; JAX prototype output; duplicated candidate binaries |

## Recommended lean policy

- Pre-commit, target under five minutes: `git diff --check`, changed native
  object compile, the smoke for the changed module, and changed-source hashes.
- Once per bundle, target under 30 minutes: one bundle-gate script containing
  the native smokes, MPI/CUDA coverage, full production link/symbol check, and
  eventually a two-rank initialized-RAMSES smoke for a few steps.
- One plan audit before and one end audit after each bundle. Time-box the plan
  audit to 15 minutes and end audit to 60 minutes. Permit at most one
  re-audit, only after a consolidated repair. A failed primary attempt should
  trigger the fallback directly.
- Add `provenance/INDEX.md`; keep one report per audit with the prompt embedded;
  put superseded artifacts under an explicitly marked archive directory.
- Target implementation >=60%, gate <=15%, audit <=25% of bundle time. If the
  audit share exceeds 25%, enlarge the next bundle.

## Recommended next action

Do not open another physics bundle yet. After operator approval, create one
short gate-consolidation bundle: commit the dirty tree, add the provenance
index, merge the native runners into one Makefile-backed bundle gate, resolve
the `phase0` mirror ownership, and build the small initialized-RAMSES SNRT
smoke with two ranks and diagnostic failure injection. Do not weaken native
conservation, MPI/CUDA, compiler/link, fail-closed, physical-model approval,
or pre-feedback integration requirements.

No deletion, archive move, code change, gate reduction, or simulation launch
was performed from this audit. In particular, the partial HDF5, JAX output,
virtual environment, generated modules, candidate binaries, physical approval
records, and existing native smokes require an explicit operator decision
before mutation.

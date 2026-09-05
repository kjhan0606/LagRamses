# Fable F-P2.6 native RT/chemistry transaction bundle closure audit

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Fable (governed fallback)
Mode: read-only source/evidence review; no jobs or files were modified by the auditor

## Verdict

`CONDITIONAL PASS`

The primary Claude Opus 5 audit did not issue a verdict: `opus-5` was rejected
by the local catalog and the `opus` alias timed out after 900 seconds. Fable
therefore performed the fallback closure audit. No blocker survives in the
current F-P2.6 native transaction boundary.

## Severity-ranked findings and disposition

### R1 — stale binary hash, medium evidence condition — closed

The evidence record had retained an earlier `bin/ramses_final3d` SHA-256 even
though the final hardening rebuild produced `b5160fb481bc33b43b7640736f0c8df837b90af11e139131b87f21517c3ba44d`.
The hash was recaptured and the evidence record was corrected. The linked
binary exports the native driver, prepared transport, transaction reduction,
and thermochemistry symbols and contains the final-hardening strings.

### R2 — complete driver fixed-point loop lacks live initialized-RAMSES evidence,
medium evidence condition, record-only

The loop in `patch/lagRamses/snrt_ramses_driver.f90` was inspected and its
frozen inputs, time-centred opacity, under-relaxation, final-trial commit,
no-photon-reuse rule, collective decisions, clean-stop behavior, and hard
32-trial bound were judged coherent. The current native smokes execute the
explicit transaction/norm/bisection controls and the source-level route check,
but do not run the full loop inside an initialized RAMSES AMR state. This is an
evidence limitation, not a found code defect, and is carried to the small
initialized-RAMSES harness/live feedback bundle.

### R3 — low diagnostic observability, record-only

The per-leaf loop retains the failure class at the level boundary rather than
preserving a structured first-failing leaf, and `trial_unassigned` is an
internal diagnostic array while the tolerance-sized residual is exposed in the
global printed ledger. This does not permit a failed transaction to commit;
structured failure/ledger diagnostics are future hardening.

### R4 — low AMR correction hardening, G5 record-only

The reviewed collective path is correct, while coarse-face/no-slot accounting
and a post-commit non-negativity guard require further distributed-AMR
coverage. These are not needed to establish the current transaction boundary
and are assigned to G5.

### R5/R6 — low configuration/initial-opacity policy records

The production environment is expected to be rank-uniform before collective
entry, and initial opacity reconstruction reports/handles errors differently
from in-loop opacity errors. Both policies are visible in source and do not
create a blocker for this bundle; they remain documented hardening records.

### R7 — governance record

F-P2.6 remains dirty-tree work alongside unrelated project work and is not yet
represented by a commit at this audit point. No commit or push was requested
in this turn; repository integration is a separate operator-controlled action.

## Original finding disposition

| Finding | Disposition |
|---|---|
| B1, H1, H2, H3 (prepared path), H5 | closed |
| H4 | closed at evidence level; route check is static, not a live injected evolution |
| M1–M4 | closed in native source/evidence |
| M5 | closed after the final binary hash refresh |
| L1–L6 | closed/documented |
| L7 | record-only G5 no-slot coarse correction accounting |
| N1 | closed by collective non-finite hydro pre-source clean-stop |

## Gate disposition

- C1 transaction/rollback: `PASS`.
- C2 unassigned receiver gate: `PASS`.
- C3 bounded fixed point: `PASS` by source and native controls, with R2
  record-only live-loop evidence gap.
- C4 native evidence/build: `CONDITIONAL` only for R2; the recorded native
  smokes, MPI path, CUDA path, GNU/`mpiifx` compilation, full link, symbols,
  hashes, and `git diff --check` pass.

Acceptance bullets 1, 2, 4, and 6 pass. Bullet 3 passes by source/native
control evidence with the R2 live-loop limitation. Bullet 5 is conditional in
the audit wording only because the same evidence limitation covers the final
receiver path inside a live RAMSES state.

## Explicit non-approvals

This audit does not approve a physical AGN or stellar SED, yield tables,
SNIa/PISN or DTD completeness, a live feedback production run, dust physics,
HDF5 restart integration, distributed AMR scaling, or publication-scale
convergence.

# SIDM project memo

Last updated: 2026-08-24 KST

## Always first

- Run `hostname` before inspecting or launching jobs.
- Do not place CPU-intensive work on `syntax`. Use `grammar` for builds and
  validation and use Lageunha only for an explicitly selected production run.
- Read `SIDM_HANDOFF.md` before changing a live namelist or executable.
- The GPFS outputs are canonical. Files below `runs/lageunha/snapshot-2026-08-24`
  are a read-only metadata snapshot and are not restart data.

## Next production task

Resume the matched Paper-I zooms after Lageunha resources are available.

1. SIDM1 must restart from `output_00010`, not the stale `nrestart=7` value.
2. SIDM3 must restart from `output_00009`, not the stale `nrestart=8` value.
3. Confirm both stored checkpoints contain the expected files for every MPI
   rank and preserve the existing executable before replacing a shared target.
4. Clear `jobcontrol.txt` only immediately before a deliberate relaunch.
5. Retain the existing affinity of 32 MPI ranks by 2 OpenMP threads for SIDM1
   and 24 MPI ranks by 2 OpenMP threads for SIDM3 unless a new benchmark
   justifies a change.
6. Continue both runs to a common epoch before making a matched profile or map.

## Canonical locations

- Project and code worktree: `/home/kjhan/BACKUP/lagRamses-SIDM`
- Branch: `sidm`
- Paper repository: `/home/kjhan/BACKUP/lagRamses-SIDM/paper`
- Paper compatibility link: `/home/kjhan/paper_sidm_overleaf`
- Live simulation root on Lageunha: `/gpfs/kjhan/Hydro/Sidm`
- Full handoff: `SIDM_HANDOFF.md`

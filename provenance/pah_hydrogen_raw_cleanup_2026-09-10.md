# PAH hydrogen-state raw-output cleanup

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
Operator's completed-test retention directive applies. Native physics and
restart evaluations finished; no listed checkpoint is needed by an active
run. Only the exact raw files below are removed, not run directories,
effective inputs, logs, executables or compact evidence.

Base: `/gpfs/kjhan/LRD_JWST/.pah-hydrogen.6fDwRw/`.

| Output directory | Disposition |
| --- | --- |
| `live/output_00001`, `live/output_00002` | Initial 8^3 general solver. Conserved physics; incomplete PAH checkpoint identity, not restart evidence. |
| `fast/output_00001`, `fast/output_00002` | Initial 8^3 optimized solver. Exact general-solver state parity; same incomplete checkpoint identity. |
| `restart.9plAuy/output_00001` | Copy of incomplete initial checkpoint; restart correctly rejected missing identity. |
| `checked.2XffQ5/output_00001`, `checked.2XffQ5/output_00002` | Corrected 4^3 live integration, complete checked identity. |
| `replay.ITn60C/output_00001`, `replay.ITn60C/output_00002` | Corrected restart; all hydro/RT/gravity/AMR/PAH identity datasets match continuous execution exactly. |

Remove `data_00001.h5` / `data_00002.h5` and the matching
`cooling_00001.out` / `cooling_00002.out` inside only these nine directories.
Preserve metadata, including effective dumped namelists, build provenance,
headers and descriptors. `COMPLETE` refers to the original output event;
this manifest records that its raw payload was subsequently removed.

Retained evidence: `evaluation.txt`, `evaluate.py`, final Intel/GNU native
logs, GUI results, all run logs/inputs/environments, build logs and binaries;
see `medium_physics_implementation_2026-09-10.md` for physics limits and
binary SHA256. Removal is permanent; regeneration requires rerunning the
retained fixtures. Byte totals are recorded after execution below.

Completed: 18 exact files removed, 470,948,520 bytes (~449 MiB). Effective
inputs, logs, metadata and compact evaluation remain in all five run folders.

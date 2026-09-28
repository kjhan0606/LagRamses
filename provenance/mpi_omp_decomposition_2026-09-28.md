# Fixed-total-core MPI/OpenMP layout comparison (2026-09-28)

Project `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`. This is a
single-run technical layout comparison on `lageunha`, not a production-scale
scaling qualification. The first comparison below reuses the 8^3, 512-cell,
two-coarse-step, source-active M5/CHIMES/dust input from the preceding OpenMP
tail experiment. A separate 16^3 follow-up is recorded below.

All layouts used the same 28 distinct physical cores (CPU IDs 1--28 on one
socket). CPU 0 was excluded because another user's process occupied it. The
effective namelist is
`.mn-scale-high-20260928/openmp-1/effective.nml` (SHA256
`61a867b5046ac14b131cc44092d7ca89811dc9ca24070b4a3e6818180552f7bd`); the
binary is `.mn-omp-tail-20260928/ramses_window3d` (SHA256
`914f0d7ac261af3f57d95cf491c19944529f15e56d4d97f3ac7a919f88564182`). The
runner explicitly exported `OMP_NUM_THREADS`, used `I_MPI_PIN_DOMAIN=omp`,
`I_MPI_PIN_ORDER=compact`, `OMP_PLACES=cores`, and pinned ranks to CPUs 1--28.
The tested layouts were 1x28, 2x14, and 4x7. Logs and effective inputs are
under `.mn-omp-tail-20260928/mpi-omp-matrix-20260928-verified/`.

| MPI ranks x OMP threads | First step (s) | Source-active second step (s) | Full executable time (s) | Max-rank photo prep (s) | Summed CHIMES photo worker time (CPU-s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1x28 | 9.33 | 36.94 | 46.268 | 21.907 | 547.092 |
| 2x14 | 4.71 | 35.97 | 40.677 | 26.466 | 546.644 |
| 4x7 | 2.57 | 30.21 | 32.772 | 23.751 | 546.698 |

At fixed total physical cores, 4x7 reduced the source-active step wall time by
18.2% relative to 1x28 and 16.0% relative to 2x14. The 2x14 result was only
2.6% faster than 1x28. These are integrated coarse-step times from the RAMSES
log, not isolated M5 kernel timings. The photo worker-time sum stayed within
0.1% across layouts, so the faster step does not indicate less chemistry
work. Rank-level preparation was more even with four ranks: max/min was 1.22
for 4x7 versus 1.36 for 2x14. The most expensive individual photo cell still
took about 13.8--14.4 s in every layout; MPI decomposition did not remove
that within-cell tail.

Physics checks: all three runs printed `Run completed`, had empty stderr, no
full `output_*` dumps, zero `NaN_CHK` counters, and the same 511 new stars in
the source-active step. Relative to 1x28, 2x14 differed by 2.3e-6 in absorbed
photons and 1.6e-6 in gas heat; 4x7 differed by 1.4e-10 and 1.3e-9,
respectively. The photon-ledger residual was about 1.2e-15 of emitted source
photons, and energy-ledger residuals remained below 1e-12 of absorbed energy.
Thus the layouts conserve well but are not bitwise identical; MPI/OpenMP
reduction/decomposition order causes small numerical drift.

Disposition: 4x7 is the best of these three layouts for this exact 512-cell
control. Keep it as a candidate, not a universal/default optimum: there was
one run per layout and the grid is too small to establish production behavior.
The next meaningful performance target remains the expensive per-cell
CHIMES photo solve/hot-cell tail, followed by confirmation on a genuinely
larger active-cell workload; do not expand this into a broad benchmark matrix.

An initial invalid runner attempt omitted the `OMP_NUM_THREADS` export and
therefore did not execute the intended 1x28 layout. It was stopped, excluded
from every result above, and its log/input were retained for provenance. The
corrected runner set and verified the thread count for all three reported
layouts.

Cleanup manifest: after evaluation, removed only
`agn_coarse_state_v1.jsonl` and `sink_log.csv` from each of the three verified
run directories. The invalid attempt produced neither file. Logs, effective
namelists, `ic_sink`, runner, binary and build identity were retained; no full
simulation snapshots were produced. No commit or push was made.

## 4096-cell follow-up and CHIMES RHS micro-optimization screen

Follow-up run on 2026-09-28, on `lageunha`, using the same 28 physical cores
(CPU IDs 36--63) and the same binary, physics input, and two coarse steps for
all three layouts. This is a uniform 16^3 = 4096-cell grid (`levelmin=levelmax=4`),
not an AMR or production scaling test. Effective input:
`.mn-omp-tail-20260928/mpi-omp-larger-control-20260928/effective.nml`
(SHA256 `10c7131aafd041c75749680378911395c7e05b7dcf920328e9011e9172cf09c8`);
`ic_sink` SHA256
`33d26c9d0b532ab897591a4d1046e36e6e9116a5c60cf1f3a63e3d3c4b832a49`. Binary
SHA256 remains
`914f0d7ac261af3f57d95cf491c19944529f15e56d4d97f3ac7a919f88564182`.
Runner: `.mn-omp-tail-20260928/run_mpi_omp_larger.sh`. The source file
`patch/lagRamses/snrt_chimes_photo.cpp` was restored after the scalar trials
and has SHA256
`492853e189f9cf45f47dd2912c99e628bb6cffc4a97e6c651816ffdcab0b5948`, matching
the pre-trial baseline.

| MPI ranks x OMP threads | First coarse step (s) | Source-active second step (s) | M5 commit (s) | RAMSES internal total (s) | Runner outer wall (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1x28 | 73.78 | 273.12 | 261.500 | 346.912 | 562 |
| 2x14 | 38.87 | 203.74 | 197.376 | 242.602 | 457.443 |
| 4x7 | 19.82 | 170.00 | 166.245 | 189.815 | 407 |

At fixed core count, 4x7 reduced the source-active coarse-step time by 37.8%,
the M5 commit time by 36.4%, and the executable's internal total by 45.3%.
The outer wall reduction was 27.6%; it includes roughly 215 seconds outside
the reported RAMSES internal timer. That timer starts on the first call to
`update_time` (`amr/update_time.f90`), after startup work. Module loading is
outside the runner's timed command and cannot explain this difference;
the remaining startup/finalization costs have not been separately measured.
The source-active step is the more useful comparison for evolution cost.
This materially reinforces 4x7 as a candidate for larger
active-cell work, but each layout was run once and this single uniform grid
does not establish a universal optimum.

Both runs created the same 2023 stars in step one and 2012 in step two (4035
total), completed normally, had empty stderr, zero reported `NaN_CHK` counters,
and produced no full `output_*` dumps. The absorbed-photon ledger differed by
2.36e-6 relative and gas heat by 2.49e-6 relative; photon residuals were
8.25e-15 and 6.69e-15 of source photons, and energy residuals were 4.94e-12
and 2.00e-12 of absorbed energy, respectively. Thus the runs agree closely
and close their ledgers, without being bitwise identical. The per-tile CHIMES
worker-time rows had different rank coverage between these two layouts, so
their summed worker time is deliberately excluded from this comparison.

Three narrow native CHIMES RHS trials were also screened against the exact
baseline object and fixture:

| Candidate | Baseline mean (s) | Candidate mean (s) | Change | Disposition |
| --- | ---: | ---: | ---: | --- |
| Precompute capture scales | 14.981 | 16.539 | +10.4% slower | Rejected |
| Precompute `c*dt*sigma` and regroup preconditioner | 14.977 | 15.634 | +4.4% slower | Rejected |
| Precompute `c*dt*sigma` in RHS only | 14.997 | 14.826 | 1.14% faster | Rejected as immaterial |

All passed the bounded native physics test; the third also left the reported
budget output unchanged. None was retained: two regressed, while the third
gain is below the project's practically meaningful ~5% threshold. No source
optimization from these trials remains in the working tree.

Cleanup manifest: after evaluation, removed only the four raw per-run files
`agn_coarse_state_v1.jsonl` and `sink_log.csv` from
`runs/mpi1-omp28/` and `runs/mpi4-omp7/` below
`.mn-omp-tail-20260928/mpi-omp-larger-control-20260928/`. Their pre-deletion
SHA256 values were, respectively, `01a11d354004fb0e7487aeb1bbd19d37d0586d2b64f741d93926862fc5321c06`,
`58ffcd9f6fcf4ef217514118a7889ba40eee09d5396052aea1471933c6ddb300`,
`a60cd8bc9eb106cd460bcde8b4793cafaa908a281a49dd2de2ae9f5b8ba96215`, and
`40756684320a8b9bfd7e8301f008984ddd9b4a120b047d6f501c937802f50301`.
Effective inputs, run logs, empty stderr captures, binary, and runner are
retained; no full simulation snapshots were created. No commit or push was
made.

### Completed middle-layout measurement (2x14)

The authorized 2x14 follow-up finished normally on Lageunha at
2026-09-28 17:01:54 KST. Its runner is
`.mn-omp-tail-20260928/run_mpi_omp_larger_2x14.sh`; effective input and logs
are under `mpi-omp-larger-control-20260928/runs/mpi2-omp14/` in that directory.
Both input hashes and the binary hash match the controls above. CPU IDs
36--63 are distinct physical cores on socket 1; both MPI ranks were observed
with 14 threads during evolution. The other long-running user process had
affinity restricted to CPU 0. The effective namelist's introductory comment
still says 8^3/MPI1/OMP4, inherited from its template; actual resolution is
set by `levelmin=levelmax=4` and layout by the runner. Its bytes were retained
to keep the three inputs identical.

The source-active step fell 25.4% from 1x28 to 2x14. Going from 2x14 to 4x7
reduced it a further 16.6%. Thus the measured ordering is 4x7, 2x14, 1x28;
4x7 remains the preferred candidate for this bounded workload. This is a
fixed-total-core decomposition comparison, not an OpenMP strong-scaling
curve, and one sample per layout does not establish statistical confidence
or a production default. No further matrix expansion is needed for this
comparison.

The 2x14 run again formed 2023 then 2012 stars, had empty stderr, zero NaN
counters, two M5 commits and no full dumps. Relative to 1x28, absorbed
photons differed by 1.157e-5 and gas heat by 1.221e-5 (about 0.0012%).
Photon residual/source was 1.021e-14; energy residual/absorbed was
1.357e-11. These compact checks support numerical agreement, not a new
scientific convergence claim. Bash timing recorded outer wall 457.443 s,
launcher-child user 3928.866 CPU-s and system 14.550 CPU-s; comparable
process CPU counters were not collected for the older two runs, so no
cross-layout CPU-time reduction is asserted. The 214.841 s difference
between outer wall and RAMSES internal total is consistent with the
earlier unseparated startup/finalization costs.

Cleanup: after evaluation, removed only `agn_coarse_state_v1.jsonl` (7166
bytes, SHA256 `fe62c61ddb41cf7fdbb8ef783fcc98bb78ca0680939735adf0f906ed2a8146d3`)
and `sink_log.csv` (268 bytes, SHA256
`0afddfe8d35fa29b70b508e08ac038c58d409a9d3cafe1168afa5a3149996d94`)
from the exact `runs/mpi2-omp14/` directory above. Inputs, run/launcher logs,
timing, executable and runner remain available. These raw files have no
separate backup; reproducing them requires rerunning the retained input.
No simulation-source changes, commit or push were made in this follow-up.

Closeout: the authorized layout comparison is complete. There is no pending
measurement required by this comparison. The expensive per-cell CHIMES RHS
remains an optional future optimization target; the rejected scalar trials
do not justify extending the current task or declaring a new completion gate.
Commit/push requested after evaluation; the statements above describe the
run-time disposition before publication of this evidence.

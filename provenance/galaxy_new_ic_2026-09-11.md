# Newly generated small-box hydrodynamic IC

Operator authorized new IC generation on2026-09-11. First128^3 input completed;
RAMSES reader/evolution and a selected-halo zoom have NOT run.

Working directory: `/gpfs/kjhan/LRD_JWST/.galaxy-ic.4BUKHc`.
Input directory: `ics_l7/level_007`. Retain these science inputs, not test-output
cleanup candidates. `l7.conf`, `transfer.py`, `transfer_metadata.json`,
`generate.sbatch`, logs and `SHA256SUMS` are preserved alongside the ICs.

## Model and intended use

- Periodic12.5 comoving Mpc (=8.425 Mpc/h), z_init99,128^3 gas cells and DM
  particles. Omega_m=.315, Omega_b=.049, H0=67.4, ns=.965, sigma8=.811.
- Massless-neutrino LCDM approximation: Neff3.046, no massive neutrino species;
  not a claim to exactly reproduce COLIBRE cosmology. MUSIC's matter+Lambda
  growth/normalization approximation is retained; CAMB includes early radiation.
- CAMB1.6.5 standard13-column z99 transfer, separate CDM/baryon density and
  velocity. As=2.018291882648546e-9 gives CAMB sigma8(z0)=.8110000000000004.
  k range1.1027338587e-5--461.6776428 h/Mpc covers even512^3 corner Nyquist.
- lagMUSIC `camb_file`, grafic2,2LPT+LLA, periodic k-space transfer.
  Seeds at levels7--12 are20260911--20260916. Preserve the hierarchy for
  later resolution variants; cross-resolution phase agreement is not yet measured.
- Initial gas mass5.752796e6 Msun/cell; DM mass3.122947e7 Msun/particle.
  Base spacing97.65625 ckpc, .9765625 pkpc atz99. These are NOT star-forming
  gas resolution or the ultimate few-pc zoom, and stars are not initialized.
  Fine DM sampling requires new refined ICs, not merely raising AMR levelmax.

This is a bounded small-box startup/target-selection input, not the full fine
tuning ensemble. Halo selection and its Lagrangian region remain future work;
no arbitrary geometric subcube is called a self-consistent halo zoom.
Primordial gas chemistry/temperature, table-domain treatment and cosmological
source/seeding applicability must be explicit in the effective hydro namelist
before execution. The IC writer's defaults are not the physics namelist.

Algorithm reference: [MUSIC author documentation](https://www-n.oca.eu/ohahn/MUSIC/).

## Generator identity and inspection

MUSIC copied from
`/home/kjhan/BACKUP/LagMUSIC/music/build-paper-Ib-lageunha-intel-extended/MUSIC`.
Binary SHA256: `345843b58c33c2555b234ecf797b829043a858031743ef4ac97370190c2a0a6d`.
External checkout HEAD35a4eb978ae0e84673d86a0e4722421118277700 has local
`output_grafic2.cc` changes (48-byte omega_b header and vector-write correction).
No external source was changed. Binary hash, rather than a clean-commit claim,
is the executable identity. Output inspection confirms the extended header.

All10 fields passed full record-marker/extent/EOF checks,128^3 dimensions,
48-byte header, finite data and matching a=.01, Omega_m=.315, Omega_L=.685,
H0=67.4, Omega_b=.049 (float32 precision). Gas overdensity range
[-.2507494092,.4026301801], mean1.5682e-13: all1+delta>0. Gas and CDM
velocity arrays differ. Every SHA256 entry verified. IC directory83,905,072
bytes; whole workspace approximately102 MiB. This is file validation, not
a RAMSES reader or physical evolution qualification.

## Resource evidence, including failed setup attempts

grammar/debug, grammar-debug,1 MPI process x8 OMP threads,32 GiB ceiling,
30minute limit. Slurm compute-step accounting below; do not add job-parent
TotalCPU to its child steps. Launcher Bash user/system times are NOT rank CPU.

| Job | Result | Allocated core-seconds (job) | Compute-step CPU seconds |
| --- | --- | ---: | ---: |
|538954|Transfer k-coverage assertion; no IC written|40|19.409|
|538955|Transfer completed; MPI PMI initialization failed before IC writing|40|26.717+.021|
|538956|MUSIC plugin name mismatch, no IC written|24|1.172|
|538958|MUSIC completed, exit0|16|2.783|

Resolved setup issues: kmax300/Mpc; explicit Slurm PMI2 and
I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so; registered plugin `camb_file`.
Successful CAMB wall3.957s; MUSIC wall1.054s, MaxRSS231088KiB (Slurm sample,
not guaranteed true peak for such a short job). Total allocated cost including
failures120core-seconds=.03333core-hours; compute-step CPU50.102s=.01392h.
These short runs do not establish MPI/OpenMP scaling, late-time memory or
galaxy-evolution cost. The campaign's fixed-input MPI/OpenMP comparison remains
attached to a sufficiently long RAMSES evolution interval, not IC generation.

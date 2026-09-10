# Gas MHD runtime

Ideal, face-centred constrained transport (CT), unsplit MUSCL-Hancock,
HLLD face fluxes and 2-D HLLD edge EMFs. Magnetic AMR interpolation and
coarse/fine EMF synchronization preserve the discrete divergence constraint.
This reuses the RAMSES MHD method; it is not divergence cleaning.

## Build and select

Use a fresh build directory directly under this repository (do not mix
hydro/MHD, CUDA/non-CUDA, or differing NVAR objects):

```sh
make -f ../bin/Makefile -j1 SOLVER=mhd SNRT=1 DUST_LIVE=1 HDF5=1 USE_FFTW=0 NVECTOR=32 EXEC=ramses_mhd ramses
```

Add `USE_CUDA=1` for optional device face batches. CUDA architectures remain
those in the Makefile unless explicitly limited to the installed GPU.
VPATH order is unchanged; MHD objects are selected explicitly. Set:

```fortran
&HYDRO_PARAMS
 mhd_enabled=.true.
 mhd_seed=1d-4,2d-4,3d-4
 scheme='muscl'
 riemann='hlld'
 riemann2d='hlld'
 pressure_fix=.false.
 mhd_omp=.true.
 mhd_gpu_faces=.false.
/
```

`mhd_seed` is **code-unit B at the initial epoch**, not gauss or today's
comoving field. Magnetic energy density is B^2/2. The cosmological
supercomoving scaling is tested separately with uniform GRAFIC initial data.
For actual cosmological initial conditions supply the proper IC fields and
units; the tiny noncosmological examples are not cosmological production ICs.

Use `gpu_hydro=.false.`: that older GPU kernel is for nonmagnetic hydro.
HDF5 output/restart is required. `mkrun.py` and the shared CLI/GUI schema
expose the same flags. Default MHD threading and CUDA flags remain false.
An example gas-only AMR input is `config/gas_mhd_ct_amr_smoke.nml`.

## Hybrid work unit and ownership

With `mhd_omp=.true.`, each CPU worker handles an active-grid batch, including
reconstruction and CT edge EMFs. Scratch arrays are private. Cell/face updates
and AMR reflux share a protected scatter region, including neighboring coarse
faces. MPI exchange occurs outside those worker regions.

With `mhd_gpu_faces=.true.`, one direction's reconstructed **HLLD face batch**
is submitted to the existing shared stream pool. A free slot executes actual
CUDA HLLD; an occupied/unavailable pool immediately returns that work to the
CPU worker. There is no calibration, latency-tuned routing, or waiting for a
busy stream. GPU copies synchronize only their acquired slot before release;
CUDA execution errors stop visibly rather than committing partial fluxes.
Faces selected by the existing high-field/density LLF safety rule still use
the CPU safety solver. Dependent dust/chemical fluxes are reconciled after
either backend. Reconstruction, edge EMFs and CT/reflux remain CPU work:
this is **partial numerical offload**, not a whole-step GPU solver.

`n_cuda_streams` is the existing runtime pool size. The pool distributes MPI
local ranks over detected devices and establishes the device in each OpenMP
thread. A single MPI rank's pool stays on one device; do not describe that
as a single-rank multi-device scheduler. Use MPI ranks for multiple GPUs.
`MHD_DISPATCH` reports rank-local cumulative real face counts and busy/absent
pool fallback counts. Reports occur at the base-level `ncontrol` cadence, not every fine
substep. Use sufficient `OMP_STACKSIZE`/`KMP_STACKSIZE` for private
MHD stencils; the NVECTOR=32 checks use 512M. Larger NVECTOR/NVAR increases
per-worker memory. Device mode currently requires NENER=0.

## Qualified scope versus physical approval

Tested: periodic 3-D gas MHD, smooth Alfven and Brio-Wu stability checks,
uniform expansion, mixed AMR/subcycling/load balancing, HDF5 restart including
MPI 2->1; and the existing stellar/RT/coadvected bulk-dust or stellar/AGN/RT
reference profiles at NENER=0. AGN thermal caps exclude unchanged B energy;
sink accretion retains gas-mesh face flux. No jet-field injection or magnetic
flux swallowing model is inferred. Thermal floors exclude B energy.

The bulk-dust mass model's existing no-sink/noncosmological limit and live IR's
`nremap=0` limit remain. Coupled IR AMR uses all levels (unset SNRT_RT_LEVEL),
HDF5, and its existing reference-control contract. These restrictions are not
lifted by MHD. The BPASS/feedback independent-population and dust material
reference-control caveats remain; numerical passes do not authorize new
physical input tables or a universal production profile.

Unsupported/not qualified here: SGS, separate dust dynamics with MHD, all
optional chemistry populations, arbitrary physical boundaries, NENER>0
combined profiles, and large cosmological galaxy-formation production runs.
Dust charging/Lorentz force, magnetic backreaction of grains, anisotropic
CR diffusion/streaming and nonideal MHD remain explicit long-term physics.

## Evidence and retention

See `../../provenance/gas_mhd_completion_plan_2026-09-10.md` and the linked
cleanup inventory. Native HLLD CPU/device comparison is buildable with the
same CUDA build options and target `mhd_hlld_device_smoke`, then run that
executable on a GPU node. It includes a deliberately occupied stream test.
Application parity uses the existing gas and coupled inputs, not a separate
Python implementation of MHD. Finished/evaluated raw dumps are deleted per
operator instruction; inputs, logs, measured summaries and build identities
are retained. Audit the output schedule before every new RAMSES launch.

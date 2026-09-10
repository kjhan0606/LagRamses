# Gas MHD: approved scope and implementation

Historical initial-bundle record. Later completion work supersedes the
CPU-serial/sink exclusions and pending-test statements below; see
`gas_mhd_completion_plan_2026-09-10.md` and `simulation/snrt/GAS_MHD.md`.

Operator approval: implement the recommended ideal gas MHD now. Dust-field
and cosmic-ray-field interactions are **long-term tasks**, not completion
gates for this bundle. No new audit stages or approval requests are added.

Numerics: existing RAMSES face-centred CT, HLLD interface fluxes and 2-D HLLD
edge EMFs, unsplit second-order MUSCL-Hancock; divergence-preserving magnetic
AMR prolongation/restriction and EMF synchronization. VPATH order is unchanged;
the Makefile explicitly selects MHD objects. Legacy hydro remains default.

The current integration uses canonical RAMSES MHD layout: rho, three momenta,
total energy at 1:5, B-left at 6:8, nonthermal energies then passive fields
from 9, B-right at NVAR+1:NVAR+3. Default profile capacity grows by three,
while explicit NVAR is a final MHD NVAR (must budget those fields). Total
energy includes B^2/2, where code B absorbs the Gaussian sqrt(4*pi).
HDF5 stores all six face components and requires identical MHD layout on
restart. A hydro executable rejects an MHD checkpoint and vice versa.

`mhd_enabled=.true.` is an explicit HYDRO_PARAMS contract requiring a
SOLVER=mhd binary, muscl/hlld/hlld and HDF5. `mhd_seed(3)` is a uniform
divergence-free face field in **code units at the starting epoch**, not gauss
or an automatically converted present-day comoving field. Cosmological
supercomoving expansion updates B and its energy using the upstream scaling.

Current build/execution restrictions: 3D, CPU MHD; no GPU hydro, separate dust
phase dynamics, SGS or sink/AGN path. These combinations must not silently run
through hydro-only kernels. Existing isotropic CR advection/pressure is not
field-aligned CR diffusion/streaming. No new magnetic star-formation subgrid
model is asserted. Cooling/RT thermal extraction excludes magnetic energy.

## Native implementation and evidence

Implemented the existing CT solver in the actual lagRamses build rather than
adding a second solver. In addition to explicit object selection, the upstream
MHD cell accesses now use `ICELL_OF` and the current Morton neighbor APIs;
the old `ncoarse + (child-1)*ngridmax + grid` layout is not valid here.
Magnetic face fields participate in initialization, AMR operations, ghost
exchange (including the variable-CPU restart refresh), HDF5 I/O and restart
layout checks. Shared thermal extraction subtracts magnetic energy.
`mkrun.py` and the shared CLI/GUI namelist schema expose the same controls.

Build, from a fresh directory immediately below the repository root:

```sh
make -f ../bin/Makefile -j1 SOLVER=mhd SNRT=1 HDF5=1 USE_FFTW=0 NVECTOR=32 EXEC=ramses_mhd ramses
```

Never reuse hydro object/module files for an MHD build, or vice versa: the
field layout differs. The validated build uses Intel ifx, 3D, NENER=0,
NVAR=21 (24 stored values including right faces), CPU MHD. SNRT is linked
but runtime-disabled in the tests. The MHD Godunov batches themselves are
currently serial within each MPI rank; OMP=2 exercises surrounding parallel
AMR/runtime work, not an OpenMP/CUDA MHD dispatch implementation.

Reproducible input: `simulation/snrt/config/gas_mhd_ct_amr_smoke.nml`.
It is a small, periodic, gas-only local shear problem, not a production
cosmological namelist. Copy it into a fresh run directory and perform the
usual output/physics audit before launching; do not overwrite the evidence.

All evolution runs below used MPI=2, OMP_NUM_THREADS=2 on lageunha,
`SNRT_RT_ENABLE=0`, four coarse steps, HDF5 every two steps, and no particles,
gravity, cooling, star formation, sinks, RT, CR or dust. All outputs are under
`/gpfs/kjhan/LRD_JWST/.gas-mhd.UDobi3/`.

| Check | Evidence directory | Result |
|---|---|---|
| Fixed 8^3 shear, nonzero field evolution | `live.6zoxCY` | Finite state; min rho 0.9999080143, min thermal energy 0.8999643266; mass 1, total energy 0.975; max discrete div B = 0; max field change 0.00726831 |
| Fixed-grid restart, step 2 to 4 | `restart.HsCNvW` | All 94 HDF5 datasets byte-identical to continuous run |
| Mixed AMR, levels 4–5, subcycling | `amr-local.BtfTWE` | 3616 level-4 and 3840 level-5 leaf cells; leaf volume 1; min rho 0.96616618, min thermal energy 0.85204087; max discrete div B 5.33e-15; mass 1.0000000000000002, total energy 0.9749999999999543 |
| Mixed-AMR restart, step 2 to 4 | `amr-restart.UFlVhH` | Same dataset names, shapes and dtypes; all 154 datasets byte-identical |

Thermal energy is measured as `E - |rho*v|^2/(2*rho) - |(B_left+B_right)/2|^2/2`.
Divergence is `sum(B_right-B_left)/dx`, and global integrals use only leaf
cells with their actual volumes. These runs started at total energy 0.975.
Final HDF5 sizes were approximately 176 KiB (fixed) and 1.8 MiB (AMR).

The numerical results above use frozen `ramses_mhd_block3d`, SHA256
`3f428bd40ef3beff1bc9299ea60cac20b6b38063a08bfb96454df80a6ab8e0f0`.
Final `ramses_mhd3d`, SHA256
`32de1c2a34005fc6c51c1ea94cfc04b8019f890fbf0f7407cc64d600045a9aba`,
also includes the subsequently noticed variable-CPU restart ghost-loop
extension to right face fields. That change rebuilt successfully; the
different-MPI-size restart branch is not exercised by the same-size results.

The baseline hydro build also succeeds in `.hydro-mhd-regression.FeJWYA/`
with SOLVER=hydro, SNRT=1, HDF5=1, NVAR=18. Its gas-only HLLC evolution
in `live.wNIWWS` completes four steps with finite fields, min rho 1,
min thermal energy 0.8999999999999998, mass 1 and total energy
0.9049999999999995 (initial 0.905); no magnetic field is present.
Shared CLI/GUI setup tests:
46 discovered, 45 passed, one display-dependent skip (`gui-tests.log`).
This includes explicit MHD controls and the gas-only wizard's required but
inactive `feedback_mode='legacy'` namelist block.

Preserved diagnostic attempts are not passes: `live.dLh8EJ` stopped at a
missing stellar namelist group; `live.JJs86C` exposed invalid legacy MHD cell
indexing and produced negative/NaN states before the blocked-layout repair.
`amr.6WeXvh` evolved correctly but did not exceed its refinement threshold,
so it is not evidence of AMR coupling. None of those outputs was deleted.

### Qualification boundary

Gas ideal-MHD integration, actual mixed-AMR CT evolution, same-size restart
and setup support are implemented and verified as above. This is not a claim
of arbitrary-physics production qualification: cosmological expansion is
wired but not live-tested here; strong-shock/Alfven convergence, physical
boundaries, changing MPI size, NENER>0, and simultaneous cooling/RT/stellar
feedback/dust chemistry require their own application evidence. Existing
HLLD degeneracy handling uses the upstream LLF fallback, not a newly
implemented HLLE solver. No full dust-carrier/chemistry MHD qualification is
implied by permitting a coadvected field layout. These limits are recorded
without making a new series of gates for the approved gas-only bundle.

## Explicit long-term scope

* Dust charging, grain Lorentz force, gyro-motion and gas/grain electromagnetic
  back reaction. Existing neutral drag is not a substitute for these.
* Field-aligned CR diffusion/streaming and associated wave-mediated heating.
  Existing isotropic CR pressure may act on the gas but does not implement
  this transport physics.
* Non-ideal gas MHD (Ohmic, ambipolar, Hall) requires a separate physical
  scope; none is implicitly supplied by an ideal-MHD switch.

Scientific basis: [RAMSES CT](https://arxiv.org/abs/astro-ph/0607230),
[HLLD](https://doi.org/10.1016/j.jcp.2005.02.017). The numerical method does
not claim universal superiority over cleaning or a resolved galactic dynamo.

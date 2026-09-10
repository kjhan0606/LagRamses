# SNRT angular refinement: bounded live implementation complete

Workspace `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
This implements the angular component of approved medium-term group 8.
It does NOT implement intragroup spectral reconstruction or close groups
1--5 or remaining PAH carbon chemistry.

## Delivered runtime choice

| Make option | Polar x azimuthal nodes | Directions | Angular state/work scaling |
| --- | ---: | ---: | ---: |
| `SNRT_ANGULAR_LEVEL=0` (default) | 8 x 10 | 80 | 1 |
| `SNRT_ANGULAR_LEVEL=1` | 16 x 20 | 320 | 4 |
| `SNRT_ANGULAR_LEVEL=2` | 24 x 30 | 720 | 9 |

These are fixed Gauss-Legendre x half-shifted uniform-azimuth product rules,
not adaptive ray splitting or a different RT method. Positive weights
integrate to 4pi. Refined Legendre roots/weights use the standard recurrence
and converged Newton solve ([NIST DLMF 3.5](https://dlmf.nist.gov/3.5)).
Default nodes, weights, direction order and arithmetic remain bitwise equal
to the original routine in HEAD with matched compiler options. Preserving
the literals alone was insufficient: a variable-bound loop changed final
roundoff bits. The fixed 80-ray calculation path is retained explicitly.

`snrt_state` owns dimensions used by transport, source deposition, photon
energy corrections, dust IR, AMR exchange/regrid and checkpoints. The driver
reports polar/azimuthal counts. Backends receive array sizes dynamically;
no new CUDA dispatch policy, spectral coefficients or physics default.
The 720 case is native-tested, not a second full MPI hydro qualification.

Build in a **fresh repo-root build directory** using the usual Make options
plus `SNRT_ANGULAR_LEVEL=1` or `2`. Never switch flags on an existing object
directory: Make does not track compiler-flag changes. VPATH was untouched.
No namelist field changed, so neither `mkrun.py` nor the NML generator needs
a corresponding field edit. Existing wizard-selected executables are NOT
silently rebuilt or replaced by the new comparison executable.

Native checkpoint direction headers and HDF5 cell-width identities already
reject incompatible resolution before state publication/payload read.
No checkpoint-version or new attribute layer was added. Changing resolution
requires fresh radiation state, not reading another level's checkpoint.

## Plan and driver evaluation

[Fable](snrt_angular_fable_2026-09-10.txt) approved Q-GOAL/Q-LEAN and recommended
keeping the change bounded. Adopted shared constants, existing smokes and
existing checkpoint guards. Used generated refined nodes rather than large
literal tables. Retained inexpensive dust/CR/stellar physics from the existing
small fixture: this checks the actual IR dimension, not just primary photons.
No new Python framework or series of extra audit gates. The driver performs
the end evaluation under the latest operator direction.

One necessary initialization fix was exposed by bounds checking:
`init_flow_fine.f90` allocated element names using the preceding table count
`dum`, then indexed them up to `nelt`. It now allocates `nelt`; the actual
input has more elements than that table count. No disabling of bounds checks,
physics suppression or yield-data edit was used to evade the error.

## Measured evidence

Evidence directory `.rt-angular.TwRcDH/`; final executable
`ramses_angular3203d`, SHA256
`ea4a72048b2718eac939a98b97e33b2c5c3fbdb5cf3ab0f251cc26986839595d`.
Intel MPI/ifx, `NVECTOR=32 NENER=1 NVAR=30 SNRT=1 DUST_LIVE=1 HDF5=1
USE_FFTW=0 FDMDEBUG=1 SNRT_ANGULAR_LEVEL=1`, CPU hydro, bounds checks.
Clean parallel build hit existing module-order prerequisites; serial build
and final incremental rebuild passed (`build-serial.log`, `build-final.log`).

- Angular smoke at each compiled level passes normalization, positive weights,
  unit directions, zero first moment, isotropic second moment, antipodal
  pairing, bad-family rejection. For `exp(20*(n.axis-1))`, axis=(1,2,3)/sqrt14:

  | Rays | Integral relative error | Vector first-moment relative error |
  | ---: | ---: | ---: |
  | 80 | 1.876730211e-3 | 5.501144076e-3 |
  | 320 | 2.343712646e-8 | 4.835711673e-8 |
  | 720 | 1.554312234e-15 | 6.600940865e-15 |

  This is an analytic angular-integral comparison, NOT spatial ray-effect,
  ionization-front, frequency, or full simulation convergence.
- `level0/default-parity-final.log`: default arrays/order bitwise unchanged.
  Invalid level3 rejects in both Make and preprocessor (logs retained).
- Native checkpoint smokes at all three compiled levels pass, including
  angular dimension rejection before mutation. The initial fixture invocation
  referenced a nonexistent candidate file; the final runs use the actual
  candidate generated from the reference NML and pass.
- Existing regrid smoke at 320 passes primary/IR refinement/restriction,
  bidirectional coarse/fine flux, subcycling, rollback and cyclic rebind.
  Both halves of its old hard-coded 40/80-ray fixture were generalized.
  Existing HDF5 bad-width case rejects with status1 before reading state.
- Native source deposition plus CUDA sparse transport on an A10 passes at
  320 and 720 directions, photon-number residuals 2.609340e-9 and 3.488253e-9.
  Current `snrt_cuda_kernels.cu` was compiled for sm86; existing shared pool/
  MG/scalar/particle CUDA objects supply its unchanged stream API. This is
  device transport evidence, NOT a multi-GPU or full GPU hydro run. The paired
  signed-energy transport's existing OpenMP-only admission is unchanged.
- Final effective NML: `.rt-angular.TwRcDH/final/physical.nml`; restart:
  `.rt-angular.TwRcDH/final-restart/physical.nml`. MPI2 x OMP2, periodic 4^3,
  noncosmological hydro+stellar/feedback+advective CR+D03/DL01 dust+SNRT;
  no sink/AGN. BPASS independent-population radiation and physical returns
  retain their reference-comparison identity; they are not a matched population.
  H/He RT chemistry plus WSS09 CIE metal cooling, not CHIMES in this build.
- `nstepmax=4`, `noutput=1`, aout2/tout1e30 outside run, foutput2,
  fbackup1000000. Two 29-MB dumps, restart from step2 to step4 with one new
  dump. Earlier evaluated run before the fixed-default-loop correction was
  also retained until final evaluation. Total generated HDF5 below .2 GB.
- `evaluation-final.log`: all **104 non-header datasets** exactly match the
  uninterrupted final result: 12 AMR, 3 coarse, 1 domain, 8 gravity, 60 hydro,
  17 particle, 3 radiation. All numeric datasets finite; SNRT attributes and
  final time/step attributes match. Step4, code time .19419469. The 64 active
  leaves have nonnegative primary/IR populations and nonzero totals: primary
  photon CODE-density sum 1.8418265577526385, IR-density sum
  1.6473918627786035e-14 (not volume-integrated universe counts/energies).
  Minimum density .0009998035582764677 and gas internal energy
  1.2031096424983634e-7 in code units, after kinetic/CR subtraction.
  Radiation cell width=49284 includes 320 directions and 136 IR nodes.

Driver verdict: angular refinement is available through the real runtime;
default stays 80 because memory/work scale by 4/9, not because higher
resolution lacks value. Select resolution with the intended science case.
Group8 spectral reconstruction/coefficient evolution remains unfinished.

## Evaluated raw output cleanup manifest

All paths below are relative to `.rt-angular.TwRcDH/`. Retain effective NMLs,
environment, build/test logs, binaries, headers and compact evaluation.
Resolved files are evaluated, not production archives or other workers' jobs.

| File | Bytes | SHA256 |
| --- | ---: | --- |
| live-fixed/output_00001/data_00001.h5 | 29012744 | 6d04ebc540f36c70e9097468d4f8a8f4caa2c9b769da3daeb094dd406b8722b4 |
| live-fixed/output_00002/data_00002.h5 | 29019912 | 007056b2781f1c20f7750ea68980a7e2b1b0f93a529199f555d563a0f06689b9 |
| restart/output_00002/data_00002.h5 | 29019912 | 785224fe6019f9009554c4f63708ee90d18330d2a4c74688657bc8bc6f3dbdc4 |
| final/output_00001/data_00001.h5 | 29012744 | 75d18da8ac298b73c388d041eb439e0d0ff02d396207b519fc23f4a78e385310 |
| final/output_00002/data_00002.h5 | 29019912 | 786cb41d7d85bb8f83743ff1c5414482aa0c8aeddcbfc55c5a3a207435085e2a |
| final-restart/output_00002/data_00002.h5 | 29019912 | 130e2c0cb42d48064299292c9ff0bb3a9c5aa722c5f94185eceafece4cf8822a |
| bad-width/radiation.h5 | 834552 | c1e998f0b0a1c432e604fa3691abbf64b33653e8eb72def13a93fd12fc96e924 |

Cleanup confirmed: all seven exact files permanently removed, 174,939,688
bytes freed; no HDF5 files remain in this bundle directory. Recovery requires
rerunning the retained inputs. No unrelated directory, raw output, source
data, checkpoint still needed for an active comparison, or log was removed.

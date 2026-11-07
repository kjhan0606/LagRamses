# Native leaf/subface QUMOND operator implementation

Source 14cdb78 adds qumond_composite.f90; 285cae1 supplies explicit Fortran module
dependencies. The module is linked to the current native RAMSES build. It walks
the actual serial periodic octree, enumerates owned leaf cells and positive
same-level/coarse–fine subfaces, and checks that each side of every leaf receives
its full physical face area. It does not construct a global auxiliary uniform
mesh or read offline connectivity.

Each subface has one oriented flux. Coarse tangential potential reconstruction
uses a fixed geometry-weighted least-squares slope, linear-exact on complete
stencils. Its normal difference uses (H+h)/2. composite_apply accumulates that
shared scalar-gradient flux as a Laplacian; composite_phantom uses the same face
vector and the simple nu family to form (nu-1)*g_N and accumulates its negative
divergence/fourpi. Both adjacent leaf cells receive identical/opposite integrated
contributions. The reconstructed interface operator is not assumed symmetric;
no CG solver or composite elliptic convergence is claimed.

The opt-in RAR_QUMOND_COMPOSITE_OPERATOR check runs after native stage-history
capture of the finest Newtonian field. It evaluates the native operators and
checks constants, interior affine stencils including coarse–fine subfaces, and
pre-mean-removal volume-integrated cancellation. It does not overwrite rho,
phi or f or publish new forces. This is the native spatial building block for
the reviewed synchronous composite operation, not another source prescription
substituted into the existing sequential Poisson solves.

## Compiled/runtime evidence

First build 1275357 failed after 11s because poisson_commons was pulled in through
a prerequisite before the new module was compiled. Preserve that failure and its
dependency-pending jobs1275358/1275359. Corrected build1275360 at285cae1 completed
in2m53s. Matched old_off/new_off/new_on execution1275361 completed in36s, each
with two snapshots. Native operator audit1275396 completed with exit0 and
history/default-off audit1275362 also completed. Their result manifests verify.

At all three captured epochs the native hierarchy has 425832 leaves and1286316
unique oriented subfaces. Saved AMR son flags independently give425832 leaves
and total leaf volume1 in both snapshots. There are70560 affine-interface
component samples. Constant-field error is exactly0; maximum affine vector error
is4.440892098500626e-16. Relative volume-weighted Laplacian cancellation residual
is1.079957708922888e-14; phantom-source cancellation residual is
2.753651518253353e-15, both measured before any mean removal. The16 old/new-off
saved gravity arrays remain bitwise identical.

## Open implementation and science requirements

No composite elliptic iteration, composite convergence criterion or synchronous
hierarchy source/force publication is implemented. Existing per-level MG and its
skipped phantom-source layer still determine the simulation's forces. The next
bundle must solve both Poisson problems with these composite residuals and
publish the full hierarchy before any consumer; native MG/FFTW may precondition
the nonsymmetric outer iteration. Smooth-source resolution convergence, particle
force consistency and independent force references remain required. MPI,
regrid/restart, common evolving timesteps and cosmological ICs remain outside
this first serial static two-level operator implementation. Simple-nu tests do
not certify other interpolation families. No galaxy RAR result is established.

Q-GOAL first: does this deliver needed simulation-ready RT/feedback/dust
capability? Q-LEAN second: excessive instrumentation/gates? Answer in that
order, then a concise verdict. Read-only; do not edit files or launch jobs.

User says proceed to next task; all medium-term groups preapproved. Group6
collision comparison done; PAH C fragmentation has recorded physical-input
gaps. Next deliverable is the angular component of approved group8 (NOT
claiming intragroup spectral reconstruction or all remaining groups done).

Existing production is S_N, 8 Gauss-Legendre mu x 10 uniform phi =80 rays,
constant in snrt_state and quadrature. Backend transport/scatter/dust kernels
accept dynamic direction counts; native and HDF5 checkpoint widths depend on
that count. Keep exact existing 80-ray node/weight/order default.

Implement compile-time SNRT_ANGULAR_LEVEL=0/1/2 -> 80/320/720 directions,
unique 8x10,16x20,24x30 product rules. No arbitrary same-count factorizations,
adaptive ray splitting, new radiation physics or NML knobs. State and
quadrature share constants. Guard unsupported level at Make and preprocessor.
Document clean separate builds to avoid stale module dimensions; preserve
VPATH order. Check checkpoint contracts reject mismatched direction counts
before payload commit; fixed unique family makes dimensions sufficient for
this change. Add missing identity checks only if needed, not a new framework.

Existing angular smoke extended: norm/positive weights, first/second moments,
antipodal pairing, analytic oriented narrow angular lobe moments at three
levels; quantify quadrature error/cost without calling this spatial ray-effect
or full RT convergence. Existing source/transport/checkpoint native smokes at
refined level, small integrated MPI2 hydro+SNRT(+existing dust if economical)
fresh/restart, minimal dumps, remove raw outputs once evaluated. CPU/GPU
dimension handling inspected and existing backend evidence extended where
device available; no GPU qualification claimed if device unavailable.

Fable once at plan; driver end evaluation per latest operator directive.
No repeated end-audit cascade or new analysis infrastructure. Runtime gain:
user can run refined discrete-ordinates angular transport through same driver.
Defaults/source spectra/chemistry coefficients unchanged. Numerical reference:
NIST DLMF3.5 Gauss-Legendre quadrature. No new physical data required.

# P0 TPU Memory Model

Status: Draft v0.1  
Date: 2026-08-28

## 1. Baseline assumptions

- Precision: float32, 4 bytes per scalar.
- Photon groups: G = 9.
- Direction batch: B_omega = 8.
- Stored radiation moments: photon number plus three flux components.
- Gas, chemistry, opacity, source, mask, and diagnostic state: 16 float32
  fields per cell before halos and temporary work arrays.
- This model is for a single globally assembled Cartesian cube. Nested blocks
  have the same total-cell accounting plus halo overhead.

The direction-batched intensity storage is:

    M_intensity = N_cell * G * B_omega * 4 bytes.

The moment storage is:

    M_moments = N_cell * G * 4 * 4 bytes.

The primary gas-state storage is:

    M_gas = N_cell * 16 * 4 bytes.

## 2. Memory table

| Cube | Cells | Intensity batch | Moments | Gas state | Explicit working target | Implicit working target |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 128^3 | 2,097,152 | 0.56 GiB | 0.28 GiB | 0.13 GiB | 2 GiB | 3 GiB |
| 256^3 | 16,777,216 | 4.50 GiB | 2.25 GiB | 1.00 GiB | 12 GiB | 24 GiB |
| 512^3 | 134,217,728 | 36.00 GiB | 18.00 GiB | 8.00 GiB | 96 GiB | 192 GiB |

The explicit target includes halos, opacity and chemistry scratch arrays,
source fields, flux temporaries, and a 1.5 safety factor. The implicit target
also reserves a second intensity batch for the Jacobi iterate and additional
diffusion-acceleration scratch arrays.

S4, S6, and S8 have identical peak intensity memory under eight-direction
batching. Their transport work scales with the number of batches: 3, 6, and
10, respectively. S6 is therefore the reference run because it doubles S4
transport work without doubling peak radiation memory.

## 3. TPU execution implications

- The minimum useful development configuration is one 128^3 cube.
- The first real snapshot target is one 256^3 static domain, sharded over TPU
  cores by spatial blocks.
- A 512^3 domain is not an initial P1 target. It requires a pod-scale run,
  proven halo exchange, and a measured rather than estimated memory budget.
- Radiation intensity is direction-batched. Photon moments and source fields
  remain resident across batches.
- Static compile shapes are required. Physical subdomains smaller than the
  allocated cube use masks; no run-time reallocation is permitted.

## 4. Mandatory measurements in P1

The first TPU prototype must measure, rather than infer:

1. compilation time for S4 and S6;
2. per-step transport time, chemistry time, and halo-exchange time;
3. high-bandwidth-memory peak allocation;
4. scaling from one device to the first multi-device mesh; and
5. the cost of X-ray chemistry relative to HI-only chemistry.

P1 cannot promote 256^3 to the science workflow until the measured peak HBM
is below 70 percent of the available aggregate HBM, leaving margin for JAX/XLA
compiler buffers and restart I/O.

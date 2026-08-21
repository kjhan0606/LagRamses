// ==========================================================================
// GPU CIC force-gather + kick/drift for move_fine / synchro_fine
// (strategy-A dynamic work sharing: each OpenMP thread that acquires a
//  stream from the shared pool becomes a GPU worker; the others run the
//  unchanged CPU path. Batches are independent, so schedule(dynamic)
//  balances the split at run time.)
//
// Faithful port of the particle arithmetic in move1/sync
// (force_fine-side tree work — get3cubefather, x0 — stays on the CPU
//  and arrives packed per superbatch):
//  - rescale to grid-cube units, 0.5..5.5 bounds check,
//  - CIC at ilevel; fall back to ilevel-1 CIC when any of the 8 cloud
//    grids is missing (ok flag), exactly as the CPU code,
//  - gather f at the 8 parent cells, kick (mode sync/move) and drift
//    (mode move), optional particle-potential store.
// Sink and tracer particles never take this path (gated in Fortran).
// ==========================================================================

#include "cuda_stream_pool.h"
#include "amr_cuda_index.cuh"
#include <atomic>
#include <cstdio>
#include <ctime>

// Accumulated cost split, used to locate the CPU/GPU crossover:
//   npart_crit = t_upload_per_call / (t_cpu_per_part - t_gpu_per_part)
static double    g_pm_t_upload = 0.0;   // seconds in mesh/son uploads
static double    g_pm_t_flush  = 0.0;   // seconds in H2D + kernel + D2H
static long long g_pm_n_upload = 0;     // upload calls
static long long g_pm_n_flush  = 0;     // flush calls
static long long g_pm_n_part   = 0;     // particles processed on the GPU
static long long g_pm_up_bytes = 0;     // bytes uploaded
static long long g_pm_mesh_uploads = 0;
static std::atomic<long long> g_pm_gather_launches{0};
static long long g_pm_rho_uploads = 0;
static std::atomic<long long> g_pm_deposit_launches{0};
static std::atomic<long long> g_pm_deposit_parts{0};
static int g_pm_block_size = 0;
static int g_pm_child_count = 0;

static bool pm_layout_valid(long long ncell, long long ncoarse,
                            int ngridmax, int block_size, int child_count)
{
    return amr_cuda_layout_valid(
        ncell, ncoarse, ngridmax, block_size, child_count, (1 << NDIM));
}

static double pm_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#define PM_MODE_MOVE 0
#define PM_MODE_SYNC 1

typedef struct {
    double scale, dx;
    double skip0, skip1, skip2;
    double cde_coeff;   // beta_cde*phip*hexp*0.5 (0 when off)
    double dt_move;     // dtnew(ilevel) for the drift (mode move)
    int    is_static;
    int    store_phi;
    int    mode;
} PmParams;

// --------------------------------------------------------------------------
// Mesh state (uploaded once per move/synchro call)
// --------------------------------------------------------------------------
static double* d_pm_f    = nullptr;   // f(1:ncell,1:3) column-major
static double* d_pm_phi  = nullptr;
static int*    d_pm_son  = nullptr;
static long long g_pm_ncell = 0;
static long long g_pm_ncell_cap = 0;
static bool g_pm_has_phi = false;
static bool g_pm_ready = false;

// --------------------------------------------------------------------------
// Per-slot staging buffers (slot = stream-pool slot index)
// --------------------------------------------------------------------------
typedef struct {
    // device
    double *d_x0;      // (3, ngcap)
    int    *d_nbf;     // (27, ngcap)
    double *d_px, *d_pv;      // (3, npcap)
    int    *d_pg;             // (npcap) 0-based superbatch grid index
    double *d_dteff;          // (npcap)
    double *d_new_v, *d_new_x;// (3, npcap)
    double *d_phi_out;        // (npcap)
    int    *d_err;
    // pinned host mirrors
    double *h_x0;
    int    *h_nbf;
    double *h_px, *h_pv;
    int    *h_pg;
    double *h_dteff;
    double *h_new_v, *h_new_x;
    double *h_phi_out;
    int ngcap, npcap;
} PmSlot;

static PmSlot g_pm_slot[MAX_CUDA_STREAMS];
static bool g_pm_slot_init = false;

static bool pm_ensure_slot(int slot, int ng, int np)
{
    if (!g_pm_slot_init) {
        for (int s = 0; s < MAX_CUDA_STREAMS; s++) {
            PmSlot* p = &g_pm_slot[s];
            p->d_x0=nullptr; p->d_nbf=nullptr; p->d_px=nullptr; p->d_pv=nullptr;
            p->d_pg=nullptr; p->d_dteff=nullptr; p->d_new_v=nullptr;
            p->d_new_x=nullptr; p->d_phi_out=nullptr; p->d_err=nullptr;
            p->h_x0=nullptr; p->h_nbf=nullptr; p->h_px=nullptr; p->h_pv=nullptr;
            p->h_pg=nullptr; p->h_dteff=nullptr; p->h_new_v=nullptr;
            p->h_new_x=nullptr; p->h_phi_out=nullptr;
            p->ngcap=0; p->npcap=0;
        }
        g_pm_slot_init = true;
    }
    PmSlot* p = &g_pm_slot[slot];
    if (ng > p->ngcap) {
        int cap = ng * 2;
        if (p->d_x0)  cudaFree(p->d_x0);
        if (p->d_nbf) cudaFree(p->d_nbf);
        if (p->h_x0)  cudaFreeHost(p->h_x0);
        if (p->h_nbf) cudaFreeHost(p->h_nbf);
        cudaError_t e1 = cudaMalloc(&p->d_x0,  (size_t)cap * 3 * sizeof(double));
        cudaError_t e2 = cudaMalloc(&p->d_nbf, (size_t)cap * 27 * sizeof(int));
        cudaError_t e3 = cudaMallocHost(&p->h_x0,  (size_t)cap * 3 * sizeof(double));
        cudaError_t e4 = cudaMallocHost(&p->h_nbf, (size_t)cap * 27 * sizeof(int));
        if (e1 || e2 || e3 || e4) { p->ngcap = 0; return false; }
        p->ngcap = cap;
    }
    if (np > p->npcap) {
        int cap = np * 2;
        if (p->d_px)      cudaFree(p->d_px);
        if (p->d_pv)      cudaFree(p->d_pv);
        if (p->d_pg)      cudaFree(p->d_pg);
        if (p->d_dteff)   cudaFree(p->d_dteff);
        if (p->d_new_v)   cudaFree(p->d_new_v);
        if (p->d_new_x)   cudaFree(p->d_new_x);
        if (p->d_phi_out) cudaFree(p->d_phi_out);
        if (p->h_px)      cudaFreeHost(p->h_px);
        if (p->h_pv)      cudaFreeHost(p->h_pv);
        if (p->h_pg)      cudaFreeHost(p->h_pg);
        if (p->h_dteff)   cudaFreeHost(p->h_dteff);
        if (p->h_new_v)   cudaFreeHost(p->h_new_v);
        if (p->h_new_x)   cudaFreeHost(p->h_new_x);
        if (p->h_phi_out) cudaFreeHost(p->h_phi_out);
        cudaError_t e = cudaSuccess;
        if (!e) e = cudaMalloc(&p->d_px,      (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMalloc(&p->d_pv,      (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMalloc(&p->d_pg,      (size_t)cap * sizeof(int));
        if (!e) e = cudaMalloc(&p->d_dteff,   (size_t)cap * sizeof(double));
        if (!e) e = cudaMalloc(&p->d_new_v,   (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMalloc(&p->d_new_x,   (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMalloc(&p->d_phi_out, (size_t)cap * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_px,      (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_pv,      (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_pg,      (size_t)cap * sizeof(int));
        if (!e) e = cudaMallocHost(&p->h_dteff,   (size_t)cap * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_new_v,   (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_new_x,   (size_t)cap * 3 * sizeof(double));
        if (!e) e = cudaMallocHost(&p->h_phi_out, (size_t)cap * sizeof(double));
        if (e) { p->npcap = 0; return false; }
        p->npcap = cap;
    }
    if (!p->d_err) cudaMalloc(&p->d_err, sizeof(int));
    return true;
}

// --------------------------------------------------------------------------
// The CIC kernel: one thread per particle
// --------------------------------------------------------------------------
__global__ void pm_cic_kernel(
    const double* __restrict__ px,
    const double* __restrict__ pv,
    const int*    __restrict__ pg,
    const double* __restrict__ dteff,
    const double* __restrict__ gx0,
    const int*    __restrict__ gnbf,
    double* __restrict__ new_v,
    double* __restrict__ new_x,
    double* __restrict__ phi_out,
    const double* __restrict__ mesh_f,
    const double* __restrict__ mesh_phi,
    const int*    __restrict__ mesh_son,
    long long ncell, int ngridmax, int ncoarse,
    int block_size, int child_count,
    int np, PmParams pp,
    int* __restrict__ err_flag)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= np) return;

    int g = pg[j]; // 0-based superbatch grid index
    const double skip[3] = {pp.skip0, pp.skip1, pp.skip2};

    double x[3];
    for (int d = 0; d < 3; d++) {
        double xx = px[j*3+d] / pp.scale + skip[d];
        xx = xx - gx0[g*3+d];
        x[d] = xx / pp.dx;
        if (x[d] < 0.5 || x[d] > 5.5) { atomicExch(err_flag, 1); return; }
    }

    // CIC at level ilevel
    double dd[3], dg[3];
    int id[3], ig[3], igg[3], igd[3];
    for (int d = 0; d < 3; d++) {
        dd[d] = x[d] + 0.5;
        id[d] = (int)dd[d];
        dd[d] = dd[d] - id[d];
        dg[d] = 1.0 - dd[d];
        ig[d] = id[d] - 1;
        igg[d] = ig[d] / 2;
        igd[d] = id[d] / 2;
    }

    int kg[8];
    kg[0] = 1 + igg[0] + 3*igg[1] + 9*igg[2];
    kg[1] = 1 + igd[0] + 3*igg[1] + 9*igg[2];
    kg[2] = 1 + igg[0] + 3*igd[1] + 9*igg[2];
    kg[3] = 1 + igd[0] + 3*igd[1] + 9*igg[2];
    kg[4] = 1 + igg[0] + 3*igg[1] + 9*igd[2];
    kg[5] = 1 + igd[0] + 3*igg[1] + 9*igd[2];
    kg[6] = 1 + igg[0] + 3*igd[1] + 9*igd[2];
    kg[7] = 1 + igd[0] + 3*igd[1] + 9*igd[2];

    int igr[8];
    bool ok = true;
    for (int ind = 0; ind < 8; ind++) {
        int fc = gnbf[g*27 + kg[ind] - 1];        // 1-based father cell
        igr[ind] = mesh_son[fc - 1];              // 1-based grid or 0
        if (igr[ind] <= 0) ok = false;
    }

    if (!ok) {
        // redo CIC at level ilevel-1
        for (int d = 0; d < 3; d++) {
            x[d] = x[d] / 2.0;
            dd[d] = x[d] + 0.5;
            id[d] = (int)dd[d];
            dd[d] = dd[d] - id[d];
            dg[d] = 1.0 - dd[d];
            ig[d] = id[d] - 1;
        }
    }

    int icg[3], icd[3];
    for (int d = 0; d < 3; d++) {
        if (ok) {
            icg[d] = ig[d] - 2*igg[d];
            icd[d] = id[d] - 2*igd[d];
        } else {
            icg[d] = ig[d];
            icd[d] = id[d];
        }
    }

    int icell[8];
    if (ok) {
        icell[0] = 1 + icg[0] + 2*icg[1] + 4*icg[2];
        icell[1] = 1 + icd[0] + 2*icg[1] + 4*icg[2];
        icell[2] = 1 + icg[0] + 2*icd[1] + 4*icg[2];
        icell[3] = 1 + icd[0] + 2*icd[1] + 4*icg[2];
        icell[4] = 1 + icg[0] + 2*icg[1] + 4*icd[2];
        icell[5] = 1 + icd[0] + 2*icg[1] + 4*icd[2];
        icell[6] = 1 + icg[0] + 2*icd[1] + 4*icd[2];
        icell[7] = 1 + icd[0] + 2*icd[1] + 4*icd[2];
    } else {
        icell[0] = 1 + icg[0] + 3*icg[1] + 9*icg[2];
        icell[1] = 1 + icd[0] + 3*icg[1] + 9*icg[2];
        icell[2] = 1 + icg[0] + 3*icd[1] + 9*icg[2];
        icell[3] = 1 + icd[0] + 3*icd[1] + 9*icg[2];
        icell[4] = 1 + icg[0] + 3*icg[1] + 9*icd[2];
        icell[5] = 1 + icd[0] + 3*icg[1] + 9*icd[2];
        icell[6] = 1 + icg[0] + 3*icd[1] + 9*icd[2];
        icell[7] = 1 + icd[0] + 3*icd[1] + 9*icd[2];
    }

    long long indp[8];
    for (int ind = 0; ind < 8; ind++) {
        if (ok) indp[ind] = amr_cuda_cell_1based(
            igr[ind], icell[ind], ncoarse, block_size, child_count);
        else    indp[ind] = gnbf[g*27 + icell[ind] - 1];
    }

    double vol[8];
    vol[0] = dg[0]*dg[1]*dg[2];
    vol[1] = dd[0]*dg[1]*dg[2];
    vol[2] = dg[0]*dd[1]*dg[2];
    vol[3] = dd[0]*dd[1]*dg[2];
    vol[4] = dg[0]*dg[1]*dd[2];
    vol[5] = dd[0]*dg[1]*dd[2];
    vol[6] = dg[0]*dd[1]*dd[2];
    vol[7] = dd[0]*dd[1]*dd[2];

    // Gather 3-force (same accumulation order as the CPU: ind outer, dim inner)
    double ff[3] = {0.0, 0.0, 0.0};
    for (int ind = 0; ind < 8; ind++) {
        long long c = (long long)(indp[ind] - 1);
        for (int d = 0; d < 3; d++) {
            ff[d] += mesh_f[(long long)d * ncell + c] * vol[ind];
        }
    }
    if (pp.store_phi && phi_out && mesh_phi) {
        phi_out[j] = mesh_phi[(long long)(indp[7] - 1)];
    }

    // Kick
    double dt = dteff[j];
    double fric = (pp.cde_coeff != 0.0) ? exp(pp.cde_coeff * dt) : 1.0;
    double nv[3];
    for (int d = 0; d < 3; d++) {
        if (pp.is_static) nv[d] = ff[d];
        else              nv[d] = pv[j*3+d] * fric + ff[d] * 0.5 * dt;
        new_v[j*3+d] = nv[d];
    }

    // Drift (move mode only)
    if (pp.mode == PM_MODE_MOVE) {
        for (int d = 0; d < 3; d++) {
            if (pp.is_static) new_x[j*3+d] = px[j*3+d];
            else              new_x[j*3+d] = px[j*3+d] + nv[d] * pp.dt_move;
        }
    }
}

// --------------------------------------------------------------------------
// C API
// --------------------------------------------------------------------------
extern "C" {

// Upload only the cells that exist: the coarse block plus, for each of the
// twotondim oct slots, the first `hw` grids (hw = highest grid index in use).
// ngridmax is an allocation ceiling that is typically ~50x the grids actually
// present, so copying the whole array wastes almost all of the transfer.
static bool pm_upload_sliced(void* dst, const void* src, size_t elemsz,
                             long long ncoarse, int ngridmax, int hw,
                             int block_size, int child_count,
                             int ncomp, long long ncell)
{
    // [RESIZABLE] Block grid-major storage makes every complete block through
    // the high-water grid a single contiguous prefix.  Copy each component's
    // prefix using the full allocated component stride.
    if (!pm_layout_valid(
            ncell, ncoarse, ngridmax, block_size, child_count) ||
        hw < 0 || hw > ngridmax)
        return false;
    const long long live_ncell = amr_cuda_live_cell_prefix(
        ncoarse, hw, block_size, child_count);
    if (live_ncell < ncoarse || live_ncell > ncell) return false;
    g_pm_up_bytes += (long long)ncomp * live_ncell * (long long)elemsz;
    for (int c = 0; c < ncomp; c++) {
        const char* hp = (const char*)src + (size_t)c * ncell * elemsz;
        char*       dp = (char*)dst       + (size_t)c * ncell * elemsz;
        if (live_ncell > 0 && cudaMemcpy(dp, hp, (size_t)live_ncell * elemsz,
                                        cudaMemcpyHostToDevice) != cudaSuccess)
            return false;
    }
    return true;
}

void cuda_pm_mesh_upload(const double* f, const int* son, const double* phi,
                         long long ncell, int with_phi,
                         long long ncoarse, int ngridmax, int hw,
                         int block_size, int child_count)
{
    g_pm_ready = false;
    if (!is_pool_initialized()) return;
    cudaGetLastError();

    if (ncell > g_pm_ncell_cap) {
        if (d_pm_f)   { cudaFree(d_pm_f);   d_pm_f   = nullptr; }
        if (d_pm_son) { cudaFree(d_pm_son); d_pm_son = nullptr; }
        if (d_pm_phi) { cudaFree(d_pm_phi); d_pm_phi = nullptr; }
        size_t free_mem = 0, total_mem = 0;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t need = (size_t)ncell * (4 * sizeof(double) + sizeof(int));
        if (need > free_mem * 9 / 10) {
            fprintf(stderr, "CUDA pm mesh: SKIP — need %.1f GB, %.1f GB free\n",
                    (double)need / 1073741824.0, (double)free_mem / 1073741824.0);
            g_pm_ncell_cap = 0;
            return;
        }
        cudaError_t e1 = cudaMalloc(&d_pm_f,   (size_t)ncell * 3 * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_pm_son, (size_t)ncell * sizeof(int));
        cudaError_t e3 = cudaMalloc(&d_pm_phi, (size_t)ncell * sizeof(double));
        if (e1 || e2 || e3) {
            fprintf(stderr, "CUDA pm mesh: allocation FAILED\n");
            if (d_pm_f)   { cudaFree(d_pm_f);   d_pm_f   = nullptr; }
            if (d_pm_son) { cudaFree(d_pm_son); d_pm_son = nullptr; }
            if (d_pm_phi) { cudaFree(d_pm_phi); d_pm_phi = nullptr; }
            g_pm_ncell_cap = 0;
            return;
        }
        g_pm_ncell_cap = ncell;
    }
    g_pm_ncell = ncell;

    double t0 = pm_now();
    bool upload_ok =
        pm_upload_sliced(d_pm_f, f, sizeof(double), ncoarse, ngridmax, hw,
                         block_size, child_count, 3, ncell) &&
        pm_upload_sliced(d_pm_son, son, sizeof(int), ncoarse, ngridmax, hw,
                         block_size, child_count, 1, ncell);
    g_pm_has_phi = (with_phi != 0);
    if (g_pm_has_phi && phi)
        upload_ok = upload_ok &&
            pm_upload_sliced(d_pm_phi, phi, sizeof(double), ncoarse, ngridmax, hw,
                             block_size, child_count, 1, ncell);
    const cudaError_t sync_error = cudaDeviceSynchronize();
    g_pm_t_upload += pm_now() - t0;
    g_pm_n_upload++;

    const cudaError_t launch_error = cudaGetLastError();
    if (!upload_ok) {
        fprintf(stderr, "CUDA pm mesh upload error: invalid block layout or copy\n");
        return;
    }
    const cudaError_t err =
        sync_error != cudaSuccess ? sync_error : launch_error;
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA pm mesh upload error: %s\n", cudaGetErrorString(err));
        return;
    }
    g_pm_ready = true;
    g_pm_mesh_uploads++;
    g_pm_block_size = block_size;
    g_pm_child_count = child_count;
}

int cuda_pm_is_ready(void) { return g_pm_ready ? 1 : 0; }

// Synchronous superbatch flush on the slot's stream:
// H2D -> kernel -> D2H. Returns 0 on success, 1 on a bounds error
// ("problem in move/sync"), -1 on CUDA failure.
int cuda_pm_flush(int slot, int ng, int np,
                  const double* h_x0, const int* h_nbf,
                  const double* h_px, const double* h_pv,
                  const int* h_pg, const double* h_dteff,
                  const double* params, int ngridmax, int ncoarse,
                  int block_size, int child_count,
                  double* out_new_v, double* out_new_x, double* out_phi)
{
    if (!g_pm_ready || slot < 0 || np <= 0 || ng <= 0) return -1;
    if (!pm_layout_valid(
            g_pm_ncell, ncoarse, ngridmax, block_size, child_count) ||
        block_size != g_pm_block_size || child_count != g_pm_child_count)
        return -1;
    if (!pm_ensure_slot(slot, ng, np)) return -1;

    double t_f0 = pm_now();
    PmSlot* p = &g_pm_slot[slot];
    cudaStream_t st = cuda_get_stream_internal(slot);

    PmParams pp;
    pp.scale     = params[0];
    pp.dx        = params[1];
    pp.skip0     = params[2];
    pp.skip1     = params[3];
    pp.skip2     = params[4];
    pp.cde_coeff = params[5];
    pp.dt_move   = params[6];
    pp.is_static = (int)params[7];
    pp.store_phi = (int)params[8] && g_pm_has_phi;
    pp.mode      = (int)params[9];

    // Stage into pinned mirrors
    memcpy(p->h_x0,  h_x0,  (size_t)ng * 3 * sizeof(double));
    memcpy(p->h_nbf, h_nbf, (size_t)ng * 27 * sizeof(int));
    memcpy(p->h_px,  h_px,  (size_t)np * 3 * sizeof(double));
    memcpy(p->h_pv,  h_pv,  (size_t)np * 3 * sizeof(double));
    memcpy(p->h_pg,  h_pg,  (size_t)np * sizeof(int));
    memcpy(p->h_dteff, h_dteff, (size_t)np * sizeof(double));

    cudaMemsetAsync(p->d_err, 0, sizeof(int), st);
    cudaMemcpyAsync(p->d_x0,  p->h_x0,  (size_t)ng * 3 * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_nbf, p->h_nbf, (size_t)ng * 27 * sizeof(int),   cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_px,  p->h_px,  (size_t)np * 3 * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_pv,  p->h_pv,  (size_t)np * 3 * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_pg,  p->h_pg,  (size_t)np * sizeof(int),        cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_dteff, p->h_dteff, (size_t)np * sizeof(double), cudaMemcpyHostToDevice, st);

    int block = 128;
    int grid = (np + block - 1) / block;
    pm_cic_kernel<<<grid, block, 0, st>>>(
        p->d_px, p->d_pv, p->d_pg, p->d_dteff,
        p->d_x0, p->d_nbf,
        p->d_new_v, p->d_new_x, p->d_phi_out,
        d_pm_f, d_pm_phi, d_pm_son,
        g_pm_ncell, ngridmax, ncoarse, block_size, child_count,
        np, pp, p->d_err);

    cudaMemcpyAsync(p->h_new_v, p->d_new_v, (size_t)np * 3 * sizeof(double), cudaMemcpyDeviceToHost, st);
    if (pp.mode == PM_MODE_MOVE)
        cudaMemcpyAsync(p->h_new_x, p->d_new_x, (size_t)np * 3 * sizeof(double), cudaMemcpyDeviceToHost, st);
    if (pp.store_phi)
        cudaMemcpyAsync(p->h_phi_out, p->d_phi_out, (size_t)np * sizeof(double), cudaMemcpyDeviceToHost, st);

    int h_err = 0;
    cudaMemcpyAsync(&h_err, p->d_err, sizeof(int), cudaMemcpyDeviceToHost, st);
    const cudaError_t sync_error = cudaStreamSynchronize(st);

    const cudaError_t launch_error = cudaGetLastError();
    const cudaError_t cerr =
        sync_error != cudaSuccess ? sync_error : launch_error;
    if (cerr != cudaSuccess) {
        fprintf(stderr, "CUDA pm flush error: %s\n", cudaGetErrorString(cerr));
        return -1;
    }
    if (h_err) return 1;

    memcpy(out_new_v, p->h_new_v, (size_t)np * 3 * sizeof(double));
    if (pp.mode == PM_MODE_MOVE)
        memcpy(out_new_x, p->h_new_x, (size_t)np * 3 * sizeof(double));
    if (pp.store_phi)
        memcpy(out_phi, p->h_phi_out, (size_t)np * sizeof(double));
    g_pm_t_flush += pm_now() - t_f0;
    g_pm_n_flush++;
    g_pm_n_part += np;
    const long long gather_count =
        g_pm_gather_launches.fetch_add(1, std::memory_order_relaxed) + 1;
    if (gather_count == 1) {
        printf("[CUDA_PM_GATHER] B=%d C=%d mesh_upload=%lld gather=%lld particles=%lld\n",
               block_size, child_count, g_pm_mesh_uploads,
               gather_count, (long long)np);
        fflush(stdout);
    }
    return 0;
}

// Cost split for tuning pm_gpu_min_part. The crossover follows from
//   npart_crit = t_upload_per_call / (t_cpu_per_part - t_gpu_per_part)
// where t_cpu_per_part comes from the same run with the GPU path disabled.
void cuda_pm_report(void)
{
    if (g_pm_n_upload == 0 && g_pm_n_flush == 0) return;
    printf(" === GPU particle path cost split ===\n");
    printf("   uploads    : %lld calls, %.3f s total, %.3f ms/call, %.1f MB/call\n",
           g_pm_n_upload, g_pm_t_upload,
           g_pm_n_upload ? 1e3 * g_pm_t_upload / (double)g_pm_n_upload : 0.0,
           g_pm_n_upload ? (double)g_pm_up_bytes / (double)g_pm_n_upload / 1048576.0 : 0.0);
    printf("   flushes    : %lld calls, %.3f s total, %lld particles\n",
           g_pm_n_flush, g_pm_t_flush, g_pm_n_part);
    printf("   per particle: %.1f ns on the GPU path\n",
           g_pm_n_part ? 1e9 * g_pm_t_flush / (double)g_pm_n_part : 0.0);
    printf("   upload amortises once npart/call exceeds t_up/(t_cpu-t_gpu)\n");
    printf("[CUDA_PM] B=%d C=%d mesh_upload=%lld gather=%lld "
           "rho_upload=%lld deposit=%lld gather_particles=%lld "
           "deposit_particles=%lld\n",
           g_pm_block_size, g_pm_child_count, g_pm_mesh_uploads,
           g_pm_gather_launches.load(std::memory_order_relaxed),
           g_pm_rho_uploads,
           g_pm_deposit_launches.load(std::memory_order_relaxed),
           g_pm_n_part,
           g_pm_deposit_parts.load(std::memory_order_relaxed));
    fflush(stdout);
}

// --------------------------------------------------------------------------
// CIC deposit (rho_fine / cic_amr): particle -> mesh via atomicAdd.
// No coarse-level fallback — cells whose cloud grid is missing are
// skipped (ok mask), exactly like the CPU cic_amr. Star/sink runs and
// the cic_levelmax special levels never take this path (Fortran gate),
// so the deposits reduce to rho += m*vol/vol_loc and the number-density
// work array phiw += vol with the static/mass_cut gates.
// --------------------------------------------------------------------------
static double* d_pm_rho  = nullptr;
static double* d_pm_phiw = nullptr;
static long long g_pm_rho_cap = 0;
static long long g_pm_rho_ncell = 0;
static bool g_pm_rho_ready = false;

__global__ void pm_deposit_kernel(
    const double* __restrict__ px,
    const double* __restrict__ pmass,
    const int*    __restrict__ pg,
    const double* __restrict__ gx0,
    const int*    __restrict__ gnbf,
    double* __restrict__ rho_acc,
    double* __restrict__ phiw_acc,
    const int* __restrict__ mesh_son,
    int ngridmax, int ncoarse, int block_size, int child_count,
    int np,
    double scale, double dx, double skip0, double skip1, double skip2,
    double vol_loc, int is_static, double mass_cut,
    int* __restrict__ err_flag)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= np) return;

    int g = pg[j];
    const double skip[3] = {skip0, skip1, skip2};

    double x[3];
    for (int d = 0; d < 3; d++) {
        double xx = px[j*3+d] / scale + skip[d];
        xx = xx - gx0[g*3+d];
        x[d] = xx / dx;
        if (x[d] < 0.5 || x[d] > 5.5) { atomicExch(err_flag, 1); return; }
    }

    double dd[3], dg[3];
    int id[3], ig[3], igg[3], igd[3];
    for (int d = 0; d < 3; d++) {
        dd[d] = x[d] + 0.5;
        id[d] = (int)dd[d];
        dd[d] = dd[d] - id[d];
        dg[d] = 1.0 - dd[d];
        ig[d] = id[d] - 1;
        igg[d] = ig[d] / 2;
        igd[d] = id[d] / 2;
    }

    double vol[8];
    vol[0] = dg[0]*dg[1]*dg[2];
    vol[1] = dd[0]*dg[1]*dg[2];
    vol[2] = dg[0]*dd[1]*dg[2];
    vol[3] = dd[0]*dd[1]*dg[2];
    vol[4] = dg[0]*dg[1]*dd[2];
    vol[5] = dd[0]*dg[1]*dd[2];
    vol[6] = dg[0]*dd[1]*dd[2];
    vol[7] = dd[0]*dd[1]*dd[2];

    int kg[8];
    kg[0] = 1 + igg[0] + 3*igg[1] + 9*igg[2];
    kg[1] = 1 + igd[0] + 3*igg[1] + 9*igg[2];
    kg[2] = 1 + igg[0] + 3*igd[1] + 9*igg[2];
    kg[3] = 1 + igd[0] + 3*igd[1] + 9*igg[2];
    kg[4] = 1 + igg[0] + 3*igg[1] + 9*igd[2];
    kg[5] = 1 + igd[0] + 3*igg[1] + 9*igd[2];
    kg[6] = 1 + igg[0] + 3*igd[1] + 9*igd[2];
    kg[7] = 1 + igd[0] + 3*igd[1] + 9*igd[2];

    int icg[3], icd[3];
    for (int d = 0; d < 3; d++) {
        icg[d] = ig[d] - 2*igg[d];
        icd[d] = id[d] - 2*igd[d];
    }
    int icell[8];
    icell[0] = 1 + icg[0] + 2*icg[1] + 4*icg[2];
    icell[1] = 1 + icd[0] + 2*icg[1] + 4*icg[2];
    icell[2] = 1 + icg[0] + 2*icd[1] + 4*icg[2];
    icell[3] = 1 + icd[0] + 2*icd[1] + 4*icg[2];
    icell[4] = 1 + icg[0] + 2*icg[1] + 4*icd[2];
    icell[5] = 1 + icd[0] + 2*icg[1] + 4*icd[2];
    icell[6] = 1 + icg[0] + 2*icd[1] + 4*icd[2];
    icell[7] = 1 + icd[0] + 2*icd[1] + 4*icd[2];

    double m = pmass[j];
    bool ok2_scalar = true;
    if (is_static && !(m > 0.0)) ok2_scalar = false;
    if (mass_cut > 0.0 && !(m < mass_cut)) ok2_scalar = false;

    for (int ind = 0; ind < 8; ind++) {
        int fc = gnbf[g*27 + kg[ind] - 1];
        int igr = mesh_son[fc - 1];
        if (igr <= 0) continue;
        long long c = (long long)amr_cuda_cell_1based(
            igr, icell[ind], ncoarse, block_size, child_count) - 1;
        atomicAdd(&rho_acc[c], m * vol[ind] / vol_loc);
        if (ok2_scalar) atomicAdd(&phiw_acc[c], vol[ind]);
    }
}

extern "C" {

void cuda_pm_rho_begin(const int* son, long long ncell,
                       long long ncoarse, int ngridmax, int hw,
                       int block_size, int child_count)
{
    g_pm_rho_ready = false;
    if (!is_pool_initialized()) return;
    cudaGetLastError();

    // son goes to the shared d_pm_son (same layout as the move/sync path)
    if (ncell > g_pm_ncell_cap) {
        if (d_pm_f)   { cudaFree(d_pm_f);   d_pm_f   = nullptr; }
        if (d_pm_son) { cudaFree(d_pm_son); d_pm_son = nullptr; }
        if (d_pm_phi) { cudaFree(d_pm_phi); d_pm_phi = nullptr; }
        cudaError_t e1 = cudaMalloc(&d_pm_f,   (size_t)ncell * 3 * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_pm_son, (size_t)ncell * sizeof(int));
        cudaError_t e3 = cudaMalloc(&d_pm_phi, (size_t)ncell * sizeof(double));
        if (e1 || e2 || e3) { g_pm_ncell_cap = 0; return; }
        g_pm_ncell_cap = ncell;
    }
    if (ncell > g_pm_rho_cap) {
        if (d_pm_rho)  { cudaFree(d_pm_rho);  d_pm_rho  = nullptr; }
        if (d_pm_phiw) { cudaFree(d_pm_phiw); d_pm_phiw = nullptr; }
        cudaError_t e1 = cudaMalloc(&d_pm_rho,  (size_t)ncell * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_pm_phiw, (size_t)ncell * sizeof(double));
        if (e1 || e2) {
            fprintf(stderr, "CUDA pm rho: allocation FAILED\n");
            if (d_pm_rho)  { cudaFree(d_pm_rho);  d_pm_rho  = nullptr; }
            if (d_pm_phiw) { cudaFree(d_pm_phiw); d_pm_phiw = nullptr; }
            g_pm_rho_cap = 0;
            return;
        }
        g_pm_rho_cap = ncell;
    }
    g_pm_rho_ncell = ncell;

    double t0 = pm_now();
    bool upload_ok = pm_upload_sliced(
        d_pm_son, son, sizeof(int), ncoarse, ngridmax, hw,
        block_size, child_count, 1, ncell);
    cudaMemset(d_pm_rho,  0, (size_t)ncell * sizeof(double));
    cudaMemset(d_pm_phiw, 0, (size_t)ncell * sizeof(double));
    const cudaError_t sync_error = cudaDeviceSynchronize();
    g_pm_t_upload += pm_now() - t0;
    g_pm_n_upload++;

    const cudaError_t launch_error = cudaGetLastError();
    if (!upload_ok || sync_error != cudaSuccess || launch_error != cudaSuccess)
        return;
    g_pm_rho_ready = true;
    g_pm_rho_uploads++;
    g_pm_block_size = block_size;
    g_pm_child_count = child_count;
}

int cuda_pm_rho_is_ready(void) { return g_pm_rho_ready ? 1 : 0; }

int cuda_pm_deposit_flush(int slot, int ng, int np,
                          const double* h_x0, const int* h_nbf,
                          const double* h_px, const double* h_mass,
                          const int* h_pg,
                          const double* params, int ngridmax, int ncoarse,
                          int block_size, int child_count)
{
    if (!g_pm_rho_ready || slot < 0 || np <= 0 || ng <= 0) return -1;
    if (!pm_layout_valid(
            g_pm_rho_ncell, ncoarse, ngridmax, block_size, child_count) ||
        block_size != g_pm_block_size || child_count != g_pm_child_count)
        return -1;
    if (!pm_ensure_slot(slot, ng, np)) return -1;

    PmSlot* p = &g_pm_slot[slot];
    cudaStream_t st = cuda_get_stream_internal(slot);

    memcpy(p->h_x0,  h_x0,  (size_t)ng * 3 * sizeof(double));
    memcpy(p->h_nbf, h_nbf, (size_t)ng * 27 * sizeof(int));
    memcpy(p->h_px,  h_px,  (size_t)np * 3 * sizeof(double));
    memcpy(p->h_dteff, h_mass, (size_t)np * sizeof(double));
    memcpy(p->h_pg,  h_pg,  (size_t)np * sizeof(int));

    cudaMemsetAsync(p->d_err, 0, sizeof(int), st);
    cudaMemcpyAsync(p->d_x0,  p->h_x0,  (size_t)ng * 3 * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_nbf, p->h_nbf, (size_t)ng * 27 * sizeof(int),   cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_px,  p->h_px,  (size_t)np * 3 * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_dteff, p->h_dteff, (size_t)np * sizeof(double), cudaMemcpyHostToDevice, st);
    cudaMemcpyAsync(p->d_pg,  p->h_pg,  (size_t)np * sizeof(int),        cudaMemcpyHostToDevice, st);

    int block = 128;
    int grid = (np + block - 1) / block;
    pm_deposit_kernel<<<grid, block, 0, st>>>(
        p->d_px, p->d_dteff, p->d_pg, p->d_x0, p->d_nbf,
        d_pm_rho, d_pm_phiw, d_pm_son,
        ngridmax, ncoarse, block_size, child_count, np,
        params[0], params[1], params[2], params[3], params[4],
        params[5], (int)params[6], params[7],
        p->d_err);

    int h_err = 0;
    cudaMemcpyAsync(&h_err, p->d_err, sizeof(int), cudaMemcpyDeviceToHost, st);
    const cudaError_t sync_error = cudaStreamSynchronize(st);

    const cudaError_t launch_error = cudaGetLastError();
    const cudaError_t cerr =
        sync_error != cudaSuccess ? sync_error : launch_error;
    if (cerr != cudaSuccess) {
        fprintf(stderr, "CUDA pm deposit error: %s\n", cudaGetErrorString(cerr));
        return -1;
    }
    if (h_err) return 1;
    const long long deposit_count =
        g_pm_deposit_launches.fetch_add(1, std::memory_order_relaxed) + 1;
    const long long deposit_parts =
        g_pm_deposit_parts.fetch_add(np, std::memory_order_relaxed) + np;
    if (deposit_count == 1) {
        printf("[CUDA_PM_DEPOSIT] B=%d C=%d rho_upload=%lld deposit=%lld particles=%lld\n",
               block_size, child_count, g_pm_rho_uploads,
               deposit_count, deposit_parts);
        fflush(stdout);
    }
    return 0;
}

// Download the accumulated GPU deposits for the host-side merge
void cuda_pm_rho_end(double* rho_add, double* phiw_add, long long ncell)
{
    g_pm_rho_ready = false;
    if (ncell != g_pm_rho_ncell || !d_pm_rho) return;
    cudaMemcpy(rho_add,  d_pm_rho,  (size_t)ncell * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(phiw_add, d_pm_phiw, (size_t)ncell * sizeof(double), cudaMemcpyDeviceToHost);
}

} // extern "C"

extern "C" void cuda_pm_finalize(void)
{
    cuda_pm_report();
    g_pm_ready = false;
    g_pm_rho_ready = false;
    if (d_pm_f)    { cudaFree(d_pm_f);    d_pm_f    = nullptr; }
    if (d_pm_son)  { cudaFree(d_pm_son);  d_pm_son  = nullptr; }
    if (d_pm_phi)  { cudaFree(d_pm_phi);  d_pm_phi  = nullptr; }
    if (d_pm_rho)  { cudaFree(d_pm_rho);  d_pm_rho  = nullptr; }
    if (d_pm_phiw) { cudaFree(d_pm_phiw); d_pm_phiw = nullptr; }
    g_pm_ncell = 0; g_pm_ncell_cap = 0;
    g_pm_rho_cap = 0; g_pm_rho_ncell = 0;
    if (g_pm_slot_init) {
        for (int s = 0; s < MAX_CUDA_STREAMS; s++) {
            PmSlot* p = &g_pm_slot[s];
            if (p->d_x0)      cudaFree(p->d_x0);
            if (p->d_nbf)     cudaFree(p->d_nbf);
            if (p->d_px)      cudaFree(p->d_px);
            if (p->d_pv)      cudaFree(p->d_pv);
            if (p->d_pg)      cudaFree(p->d_pg);
            if (p->d_dteff)   cudaFree(p->d_dteff);
            if (p->d_new_v)   cudaFree(p->d_new_v);
            if (p->d_new_x)   cudaFree(p->d_new_x);
            if (p->d_phi_out) cudaFree(p->d_phi_out);
            if (p->d_err)     cudaFree(p->d_err);
            if (p->h_x0)      cudaFreeHost(p->h_x0);
            if (p->h_nbf)     cudaFreeHost(p->h_nbf);
            if (p->h_px)      cudaFreeHost(p->h_px);
            if (p->h_pv)      cudaFreeHost(p->h_pv);
            if (p->h_pg)      cudaFreeHost(p->h_pg);
            if (p->h_dteff)   cudaFreeHost(p->h_dteff);
            if (p->h_new_v)   cudaFreeHost(p->h_new_v);
            if (p->h_new_x)   cudaFreeHost(p->h_new_x);
            if (p->h_phi_out) cudaFreeHost(p->h_phi_out);
        }
        g_pm_slot_init = false;
    }
}

} // extern "C"

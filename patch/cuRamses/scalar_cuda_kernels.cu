// ==========================================================================
// GPU Newton-Gauss-Seidel sweeps for the nGR scalar-field solvers
// (f(R) Hu-Sawicki, symmetron, dilaton, nDGP, cubic Galileon).
//
// Mirrors the CPU sweeps in force_fine.kjhan.f90:
//  - red-black coloring by popcount(ind-1) parity, red then black per call,
//    halo exchange only after both colors (same as the CPU solve loop);
//  - same-level neighbors resolved through precomputed face (6) and
//    edge-diagonal (12) grid tables (vain_face_grid / vain_{xy,xz,yz}_grid);
//  - coarse-fine Dirichlet closure via values precomputed on the CPU once
//    per level solve (the parent level is frozen during the solve), with a
//    live-cell escape for same-level cells the grid tables cannot reach;
//  - res_max/src_max via block reduction + atomicMax (exact, order-free).
//
// For the Vainshtein solvers the CPU updates same-color cells sequentially,
// so diagonal (same-color) neighbor reads differ in order from the GPU's
// parallel sweep. Both orderings converge to the same tolerance-limited
// solution; bitwise sweep-for-sweep parity holds only for the 6-point
// solvers (f(R)/symmetron/dilaton).
// ==========================================================================

#include "cuda_stream_pool.h"
#include "amr_cuda_index.cuh"
#include <cstdio>
#include <cstring>
#include <cstdlib>

#define SCAL_SENTINEL 1.0e299

static double scal_ull_to_double(unsigned long long v)
{
    double d;
    memcpy(&d, &v, sizeof(d));
    return d;
}

// Model ids (must match scalar_cuda_interface.f90)
#define SCAL_MODEL_FR        0
#define SCAL_MODEL_SYMMETRON 1
#define SCAL_MODEL_DILATON   2
#define SCAL_MODEL_NDGP      3
#define SCAL_MODEL_GALILEON  4

typedef struct { double p[12]; } ScalParams;

// --------------------------------------------------------------------------
// Module state
// --------------------------------------------------------------------------
static cudaStream_t g_sc_stream = nullptr;

static double* d_sc_field = nullptr;   // scalar_gr, full ncell array
static double* d_sc_rho   = nullptr;   // rho, full ncell array
static long long g_sc_ncell = 0;
static long long g_sc_ncell_cap = 0;

static int* d_sc_igrid = nullptr;      // active grid list (1-based)
static int* d_sc_face  = nullptr;      // (6, ngrid) face-neighbor grids
static int* d_sc_edge  = nullptr;      // (12, ngrid) edge-diagonal grids
static int* d_sc_bnd_slot = nullptr;   // per grid: -1 or slot into bnd blocks
static int g_sc_ngrid = 0;
static int g_sc_grid_cap = 0;

static int*    d_sc_bnd_live = nullptr; // (noff,8,nbnd) live cell idx or 0
static double* d_sc_bnd_val  = nullptr; // (noff,8,nbnd) frozen value / sentinel
static int g_sc_bnd_cap = 0;            // capacity in doubles/ints
static int g_sc_noff = 6;

static unsigned long long* d_sc_maxred = nullptr; // [0]=res_max, [1]=src_max

// Halo exchange (same pattern as the MG phi halo)
static int*    d_sc_emit_cells = nullptr;
static int*    d_sc_recv_cells = nullptr;
static double* d_sc_emit_buf   = nullptr;
static double* d_sc_recv_buf   = nullptr;
static int g_sc_n_emit = 0, g_sc_n_recv = 0;
static int g_sc_emit_cap = 0, g_sc_recv_cap = 0;

static bool g_sc_ready = false;
static long long g_sc_uploads = 0;
static long long g_sc_sweeps = 0;
static int g_sc_block_size = 0;
static int g_sc_child_count = 0;

// --------------------------------------------------------------------------
// Device helpers
// --------------------------------------------------------------------------

// Offset table: ids 0..5 faces (-x,+x,-y,+y,-z,+z), 6..9 xy, 10..13 xz,
// 14..17 yz diagonals. Must match scalar_gpu_precompute_bnd (Fortran).
__constant__ int c_sc_off[18][3] = {
    {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1},
    {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},
    {1,0,1},{1,0,-1},{-1,0,1},{-1,0,-1},
    {0,1,1},{0,1,-1},{0,-1,1},{0,-1,-1}
};

// Resolve the same-level neighbor value of cell (my_grid, ind0) at offset
// off_id. Falls back to the precomputed boundary entry when the grid tables
// have no path; SCAL_SENTINEL in the table means "use u_c" (zero-gradient).
__device__ __forceinline__ double scal_nb_val(
    const double* __restrict__ field,
    int my_grid, int ind0, int off_id, double u_c,
    const int* __restrict__ face6, const int* __restrict__ edge12,
    int bnd_base, const int* __restrict__ bnd_live,
    const double* __restrict__ bnd_val,
    int ngridmax, int ncoarse, int block_size, int child_count, int noff)
{
    int ox = c_sc_off[off_id][0];
    int oy = c_sc_off[off_id][1];
    int oz = c_sc_off[off_id][2];

    int tx = (ind0 & 1) + ox;
    int ty = ((ind0 >> 1) & 1) + oy;
    int tz = ((ind0 >> 2) & 1) + oz;

    int dxs = (tx < 0) ? -1 : ((tx > 1) ? 1 : 0);
    int dys = (ty < 0) ? -1 : ((ty > 1) ? 1 : 0);
    int dzs = (tz < 0) ? -1 : ((tz > 1) ? 1 : 0);

    int nb_ind0 = (tx & 1) + 2 * (ty & 1) + 4 * (tz & 1);

    int g = 0;
    int nnz = (dxs != 0) + (dys != 0) + (dzs != 0);
    if (nnz == 0) {
        g = my_grid;
    } else if (nnz == 1) {
        int slot;
        if (dxs != 0)      slot = (dxs < 0) ? 0 : 1;
        else if (dys != 0) slot = (dys < 0) ? 2 : 3;
        else               slot = (dzs < 0) ? 4 : 5;
        g = face6[slot];
    } else {
        // Two axes: edge diagonal (offsets never use three axes)
        int idx;
        if (dzs == 0) {          // xy block
            idx = ((dxs + 1) / 2) * 2 + ((dys + 1) / 2);
        } else if (dys == 0) {   // xz block
            idx = 4 + ((dxs + 1) / 2) * 2 + ((dzs + 1) / 2);
        } else {                 // yz block
            idx = 8 + ((dys + 1) / 2) * 2 + ((dzs + 1) / 2);
        }
        g = edge12[idx];
    }

    if (g > 0) {
        const long long cell = amr_cuda_cell_1based(
            g, nb_ind0 + 1, ncoarse, block_size, child_count);
        return field[cell - 1];
    }

    // Boundary closure (precomputed on CPU)
    if (bnd_base >= 0) {
        int k = bnd_base + ind0 * noff + off_id;
        int live = bnd_live[k];
        if (live > 0) return field[live - 1];
        double v = bnd_val[k];
        if (v < SCAL_SENTINEL) return v;
    }
    return u_c; // zero-gradient fallback (matches fallback=u_c on CPU)
}

// Exact double atomicMax for non-negative values (IEEE bit pattern of
// non-negative doubles preserves ordering as unsigned integers).
__device__ __forceinline__ void scal_atomic_max_nonneg(
    unsigned long long* addr, double val)
{
    atomicMax(addr, __double_as_longlong(val));
}

// --------------------------------------------------------------------------
// The Newton-GS sweep kernel: one thread per active grid, 4 cells per color
// --------------------------------------------------------------------------
__global__ void scal_gs_kernel(
    double* __restrict__ field,
    const double* __restrict__ rho_arr,
    const int* __restrict__ igrid_arr,
    const int* __restrict__ face_arr,
    const int* __restrict__ edge_arr,
    const int* __restrict__ bnd_slot,
    const int* __restrict__ bnd_live,
    const double* __restrict__ bnd_val,
    int ngrid, int ngridmax, int ncoarse, int block_size, int child_count,
    int noff,
    int model, int color, int tracker,
    ScalParams sp,
    unsigned long long* __restrict__ maxred)
{
    extern __shared__ double sdata[]; // [blockDim]: res | [blockDim]: src

    int gidx = blockIdx.x * blockDim.x + threadIdx.x;
    double my_res = 0.0, my_src = 0.0;

    if (gidx < ngrid) {
        int my_grid = igrid_arr[gidx]; // 1-based
        int face6[6];
        for (int j = 0; j < 6; j++) face6[j] = face_arr[gidx * 6 + j];
        int edge12[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        if (noff > 6)
            for (int j = 0; j < 12; j++) edge12[j] = edge_arr[gidx * 12 + j];
        int slot = bnd_slot[gidx];
        int bnd_base = (slot >= 0) ? slot * 8 * noff : -1;

        const double dx2_inv = sp.p[0];

        for (int ind0 = 0; ind0 < 8; ind0++) {
            if ((__popc(ind0) & 1) != color) continue;

            long long icell = amr_cuda_cell_1based(
                my_grid, ind0 + 1, ncoarse, block_size, child_count);
            double u_c = field[icell - 1];

#define NBV(oid) scal_nb_val(field, my_grid, ind0, (oid), u_c, face6, edge12, \
                             bnd_base, bnd_live, bnd_val, ngridmax, ncoarse, \
                             block_size, child_count, noff)

            if (model == SCAL_MODEL_FR) {
                // p: 0=dx2_inv 1=a2_over_3 2=rho_coeff 3=R_bar 4=R_bar0
                //    5=fR0_abs 6=small_fR 7=inv_np1 8=np1 9=rho_tot
                double u_abs = fabs(u_c);
                double lapl = 0.0;
                for (int idim = 0; idim < 3; idim++) {
                    double u_l = NBV(2 * idim);
                    double u_r = NBV(2 * idim + 1);
                    lapl += (u_l + u_r - 2.0 * u_c) * dx2_inv;
                }
                double R_of_u;
                if (u_abs > sp.p[6]) {
                    int np1 = (int)sp.p[8];
                    if (np1 == 1)      R_of_u = sp.p[4] * sp.p[5] / u_abs;
                    else if (np1 == 2) R_of_u = sp.p[4] * sqrt(sp.p[5] / u_abs);
                    else               R_of_u = sp.p[4] * pow(sp.p[5] / u_abs, sp.p[7]);
                } else {
                    R_of_u = sp.p[3];
                }
                double source = sp.p[1] * (R_of_u - sp.p[3])
                              - sp.p[2] * (rho_arr[icell - 1] - sp.p[9]);
                double residual = lapl - source;
                double dR_du = (u_abs > sp.p[6]) ? (-R_of_u / (sp.p[8] * u_c)) : 0.0;
                double jacobian = -6.0 * dx2_inv - sp.p[1] * dR_du;
                if (fabs(jacobian) > 1e-30) {
                    double delta_u = -residual / jacobian;
                    if (fabs(delta_u) > 0.5 * u_abs && u_abs > sp.p[6])
                        delta_u = copysign(0.5 * u_abs, delta_u);
                    double u_new = u_c + delta_u;
                    if (u_new > 0.0) u_new = -0.5 * fmax(u_abs, sp.p[6]);
                    field[icell - 1] = u_new;
                }
                my_res = fmax(my_res, fabs(residual));
                my_src = fmax(my_src, fabs(source));

            } else if (model == SCAL_MODEL_SYMMETRON) {
                // p: 0=dx2_inv 1=a2_over_2L2 2=(a_ssb/a)^3
                double lapl = 0.0;
                for (int idim = 0; idim < 3; idim++) {
                    lapl += NBV(2 * idim) + NBV(2 * idim + 1);
                }
                lapl = (lapl - 6.0 * u_c) * dx2_inv;
                double rho_ratio = rho_arr[icell - 1] * sp.p[2];
                double mass_term = sp.p[1] * (rho_ratio - 1.0);
                double residual = lapl - mass_term * u_c - sp.p[1] * u_c * u_c * u_c;
                double jacobian = -6.0 * dx2_inv - mass_term
                                - 3.0 * sp.p[1] * u_c * u_c;
                if (fabs(jacobian) > 1e-30) {
                    field[icell - 1] = u_c - residual / jacobian;
                }
                my_src = fmax(my_src, fabs(mass_term * u_c)
                                    + fabs(sp.p[1] * u_c * u_c * u_c));
                my_res = fmax(my_res, fabs(residual));

            } else if (model == SCAL_MODEL_DILATON) {
                // p: 0=dx2_inv 1=cA 2=cV 3=pexp 4=vbar 5=chibar_d
                //    6=A2/beta0 7=-3*Om*beta0 8=-3*Om*A2*pexp 9=1e-30*chibar
                double lapl = 0.0;
                for (int idim = 0; idim < 3; idim++) {
                    lapl += NBV(2 * idim) + NBV(2 * idim + 1);
                }
                lapl = (lapl - 6.0 * u_c) * dx2_inv;
                double rho_c = rho_arr[icell - 1];
                double wfac = sp.p[6] * fmax(u_c, 1e-30);
                double vphi  = sp.p[7] * pow(wfac, sp.p[3]);
                double dvphi = sp.p[8] * pow(wfac, sp.p[3] - 1.0);
                double source = sp.p[1] * (rho_c * u_c - sp.p[5])
                              + sp.p[2] * (vphi - sp.p[4]);
                double residual = lapl - source;
                double jacobian = -6.0 * dx2_inv - sp.p[1] * rho_c - sp.p[2] * dvphi;
                if (fabs(jacobian) > 1e-30) {
                    double delta_u = -residual / jacobian;
                    if (fabs(delta_u) > 0.5 * fabs(u_c))
                        delta_u = copysign(0.5 * fabs(u_c), delta_u);
                    double u_new = u_c + delta_u;
                    if (u_new <= 0.0) u_new = 0.5 * fmax(fabs(u_c), sp.p[9]);
                    field[icell - 1] = u_new;
                }
                my_src = fmax(my_src, fabs(source));
                my_res = fmax(my_res, fabs(residual));

            } else {
                // nDGP / Galileon Vainshtein operator
                // p: 0=dx2_inv 1=coeff 2=src_coeff(Om*a/beta) 3=rho_tot
                //    4=sclamp_floor(1e-2*Om*a/|beta|) 5=dx2
                double phi_xm = NBV(0), phi_xp = NBV(1);
                double phi_ym = NBV(2), phi_yp = NBV(3);
                double phi_zm = NBV(4), phi_zp = NBV(5);

                double lapl = (phi_xp + phi_xm + phi_yp + phi_ym
                             + phi_zp + phi_zm - 6.0 * u_c) * dx2_inv;
                double phi_xx = (phi_xp + phi_xm - 2.0 * u_c) * dx2_inv;
                double phi_yy = (phi_yp + phi_ym - 2.0 * u_c) * dx2_inv;
                double phi_zz = (phi_zp + phi_zm - 2.0 * u_c) * dx2_inv;

                double dpp, dpm, dmp, dmm, dmix_du;
                double mix_xy2, mix_xz2, mix_yz2;
                double phi_pp, phi_pm, phi_mp, phi_mm;

                // xy
                phi_pp = NBV(6); phi_pm = NBV(7); phi_mp = NBV(8); phi_mm = NBV(9);
                dpp = (phi_pp - phi_xp - phi_yp + u_c) * dx2_inv;
                dpm = (phi_xp - phi_pm - u_c + phi_ym) * dx2_inv;
                dmp = (phi_yp - phi_mp - u_c + phi_xm) * dx2_inv;
                dmm = (u_c - phi_ym - phi_xm + phi_mm) * dx2_inv;
                mix_xy2 = 0.25 * (dpp*dpp + dpm*dpm + dmp*dmp + dmm*dmm);
                if (model == SCAL_MODEL_NDGP && tracker) {
                    double cmix = 0.25 * (phi_pp - phi_pm - phi_mp + phi_mm) * dx2_inv;
                    mix_xy2 = cmix * cmix;
                }
                dmix_du = 0.5 * dx2_inv * (dpp - dpm - dmp + dmm);

                // xz
                phi_pp = NBV(10); phi_pm = NBV(11); phi_mp = NBV(12); phi_mm = NBV(13);
                dpp = (phi_pp - phi_xp - phi_zp + u_c) * dx2_inv;
                dpm = (phi_xp - phi_pm - u_c + phi_zm) * dx2_inv;
                dmp = (phi_zp - phi_mp - u_c + phi_xm) * dx2_inv;
                dmm = (u_c - phi_zm - phi_xm + phi_mm) * dx2_inv;
                mix_xz2 = 0.25 * (dpp*dpp + dpm*dpm + dmp*dmp + dmm*dmm);
                if (model == SCAL_MODEL_NDGP && tracker) {
                    double cmix = 0.25 * (phi_pp - phi_pm - phi_mp + phi_mm) * dx2_inv;
                    mix_xz2 = cmix * cmix;
                }
                dmix_du += 0.5 * dx2_inv * (dpp - dpm - dmp + dmm);

                // yz
                phi_pp = NBV(14); phi_pm = NBV(15); phi_mp = NBV(16); phi_mm = NBV(17);
                dpp = (phi_pp - phi_yp - phi_zp + u_c) * dx2_inv;
                dpm = (phi_yp - phi_pm - u_c + phi_zm) * dx2_inv;
                dmp = (phi_zp - phi_mp - u_c + phi_ym) * dx2_inv;
                dmm = (u_c - phi_zm - phi_ym + phi_mm) * dx2_inv;
                mix_yz2 = 0.25 * (dpp*dpp + dpm*dpm + dmp*dmp + dmm*dmm);
                if (model == SCAL_MODEL_NDGP && tracker) {
                    double cmix = 0.25 * (phi_pp - phi_pm - phi_mp + phi_mm) * dx2_inv;
                    mix_yz2 = cmix * cmix;
                }
                dmix_du += 0.5 * dx2_inv * (dpp - dpm - dmp + dmm);

                double lapl2 = lapl * lapl;
                double trace_ij2 = phi_xx*phi_xx + phi_yy*phi_yy + phi_zz*phi_zz
                                 + 2.0 * (mix_xy2 + mix_xz2 + mix_yz2);
                double vain_term = lapl2 - trace_ij2;

                double source = sp.p[2] * (rho_arr[icell - 1] - sp.p[3]);

                double residual, jacobian, a_tgt = source;
                if (model == SCAL_MODEL_GALILEON) {
                    double coeff_G = sp.p[1];
                    double tbar_ij = fmax(trace_ij2 - lapl2 / 3.0, 0.0);
                    double qcoeff = 2.0 * coeff_G / 3.0;
                    double disc = 1.0 + 4.0 * qcoeff * (source + coeff_G * tbar_ij);
                    if (fabs(qcoeff) > 1e-30) {
                        a_tgt = (disc > 0.0) ? (-1.0 + sqrt(disc)) / (2.0 * qcoeff)
                                             : -1.0 / (2.0 * qcoeff);
                    }
                    if (tracker && fabs(qcoeff) > 1e-30) {
                        residual = lapl - a_tgt;
                        jacobian = -6.0 * dx2_inv;
                    } else {
                        residual = lapl + coeff_G * vain_term - source;
                        jacobian = -6.0 * dx2_inv
                                 + coeff_G * (-8.0 * lapl * dx2_inv - 2.0 * dmix_du);
                    }
                } else { // nDGP
                    residual = lapl + sp.p[1] * vain_term - source;
                    jacobian = -6.0 * dx2_inv
                             + sp.p[1] * (-8.0 * lapl * dx2_inv - 2.0 * dmix_du);
                }

                if (fabs(jacobian) > 1e-30) {
                    double delta_u = -residual / jacobian;
                    double sclamp;
                    if (model == SCAL_MODEL_GALILEON && tracker)
                        sclamp = 0.5 * sp.p[5] * fmax(fmax(fabs(a_tgt), fabs(source)), sp.p[4]);
                    else
                        sclamp = 0.5 * sp.p[5] * fmax(fabs(source), sp.p[4]);
                    if (fabs(delta_u) > sclamp)
                        delta_u = copysign(sclamp, delta_u);
                    field[icell - 1] = u_c + delta_u;
                }
                my_res = fmax(my_res, fabs(residual));
                my_src = fmax(my_src, fabs(source));
            }
#undef NBV
        }
    }

    // Block reduction of the two maxima, then one atomic per block
    double* sres = sdata;
    double* ssrc = sdata + blockDim.x;
    sres[threadIdx.x] = my_res;
    ssrc[threadIdx.x] = my_src;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            sres[threadIdx.x] = fmax(sres[threadIdx.x], sres[threadIdx.x + s]);
            ssrc[threadIdx.x] = fmax(ssrc[threadIdx.x], ssrc[threadIdx.x + s]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        scal_atomic_max_nonneg(&maxred[0], sres[0]);
        scal_atomic_max_nonneg(&maxred[1], ssrc[0]);
    }
}

// --------------------------------------------------------------------------
// Halo kernels (same as MG phi halo, on the scalar field)
// --------------------------------------------------------------------------
__global__ void scal_halo_gather_kernel(
    const double* __restrict__ field, double* __restrict__ buf,
    const int* __restrict__ cells, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) buf[i] = field[cells[i] - 1];
}

__global__ void scal_halo_scatter_kernel(
    double* __restrict__ field, const double* __restrict__ buf,
    const int* __restrict__ cells, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) field[cells[i] - 1] = buf[i];
}

// --------------------------------------------------------------------------
// C API
// --------------------------------------------------------------------------
extern "C" {

void cuda_scal_upload(const double* field, const double* rho,
                      long long ncell,
                      const int* igrid, const int* face, const int* edge,
                      const int* bnd_slot,
                      const int* bnd_live, const double* bnd_val,
                      int ngrid, int nbnd, int noff)
{
    g_sc_ready = false;
    g_sc_uploads = 0;
    g_sc_sweeps = 0;
    g_sc_block_size = 0;
    g_sc_child_count = 0;
    if (!is_pool_initialized()) return;
    cudaGetLastError();

    if (!g_sc_stream)
        cudaStreamCreateWithFlags(&g_sc_stream, cudaStreamNonBlocking);

    // Full-ncell arrays
    if (ncell > g_sc_ncell_cap) {
        if (d_sc_field) { cudaFree(d_sc_field); d_sc_field = nullptr; }
        if (d_sc_rho)   { cudaFree(d_sc_rho);   d_sc_rho   = nullptr; }
        size_t free_mem = 0, total_mem = 0;
        cudaMemGetInfo(&free_mem, &total_mem);
        size_t need = (size_t)ncell * 2 * sizeof(double);
        if (need > free_mem * 9 / 10) {
            fprintf(stderr, "CUDA scalar: SKIP — need %.1f GB, %.1f GB free\n",
                    (double)need / 1073741824.0, (double)free_mem / 1073741824.0);
            g_sc_ncell_cap = 0;
            return;
        }
        cudaError_t e1 = cudaMalloc(&d_sc_field, (size_t)ncell * sizeof(double));
        cudaError_t e2 = cudaMalloc(&d_sc_rho,   (size_t)ncell * sizeof(double));
        if (e1 || e2) {
            fprintf(stderr, "CUDA scalar: cell alloc FAILED\n");
            if (d_sc_field) { cudaFree(d_sc_field); d_sc_field = nullptr; }
            if (d_sc_rho)   { cudaFree(d_sc_rho);   d_sc_rho   = nullptr; }
            g_sc_ncell_cap = 0;
            return;
        }
        g_sc_ncell_cap = ncell;
    }
    g_sc_ncell = ncell;

    // Grid tables
    if (ngrid > g_sc_grid_cap) {
        if (d_sc_igrid)    { cudaFree(d_sc_igrid);    d_sc_igrid    = nullptr; }
        if (d_sc_face)     { cudaFree(d_sc_face);     d_sc_face     = nullptr; }
        if (d_sc_edge)     { cudaFree(d_sc_edge);     d_sc_edge     = nullptr; }
        if (d_sc_bnd_slot) { cudaFree(d_sc_bnd_slot); d_sc_bnd_slot = nullptr; }
        int cap = ngrid * 2;
        cudaError_t e1 = cudaMalloc(&d_sc_igrid,    (size_t)cap * sizeof(int));
        cudaError_t e2 = cudaMalloc(&d_sc_face,     (size_t)cap * 6 * sizeof(int));
        cudaError_t e3 = cudaMalloc(&d_sc_edge,     (size_t)cap * 12 * sizeof(int));
        cudaError_t e4 = cudaMalloc(&d_sc_bnd_slot, (size_t)cap * sizeof(int));
        if (e1 || e2 || e3 || e4) {
            fprintf(stderr, "CUDA scalar: grid alloc FAILED\n");
            g_sc_grid_cap = 0;
            return;
        }
        g_sc_grid_cap = cap;
    }
    g_sc_ngrid = ngrid;
    g_sc_noff = noff;

    // Boundary blocks
    int bnd_n = nbnd * 8 * noff;
    if (bnd_n > g_sc_bnd_cap) {
        if (d_sc_bnd_live) { cudaFree(d_sc_bnd_live); d_sc_bnd_live = nullptr; }
        if (d_sc_bnd_val)  { cudaFree(d_sc_bnd_val);  d_sc_bnd_val  = nullptr; }
        int cap = bnd_n * 2;
        cudaError_t e1 = cudaMalloc(&d_sc_bnd_live, (size_t)cap * sizeof(int));
        cudaError_t e2 = cudaMalloc(&d_sc_bnd_val,  (size_t)cap * sizeof(double));
        if (e1 || e2) {
            fprintf(stderr, "CUDA scalar: bnd alloc FAILED\n");
            g_sc_bnd_cap = 0;
            return;
        }
        g_sc_bnd_cap = cap;
    }

    if (!d_sc_maxred) cudaMalloc(&d_sc_maxred, 2 * sizeof(unsigned long long));

    cudaMemcpy(d_sc_field, field, (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sc_rho,   rho,   (size_t)ncell * sizeof(double), cudaMemcpyHostToDevice);
    if (ngrid > 0) {
        cudaMemcpy(d_sc_igrid, igrid, (size_t)ngrid * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sc_face,  face,  (size_t)ngrid * 6 * sizeof(int), cudaMemcpyHostToDevice);
        if (noff > 6)
            cudaMemcpy(d_sc_edge, edge, (size_t)ngrid * 12 * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sc_bnd_slot, bnd_slot, (size_t)ngrid * sizeof(int), cudaMemcpyHostToDevice);
    }
    if (bnd_n > 0) {
        cudaMemcpy(d_sc_bnd_live, bnd_live, (size_t)bnd_n * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_sc_bnd_val,  bnd_val,  (size_t)bnd_n * sizeof(double), cudaMemcpyHostToDevice);
    }

    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA scalar upload error: %s\n", cudaGetErrorString(err));
        return;
    }
    g_sc_ready = true;
    g_sc_uploads = 1;
}

int cuda_scal_is_ready(void) { return g_sc_ready ? 1 : 0; }

// One full sweep: red then black (halo exchange happens after both,
// mirroring the CPU solve loop).
void cuda_scal_sweep(int model, const double* params,
                     int ngridmax, int ncoarse, int block_size, int child_count,
                     int tracker,
                     double* res_max, double* src_max)
{
    *res_max = 0.0;
    *src_max = 0.0;
    if (!g_sc_ready || g_sc_ngrid <= 0) return;
    const long long expected =
        (long long)ncoarse + (long long)child_count * ngridmax;
    if (!amr_cuda_layout_valid(g_sc_ncell, ncoarse, ngridmax,
                               block_size, child_count, (1 << NDIM)) ||
        (g_sc_block_size != 0 && g_sc_block_size != block_size) ||
        (g_sc_child_count != 0 && g_sc_child_count != child_count)) {
        fprintf(stderr,
                "CUDA scalar invalid block layout B=%d C=%d "
                "ngm=%d ncell=%lld expected=%lld\n",
                block_size, child_count, ngridmax, g_sc_ncell, expected);
        fflush(stderr);
        std::abort();
    }
    g_sc_block_size = block_size;
    g_sc_child_count = child_count;

    ScalParams sp;
    for (int i = 0; i < 12; i++) sp.p[i] = params[i];

    cudaMemsetAsync(d_sc_maxred, 0, 2 * sizeof(unsigned long long), g_sc_stream);

    int block = 128;
    int grid = (g_sc_ngrid + block - 1) / block;
    size_t smem = 2 * block * sizeof(double);

    for (int color = 0; color < 2; color++) {
        scal_gs_kernel<<<grid, block, smem, g_sc_stream>>>(
            d_sc_field, d_sc_rho, d_sc_igrid, d_sc_face, d_sc_edge,
            d_sc_bnd_slot, d_sc_bnd_live, d_sc_bnd_val,
            g_sc_ngrid, ngridmax, ncoarse, block_size, child_count, g_sc_noff,
            model, color, tracker, sp, d_sc_maxred);
    }
    const cudaError_t sync_error = cudaStreamSynchronize(g_sc_stream);

    unsigned long long h_red[2];
    const cudaError_t copy_error = cudaMemcpy(
        h_red, d_sc_maxred, 2 * sizeof(unsigned long long),
        cudaMemcpyDeviceToHost);

    const cudaError_t launch_error = cudaGetLastError();
    const cudaError_t err = sync_error != cudaSuccess ? sync_error :
                            (copy_error != cudaSuccess ? copy_error : launch_error);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA scalar sweep error: %s\n", cudaGetErrorString(err));
        g_sc_ready = false;
    } else {
        *res_max = scal_ull_to_double(h_red[0]);
        *src_max = scal_ull_to_double(h_red[1]);
        g_sc_sweeps++;
    }
}

void cuda_scal_download(double* field, long long ncell)
{
    if (!g_sc_ready || ncell != g_sc_ncell) return;
    cudaMemcpy(field, d_sc_field, (size_t)ncell * sizeof(double),
               cudaMemcpyDeviceToHost);
}

void cuda_scal_halo_setup(const int* emit, int n_emit,
                          const int* recv, int n_recv)
{
    if (!is_pool_initialized()) return;
    if (n_emit > g_sc_emit_cap) {
        if (d_sc_emit_cells) cudaFree(d_sc_emit_cells);
        if (d_sc_emit_buf)   cudaFree(d_sc_emit_buf);
        int cap = n_emit * 2;
        cudaMalloc(&d_sc_emit_cells, (size_t)cap * sizeof(int));
        cudaMalloc(&d_sc_emit_buf,   (size_t)cap * sizeof(double));
        g_sc_emit_cap = cap;
    }
    if (n_recv > g_sc_recv_cap) {
        if (d_sc_recv_cells) cudaFree(d_sc_recv_cells);
        if (d_sc_recv_buf)   cudaFree(d_sc_recv_buf);
        int cap = n_recv * 2;
        cudaMalloc(&d_sc_recv_cells, (size_t)cap * sizeof(int));
        cudaMalloc(&d_sc_recv_buf,   (size_t)cap * sizeof(double));
        g_sc_recv_cap = cap;
    }
    if (n_emit > 0)
        cudaMemcpy(d_sc_emit_cells, emit, (size_t)n_emit * sizeof(int),
                   cudaMemcpyHostToDevice);
    if (n_recv > 0)
        cudaMemcpy(d_sc_recv_cells, recv, (size_t)n_recv * sizeof(int),
                   cudaMemcpyHostToDevice);
    g_sc_n_emit = n_emit;
    g_sc_n_recv = n_recv;
}

void cuda_scal_halo_gather(double* buf, int n)
{
    if (!g_sc_ready || n <= 0 || n != g_sc_n_emit) return;
    int block = 256;
    int grid = (n + block - 1) / block;
    scal_halo_gather_kernel<<<grid, block, 0, g_sc_stream>>>(
        d_sc_field, d_sc_emit_buf, d_sc_emit_cells, n);
    cudaStreamSynchronize(g_sc_stream);
    cudaMemcpy(buf, d_sc_emit_buf, (size_t)n * sizeof(double),
               cudaMemcpyDeviceToHost);
}

void cuda_scal_halo_scatter(const double* buf, int n)
{
    if (!g_sc_ready || n <= 0 || n != g_sc_n_recv) return;
    cudaMemcpy(d_sc_recv_buf, buf, (size_t)n * sizeof(double),
               cudaMemcpyHostToDevice);
    int block = 256;
    int grid = (n + block - 1) / block;
    scal_halo_scatter_kernel<<<grid, block, 0, g_sc_stream>>>(
        d_sc_field, d_sc_recv_buf, d_sc_recv_cells, n);
    cudaStreamSynchronize(g_sc_stream);
}

// Free the big per-solve arrays (grid/bnd/halo tables keep capacity)
void cuda_scal_release(void)
{
    if (g_sc_uploads + g_sc_sweeps > 0) {
        printf("[CUDA_NGR] B=%d C=%d uploads=%lld scalar_sweeps=%lld\n",
               g_sc_block_size, g_sc_child_count, g_sc_uploads, g_sc_sweeps);
        fflush(stdout);
    }
    g_sc_uploads = 0;
    g_sc_sweeps = 0;
    g_sc_ready = false;
    if (d_sc_field) { cudaFree(d_sc_field); d_sc_field = nullptr; }
    if (d_sc_rho)   { cudaFree(d_sc_rho);   d_sc_rho   = nullptr; }
    g_sc_ncell = 0;
    g_sc_ncell_cap = 0;
}

void cuda_scal_finalize(void)
{
    cuda_scal_release();
    if (d_sc_igrid)      { cudaFree(d_sc_igrid);      d_sc_igrid      = nullptr; }
    if (d_sc_face)       { cudaFree(d_sc_face);       d_sc_face       = nullptr; }
    if (d_sc_edge)       { cudaFree(d_sc_edge);       d_sc_edge       = nullptr; }
    if (d_sc_bnd_slot)   { cudaFree(d_sc_bnd_slot);   d_sc_bnd_slot   = nullptr; }
    if (d_sc_bnd_live)   { cudaFree(d_sc_bnd_live);   d_sc_bnd_live   = nullptr; }
    if (d_sc_bnd_val)    { cudaFree(d_sc_bnd_val);    d_sc_bnd_val    = nullptr; }
    if (d_sc_maxred)     { cudaFree(d_sc_maxred);     d_sc_maxred     = nullptr; }
    if (d_sc_emit_cells) { cudaFree(d_sc_emit_cells); d_sc_emit_cells = nullptr; }
    if (d_sc_recv_cells) { cudaFree(d_sc_recv_cells); d_sc_recv_cells = nullptr; }
    if (d_sc_emit_buf)   { cudaFree(d_sc_emit_buf);   d_sc_emit_buf   = nullptr; }
    if (d_sc_recv_buf)   { cudaFree(d_sc_recv_buf);   d_sc_recv_buf   = nullptr; }
    g_sc_grid_cap = 0; g_sc_bnd_cap = 0;
    g_sc_emit_cap = 0; g_sc_recv_cap = 0;
    if (g_sc_stream) { cudaStreamDestroy(g_sc_stream); g_sc_stream = nullptr; }
}

} // extern "C"

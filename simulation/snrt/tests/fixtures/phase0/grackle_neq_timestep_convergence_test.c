#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,
                         double,double,double,double,double);

enum { NCELL = 16, NSPEC = 9, NLEVEL = 18 };
static const int subdivisions[NLEVEL] = {
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
    16384, 32768, 65536, 131072
};
static const double density_grid[4] = {1e-4, 1e-2, 1.0, 1e2};
static const double energy_grid[4] = {1e10, 1e12, 1e14, 1e16};
/* With density_unit=m_H/cc, Grackle's comoving conversion preserves the
 * physical H number density at 0.76*density_grid (not /a^3). With
 * length_unit=time_unit=1 and a=0.1, its specific-energy unit is
 * [LU/(a TU)]^2=100 erg/g, so these inputs span 1e12..1e18 erg/g.
 * The qualification trajectory reaches physical intervals of order 7e14 s
 * near a=0.12; this matrix brackets that interval and a margin. */
static const double dts[] = {1e10, 1e12, 1e14, 7e14, 1e15};
static const double density_unit = 1.6726219e-24;

struct state {
    double rho[NCELL], energy[NCELL], species[NSPEC * NCELL];
    double metal[NCELL], length[NCELL];
    int cooling[NCELL];
};

#define SP(s,k,i) ((s)->species[(k) * NCELL + (i)])

static void initial_state(struct state *s)
{
    memset(s, 0, sizeof(*s));
    for (int i = 0; i < NCELL; ++i) {
        const double rho = density_grid[i / 4];
        const double hfrac = 0.76;
        const double h = hfrac * rho * (1.0 - 1e-3);
        const double he = (1.0 - hfrac) * rho * (1.0 - 1e-3);
        s->rho[i] = rho;
        s->energy[i] = energy_grid[i % 4];
        s->metal[i] = 1e-3 * rho;
        s->length[i] = 1.0;
        s->cooling[i] = 1;
        SP(s, 0, i) = h * (1.0 - 2e-4 - 2e-6);
        SP(s, 1, i) = h * 2e-4;
        SP(s, 3, i) = he;
        SP(s, 6, i) = h * 2e-6;
        SP(s, 8, i) = SP(s, 1, i);
    }
}

static int evolve(const struct state *initial, double dt, int nstep,
                  struct state *result)
{
    *result = *initial;
    for (int step = 0; step < nstep; ++step) {
        int status = lag_grackle_neq_step(
            NCELL, result->rho, result->energy, result->species,
            result->metal, result->length, result->cooling,
            density_unit, 1.0, 1.0, 0.1, dt / nstep);
        if (status != 0) return status;
    }
    return 0;
}

static int valid_state(const struct state *s)
{
    for (int i = 0; i < NCELL; ++i) {
        double h = 0.0, he = 0.0;
        double charge = SP(s, 1, i) - SP(s, 2, i) + SP(s, 7, i) / 2.0
                      + SP(s, 4, i) / 4.0 + SP(s, 5, i) / 2.0;
        for (int k = 0; k < NSPEC; ++k) {
            if (!isfinite(SP(s, k, i)) || SP(s, k, i) < 0.0) return 0;
        }
        h = SP(s, 0, i) + SP(s, 1, i) + SP(s, 2, i)
          + SP(s, 6, i) + SP(s, 7, i);
        he = SP(s, 3, i) + SP(s, 4, i) + SP(s, 5, i);
        if (!isfinite(s->energy[i]) || s->energy[i] <= 0.0) return 0;
        if (fabs(h + he + s->metal[i] - s->rho[i]) > 1e-6 * s->rho[i]) return 0;
        if (fabs(SP(s, 8, i) - charge) > 1e-6 * s->rho[i]) return 0;
    }
    return 1;
}

static void compare_states(const struct state *a, const struct state *b,
                           double *max_energy_rel, double *max_species_massfrac,
                           double *max_h_l1, double *max_he_l1,
                           int *max_energy_cell, int *max_species_cell)
{
    *max_energy_rel = *max_species_massfrac = *max_h_l1 = *max_he_l1 = 0.0;
    *max_energy_cell = *max_species_cell = 0;
    for (int i = 0; i < NCELL; ++i) {
        double h0 = 0.0, he0 = 0.0, h_l1 = 0.0, he_l1 = 0.0;
        double de = fabs(a->energy[i] - b->energy[i]) / fmax(fabs(b->energy[i]), 1e-300);
        if (de > *max_energy_rel) {
            *max_energy_rel = de;
            *max_energy_cell = i;
        }
        for (int k = 0; k < NSPEC; ++k) {
            double ds = fabs(SP(a, k, i) - SP(b, k, i)) / a->rho[i];
            if (ds > *max_species_massfrac) {
                *max_species_massfrac = ds;
                *max_species_cell = i;
            }
            if (k == 0 || k == 1 || k == 2 || k == 6 || k == 7) {
                h0 += SP(b, k, i);
                h_l1 += fabs(SP(a, k, i) - SP(b, k, i));
            } else if (k == 3 || k == 4 || k == 5) {
                he0 += SP(b, k, i);
                he_l1 += fabs(SP(a, k, i) - SP(b, k, i));
            }
        }
        h_l1 /= fmax(h0, 1e-300);
        he_l1 /= fmax(he0, 1e-300);
        if (h_l1 > *max_h_l1) *max_h_l1 = h_l1;
        if (he_l1 > *max_he_l1) *max_he_l1 = he_l1;
    }
}

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    omp_set_num_threads(1);
    if (lag_grackle_neq_init(argv[1], 1, 0, 0, 1,
                             density_unit, 1.0, 1.0, 0.1, 10.0) != 0) return 3;

    struct state initial, states[NLEVEL];
    initial_state(&initial);
    for (int it = 0; it < (int)(sizeof(dts) / sizeof(dts[0])); ++it) {
        for (int level = 0; level < NLEVEL; ++level) {
            int status = evolve(&initial, dts[it], subdivisions[level], &states[level]);
            if (status != 0 || !valid_state(&states[level])) {
                fprintf(stderr, "FAIL dt=%.1e subdivisions=%d solver_status=%d\n",
                        dts[it], subdivisions[level], status);
                return 4;
            }
        }
        printf("convergence dt=%.1e s (UVB=HM2012, a=0.1):\n", dts[it]);
        printf(" subdivisions  max_dE/E(ref131072)  cell  max_dspecies/rho  cell  max_H_L1  max_He_L1\n");
        for (int level = 0; level < NLEVEL - 1; ++level) {
            double de, ds, dh, dhe;
            int ecell, scall;
            compare_states(&states[level], &states[NLEVEL - 1], &de, &ds, &dh, &dhe,
                           &ecell, &scall);
            printf(" %12d  %17.6e %5d  %16.6e %5d  %9.3e  %10.3e\n",
                   subdivisions[level], de, ecell, ds, scall, dh, dhe);
        }
        double de, ds, dh, dhe;
        int ecell, scall;
        compare_states(&states[NLEVEL - 2], &states[NLEVEL - 1], &de, &ds, &dh, &dhe,
                       &ecell, &scall);
        printf(" adjacent(%d,%d): dE/E=%.6e cell=%d(rho_code=%.3e,u0=%.3e) "
               "dspecies/rho=%.6e cell=%d(rho_code=%.3e,u0=%.3e) H_L1=%.6e He_L1=%.6e\n",
               subdivisions[NLEVEL - 2], subdivisions[NLEVEL - 1], de, ecell,
               density_grid[ecell/4], energy_grid[ecell%4], ds, scall,
               density_grid[scall/4], energy_grid[scall%4], dh, dhe);
    }
    return 0;
}

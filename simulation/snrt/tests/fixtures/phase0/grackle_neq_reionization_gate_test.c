#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,
                         double,double,double,double,double);

int main(int argc, char **argv)
{
    if (argc != 4) return 2;
    const int uv = 1;
    const double a = strtod(argv[2], NULL);
    const double z_reion = strtod(argv[3], NULL);
    const double rho0 = 1.0e-4;
    const double density_unit = 1.6726219e-24;
    const double length_unit = 3.0856776e21;
    const double time_unit = 3.15576e13;
    double rho[1] = {rho0}, energy[1] = {1.0e12};
    double metal[1] = {1.0e-8}, length[1] = {1.0e-5};
    double species[9] = {0};
    int cooling[1] = {1};
    const double hydrogen = 0.76 * (rho0 - metal[0]);
    const double helium = 0.24 * (rho0 - metal[0]);

    species[0] = hydrogen * (1.0 - 2.0e-4 - 2.0e-6);
    species[1] = hydrogen * 2.0e-4;
    species[3] = helium;
    species[6] = hydrogen * 2.0e-6;
    species[8] = species[1];

    if (lag_grackle_neq_init(argv[1], uv, 0, 0, 1, density_unit,
                             length_unit, time_unit, a, z_reion) != 0) return 3;
    if (lag_grackle_neq_step(1, rho, energy, species, metal, length, cooling,
                             density_unit, length_unit, time_unit, a, 1.0e-2) != 0) return 4;

    printf("RESULT a=%.9g z_reion=%.9g u=%.17g", a, z_reion, energy[0]);
    for (int k = 0; k < 9; ++k) printf(" s%d=%.17g", k, species[k]);
    putchar('\n');
    return 0;
}

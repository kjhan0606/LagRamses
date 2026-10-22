/* Analysis-only adapter: evaluate Grackle temperature, never evolve a state.
 * Columns are rho, Grackle specific energy, nine species, metal density.
 * Energy/unit conversion is performed by the snapshot reader. */
#include <grackle.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

int snapshot_temperature(int n, double *p, double density_unit,
                         double length_unit, double time_unit, double a,
                         const char *data_file, double *temperature)
{
    if(n<1 || n>INT_MAX/12 || sizeof(gr_float)!=sizeof(double) ||
       !gr_check_consistency() || !isfinite(a) || a<=0)return 1;
    chemistry_data c;
    chemistry_data_storage rates={0};
    code_units units={1,density_unit,length_unit,time_unit,
                      length_unit/(a*time_unit),1.0,a};
    if(!local_initialize_chemistry_parameters(&c))return 2;
    c.use_grackle=1; c.primordial_chemistry=2;
    c.Gamma=5.0/3.0; c.metal_cooling=1; c.UVbackground=0;
    c.grackle_data_file=(char *)data_file; c.omp_nthreads=1;
    c.h2_on_dust=0; c.dust_chemistry=0; c.photoelectric_heating=0;
    c.use_radiative_transfer=0; c.cmb_temperature_floor=0;
    if(!local_initialize_chemistry_data(&c,&rates,&units))return 3;
    int dims[3]={n,1,1},start[3]={0,0,0},end[3]={n-1,0,0};
    grackle_field_data f;
    gr_initialize_field_data(&f);
    f.grid_rank=1; f.grid_dimension=dims; f.grid_start=start; f.grid_end=end;
    f.density=p; f.internal_energy=p+n;
    f.HI_density=p+2*n; f.HII_density=p+3*n; f.HM_density=p+4*n;
    f.HeI_density=p+5*n; f.HeII_density=p+6*n; f.HeIII_density=p+7*n;
    f.H2I_density=p+8*n; f.H2II_density=p+9*n; f.e_density=p+10*n;
    f.metal_density=p+11*n;
    int ok=local_calculate_temperature(&c,&rates,&units,&f,temperature);
    local_free_chemistry_data(&c,&rates);
    return ok?0:4;
}

/* Native-double stdin/stdout protocol avoids changing/rebuilding the existing
 * non-PIC static Grackle installation merely for an analysis adapter. */
int main(int argc,char **argv)
{
    if(argc!=7)return 10;
    int n=atoi(argv[1]);
    if(n<1 || n>INT_MAX/12)return 11;
    double *p=malloc((size_t)n*12*sizeof(double));
    double *temperature=malloc((size_t)n*sizeof(double));
    if(!p || !temperature)return 12;
    if(fread(p,sizeof(double),(size_t)n*12,stdin)!=(size_t)n*12)return 13;
    int status=snapshot_temperature(n,p,strtod(argv[2],NULL),strtod(argv[3],NULL),
                                    strtod(argv[4],NULL),strtod(argv[5],NULL),argv[6],temperature);
    if(!status && fwrite(temperature,sizeof(double),n,stdout)!=(size_t)n)status=14;
    free(p);free(temperature);
    return status;
}

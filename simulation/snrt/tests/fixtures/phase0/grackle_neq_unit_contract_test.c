#include <grackle.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,
                         double,double,double,double,double);

struct state {
    double rho, energy, species[9], metal, length;
    int cooling;
};

static void set_state(struct state *s,double a,double rho_proper,double hfrac,
                      double xhii,double temperature_over_mu)
{
    const double mh=1.6726219e-24, kb=1.380649e-16;
    const double rho_base=mh;             /* one H-mass per cc at a=1 */
    const double length_base=1.0e21;      /* comoving cm */
    const double time_unit=1.0e15;        /* fixed Grackle seconds */
    const double density_unit=rho_base/(a*a*a);
    const double length_unit=length_base*a;
    const double velocity_unit=length_unit/(a*time_unit);
    const double t2_unit=mh/kb*velocity_unit*velocity_unit;
    const double h=hfrac*rho_proper/density_unit;
    const double he=(1.0-hfrac)*rho_proper/density_unit;

    memset(s,0,sizeof(*s));
    s->rho=rho_proper/density_unit;
    s->metal=0.0;
    s->length=1.0;
    /* Keep radiative cooling off in this adapter-vs-direct check: Compton
     * redshift-dependent physics, not a unit-conversion error. Chemistry still
     * uses the mapped specific energy to evaluate its temperature-dependent rates. */
    s->cooling=0;
    s->species[0]=h*(1.0-xhii);
    s->species[1]=h*xhii;
    s->species[3]=he;
    s->species[8]=s->species[1];
    s->energy=temperature_over_mu/(t2_unit*(2.0/3.0));
}

static int evolve(struct state *s,double a,double rho0,double l0,double dt)
{
    const double dunit=rho0/(a*a*a),lunit=l0*a,tunit=1.0e15;
    return lag_grackle_neq_step(1,&s->rho,&s->energy,s->species,&s->metal,
                                &s->length,&s->cooling,dunit,lunit,tunit,a,dt);
}

static double max_fraction_difference(const struct state *x,const struct state *y)
{
    double result=fabs(x->energy-y->energy)/fmax(fabs(y->energy),1e-300);
    for(int k=0;k<9;k++)
        result=fmax(result,fabs(x->species[k]/x->rho-y->species[k]/y->rho));
    return result;
}

/* Per-epoch direct Grackle reference, bypassing the adapter's pack/sort path. */
static int evolve_direct(const char *path,struct state *s,double a,
                         double rho0,double l0,double dt)
{
    chemistry_data chemistry;
    chemistry_data_storage rates={0};
    const double dunit=rho0/(a*a*a),lunit=l0*a,tunit=1.0e15;
    code_units units={1,dunit,lunit,tunit,lunit/(a*tunit),1.0,a};
    int dimensions[3]={1,1,1},start[3]={0,0,0},end[3]={0,0,0};
    grackle_field_data fields;

    if(!local_initialize_chemistry_parameters(&chemistry))return 1;
    chemistry.use_grackle=1;
    chemistry.primordial_chemistry=2;
    chemistry.with_radiative_cooling=s->cooling;
    chemistry.metal_cooling=1;
    chemistry.UVbackground=0;
    chemistry.Gamma=5.0/3.0;
    chemistry.grackle_data_file=(char *)path;
    chemistry.omp_nthreads=1;
    chemistry.h2_on_dust=0;
    chemistry.dust_chemistry=0;
    chemistry.photoelectric_heating=0;
    chemistry.use_dust_density_field=0;
    chemistry.use_isrf_field=0;
    chemistry.use_radiative_transfer=0;
    chemistry.cmb_temperature_floor=0;
    chemistry.self_shielding_method=0;
    chemistry.H2_self_shielding=0;
    chemistry.exit_after_iterations_exceeded=1;
    chemistry.HydrogenFractionByMass=0.76;
    if(!local_initialize_chemistry_data(&chemistry,&rates,&units))return 2;

    gr_initialize_field_data(&fields);
    fields.grid_rank=1;fields.grid_dimension=dimensions;
    fields.grid_start=start;fields.grid_end=end;
    fields.density=&s->rho;fields.internal_energy=&s->energy;
    fields.HI_density=&s->species[0];fields.HII_density=&s->species[1];
    fields.HM_density=&s->species[2];fields.HeI_density=&s->species[3];
    fields.HeII_density=&s->species[4];fields.HeIII_density=&s->species[5];
    fields.H2I_density=&s->species[6];fields.H2II_density=&s->species[7];
    fields.e_density=&s->species[8];fields.metal_density=&s->metal;
    fields.H2_self_shielding_length=&s->length;
    const int ok=local_solve_chemistry(&chemistry,&rates,&units,&fields,dt);
    local_free_chemistry_data(&chemistry,&rates);
    return ok?0:3;
}

int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    const double a0=0.1,a1=0.2,rho0=1.6726219e-24,l0=1.0e21;
    const double rho_proper=1.0e-2*1.6726219e-24/0.76;
    const double time_unit=1.0e15,dt=2.0e-2;
    struct state bridge[2],direct[2];
    const double epochs[2]={a0,a1};
    double max_direct_difference=0.0;

    omp_set_num_threads(1);
    for(int j=0;j<2;j++)set_state(&bridge[j],epochs[j],rho_proper,0.76,0.5,1.0e5);
    memcpy(direct,bridge,sizeof(bridge));
    if(lag_grackle_neq_init(argv[1],0,0,0,1,rho0/(a0*a0*a0),l0*a0,
                            time_unit,a0,8.5)!=0)return 3;
    double vu[2],tu[2],proper_density[2];
    for(int j=0;j<2;j++) {
        const double a=epochs[j],dunit=rho0/(a*a*a),lunit=l0*a;
        code_units units={1,dunit,lunit,time_unit,lunit/(a*time_unit),1.0,a};
        vu[j]=get_velocity_units(&units);
        tu[j]=get_temperature_units(&units);
        proper_density[j]=bridge[j].rho*dunit;
        if(evolve(&bridge[j],a,rho0,l0,dt)!=0)return 4+j;
        if(evolve_direct(argv[1],&direct[j],a,rho0,l0,dt)!=0)return 6+j;
        max_direct_difference=fmax(max_direct_difference,
                                   max_fraction_difference(&bridge[j],&direct[j]));
    }
    const double rho_error=fabs(proper_density[0]-proper_density[1])/rho_proper;
    const double vu_error=fabs(vu[0]-vu[1])/fmax(vu[0],vu[1]);
    const double tu_error=fabs(tu[0]-tu[1])/fmax(tu[0],tu[1]);
    printf("cosmo_unit_contract epochs=%.3g/%.3g rho_proper_rel_diff=%.3e "
           "velocity_unit_rel_diff=%.3e temperature_unit_rel_diff=%.3e "
           "adapter_vs_direct_max_diff=%.3e cross_epoch_state_diff=%.3e\n",
           a0,a1,rho_error,vu_error,tu_error,max_direct_difference,
           max_fraction_difference(&bridge[0],&bridge[1]));
    if(rho_error>1e-12||vu_error>1e-12||tu_error>1e-12||
       max_direct_difference>1e-10)return 8;

    /* The adapter rejects a changing Grackle time unit after rates initialize. */
    struct state invalid;
    set_state(&invalid,a0,rho_proper,0.76,0.5,1.0e5);
    if(lag_grackle_neq_step(1,&invalid.rho,&invalid.energy,invalid.species,
                            &invalid.metal,&invalid.length,&invalid.cooling,
                            rho0/(a0*a0*a0),l0*a0,time_unit*1.01,a0,dt)==0)return 7;
    puts("Grackle cosmological unit contract: PASS");
    return 0;
}

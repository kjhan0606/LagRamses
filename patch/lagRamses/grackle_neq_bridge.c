/* Modern local Grackle API. No legacy global-state Grackle interface.
 * One initialized rate/context per OpenMP owner, transactional tile updates. */
#include <grackle.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

struct owner { chemistry_data chemistry; chemistry_data_storage rates;
    double *buffer; int capacity; };
static struct owner *owners;
static int nowners;
static code_units initial_units;
static char *data_path;
static int uv_requested;
static double uv_reion_redshift;

int lag_grackle_neq_init(const char *path,int uv,int shielding,int h2shield,
                        int comoving,double density,double length,double time,double a,
                        double z_reion)
{
    if(owners || sizeof(gr_float)!=sizeof(double) || !gr_check_consistency() ||
       !isfinite(z_reion) || z_reion<=-1.0)return 1;
    /* Grackle's comoving velocity unit is LU/(a TU), not LU/TU. */
    initial_units=(code_units){comoving,density,length,time,length/(a*time),1.0,a};
    uv_requested=uv!=0;
    uv_reion_redshift=z_reion;
    nowners=omp_get_max_threads();
    owners=calloc(nowners,sizeof(*owners));
    if(!owners)return 1;
    data_path=strdup(path);
    if(!data_path)return 1;
    for(int i=0;i<nowners;i++) {
        chemistry_data *c=&owners[i].chemistry;
        if(!local_initialize_chemistry_parameters(c))return 2;
        c->use_grackle=1;c->primordial_chemistry=2;
        c->with_radiative_cooling=1;c->metal_cooling=1;c->UVbackground=uv;
        c->grackle_data_file=data_path;c->Gamma=5.0/3.0;
        c->h2_on_dust=0;c->dust_chemistry=0;c->photoelectric_heating=0;
        c->use_dust_density_field=0;c->use_isrf_field=0;
        c->use_radiative_transfer=0;c->cmb_temperature_floor=0;
        c->self_shielding_method=shielding;c->H2_self_shielding=h2shield;
        c->omp_nthreads=1;
        c->exit_after_iterations_exceeded=1;
        if(!local_initialize_chemistry_data(c,&owners[i].rates,&initial_units))return 3;
    }
    fprintf(stderr,"Grackle NEQ: nine species, equilibrium metals; %d private owners; no dust/RT; UV=%d z_reion=%.8g\n",
            nowners,uv_requested,uv_reion_redshift);
    return 0;
}

static int budget(const double *s,double *h,double *he,double *e)
{
    for(int k=0;k<9;k++)if(!isfinite(s[k])||s[k]<0)return 1;
    *h=s[0]+s[1]+s[2]+s[6]+s[7];
    *he=s[3]+s[4]+s[5];
    *e=s[1]-s[2]+s[7]/2+s[4]/4+s[5]/2;
    return 0;
}

struct composition { double hydrogen_fraction; int cell,cooling; };
static int composition_order(const void *a,const void *b)
{
    const struct composition *x=a,*y=b;
    if(x->cooling!=y->cooling)return x->cooling-y->cooling;
    return (x->hydrogen_fraction>y->hydrogen_fraction)-(x->hydrogen_fraction<y->hydrogen_fraction);
}
int lag_grackle_neq_step(int n,double *density,double *energy,double *species,
                        double *metals,double *shield_length,int *cooling,
                        double dunit,double lunit,double tunit,double a,double dt)
{
    int id=omp_in_parallel()?omp_get_thread_num():0;
    if(!owners||id>=nowners||n<1||!isfinite(dt)||dt<0)return 1;
    if(!isfinite(dunit)||dunit<=0||!isfinite(lunit)||lunit<=0||
       !isfinite(tunit)||tunit<=0||!isfinite(a)||a<=0)return 1;
    if(fabs(tunit-initial_units.time_units)>
       1e-12*fmax(tunit,initial_units.time_units))return 1;
    struct owner *o=&owners[id];
    /* Match RAMSES' instantaneous z_reion switch. Keep the Grackle UV tables
     * initialized when requested, but suppress their rates until the configured
     * cosmological transition. Grackle recomputes UV rates from `a` at each
     * local_solve_chemistry call, so this per-owner toggle is safe for OMP. */
    o->chemistry.UVbackground = uv_requested &&
        (!initial_units.comoving_coordinates || (1.0/a-1.0)<uv_reion_redshift);
    if(o->capacity<n) {
        double *next=realloc(o->buffer,(size_t)n*24*sizeof(double));
        if(!next)return 1;
        o->buffer=next;o->capacity=n;
    }
    struct composition *order=malloc(n*sizeof(*order));
    if(!order)return 1;
    int status=0;
    double *b=o->buffer,*saved=b+13*n;
    memcpy(saved,energy,n*sizeof(double));
    memcpy(saved+n,species,9*n*sizeof(double));
    /* `tunit` is the fixed Grackle time unit; dt is already converted to it
     * by the RAMSES caller.  Recompute the epoch-dependent proper conversion
     * factors from the current scale factor. */
    code_units units={initial_units.comoving_coordinates,dunit,lunit,tunit,
                      lunit/(a*tunit),1.0,a};
    for(int i=0;i<n;i++) {
        double s[9],h,he,e;
        for(int k=0;k<9;k++)s[k]=species[k*n+i];
        if((cooling[i]!=0&&cooling[i]!=1)||!isfinite(shield_length[i])||shield_length[i]<=0||
           !isfinite(density[i])||density[i]<=0||!isfinite(energy[i])||energy[i]<=0||
           !isfinite(metals[i])||metals[i]<0||metals[i]>=density[i]||
           budget(s,&h,&he,&e)||h<=0||he<=0||fabs(s[8]-e)>1e-8*density[i]||
           fabs(h+he+metals[i]-density[i])>1e-8*density[i]){
            /* Failure-only diagnostics: never change admission tolerances or
             * repair a rejected state. Recompute even after short-circuiting. */
            h=he=e=NAN;
            int invalid_species=budget(s,&h,&he,&e);
            #pragma omp critical(grackle_input_failure_report)
            {
                fprintf(stderr,"Grackle NEQ input rejection: owner=%d cell=%d "
                        "rho=%.17g u=%.17g metal=%.17g length=%.17g cooling=%d "
                        "invalid_species=%d H=%.17g He=%.17g electron=%.17g "
                        "expected_electron=%.17g mass_residual=%.17g a=%.17g dt=%.17g\n",
                        id,i,density[i],energy[i],metals[i],shield_length[i],cooling[i],
                        invalid_species,h,he,s[8],e,h+he+metals[i]-density[i],a,dt);
                fprintf(stderr,"Grackle NEQ rejected species:");
                for(int k=0;k<9;k++)fprintf(stderr," %.17g",s[k]);
                fprintf(stderr,"\n");
                fflush(stderr);
            }
            status=2;goto finish;
        }
        order[i]=(struct composition){h/(h+he),i,cooling[i]};
    }
    qsort(order,n,sizeof(*order),composition_order);
    for(int first=0;first<n;) {
        int last=first+1;
        while(last<n && order[last].cooling==order[first].cooling &&
              order[last].hydrogen_fraction==order[first].hydrogen_fraction)last++;
        int count=last-first;
        for(int j=0;j<count;j++) {
            int i=order[first+j].cell;
            b[j]=density[i];b[n+j]=energy[i];
            for(int k=0;k<9;k++)b[(k+2)*n+j]=species[k*n+i];
            b[11*n+j]=metals[i];b[12*n+j]=shield_length[i];
        }
        int dimensions[3]={count,1,1},start[3]={0,0,0},end[3]={count-1,0,0};
        grackle_field_data f;
        gr_initialize_field_data(&f);
        f.grid_rank=1;f.grid_dimension=dimensions;f.grid_start=start;f.grid_end=end;
        f.density=b;f.internal_energy=b+n;
        f.HI_density=b+2*n;f.HII_density=b+3*n;f.HM_density=b+4*n;
        f.HeI_density=b+5*n;f.HeII_density=b+6*n;f.HeIII_density=b+7*n;
        f.H2I_density=b+8*n;f.H2II_density=b+9*n;f.e_density=b+10*n;
        f.metal_density=b+11*n;f.H2_self_shielding_length=b+12*n;
        o->chemistry.with_radiative_cooling=order[first].cooling;
        /* Grackle uses this scalar for nucleus normalization. Never force
         * enriched gas back to the primordial mixture or quantize its ratio. */
        o->chemistry.HydrogenFractionByMass=order[first].hydrogen_fraction;
        if(!local_solve_chemistry(&o->chemistry,&o->rates,&units,&f,dt)){status=3;goto finish;}
        for(int j=0;j<count;j++) {
            int i=order[first+j].cell;
            double s[9],old[9],h,he,e,oh,ohe,oe;
            for(int k=0;k<9;k++){s[k]=b[(k+2)*n+j];old[k]=species[k*n+i];}
            if(budget(s,&h,&he,&e)||budget(old,&oh,&ohe,&oe)||
               !isfinite(b[n+j])||b[n+j]<=0||
               fabs(h-oh)>1e-6*density[i]||fabs(he-ohe)>1e-6*density[i]||
               fabs(s[8]-e)>1e-6*density[i]){status=4;goto finish;}
            saved[i]=b[n+j];
            for(int k=0;k<9;k++)saved[(k+1)*n+i]=s[k];
        }
        first=last;
    }
    memcpy(energy,saved,n*sizeof(double));
    memcpy(species,saved+n,9*n*sizeof(double));
finish:
    free(order);
    return status;
}

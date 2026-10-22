#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <grackle.h>
int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,double,double,double,double,double);
/* Independent one-cell local API reference, with no adapter packing/sorting. */
static int reference(const char *path,double rho,double *energy,double *s,
                     double metals,double length,int cool,double dt)
{
    chemistry_data c;
    chemistry_data_storage rates={0};
    code_units units={0,1.6726219e-24,1,1,1,1,1};
    if(!local_initialize_chemistry_parameters(&c))return 1;
    c.use_grackle=1;c.primordial_chemistry=2;c.metal_cooling=1;
    c.with_radiative_cooling=cool;c.UVbackground=0;c.Gamma=5.0/3.0;
    c.grackle_data_file=(char *)path;c.omp_nthreads=1;
    c.h2_on_dust=0;c.dust_chemistry=0;c.photoelectric_heating=0;
    c.use_dust_density_field=0;c.use_isrf_field=0;c.use_radiative_transfer=0;
    c.cmb_temperature_floor=0;c.self_shielding_method=0;c.H2_self_shielding=0;
    c.exit_after_iterations_exceeded=1;
    double h=s[0]+s[1]+s[2]+s[6]+s[7],he=s[3]+s[4]+s[5];
    c.HydrogenFractionByMass=h/(h+he);
    /* The direct-reference contexts read HDF5 tables. Serialize that setup
     * just as the production adapter initializes its owners before the OMP
     * solve; the linked HDF5 library need not be thread-safe. */
    int initialized;
    #pragma omp critical(grackle_reference_table_init)
    initialized=local_initialize_chemistry_data(&c,&rates,&units);
    if(!initialized)return 2;
    int dims[3]={1,1,1},start[3]={0,0,0},end[3]={0,0,0};
    grackle_field_data f;
    gr_initialize_field_data(&f);
    f.grid_rank=1;f.grid_dimension=dims;f.grid_start=start;f.grid_end=end;
    f.density=&rho;f.internal_energy=energy;f.metal_density=&metals;
    f.HI_density=s;f.HII_density=s+1;f.HM_density=s+2;
    f.HeI_density=s+3;f.HeII_density=s+4;f.HeIII_density=s+5;
    f.H2I_density=s+6;f.H2II_density=s+7;f.e_density=s+8;
    f.H2_self_shielding_length=&length;
    int status=!local_solve_chemistry(&c,&rates,&units,&f,dt);
    local_free_chemistry_data(&c,&rates);
    return status;
}
int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    omp_set_num_threads(2);
    if(lag_grackle_neq_init(argv[1],0,0,0,0,1.6726219e-24,1,1,1,8.5))return 3;
    int failed=0;
    #pragma omp parallel for reduction(+:failed)
    for(int tile=0;tile<4;tile++) {
        double density[4],energy[4],s[36]={0},z[4],len[4],initial[4],before[36];int cooling[4];
        for(int i=0;i<4;i++) {
            density[i]=pow(10,i-2);z[i]=.001*density[i];len[i]=3e18;
            energy[i]=initial[i]=1e12;
            double x=(tile+i)%2?.68:.76;
            double h=x*(density[i]-z[i]),he=(1-x)*(density[i]-z[i]);
            s[i]=h*(1-2e-4-2e-6);s[4+i]=h*2e-4;s[12+i]=he;s[24+i]=h*2e-6;s[32+i]=s[4+i];
            cooling[i]=i%2;
        }
        memcpy(before,s,sizeof(s));
        int status=lag_grackle_neq_step(4,density,energy,s,z,len,cooling,1.6726219e-24,1,1,1,1e8);
        if(status){fprintf(stderr,"tile=%d solve_status=%d\n",tile,status);failed++;}
        for(int i=0;i<4;i++)if(!isfinite(energy[i])||energy[i]<=0||(!cooling[i]&&energy[i]!=initial[i]))failed++;
        for(int i=0;i<4;i++) {
            double direct[9],u=initial[i];
            for(int k=0;k<9;k++)direct[k]=before[4*k+i];
            if(reference(argv[1],density[i],&u,direct,z[i],len[i],cooling[i],1e8))failed++;
            if(fabs(u-energy[i])>1e-10*fmax(u,energy[i]))failed++;
            for(int k=0;k<9;k++)if(fabs(direct[k]-s[4*k+i])>1e-10*density[i])failed++;
            /* The split-step discrepancy is reported, not silently assumed zero. */
            double split[9],half=initial[i];
            for(int k=0;k<9;k++)split[k]=before[4*k+i];
            if(reference(argv[1],density[i],&half,split,z[i],len[i],cooling[i],5e7)||
               reference(argv[1],density[i],&half,split,z[i],len[i],cooling[i],5e7))failed++;
            #pragma omp critical
            printf("tile=%d cell=%d split_energy_relative=%.6e\n",tile,i,fabs(half-u)/u);
        }
        double saved[36];memcpy(saved,s,sizeof(s));double old[4];memcpy(old,energy,sizeof(old));
        s[0]=-1;double invalid[36];memcpy(invalid,s,sizeof(s));
        if(!lag_grackle_neq_step(4,density,energy,s,z,len,cooling,1.6726219e-24,1,1,1,1e8))failed++;
        if(memcmp(invalid,s,sizeof(s))||memcmp(old,energy,sizeof(old)))failed++;
    }
    printf("Grackle NEQ private-owner, mixed delayed-cooling and transactional tests: %s (%d)\n",failed?"FAIL":"PASS",failed);
    return failed?1:0;
}

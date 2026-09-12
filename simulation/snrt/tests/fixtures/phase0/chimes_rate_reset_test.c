/* Native regression for sparse active-species accumulator initialization.
 * No physical tables are loaded: reaction counts are zero. Pass a native
 * library path to compare the old and corrected accumulation contracts. */
#include <stdio.h>
#include <math.h>
#include <dlfcn.h>
#include "chimes_proto.h"
#include "chimes_vars.h"

int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    void *lib=dlopen(argv[1],RTLD_NOW|RTLD_LOCAL);
    if(!lib){fprintf(stderr,"%s\n",dlerror());return 2;}
    typedef void (*rate_fn)(struct Species_Structure *,struct gasVariables *,
                           const struct globalVariables *,struct UserData);
    rate_fn evaluate=(rate_fn)dlsym(lib,"update_rate_vector");
    if(!evaluate)return 2;
#define EMPTY_TABLE(name) do { \
    __typeof__(name) *table=dlsym(lib,#name); \
    if(!table)return 2; \
    table->N_reactions[0]=table->N_reactions[1]=0; \
} while(0)
    EMPTY_TABLE(chimes_table_T_dependent);
    EMPTY_TABLE(chimes_table_constant);
    EMPTY_TABLE(chimes_table_recombination_AB);
    EMPTY_TABLE(chimes_table_grain_recombination);
#undef EMPTY_TABLE
    struct globalVariables cfg={0};
    struct gasVariables gas={0};
    struct chimes_current_rates_struct rates={0};
    struct Species_Structure species[3]={0};
    struct UserData data={0};
    cfg.totalNumberOfSpecies=3;
    data.myGlobalVars=&cfg;data.myGasVars=&gas;
    data.species=species;data.chimes_current_rates=&rates;data.network_size=2;
    species[0].include_species=species[2].include_species=1;
    for(int repeat=0;repeat<2;repeat++){
        for(int s=0;s<3;s++){
            species[s].creation_rate=NAN;species[s].destruction_rate=NAN;
        }
        /* Old explicit caller cleared compact i=0,1, not active j=0,2. */
        for(int i=0;i<2;i++){
            species[i].creation_rate=0;species[i].destruction_rate=0;
        }
        evaluate(species,&gas,&cfg,data);
        for(int s=0;s<3;s++)if(species[s].creation_rate!=0 || species[s].destruction_rate!=0){
            fprintf(stderr,"RATE_RESET_FAIL species=%d repeat=%d\n",s,repeat);
            return 1;
        }
    }
    puts("CHIMES_RATE_RESET_PASS");
    dlclose(lib);return 0;
}

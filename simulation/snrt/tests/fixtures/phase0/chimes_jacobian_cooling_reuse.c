/* Dependency-exact cooling-coefficient and reaction-row reuse inside ONE
 * dense DQ Jacobian.
 * Experimental interposer, pinned to SUNDIALS 5.8. No rounded keys, altered
 * increments, omitted species, stale Jacobians or changed solver tolerances.
 * Cooling/heating summation remains evaluated for every perturbed abundance.
 * SNRT_CHIMES_JAC_REACTION_REUSE enables fixed-order reaction-row gathering;
 * unsupported secondary hooks retain the original complete reaction vector.
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include "cvode_impl.h"
#include "cvode_ls_impl.h"
#include <sunmatrix/sunmatrix_dense.h>
#include "chimes_proto.h"
#include "chimes_vars.h"
#if SUNDIALS_VERSION_MAJOR != 5 || SUNDIALS_VERSION_MINOR != 8
#error "This private-ABI experiment requires SUNDIALS 5.8"
#endif
typedef int (*jac_fn)(realtype,N_Vector,N_Vector,SUNMatrix,CVodeMem,N_Vector);
typedef void (*cool_fn)(struct gasVariables *,const struct globalVariables *,struct UserData);
typedef void (*vector_fn)(struct Species_Structure *,struct gasVariables *,const struct globalVariables *,struct UserData);
static _Thread_local jac_fn original_jac;
static _Thread_local cool_fn original_cool;
static _Thread_local vector_fn original_vector;
static _Thread_local int active,valid;
static _Thread_local struct chimes_current_rates_struct *owner;
static _Thread_local ChimesFloat temperature,ne,nHI,nHII;
static _Thread_local long calls,hits,misses,checked;
struct reaction {ChimesFloat *rate,previous;int r[3],p[3];};
static _Thread_local struct reaction *reactions;
static _Thread_local ChimesFloat **terms;
static _Thread_local int starts[2*CHIMES_TOTSIZE+1],nr,vector_valid;
static _Thread_local long rows_reused,rows_recomputed;
extern int snrt_chimes_jac_secondary_is_dark_noop(const struct globalVariables *) __attribute__((weak));
static void *symbol(const char *name)
{
    void *p=dlsym(RTLD_NEXT,name);
    if(!p){fprintf(stderr,"missing %s\n",name);abort();}
    return p;
}

static void append(int n,int arity,int products,const int *r,const int *p,ChimesFloat *rates)
{
    for(int i=0;i<n;i++) {
        struct reaction *q=&reactions[nr++];q->rate=&rates[i];
        for(int j=0;j<3;j++) {
            q->r[j]=j<arity?r[i*arity+j]:-1;
            q->p[j]=j<products?p[i*products+j]:-1;
        }
    }
}
static void build_dependencies(struct UserData *d)
{
    int m=d->mol_flag_index;
    int capacity=chimes_table_T_dependent.N_reactions[m]+chimes_table_constant.N_reactions[m]
        +chimes_table_recombination_AB.N_reactions[m]+chimes_table_grain_recombination.N_reactions[m]
        +chimes_table_cosmic_ray.N_reactions[m]+chimes_table_H2_collis_dissoc.N_reactions[m]
        +chimes_table_CO_cosmic_ray.N_reactions[m]+1;
    reactions=calloc(capacity,sizeof(*reactions));terms=malloc(6*capacity*sizeof(*terms));
    if(!reactions || !terms)abort();
    nr=0;vector_valid=0;memset(starts,0,sizeof(starts));
    struct chimes_current_rates_struct *c=d->chimes_current_rates;
#define ADD(table,a,b,rate) append(table.N_reactions[m],a,b,table.reactants,table.products,c->rate)
    ADD(chimes_table_T_dependent,3,3,T_dependent_rate);
    ADD(chimes_table_constant,2,3,constant_rate);
    ADD(chimes_table_recombination_AB,2,1,recombination_AB_rate);
    ADD(chimes_table_grain_recombination,2,1,grain_recombination_rate);
    if(d->myGasVars->cr_rate>0) {ADD(chimes_table_cosmic_ray,1,3,cosmic_ray_rate);}
    if(m) {
        append(1,2,1,chimes_table_H2_dust_formation.reactants,chimes_table_H2_dust_formation.products,
               &c->H2_dust_formation_rate);
        ADD(chimes_table_H2_collis_dissoc,2,3,H2_collis_dissoc_rate);
        if(d->myGasVars->cr_rate>0) {ADD(chimes_table_CO_cosmic_ray,1,2,CO_cosmic_ray_rate);}
    }
#undef ADD
    /* Preserve the original reaction/addition order separately for each
     * species' creation and destruction sum. No delta subtraction, regrouping
     * or cancellation of repeated reactants is allowed. */
    for(int k=0;k<nr;k++)for(int side=0;side<2;side++)for(int j=0;j<3;j++) {
        int species=side?reactions[k].p[j]:reactions[k].r[j];
        if(species<0)break;
        if(species>=d->myGlobalVars->totalNumberOfSpecies)abort();
        starts[side*CHIMES_TOTSIZE+species+1]++;
    }
    for(int i=1;i<=2*CHIMES_TOTSIZE;i++)starts[i]+=starts[i-1];
    int cursor[2*CHIMES_TOTSIZE];memcpy(cursor,starts,sizeof(cursor));
    for(int k=0;k<nr;k++)for(int side=0;side<2;side++)for(int j=0;j<3;j++) {
        int species=side?reactions[k].p[j]:reactions[k].r[j];
        if(species<0)break;
        terms[cursor[side*CHIMES_TOTSIZE+species]++]=reactions[k].rate;
    }
}

void update_rate_vector(struct Species_Structure *species,struct gasVariables *g,
                        const struct globalVariables *c,struct UserData d)
{
    if(!original_vector)original_vector=(vector_fn)symbol("update_rate_vector");
    if(!active || !reactions){original_vector(species,g,c,d);return;}
    unsigned char dirty[CHIMES_TOTSIZE]={0};
    for(int k=0;k<nr;k++) {
        struct reaction *q=&reactions[k];
        if(!vector_valid || memcmp(q->rate,&q->previous,sizeof(q->previous))) {
            for(int side=0;side<2;side++)for(int j=0;j<3;j++) {
                int s=side?q->p[j]:q->r[j];if(s<0)break;dirty[s]=1;
            }
            q->previous=*q->rate;
        }
    }
    if(!vector_valid) {
        original_vector(species,g,c,d);vector_valid=1;return;
    }
    for(int s=0;s<c->totalNumberOfSpecies;s++) {
        if(!dirty[s]){rows_reused++;continue;}
        rows_recomputed++;
        ChimesFloat destruction=0,creation=0;
        for(int k=starts[s];k<starts[s+1];k++)destruction+=*terms[k];
        for(int k=starts[CHIMES_TOTSIZE+s];k<starts[CHIMES_TOTSIZE+s+1];k++)creation+=*terms[k];
        species[s].creation_rate=creation;species[s].destruction_rate=destruction;
    }
}

void update_cooling_rates(struct gasVariables *g,const struct globalVariables *c,struct UserData d)
{
    if(!original_cool)original_cool=(cool_fn)symbol("update_cooling_rates");
    if(!active){original_cool(g,c,d);return;}
    /* These are ALL evolving inputs to upstream update_cooling_rates:
     * its 1D, 2D and 4D interpolation branches, including the high-T ne
     * multiplier. Tables and configuration are immutable within this call.
     * The live coefficient buffer is written only by this routine on the
     * admitted CPU path. No reuse across solves or buffer allocations. */
    const ChimesFloat e=g->abundances[c->speciesIndices[sp_elec]]*g->nH_tot;
    const ChimesFloat h=g->abundances[c->speciesIndices[sp_HI]]*g->nH_tot;
    const ChimesFloat hp=g->abundances[c->speciesIndices[sp_HII]]*g->nH_tot;
    if(valid && owner==d.chimes_current_rates && temperature==g->temperature &&
       ne==e && nHI==h && nHII==hp) {hits++;return;}
    original_cool(g,c,d);misses++;
    owner=d.chimes_current_rates;temperature=g->temperature;ne=e;nHI=h;nHII=hp;valid=1;
}

int cvLsDenseDQJac(realtype t,N_Vector y,N_Vector fy,SUNMatrix J,CVodeMem mem,N_Vector tmp)
{
    if(!original_jac)original_jac=(jac_fn)symbol("cvLsDenseDQJac");
    struct UserData *d=mem->cv_user_data;
    if(d->myGlobalVars->N_spectra || getenv("SNRT_CHIMES_RHS_BACKEND") ||
       getenv("SNRT_CHIMES_CVODE_BACKEND"))return original_jac(t,y,fy,J,mem,tmp);
    calls++;active=1;valid=0;
    if(getenv("SNRT_CHIMES_JAC_REACTION_REUSE") &&
       (!snrt_chimes_secondary_rates || !d->myGlobalVars->rt_update_flux ||
        (snrt_chimes_jac_secondary_is_dark_noop &&
         snrt_chimes_jac_secondary_is_dark_noop(d->myGlobalVars))))build_dependencies(d);
    int rc=original_jac(t,y,fy,J,mem,tmp);
    active=0;valid=0;
    free(reactions);free(terms);reactions=NULL;terms=NULL;
    if(!rc && getenv("SNRT_CHIMES_JAC_VERIFY")) {
        SUNMatrix ref=SUNMatClone(J);if(!ref)abort();
        int reference_status=original_jac(t,y,fy,ref,mem,tmp);
        sunindextype n=SUNDenseMatrix_Columns(J);
        if(reference_status || memcmp(SUNDenseMatrix_Data(J),SUNDenseMatrix_Data(ref),sizeof(realtype)*n*n)) {
            fprintf(stderr,"JACOBIAN_COOLING_REUSE_MATRIX_MISMATCH\n");abort();
        }
        SUNMatDestroy(ref);checked++;
    }
    return rc;
}
/* Main-thread counters only; OMP endpoint parity is checked by the fixture. */
__attribute__((destructor)) static void report(void)
{
    if(calls)fprintf(stderr,"JACOBIAN_COOLING_REUSE calls=%ld hits=%ld misses=%ld matrices_checked=%ld rows_reused=%ld rows_recomputed=%ld\n",
                     calls,hits,misses,checked,rows_reused,rows_recomputed);
}

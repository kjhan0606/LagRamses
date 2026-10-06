/* Experimental LD_PRELOAD adapter, pinned to SUNDIALS 5.8 internals.
 * Not a production ABI. Retains upstream DQ increments, dense LU, tolerances
 * and Jacobian refresh policy. Never combine with the RHS cost interposer.
 * Parallelism is within one cell only when no outer OMP team is active.
 * Reference increment formula: SUNDIALS 5.8 src/cvode/cvode_ls.c (BSD-3-Clause).
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <dlfcn.h>
#include <omp.h>
#include "cvode_impl.h"
#include "cvode_ls_impl.h"
#include <sunmatrix/sunmatrix_dense.h>
#include "chimes_proto.h"
#include "chimes_vars.h"
#if SUNDIALS_VERSION_MAJOR != 5 || SUNDIALS_VERSION_MINOR != 8
#error "This private-ABI experiment requires SUNDIALS 5.8"
#endif

extern void snrt_chimes_rate_cache_invalidate(void);
extern int snrt_chimes_jac_context_clone(const struct globalVariables *,struct globalVariables *);
extern void snrt_chimes_jac_context_release(const struct globalVariables *,struct globalVariables *);
extern int snrt_chimes_jac_context_flag(const struct globalVariables *);
typedef int (*original_fn)(realtype,N_Vector,N_Vector,SUNMatrix,CVodeMem,N_Vector);
static _Thread_local original_fn original;
static _Thread_local long parallel_calls, fallback_calls, checked_calls;
static _Thread_local int max_workers;

int cvLsDenseDQJac(realtype t,N_Vector y,N_Vector fy,SUNMatrix J,
                  CVodeMem mem,N_Vector tmp)
{
    if(!original){original=(original_fn)dlsym(RTLD_NEXT,"cvLsDenseDQJac");if(!original)abort();}
    struct UserData *d=mem->cv_user_data;
    const char *choice=getenv("SNRT_CHIMES_JAC_THREADS");
    char *end=NULL;long requested=choice?strtol(choice,&end,10):0;
    if(choice && (!*choice || *end || requested<1 || requested>8))abort();
    if(!requested || omp_in_parallel() || d->myGlobalVars->N_spectra ||
       d->myGasVars->hybrid_data ||
       getenv("SNRT_CHIMES_RHS_BACKEND") || getenv("SNRT_CHIMES_CVODE_BACKEND")) {
        fallback_calls++;return original(t,y,fy,J,mem,tmp);
    }
    const sunindextype n=SUNDenseMatrix_Columns(J);
    if(n>CHIMES_TOTSIZE+1 || d->myGlobalVars->totalNumberOfSpecies>CHIMES_TOTSIZE)abort();
    realtype inc[CHIMES_TOTSIZE+1];
    int status[CHIMES_TOTSIZE+1]={0},flags[8]={0};
    const realtype norm=N_VWrmsNorm(fy,mem->cv_ewt);
    const realtype floor=norm!=0 ? 1000.0*fabs(mem->cv_h)*mem->cv_uround*n*norm : 1.0;
    const realtype srur=sqrt(mem->cv_uround);
    for(sunindextype j=0;j<n;j++) {
        inc[j]=fmax(srur*fabs(NV_Ith_S(y,j)),floor/NV_Ith_S(mem->cv_ewt,j));
        if(mem->cv_constraintsSet) {
            realtype c=NV_Ith_S(mem->cv_constraints,j),v=NV_Ith_S(y,j)+inc[j];
            if((fabs(c)==1 && v*c<0)||(fabs(c)==2 && v*c<=0))inc[j]=-inc[j];
        }
    }
    int buffer_size;determine_current_rates_buffer_size(&buffer_size,d->myGlobalVars);
    parallel_calls++;
    int actual_workers=0;
#pragma omp parallel num_threads(requested)
    {
#pragma omp single
        actual_workers=omp_get_num_threads();
        struct UserData local=*d;
        struct gasVariables gas=*d->myGasVars;
        struct globalVariables global;
        struct chimes_current_rates_struct rates=*d->chimes_current_rates;
        struct Species_Structure species[CHIMES_TOTSIZE];
        ChimesFloat abundances[CHIMES_TOTSIZE];
        memcpy(species,d->species,sizeof(*species)*d->myGlobalVars->totalNumberOfSpecies);
        memcpy(abundances,gas.abundances,sizeof(*abundances)*d->myGlobalVars->totalNumberOfSpecies);
        if(snrt_chimes_jac_context_clone(d->myGlobalVars,&global))abort();
        allocate_current_rates_memory(&rates,&global);
        if(!rates.data_buffer)abort();
        memcpy(rates.data_buffer,d->chimes_current_rates->data_buffer,sizeof(ChimesFloat)*buffer_size);
        gas.abundances=abundances;local.myGasVars=&gas;local.myGlobalVars=&global;
        local.species=species;local.chimes_current_rates=&rates;
        N_Vector yp=N_VClone(y),fp=N_VClone(y),column=N_VCloneEmpty(y);
        if(!yp || !fp || !column)abort();
        N_VScale(1,y,yp);
        snrt_chimes_rate_cache_invalidate();
#pragma omp for schedule(static)
        for(sunindextype j=0;j<n;j++) {
            NV_Ith_S(yp,j)=NV_Ith_S(y,j)+inc[j];
            status[j]=mem->cv_f(t,yp,fp,&local);
            NV_Ith_S(yp,j)=NV_Ith_S(y,j);
            if(!status[j]) {
                N_VSetArrayPointer(SUNDenseMatrix_Column(J,j),column);
                const realtype inv=1.0/inc[j];
                N_VLinearSum(inv,fp,-inv,fy,column);
            }
        }
        flags[omp_get_thread_num()]=snrt_chimes_jac_context_flag(&global);
        snrt_chimes_rate_cache_invalidate();
        N_VSetArrayPointer(NULL,column);N_VDestroy(column);N_VDestroy(fp);N_VDestroy(yp);
        free_current_rates_memory(&rates,&global);
        snrt_chimes_jac_context_release(d->myGlobalVars,&global);
    }
    if(actual_workers>max_workers)max_workers=actual_workers;
    ((CVLsMem)mem->cv_lmem)->nfeDQ+=n;
    int retry=0;for(int j=0;j<n;j++)retry|=status[j]!=0;
    for(int j=0;j<requested;j++)retry|=flags[j];
    if(retry){fallback_calls++;return original(t,y,fy,J,mem,tmp);}
    if(getenv("SNRT_CHIMES_JAC_VERIFY")) {
        SUNMatrix ref=SUNMatClone(J);if(!ref)abort();
        int rc=original(t,y,fy,ref,mem,tmp);
        if(rc || memcmp(SUNDenseMatrix_Data(ref),SUNDenseMatrix_Data(J),sizeof(realtype)*n*n)) {
            fprintf(stderr,"PARALLEL_JACOBIAN_MATRIX_MISMATCH\n");abort();
        }
        SUNMatDestroy(ref);checked_calls++;return 0;
    }
    /* Reproduce the upstream final-column mutable gas/rate state. */
    const realtype saved=NV_Ith_S(y,n-1);
    NV_Ith_S(y,n-1)=saved+inc[n-1];
    int rc=mem->cv_f(t,y,tmp,d);
    NV_Ith_S(y,n-1)=saved;
    ((CVLsMem)mem->cv_lmem)->nfeDQ++;
    return rc;
}
__attribute__((destructor)) static void report(void)
{
    if(parallel_calls || fallback_calls)
        fprintf(stderr,"PARALLEL_JACOBIAN calls=%ld fallback=%ld matrices_checked=%ld max_workers=%d\n",
                parallel_calls,fallback_calls,checked_calls,max_workers);
}

/* Diagnostic LD_PRELOAD only: no production counters or solver changes.
 * Run a serial native probe. Aggregate per-thread counters are printed only
 * for the main thread at exit; do not interpret this as an OMP profiler. */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <dlfcn.h>
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "chimes_proto.h"
#include "chimes_vars.h"
static _Thread_local CVRhsFn original_rhs;
static _Thread_local double network_s,cvode_s,rhs_s,lu_s,solve_s;
static _Thread_local long calls,integrations,rhs_calls,lu_calls,solve_calls;
static _Thread_local long steps,nfe,nfeLS,nje,netfails,nconvfails;
static _Thread_local double coefficients_s,rates_s,vector_s,cooling_s;
/* Optional exact-value cache experiment; never enabled in production. */
static _Thread_local int cache_active,cache_valid;
static _Thread_local double cache_temperature;
static _Thread_local long cache_hits;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void *symbol(const char *name){void *p=dlsym(RTLD_NEXT,name);if(!p){fprintf(stderr,"missing %s\n",name);abort();}return p;}
static int timed_rhs(realtype t,N_Vector y,N_Vector dy,void *data){
    double start=now();int rc=original_rhs(t,y,dy,data);rhs_s+=now()-start;rhs_calls++;return rc;
}
int CVodeInit(void *mem,CVRhsFn rhs,realtype t,N_Vector y){
    static _Thread_local int (*fn)(void *,CVRhsFn,realtype,N_Vector);
    if(!fn)fn=symbol("CVodeInit");original_rhs=rhs;return fn(mem,timed_rhs,t,y);
}
int CVode(void *mem,realtype end,N_Vector y,realtype *t,int task){
    static _Thread_local int (*fn)(void *,realtype,N_Vector,realtype *,int);
    if(!fn)fn=symbol("CVode");double start=now();int rc=fn(mem,end,y,t,task);cvode_s+=now()-start;integrations++;
    long v=0;
#define COUNT(api,dest) do {v=0;if(api(mem,&v)==0)dest+=v;}while(0)
    COUNT(CVodeGetNumSteps,steps);COUNT(CVodeGetNumRhsEvals,nfe);
    COUNT(CVodeGetNumLinRhsEvals,nfeLS);COUNT(CVodeGetNumJacEvals,nje);
    COUNT(CVodeGetNumErrTestFails,netfails);COUNT(CVodeGetNumNonlinSolvConvFails,nconvfails);
#undef COUNT
    return rc;
}
int SUNLinSolSetup_Dense(SUNLinearSolver s,SUNMatrix a){
    static _Thread_local int (*fn)(SUNLinearSolver,SUNMatrix);
    if(!fn)fn=symbol("SUNLinSolSetup_Dense");double start=now();int rc=fn(s,a);lu_s+=now()-start;lu_calls++;return rc;
}
int SUNLinSolSolve_Dense(SUNLinearSolver s,SUNMatrix a,N_Vector x,N_Vector b,realtype tol){
    static _Thread_local int (*fn)(SUNLinearSolver,SUNMatrix,N_Vector,N_Vector,realtype);
    if(!fn)fn=symbol("SUNLinSolSolve_Dense");double start=now();int rc=fn(s,a,x,b,tol);solve_s+=now()-start;solve_calls++;return rc;
}
int chimes_network(struct gasVariables *g,const struct globalVariables *c){
    static _Thread_local int (*fn)(struct gasVariables *,const struct globalVariables *);
    if(!fn)fn=symbol("chimes_network");
    cache_valid=0;cache_active=getenv("SNRT_COST_EXACT_CACHE")!=NULL;
    double start=now();int rc=fn(g,c);network_s+=now()-start;calls++;
    cache_active=0;return rc;
}
__attribute__((destructor)) static void report(void){
    fprintf(stderr,"CHIMES_EXACT_CACHE hits=%ld\n",cache_hits);
    fprintf(stderr,"CHIMES_RHS_COST coefficients_s=%.9f rates_s=%.9f vector_s=%.9f cooling_s=%.9f\n",
        coefficients_s,rates_s,vector_s,cooling_s);
    fprintf(stderr,"CHIMES_COST calls=%ld integrations=%ld network_s=%.9f cvode_s=%.9f rhs_s=%.9f lu_s=%.9f solve_s=%.9f rhs_calls=%ld lu_calls=%ld solve_calls=%ld steps=%ld nfe=%ld nfeLS=%ld nje=%ld errorfails=%ld convfails=%ld\n",
        calls,integrations,network_s,cvode_s,rhs_s,lu_s,solve_s,rhs_calls,lu_calls,solve_calls,steps,nfe,nfeLS,nje,netfails,nconvfails);
}

void update_rate_coefficients(struct gasVariables *g,const struct globalVariables *c,struct UserData d,int mode){
    static _Thread_local void (*fn)(struct gasVariables *,const struct globalVariables *,struct UserData,int);
    if(!fn)fn=symbol("update_rate_coefficients");
    int reuse=cache_active && c->N_spectra==0 && mode==1 && cache_valid && g->temperature==cache_temperature;
    double start=now();fn(g,c,d,reuse?0:mode);coefficients_s+=now()-start;
    if(reuse)cache_hits++;
    if(cache_active && c->N_spectra==0 && mode==1){cache_temperature=g->temperature;cache_valid=1;}
}
void update_rates(struct gasVariables *g,const struct globalVariables *c,struct UserData d){
    static _Thread_local void (*fn)(struct gasVariables *,const struct globalVariables *,struct UserData);
    if(!fn)fn=symbol("update_rates");double start=now();fn(g,c,d);rates_s+=now()-start;
}
void update_rate_vector(struct Species_Structure *s,struct gasVariables *g,const struct globalVariables *c,struct UserData d){
    static _Thread_local void (*fn)(struct Species_Structure *,struct gasVariables *,const struct globalVariables *,struct UserData);
    if(!fn)fn=symbol("update_rate_vector");double start=now();fn(s,g,c,d);vector_s+=now()-start;
}
ChimesFloat calculate_total_cooling_rate(struct gasVariables *g,const struct globalVariables *c,struct UserData d,int mode){
    static _Thread_local ChimesFloat (*fn)(struct gasVariables *,const struct globalVariables *,struct UserData,int);
    if(!fn)fn=symbol("calculate_total_cooling_rate");double start=now();ChimesFloat result=fn(g,c,d,mode);cooling_s+=now()-start;return result;
}

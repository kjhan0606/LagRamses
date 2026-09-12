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
    if(!fn)fn=symbol("chimes_network");double start=now();int rc=fn(g,c);network_s+=now()-start;calls++;return rc;
}
__attribute__((destructor)) static void report(void){
    fprintf(stderr,"CHIMES_COST calls=%ld integrations=%ld network_s=%.9f cvode_s=%.9f rhs_s=%.9f lu_s=%.9f solve_s=%.9f rhs_calls=%ld lu_calls=%ld solve_calls=%ld steps=%ld nfe=%ld nfeLS=%ld nje=%ld errorfails=%ld convfails=%ld\n",
        calls,integrations,network_s,cvode_s,rhs_s,lu_s,solve_s,rhs_calls,lu_calls,solve_calls,steps,nfe,nfeLS,nje,netfails,nconvfails);
}

/* Diagnostic LD_PRELOAD only: no production counters or solver changes.
 * Run a serial native probe. Aggregate per-thread counters are printed only
 * for the main thread at exit; do not interpret this as an OMP profiler. */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <dlfcn.h>
#include <math.h>
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "chimes_proto.h"
#include "chimes_vars.h"
static _Thread_local CVRhsFn original_rhs;
static _Thread_local struct UserData *current_data;
static _Thread_local double network_s,cvode_s,rhs_s,lu_s,solve_s;
static _Thread_local long calls,integrations,rhs_calls,lu_calls,solve_calls;
static _Thread_local long steps,nfe,nfeLS,nje,netfails,nconvfails;
static _Thread_local double coefficients_s,rates_s,vector_s,cooling_s;
/* Optional exact-value cache experiment; never enabled in production. */
static _Thread_local int cache_active,cache_valid;
static _Thread_local double cache_temperature;
static _Thread_local long cache_hits;
static _Thread_local int jac_age_active;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void *symbol(const char *name){void *p=dlsym(RTLD_NEXT,name);if(!p){fprintf(stderr,"missing %s\n",name);abort();}return p;}
/* Diagnostic convergence check only: may tighten, never loosen, tolerances.
 * Disabled unless explicitly requested; no production configuration change. */
static double tolerance_scale(void){
    const char *value=getenv("SNRT_COST_TOLERANCE_SCALE");
    if(!value)return 1;
    char *end=NULL;double scale=strtod(value,&end);
    if(!*value || *end || !isfinite(scale) || scale<=0 || scale>1)abort();
    return scale;
}
int CVodeSStolerances(void *mem,realtype rel,realtype abs){
    static _Thread_local int (*fn)(void *,realtype,realtype);
    if(!fn)fn=symbol("CVodeSStolerances");
    double scale=tolerance_scale();return fn(mem,rel*scale,abs*scale);
}
int CVodeSVtolerances(void *mem,realtype rel,N_Vector abs){
    static _Thread_local int (*fn)(void *,realtype,N_Vector);
    if(!fn)fn=symbol("CVodeSVtolerances");
    double scale=tolerance_scale();
    const char *molecular=getenv("SNRT_COST_MOLECULAR_ATOL_SCALE");
    if(scale==1 && !molecular)return fn(mem,rel,abs);
    N_Vector tighter=N_VClone(abs);if(!tighter)abort();
    N_VScale(scale,abs,tighter);
    if(molecular) {
        char *end=NULL;double factor=strtod(molecular,&end);
        if(!*molecular || *end || !isfinite(factor) || factor<=0 || factor>1 || !current_data)abort();
        int j=0;
        for(int s=0;s<current_data->myGlobalVars->totalNumberOfSpecies;s++)if(current_data->species[s].include_species) {
            if(s>=sp_H2)N_VGetArrayPointer(tighter)[j]*=factor;
            j++;
        }
    }
    int rc=fn(mem,rel*scale,tighter);
    N_VDestroy(tighter);return rc;
}
/* Bounded solver experiment: retain CVODE error/convergence controls while
 * varying only the maximum Jacobian age. No setting is changed by default. */
int CVodeSetLinearSolver(void *mem,SUNLinearSolver solver,SUNMatrix matrix){
    static _Thread_local int (*fn)(void *,SUNLinearSolver,SUNMatrix);
    if(!fn)fn=symbol("CVodeSetLinearSolver");
    int rc=fn(mem,solver,matrix);
    const char *age=getenv("SNRT_COST_JACOBIAN_AGE");
    if(!rc && age && jac_age_active){
        char *end=NULL;long value=strtol(age,&end,10);
        if(!*age || *end || value<1 || value>1000)abort();
        rc=CVodeSetJacEvalFrequency(mem,value);
    }
    return rc;
}
static int timed_rhs(realtype t,N_Vector y,N_Vector dy,void *data){
    current_data=data;
    double start=now();int rc=original_rhs(t,y,dy,data);rhs_s+=now()-start;rhs_calls++;return rc;
}
int CVodeInit(void *mem,CVRhsFn rhs,realtype t,N_Vector y){
    static _Thread_local int (*fn)(void *,CVRhsFn,realtype,N_Vector);
    if(!fn)fn=symbol("CVodeInit");original_rhs=rhs;int rc=fn(mem,timed_rhs,t,y);
    if(!rc && getenv("SNRT_COST_NONNEGATIVE") && current_data) {
        N_Vector constraints=N_VClone(y);if(!constraints)abort();N_VConst(0,constraints);
        for(int i=0;i<current_data->network_size;i++)N_VGetArrayPointer(constraints)[i]=1;
        rc=CVodeSetConstraints(mem,constraints);N_VDestroy(constraints);
    }
    return rc;
}
int CVodeSetUserData(void *mem,void *data){
    static _Thread_local int (*fn)(void *,void *);
    if(!fn)fn=symbol("CVodeSetUserData");current_data=data;return fn(mem,data);
}
int CVode(void *mem,realtype end,N_Vector y,realtype *t,int task){
    static _Thread_local int (*fn)(void *,realtype,N_Vector,realtype *,int);
    if(!fn)fn=symbol("CVode");double start=now();int rc=fn(mem,end,y,t,task);cvode_s+=now()-start;integrations++;
    if(getenv("SNRT_COST_RESTART_ERROR")) {
        for(int retry=0;retry<4 && rc==CV_ERR_FAILURE && *t<end;retry++) {
            realtype previous=*t;
            /* CVODE returns its last accepted state on this failure. Reset
             * multistep history only, not abundances, rates or tolerances. */
            if(CVodeReInit(mem,previous,y))break;
            double restart=now();rc=fn(mem,end,y,t,task);cvode_s+=now()-restart;
            fprintf(stderr,"CHIMES_HISTORY_RESTART retry=%d from=%.17g to=%.17g flag=%d\n",retry+1,(double)previous,(double)*t,rc);
            if(*t<=previous)break;
        }
    }
    if(rc<0 && getenv("SNRT_COST_FAILURE_DETAIL")) {
        N_Vector err=N_VClone(y),weight=N_VClone(y);
        if(!err || !weight)abort();
        realtype h=0;int worst=-1;double maximum=0;
        CVodeGetCurrentStep(mem,&h);
        if(!CVodeGetEstLocalErrors(mem,err) && !CVodeGetErrWeights(mem,weight)) {
            int j=0;
            for(int s=0;s<current_data->myGlobalVars->totalNumberOfSpecies;s++)if(current_data->species[s].include_species) {
                double v=fabs(N_VGetArrayPointer(err)[j]*N_VGetArrayPointer(weight)[j]);j++;
                if(v>maximum){maximum=v;worst=s;}
            }
        }
        fprintf(stderr,"CHIMES_FAILURE_DETAIL flag=%d t=%.17g end=%.17g h=%.17g worst_species=%d weighted_error=%.17g\n",rc,(double)*t,(double)end,(double)h,worst,maximum);
        N_VDestroy(err);N_VDestroy(weight);
    }
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
    int previous_jac_age_active=jac_age_active;
    jac_age_active=!getenv("SNRT_COST_JACOBIAN_COLD_ONLY") ||
        (c->N_spectra==0 && g->temperature<=200 && g->nH_tot<=0.2);
    double start=now();int rc=fn(g,c);network_s+=now()-start;calls++;
    jac_age_active=previous_jac_age_active;
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

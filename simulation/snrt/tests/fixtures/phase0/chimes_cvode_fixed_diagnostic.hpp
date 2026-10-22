// Diagnostic-only inclusion after Cell, tile_cell and compare are defined.
// The executable exports CVode to interpose the dynamically linked CHIMES
// call. It delegates with unchanged arguments and never advances the solver
// itself. This reports terminal statistics, not the first divergent step.
#include <cvode/cvode.h>
static bool fixed_diag_active=false;
static int fixed_backend=0, fixed_repeat=0;
static double fixed_tolerance_factor=1;
static unsigned long fixed_solver_calls=0;

extern "C" int CVode(void *mem,realtype tout,N_Vector y,realtype *tret,int task) {
    using Fn=int (*)(void*,realtype,N_Vector,realtype*,int);
    static const Fn real_cvode=reinterpret_cast<Fn>(dlsym(RTLD_NEXT,"CVode"));
    require(real_cvode!=nullptr,"resolve actual CVode");
    const int flag=real_cvode(mem,tout,y,tret,task);
    if(!fixed_diag_active)return flag;
    long steps=-1,rhs=-1,errfails=-1,convfails=-1;
    int order=-1;
    realtype last_h=NAN,next_h=NAN,current_t=NAN;
    int queries=CVodeGetNumSteps(mem,&steps);
    queries|=CVodeGetNumRhsEvals(mem,&rhs);
    queries|=CVodeGetNumErrTestFails(mem,&errfails);
    queries|=CVodeGetNumNonlinSolvConvFails(mem,&convfails);
    queries|=CVodeGetLastOrder(mem,&order);
    queries|=CVodeGetLastStep(mem,&last_h);
    queries|=CVodeGetCurrentStep(mem,&next_h);
    queries|=CVodeGetCurrentTime(mem,&current_t);
    ++fixed_solver_calls; // Only the selected worker runs cells in this mode.
    std::printf("FIXED_CVODE worker=%d gpu=%d repeat=%d tol_factor=%.17g flag=%d requested_t=%.17g returned_t=%.17g internal_t=%.17g steps=%ld rhs=%ld error_failures=%ld convergence_failures=%ld last_order=%d last_h=%.17g next_h=%.17g stats_status=%d\n",
        omp_get_thread_num(),fixed_backend,fixed_repeat,fixed_tolerance_factor,
        flag,(double)tout,(double)*tret,(double)current_t,steps,rhs,errfails,
        convfails,order,(double)last_h,(double)next_h,queries);
    return flag;
}

static void fixed_route_diagnostic(){
    require(diagnostic_mode,"fixed route requires diagnostic mode");
    require(nc>0 && nc<=128,"fixed route bounded tile");
    make_cells();
    for(int i=0;i<nc;i++){
        require(fixture_cell_index[i]%12==5,"fixed route logical input");
        require(!std::memcmp(input[i].controls,input[0].controls,sizeof(input[0].controls)) &&
                !std::memcmp(input[i].elements,input[0].elements,sizeof(input[0].elements)) &&
                !std::memcmp(input[i].old,input[0].old,sizeof(input[0].old)),
                "fixed route repeated identical state");
    }
    mode("cpu"); // CPU CVODE + dense LU in both arms; only RHS differs.
    setenv("SNRT_CHIMES_GPU_WORKER","any",1);
    require(!std::getenv("SNRT_CHIMES_RHS_VERIFY"),"fixed route uses actual backend results");
    cuda_pool_init(gpu_local_rank(),3);
    fixed_diag_active=true;
    Cell baseline[2][nc],previous[2][nc],trial[nc];
    int factor_index=0;
    for(const char *factor:{"1","0.1","0.01"}){
        setenv("SNRT_CHIMES_TEST_TOLERANCE_FACTOR",factor,1);
        fixed_tolerance_factor=std::strtod(factor,nullptr);
        for(int owner:{0,3})for(int repeat=0;repeat<2;repeat++)for(int backend:{0,1}){
            fixed_backend=backend;fixed_repeat=repeat;
            require(!snrt_chimes_rhs_configure_mode(backend),"fixed backend configure");
            std::memcpy(trial,input,sizeof(input));
            unsigned long long before[7],after[7];snrt_chimes_rhs_counts(before);
            const unsigned long calls_before=fixed_solver_calls;
            int dispatch=-1,team=0;
#pragma omp parallel num_threads(8)
            {
                if(omp_get_thread_num()==owner){
                    team=omp_get_num_threads();
                    dispatch=snrt_chimes_rhs_tile_run(nc,tile_cell,trial);
                }
            }
            require(team==8 && dispatch==0,"fixed worker tile completed");
            snrt_chimes_rhs_counts(after);
            require(fixed_solver_calls>calls_before,"CVode interposition exercised");
            require(after[3]==before[3],"fixed route device errors");
            require(!backend || after[0]>before[0],"fixed GPU route performed GPU work");
            require(backend || after[0]==before[0],"fixed CPU route stayed CPU");
            int accepted=0;
            for(int i=0;i<nc;i++){
                const auto &c=trial[i];
                require(c.tile_worker==owner && c.tile_gpu_lease==backend,"fixed route ownership");
                if(c.status){
                    require(c.temperature==input[i].controls[1] &&
                        !std::memcmp(c.next,input[i].old,sizeof(c.next)),
                        "rejected state restored");
                    continue;
                }
                ++accepted;
                double elements[11],q=0;
                require(std::isfinite(c.temperature) && c.temperature>=10,
                    "accepted finite temperature");
                require(!snrt_chimes_budget(c.next,elements,&q),"accepted finite nonnegative species");
                require(std::abs(q)<1e-10,"accepted charge");
                for(int e=0;e<11;e++)require(
                    std::abs(elements[e]-input[i].elements[e])<
                    1e-8*std::max(input[i].elements[e],1e-20),"accepted element budget");
            }
            char label[128];
            std::snprintf(label,sizeof(label),"fixed-tol%s-gpu%d-worker%d-repeat%d",factor,backend,owner,repeat);
            if(owner==0 && repeat==0){
                if(factor_index)compare(previous[backend],trial,label,false);
                std::memcpy(baseline[backend],trial,sizeof(trial));
                std::memcpy(previous[backend],trial,sizeof(trial));
            }else compare(baseline[backend],trial,label,false);
            std::printf("FIXED_RESULT %s accepted=%d rejected=%d gpu_calls=%llu cpu_fallback=%llu cvode_calls=%lu\n",
                label,accepted,nc-accepted,after[0]-before[0],after[1]-before[1],
                fixed_solver_calls-calls_before);
            // Save one full endpoint per repeated-state tile for convergence
            // analysis. A rejected cell exports no candidate physical state.
            for(int i=0;i<nc;i++)
                std::printf("FIXED_STATUS %s cell=%d status=%d\n",label,i,trial[i].status);
            if(trial[0].status==0){
                std::printf("FIXED_ENDPOINT %s T=%.17g species=",label,trial[0].temperature);
                for(int j=0;j<ns;j++)std::printf("%s%.17g",j?",":"",trial[0].next[j]);
                std::putchar('\n');
            }
        }
        char label[80];std::snprintf(label,sizeof(label),"fixed-cpu-gpu-tol%s",factor);
        compare(baseline[0],baseline[1],label,false);
        ++factor_index;
    }
    fixed_diag_active=false;
    unsetenv("SNRT_CHIMES_TEST_TOLERANCE_FACTOR");
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_FIXED_ROUTE_DIAGNOSTIC_COMPLETE"); // Not an accuracy PASS.
}

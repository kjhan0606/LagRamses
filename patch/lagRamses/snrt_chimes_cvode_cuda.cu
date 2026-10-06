/* Experimental batched dense-direct backend for SUNDIALS 5.8 / CHIMES.
 * Independent CVODE BDF states, tolerances, roots and CPU RHS are preserved.
 * Simultaneously pending dense SETUP and SOLVE calls are coalesced separately.
 * No MPI calls, background threads, global device sync, or published-state changes.
 */
#include "snrt_chimes_cvode_cuda.h"
#include "snrt_chimes_rhs_cuda.h"
#include "../cuRamses/cuda_stream_pool.h"
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_serial.h>
#include <nvector/nvector_cuda.h>
#include <sundials/sundials_config.h>
#include <cublas_v2.h>
#include <omp.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <vector>
#include <new>
#if SUNDIALS_VERSION_MAJOR != 5 || SUNDIALS_VERSION_MINOR != 8
#error This adapter requires the validated SUNDIALS 5.8 ABI
#endif
static_assert(sizeof(realtype)==sizeof(double), "CHIMES CUDA requires FP64");
namespace {
constexpr int max_batch=64, max_n=256, min_batch=4, fallback=-9999;
constexpr auto gather_time=std::chrono::microseconds(50);
bool enabled=false; // configured on serial startup, then read-only
bool solve_enabled=false;
std::atomic<unsigned long long> cpu{0}, gpu{0}, batches{0}, busy{0}, errors{0}, largest{0};
std::atomic<unsigned long long> gpu_solves{0},solve_batches{0};
std::atomic<unsigned long long> active_batches{0},peak_active_batches{0},used_stream_mask{0};
struct Request { SUNLinearSolver s; SUNMatrix a; int n; bool claimed=false,done=false; int status=fallback; };
struct SolveRequest { SUNLinearSolver s; SUNMatrix a; N_Vector x,b; int n; bool claimed=false,done=false; int status=fallback; };
std::mutex queue_mutex;
std::condition_variable changed;
std::vector<Request*> pending;
std::mutex solve_queue_mutex;
std::condition_variable solve_changed;
std::vector<SolveRequest*> solve_pending;

struct Buffers {
    int device=-1, ncap=0, bcap=0;
    double *host=nullptr,*device_a=nullptr,**device_ptr=nullptr;
    double *host_b=nullptr,*device_b=nullptr,**device_bptr=nullptr;
    int *host_ip=nullptr,*host_info=nullptr,*device_ip=nullptr,*device_info=nullptr;
    cublasHandle_t handle=nullptr;
    void clear() {
        if(device>=0)cudaSetDevice(device);
        if(handle)cublasDestroy(handle);
        cudaFree(device_a);cudaFree(device_ptr);cudaFree(device_b);cudaFree(device_bptr);
        cudaFree(device_ip);cudaFree(device_info);
        cudaFreeHost(host);cudaFreeHost(host_b);cudaFreeHost(host_ip);cudaFreeHost(host_info);
        *this=Buffers{};
    }
    bool reserve(int n,int b,cudaStream_t stream) {
        int dev=-1;
        if(cudaGetDevice(&dev)!=cudaSuccess)return false;
        if(device==dev && n<=ncap && b<=bcap && handle)
            return cublasSetStream(handle,stream)==CUBLAS_STATUS_SUCCESS;
        clear();device=dev;
        if(cudaSetDevice(dev)!=cudaSuccess)return false;
        ncap=n;bcap=b;
        size_t na=(size_t)n*n*b*sizeof(double), nb=(size_t)n*b*sizeof(double), ni=(size_t)n*b*sizeof(int);
        if(cudaMallocHost(&host,na)!=cudaSuccess || cudaMallocHost(&host_ip,ni)!=cudaSuccess ||
           cudaMallocHost(&host_b,nb)!=cudaSuccess || cudaMallocHost(&host_info,b*sizeof(int))!=cudaSuccess ||
           cudaMalloc(&device_a,na)!=cudaSuccess || cudaMalloc(&device_ip,ni)!=cudaSuccess ||
           cudaMalloc(&device_b,nb)!=cudaSuccess || cudaMalloc(&device_bptr,b*sizeof(double*))!=cudaSuccess ||
           cudaMalloc(&device_info,b*sizeof(int))!=cudaSuccess ||
           cudaMalloc(&device_ptr,b*sizeof(double*))!=cudaSuccess ||
           cublasCreate(&handle)!=CUBLAS_STATUS_SUCCESS ||
           cublasSetStream(handle,stream)!=CUBLAS_STATUS_SUCCESS){clear();return false;}
        return true;
    }
};
Buffers scratch[MAX_CUDA_STREAMS];

void factor_batch(std::vector<Request*>& rs,int slot) {
    cudaStream_t stream=cuda_get_stream_internal(slot);
    Buffers &w=scratch[slot];
    const int n=rs[0]->n, count=(int)rs.size();
    bool ok=w.reserve(n,count,stream);
    bool tracking=false;
    if(ok) {
        const auto active=active_batches.fetch_add(1)+1;
        auto peak=peak_active_batches.load();
        while(peak<active && !peak_active_batches.compare_exchange_weak(peak,active)){}
        if(slot>=0 && slot<64)used_stream_mask.fetch_or(1ull<<slot);
        tracking=true;
    }
    double* pointers[max_batch];
    const size_t one=(size_t)n*n, bytes=one*count*sizeof(double);
    if(ok) {
        for(int i=0;i<count;i++) {
            std::memcpy(w.host+i*one,SM_DATA_D(rs[i]->a),one*sizeof(double));
            pointers[i]=w.device_a+i*one;
        }
        ok=cudaMemcpyAsync(w.device_a,w.host,bytes,cudaMemcpyHostToDevice,stream)==cudaSuccess &&
           cudaMemcpyAsync(w.device_ptr,pointers,count*sizeof(double*),cudaMemcpyHostToDevice,stream)==cudaSuccess;
    }
    if(ok)ok=cublasDgetrfBatched(w.handle,n,w.device_ptr,n,w.device_ip,w.device_info,count)==CUBLAS_STATUS_SUCCESS;
    if(ok)ok=cudaMemcpyAsync(w.host,w.device_a,bytes,cudaMemcpyDeviceToHost,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.host_ip,w.device_ip,n*count*sizeof(int),cudaMemcpyDeviceToHost,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.host_info,w.device_info,count*sizeof(int),cudaMemcpyDeviceToHost,stream)==cudaSuccess;
    // Drain even on a failed enqueue: never release a lease with live DMA.
    if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
    if(tracking)active_batches.fetch_sub(1);
    if(ok) {
        batches++;
        auto high=largest.load();
        while(high<(unsigned long long)count &&
              !largest.compare_exchange_weak(high,(unsigned long long)count)){}
        for(int i=0;i<count;i++) {
            bool valid=w.host_info[i]==0;
            for(int j=0;j<n && valid;j++)valid=w.host_ip[i*n+j]>=1 && w.host_ip[i*n+j]<=n;
            // Do not publish partial LU on singular/error: CPU recomputes
            // from the untouched input and returns its native failure code.
            if(!valid)continue;
            std::memcpy(SM_DATA_D(rs[i]->a),w.host+i*one,one*sizeof(double));
            auto *content=(SUNLinearSolverContent_Dense)rs[i]->s->content;
            for(int j=0;j<n;j++)content->pivots[j]=w.host_ip[i*n+j]-1;
            content->last_flag=SUNLS_SUCCESS;rs[i]->status=SUNLS_SUCCESS;gpu++;
        }
    } else {errors+=rs.size();w.clear();}
}

void solve_batch(std::vector<SolveRequest*>& rs,int slot) {
    cudaStream_t stream=cuda_get_stream_internal(slot);
    Buffers &w=scratch[slot];
    const int n=rs[0]->n,count=(int)rs.size();
    bool ok=w.reserve(n,count,stream),tracking=false;
    if(ok) {
        const auto active=active_batches.fetch_add(1)+1;
        auto peak=peak_active_batches.load();
        while(peak<active && !peak_active_batches.compare_exchange_weak(peak,active)){}
        if(slot>=0 && slot<64)used_stream_mask.fetch_or(1ull<<slot);
        tracking=true;
    }
    double *a_ptrs[max_batch],*b_ptrs[max_batch],*x_ptrs[max_batch];
    const size_t one=(size_t)n*n,abytes=one*count*sizeof(double),
                 bbytes=(size_t)n*count*sizeof(double),ipbytes=(size_t)n*count*sizeof(int);
    if(ok)for(int i=0;i<count && ok;i++) {
        auto *content=(SUNLinearSolverContent_Dense)rs[i]->s->content;
        double *rhs=N_VGetArrayPointer(rs[i]->b),*solution=N_VGetArrayPointer(rs[i]->x);
        ok=content && content->pivots && content->last_flag==SUNLS_SUCCESS &&
           SM_DATA_D(rs[i]->a) && rhs && solution;
        if(!ok)break;
        for(int j=0;j<n;j++) {
            const sunindextype pivot=content->pivots[j];
            if(pivot<0 || pivot>=n){ok=false;break;}
            w.host_ip[(size_t)i*n+j]=(int)pivot+1; // cuBLAS pivots are 1-based
        }
        if(!ok)break;
        std::memcpy(w.host+(size_t)i*one,SM_DATA_D(rs[i]->a),one*sizeof(double));
        std::memcpy(w.host_b+(size_t)i*n,rhs,(size_t)n*sizeof(double));
        a_ptrs[i]=w.device_a+(size_t)i*one;
        b_ptrs[i]=w.device_b+(size_t)i*n;
        x_ptrs[i]=solution;
    }
    int info=-1;
    if(ok)ok=cudaMemcpyAsync(w.device_a,w.host,abytes,cudaMemcpyHostToDevice,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.device_ptr,a_ptrs,count*sizeof(double*),cudaMemcpyHostToDevice,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.device_ip,w.host_ip,ipbytes,cudaMemcpyHostToDevice,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.device_b,w.host_b,bbytes,cudaMemcpyHostToDevice,stream)==cudaSuccess &&
        cudaMemcpyAsync(w.device_bptr,b_ptrs,count*sizeof(double*),cudaMemcpyHostToDevice,stream)==cudaSuccess;
    if(ok)ok=cublasDgetrsBatched(w.handle,CUBLAS_OP_N,n,1,
        reinterpret_cast<const double *const *>(w.device_ptr),n,w.device_ip,
        w.device_bptr,n,&info,count)==CUBLAS_STATUS_SUCCESS;
    if(ok)ok=cudaMemcpyAsync(w.host_b,w.device_b,bbytes,cudaMemcpyDeviceToHost,stream)==cudaSuccess;
    // Drain every submitted operation before reusing this stream's staging buffers.
    if(cudaStreamSynchronize(stream)!=cudaSuccess || info!=0)ok=false;
    if(tracking)active_batches.fetch_sub(1);
    if(ok) {
        for(int i=0;i<count;i++) {
            std::memcpy(x_ptrs[i],w.host_b+(size_t)i*n,(size_t)n*sizeof(double));
            auto *content=(SUNLinearSolverContent_Dense)rs[i]->s->content;
            content->last_flag=SUNLS_SUCCESS;
            rs[i]->status=SUNLS_SUCCESS;gpu_solves++;
        }
        solve_batches++;
    } else {errors+=rs.size();w.clear();}
}

bool serve_setup(Request &r,int team) {
    const int nstreams=std::max(1,cuda_get_n_streams());
    const int target=std::max(min_batch,std::min(max_batch,(team+nstreams-1)/nstreams));
    const auto deadline=std::chrono::steady_clock::now()+gather_time;
    std::unique_lock<std::mutex> lock(queue_mutex);
    pending.push_back(&r);changed.notify_all();
    auto compatible_count=[&]() {
        return std::count_if(pending.begin(),pending.end(),[&](Request *q){
            return !q->claimed && !q->done && q->n==r.n;
        });
    };
    for(;;) {
        if(r.done)return r.status!=fallback;
        if(r.claimed){changed.wait(lock,[&]{return r.done;});continue;}
        if(compatible_count()<target && std::chrono::steady_clock::now()<deadline) {
            changed.wait_until(lock,deadline,[&]{
                return r.done||r.claimed||compatible_count()>=target;
            });
            continue;
        }

        // Claim a homogeneous matrix-size batch before lease acquisition.
        // This prevents every waiter from racing the pool and abandoning work
        // before one worker has reserved the requests it can serve.
        std::vector<Request*> group;group.reserve((size_t)target);
        auto own=std::find(pending.begin(),pending.end(),&r);
        if(own!=pending.end())pending.erase(own);
        r.claimed=true;group.push_back(&r);
        while(group.size()<(size_t)target) {
            auto it=std::find_if(pending.begin(),pending.end(),[&](Request *q){
                return !q->claimed && !q->done && q->n==r.n;
            });
            if(it==pending.end())break;
            Request *q=*it;pending.erase(it);q->claimed=true;group.push_back(q);
        }
        if(group.size()<(size_t)min_batch) {
            for(auto *q:group)q->done=true; // insufficient coalescing: native CPU LU
            changed.notify_all();
            return false;
        }

        lock.unlock();
        const int slot=cuda_acquire_stream();
        if(slot<0) {
            busy+=group.size();
            lock.lock();for(auto *q:group)q->done=true;changed.notify_all();
            return false;
        }
        factor_batch(group,slot);
        cuda_release_stream(slot);
        lock.lock();for(auto *q:group)q->done=true;changed.notify_all();
    }
}

bool serve_solve(SolveRequest &r,int team) {
    const int nstreams=std::max(1,cuda_get_n_streams());
    const int target=std::max(min_batch,std::min(max_batch,(team+nstreams-1)/nstreams));
    const auto deadline=std::chrono::steady_clock::now()+gather_time;
    std::unique_lock<std::mutex> lock(solve_queue_mutex);
    solve_pending.push_back(&r);solve_changed.notify_all();
    auto compatible_count=[&]() {
        return std::count_if(solve_pending.begin(),solve_pending.end(),[&](SolveRequest *q){
            return !q->claimed && !q->done && q->n==r.n;
        });
    };
    for(;;) {
        if(r.done)return r.status!=fallback;
        if(r.claimed){solve_changed.wait(lock,[&]{return r.done;});continue;}
        if(compatible_count()<target && std::chrono::steady_clock::now()<deadline) {
            solve_changed.wait_until(lock,deadline,[&]{
                return r.done||r.claimed||compatible_count()>=target;
            });
            continue;
        }
        std::vector<SolveRequest*> group;group.reserve((size_t)target);
        auto own=std::find(solve_pending.begin(),solve_pending.end(),&r);
        if(own!=solve_pending.end())solve_pending.erase(own);
        r.claimed=true;group.push_back(&r);
        while(group.size()<(size_t)target) {
            auto it=std::find_if(solve_pending.begin(),solve_pending.end(),[&](SolveRequest *q){
                return !q->claimed && !q->done && q->n==r.n;
            });
            if(it==solve_pending.end())break;
            SolveRequest *q=*it;solve_pending.erase(it);q->claimed=true;group.push_back(q);
        }
        if(group.size()<(size_t)min_batch) {
            for(auto *q:group)q->done=true;
            solve_changed.notify_all();
            return false;
        }
        lock.unlock();
        const int slot=cuda_acquire_stream();
        if(slot<0) {
            busy+=group.size();
            lock.lock();for(auto *q:group)q->done=true;solve_changed.notify_all();
            return false;
        }
        solve_batch(group,slot);
        cuda_release_stream(slot);
        lock.lock();for(auto *q:group)q->done=true;solve_changed.notify_all();
    }
}

int setup(SUNLinearSolver s,SUNMatrix a) {
    if(!s || !a)return SUNLS_MEM_NULL;
    const int n=(int)SM_ROWS_D(a), team=omp_in_parallel()?omp_get_num_threads():1;
    if(!enabled || snrt_chimes_rhs_tile_active() || snrt_chimes_rhs_cpu_direct_active() ||
       team<min_batch || n<1 || n>max_n || SM_COLUMNS_D(a)!=n) {
        cpu++;return SUNLinSolSetup_Dense(s,a);
    }
    Request r{s,a,n};
    serve_setup(r,team);
    if(r.status==fallback){cpu++;return SUNLinSolSetup_Dense(s,a);}
    return r.status;
}
int solve_linear(SUNLinearSolver s,SUNMatrix a,N_Vector x,N_Vector b,realtype tol) {
    if(!s || !a || !x || !b)return SUNLS_MEM_NULL;
    const int n=(int)SM_ROWS_D(a),team=omp_in_parallel()?omp_get_num_threads():1;
    if(!enabled || !solve_enabled || snrt_chimes_rhs_tile_active() ||
       snrt_chimes_rhs_cpu_direct_active() || team<min_batch ||
       n<1 || n>max_n || SM_COLUMNS_D(a)!=n) {
        cpu++;return SUNLinSolSolve_Dense(s,a,x,b,tol);
    }
    SolveRequest r{s,a,x,b,n};
    serve_solve(r,team);
    if(r.status==fallback){cpu++;return SUNLinSolSolve_Dense(s,a,x,b,tol);}
    return r.status;
}
SUNLinearSolver_ID custom_id(SUNLinearSolver) {return SUNLINEARSOLVER_CUSTOM;}
}
extern "C" int snrt_chimes_cvode_configure_mode(int use_cuda) {
    if(use_cuda!=0 && use_cuda!=1)return 1;
    enabled=use_cuda!=0;
    solve_enabled=false;
    if(enabled)
        std::fprintf(stderr,"CHIMES CVODE experimental CUDA batched LU; CPU RHS/control; min_batch=%d wait_us=50; CPU fallback on unavailable lease\n",min_batch);
    return 0;
}
extern "C" int snrt_chimes_cvode_is_cuda_vector(N_Vector y) {
    return y && N_VGetVectorID(y)==SUNDIALS_NVEC_CUDA;
}
extern "C" int snrt_chimes_cvode_configure(void) {
    const char *mode=std::getenv("SNRT_CHIMES_CVODE_BACKEND");
    if(!mode || !*mode || std::strcmp(mode,"cpu")==0)return snrt_chimes_cvode_configure_mode(0);
    if(std::strcmp(mode,"cuda_batched_lu")!=0 && std::strcmp(mode,"cuda_batched_lu_solve")!=0) {
        std::fprintf(stderr,"Unknown SNRT_CHIMES_CVODE_BACKEND=%s\n",mode);return 1;
    }
    const int status=snrt_chimes_cvode_configure_mode(1);
    if(!status && std::strcmp(mode,"cuda_batched_lu_solve")==0) {
        solve_enabled=true;
        std::fprintf(stderr,"CHIMES CVODE experimental CUDA batched dense solve enabled\n");
    }
    return status;
}
extern "C" SUNLinearSolver snrt_chimes_cvode_linear_solver(N_Vector y,SUNMatrix a) {
    if(N_VGetVectorID(y)==SUNDIALS_NVEC_CUDA)
        return SUNLinSol_SPGMR(y,PREC_NONE,0);
    SUNLinearSolver s=SUNLinSol_Dense(y,a);
    if(s && enabled && !snrt_chimes_rhs_cpu_direct_active() &&
       !snrt_chimes_rhs_tile_active()){
        s->ops->setup=setup;s->ops->getid=custom_id;
        if(solve_enabled)s->ops->solve=solve_linear;
    }
    return s;
}
extern "C" void snrt_chimes_cvode_counts(unsigned long long out[10]) {
    out[0]=cpu.load();out[1]=gpu.load();out[2]=batches.load();out[3]=busy.load();
    out[4]=errors.load();out[5]=largest.load();out[6]=peak_active_batches.load();
    out[7]=(unsigned long long)__builtin_popcountll(used_stream_mask.load());
    out[8]=gpu_solves.load();out[9]=solve_batches.load();
}
extern "C" void snrt_chimes_cvode_cuda_finalize(void) {
    // Caller has joined all OMP workers; call before cuda_pool_finalize.
    for(auto &w:scratch)w.clear();
}

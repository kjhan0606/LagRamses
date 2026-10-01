/* Experimental batched LU backend for SUNDIALS 5.8 / CHIMES.
 * Independent CVODE BDF states, tolerances, roots and CPU RHS are preserved.
 * Only simultaneously pending linear SETUPs are coalesced. No MPI calls,
 * background threads, global device sync, or changes to published cell state.
 */
#include "snrt_chimes_cvode_cuda.h"
#include "../cuRamses/cuda_stream_pool.h"
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
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
std::atomic<unsigned long long> cpu{0}, gpu{0}, batches{0}, busy{0}, errors{0}, largest{0};
struct Request { SUNLinearSolver s; SUNMatrix a; int n; bool done=false; int status=fallback; };
std::mutex queue_mutex;
std::condition_variable changed;
std::vector<Request*> pending;
bool processing=false;

struct Buffers {
    int device=-1, ncap=0, bcap=0;
    double *host=nullptr,*device_a=nullptr,**device_ptr=nullptr;
    int *host_ip=nullptr,*host_info=nullptr,*device_ip=nullptr,*device_info=nullptr;
    cublasHandle_t handle=nullptr;
    void clear() {
        if(device>=0)cudaSetDevice(device);
        if(handle)cublasDestroy(handle);
        cudaFree(device_a);cudaFree(device_ptr);cudaFree(device_ip);cudaFree(device_info);
        cudaFreeHost(host);cudaFreeHost(host_ip);cudaFreeHost(host_info);
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
        size_t na=(size_t)n*n*b*sizeof(double), ni=(size_t)n*b*sizeof(int);
        if(cudaMallocHost(&host,na)!=cudaSuccess || cudaMallocHost(&host_ip,ni)!=cudaSuccess ||
           cudaMallocHost(&host_info,b*sizeof(int))!=cudaSuccess ||
           cudaMalloc(&device_a,na)!=cudaSuccess || cudaMalloc(&device_ip,ni)!=cudaSuccess ||
           cudaMalloc(&device_info,b*sizeof(int))!=cudaSuccess ||
           cudaMalloc(&device_ptr,b*sizeof(double*))!=cudaSuccess ||
           cublasCreate(&handle)!=CUBLAS_STATUS_SUCCESS ||
           cublasSetStream(handle,stream)!=CUBLAS_STATUS_SUCCESS){clear();return false;}
        return true;
    }
};
Buffers scratch[MAX_CUDA_STREAMS];

void factor_batch(std::vector<Request*>& rs) {
    if(rs.size()<min_batch)return; // small batches: original CPU factorization
    const int slot=cuda_acquire_stream();
    if(slot<0){busy+=rs.size();return;}
    cudaStream_t stream=cuda_get_stream_internal(slot);
    Buffers &w=scratch[slot];
    const int n=rs[0]->n, count=(int)rs.size();
    bool ok=w.reserve(n,count,stream);
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
    if(ok) {
        batches++;
        largest.store(std::max(largest.load(),(unsigned long long)count));
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
    cuda_release_stream(slot);
}

int setup(SUNLinearSolver s,SUNMatrix a) {
    if(!s || !a)return SUNLS_MEM_NULL;
    const int n=(int)SM_ROWS_D(a), team=omp_in_parallel()?omp_get_num_threads():1;
    if(!enabled || team<min_batch || n<1 || n>max_n || SM_COLUMNS_D(a)!=n) {
        cpu++;return SUNLinSolSetup_Dense(s,a);
    }
    Request r{s,a,n};
    std::unique_lock<std::mutex> lock(queue_mutex);
    pending.push_back(&r);changed.notify_all();
    while(!r.done) {
        if(processing){changed.wait(lock,[&]{return r.done || !processing;});continue;}
        processing=true;
        const int target=std::min(team,max_batch);
        changed.wait_until(lock,std::chrono::steady_clock::now()+gather_time,
                           [&]{return (int)pending.size()>=target;});
        std::vector<Request*> group;group.swap(pending);
        lock.unlock();
        // A network switch changes matrix size. Never mix different sizes.
        std::sort(group.begin(),group.end(),[](Request *x,Request *y){return x->n<y->n;});
        for(size_t first=0;first<group.size();) {
            size_t last=first+1;
            while(last<group.size() && group[last]->n==group[first]->n && last-first<max_batch)last++;
            std::vector<Request*> compatible(group.begin()+first,group.begin()+last);
            factor_batch(compatible);first=last;
        }
        lock.lock();
        for(auto q:group)q->done=true;
        processing=false;changed.notify_all();
    }
    lock.unlock();
    if(r.status==fallback){cpu++;return SUNLinSolSetup_Dense(s,a);}
    return r.status;
}
SUNLinearSolver_ID custom_id(SUNLinearSolver) {return SUNLINEARSOLVER_CUSTOM;}
}
extern "C" int snrt_chimes_cvode_configure(void) {
    const char *mode=std::getenv("SNRT_CHIMES_CVODE_BACKEND");
    if(!mode || !*mode || std::strcmp(mode,"cpu")==0){enabled=false;return 0;}
    if(std::strcmp(mode,"cuda_batched_lu")!=0) {
        std::fprintf(stderr,"Unknown SNRT_CHIMES_CVODE_BACKEND=%s\n",mode);return 1;
    }
    enabled=true;
    std::fprintf(stderr,"CHIMES CVODE experimental CUDA batched LU; CPU RHS/triangular solve; min_batch=%d wait_us=50; CPU fallback on unavailable lease\n",min_batch);
    return 0;
}
extern "C" SUNLinearSolver snrt_chimes_cvode_linear_solver(N_Vector y,SUNMatrix a) {
    auto s=SUNLinSol_Dense(y,a);
    if(s && enabled){s->ops->setup=setup;s->ops->getid=custom_id;}
    return s;
}
extern "C" void snrt_chimes_cvode_counts(unsigned long long out[6]) {
    out[0]=cpu;out[1]=gpu;out[2]=batches;out[3]=busy;out[4]=errors;out[5]=largest;
}
extern "C" void snrt_chimes_cvode_cuda_finalize(void) {
    // Caller has joined all OMP workers; call before cuda_pool_finalize.
    for(auto &w:scratch)w.clear();
}

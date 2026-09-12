#include "../cuRamses/cuda_stream_pool.h"
#include "snrt_ir_cell.h"
#include "snrt_dust_material_cell.h"

__global__ void exchange_kernel(const double *input,const double *table,double *out,int *error,
    int n,int nt,double dt,double floor_t) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<n){int rc=dust_exchange_cell(input,table,out,nt,dt,floor_t,i);if(rc)atomicMax(error,rc);}
}
extern "C" int snrt_exchange_batch_c(const double *input,const double *table,double *output,
    int n,int nt,double dt,double floor_t,int slot) {
  cudaStream_t stream=cuda_get_stream_internal(slot);
  double *buf=nullptr;int *error=nullptr;int status=7,kernel_status=0;
  const size_t ni=5*size_t(n),nk=2*size_t(nt),no=4*size_t(n);
  if(cudaMallocAsync(&buf,(ni+nk+no)*sizeof(double),stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&error,sizeof(int),stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(buf,input,ni*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(buf+ni,table,nk*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemsetAsync(error,0,sizeof(int),stream)!=cudaSuccess)goto done;
  exchange_kernel<<<(n+127)/128,128,0,stream>>>(buf,buf+ni,buf+ni+nk,error,n,nt,dt,floor_t);
  if(cudaGetLastError()!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(&kernel_status,error,sizeof(int),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  if(kernel_status){status=kernel_status;goto done;}
  if(cudaMemcpyAsync(output,buf+ni+nk,no*sizeof(double),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  status=0;
done:
  if(error&&cudaFreeAsync(error,stream)!=cudaSuccess)status=7;
  if(buf&&cudaFreeAsync(buf,stream)!=cudaSuccess)status=7;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)status=7;
  return status;
}

__global__ void scatter_kernel(const double *input,const double *weight,double *out,int *error,
    int n,int ng,int nd,double sum_w) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<n){int rc=snrt_isotropic_scatter_cell(input,weight,out,n,ng,nd,sum_w,i);if(rc)atomicMax(error,rc);}
}
extern "C" int snrt_scatter_batch_c(const double *input,const double *weight,double *output,
    int n,int ng,int nd,double sum_w,int slot) {
  const size_t no=size_t(n)*ng*nd,ni=no+size_t(n)*ng;
  cudaStream_t stream=cuda_get_stream_internal(slot);
  double *buf=nullptr;int *error=nullptr;int status=7,kernel_status=0;
  if(cudaMallocAsync(&buf,(ni+nd+no)*sizeof(double),stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&error,sizeof(int),stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(buf,input,ni*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(buf+ni,weight,nd*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemsetAsync(error,0,sizeof(int),stream)!=cudaSuccess)goto done;
  scatter_kernel<<<(n+127)/128,128,0,stream>>>(buf,buf+ni,buf+ni+nd,error,n,ng,nd,sum_w);
  if(cudaGetLastError()!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(&kernel_status,error,sizeof(int),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  if(kernel_status){status=kernel_status;goto done;}
  if(cudaMemcpyAsync(output,buf+ni+nd,no*sizeof(double),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  status=0;
done:
  if(error && cudaFreeAsync(error,stream)!=cudaSuccess)status=7;
  if(buf && cudaFreeAsync(buf,stream)!=cudaSuccess)status=7;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)status=7;
  return status;
}

__global__ void ir_transport_kernel(const double *q,const double *rho,const int *blocked,
    const double *dir,const double *sigma,double *out,int *error,int n,int ng,int nd,double cdt,double ratio,int cell_sigma) {
  // A block owns one cell: adjacent lanes read/write adjacent frequencies.
  // Keep the original axis order and FP64 expressions; no reduced precision.
  const int i=blockIdx.x;
  const size_t rays=size_t(ng)*nd,total=rays*n,groups=size_t(ng)*n;
  for(size_t ray=threadIdx.x;ray<rays;ray+=blockDim.x) {
    const int d=ray/ng;
    const double old=q[size_t(i)*7*rays+ray];double value=old;
    for(int axis=0;axis<3;++axis) {
      int face=2*axis,outgoing=face+1;
      if(dir[3*d+axis]<0){face=2*axis+1;outgoing=face-1;}
      const double factor=ratio*fabs(dir[3*d+axis]);
      if(!blocked[6*i+outgoing])value-=factor*old;
      value+=factor*q[(size_t(i)*7+face+1)*rays+ray];
    }
    if(!isfinite(value))atomicMax(error,2);
    out[size_t(i)*rays+ray]=value;
  }
  for(int g=threadIdx.x;g<ng;g+=blockDim.x) {
    const double tau=cell_sigma?cdt*rho[size_t(i)*ng+g]:cdt*sigma[g]*rho[i];
    if(!isfinite(tau)||tau<0){atomicMax(error,2);continue;}
    const double transmit=exp(-tau);
    const double loss=tau<1e-4?tau*(1-tau/2+tau*tau/6-tau*tau*tau/24):1-transmit;
    const double response=tau<1e-4?1-tau/2+tau*tau/6:loss/fmax(tau,DBL_MIN);
    out[total+size_t(i)*ng+g]=transmit;
    out[total+groups+size_t(i)*ng+g]=loss;
    out[total+2*groups+size_t(i)*ng+g]=response;
  }
}
__global__ void ir_absorb_kernel(const double *input,const double *weight,double *out,
    int *error,int n,int ng,int nd,double dt,double sum_w) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<n){int rc=snrt_ir_absorb_cell(input,weight,out,n,ng,nd,dt,sum_w,i);if(rc)atomicMax(error,rc);}
}
// One stream-local allocation packs all real inputs and output; caller owns
// the lease until cleanup completes. No device-wide sync or error replay.
extern "C" int snrt_ir_batch_c(const double *input,const int *blocked,const double *coeff,
    double *output,int n,int ng,int nd,double a,double b,int slot,int op) {
  const size_t rays=size_t(ng)*nd,groups=size_t(ng)*n;
  const size_t ni=op!=1?7*rays*n+(op==2?groups:n):rays*n+4*groups;
  const size_t nk=op!=1?3*size_t(nd)+ng:size_t(nd);
  const size_t no=op!=1?rays*n+3*groups:rays*n+n;
  cudaStream_t stream=cuda_get_stream_internal(slot);
  double *buf=nullptr;int *flags=nullptr,*error=nullptr;
  int status=7,kernel_status=0;
  if(cudaMallocAsync(&buf,(ni+nk+no)*sizeof(double),stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&error,sizeof(int),stream)!=cudaSuccess)goto done;
  if(op!=1) {
    if(cudaMallocAsync(&flags,6*size_t(n)*sizeof(int),stream)!=cudaSuccess)goto done;
    if(cudaMemcpyAsync(flags,blocked,6*size_t(n)*sizeof(int),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  }
  if(cudaMemcpyAsync(buf,input,ni*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(buf+ni,coeff,nk*sizeof(double),cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemsetAsync(error,0,sizeof(int),stream)!=cudaSuccess)goto done;
  if(op!=1)ir_transport_kernel<<<n,128,0,stream>>>(buf,buf+7*rays*n,flags,
      buf+ni,buf+ni+3*nd,buf+ni+nk,error,n,ng,nd,a,b,op==2);
  else ir_absorb_kernel<<<(n+127)/128,128,0,stream>>>(buf,buf+ni,buf+ni+nk,error,n,ng,nd,a,b);
  if(cudaGetLastError()!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(&kernel_status,error,sizeof(int),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  if(kernel_status){status=kernel_status;goto done;}
  if(cudaMemcpyAsync(output,buf+ni+nk,no*sizeof(double),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  status=0;
done:
  if(error && cudaFreeAsync(error,stream)!=cudaSuccess)status=7;
  if(flags && cudaFreeAsync(flags,stream)!=cudaSuccess)status=7;
  if(buf && cudaFreeAsync(buf,stream)!=cudaSuccess)status=7;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)status=7;
  return status;
}

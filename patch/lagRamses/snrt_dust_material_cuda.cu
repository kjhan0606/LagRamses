#include <cuda_runtime.h>
#include "snrt_dust_material_cell.h"
#include "../cuRamses/cuda_stream_pool.h"

__global__ void dust_material_kernel(const double *input,const double *table,double *output,int *error,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<nc) {
    const int status=dust_material_cell(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,i);
    if(status)atomicMax(error,status);
  }
}

static int dust_material_stream(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,cudaStream_t stream) {
  if(nc<1||ng<1||nt<2||!isfinite(dt)||dt<=0)return 7;
  double *in_d=nullptr,*table_d=nullptr,*out_d=nullptr;int *error_d=nullptr;
  int status=7,kernel_status=0;
  const size_t in_bytes=sizeof(double)*(use_u==2?7:4)*size_t(nc),table_bytes=sizeof(double)*(ng+3)*size_t(nt);
  const size_t out_bytes=sizeof(double)*(ng+2+(use_u==2))*size_t(nc);
  if(cudaMallocAsync(&in_d,in_bytes,stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&table_d,table_bytes,stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&out_d,out_bytes,stream)!=cudaSuccess)goto done;
  if(cudaMallocAsync(&error_d,sizeof(int),stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(in_d,input,in_bytes,cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(table_d,table,table_bytes,cudaMemcpyHostToDevice,stream)!=cudaSuccess)goto done;
  if(cudaMemsetAsync(error_d,0,sizeof(int),stream)!=cudaSuccess)goto done;
  dust_material_kernel<<<(nc+127)/128,128,0,stream>>>(in_d,table_d,out_d,error_d,nc,ng,nt,use_u,dt,background,bath,tolerance);
  if(cudaGetLastError()!=cudaSuccess||cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  if(cudaMemcpyAsync(&kernel_status,error_d,sizeof(int),cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  if(kernel_status){status=kernel_status;goto done;}
  // Publish only a successful trial; a copy failure is left to the enclosing
  // native radiation/material transaction, never replayed on CPU here.
  if(cudaMemcpyAsync(output,out_d,out_bytes,cudaMemcpyDeviceToHost,stream)!=cudaSuccess)goto done;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)goto done;
  status=0;
done:
  if(error_d && cudaFreeAsync(error_d,stream)!=cudaSuccess)status=7;
  if(out_d && cudaFreeAsync(out_d,stream)!=cudaSuccess)status=7;
  if(table_d && cudaFreeAsync(table_d,stream)!=cudaSuccess)status=7;
  if(in_d && cudaFreeAsync(in_d,stream)!=cudaSuccess)status=7;
  if(cudaStreamSynchronize(stream)!=cudaSuccess)status=7;
  return status;
}
extern "C" int snrt_dust_material_cuda_c(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int) {
  return dust_material_stream(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,nullptr);
}
extern "C" int snrt_dust_material_batch_c(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int slot) {
  return dust_material_stream(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,cuda_get_stream_internal(slot));
}
extern "C" int snrt_hybrid_try_acquire_c(long long bytes,int sharers) {
  const int slot=cuda_acquire_stream();
  if(slot<0)return -1;
  size_t free=0,total=0;
  if(cudaMemGetInfo(&free,&total)!=cudaSuccess || bytes>0.8*double(free)/sharers) {
    cuda_release_stream(slot);
    return -1; // before allocation or kernel launch only
  }
  return slot;
}

#include <cuda_runtime.h>
#include "snrt_dust_material_cell.h"

__global__ void dust_material_kernel(const double *input,const double *table,double *output,int *error,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance) {
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<nc) {
    const int status=dust_material_cell(input,table,output,nc,ng,nt,use_u,dt,background,bath,tolerance,i);
    if(status)atomicMax(error,status);
  }
}

extern "C" int snrt_dust_material_cuda_c(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int) {
  if(nc<1||ng<1||nt<2||!isfinite(dt)||dt<=0)return 7;
  double *in_d=nullptr,*table_d=nullptr,*out_d=nullptr;int *error_d=nullptr;
  int status=7,kernel_status=0;
  const size_t in_bytes=sizeof(double)*4*size_t(nc),table_bytes=sizeof(double)*(ng+3)*size_t(nt);
  const size_t out_bytes=sizeof(double)*(ng+2)*size_t(nc);
  if(cudaMalloc(&in_d,in_bytes)!=cudaSuccess)goto done;
  if(cudaMalloc(&table_d,table_bytes)!=cudaSuccess)goto done;
  if(cudaMalloc(&out_d,out_bytes)!=cudaSuccess)goto done;
  if(cudaMalloc(&error_d,sizeof(int))!=cudaSuccess)goto done;
  if(cudaMemcpy(in_d,input,in_bytes,cudaMemcpyHostToDevice)!=cudaSuccess)goto done;
  if(cudaMemcpy(table_d,table,table_bytes,cudaMemcpyHostToDevice)!=cudaSuccess)goto done;
  if(cudaMemset(error_d,0,sizeof(int))!=cudaSuccess)goto done;
  dust_material_kernel<<<(nc+127)/128,128>>>(in_d,table_d,out_d,error_d,nc,ng,nt,use_u,dt,background,bath,tolerance);
  if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)goto done;
  if(cudaMemcpy(&kernel_status,error_d,sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess)goto done;
  if(kernel_status){status=kernel_status;goto done;}
  // Publish only a successful trial; a copy failure is left to the enclosing
  // native radiation/material transaction, never replayed on CPU here.
  if(cudaMemcpy(output,out_d,out_bytes,cudaMemcpyDeviceToHost)!=cudaSuccess)goto done;
  status=0;
done:
  cudaFree(error_d);cudaFree(out_d);cudaFree(table_d);cudaFree(in_d);
  return status;
}

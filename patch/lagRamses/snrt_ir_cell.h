#pragma once
#include <cmath>
#include <cfloat>
#include <cstddef>
#ifdef __CUDACC__
#define IR_HD __host__ __device__
#else
#define IR_HD
#endif
// FP64 scalar kernels shared by CPU and GPU. Fortran retains independent
// reference loops and all domain/interface/energy-balance reductions.
IR_HD inline int snrt_ir_transport_cell(const double *q,const double *rho,const int *blocked,
    const double *direction,const double *sigma,double *out,int n,int ng,int nd,
    double cdt,double ratio,int i) {
  const size_t rays=size_t(ng)*nd,total=rays*n,groups=size_t(ng)*n;
  for(int d=0;d<nd;++d)for(int g=0;g<ng;++g) {
    const size_t ray=size_t(d)*ng+g;
    const double old=q[size_t(i)*7*rays+ray];
    double value=old;
    for(int axis=0;axis<3;++axis) {
      int face=2*axis,outgoing=face+1;
      if(direction[3*d+axis]<0){face=2*axis+1;outgoing=face-1;}
      const double factor=ratio*fabs(direction[3*d+axis]);
      if(!blocked[6*i+outgoing])value-=factor*old;
      value+=factor*q[(size_t(i)*7+face+1)*rays+ray];
    }
    if(!isfinite(value))return 2;
    out[size_t(i)*rays+ray]=value;
  }
  for(int g=0;g<ng;++g) {
    const double tau=cdt*sigma[g]*rho[i];
    if(!isfinite(tau)||tau<0)return 2;
    const double transmit=exp(-tau);
    double loss,response;
    if(tau<1e-4) {
      loss=tau*(1-tau/2+tau*tau/6-tau*tau*tau/24);
      response=1-tau/2+tau*tau/6;
    } else {
      loss=1-transmit;response=loss/fmax(tau,DBL_MIN);
    }
    out[total+size_t(i)*ng+g]=transmit;
    out[total+groups+size_t(i)*ng+g]=loss;
    out[total+2*groups+size_t(i)*ng+g]=response;
  }
  return 0;
}
IR_HD inline int snrt_ir_absorb_cell(const double *input,const double *weight,double *out,
    int n,int ng,int nd,double dt,double sum_w,int i) {
  const size_t rays=size_t(ng)*nd,total=rays*n,groups=size_t(ng)*n;
  double absorbed=0;
  for(int d=0;d<nd;++d)for(int g=0;g<ng;++g) {
    const size_t k=size_t(i)*ng+g,ray=size_t(i)*rays+size_t(d)*ng+g;
    const double source=dt*input[total+3*groups+k]/sum_w;
    const double response=input[total+2*groups+k];
    const double candidate=input[ray]*input[total+k]+source*response;
    absorbed+=weight[d]*(input[ray]*input[total+groups+k]+source*(1-response));
    if(!isfinite(candidate)||candidate<0)return 2;
    out[ray]=candidate;
  }
  if(!isfinite(absorbed)||absorbed<0)return 2;
  out[total+i]=absorbed;
  return 0;
}
// Exact local isotropic elastic-scattering solution for direction-integrated
// photon bins. No absorption/heating; total photons in EACH group is invariant.
// First-order splitting from spatial transport is owned by the live caller.
IR_HD inline int snrt_isotropic_scatter_cell(const double *input,const double *weight,double *out,
    int n,int ng,int nd,double sum_w,int i) {
  const size_t rays=size_t(ng)*nd,total=rays*n;
  for(int g=0;g<ng;++g) {
    const double tau=input[total+size_t(i)*ng+g];
    if(!isfinite(tau)||tau<0)return 2;
    double photons=0;
    for(int d=0;d<nd;++d) {
      const double q=input[size_t(i)*rays+size_t(g)*nd+d];
      if(!isfinite(q)||q<0)return 2;
      photons+=q;
    }
    if(!isfinite(photons))return 2;
    const double mixed=-expm1(-tau);
    for(int d=0;d<nd;++d) {
      const size_t k=size_t(i)*rays+size_t(g)*nd+d;
      const double q=input[k]*(1-mixed)+photons*(weight[d]/sum_w)*mixed;
      if(!isfinite(q)||q<0)return 2;
      out[k]=q;
    }
  }
  return 0;
}
#undef IR_HD

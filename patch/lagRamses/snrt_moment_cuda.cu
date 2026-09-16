// Fixed-grid benchmark backend: FP64 CUDA S_N and M_N, optional FP64 WMMA
// covariance products. Kinetic MUSCL/SSPRK2 with an explicit moment-projection
// receipt: finite closure residuals are not carried as an unrealizable state.
// Mixed precision is explicit. No coupled RAMSES claim, isotropic mixing or CPU fallback.
// Build with --fmad=false: this is the established compensated-sum baseline.
// Do not silently change its contraction policy in a transport comparison.
#include <cuda_runtime.h>
#include <mma.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <new>
#include <cfloat>
#include <algorithm>
#include "snrt_cweno3.h"

namespace {
constexpr int Q=8192, INLINE_Q=384, D=35, P=48, K=32, TILE=16, HALO=20, THREADS=128;
struct Basis { int nm,nq,limiter; const double *y,*w,*dir; };
template<class Real,bool Global> struct AngularWorkspace;
template<class Real> struct AngularWorkspace<Real,false> {
  Real eta[INLINE_Q],nexteta[INLINE_Q],inc[INLINE_Q],p[INLINE_Q],np[INLINE_Q];
};
template<class Real> struct AngularWorkspace<Real,true> {
  Real *eta,*nexteta,*inc,*p,*np;
};
template<class Real,bool Global> struct Workspace : AngularWorkspace<Real,Global> {
  Real target[D], mean[D], nextmean[D], g[D], step[D], v[D];
  __align__(32) Real h[P*P];
  Real l[D*D];
  __align__(32) Real a[P*K];
  Real objective,nextobjective,err,radius,dec,scale,merit,slope,change;
  double reduce[4];
  double norm; // never cast tiny physical photon densities to FP32
  double alpha[D]; // compact warm-start coefficients; never physical state
  int status,accept,warm;
};
template<class Real> __device__ Real dotrow(const Real *x,const Real *y,int n) {
  Real s=0; for(int j=0;j<n;++j)s+=x[j]*y[j]; return s;
}
// Warp/block reductions avoid repeating the same 384-term sum in 128 lanes.
__device__ double block_reduce(double v,double *scratch,bool maximum){
  const int lane=threadIdx.x%32,warp=threadIdx.x/32;
  for(int offset=16;offset;offset/=2){
    double other=__shfl_down_sync(0xffffffff,v,offset);
    v=maximum?fmax(v,other):v+other;
  }
  if(lane==0)scratch[warp]=v;
  __syncthreads();
  if(warp==0){
    v=lane<blockDim.x/32?scratch[lane]:(maximum?-DBL_MAX:0);
    for(int offset=16;offset;offset/=2){
      double other=__shfl_down_sync(0xffffffff,v,offset);
      v=maximum?fmax(v,other):v+other;
    }
    if(lane==0)scratch[0]=v;
  }
  __syncthreads();double result=scratch[0];__syncthreads();return result;
}
template<class Real,bool Global> __device__ void dual(Basis b,Workspace<Real,Global> &s,bool trial) {
  Real *eta=trial?s.nexteta:s.eta,*p=trial?s.np:s.p,*mean=trial?s.nextmean:s.mean;
  double top=-DBL_MAX;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)top=fmax(top,double(eta[q]));
  top=block_reduce(top,s.reduce,true);
  double z=0;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
    p[q]=Real(b.w[q])*exp(eta[q]-Real(top));z+=double(p[q]);
  }
  z=block_reduce(z,s.reduce,false);
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)p[q]/=Real(z);
  if(threadIdx.x==0){if(trial)s.nextobjective=Real(top+log(z));else s.objective=Real(top+log(z));}
  __syncthreads();
  for(int j=threadIdx.x;j<b.nm-1;j+=blockDim.x){
    Real t=0,comp=0;
    for(int q=0;q<b.nq;++q){
      Real value=Real(b.y[j+1+b.nm*q])*p[q];
      if constexpr(sizeof(Real)==4){Real y=value-comp,next=t+y;comp=(next-t)-y;t=next;}
      else t+=value;
    }
    mean[j]=t;
  }
  __syncthreads();
}
template<class Real,bool Global> __device__ void covariance(Basis b,Workspace<Real,Global> &s,int tensor) {
  const int d=b.nm-1,edge=sizeof(Real)==8?8:16,nt=(d+edge-1)/edge,warp=threadIdx.x/32;
  for(int j=threadIdx.x;j<P*P;j+=blockDim.x)s.h[j]=0;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)s.inc[q]=sqrt(s.p[q]);
  __syncthreads();
  for(int start=0;start<b.nq;start+=K){
    for(int k=threadIdx.x;k<P*K;k+=blockDim.x){
      int j=k/K,q=start+k%K;
      s.a[k]=(j<d&&q<b.nq)?(Real(b.y[j+1+b.nm*q])-s.mean[j])*s.inc[q]:0;
    }
    __syncthreads();
    if(tensor){
#if __CUDA_ARCH__ >= 800
      using namespace nvcuda;
      if constexpr(sizeof(Real)==8){
      for(int t=warp;t<nt*nt;t+=blockDim.x/32){
        int i=(t/nt)*8,j=(t%nt)*8;
        wmma::fragment<wmma::matrix_a,8,8,4,Real,wmma::row_major> a;
        wmma::fragment<wmma::matrix_b,8,8,4,Real,wmma::col_major> bb;
        wmma::fragment<wmma::accumulator,8,8,4,Real> c;
        wmma::load_matrix_sync(c,s.h+i*P+j,P,wmma::mem_row_major);
        for(int k=0;k<K;k+=4){
          wmma::load_matrix_sync(a,s.a+i*K+k,K);
          wmma::load_matrix_sync(bb,s.a+j*K+k,K);
          wmma::mma_sync(c,a,bb,c);
        }
        wmma::store_matrix_sync(s.h+i*P+j,c,P,wmma::mem_row_major);
      }
      }else{
        for(int t=warp;t<nt*nt;t+=blockDim.x/32){
          int i=(t/nt)*16,j=(t%nt)*16;
          wmma::fragment<wmma::matrix_a,16,16,8,wmma::precision::tf32,wmma::row_major> a;
          wmma::fragment<wmma::matrix_b,16,16,8,wmma::precision::tf32,wmma::col_major> bb;
          wmma::fragment<wmma::matrix_a,16,16,8,wmma::precision::tf32,wmma::row_major> al;
          wmma::fragment<wmma::matrix_b,16,16,8,wmma::precision::tf32,wmma::col_major> bl;
          wmma::fragment<wmma::accumulator,16,16,8,float> c;
          wmma::load_matrix_sync(c,s.h+i*P+j,P,wmma::mem_row_major);
          for(int k=0;k<K;k+=8){
            wmma::load_matrix_sync(a,s.a+i*K+k,K);
            wmma::load_matrix_sync(bb,s.a+j*K+k,K);
            // Three TF32 products recover the leading FP32 input residuals.
            // Plain TF32 covariance lost convergence on M3-M5; do not hide
            // that by accepting its inaccurate original-moment residual.
            for(int e=0;e<a.num_elements;++e){float v=a.x[e];a.x[e]=wmma::__float_to_tf32(v);al.x[e]=wmma::__float_to_tf32(v-a.x[e]);}
            for(int e=0;e<bb.num_elements;++e){float v=bb.x[e];bb.x[e]=wmma::__float_to_tf32(v);bl.x[e]=wmma::__float_to_tf32(v-bb.x[e]);}
            wmma::mma_sync(c,a,bb,c);
            wmma::mma_sync(c,a,bl,c);
            wmma::mma_sync(c,al,bb,c);
          }
          wmma::store_matrix_sync(s.h+i*P+j,c,P,wmma::mem_row_major);
        }
      }
#endif
    }else{
      for(int t=threadIdx.x;t<d*d;t+=blockDim.x){
        int i=t/d,j=t%d;Real c=s.h[i*P+j];
        for(int k=0;k<K;++k)c+=s.a[i*K+k]*s.a[j*K+k];
        s.h[i*P+j]=c;
      }
    }
    __syncthreads();
  }
}
// Cooperative left-looking Cholesky: each column's off-diagonal entries
// are independent. Preserve each dot-product's summation order.
template<class Real,bool Global> __device__ bool newton(Workspace<Real,Global> &s,int d){
  double scale=1;
  for(int t=threadIdx.x;t<d*d;t+=blockDim.x)scale=fmax(scale,fabs(double(s.h[(t/d)*P+t%d])));
  Real matrix_scale=Real(block_reduce(scale,s.reduce,true)),shift=0;
  for(int attempt=0;attempt<10;++attempt){
    for(int t=threadIdx.x;t<d*d;t+=blockDim.x)s.l[t]=0;
    if(threadIdx.x==0)s.accept=1;
    __syncthreads();
    for(int j=0;j<d;++j){
      if(threadIdx.x==0){
        Real v=s.h[j*P+j]-dotrow(s.l+j*d,s.l+j*d,j);v+=shift;
        if(!isfinite(v)||v<=(sizeof(Real)==8?DBL_EPSILON:FLT_EPSILON)*matrix_scale)s.accept=0;
        else s.l[j*d+j]=sqrt(v);
      }
      __syncthreads();
      bool column_ok=s.accept!=0;
      __syncthreads(); // all warps consume this flag before another write
      if(!column_ok)break;
      for(int i=j+1+threadIdx.x;i<d;i+=blockDim.x)
        s.l[i*d+j]=(s.h[i*P+j]-dotrow(s.l+i*d,s.l+j*d,j))/s.l[j*d+j];
      __syncthreads();
    }
    bool factor_ok=s.accept!=0;
    __syncthreads(); // next attempt resets both the flag and factor storage
    if(factor_ok)break;
    shift=fmax(Real((sizeof(Real)==8?1e-15:1e-7)*matrix_scale),Real(10)*shift);
  }
  bool factor_ok=s.accept!=0;
  __syncthreads();
  if(!factor_ok)return false;
  // Small triangular solves retain their original arithmetic order.
  if(threadIdx.x==0){
    for(int i=0;i<d;++i)s.v[i]=(s.g[i]-dotrow(s.l+i*d,s.v,i))/s.l[i*d+i];
    for(int i=d-1;i>=0;--i){
      Real v=s.v[i];for(int j=i+1;j<d;++j)v-=s.l[j*d+i]*s.step[j];
      s.step[i]=v/s.l[i*d+i];if(!isfinite(s.step[i]))s.accept=0;
    }
  }
  __syncthreads();
  bool solved=s.accept!=0;
  __syncthreads(); // caller reuses accept: do not leave any delayed readers
  return solved;
}
template<class Real,bool Global> __device__ Real local_change(Basis b,Workspace<Real,Global> &s){
  double partial=0;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
    Real inc=s.scale*s.inc[q],term;
    if(inc>50){
      Real post=log(b.w[q])+s.eta[q]-s.objective+inc;
      term=post>50?Real(1e30):exp(post)-s.p[q];
    }else if(inc< -50)term=-s.p[q];
    else {Real v=exp(inc);term=v==1?s.p[q]*inc:(v==0?-s.p[q]:s.p[q]*(v-1)*(inc/log(v)));}
    partial+=double(term);
  }
  Real total=Real(block_reduce(partial,s.reduce,false));
  if(!isfinite(total)||total>Real(1e20))return Real(DBL_MAX);
  if(total<=-.5)return s.nextobjective-s.objective;
  Real v=1+total;return v==1?total:log(v)*(total/(v-1));
}
__device__ int cell_index(int i,int j,int k,int n){
  return (i<0||j<0||k<0||i>=n||j>=n||k>=n)?-1:i+n*(j+n*k);
}
__device__ bool record_failure(int cell,int reason,double residual,int *failure,
    int *retry_mask,int retrying,unsigned long long *stats){
  if(retry_mask&&!retrying){retry_mask[cell]=1;atomicAdd(stats,1ULL);return false;}
  if(retrying)atomicAdd(stats+1,1ULL);
  if(atomicCAS(failure,0,cell+1)==0){
    printf("GPU_CLOSURE_FAIL cell=%d reason=%d residual=%.17g retry=%d\n",cell+1,reason,residual,retrying);
    return true;
  }
  return false;
}
// One block cooperates on one independent closure. Tile halos bound angular
// storage. Final checks use contiguous cell batches and discard intensities.
template<class Real,bool Global> __global__ void reconstruct(Basis b,int n,const double *u,double *angular,
                            int ox,int oy,int oz,int nx,int ny,int nz,
                            int offset,int count,int tensor,int *failure,double tolerance,int *retry_mask,int retrying,unsigned long long *stats,
                            const double *hint,const double *hint_target,const int *hint_valid,
                            double *next_hint,double *next_target,int *next_valid,int pass,double *scratch,double *compact_cache=nullptr){
  extern __shared__ __align__(32) unsigned char shared[];
  Workspace<Real,Global> &s=*reinterpret_cast<Workspace<Real,Global>*>(shared);
  int t=blockIdx.x,cell;
  if(offset>=0){if(t>=count)return;cell=offset+t;}
  else {int i=t%nx,j=t/nx%ny,k=t/(nx*ny);cell=cell_index(ox+i-2,oy+j-2,oz+k-2,n);}
  if(cell<0){if(!retrying)for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)angular[t*b.nq+q]=0;return;}
  if(retrying&&retry_mask[cell]==0)return;
  if(retry_mask&&!retrying&&threadIdx.x==0)retry_mask[cell]=0;
  if(threadIdx.x==0){
    if constexpr(Global){
      Real *storage=reinterpret_cast<Real*>(scratch)+size_t(t)*5*b.nq;
      s.eta=storage;s.nexteta=storage+b.nq;s.inc=storage+2*b.nq;s.p=storage+3*b.nq;s.np=storage+4*b.nq;
    }
    s.norm=u[cell*b.nm];s.status=0;s.radius=1;
    for(int j=0;j<b.nm;++j)if(!isfinite(u[cell*b.nm+j]))s.status=1;
    if(s.norm<0)s.status=1;
    if(s.norm==0)for(int j=1;j<b.nm;++j)if(u[cell*b.nm+j]!=0)s.status=1;
    if(s.status&&atomicCAS(failure,0,cell+1)==0)
      printf("GPU_INPUT_FAIL cell=%d norm=%.17g moment1=%.17g\n",cell+1,s.norm,u[cell*b.nm+1]);
  }
  __syncthreads();
  if(s.status)return;
  if(s.norm==0){
    if(next_valid&&threadIdx.x==0)next_valid[cell]=0;
    if(angular)for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)angular[t*b.nq+q]=0;
    return;
  }
  if(!retrying&&threadIdx.x==0)atomicAdd(stats+2,1ULL);
  for(int j=threadIdx.x;j<b.nm-1;j+=blockDim.x)s.target[j]=u[cell*b.nm+j+1]/s.norm;
  __syncthreads();
  const int d=b.nm-1;
  if(threadIdx.x==0)for(int j=0;j<d;++j)if(!isfinite(s.target[j]))s.status=1;
  __syncthreads();
  if(s.status){if(threadIdx.x==0)atomicCAS(failure,0,cell+1);return;}
  __syncthreads(); // seed initialization below reuses status
  // One cache record per cell, updated only after full original-moment
  // verification. Kernel launches in the single stream serialize tile/stage
  // reuse; each cell occurs at most once within a launch.
  for(int seed_attempt=0;seed_attempt<2;++seed_attempt){
    if(threadIdx.x==0){
      s.status=0;s.radius=1;
      s.warm=seed_attempt==0&&!retrying&&hint_valid&&hint_valid[cell];
      if(s.warm)for(int j=0;j<d;++j){
        double a=hint[size_t(cell)*d+j],target=hint_target[size_t(cell)*d+j];
        if(!isfinite(a)||fabs(a)>1e6||!isfinite(target)||fabs(target-u[cell*b.nm+j+1]/s.norm)>.25)s.warm=0;
      }
      if(s.warm)atomicAdd(stats+4,1ULL);
    }
    __syncthreads();
    for(int j=threadIdx.x;j<d;j+=blockDim.x)s.alpha[j]=s.warm?hint[size_t(cell)*d+j]:0;
    __syncthreads();
    for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
      double eta=0;for(int j=0;j<d;++j)eta+=s.alpha[j]*(b.y[j+1+b.nm*q]-double(s.target[j]));
      s.eta[q]=Real(eta);
    }
    __syncthreads();
  // Match the native solver budget; the residual tolerance is unchanged.
  for(int it=0;it<1024;++it){
    if(threadIdx.x==0)atomicAdd(stats+7,1ULL);
    dual(b,s,false);
    if(threadIdx.x==0){
      s.err=0;s.merit=0;
      for(int j=0;j<d;++j){s.g[j]=s.mean[j]-s.target[j];if(!isfinite(s.g[j]))s.status=1;s.err=fmax(s.err,fabs(s.g[j]));s.merit+=s.g[j]*s.g[j];}
      if(!isfinite(s.err))s.status=1;
    }
    __syncthreads();
    if(s.status)break;
    if(s.err<=tolerance){
      // FP32 closure is explicitly approximate; normalize its probabilities
      // in FP64 so monopole/photon conservation is not reduced to FP32.
      double partition=0;for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)partition+=double(s.p[q]);
      partition=block_reduce(partition,s.reduce,false);
      if(threadIdx.x==0)s.accept=1;
      __syncthreads();
      for(int j=threadIdx.x;j<d;j+=blockDim.x){
        double moment=0;for(int q=0;q<b.nq;++q)moment+=b.y[j+1+b.nm*q]*(double(s.p[q])/partition);
        if(!isfinite(moment)||fabs(moment-u[cell*b.nm+j+1]/s.norm)>(1.5*tolerance))atomicExch(&s.accept,0);
      }
      for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
        double value=u[cell*b.nm]*(double(s.p[q])/partition)/b.w[q];
        if(!isfinite(value)||value<0)atomicExch(&s.accept,0);
        if(angular)angular[t*b.nq+q]=value;
      }
      __syncthreads();
      if(s.accept){
        if(compact_cache){
          // Cache the FP64 exponential itself, not a lossy angular field.
          // Normalization is shifted before exponentiation, as on the CPU.
          double top=-DBL_MAX;
          for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
            double eta=0;for(int j=0;j<d;++j)eta+=s.alpha[j]*(b.y[j+1+b.nm*q]-u[cell*b.nm+j+1]/s.norm);
            top=fmax(top,eta);
          }
          top=block_reduce(top,s.reduce,true);
          double z=0;
          for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
            double eta=0;for(int j=0;j<d;++j)eta+=s.alpha[j]*(b.y[j+1+b.nm*q]-u[cell*b.nm+j+1]/s.norm);
            z+=b.w[q]*exp(eta-top);
          }
          z=block_reduce(z,s.reduce,false);
          for(int j=threadIdx.x;j<d;j+=blockDim.x)compact_cache[t*(b.nm+1)+j]=s.alpha[j];
          if(threadIdx.x==0){compact_cache[t*(b.nm+1)+d]=top;compact_cache[t*(b.nm+1)+d+1]=z;}
          for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
            double eta=0;for(int j=0;j<d;++j)eta+=s.alpha[j]*(b.y[j+1+b.nm*q]-u[cell*b.nm+j+1]/s.norm);
            angular[t*b.nq+q]=s.norm*exp(eta-top)/z;
          }
        }
        if(hint){
          for(int j=threadIdx.x;j<d;j+=blockDim.x){
            next_hint[size_t(cell)*d+j]=s.alpha[j];
            next_target[size_t(cell)*d+j]=u[cell*b.nm+j+1]/s.norm;
          }
          __syncthreads();if(threadIdx.x==0)next_valid[cell]=1;
        }
        if(s.warm&&threadIdx.x==0)atomicAdd(stats+5,1ULL);
        return;
      }
      if(threadIdx.x==0)s.status=5;
      __syncthreads();break;
    }
    if(tensor&&it==0&&threadIdx.x==0)atomicAdd(stats+3,1ULL);
    covariance(b,s,tensor);
    bool solved=newton(s,d);
    if(!solved){if(threadIdx.x==0)s.status=2;__syncthreads();break;}
    for(int j=threadIdx.x;j<d;j+=blockDim.x)s.v[j]=dotrow(s.h+j*P,s.step,d);
    __syncthreads();
    if(threadIdx.x==0){
      s.dec=dotrow(s.g,s.step,d);
      if(!isfinite(s.dec)||s.dec<=0)s.status=2;
      s.scale=fmin(Real(1),Real(s.radius/sqrt(s.dec)));s.slope=0;
      for(int j=0;j<d;++j)s.slope+=s.g[j]*s.v[j];
    }
    __syncthreads();
    if(s.status)break;
    for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){Real v=0;for(int j=0;j<d;++j)v-=s.step[j]*(Real(b.y[j+1+b.nm*q])-s.target[j]);s.inc[q]=v;}
    __syncthreads();
    int ls;
    for(ls=0;ls<48;++ls){
      for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)s.nexteta[q]=s.eta[q]+s.scale*s.inc[q];
      __syncthreads();dual(b,s,true);
      Real change=local_change(b,s);
      if(threadIdx.x==0){
        Real err=0,merit=0;
        bool finite=true;
        for(int j=0;j<d;++j){Real g=s.nextmean[j]-s.target[j];finite=finite&&isfinite(g);err=fmax(err,fabs(g));merit+=g*g;}
        s.change=change;
        s.accept=(err<=tolerance)||(s.change<0&&s.change<=-1e-4*s.scale*s.dec);
        if(s.err<(sizeof(Real)==8?1e-7:1e-3)&&merit<s.merit&&merit<=s.merit-1e-4*s.scale*fmax(Real(0),s.slope))s.accept=1;
        if(!finite)s.accept=0;
        if(!s.accept)s.scale*=.5;
      }
      __syncthreads();if(s.accept)break;
    }
    if(ls==48){if(threadIdx.x==0)s.status=3;__syncthreads();break;}
    if(threadIdx.x==0&&s.err>=(sizeof(Real)==8?1e-7:1e-3)){
      Real curvature=0;for(int j=0;j<d;++j)curvature+=s.step[j]*s.v[j];
      Real predicted=s.scale*s.dec-.5*s.scale*s.scale*curvature;
      if(predicted>0){Real ratio=-s.change/predicted;
        if(ratio<.25)s.radius=fmax(1e-8,.25*s.radius);
        if(ratio>.75&&s.scale*sqrt(s.dec)>.8*s.radius)s.radius=fmin(Real(16),Real(2)*s.radius);
      }
    }
    for(int q=threadIdx.x;q<b.nq;q+=blockDim.x)s.eta[q]=s.nexteta[q];
    for(int j=threadIdx.x;j<d;j+=blockDim.x)s.alpha[j]-=double(s.scale)*double(s.step[j]);
    __syncthreads();
  }
    if(seed_attempt==0&&s.warm){
      if(threadIdx.x==0)atomicAdd(stats+6,1ULL);
      __syncthreads();continue;
    }
    break;
  }
  if(threadIdx.x==0&&record_failure(cell,s.status?s.status:4,double(s.err),failure,retry_mask,retrying,stats)){
    // Capture only the first failed target for a cheap isolated reproduction.
    // No accepted state, stopping tolerance or iteration path is changed.
    printf("GPU_FAIL_LOCATION pass=%d offset=%d tile=%d,%d,%d norm=%.17g radius=%.17g\n",pass,offset,ox,oy,oz,s.norm,double(s.radius));
    for(int j=0;j<b.nm;++j)printf("GPU_FAIL_U %d %.17g\n",j,u[cell*b.nm+j]);
    for(int j=0;j<d;++j)printf("GPU_FAIL_ALPHA %d %.17g\n",j,s.alpha[j]);
  }
}

__device__ double minmod(double a,double b){return (a>0&&b>0)?fmin(a,b):((a<0&&b<0)?fmax(a,b):0);}
__device__ double limited_slope(double a,double b,int limiter){
  double c=minmod(a,b);
  if(limiter==2)return minmod(.5*a+.5*b,2*c);
  if(limiter==3)return minmod(.5*a+.5*b,1.5*c);
  return c;
}
__device__ double intensity(const double *a,int q,int i,int j,int k,int nx,int ny,int nz,int nq){
  if(i<0||j<0||k<0||i>=nx||j>=ny||k>=nz)return 0;
  return a[q+nq*(i+nx*(j+ny*k))];
}
// One gradient per tile cell/angle, reused by both faces and all moments.
// The one-cell gradient halo reads the existing two-cell intensity halo.
__global__ void vertex_gradients(Basis b,int n,int mn,const double *a,double *gradient,
                                 int ox,int oy,int oz,int tx,int ty,int tz){
  const int gx=tx+2,gy=ty+2,t=blockIdx.x;
  const int i=t%gx,j=t/gx%gy,k=t/(gx*gy);
  const int x=mn?i+1:ox+i-1,y=mn?j+1:oy+j-1,z=mn?k+1:oz+k-1;
  const int nx=mn?tx+4:n,ny=mn?ty+4:n,nz=mn?tz+4:n;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
    double g[3],v=intensity(a,q,x,y,z,nx,ny,nz,b.nq),theta=1;
    if(v==0){for(int ax=0;ax<3;++ax)gradient[q+b.nq*(ax+3*t)]=0;continue;}
    g[0]=.5*(intensity(a,q,x+1,y,z,nx,ny,nz,b.nq)-intensity(a,q,x-1,y,z,nx,ny,nz,b.nq));
    g[1]=.5*(intensity(a,q,x,y+1,z,nx,ny,nz,b.nq)-intensity(a,q,x,y-1,z,nx,ny,nz,b.nq));
    g[2]=.5*(intensity(a,q,x,y,z+1,nx,ny,nz,b.nq)-intensity(a,q,x,y,z-1,nx,ny,nz,b.nq));
    if(g[0]==0&&g[1]==0&&g[2]==0){for(int ax=0;ax<3;++ax)gradient[q+b.nq*(ax+3*t)]=0;continue;}
    for(int sz=-1;sz<=1;sz+=2)for(int sy=-1;sy<=1;sy+=2)for(int sx=-1;sx<=1;sx+=2){
      double lo=v,hi=v;
      for(int iz=0;iz<=1;++iz)for(int iy=0;iy<=1;++iy)for(int ix=0;ix<=1;++ix){
        const double vj=intensity(a,q,x+ix*sx,y+iy*sy,z+iz*sz,nx,ny,nz,b.nq);
        lo=fmin(lo,vj);hi=fmax(hi,vj);
      }
      const double change=.5*(sx*g[0]+sy*g[1]+sz*g[2]);
      if(change>0)theta=fmin(theta,(hi-v)/change);
      if(change<0)theta=fmin(theta,(lo-v)/change);
    }
    const double excursion=.5*(fabs(g[0])+fabs(g[1])+fabs(g[2]));
    if(excursion>0)theta=fmin(theta,.75*v/excursion);
    for(int ax=0;ax<3;++ax)gradient[q+b.nq*(ax+3*t)]=fmax(0.,theta)*g[ax];
  }
}
__global__ void cweno_faces(Basis b,int n,int mn,const double *a,double *faces,int ox,int oy,int oz,int tx,int ty,int tz){
  const int gx=tx+2,gy=ty+2,t=blockIdx.x,i=t%gx,j=t/gx%gy,k=t/(gx*gy);
  const int x=mn?i+1:ox+i-1,y=mn?j+1:oy+j-1,z=mn?k+1:oz+k-1;
  const int nx=mn?tx+4:n,ny=mn?ty+4:n,nz=mn?tz+4:n;
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
    double s[27],f[6];
    for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx)
      s[dx+1+3*(dy+1)+9*(dz+1)]=intensity(a,q,x+dx,y+dy,z+dz,nx,ny,nz,b.nq);
    snrt_cweno3_faces(s,f);
    for(int side=0;side<6;++side)faces[q+b.nq*(side+6*t)]=f[side];
  }
}
__global__ void flux(Basis b,int n,int mn,const double *input,const double *a,double *out,
                     double *escaped,double *receipt,double *firstbase,const double *gradient,double cdx,int ox,int oy,int oz,int nx,int ny,int nz,int *failure,
                     double rkweight,double projectweight,double baseextra){
  __shared__ double net[THREADS],loss[THREADS],base[THREADS];
  int tilecell=blockIdx.x,tx=min(TILE,n-ox),ty=min(TILE,n-oy);
  int i=tilecell%tx,j=tilecell/tx%ty,k=tilecell/(tx*ty);
  int gi=ox+i,gj=oy+j,gk=oz+k,cell=gi+n*(gj+n*gk);
  int pos[3]={mn?i+2:gi,mn?j+2:gj,mn?k+2:gk};
  double value=0,projected=0,boundary=0,total_escape=0;
  // Each moment and escape receipt accumulates in the original q order.
  // Shared storage is independent of the total quadrature size.
  for(int q0=0;q0<b.nq;q0+=THREADS){
    const int lane=threadIdx.x,q=q0+lane;
    if(q<b.nq){
    double v=0,e=0;
    for(int ax=0;ax<3;++ax)for(int side=0;side<2;++side){
      int l[3]={pos[0],pos[1],pos[2]},r[3],lm[3],rp[3];l[ax]+=-1+side;
      for(int d=0;d<3;++d){r[d]=l[d];lm[d]=l[d];rp[d]=l[d];}r[ax]++;lm[ax]--;rp[ax]+=2;
      double il=intensity(a,q,l[0],l[1],l[2],nx,ny,nz,b.nq),ir=intensity(a,q,r[0],r[1],r[2],nx,ny,nz,b.nq);
      double dl,dr;
      if(b.limiter==5){
        int gl[3]={i+1,j+1,k+1};gl[ax]+=-1+side;
        const int stride=ax==0?1:(ax==1?tx+2:(tx+2)*(ty+2));
        const int lc=gl[0]+(tx+2)*(gl[1]+(ty+2)*gl[2]),rc=lc+stride;
        il=gradient[q+b.nq*(2*ax+1+6*lc)];ir=gradient[q+b.nq*(2*ax+6*rc)];dl=0;dr=0;
      }else if(b.limiter==4){
        int gl[3]={i+1,j+1,k+1};gl[ax]+=-1+side;
        const int stride=ax==0?1:(ax==1?tx+2:(tx+2)*(ty+2));
        const int lc=gl[0]+(tx+2)*(gl[1]+(ty+2)*gl[2]),rc=lc+stride;
        dl=gradient[q+b.nq*(ax+3*lc)];dr=gradient[q+b.nq*(ax+3*rc)];
      }else{
        dl=limited_slope(il-intensity(a,q,lm[0],lm[1],lm[2],nx,ny,nz,b.nq),ir-il,b.limiter);
        dr=limited_slope(ir-il,intensity(a,q,rp[0],rp[1],rp[2],nx,ny,nz,b.nq)-ir,b.limiter);
      }
      double mu=b.dir[ax+3*q],f=cdx*(fmax(mu,0.)*(il+.5*dl)+fmin(mu,0.)*(ir-.5*dr));
      v+=(1-2*side)*f;int g=ax==0?gi:(ax==1?gj:gk);
      if(side==0&&g==0)e-=b.w[q]*f;if(side==1&&g==n-1)e+=b.w[q]*f;
    }
    net[lane]=v;loss[lane]=e;
    if(mn){
      base[lane]=intensity(a,q,pos[0],pos[1],pos[2],nx,ny,nz,b.nq);
      // Do not clip: at the admitted CFL the angular Euler update is positive.
      if(!isfinite(base[lane]+v)||base[lane]+v<0)atomicCAS(failure,0,cell+1);
    }
    if(!mn)out[q+b.nq*cell]=input[q+b.nq*cell]+v;
    }
    __syncthreads();
    if(mn&&threadIdx.x<b.nm){
      for(int l=0;l<min(THREADS,b.nq-q0);++l){
        const int q=q0+l;const double y=b.y[threadIdx.x+b.nm*q];
        value+=y*(b.w[q]*(base[l]+net[l]));
        projected+=y*(b.w[q]*base[l]);boundary+=y*loss[l];
      }
    }
    if(threadIdx.x==0)for(int l=0;l<min(THREADS,b.nq-q0);++l)total_escape+=loss[l];
    __syncthreads();
  }
  if(mn&&threadIdx.x<b.nm){
    const int j=threadIdx.x;
    const double old=input[j+b.nm*cell],correction=projected-old;
    out[j+b.nm*cell]=value;
    // Also use the first positive projection as SSPRK's retained base state.
    // Its correction has total weight 1; the second Euler correction has 1/2.
    if(firstbase)firstbase[j+b.nm*cell]=projected;
    receipt[j+b.nm*(0+4*tilecell)]=projectweight*correction;
    receipt[j+b.nm*(1+4*tilecell)]=projectweight*fabs(correction);
    receipt[j+b.nm*(2+4*tilecell)]=rkweight*boundary;
    receipt[j+b.nm*(3+4*tilecell)]=rkweight*(value-old)+baseextra*correction;
  }
  if(threadIdx.x==0)escaped[cell]+=total_escape*rkweight;
}
// Original <=384/minmod flux kept as the numerical compatibility path.
__global__ void flux_legacy(Basis b,int n,int mn,const double *input,const double *a,double *out,
                     double *escaped,double *receipt,double *firstbase,double cdx,int ox,int oy,int oz,int nx,int ny,int nz,int *failure){
  __shared__ double net[INLINE_Q],loss[INLINE_Q],base[INLINE_Q];
  int tilecell=blockIdx.x,tx=min(TILE,n-ox),ty=min(TILE,n-oy);
  int i=tilecell%tx,j=tilecell/tx%ty,k=tilecell/(tx*ty);
  int gi=ox+i,gj=oy+j,gk=oz+k,cell=gi+n*(gj+n*gk);
  int pos[3]={mn?i+2:gi,mn?j+2:gj,mn?k+2:gk};
  for(int q=threadIdx.x;q<b.nq;q+=blockDim.x){
    double v=0,e=0;
    for(int ax=0;ax<3;++ax)for(int side=0;side<2;++side){
      int l[3]={pos[0],pos[1],pos[2]},r[3],lm[3],rp[3];l[ax]+=-1+side;
      for(int d=0;d<3;++d){r[d]=l[d];lm[d]=l[d];rp[d]=l[d];}r[ax]++;lm[ax]--;rp[ax]+=2;
      double il=intensity(a,q,l[0],l[1],l[2],nx,ny,nz,b.nq),ir=intensity(a,q,r[0],r[1],r[2],nx,ny,nz,b.nq);
      double dl=minmod(il-intensity(a,q,lm[0],lm[1],lm[2],nx,ny,nz,b.nq),ir-il);
      double dr=minmod(ir-il,intensity(a,q,rp[0],rp[1],rp[2],nx,ny,nz,b.nq)-ir);
      double mu=b.dir[ax+3*q],f=cdx*(fmax(mu,0.)*(il+.5*dl)+fmin(mu,0.)*(ir-.5*dr));
      v+=(1-2*side)*f;int g=ax==0?gi:(ax==1?gj:gk);
      if(side==0&&g==0)e-=b.w[q]*f;if(side==1&&g==n-1)e+=b.w[q]*f;
    }
    net[q]=v;loss[q]=e;
    if(mn){
      base[q]=intensity(a,q,pos[0],pos[1],pos[2],nx,ny,nz,b.nq);
      if(!isfinite(base[q]+v)||base[q]+v<0)atomicCAS(failure,0,cell+1);
    }
    if(!mn)out[q+b.nq*cell]=input[q+b.nq*cell]+v;
  }
  __syncthreads();
  if(mn)for(int j=threadIdx.x;j<b.nm;j+=blockDim.x){
    double value=0,projected=0,boundary=0;
    for(int q=0;q<b.nq;++q){
      const double y=b.y[j+b.nm*q];
      value+=y*(b.w[q]*(base[q]+net[q]));
      projected+=y*(b.w[q]*base[q]);boundary+=y*loss[q];
    }
    const double old=input[j+b.nm*cell],correction=projected-old;
    out[j+b.nm*cell]=value;
    if(firstbase)firstbase[j+b.nm*cell]=projected;
    receipt[j+b.nm*(0+4*tilecell)]=(firstbase?1.:.5)*correction;
    receipt[j+b.nm*(1+4*tilecell)]=(firstbase?1.:.5)*fabs(correction);
    receipt[j+b.nm*(2+4*tilecell)]=.5*boundary;
    receipt[j+b.nm*(3+4*tilecell)]=.5*(value-old)+(firstbase?.5*correction:0.);
  }
  if(threadIdx.x==0){double v=0;for(int q=0;q<b.nq;++q)v+=loss[q];escaped[cell]+=v*.5;}
}
// Bounded tile scratch and deterministic reduction; no per-grid moment ledger
// arrays, atomics on floating sums, or saved angular field are introduced.
__global__ void reduce_receipt(int nm,int cells,const double *tile,double *total){
  __shared__ double scratch[4];
  const int j=blockIdx.x;
  for(int kind=0;kind<4;++kind){
    double value=0;
    for(int cell=threadIdx.x;cell<cells;cell+=blockDim.x)value+=tile[j+nm*(kind+4*cell)];
    value=block_reduce(value,scratch,false);
    if(threadIdx.x==0)total[j+nm*kind]+=value;
    __syncthreads();
  }
}
__global__ void average(const double *a,const double *b,double *out,size_t count,int mn,int *failure){
  size_t i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<count){double v=.5*(a[i]+b[i]);out[i]=v;if(!isfinite(v)||(!mn&&v<0))atomicCAS(failure,0,1);}
}
__global__ void rk_combine(const double *a,const double *b,double *out,size_t count,double wa,double wb,int mn,int *failure){
  size_t i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i<count){double v=wa*a[i]+wb*b[i];out[i]=v;if(!isfinite(v)||(!mn&&v<0))atomicCAS(failure,0,1);}
}
struct Context {
  int n,nc,nv,mn,tensor,mixed,pass=0;Basis b{};
  double *state=nullptr,*stage=nullptr,*work=nullptr,*result=nullptr,*angular=nullptr,*escape=nullptr;
  double *closure_scratch=nullptr,*gradient=nullptr;
  size_t closure_scratch_count=0,gradient_count=0;
  double *hint=nullptr,*hint_target=nullptr;int *hint_valid=nullptr;
  double *next_hint=nullptr,*next_target=nullptr;int *next_valid=nullptr;
  double *tile_receipt=nullptr,*receipt=nullptr;
  // Freeze seeds for a complete Euler pass: overlapping tile halos must
  // reconstruct identical input states identically for face conservation.
  void swap_hints(){std::swap(hint,next_hint);std::swap(hint_target,next_target);std::swap(hint_valid,next_valid);}
  int *failure=nullptr,*retry_mask=nullptr;unsigned long long *stats=nullptr;size_t bytes=0;cudaStream_t stream=nullptr;
  ~Context(){if(stats){unsigned long long counts[8]={};
    if(cudaMemcpy(counts,stats,sizeof(counts),cudaMemcpyDeviceToHost)==cudaSuccess)
      printf("RT_GPU_COUNTERS primary_nonvacuum=%llu fp64_retry_attempts=%llu failed_retries=%llu tensor_closures=%llu\n",counts[2],counts[0],counts[1],counts[3]);
    printf("RT_GPU_WARM attempts=%llu accepted=%llu cold_retries=%llu iterations=%llu\n",counts[4],counts[5],counts[6],counts[7]);
    cudaFree(stats);}
    cudaFree(hint);cudaFree(hint_target);cudaFree(hint_valid);
    cudaFree(next_hint);cudaFree(next_target);cudaFree(next_valid);
    cudaFree(tile_receipt);cudaFree(receipt);
    cudaFree(closure_scratch);
    cudaFree(gradient);
    cudaFree(retry_mask);cudaFree(const_cast<double*>(b.y));cudaFree(const_cast<double*>(b.w));cudaFree(const_cast<double*>(b.dir));cudaFree(state);cudaFree(stage);cudaFree(work);cudaFree(result);cudaFree(angular);cudaFree(escape);cudaFree(failure);if(stream)cudaStreamDestroy(stream);}
};
bool checked(cudaError_t e,const char *where){if(e==cudaSuccess)return true;fprintf(stderr,"RT CUDA %s: %s\n",where,cudaGetErrorString(e));return false;}
bool allocate(double *&p,size_t count){return checked(cudaMalloc(&p,count*sizeof(double)),"allocate");}
template<bool Global> void launch_closure_impl(Context &c,const double *in,double *angular,int ox,int oy,int oz,
                    int nx,int ny,int nz,int offset,int count){
  int blocks=offset>=0?count:nx*ny*nz;
  if(c.mixed){
    reconstruct<float,Global><<<blocks,THREADS,sizeof(Workspace<float,Global>),c.stream>>>
      (c.b,c.n,in,angular,ox,oy,oz,nx,ny,nz,offset,count,c.tensor,c.failure,2e-6,c.retry_mask,0,c.stats,c.hint,c.hint_target,c.hint_valid,c.next_hint,c.next_target,c.next_valid,c.pass,c.closure_scratch);
    // GPU-only adaptive retry. Same declared mixed tolerance and physical
    // target; only failed cells use FP64 arithmetic, never CPU fallback.
    reconstruct<double,Global><<<blocks,THREADS,sizeof(Workspace<double,Global>),c.stream>>>
      (c.b,c.n,in,angular,ox,oy,oz,nx,ny,nz,offset,count,0,c.failure,2e-6,c.retry_mask,1,c.stats,c.hint,c.hint_target,c.hint_valid,c.next_hint,c.next_target,c.next_valid,c.pass,c.closure_scratch);
  }else{
    reconstruct<double,Global><<<blocks,THREADS,sizeof(Workspace<double,Global>),c.stream>>>
      (c.b,c.n,in,angular,ox,oy,oz,nx,ny,nz,offset,count,c.tensor,c.failure,2e-12,nullptr,0,c.stats,c.hint,c.hint_target,c.hint_valid,c.next_hint,c.next_target,c.next_valid,c.pass,c.closure_scratch);
  }
}
void launch_closure(Context &c,const double *in,double *angular,int ox,int oy,int oz,
                    int nx,int ny,int nz,int offset,int count){
  if(c.b.nq<=INLINE_Q)launch_closure_impl<false>(c,in,angular,ox,oy,oz,nx,ny,nz,offset,count);
  else launch_closure_impl<true>(c,in,angular,ox,oy,oz,nx,ny,nz,offset,count);
}
int euler(Context &c,const double *in,double *out,double cdx){
  const double rkweight=c.b.limiter==5?(c.pass==3?2./3:1./6):.5;
  const double projectweight=c.pass==1?1.:rkweight;
  const double baseextra=c.pass==1?1.-rkweight:0.;
  for(int z=0;z<c.n;z+=TILE)for(int y=0;y<c.n;y+=TILE)for(int x=0;x<c.n;x+=TILE){
    int tx=std::min(TILE,c.n-x),ty=std::min(TILE,c.n-y),tz=std::min(TILE,c.n-z);
    if(c.mn)launch_closure(c,in,c.angular,x,y,z,tx+4,ty+4,tz+4,-1,0);
    if(c.b.limiter==4)vertex_gradients<<<(tx+2)*(ty+2)*(tz+2),THREADS,0,c.stream>>>
      (c.b,c.n,c.mn,c.mn?c.angular:in,c.gradient,x,y,z,tx,ty,tz);
    if(c.b.limiter==5)cweno_faces<<<(tx+2)*(ty+2)*(tz+2),THREADS,0,c.stream>>>
      (c.b,c.n,c.mn,c.mn?c.angular:in,c.gradient,x,y,z,tx,ty,tz);
    if(c.b.nq<=INLINE_Q&&c.b.limiter==1)
      flux_legacy<<<tx*ty*tz,THREADS,0,c.stream>>>(c.b,c.n,c.mn,in,c.mn?c.angular:in,out,c.escape,c.tile_receipt,(c.mn&&c.pass==1)?c.result:nullptr,cdx,x,y,z,c.mn?tx+4:c.n,c.mn?ty+4:c.n,c.mn?tz+4:c.n,c.failure);
    else
      flux<<<tx*ty*tz,THREADS,0,c.stream>>>(c.b,c.n,c.mn,in,c.mn?c.angular:in,out,c.escape,c.tile_receipt,(c.mn&&c.pass==1)?c.result:nullptr,c.gradient,cdx,x,y,z,c.mn?tx+4:c.n,c.mn?ty+4:c.n,c.mn?tz+4:c.n,c.failure,rkweight,projectweight,baseextra);
    if(c.mn)reduce_receipt<<<c.b.nm,THREADS,0,c.stream>>>(c.b.nm,tx*ty*tz,c.tile_receipt,c.receipt);
  }
  if(c.mn)c.swap_hints();
  return checked(cudaGetLastError(),"Euler kernels")?0:90;
}
} // namespace

extern "C" void *rt_gpu_create(int n,int nm,int nq,int mn,int mode,int limiter,const double *y,const double *w,const double *dir){
  if(n<2||n>128||nm<4||nm>36||nq<1||nq>Q||mode<1||mode>4||(!mn&&mode!=1)||limiter<1||limiter>5)return nullptr;
  cudaDeviceProp prop{};int device=0;
  if(!checked(cudaGetDevice(&device),"device")||!checked(cudaGetDeviceProperties(&prop,device),"properties"))return nullptr;
  if(prop.major<8){fprintf(stderr,"RT GPU benchmark requires sm80 or newer\n");return nullptr;}
  // Explicit supported FP64 tensor architectures, not all Ampere/Ada GPUs.
  if(mode==2&&!((prop.major==8&&prop.minor==0)||prop.major==9)){
    fprintf(stderr,"RT Tensor backend requires FP64 WMMA (validated targets sm80/sm90), device=%s\n",prop.name);return nullptr;
  }
  if(mode==4&&prop.major<8)return nullptr;
  Context *c=new(std::nothrow)Context;if(!c)return nullptr;
  c->n=n;c->nc=n*n*n;c->nv=mn?nm:nq;c->mn=mn;c->tensor=(mode==2||mode==4);c->mixed=mode>=3;c->b.nm=nm;c->b.nq=nq;c->b.limiter=limiter;
  c->bytes=size_t(c->nc)*c->nv*sizeof(double);
  bool ok=checked(cudaStreamCreate(&c->stream),"stream")&&
    checked(cudaFuncSetAttribute(reconstruct<double,false>,cudaFuncAttributeMaxDynamicSharedMemorySize,sizeof(Workspace<double,false>)),"small closure shared")&&
    checked(cudaFuncSetAttribute(reconstruct<float,false>,cudaFuncAttributeMaxDynamicSharedMemorySize,sizeof(Workspace<float,false>)),"small mixed shared")&&
    checked(cudaFuncSetAttribute(reconstruct<double,true>,cudaFuncAttributeMaxDynamicSharedMemorySize,sizeof(Workspace<double,true>)),"large closure shared")&&
    checked(cudaFuncSetAttribute(reconstruct<float,true>,cudaFuncAttributeMaxDynamicSharedMemorySize,sizeof(Workspace<float,true>)),"large mixed shared")&&
    allocate(c->state,size_t(c->nc)*c->nv)&&allocate(c->stage,size_t(c->nc)*c->nv)&&
    allocate(c->work,size_t(c->nc)*c->nv)&&allocate(c->result,size_t(c->nc)*c->nv)&&
    allocate(c->escape,c->nc)&&checked(cudaMalloc(&c->failure,sizeof(int)),"failure allocation")&&
    checked(cudaMalloc(&c->retry_mask,c->nc*sizeof(int)),"retry mask")&&
    checked(cudaMalloc(&c->stats,8*sizeof(unsigned long long)),"counters")&&
    checked(cudaMemset(c->stats,0,8*sizeof(unsigned long long)),"clear counters");
  double *dy=nullptr,*dw=nullptr,*dd=nullptr;
  ok=ok&&allocate(dy,nm*nq)&&allocate(dw,nq)&&allocate(dd,3*nq);
  c->b.y=dy;c->b.w=dw;c->b.dir=dd;
  if(mn)ok=ok&&allocate(c->tile_receipt,size_t(TILE)*TILE*TILE*4*nm)&&allocate(c->receipt,4*nm)&&
    allocate(c->angular,size_t(HALO)*HALO*HALO*nq)&&
    allocate(c->hint,size_t(c->nc)*(nm-1))&&allocate(c->hint_target,size_t(c->nc)*(nm-1))&&
    checked(cudaMalloc(&c->hint_valid,c->nc*sizeof(int)),"warm validity")&&
    checked(cudaMemset(c->hint_valid,0,c->nc*sizeof(int)),"clear warm validity")&&
    allocate(c->next_hint,size_t(c->nc)*(nm-1))&&allocate(c->next_target,size_t(c->nc)*(nm-1))&&
    checked(cudaMalloc(&c->next_valid,c->nc*sizeof(int)),"next warm validity")&&
    checked(cudaMemset(c->next_valid,0,c->nc*sizeof(int)),"clear next warm validity");
  if(mn&&nq>INLINE_Q){
    const size_t edge=std::min(n,TILE)+4;
    // Also covers the final closure's <=min(n^3,8000) contiguous-cell batch.
    c->closure_scratch_count=edge*edge*edge*5*nq;
    ok=ok&&allocate(c->closure_scratch,c->closure_scratch_count);
  }
  if(limiter==4||limiter==5){
    const size_t edge=std::min(n,TILE)+2;
    c->gradient_count=edge*edge*edge*(limiter==5?6:3)*nq;
    ok=ok&&allocate(c->gradient,c->gradient_count);
  }
  if(ok)ok=checked(cudaMemcpy(dy,y,nm*nq*sizeof(double),cudaMemcpyHostToDevice),"basis")&&
    checked(cudaMemcpy(dw,w,nq*sizeof(double),cudaMemcpyHostToDevice),"weights")&&
    checked(cudaMemcpy(dd,dir,3*nq*sizeof(double),cudaMemcpyHostToDevice),"directions");
  if(!ok){delete c;return nullptr;}
  const char *names[]={"","cuda_fp64","tensor_fp64","cuda_mixed","tensor_tf32x3_mixed"};
  const size_t shared=nq<=INLINE_Q?(c->mixed?sizeof(Workspace<float,false>):sizeof(Workspace<double,false>)):
    (c->mixed?sizeof(Workspace<float,true>):sizeof(Workspace<double,true>));
  printf("RT_GPU device=%s cc=%d.%d mode=%s closure_tolerance=%.3g transport=FP64 shared_per_closure=%zu state_bytes=%zu angular_scratch_bytes=%zu\n",prop.name,prop.major,prop.minor,names[mode],c->mixed?2e-6:2e-12,shared,c->bytes,mn?size_t(HALO)*HALO*HALO*nq*8:0);
  if(mn)printf("RT_GPU compact_warm_cache_bytes=%zu frozen_stage_seeds=1 cooperative_cholesky=1\n",2*size_t(c->nc)*(2*(nm-1)*sizeof(double)+sizeof(int)));
  if(mn)printf("RT_MN_UPDATE kinetic_projection=1 explicit_projection_receipt=1 receipt_scratch_bytes=%zu\n",(size_t(TILE)*TILE*TILE+1)*4*nm*sizeof(double));
  printf("RT_ANGULAR nodes=%d limiter=%s closure_tile_scratch_bytes=%zu streamed_face_angles=%d\n",nq,limiter==5?"cweno3":(limiter==4?"mlp":(limiter==3?"mc15":(limiter==2?"mc":"minmod"))),c->closure_scratch_count*sizeof(double),c->b.nq<=INLINE_Q&&limiter==1?0:THREADS);
  if(limiter==4)printf("RT_SPATIAL vertex_mlp=1 positive_corner_fraction=0.25 gradient_tile_scratch_bytes=%zu\n",c->gradient_count*sizeof(double));
  if(limiter==5)printf("RT_SPATIAL cweno3=1 bernstein_positive=1 exact_quadratic_face_average=1 ssprk=3 face_tile_scratch_bytes=%zu\n",c->gradient_count*sizeof(double));
  return c;
}
extern "C" int rt_gpu_transport(void *handle,const double *state,double cdx,double *out,double *escaped,double *receipt){
  auto &c=*static_cast<Context*>(handle);
  if(!std::isfinite(cdx)||cdx<0||cdx>1./12)return 1;
  if(cdx==0){std::copy(state,state+size_t(c.nc)*c.nv,out);*escaped=0;std::fill(receipt,receipt+4*c.b.nm,0.);return 0;}
  if(!checked(cudaMemcpyAsync(c.state,state,c.bytes,cudaMemcpyHostToDevice,c.stream),"input")||
     !checked(cudaMemsetAsync(c.failure,0,sizeof(int),c.stream),"clear failure")||
     !checked(cudaMemsetAsync(c.escape,0,c.nc*sizeof(double),c.stream),"clear escape"))return 90;
  if(c.mn&&!checked(cudaMemsetAsync(c.receipt,0,4*c.b.nm*sizeof(double),c.stream),"clear receipt"))return 90;
  c.pass=1;int status=euler(c,c.state,c.stage,cdx);if(status)return status;
  c.pass=2;status=euler(c,c.stage,c.work,cdx);if(status)return status;
  size_t count=size_t(c.nc)*c.nv;
  if(c.b.limiter==5){
    rk_combine<<<(count+255)/256,256,0,c.stream>>>(c.mn?c.result:c.state,c.work,c.stage,count,.75,.25,c.mn,c.failure);
    c.pass=3;status=euler(c,c.stage,c.work,cdx);if(status)return status;
    rk_combine<<<(count+255)/256,256,0,c.stream>>>(c.mn?c.result:c.state,c.work,c.result,count,1./3,2./3,c.mn,c.failure);
  }else average<<<(count+255)/256,256,0,c.stream>>>(c.mn?c.result:c.state,c.work,c.result,count,c.mn,c.failure);
  c.pass=c.b.limiter==5?4:3;
  if(c.mn)for(int offset=0;offset<c.nc;offset+=8000){
    launch_closure(c,c.result,nullptr,0,0,0,0,0,0,offset,std::min(8000,c.nc-offset));
  }
  if(c.mn)c.swap_hints();
  int failed=0;
  if(!checked(cudaMemcpyAsync(&failed,c.failure,sizeof(int),cudaMemcpyDeviceToHost,c.stream),"failure copy")||
     !checked(cudaStreamSynchronize(c.stream),"transport synchronize"))return 90;
  if(failed){fprintf(stderr,"RT_GPU failed cell/index=%d; no CPU fallback, no output commit\n",failed);return 2;}
  std::vector<double> loss(c.nc);
  if(!checked(cudaMemcpy(out,c.result,c.bytes,cudaMemcpyDeviceToHost),"output")||
     !checked(cudaMemcpy(loss.data(),c.escape,c.nc*sizeof(double),cudaMemcpyDeviceToHost),"escape"))return 90;
  if(c.mn){if(!checked(cudaMemcpy(receipt,c.receipt,4*c.b.nm*sizeof(double),cudaMemcpyDeviceToHost),"projection receipt"))return 90;}
  else std::fill(receipt,receipt+4*c.b.nm,0.);
  *escaped=0;for(double v:loss)*escaped+=v;return 0;
}
extern "C" void rt_gpu_destroy(void *handle){delete static_cast<Context*>(handle);}

#ifdef MN_LIVE_POOL
#include "../cuRamses/cuda_stream_pool.h"
// Live AMR adapters use the shared multigpu stream lease, not the benchmark
// context or its uniform mesh. MPI/face ownership remains on the host.
namespace {
struct LiveBuffers {
  cudaStream_t stream;std::vector<void*> pointers;
  explicit LiveBuffers(int slot):stream(cuda_get_stream_internal(slot)){}
  template<class T> bool alloc(T *&p,size_t n){
    if(cudaMallocAsync(&p,n*sizeof(T),stream)!=cudaSuccess)return false;
    pointers.push_back(p);return true;
  }
  ~LiveBuffers(){for(void *p:pointers)cudaFreeAsync(p,stream);cudaStreamSynchronize(stream);}
};
__global__ void live_project(Basis b,int nc,const double *angular,double *moments){
  int c=blockIdx.x,j=threadIdx.x;
  if(c>=nc||j>=b.nm)return;
  double v=0;for(int q=0;q<b.nq;++q)v+=b.y[j+b.nm*q]*b.w[q]*angular[c*b.nq+q];
  moments[c*b.nm+j]=v;
}
__global__ void live_face_flux(int nf,int nm,int nq,const double *y,const double *w,const double *direction,
    const double *left,const double *right,const double *geometry,double *flux){
  int f=blockIdx.x,j=threadIdx.x;
  if(f>=nf||j>=nm)return;
  int axis=int(geometry[5*f])-1;double sign=geometry[5*f+1],value=0;
  for(int q=0;q<nq;++q){
    const int k=4*(q+nq*f);double mu=sign*direction[axis+3*q];
    double il=left[k]+sign*.5*geometry[5*f+2]*left[k+axis+1];
    double ir=right[k]-sign*.5*geometry[5*f+3]*right[k+axis+1];
    double flow=geometry[5*f+4]*(fmax(mu,0.)*il+fmin(mu,0.)*ir);
    value+=y[j+nm*q]*(w[q]*flow);
  }
  flux[j+nm*f]=value;
}
}
static int snrt_mn_cuda_closure_impl(int slot,int nc,int nm,int nq,const double *y,const double *w,
    const double *direction,const double *input,double *projected,double *cache,double *angular_host){
  if(slot<0||nc<1||nc>256||nm<9||nm>36||nq!=384)return 7;
  LiveBuffers mem(slot);double *dy=nullptr,*dw=nullptr,*dd=nullptr,*du=nullptr,*da=nullptr,*dp=nullptr,*dc=nullptr;
  int *failure=nullptr;unsigned long long *stats=nullptr;int failed=0;
  if(!mem.alloc(dy,nm*nq)||!mem.alloc(dw,nq)||!mem.alloc(dd,3*nq)||!mem.alloc(du,nm*nc)||
     !mem.alloc(da,nq*nc)||!mem.alloc(dp,nm*nc)||!mem.alloc(dc,(nm+1)*nc)||!mem.alloc(failure,1)||!mem.alloc(stats,8))return 7;
  auto copy=[&](double *d,const double *h,size_t n){return cudaMemcpyAsync(d,h,8*n,cudaMemcpyHostToDevice,mem.stream)==cudaSuccess;};
  if(!copy(dy,y,nm*nq)||!copy(dw,w,nq)||!copy(dd,direction,3*nq)||!copy(du,input,nm*nc))return 7;
  if(cudaMemsetAsync(failure,0,sizeof(int),mem.stream)!=cudaSuccess||
     cudaMemsetAsync(stats,0,8*sizeof(unsigned long long),mem.stream)!=cudaSuccess||
     cudaMemsetAsync(dc,0,8*(nm+1)*nc,mem.stream)!=cudaSuccess)return 7;
  if(cudaFuncSetAttribute(reconstruct<double,false>,cudaFuncAttributeMaxDynamicSharedMemorySize,sizeof(Workspace<double,false>))!=cudaSuccess)return 7;
  Basis b{nm,nq,3,dy,dw,dd};
  reconstruct<double,false><<<nc,THREADS,sizeof(Workspace<double,false>),mem.stream>>>(b,1,du,da,0,0,0,0,0,0,
      0,nc,0,failure,2e-12,nullptr,0,stats,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,0,nullptr,dc);
  if(cudaGetLastError()!=cudaSuccess)return 7;
  if(cudaMemcpyAsync(&failed,failure,sizeof(int),cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess||cudaStreamSynchronize(mem.stream)!=cudaSuccess)return 7;
  if(failed)return 2;
  live_project<<<nc,64,0,mem.stream>>>(b,nc,da,dp);
  if(cudaGetLastError()!=cudaSuccess)return 7;
  if(cudaMemcpyAsync(projected,dp,8*nm*nc,cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess||
     (angular_host&&cudaMemcpyAsync(angular_host,da,8*nq*nc,cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess)||
     cudaMemcpyAsync(cache,dc,8*(nm+1)*nc,cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess||cudaStreamSynchronize(mem.stream)!=cudaSuccess)return 7;
  return 0;
}
extern "C" int snrt_mn_cuda_closure_c(int slot,int nc,int nm,int nq,const double *y,const double *w,
    const double *direction,const double *input,double *projected,double *cache){
  return snrt_mn_cuda_closure_impl(slot,nc,nm,nq,y,w,direction,input,projected,cache,nullptr);
}
extern "C" int snrt_mn_cuda_closure_angular_c(int slot,int nc,int nm,int nq,const double *y,const double *w,
    const double *direction,const double *input,double *projected,double *cache,double *angular){
  return snrt_mn_cuda_closure_impl(slot,nc,nm,nq,y,w,direction,input,projected,cache,angular);
}
extern "C" int snrt_mn_cuda_project_c(int slot,int nc,int nm,int nq,const double *y,const double *w,
    const double *angular,double *moments){
  if(slot<0||nc<1||nc>256||nm<1||nm>36||nq!=384)return 7;
  LiveBuffers mem(slot);double *dy=nullptr,*dw=nullptr,*da=nullptr,*dm=nullptr;
  if(!mem.alloc(dy,nm*nq)||!mem.alloc(dw,nq)||!mem.alloc(da,nq*nc)||!mem.alloc(dm,nm*nc))return 7;
  if(cudaMemcpyAsync(dy,y,8*nm*nq,cudaMemcpyHostToDevice,mem.stream)!=cudaSuccess||
     cudaMemcpyAsync(dw,w,8*nq,cudaMemcpyHostToDevice,mem.stream)!=cudaSuccess||
     cudaMemcpyAsync(da,angular,8*nq*nc,cudaMemcpyHostToDevice,mem.stream)!=cudaSuccess)return 7;
  Basis b{nm,nq,3,dy,dw,nullptr};
  live_project<<<nc,64,0,mem.stream>>>(b,nc,da,dm);
  if(cudaGetLastError()!=cudaSuccess||
     cudaMemcpyAsync(moments,dm,8*nm*nc,cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess||
     cudaStreamSynchronize(mem.stream)!=cudaSuccess)return 7;
  return 0;
}
extern "C" int snrt_mn_cuda_flux_c(int slot,int nf,int nm,int nq,const double *y,const double *w,
    const double *direction,const double *left,const double *right,const double *geometry,double *flux){
  if(slot<0||nf<1||nf>256||nm<9||nm>36||nq<1||nq>8)return 7;
  LiveBuffers mem(slot);double *dy=nullptr,*dw=nullptr,*dd=nullptr,*dl=nullptr,*dr=nullptr,*dg=nullptr,*df=nullptr;
  if(!mem.alloc(dy,nm*nq)||!mem.alloc(dw,nq)||!mem.alloc(dd,3*nq)||!mem.alloc(dl,4*nq*nf)||
     !mem.alloc(dr,4*nq*nf)||!mem.alloc(dg,5*nf)||!mem.alloc(df,nm*nf))return 7;
  auto copy=[&](double *d,const double *h,size_t n){return cudaMemcpyAsync(d,h,8*n,cudaMemcpyHostToDevice,mem.stream)==cudaSuccess;};
  if(!copy(dy,y,nm*nq)||!copy(dw,w,nq)||!copy(dd,direction,3*nq)||!copy(dl,left,4*nq*nf)||
     !copy(dr,right,4*nq*nf)||!copy(dg,geometry,5*nf))return 7;
  live_face_flux<<<nf,64,0,mem.stream>>>(nf,nm,nq,dy,dw,dd,dl,dr,dg,df);
  if(cudaGetLastError()!=cudaSuccess)return 7;
  if(cudaMemcpyAsync(flux,df,8*nm*nf,cudaMemcpyDeviceToHost,mem.stream)!=cudaSuccess||cudaStreamSynchronize(mem.stream)!=cudaSuccess)return 7;
  return 0;
}
#endif

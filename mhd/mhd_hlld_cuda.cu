// NENER=0 HLLD: numerical counterpart of mhd/godunov_utils.f90::hlld.
// Miyoshi & Kusano (2005). Reconstruction/edge EMF/CT remain on the CPU.
#include "../patch/cuRamses/cuda_stream_pool.h"
#include <cmath>
#include <cstdio>

struct State { double r,u,v,w,b,c,p,e,vb,ei; };
__device__ static State primitive(const double* q,double a,double gamma) {
    State s;
    s.r=q[0]; s.u=q[2]; s.v=q[4]; s.w=q[6]; s.b=q[5]; s.c=q[7];
    const double em=.5*(a*a+s.b*s.b+s.c*s.c);
    s.ei=q[1]/(gamma-1.);
    s.e=s.ei+.5*(s.u*s.u+s.v*s.v+s.w*s.w)*s.r+em;
    s.p=q[1]+em; s.vb=s.u*a+s.v*s.b+s.w*s.c;
    return s;
}
__device__ static double fast(const double* q,double a,double gamma) {
    const double b2=a*a+q[5]*q[5]+q[7]*q[7],c2=gamma*q[1]/q[0];
    const double d2=.5*(b2/q[0]+c2);
    return sqrt(d2+sqrt(d2*d2-c2*a*a/q[0]));
}
__device__ static State star(State s,double speed,double us,double ps,double a) {
    State t=s;
    t.r=s.r*(speed-s.u)/(speed-us); t.u=us; t.p=ps;
    const double denom=s.r*(speed-s.u)*(speed-us)-a*a;
    const double numer=s.r*(speed-s.u)*(speed-s.u)-a*a;
    if(!(fabs(denom)<1.e-4*a*a)) {
        t.v=s.v-a*s.b*(us-s.u)/denom; t.b=s.b*numer/denom;
        t.w=s.w-a*s.c*(us-s.u)/denom; t.c=s.c*numer/denom;
    }
    t.vb=us*a+t.v*t.b+t.w*t.c;
    t.e=((speed-s.u)*s.e-s.p*s.u+ps*us+a*(s.vb-t.vb))/(speed-us);
    t.ei=s.ei*(speed-s.u)/(speed-us);
    return t;
}
__global__ static void face_kernel(const double* l,const double* r,const int* mask,
                                   double* flux,int n,int nv,double gamma) {
    const int j=blockIdx.x*blockDim.x+threadIdx.x;
    if(j>=n || !mask[j])return;
    l+=static_cast<size_t>(j)*nv; r+=static_cast<size_t>(j)*nv;
    double* f=flux+static_cast<size_t>(j)*(nv+1);
    const double a=.5*(l[3]+r[3]),sign= a<0 ? -1. : 1.;
    const State L=primitive(l,a,gamma),R=primitive(r,a,gamma);
    const double cf=fmax(fast(l,a,gamma),fast(r,a,gamma));
    const double sl=fmin(L.u,R.u)-cf,sr=fmax(L.u,R.u)+cf;
    const double rl=L.r*(L.u-sl),rr=R.r*(sr-R.u);
    const double us=(rr*R.u+rl*L.u+(L.p-R.p))/(rr+rl);
    const double ps=(rr*L.p+rl*R.p+rl*rr*(L.u-R.u))/(rr+rl);
    const State LS=star(L,sl,us,ps,a),RS=star(R,sr,us,ps,a);
    const double ql=sqrt(LS.r),qr=sqrt(RS.r);
    const double sal=us-fabs(a)/ql,sar=us+fabs(a)/qr;
    const double vss=(ql*LS.v+qr*RS.v+sign*(RS.b-LS.b))/(ql+qr);
    const double wss=(ql*LS.w+qr*RS.w+sign*(RS.c-LS.c))/(ql+qr);
    const double bss=(ql*RS.b+qr*LS.b+sign*ql*qr*(RS.v-LS.v))/(ql+qr);
    const double css=(ql*RS.c+qr*LS.c+sign*ql*qr*(RS.w-LS.w))/(ql+qr);
    const double vbss=us*a+vss*bss+wss*css;
    State S;
    if(sl>0.) S=L;
    else if(sal>0.) S=LS;
    else if(us>0.) {
        S=LS;S.v=vss;S.w=wss;S.b=bss;S.c=css;S.vb=vbss;
        S.e=LS.e-sign*ql*(LS.vb-vbss);
    } else if(sar>0.) {
        S=RS;S.v=vss;S.w=wss;S.b=bss;S.c=css;S.vb=vbss;
        S.e=RS.e+sign*qr*(RS.vb-vbss);
    } else if(sr>0.) S=RS;
    else S=R;
    f[0]=S.r*S.u;f[1]=(S.e+S.p)*S.u-a*S.vb;
    f[2]=S.r*S.u*S.u+S.p-a*a;f[3]=0.;
    f[4]=S.r*S.u*S.v-a*S.b;f[5]=S.b*S.u-a*S.v;
    f[6]=S.r*S.u*S.w-a*S.c;f[7]=S.c*S.u-a*S.w;
    for(int k=8;k<nv;++k)f[k]=f[0]*(f[0]>0.?l[k]:r[k]);
    f[nv]=S.u*S.ei;
}

// 1=completed device work; 0=pool absent/busy (caller computes on CPU);
// -1=actual CUDA error (caller must stop, never consume partial output).
extern "C" int mhd_hlld_batch(const double* l,const double* r,const int* mask,
                               double* flux,int n,int nv,double gamma) {
    if(n<=0 || nv<8 || nv!=NVAR)return -1;
    const int slot=cuda_acquire_stream();
    if(slot<0)return 0;
    StreamSlot& s=get_pool()[slot];
    cudaError_t err=cudaSuccess;
    auto check=[&](cudaError_t e){if(err==cudaSuccess)err=e;return e==cudaSuccess;};
    if(n>s.mhd_face_cap || nv!=s.mhd_nvar) {
        if(s.d_mhd_left)check(cudaFree(s.d_mhd_left));
        if(s.d_mhd_right)check(cudaFree(s.d_mhd_right));
        if(s.d_mhd_flux)check(cudaFree(s.d_mhd_flux));
        if(s.d_mhd_mask)check(cudaFree(s.d_mhd_mask));
        s.d_mhd_left=s.d_mhd_right=s.d_mhd_flux=nullptr;s.d_mhd_mask=nullptr;
        s.mhd_face_cap=0;s.mhd_nvar=0;
        const size_t count=static_cast<size_t>(n)*nv;
        check(cudaMalloc(&s.d_mhd_left,count*sizeof(double)));
        check(cudaMalloc(&s.d_mhd_right,count*sizeof(double)));
        check(cudaMalloc(&s.d_mhd_flux,static_cast<size_t>(n)*(nv+1)*sizeof(double)));
        check(cudaMalloc(&s.d_mhd_mask,static_cast<size_t>(n)*sizeof(int)));
        if(err==cudaSuccess){s.mhd_face_cap=n;s.mhd_nvar=nv;}
    }
    const size_t qbytes=static_cast<size_t>(n)*nv*sizeof(double);
    const size_t fbytes=static_cast<size_t>(n)*(nv+1)*sizeof(double);
    if(err==cudaSuccess)check(cudaMemcpyAsync(s.d_mhd_left,l,qbytes,cudaMemcpyHostToDevice,s.stream));
    if(err==cudaSuccess)check(cudaMemcpyAsync(s.d_mhd_right,r,qbytes,cudaMemcpyHostToDevice,s.stream));
    if(err==cudaSuccess)check(cudaMemcpyAsync(s.d_mhd_mask,mask,n*sizeof(int),cudaMemcpyHostToDevice,s.stream));
    if(err==cudaSuccess) {
        face_kernel<<<(n+127)/128,128,0,s.stream>>>(s.d_mhd_left,s.d_mhd_right,s.d_mhd_mask,s.d_mhd_flux,n,nv,gamma);
        check(cudaGetLastError());
    }
    if(err==cudaSuccess)check(cudaMemcpyAsync(flux,s.d_mhd_flux,fbytes,cudaMemcpyDeviceToHost,s.stream));
    check(cudaStreamSynchronize(s.stream));
    cuda_release_stream(slot);
    if(err!=cudaSuccess){fprintf(stderr,"MHD CUDA error: %s\n",cudaGetErrorString(err));return -1;}
    return 1;
}

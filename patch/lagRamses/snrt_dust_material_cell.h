#pragma once
#include <cmath>
#include <cfloat>
#ifdef __CUDACC__
#define DUST_HD __host__ __device__
#else
#define DUST_HD
#endif

// Same FP64 scalar material solve on both processors. The Fortran solver
// remains the independent reference. No gas, radiation or persistent writes.
DUST_HD inline double dust_material_u(const double *a,int nt,int use_u,
    double temperature,double density,double capacity) {
  if(!use_u) return capacity*temperature;
  const double x=log(temperature); int k=0;
  while(k<nt-2 && x>a[k+1]) ++k;
  const double w=fmax(0.,fmin(1.,(x-a[k])/(a[k+1]-a[k])));
  return density*(a[2*nt+k]+w*(a[2*nt+k+1]-a[2*nt+k]));
}

DUST_HD inline int dust_material_cell(const double *input,const double *a,double *out,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int i) {
  const double heating=input[i],density=input[nc+i],old=input[2*nc+i],capacity=input[3*nc+i];
  if(!isfinite(heating)||!isfinite(density)||!isfinite(old)||!isfinite(capacity)||
      heating<0||density<0||old<0||capacity<=0) return 2;
  const double *p=a+nt,*band=a+3*nt;
  double *rate=out+size_t(i)*ng;
  for(int g=0;g<ng;++g)rate[g]=0;
  if(density==0) {
    if(heating!=0 || (use_u && old!=0)) return 5;
    out[size_t(ng)*nc+i]=use_u?0:old/capacity;
    out[size_t(ng+1)*nc+i]=old;
    return 0;
  }
  const double target=old+dt*heating;
  const double floor=dust_material_u(a,nt,use_u,bath,density,capacity);
  const double top=dust_material_u(a,nt,use_u,exp(a[nt-1]),density,capacity);
  if(!isfinite(target)||target<floor*(1-tolerance)||
      target>(top+dt*density*(p[nt-1]-background))*(1+64*DBL_EPSILON))return 5;
  double lower=0,upper=fmin(density*(p[nt-1]-background),fmax(heating+(old-floor)/dt,0.));
  for(int it=0;it<80;++it) {
    const double mid=lower+0.5*(upper-lower),power=background+mid/density;
    int k=0;while(k<nt-2 && power>p[k+1])++k;
    const double w=(power-p[k])/(p[k+1]-p[k]);
    const double temperature=exp(a[k]+w*(a[k+1]-a[k]));
    const double residual=(dust_material_u(a,nt,use_u,temperature,density,capacity)-old)/dt+mid-heating;
    if(residual>0)upper=mid;else lower=mid;
  }
  const double emitted=lower+0.5*(upper-lower),increment=emitted/density;
  const double power=background+increment;
  int k=0;while(k<nt-2 && power>p[k+1])++k;
  const double w=(power-p[k])/(p[k+1]-p[k]);
  const double temperature=exp(a[k]+w*(a[k+1]-a[k]));
  out[size_t(ng)*nc+i]=temperature;
  out[size_t(ng+1)*nc+i]=dust_material_u(a,nt,use_u,temperature,density,capacity);
  if(!isfinite(increment)||increment>p[nt-1]-background)return 5;
  for(int t=0;t<nt-1;++t) {
    const double start=fmax(p[t]-background,0.),finish=fmax(p[t+1]-background,0.);
    const double width=fmax(fmin(increment,finish)-start,0.);
    for(int g=0;g<ng;++g)
      rate[g]+=width*(band[size_t(t+1)*ng+g]-band[size_t(t)*ng+g])/(p[t+1]-p[t])*density;
  }
  for(int g=0;g<ng;++g)if(!isfinite(rate[g])||rate[g]<0)return 2;
  return isfinite(temperature)&&isfinite(out[size_t(ng+1)*nc+i])?0:2;
}
#undef DUST_HD

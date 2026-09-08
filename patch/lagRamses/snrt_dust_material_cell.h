#pragma once
#include <cmath>
#include <cfloat>
#ifdef __CUDACC__
#define DUST_HD __host__ __device__
#else
#define DUST_HD
#endif

// Same FP64 scalar material solve on both processors. The Fortran solver
// remains the no-exchange independent reference. No persistent writes.
DUST_HD inline double dust_material_u(const double *a,int nt,int use_u,
    double temperature,double density,double capacity) {
  if(!use_u) return capacity*temperature;
  const double x=log(temperature); int k=0;
  while(k<nt-2 && x>a[k+1]) ++k;
  const double w=fmax(0.,fmin(1.,(x-a[k])/(a[k+1]-a[k])));
  return density*(a[2*nt+k]+w*(a[2*nt+k+1]-a[2*nt+k]));
}

// Eliminate the BE gas equation with K(T)=K0*sqrt(T/T0), not a frozen
// thermal speed. Cv and the geometric/accommodation factors stay fixed.
// In y=sqrt(T/max(T0,Td)), the scaled cubic is
// b*(y*y-t) + a*y*(y*y-d)=0. It has one positive root for T0,Td>0.
// Safeguarded Newton starts at the upper bracket (the cubic is convex for
// y>=0). Scaling a,b separately keeps both weak and stiff limits finite.
DUST_HD inline double dust_gas_transfer(double gas,double cv,double kd,double td) {
  if(kd==0 || gas==cv*td)return 0;
  if(gas<=0)return NAN; // nonzero K0 at zero thermal speed is inconsistent
  const double t0=gas/cv,scale=fmax(t0,td),speed=sqrt(scale)/sqrt(t0);
  if(!isfinite(t0)||!isfinite(speed))return NAN;
  const double divider=cv/speed;
  double a,b;
  if(kd<=divider){const double s=kd/divider;b=1/(1+s);a=s*b;}
  else {const double s=divider/kd;a=1/(1+s);b=s*a;}
  const double t=t0/scale,d=td/scale;
  double lo=sqrt(fmin(t,d)),hi=1,y=hi;
  bool solved=false;
  for(int it=0;it<96;++it) {
    const double yy=y*y, f=b*(yy-t)+a*y*(yy-d);
    const double norm=b*(yy+t)+a*y*(yy+d);
    if(fabs(f)<=8*DBL_EPSILON*norm){solved=true;break;}
    if(f>0)hi=y;else lo=y;
    if(hi-lo<=8*DBL_EPSILON*hi){y=lo+.5*(hi-lo);solved=true;break;}
    const double derivative=2*b*y+a*(3*yy-d);
    const double next=y-f/derivative;
    y=(derivative>0&&isfinite(next)&&next>lo&&next<hi)?next:lo+.5*(hi-lo);
  }
  if(!solved)return NAN;
  const double r=(a*y)/(b+a*y);
  // Avoid Eg-Eg' cancellation for arbitrarily weak exchange.
  return r*(gas-cv*td);
}

DUST_HD inline int dust_material_cell(const double *input,const double *a,double *out,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int i) {
  const double heating=input[i],density=input[nc+i],old=input[2*nc+i],capacity=input[3*nc+i];
  if(!isfinite(heating)||!isfinite(density)||!isfinite(old)||!isfinite(capacity)||
      heating<0||density<0||old<0||capacity<=0) return 2;
  const double *p=a+nt,*band=a+3*nt;
  // Modes 2/3 share gas ABI; 2 freezes K, 3 solves its thermal speed.
  double gas=0,gas_cv=0,ratio=0,kd=0;
  if(use_u>=2) {
    gas=input[4*size_t(nc)+i];gas_cv=input[5*size_t(nc)+i];
    const double conductance=input[6*size_t(nc)+i];kd=dt*conductance;
    if(!isfinite(gas)||gas<0||!isfinite(gas_cv)||gas_cv<=0||
        !isfinite(conductance)||conductance<0||!isfinite(kd))return 2;
    if(use_u==2)ratio=kd<=gas_cv?(kd/gas_cv)/(1+kd/gas_cv):1/(1+gas_cv/kd);
    out[size_t(ng+2)*nc+i]=0;
  }
  double *rate=out+size_t(i)*ng;
  for(int g=0;g<ng;++g)rate[g]=0;
  if(density==0) {
    if(heating!=0 || (use_u && old!=0)) return 5;
    out[size_t(ng)*nc+i]=use_u?0:old/capacity;
    out[size_t(ng+1)*nc+i]=old;
    return 0;
  }
  const double target=old+dt*heating+ratio*gas;
  const double floor=dust_material_u(a,nt,use_u,bath,density,capacity);
  const double top=dust_material_u(a,nt,use_u,exp(a[nt-1]),density,capacity);
  const double qfloor=use_u==3?dust_gas_transfer(gas,gas_cv,kd,bath):0;
  const double qtop=use_u==3?dust_gas_transfer(gas,gas_cv,kd,exp(a[nt-1])):0;
  if(!isfinite(qfloor)||!isfinite(qtop))return 6;
  if(!isfinite(target)||target<(floor+ratio*gas_cv*bath-qfloor)*(1-tolerance)||
      target>(top+dt*density*(p[nt-1]-background)+ratio*gas_cv*exp(a[nt-1])-qtop)*(1+64*DBL_EPSILON))return 5;
  double lower=0,upper=fmin(density*(p[nt-1]-background),
      fmax(heating+(old-floor)/dt+(ratio*(gas-gas_cv*bath)+qfloor)/dt,0.));
  for(int it=0;it<80;++it) {
    const double mid=lower+0.5*(upper-lower),power=background+mid/density;
    int k=0;while(k<nt-2 && power>p[k+1])++k;
    const double w=(power-p[k])/(p[k+1]-p[k]);
    const double temperature=exp(a[k]+w*(a[k+1]-a[k]));
    const double q=use_u==3?dust_gas_transfer(gas,gas_cv,kd,temperature):ratio*(gas-gas_cv*temperature);
    if(!isfinite(q))return 6;
    const double residual=(dust_material_u(a,nt,use_u,temperature,density,capacity)-old)/dt+mid-heating-q/dt;
    if(residual>0)upper=mid;else lower=mid;
  }
  const double emitted=lower+0.5*(upper-lower),increment=emitted/density;
  const double power=background+increment;
  int k=0;while(k<nt-2 && power>p[k+1])++k;
  const double w=(power-p[k])/(p[k+1]-p[k]);
  const double temperature=exp(a[k]+w*(a[k+1]-a[k]));
  out[size_t(ng)*nc+i]=temperature;
  out[size_t(ng+1)*nc+i]=dust_material_u(a,nt,use_u,temperature,density,capacity);
  if(use_u>=2) {
    const double q=use_u==3?dust_gas_transfer(gas,gas_cv,kd,temperature):ratio*(gas-gas_cv*temperature);
    if(!isfinite(q)||gas-q<0)return 5;
    out[size_t(ng+2)*nc+i]=q;
  }
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
// Conservative backward-Euler gas/dust exchange, frozen gas heat capacity
// and collision conductance. Input is cell-major [Eg,Cg,Ed,n_dust_ref,K].
// Table is [log(T), U/H]; output [Eg',Ed',Td',Qgas_to_dust]. No radiation.
DUST_HD inline int dust_exchange_cell(const double *input,const double *a,double *out,
    int nt,double dt,double floor_t,int i) {
  const double *v=input+5*size_t(i);
  const double eg=v[0],cg=v[1],ed=v[2],rho=v[3],conductance=v[4];
  for(int k=0;k<5;++k)if(!isfinite(v[k])||v[k]<0)return 2;
  if(cg<=0)return 2;
  double *o=out+4*size_t(i);
  if(rho==0){if(ed!=0)return 5;o[0]=eg;o[1]=ed;o[2]=floor_t;o[3]=0;return 0;}
  const double bottom=rho*a[nt],top=rho*a[2*nt-1];
  if(ed<bottom||ed>top||!isfinite(top))return 5;
  int old_k=0;while(old_k<nt-2&&ed/rho>a[nt+old_k+1])++old_k;
  const double old_t=exp(a[old_k]+(ed/rho-a[nt+old_k])/(a[nt+old_k+1]-a[nt+old_k]) *
       (a[old_k+1]-a[old_k]));
  if(dt==0||conductance==0){o[0]=eg;o[1]=ed;o[2]=old_t;o[3]=0;return 0;}
  const double kd=conductance*dt;
  if(!isfinite(kd))return 2;
  const double ratio=kd<=cg?(kd/cg)/(1+kd/cg):1/(1+cg/kd);
  if(!isfinite(ratio))return 2;
  double lo=log(floor_t),hi=a[nt-1],value=0,temperature=0,residual=0;
  // Residual Q - ratio*(Eg-Cg*Td), monotone in Td, bounded by supplied U(T).
  for(int it=-2;it<96;++it) {
    const double x=it==-2?lo:it==-1?hi:lo+.5*(hi-lo);
    int k=0;while(k<nt-2&&x>a[k+1])++k;
    value=rho*(a[nt+k]+(x-a[k])/(a[k+1]-a[k])*(a[nt+k+1]-a[nt+k]));
    temperature=exp(x);
    residual=(value-ed)-ratio*(eg-cg*temperature);
    if(!isfinite(value)||!isfinite(residual))return 2;
    const double tol=64*DBL_EPSILON*fmax(fmax(value,ed),ratio*fmax(eg,cg*temperature));
    if(it==-2){if(residual>tol)return 5;continue;}
    if(it==-1){if(residual < -tol)return 5;continue;}
    if(residual>0)hi=x;else lo=x;
  }
  const double q=value-ed,next_gas=eg-q;
  if(!isfinite(next_gas)||next_gas<0)return 5;
  if(fabs(residual)>256*DBL_EPSILON*fmax(fmax(value,ed),ratio*fmax(eg,cg*temperature)))return 6;
  o[0]=next_gas;o[1]=value;o[2]=temperature;o[3]=q;
  return 0;
}
#undef DUST_HD

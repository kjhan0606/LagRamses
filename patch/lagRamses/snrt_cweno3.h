// Quadratic 3-D CWENO reconstruction; mirror of mn_cweno_faces in Fortran.
// Cell averages in x-fastest [-1,1]^3 order; output -x,+x,-y,+y,-z,+z.
// No allocation, global state, floors or post-update intensity clipping.
#ifndef SNRT_CWENO3_H
#define SNRT_CWENO3_H
#include <cmath>
#include <cfloat>
#ifdef __CUDACC__
#define SNRT_CWENO_HD __host__ __device__
#else
#define SNRT_CWENO_HD
#endif
SNRT_CWENO_HD inline void snrt_cweno3_faces(const double *stencil,double *face) {
  double scale=0;
  for(int j=0;j<27;++j)scale=fmax(scale,stencil[j]);
  for(int j=0;j<6;++j)face[j]=0;
  if(scale<=0||stencil[13]<=0)return;
  double s[27],p[9][9]={},w[9],c[9]={};
  for(int j=0;j<27;++j)s[j]=stencil[j]/scale;
  const double v=s[13];
  p[0][0]=.5*(s[14]-s[12]);p[0][1]=.5*(s[16]-s[10]);p[0][2]=.5*(s[22]-s[4]);
  p[0][3]=s[14]-2*v+s[12];p[0][4]=s[16]-2*v+s[10];p[0][5]=s[22]-2*v+s[4];
  p[0][6]=.5*(s[17]-s[11]-s[15]+s[9]);
  p[0][7]=.5*(s[23]-s[5]-s[21]+s[3]);
  p[0][8]=.5*(s[25]-s[7]-s[19]+s[1]);
  int k=0;
  for(int sz=-1;sz<=1;sz+=2)for(int sy=-1;sy<=1;sy+=2)for(int sx=-1;sx<=1;sx+=2){
    ++k;p[k][0]=sx*(s[13+sx]-v);p[k][1]=sy*(s[13+3*sy]-v);p[k][2]=sz*(s[13+9*sz]-v);
  }
  double sum=0;
  for(k=0;k<9;++k){
    double linear=0,diagonal=0,cross=0;
    for(int a=0;a<3;++a){linear+=p[k][a]*p[k][a];diagonal+=p[k][a+3]*p[k][a+3];cross+=p[k][a+6]*p[k][a+6];}
    const double beta=linear+(13./3)*diagonal+(7./6)*cross;
    w[k]=(k==0?.5:1./16)/((1e-12+beta)*(1e-12+beta));sum+=w[k];
  }
  for(k=0;k<9;++k){w[k]/=sum;for(int a=0;a<9;++a)c[a]+=p[k][a]*w[k];}
  const double bx[3]={-.5,0,.5},bq[3]={1./6,-1./3,1./6};
  double lo=v,hi=v;
  for(int z=0;z<3;++z)for(int y=0;y<3;++y)for(int x=0;x<3;++x){
    const double value=v+c[0]*bx[x]+c[1]*bx[y]+c[2]*bx[z]+c[3]*bq[x]+c[4]*bq[y]+c[5]*bq[z]+
      c[6]*bx[x]*bx[y]+c[7]*bx[x]*bx[z]+c[8]*bx[y]*bx[z];
    lo=fmin(lo,value);hi=fmax(hi,value);
  }
  double theta=1;
  if(lo<0)theta=fmin(theta,v/(v-lo));
  if(hi>2*v)theta=fmin(theta,v/(hi-v));
  if(theta<1)theta*=1-16*DBL_EPSILON;
  for(int a=0;a<9;++a)c[a]*=theta;
  for(int a=0;a<3;++a){face[2*a]=scale*(v-.5*c[a]+c[a+3]/6);face[2*a+1]=scale*(v+.5*c[a]+c[a+3]/6);}
}
#undef SNRT_CWENO_HD
#endif

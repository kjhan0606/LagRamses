// The existing spectral material law, without S_N's FP32 number encoding
// or spatial transport. All angular arrays are one-cell transaction scratch.
#include "snrt_band_spectrum.h"

template<int K> static int material(int nd,int ng,double *number,double *energy,
    const double *edges,const double *columns,double *budget,double *an,double *ae,
    double *returned,snrt_band::Secondary secondary,double xi,double *deposition,
    const snrt_band::Dust *dust,bool scatter) {
  std::vector<snrt_band::Grid<K>> grid;
  for(int g=0;g<ng;++g)grid.emplace_back(edges[g],edges[g+1],g==ng-1,dust?dust->energy+g*K:nullptr);
  std::vector<double> oldn(number,number+nd*ng),olde(energy,energy+nd*ng),shift(olde);
  std::vector<double> zero(ng,0),hhe(3*ng,0),hhee(3*ng,0),dn(ng,0),de(ng,0),raw(ng,0),used(ng,0);
  double total=0;
  int rc=snrt_band::absorb_cell<K,double>(number,shift.data(),oldn.data(),olde.data(),zero.data(),columns,
      grid,budget,hhe.data(),hhee.data(),returned,raw.data(),used.data(),&total,1,1,nd,ng,0,
      secondary,&xi,deposition,dust,dn.data(),de.data());
  if(rc)return rc;
  if(dust && scatter){rc=snrt_band::scatter_cell<K,double>(number,shift.data(),zero.data(),grid,*dust,1,1,nd,ng,0);if(rc)return rc;}
  std::copy(shift.begin(),shift.end(),energy);
  for(int g=0;g<ng;++g){
    for(int s=0;s<3;++s){an[4*g+s]=hhe[s*ng+g];ae[4*g+s]=hhee[s*ng+g];}
    an[4*g+3]=dn[g];ae[4*g+3]=de[g];
  }
  return 0;
}
extern "C" int snrt_mn_band_material_c(int nd,int ng,int nb,int scatter,double *number,double *energy,
    const double *edges,const double *columns,double *budget,double *an,double *ae,double *returned,
    snrt_band::Secondary secondary,double xi,double *deposition,const double *grains,
    const double *abs,const double *transport,const double *nodes,const double *weights) {
  if(nd<1||ng!=9||(nb!=0&&nb!=4&&nb!=6))return 1;
  if(nb==0)return material<64>(nd,ng,number,energy,edges,columns,budget,an,ae,returned,secondary,xi,deposition,nullptr,false);
  snrt_band::Dust dust{grains,abs,transport,nodes,weights,nb};
  return material<128>(nd,ng,number,energy,edges,columns,budget,an,ae,returned,secondary,xi,deposition,&dust,scatter!=0);
}

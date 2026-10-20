#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

extern "C" int snrt_mn_band_material_c(int,int,int,int,double*,double*,const double*,const double*,double*,
    double*,double*,double*,int(*)(double,double,double*),double,double*,const double*,double*,const double*,
    const double*,const double*,const double*,const double*);

int main() {
  constexpr int nd=2,ng=9,nb=4,K=128;
  const double edges[ng+1]={13.6,24.59,54.42,100,200,400,800,1500,3000,10000};
  double number[nd*ng]={},energy[nd*ng]={},columns[3*nd]={1e18,0,0,1e18,0,0},
      budget[3*nd]={1e40,1e40,1e40,1e40,1e40,1e40};
  double an[4*ng]={},ae[4*ng]={},returned[ng]={},deposition[8*nd]={};
  // Arrays use group-major, then direction order, matching angular_read's
  // quadrature-weighted ray payload.  The same spectrum in each ray isolates
  // directional integration from energy-reconstruction details.
  constexpr double ray_share[nd]={0.35,0.65};
  const double ray_energy=20.0;
  for(int d=0;d<nd;++d) {
    number[d]=ray_share[d]*1e8;
    energy[d]=ray_energy*number[d];
  }
  const double direction[3*nd]={0.6,0.8,0.0,-0.8,0.6,0.0};
  double first[3*4*ng]={},grains[nb*nd]={5e-4,0,0,0};
  std::vector<double> absorption(nb*ng*K,0.0),scattering(nb*ng*K,0.0),nodes(ng*K),weights(nd,1.0);
  for(int g=0;g<ng;++g)for(int j=0;j<K;++j) {
    nodes[g*K+j]=edges[g]*std::exp(std::log(edges[g+1]/edges[g])*j/(K-1));
    if(g==0)absorption[g*K+j]=1e3;
  }
  const int status=snrt_mn_band_material_c(nd,ng,nb,0,number,energy,edges,columns,budget,an,ae,returned,
      nullptr,0.0,deposition,direction,first,grains,absorption.data(),scattering.data(),nodes.data(),weights.data());
  if(status!=0||!(ae[3]>0)||!(an[3]>0)) {
    std::fprintf(stderr,"band momentum smoke: absorber rejected or dust channel empty (status=%d, E=%.17g)\n",status,ae[3]);
    return 1;
  }
  for(int s=0;s<4;++s) {
    const double captured=ae[s];
    const double px=first[3*s],py=first[3*s+1],pz=first[3*s+2];
    const double scale=std::max(captured,1.0);
    const double expected_x=captured*(ray_share[0]*direction[0]+ray_share[1]*direction[3]);
    const double expected_y=captured*(ray_share[0]*direction[1]+ray_share[1]*direction[4]);
    if(std::abs(px-expected_x)>2e-12*scale||
       std::abs(py-expected_y)>2e-12*scale||std::abs(pz)>2e-12*scale) {
      std::fprintf(stderr,"band momentum smoke: receiver %d E=%.17g J=(%.17g,%.17g,%.17g), expected=(%.17g,%.17g,0)\n",
          s,captured,px,py,pz,expected_x,expected_y);
      return 2;
    }
  }
  std::puts("SNRT_MOMENT_BAND_MOMENTUM_OK");
  return 0;
}

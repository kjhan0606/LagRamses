#ifndef SNRT_BAND_SPECTRUM_H
#define SNRT_BAND_SPECTRUM_H
#include <array>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <vector>
#include <cstdio>

// hhe_maxent64_v1: a two-moment closure, NOT recovery of a unique SED.
// Positive trapezoidal dE prior on logarithmic nodes, including endpoints.
// H/He ground-state fits: Verner et al. 1996, ApJ 465, 487, Table 1.
// https://www.pa.uky.edu/~verner/photo.html (same fits as primordial.py).
namespace snrt_band {
constexpr int nodes=64;
constexpr int dust_nodes=128;
struct Dust {
  // columns[bin][cell] in g/cm2; immutable kappa[bin][group][node]
  // in cm2/g. Angular weights are normalized direction-integrated shares.
  const double *columns,*abs,*transport,*energy,*weights;
  int bins=4;
};
constexpr double threshold[3]={13.60,24.59,54.42};
constexpr double moment_tolerance=8*FLT_EPSILON;
// Immutable, already-loaded physical-table callback. Output order is heat,
// HI/HeI/HeII ionization, excitation; no table loading inside OpenMP.
using Secondary = int (*)(double,double,double*);
inline double sigma(int s,double e) {
  constexpr double fit[3][7]={{.4298,5.475e4,32.88,2.963,0,0,0},
      {13.61,949.2,1.469,3.188,2.039,.4434,2.136},
      {1.720,1.369e4,32.88,2.963,0,0,0}};
  if(s<0 || s>2 || e<threshold[s] || e>5e4)return 0;
  const double *f=fit[s],x=e/f[0]-f[5],y=std::hypot(x,f[6]);
  return 1e-18*f[1]*((x-1)*(x-1)+f[4]*f[4])*std::pow(y,.5*f[3]-5.5)*
      std::pow(1+std::sqrt(y/f[2]),-f[3]);
}
template<int K=nodes> struct Grid {
  std::array<double,K> e,x,w;
  double cross[3][K];
  Grid(double lo,double hi,bool final_band=false,const double *bound=nullptr) {
    for(int j=0;j<K;++j) {
      e[j]=lo*std::exp(std::log(hi/lo)*j/(K-1));
      if(j==0)e[j]=lo;
      if(j==K-1)e[j]=hi;
      if(bound)e[j]=bound[j]; // caller has checked exact source grid support
      x[j]=(e[j]-lo)/(hi-lo);
      // Upper endpoint is the band's left limit: a threshold does not
      // acquire opacity in the band immediately BELOW that threshold.
      const double eval=(j==K-1 && !final_band)?std::nextafter(hi,lo):e[j];
      for(int s=0;s<3;++s)cross[s][j]=sigma(s,eval);
    }
    for(int j=0;j<K;++j)w[j]=(e[std::min(j+1,K-1)]-e[std::max(j-1,0)])/(2*(hi-lo));
  }
  bool reconstruct(double mean,std::array<double,K>&p) const {
    if(!std::isfinite(mean) || mean<e[0]*(1-moment_tolerance) ||
        mean>e[K-1]*(1+moment_tolerance))return false;
    p.fill(0);
    if(mean<=e[0]){p[0]=1;return true;}
    if(mean>=e[K-1]){p[K-1]=1;return true;}
    const double target=(mean-e[0])/(e[K-1]-e[0]);
    auto evaluate=[&](double b) {
      double norm=0,m=0;
      for(int j=0;j<K;++j){p[j]=w[j]*std::exp(b*(x[j]-(b>0?1:0)));norm+=p[j];m+=p[j]*x[j];}
      for(double &v:p)v/=norm;
      return m/norm;
    };
    double left=-1,right=1;
    while(evaluate(left)>target && left>-1e12)left*=2;
    while(evaluate(right)<target && right<1e12)right*=2;
    for(int it=0;it<80;++it) {
      const double b=.5*(left+right),m=evaluate(b);
      if(std::abs(m-target)<2e-15)return true;
      if(m<target)left=b;else right=b;
    }
    return std::abs(evaluate(.5*(left+right))-target)<2e-13;
  }
};

// Cell-local finite-inventory operator, after the existing upwind stencil.
// Band order matches the fixed operator. Each species cap is constant across
// directions and nodes within a band; rejected packets retain E/N. A depleted
// species must not veto absorption by another species with remaining atoms.
// Every pointer here addresses transaction scratch, never published state.
template<int K=nodes> inline int absorb_cell(float *next,double *shift,const float *number,const double *energy,
    const double *reference,const double *columns,const std::vector<Grid<K>>&grids,
    float *budget,float *hhe,double *hhe_e,float *returned,float *raw,float *absorbed,float *total,
    int no,int nw,int nd,int ng,int cell,Secondary secondary=nullptr,
    const double *xi=nullptr,double *deposition=nullptr,const Dust *grains=nullptr,
    float *dust=nullptr,double *dust_e=nullptr) {
  struct Ray {double n[4]={},e[4]={},survive=0,survive_e=0;};
  std::vector<Ray> rays(nd);
  const size_t groups=size_t(no)*ng;
  bool occupied=false;
  if(grains)for(int b=0;b<grains->bins;++b)occupied=occupied || grains->columns[b*no+cell]>0;
  const bool iron=grains && grains->bins==6 &&
      (grains->columns[4*no+cell]>0 || grains->columns[5*no+cell]>0);
  const int channels=occupied?4:3;
  for(int g=0;g<ng;++g) {
    const auto &grid=grids[g];
    const size_t out=size_t(g)*no+cell;
    if(!occupied && grid.e.back()<=threshold[0])continue; // exact zero-opacity transport
    if(!occupied && columns[cell]==0 && columns[no+cell]==0 && columns[2*no+cell]==0)continue;
    double fraction[4][K],transmit[K];
    for(int j=0;j<K;++j) {
      double tau=0;
      for(int s=0;s<3;++s){fraction[s][j]=columns[s*no+cell]*grid.cross[s][j];tau+=fraction[s][j];}
      if(occupied){
        fraction[3][j]=0;
        for(int b=0;b<grains->bins;++b)fraction[3][j]+=grains->columns[b*no+cell]*grains->abs[(b*ng+g)*K+j];
        tau+=fraction[3][j];
      }
      if(!std::isfinite(tau))return 3;
      const double loss=-std::expm1(-tau);
      // Use complementary fractions from the SAME rounded absorption.
      // Once loss rounds to one, a separate exp(-tau) would create an
      // energy-only extinction tail even though the ledger absorbed it all.
      transmit[j]=1-loss;
      for(int s=0;s<channels;++s)fraction[s][j]=tau>0?loss*fraction[s][j]/tau:0;
    }
    double wanted[4]={},wanted_e[4]={},node_count[3][K]={};
    for(int d=0;d<nd;++d) {
      const size_t k=(size_t(g)*nd+d)*nw+cell;
      Ray &r=rays[d];r=Ray{};
      if(number[k]==0)continue;
      std::array<double,K> p;
      const double mean=energy[k]/number[k];
      if(!grid.reconstruct(mean,p)){
        std::fprintf(stderr,"SNRT band inadmissible transported mean cell=%d group=%d N=%.17g E=%.17g mean=%.17g\n",cell,g, double(number[k]),energy[k],mean);
        return 3;
      }
      // Endpoint projection is only a FP32 moment-roundoff operation. Scale
      // represented number (within that tolerance) to preserve actual E.
      const double bounded=std::max(grid.e[0],std::min(grid.e.back(),mean));
      const double count=energy[k]/bounded;
      for(int j=0;j<K;++j) {
        const double q=count*p[j];
        // Preserve the live Fe comparison's <=4 eV absorption restriction.
        // A low GROUP MEAN is insufficient: do not discard a hard tail or
        // assign missing photoelectron escape energy to grain heat.
        if(iron && grid.e[j]>4.0 && q>0 && fraction[3][j]>0)return 10;
        r.survive+=q*transmit[j];r.survive_e+=q*transmit[j]*grid.e[j];
        for(int s=0;s<channels;++s){
          const double captured=q*fraction[s][j];
          r.n[s]+=captured;r.e[s]+=captured*grid.e[j];
          if(secondary && s<3)node_count[s][j]+=captured;
        }
      }
      for(int s=0;s<channels;++s){wanted[s]+=r.n[s];wanted_e[s]+=r.e[s];}
    }
    double cap[3]={1,1,1};
    for(int s=0;s<3;++s)if(wanted[s]>0) {
      // Leave only rounding headroom, not a physical residual reservoir.
      const double available=budget[s*no+cell];
      if(wanted[s]>=available*(1-4*FLT_EPSILON))
        cap[s]=std::min(1.,available*(1-4*FLT_EPSILON)/wanted[s]);
    }
    double requested=0,accepted=0;
    for(int s=0;s<3;++s) {
      const float used=float(cap[s]*wanted[s]);
      hhe[s*groups+out]=used;hhe_e[s*groups+out]=cap[s]*wanted_e[s];
      budget[s*no+cell]-=used;
      if(budget[s*no+cell]<0)return 3;
      requested+=wanted[s];accepted+=cap[s]*wanted[s];
      if(secondary) {
        deposition[s*no+cell]+=cap[s]*wanted[s];
        for(int j=0;j<K;++j) {
          const double excess=std::max(0.,grid.e[j]-threshold[s]);
          const double electron_energy=cap[s]*node_count[s][j]*excess;
          if(electron_energy==0)continue;
          double f[5],norm=0;
          if(secondary(excess,xi[cell],f)!=0)return 9;
          for(double v:f){if(!std::isfinite(v)||v<0)return 9;norm+=v;}
          if(std::abs(norm-1)>2e-12)return 9;
          for(int k=0;k<5;++k)deposition[(3+k)*no+cell]+=electron_energy*f[k];
        }
      }
    }
    // Grain captures are not atom-limited. Rejected H/He packets are
    // returned with their energy, not passed through a second dust sink.
    if(occupied){
      dust[out]=float(wanted[3]);dust_e[out]=wanted_e[3];
      requested+=wanted[3];accepted+=wanted[3];
    }
    raw[out]=float(requested);returned[out]=float(requested-accepted);absorbed[out]=float(accepted);
    total[cell]+=absorbed[out];
    for(int d=0;d<nd;++d) {
      const size_t k=(size_t(g)*nd+d)*nw+cell;
      const Ray &r=rays[d];double n=r.survive,e=r.survive_e;
      for(int s=0;s<3;++s){n+=(1-cap[s])*r.n[s];e+=(1-cap[s])*r.e[s];}
      // Stable positive survivor sums avoid subtracting nearly equal E's
      // in optically thick cells. Moment solve conserves E to FP64 tolerance.
      next[k]=float(n);shift[k]=e-reference[g]*next[k];
      if(!std::isfinite(n)||!std::isfinite(e)||n<0||e<0 || (next[k]==0 && e!=0)){
        std::fprintf(stderr,"SNRT band unrepresentable survivor cell=%d group=%d N=%.17g E=%.17g\n",cell,g,n,e);
        return 3;
      }
    }
  }
  return 0;
}

template<int K> inline int scatter_cell(float *number,double *shift,const double *reference,
    const std::vector<Grid<K>>&grids,const Dust &grains,int no,int nw,int nd,int ng,int cell) {
  bool occupied=false;
  for(int b=0;b<grains.bins;++b)occupied=occupied || grains.columns[b*no+cell]>0;
  if(!occupied)return 0;
  std::vector<double> keep_n(nd),keep_e(nd);
  for(int g=0;g<ng;++g){
    const auto &grid=grids[g];
    double fraction[K],sum_n=0,sum_e=0;
    bool scatters=false;
    for(int j=0;j<K;++j){
      double tau=0;
      for(int b=0;b<grains.bins;++b)tau+=grains.columns[b*no+cell]*grains.transport[(b*ng+g)*K+j];
      if(!std::isfinite(tau))return 3;
      fraction[j]=-std::expm1(-tau);scatters=scatters || tau>0;
    }
    if(!scatters)continue;
    std::fill(keep_n.begin(),keep_n.end(),0);std::fill(keep_e.begin(),keep_e.end(),0);
    for(int d=0;d<nd;++d){
      const size_t k=(size_t(g)*nd+d)*nw+cell;
      if(number[k]==0)continue;
      const double e=reference[g]*number[k]+shift[k],mean=e/number[k];
      std::array<double,K> p;
      if(!grid.reconstruct(mean,p))return 3;
      const double count=e/std::max(grid.e[0],std::min(grid.e.back(),mean));
      for(int j=0;j<K;++j){
        const double scattered=count*p[j]*fraction[j],kept=count*p[j]*(1-fraction[j]);
        keep_n[d]+=kept;keep_e[d]+=kept*grid.e[j];
        sum_n+=scattered;sum_e+=scattered*grid.e[j];
      }
    }
    // One elastic angular relaxation of the surviving reconstructed
    // spectrum. No material force or heating in this stationary comparison.
    for(int d=0;d<nd;++d){
      const size_t k=(size_t(g)*nd+d)*nw+cell;
      const double n=keep_n[d]+grains.weights[d]*sum_n,e=keep_e[d]+grains.weights[d]*sum_e;
      number[k]=float(n);shift[k]=e-reference[g]*number[k];
      if(!std::isfinite(n)||!std::isfinite(e)||n<0||e<0 || (number[k]==0 && e!=0))return 3;
    }
  }
  return 0;
}
} // namespace snrt_band
#endif

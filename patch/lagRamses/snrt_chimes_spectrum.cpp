#include "snrt_chimes_spectrum.h"
#include "snrt_chimes_spectrum_internal.h"
#include <hdf5.h>
#include <openssl/evp.h>
#include <fstream>
#include <memory>
#include <string>
#include <limits>
#include <stdexcept>
#include <cstdint>

namespace {
constexpr int ng=9,K=128,nr=311;
constexpr double edges[10]={.01,1,5.6,11.2,13.6,24.59,54.42,500,2000,10000};
using snrt_chimes_detail::Bank;
struct H5Object {
  hid_t id; herr_t (*close)(hid_t);
  ~H5Object(){if(id>=0)close(id);}
};
void require(bool valid){if(!valid)throw std::runtime_error("invalid atomic spectral bank");}
template<class T> std::vector<T> read(hid_t file,const char *name,hid_t type,std::vector<hsize_t> shape) {
  H5Object data{H5Dopen2(file,name,H5P_DEFAULT),H5Dclose};require(data.id>=0);
  H5Object space{H5Dget_space(data.id),H5Sclose};require(space.id>=0);
  require(H5Sget_simple_extent_ndims(space.id)==int(shape.size()));
  std::vector<hsize_t> dims(shape.size());
  require(H5Sget_simple_extent_dims(space.id,dims.data(),nullptr)>=0 && dims==shape);
  size_t count=1;for(auto d:dims){require(d>0 && d<=10000000/count);count*=d;}
  std::vector<T> v(count);require(H5Dread(data.id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,v.data())>=0);return v;
}
std::string attribute(hid_t file,const char *name) {
  H5Object attr{H5Aopen(file,name,H5P_DEFAULT),H5Aclose};require(attr.id>=0);
  H5Object type{H5Aget_type(attr.id),H5Tclose};require(type.id>=0 && H5Tget_class(type.id)==H5T_STRING);
  H5Object space{H5Aget_space(attr.id),H5Sclose};require(space.id>=0 && H5Sget_simple_extent_npoints(space.id)==1);
  if(H5Tis_variable_str(type.id)>0){
    char *value=nullptr;require(H5Aread(attr.id,type.id,&value)>=0 && value);
    std::string result(value);H5free_memory(value);return result;
  }
  const size_t n=H5Tget_size(type.id);require(n>0 && n<1024);
  std::vector<char> value(n+1,0);require(H5Aread(attr.id,type.id,value.data())>=0);return value.data();
}
std::array<unsigned char,32> digest(const char *path) {
  std::ifstream file(path,std::ios::binary);require(bool(file));
  std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> c(EVP_MD_CTX_new(),EVP_MD_CTX_free);
  require(bool(c) && EVP_DigestInit_ex(c.get(),EVP_sha256(),nullptr)==1);
  std::array<char,65536> buf;
  while(file){file.read(buf.data(),buf.size());require(EVP_DigestUpdate(c.get(),buf.data(),file.gcount())==1);}
  require(file.eof());std::array<unsigned char,32> out;unsigned len=0;
  require(EVP_DigestFinal_ex(c.get(),out.data(),&len)==1 && len==32);return out;
}
void verify_reaction_order(const std::vector<int> &map) {
  // Exact ordered mapping extracted from the pinned main table. Shape and
  // index-range checks alone do not detect a valid-species miswire.
  std::vector<unsigned char> bytes(map.size()*4);
  for(size_t j=0;j<map.size();++j)for(int k=0;k<4;++k)
    bytes[4*j+k]=(uint32_t(map[j])>>(8*k))&255;
  unsigned char hash[32];unsigned n=0;
  require(EVP_Digest(bytes.data(),bytes.size(),hash,&n,EVP_sha256(),nullptr)==1 && n==32);
  const char *hex="0123456789abcdef";std::string encoded;
  for(auto v:hash){encoded+=hex[v>>4];encoded+=hex[v&15];}
  require(encoded=="8cf424068bd05d205155e9359d07a84c8abbfb5f880b268cf33becf94148d351");
}
}

extern "C" int snrt_chimes_band_load(const char *path,void **handle,int *reactions,int *shells,double *identity) {
  if(!path || !handle || *handle || !reactions || !shells || !identity)return 1;
  try {
    const auto before=digest(path);
    H5Object file{H5Fopen(path,H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose};require(file.id>=0);
    require(attribute(file.id,"schema")=="snrt_chimes_atomic_shell_nodes_v1");
    require(attribute(file.id,"main_sha256")=="8bde78faacd59249fad5e810cae43311ed03ef09131c62b0a0a5c39e8407cb0a");
    auto b=std::make_unique<Bank>();
    const auto e=read<double>(file.id,"edges_ev",H5T_NATIVE_DOUBLE,{10});
    require(std::equal(e.begin(),e.end(),edges));
    b->energy=read<double>(file.id,"node_ev",H5T_NATIVE_DOUBLE,{ng,K});
    b->reaction=read<int>(file.id,"reactions",H5T_NATIVE_INT,{nr,5});
    verify_reaction_order(b->reaction);
    H5Object data{H5Dopen2(file.id,"binding_ev",H5P_DEFAULT),H5Dclose};require(data.id>=0);
    H5Object space{H5Dget_space(data.id),H5Sclose};require(space.id>=0);
    const auto ns=H5Sget_simple_extent_npoints(space.id);require(ns>0 && ns<=4096);
    b->binding=read<double>(file.id,"binding_ev",H5T_NATIVE_DOUBLE,{hsize_t(ns)});
    b->shell=read<int>(file.id,"shell_reaction",H5T_NATIVE_INT,{hsize_t(ns)});
    b->sigma=read<double>(file.id,"sigma_cm2",H5T_NATIVE_DOUBLE,{hsize_t(ns),ng,K});
    const int starts[5]={0,8,124,159,311};
    for(int category=0;category<4;++category)for(int r=starts[category];r<starts[category+1];++r){
      const int *m=b->reaction.data()+r*5;
      require(m[0]==category && m[1]==r-starts[category] && m[2]>0 && m[2]<157 && m[3]>0 && m[3]<157);
      require(m[4]>=1 && m[4]<=10 && (category<2 ? m[4]==1 : m[3]-m[2]==m[4]));
    }
    for(int g=0;g<ng;++g){
      const snrt_band::Grid<K> reference(edges[g],edges[g+1],g==ng-1);
      for(int j=0;j<K;++j){const double x=b->energy[g*K+j];
        require(std::isfinite(x) && std::abs(x/reference.e[j]-1)<2e-13);
      }
      require(b->energy[g*K]==edges[g] && b->energy[(g+1)*K-1]==edges[g+1]);
      b->grids.emplace_back(edges[g],edges[g+1],g==ng-1,b->energy.data()+g*K);
    }
    std::vector<int> present(nr,0);
    for(size_t s=0;s<size_t(ns);++s){
      require(b->shell[s]>=0 && b->shell[s]<nr);present[b->shell[s]]++;
      require(std::isfinite(b->binding[s]) && b->binding[s]>0);
      for(int g=0;g<ng;++g)for(int j=0;j<K;++j){
        const double v=b->sigma[(s*ng+g)*K+j],ev=b->energy[g*K+j];
        const double eval=j==K-1 && g<ng-1 ? std::nextafter(ev,edges[g]) : ev;
        require(std::isfinite(v) && v>=0 && (eval>=b->binding[s] || v==0));
      }
    }
    for(int p:present)require(p>0);
    require(before==digest(path)); // no torn replacement during admission
    for(int i=0;i<32;++i)identity[i]=before[i];
    *reactions=nr;*shells=int(ns);*handle=b.release();return 0;
  } catch(...) {return 2;}
}
extern "C" void snrt_chimes_band_free(void *handle){delete static_cast<Bank*>(handle);}
extern "C" int snrt_chimes_molecular_load(const char *path,void **handle,double *identity) {
  if(!path || !handle || *handle || !identity)return 1;
  try {
    const auto before=digest(path);
    H5Object file{H5Fopen(path,H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose};require(file.id>=0);
    require(attribute(file.id,"schema")=="snrt_chimes_molecular_nodes_v1");
    require(attribute(file.id,"main_sha256")=="8bde78faacd59249fad5e810cae43311ed03ef09131c62b0a0a5c39e8407cb0a");
    const auto e=read<double>(file.id,"edges_ev",H5T_NATIVE_DOUBLE,{10});
    require(std::equal(e.begin(),e.end(),edges));
    const auto nodes=read<double>(file.id,"node_ev",H5T_NATIVE_DOUBLE,{ng,K});
    for(int g=0;g<ng;++g){
      const snrt_band::Grid<K> grid(edges[g],edges[g+1],g==ng-1);
      for(int j=0;j<K;++j)require(std::isfinite(nodes[g*K+j]) && std::abs(nodes[g*K+j]/grid.e[j]-1)<2e-13);
    }
    auto b=std::make_unique<snrt_chimes_detail::MoleculeBank>();
    b->reaction=read<int>(file.id,"reactions",H5T_NATIVE_INT,{32,3});
    b->absorption=read<double>(file.id,"absorption_cm2",H5T_NATIVE_DOUBLE,{32,ng,K});
    b->dissociation=read<double>(file.id,"dissociation_cm2",H5T_NATIVE_DOUBLE,{32,ng,K});
    // Exact pinned CHIMES empirical channels followed by H2 and CO.
    constexpr int mapping[96]={14,7,0,32,23,0,138,1,2,141,152,0,141,140,1,151,23,2,
      140,151,0,140,23,1,142,7,7,147,149,137,147,150,1,150,149,1,146,150,0,146,145,1,
      149,7,2,145,149,0,145,7,1,139,138,1,139,137,2,144,148,2,143,156,0,143,23,23,
      152,138,23,152,2,140,152,24,137,152,151,1,153,2,141,153,138,140,153,152,1,153,151,137,
      137,1,1,148,7,23};
    require(std::equal(b->reaction.begin(),b->reaction.end(),mapping));
    for(size_t i=0;i<b->absorption.size();++i){
      const double a=b->absorption[i],d=b->dissociation[i];
      require(std::isfinite(a) && std::isfinite(d) && d>=0 && a>=d);
    }
    require(before==digest(path));
    for(int i=0;i<32;++i)identity[i]=before[i];
    *handle=b.release();return 0;
  }catch(...){return 2;}
}
extern "C" void snrt_chimes_molecular_free(void *handle){delete static_cast<snrt_chimes_detail::MoleculeBank*>(handle);}
extern "C" int snrt_chimes_band_reactions(void *handle,int reactions,int *mapping){
  if(!handle || reactions!=nr || !mapping)return 1;
  const auto &b=*static_cast<const Bank*>(handle);std::copy(b.reaction.begin(),b.reaction.end(),mapping);return 0;
}
extern "C" int snrt_chimes_band_nodes(void *handle,int nd,const double *number,const double *energy,
    double *node_number,double *node_energy) {
  // Private moving-grain staging: reconstruct every ray using the SAME bank
  // as photo chemistry. Output layout is Fortran (nd,128,9); publish atomically.
  if(!handle || nd<1 || nd>720 || !number || !energy || !node_number || !node_energy)return 1;
  const auto &b=*static_cast<const Bank*>(handle);
  try {
    bool dark=true;
    for(int ray=0;ray<ng*nd;++ray){
      const double n=number[ray],e=energy[ray];
      if(!std::isfinite(n)||!std::isfinite(e)||n<0||e<0||(n==0&&e!=0))return 2;
      if(n!=0)dark=false;
    }
    if(dark){
      // No failing work remains: publish zeros without private node copies.
      std::fill_n(node_number,size_t(nd)*K*ng,0.);
      std::fill_n(node_energy,size_t(nd)*K*ng,0.);
      return 0;
    }
    std::vector<double> nn(size_t(nd)*K*ng,0),ee(nn.size(),0);
    for(int g=0;g<ng;++g)for(int d=0;d<nd;++d){
      const int ray=g*nd+d;const double n=number[ray],e=energy[ray];
      if(!std::isfinite(n)||!std::isfinite(e)||n<0||e<0||(n==0&&e!=0))return 2;
      if(n==0)continue;
      std::array<double,K> w;const double mean=e/n;
      if(mean<b.grids[g].e.front()*(1-2e-13)||mean>b.grids[g].e.back()*(1+2e-13)||
          !b.grids[g].reconstruct(mean,w))return 2;
      double sum_n=0,sum_e=0;
      for(int k=0;k<K;++k){sum_n+=w[k];sum_e+=w[k]*b.grids[g].e[k];}
      if(!(sum_n>0&&sum_e>0))return 2;
      for(int k=0;k<K;++k){
        const size_t j=d+size_t(nd)*(k+K*g);
        nn[j]=n*(w[k]/sum_n);ee[j]=e*(w[k]*b.grids[g].e[k]/sum_e);
      }
    }
    std::copy(nn.begin(),nn.end(),node_number);std::copy(ee.begin(),ee.end(),node_energy);return 0;
  }catch(...){return 3;}
}

extern "C" int snrt_chimes_band_moments(void *handle,int directions,int reactions,
    const double *number,const double *energy,double *moments) {
  if(!handle || directions<=0 || directions>100000 || reactions!=nr || !number || !energy || !moments)return 1;
  const auto &b=*static_cast<const Bank*>(handle);
  try {
    std::vector<double> q(ng*K,0),out(ng*3*nr,0);
    // Reconstruct EACH ray; averaging E/N before this nonlinear closure
    // would erase the opacity difference between independently hardened rays.
    for(int g=0;g<ng;++g)for(int d=0;d<directions;++d){
      const size_t k=size_t(g)*directions+d;const double n=number[k],e=energy[k];
      if(!std::isfinite(n) || !std::isfinite(e) || n<0 || e<0 || (n==0 && e!=0))return 2;
      if(n==0)continue;
      const double mean=e/n;std::array<double,K> p;
      if(!b.grids[g].reconstruct(mean,p))return 2;
      const double count=e/std::max(edges[g],std::min(edges[g+1],mean));
      for(int j=0;j<K;++j)q[g*K+j]+=count*p[j];
    }
    for(int g=0;g<ng;++g){
      if(std::none_of(q.begin()+g*K,q.begin()+(g+1)*K,[](double v){return v>0;}))continue;
      for(size_t s=0;s<b.shell.size();++s)for(int j=0;j<K;++j){
        const double captures=q[g*K+j]*b.sigma[(s*ng+g)*K+j];
        const double ev=b.energy[g*K+j];const int r=b.shell[s];
        out[(g*3)*nr+r]+=captures;
        out[(g*3+1)*nr+r]+=captures*ev;
        out[(g*3+2)*nr+r]+=captures*std::max(0.,ev-b.binding[s]);
      }
    }
    for(double v:out)if(!std::isfinite(v) || v<0)return 3;
    std::copy(out.begin(),out.end(),moments);return 0;
  } catch(...) {return 3;}
}

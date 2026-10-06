/* Dark net cooling for N_spectra==0. Operation order follows
 * calculate_total_cooling_rate / update_cooling_rates (Richings 2020).
 * Metal terms are computed in parallel and summed in table order.
 * H2O vibrational tables are indexed with the H2O temperature bin, matching
 * the host, even though those arrays are allocated on the molecular axis. */
__device__ double dmaxd(double a,double b){return a>b?a:b;}
__device__ double at2(const double *t,int x,int y,int Ny){return t[x*Ny+y];}
__device__ double at3(const double *t,int x,int y,int z,int Ny,int Nz){return t[(x*Ny+y)*Nz+z];}
__device__ double at5(const double *t,int x,int y,int z,int v,int w,int Ny,int Nz,int Nv,int Nw){
 return t[(((x*Ny+y)*Nz+z)*Nv+v)*Nw+w];
}
__device__ double lin_fix(const double *t,int x,int y,double dy,int Ny){
 double my=1.0-dy;return my*at2(t,x,y,Ny)+dy*at2(t,x,y+1,Ny);
}
__device__ double bilin(const double *t,int x,int y,double dx,double dy,int Ny){
 double mx=1.0-dx,my=1.0-dy;
 return mx*my*at2(t,x,y,Ny)+mx*dy*at2(t,x,y+1,Ny)+dx*my*at2(t,x+1,y,Ny)+dx*dy*at2(t,x+1,y+1,Ny);
}
__device__ double bilin_fix(const double *t,int x,int y,int z,double dy,double dz,int Ny,int Nz){
 double my=1.0-dy,mz=1.0-dz;
 return my*mz*at3(t,x,y,z,Ny,Nz)+my*dz*at3(t,x,y,z+1,Ny,Nz)
   +dy*mz*at3(t,x,y+1,z,Ny,Nz)+dy*dz*at3(t,x,y+1,z+1,Ny,Nz);
}
__device__ double quad_fix(const double *t,int x,int y,int z,int v,int w,
    double dy,double dz,double dv,double dw,int Ny,int Nz,int Nv,int Nw){
 const double my=1.0-dy,mz=1.0-dz,mv=1.0-dv,mw=1.0-dw;
 double o=my*mz*mv*mw*at5(t,x,y,z,v,w,Ny,Nz,Nv,Nw);
 o+=my*mz*dv*mw*at5(t,x,y,z,v+1,w,Ny,Nz,Nv,Nw);
 o+=my*dz*mv*mw*at5(t,x,y,z+1,v,w,Ny,Nz,Nv,Nw);
 o+=dy*mz*mv*mw*at5(t,x,y+1,z,v,w,Ny,Nz,Nv,Nw);
 o+=dy*dz*mv*mw*at5(t,x,y+1,z+1,v,w,Ny,Nz,Nv,Nw);
 o+=dy*mz*dv*mw*at5(t,x,y+1,z,v+1,w,Ny,Nz,Nv,Nw);
 o+=my*dz*dv*mw*at5(t,x,y,z+1,v+1,w,Ny,Nz,Nv,Nw);
 o+=dy*dz*dv*mw*at5(t,x,y+1,z+1,v+1,w,Ny,Nz,Nv,Nw);
 o+=my*mz*mv*dw*at5(t,x,y,z,v,w+1,Ny,Nz,Nv,Nw);
 o+=my*mz*dv*dw*at5(t,x,y,z,v+1,w+1,Ny,Nz,Nv,Nw);
 o+=my*dz*mv*dw*at5(t,x,y,z+1,v,w+1,Ny,Nz,Nv,Nw);
 o+=dy*mz*mv*dw*at5(t,x,y+1,z,v,w+1,Ny,Nz,Nv,Nw);
 o+=dy*dz*mv*dw*at5(t,x,y+1,z+1,v,w+1,Ny,Nz,Nv,Nw);
 o+=dy*mz*dv*dw*at5(t,x,y+1,z,v+1,w+1,Ny,Nz,Nv,Nw);
 o+=my*dz*dv*dw*at5(t,x,y,z+1,v+1,w+1,Ny,Nz,Nv,Nw);
 o+=dy*dz*dv*dw*at5(t,x,y+1,z+1,v+1,w+1,Ny,Nz,Nv,Nw);
 return o;
}
__device__ double mean_mu(const double *x,const double *elem){
 double den=0;for(int i=0;i<ns;i++)den+=x[i];
 return (1.0+elem[0]*4.0f+elem[1]*12.0f+elem[2]*14.0f+elem[3]*16.0f+elem[4]*20.0f
   +elem[5]*24.0f+elem[6]*28.0f+elem[7]*32.0f+elem[8]*40.0f+elem[9]*56.0f)/den;
}
__device__ double thermal_dv(double T,double mu){
 return sqrt((3.0f*BOLTZMANNCGS)*T/(PROTON_MASS*mu));
}
__device__ double rot_factor(const double *L0t,const double *Llte,const double *nhalf,const double *a,
    int iT,double dT,int iN,double dN,int nN,double neff){
 double L0=powten(interp(L0t,iT,dT));
 double Ll=powten(bilin(Llte,iT,iN,dT,dN,nN));
 double nhf=powten(bilin(nhalf,iT,iN,dT,dN,nN));
 double aa=powten(bilin(a,iT,iN,dT,dN,nN));
 return 1.0/((1.0f/L0)+(neff/Ll)+(1.0f/L0)*pow(neff/nhf,aa)*(1.0f-(nhf*L0/Ll)));
}
__device__ double oh_rot(const double *x,const double *elem,double T,double nh,double column,double extinction){
 double dv=thermal_dv(T,mean_mu(x,elem));
 double N_tau=1.485e11f*dv;
 double tau_T=4.0f*column/(10.0f*(T/27.0f)*6.8e-4f*N_tau);
 double c_tau=tau_T*sqrt(2.0f*PI*log(2.13f+((tau_T/EULERS_CONST)*(tau_T/EULERS_CONST))))
   /((exp(-extinction)/(1.0f+(extinction*extinction)))
     +2.0f*extinction*sqrt(log(1.0f+(tau_T/EULERS_CONST)))*sqrt(log(tau_T/(extinction*EULERS_CONST))));
 if(isnan(c_tau))c_tau=0.0f;
 double n_cr=1.5e10f*sqrt(T/1.0e3f);
 double ym=log(1.0f+(c_tau/(1.0f+10.0f*(n_cr/nh))));
 return x[sp_OH]*(2.0f*BOLTZMANNCGS*T*T*2.3e-2f/(nh*27.0f))
   *((2.0f+ym+0.6f*ym*ym)/(1.0f+c_tau+(n_cr/nh)+1.5f*sqrt(n_cr/nh)));
}

__global__ void dark_cooling_metal(Tables t,CoolTables c,const double *v,const int *entry,const Input *in,Output *out,CompactOutput *packed){
 const int tid=threadIdx.x;
 const Input &x=in[blockIdx.x];
 __shared__ double term[128],cool,logT;
 __shared__ int low2,low4,iT,iTmol,iT2,iNe2,iThi2,iT4,inHI,iNe4,inHII,iThi4;
 __shared__ double dT,dTmol,dT2,dNe2,dThi2,dT4,dnHI,dNe4,dnHII,dThi4;
 if(tid==0)cool=0;
 __syncthreads();
 if(!c.enabled || blockDim.x!=128 || (x.hybrid_mode!=0 && x.hybrid_mode!=1)){
   if(tid==0){
     out[blockIdx.x].net_cooling=0;out[blockIdx.x].cooling_valid=0;
     if(packed){packed[blockIdx.x].net_cooling=0;packed[blockIdx.x].cooling_valid=0;}
   }
   return;
 }
 if(tid==0){
   logT=log10(x.t);
   index(v+t.offset[AXT],t.nt,logT,iT,dT);
   index(v+c.o_tmol,c.nTmol,logT,iTmol,dTmol);
 }
 __syncthreads();
 for(int base=0;base<c.n1;base+=128){
   int i=base+tid;double z=0;
   if(i<c.n1)z=powten(lin_fix(v+c.o_rates,i,iT,dT,t.nt))*x.x[entry[c.base1+i]]*x.x[sp_elec];
   term[tid]=z;__syncthreads();
   if(tid==0)for(int j=0;j<128 && base+j<c.n1;j++)cool+=term[j];
   __syncthreads();
 }
 if(c.n2){
   if(tid==0){
     low2=logT<v[c.o_T2+c.nT2-1];
     if(low2){
       index(v+c.o_T2,c.nT2,logT,iT2,dT2);
       index(v+c.o_ne2,c.nNe2,log10(dmaxd(x.x[sp_elec]*x.nh,DBL_MIN)),iNe2,dNe2);
     }else index(v+c.o_Thi2,c.nThi2,logT,iThi2,dThi2);
   }
   __syncthreads();
   for(int base=0;base<c.n2;base+=128){
     int i=base+tid;double z=0;
     if(i<c.n2){
       double rate=low2?powten(bilin_fix(v+c.o_rates2,i,iT2,iNe2,dT2,dNe2,c.nT2,c.nNe2))
                       :powten(lin_fix(v+c.o_hi2,i,iThi2,dThi2,c.nThi2));
       z=rate*x.x[entry[c.base2+i]]*x.x[sp_elec];
     }
     term[tid]=z;__syncthreads();
     if(tid==0)for(int j=0;j<128 && base+j<c.n2;j++)cool+=term[j];
     __syncthreads();
   }
 }
 if(c.n4){
   if(tid==0){
     low4=logT<v[c.o_T4+c.nT4-1];
     if(low4){
       index(v+c.o_T4,c.nT4,logT,iT4,dT4);
       index(v+c.o_nHI,c.nHI4,log10(dmaxd(x.x[sp_HI]*x.nh,DBL_MIN)),inHI,dnHI);
       index(v+c.o_ne4,c.nNe4,log10(dmaxd(x.x[sp_elec]*x.nh,DBL_MIN)),iNe4,dNe4);
       index(v+c.o_nHII,c.nHII4,log10(dmaxd(x.x[sp_HII]*x.nh,DBL_MIN)),inHII,dnHII);
     }else index(v+c.o_Thi4,c.nThi4,logT,iThi4,dThi4);
   }
   __syncthreads();
   for(int base=0;base<c.n4;base+=128){
     int i=base+tid;double z=0;
     if(i<c.n4){
       double rate;
       if(low4)rate=powten(quad_fix(v+c.o_rates4,i,iT4,inHI,iNe4,inHII,dT4,dnHI,dNe4,dnHII,c.nT4,c.nHI4,c.nNe4,c.nHII4));
       else rate=powten(lin_fix(v+c.o_hi4,i,iThi4,dThi4,c.nThi4))*(x.x[sp_elec]*x.nh);
       z=rate*x.x[entry[c.base4+i]]/x.nh;
     }
     term[tid]=z;__syncthreads();
     if(tid==0)for(int j=0;j<128 && base+j<c.n4;j++)cool+=term[j];
     __syncthreads();
   }
 }
 if(tid==0){out[blockIdx.x].net_cooling=cool;out[blockIdx.x].cooling_valid=2;}
}

/* One divide by nH, then a left-to-right product sum. A per-term divide
 * matches the host expression but takes about 100 us on A10; this form
 * stays within a couple of ulps on the heating sum. */
__global__ void dark_cooling_cr(int base,int ncr0,int ncr1,const Input *in,const Output *out,double *heat){
 const int cell=blockIdx.x;
 if(out[cell].cooling_valid!=2){heat[cell]=0;return;}
 const int ncr=in[cell].mol?ncr1:ncr0;
 const double scale=3.2e-11f/in[cell].nh;
 const double *rate=out[cell].rate+base;
 double h=0;
 for(int i=0;i<ncr;i++)h+=scale*rate[i];
 heat[cell]=h;
}

__global__ void dark_cooling_mol(Tables t,CoolTables c,const double *v,const int *entry,const Input *in,Output *out,CompactOutput *packed,const double *cr_heat){
 const int cell=blockIdx.x;
 if(out[cell].cooling_valid!=2)return;
 const Input &x=in[cell];
 const Output &y=out[cell];
 double cool=y.net_cooling,heat=cr_heat[cell],logT=log10(x.t);
 int ok=1,iT,iTmol;double dT,dTmol;
 index(v+t.offset[AXT],t.nt,logT,iT,dT);
 index(v+c.o_tmol,c.nTmol,logT,iTmol,dTmol);
 int xi;double xw;
 index(v+t.offset[AXX],t.nx,log10(dmaxd(x.x[sp_HII],DBL_MIN)),xi,xw);
 for(int i=0;i<2;i++){
   double cr=powten(interp(v+t.offset[CRSEC]+i*t.nx,xi,xw));
   heat-=3.2e-11f*y.rate[t.base[CR]+t.secondary[i]]*cr/(x.nh*(1.0f+cr));
 }
 heat+=0;
 cool+=1.017e-37*x.cmb*x.cmb*x.cmb*x.cmb*(x.t-x.cmb)*x.x[sp_elec]/x.nh;
 if(x.mol==1){
   const double xHI=x.x[sp_HI],xH2=x.x[sp_H2],xHII=x.x[sp_HII],xHe=x.x[sp_HeI],xe=x.x[sp_elec];
   double H2_low=powten(interp(v+c.o_h2h2,iTmol,dTmol))*xH2;
   H2_low+=powten(interp(v+c.o_h2hi,iTmol,dTmol))*xHI;
   H2_low+=powten(interp(v+c.o_h2hii,iTmol,dTmol))*xHII;
   H2_low+=powten(interp(v+c.o_h2hei,iTmol,dTmol))*xHe;
   H2_low+=powten(interp(v+c.o_h2e,iTmol,dTmol))*xe;
   H2_low*=xH2;
   if(H2_low>0.0f){
     double lte=powten(interp(v+c.o_h2lte,iTmol,dTmol))*xH2/x.nh;
     cool+=lte/(1.0f+(lte/H2_low));
   }
   cool+=7.2e-12f*y.coefficient[t.base[H2C]+t.cooling_index[1]]*xHI*xH2;
   cool+=7.2e-12f*y.coefficient[t.base[TD]+t.cooling_index[0]]*xH2*xH2;
   double H2_crit=(xHI+xH2==0)?0.0:(xHI+xH2)/((xHI/y.crit[0])+(xH2/y.crit[1]));
   heat+=((2.93e-12f*y.rate[t.base[TD]+t.cooling_index[2]])+(5.65e-12f*y.rate[t.base[CON]+t.cooling_index[3]]))*(1.0f/(x.nh+H2_crit));
   heat+=7.16e-12f*(y.rate[t.base[H2D]]/x.nh)*(x.nh/(x.nh+H2_crit));
   double mu=mean_mu(x.x,x.elem);
   double co_col=x.co_column,h2o_col=x.h2o_column,oh_col=x.oh_column;
   if(x.rewrite_columns){
     double column=x.cell_size*x.nh;
     co_col=x.x[sp_CO]*column;h2o_col=x.x[sp_H2O]*column;oh_col=x.x[sp_OH]*column;
   }
   if((x.elem_mask&5)==5){
     if(!c.have_co)ok=0;
     else{
       double logN=x.static_mol==1?log10(dmaxd(1.0e5f*co_col/thermal_dv(x.t,mu),DBL_MIN))
         :log10(dmaxd(1.0e5f*x.x[sp_CO]*x.nh/dmaxd(fabs(x.div_vel),DBL_MIN),DBL_MIN));
       int iR,iV;double dR,dV;
       index(v+c.o_corN,c.nCOrot,logN,iR,dR);
       index(v+c.o_covN,c.nCOvib,logN,iV,dV);
       double neff=x.nh*(xH2+9.857f*pow(x.t/1.0e3f,0.25f)*xHI+680.13f*pow(x.t,-0.25f)*xe);
       cool+=xH2*x.x[sp_CO]*rot_factor(v+c.o_corL0,v+c.o_corLlte,v+c.o_corNh,v+c.o_corA,iTmol,dTmol,iR,dR,c.nCOrot,neff);
       double vL0=powten(interp(v+c.o_covL0,iTmol,dTmol));
       double vLl=powten(bilin(v+c.o_covLlte,iTmol,iV,dTmol,dV,c.nCOvib));
       double vneff=x.nh*(xH2+50.0f*xHI+9035.09f*exp(68.0f/pow(x.t,1.0f/3.0f))*pow(x.t/300.0f,0.938f)*xe);
       cool+=xH2*x.x[sp_CO]/((1.0f/vL0)+(vneff/vLl));
     }
   }
   if(x.elem_mask&4){
     if(!c.have_h2o)ok=0;
     else{
       double logN=x.static_mol==1?log10(dmaxd(1.0e5f*h2o_col/thermal_dv(x.t,mu),DBL_MIN))
         :log10(dmaxd(1.0e5f*x.x[sp_H2O]*x.nh/dmaxd(fabs(x.div_vel),DBL_MIN),DBL_MIN));
       int iN,iH,iV;double dN,dH,dV;
       index(v+c.o_rotN,c.nH2Orot,logN,iN,dN);
       double neff=x.nh*(xH2+10.0f*xHI+powten(-8.02f+(15.749f/pow(x.t,1.0f/6.0f))-(47.137f/pow(x.t,1.0f/3.0f))+(76.648f/sqrt(x.t))-(60.191f/pow(x.t,2.0f/3.0f)))*xe/(7.4e-12f*sqrt(x.t)));
       if(logT>=2.0f){
         index(v+c.o_hiT,c.nH2Ohi,logT,iH,dH);
         cool+=x.x[sp_H2O]*xH2*rot_factor(v+c.o_hiL0,v+c.o_hiLlte,v+c.o_hiNh,v+c.o_hiA,iH,dH,iN,dN,c.nH2Orot,neff);
       }else{
         index(v+c.o_loT,c.nH2Olo,logT,iH,dH);
         cool+=0.75f*x.x[sp_H2O]*xH2*rot_factor(v+c.o_orL0,v+c.o_orLlte,v+c.o_orNh,v+c.o_orA,iH,dH,iN,dN,c.nH2Orot,neff);
         cool+=0.25f*x.x[sp_H2O]*xH2*rot_factor(v+c.o_paL0,v+c.o_paLlte,v+c.o_paNh,v+c.o_paA,iH,dH,iN,dN,c.nH2Orot,neff);
       }
       if(iH+1>=c.nTmol)ok=0;
       else{
         index(v+c.o_vibN,c.nH2Ovib,logN,iV,dV);
         double L0=powten(interp(v+c.o_vibL0,iH,dH));
         double Ll=powten(bilin(v+c.o_vibLlte,iH,iV,dH,dV,c.nH2Ovib));
         double vneff=x.nh*(xH2+10.0f*xHI+4.0625e8f*exp(47.5f/pow(x.t,1.0f/3.0f))*xe/sqrt(x.t));
         cool+=x.x[sp_H2O]*xH2/((1.0f/L0)+(vneff/Ll));
       }
       cool+=oh_rot(x.x,x.elem,x.t,x.nh,oh_col,x.extinction);
     }
   }
   cool+=powten(interp(v+c.o_gg,iT,dT))*x.dust*(x.t-x.td);
 }
 if(x.hybrid_mode==1 && !(x.t>x.tmol || x.dust==0)){
   int ig;double dg;
   index(v+t.offset[AXT],t.nt,log10(x.t),ig,dg);
   heat+=powten(interp(v+c.o_gg,ig,dg))*x.dust*(x.t-x.td);
 }
 double cool_s=cool*x.nh*x.nh;
 double heat_s=heat*x.nh*x.nh+x.constant_heating;
 double net=cool_s-heat_s;
 if(!isfinite(net))ok=0;
 out[blockIdx.x].net_cooling=net;
 out[blockIdx.x].cooling_valid=ok;
 if(packed){
   packed[blockIdx.x].net_cooling=net;
   packed[blockIdx.x].cooling_valid=ok && packed[blockIdx.x].valid;
 }
}

inline bool launch_dark_cooling(const Tables &t,const CoolTables &c,const double *v,const int *entry,const Input *in,Output *out,CompactOutput *packed,double *cr_heat,int n,cudaStream_t stream){
 dark_cooling_metal<<<n,128,0,stream>>>(t,c,v,entry,in,out,packed);
 if(cudaGetLastError()!=cudaSuccess)return false;
 dark_cooling_cr<<<n,1,0,stream>>>(t.base[CR],t.count[CR][0],t.count[CR][1],in,out,cr_heat);
 if(cudaGetLastError()!=cudaSuccess)return false;
 dark_cooling_mol<<<n,1,0,stream>>>(t,c,v,entry,in,out,packed,cr_heat);
 return cudaGetLastError()==cudaSuccess;
}

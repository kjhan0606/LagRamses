#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,double,double,double,double,double);
int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    omp_set_num_threads(2);
    if(lag_grackle_neq_init(argv[1],0,0,0,0,1.6726219e-24,1,1,1,8.5))return 3;
    enum { N=16 };
    const double densities[4]={1e-4,1e-2,1,1e2};
    const double energies[4]={1e10,1e12,1e14,1e16};
    double rho[N],u[N],s[9*N]={0},z[N],length[N];int cool[N];
    for(int i=0;i<N;i++) {
        double hfrac=.76, d=densities[i/4];
        double h=hfrac*(d-1e-3*d),he=(1-hfrac)*(d-1e-3*d);
        rho[i]=d;z[i]=1e-3*d;length[i]=1;u[i]=energies[i%4];cool[i]=1;
        s[0*N+i]=h*(1-2e-4-2e-6);s[1*N+i]=h*2e-4;s[3*N+i]=he;
        s[6*N+i]=h*2e-6;s[8*N+i]=s[1*N+i];
    }
    const double dts[]={1e10,1e12,1e14};
    for(int t=0;t<3;t++) {
        double before_u[N],before_s[9*N];memcpy(before_u,u,sizeof(u));memcpy(before_s,s,sizeof(s));
        int status=lag_grackle_neq_step(N,rho,u,s,z,length,cool,1.6726219e-24,1,1,1,dts[t]);
        if(status) {
            memcpy(u,before_u,sizeof(u));memcpy(s,before_s,sizeof(s));
            printf("matrix dt=%.1e: controlled solver rejection status=%d\n",dts[t],status);
            continue;
        }
        double maxdrift=0,max_energy_change=0;
        for(int i=0;i<N;i++) {
            double h=s[i]+s[N+i]+s[2*N+i]+s[6*N+i]+s[7*N+i];
            double he=s[3*N+i]+s[4*N+i]+s[5*N+i];
            double charge=s[N+i]-s[2*N+i]+s[7*N+i]/2+s[4*N+i]/4+s[5*N+i]/2;
            double drift=fmax(fabs(h+he+z[i]-rho[i]),fabs(s[8*N+i]-charge))/rho[i];
            if(!isfinite(u[i])||u[i]<=0||drift>1e-6)return 4;
            max_energy_change=fmax(max_energy_change,fabs(u[i]-before_u[i])/before_u[i]);
            for(int k=0;k<9;k++)if(!isfinite(s[k*N+i])||s[k*N+i]<0)return 5;
            maxdrift=fmax(maxdrift,drift);
        }
        printf("matrix dt=%.1e: PASS max_budget_drift=%.3e max_relative_energy_change=%.3e\n",
               dts[t],maxdrift,max_energy_change);
    }
    return 0;
}

#include <stdio.h>
#include <math.h>
int lag_grackle_neq_init(const char *,int,int,int,int,double,double,double,double,double);
int lag_grackle_neq_step(int,double *,double *,double *,double *,double *,int *,double,double,double,double,double);
int main(int argc,char **argv)
{
    if(argc!=2)return 2;
    if(lag_grackle_neq_init(argv[1],1,0,0,1,1.6726219e-24,3.0856776e21,
                            3.15576e13,0.1,8.5))return 3;
    double rho[1]={1e-3},z[1]={1e-8},length[1]={1e-5},cool[1]={1};
    double energy[2]={1e12,1e12},species[18]={0};
    for(int j=0;j<2;j++) {
        double h=.76*(rho[0]-z[0]),he=.24*(rho[0]-z[0]);
        double *s=species+9*j;
        s[0]=h*(1-2e-4-2e-6);s[1]=h*2e-4;s[3]=he;
        s[6]=h*2e-6;s[8]=s[1];
    }
    int c[1]={1};
    if(lag_grackle_neq_step(1,rho,energy,species,z,length,c,1.6726219e-24,
                            3.0856776e21,3.15576e13,.1,1e-4))return 4;
    if(lag_grackle_neq_step(1,rho,energy+1,species+9,z,length,c,1.6726219e-24,
                            3.0856776e21,3.15576e13,.5,1e-4))return 5;
    double df=0;
    for(int k=0;k<9;k++)df=fmax(df,fabs(species[k]-species[9+k])/rho[0]);
    double de=fabs(energy[0]-energy[1])/energy[0];
    printf("Grackle UV scale-factor comparison: a=0.1 versus 0.5; relative dE=%.9e dspecies_max=%.9e\n",de,df);
    if(!isfinite(de)||!isfinite(df)||de+df<1e-12)return 6;
    return 0;
}

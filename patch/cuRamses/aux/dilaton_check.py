#!/usr/bin/env python3
"""Validate the Brax+12 dilaton as coded: background fixed point,
linear response mu(k), unscreened limit 2*beta(a)^2, screening."""
import numpy as np
N=32; dx=1.0/N
i=np.arange(N); kd1=(2-2*np.cos(2*np.pi*i/N))/dx**2
Om=0.3089; B=100.0; B2=(B/2997.92458)**2
beta0=0.5; L=5.0; xi=L/2997.92458; A2=1/(3*xi**2); s=3*Om; pexp=1-3/s

def S(chi,rho,a,chibar):
    cA=3*Om*A2*B2/a; cV=a*a*B2
    w=np.maximum(A2*chi/beta0,1e-30)
    v=-3*Om*beta0*w**pexp
    vbar=-3*Om*beta0*(A2*chibar/beta0)**pexp
    dS=cA*rho - 3*Om*A2*pexp*w**(pexp-1)*cV
    return cA*(rho*chi-chibar)+cV*(v-vbar), dS
def lap(u): return (np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)+np.roll(u,1,2)+np.roll(u,-1,2)-6*u)/dx**2

for a in [1.0,0.5,0.25]:
    chibar=beta0*a**s/A2
    # 1) background fixed point
    chi=np.full((N,N,N),chibar); rho=np.ones((N,N,N))
    r0,_=S(chi,rho,a,chibar); print(f"a={a}: bg residual={np.max(np.abs(lap(chi)-r0)):.2e}", end='')
    # 2) linear response
    x=(np.arange(N)+0.5)*dx; A=1e-6
    for km in [2,8]:
        delta=A*np.sin(2*np.pi*km*x)[:,None,None]*np.ones((N,N,N)); delta-=delta.mean()
        rho=1+delta; chi=np.full((N,N,N),chibar)
        for it in range(3000):
            Sv,dS=S(chi,rho,a,chibar)
            res=lap(chi)-Sv
            chi=chi-res/(-6/dx**2-dS)
        k2d=(2-2*np.cos(2*np.pi*km*dx))/dx**2
        # F5 = -c~^2 a^2 A2 chi grad chi ; c~^2 = 1/B2
        gchi=np.gradient(chi,dx,axis=0)
        phik=-1.5*Om*a*np.fft.rfftn(delta)/np.where(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1]>0,kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1],1)
        phi=np.fft.irfftn(phik,s=(N,N,N),axes=(0,1,2))
        gphi=np.gradient(phi,dx,axis=0)
        m=np.abs(gphi)>0.3*np.max(np.abs(gphi))
        F5=-(1/B2)*a*a*A2*chi[m]*gchi[m]
        meas=np.median(F5/(-gphi[m]))
        m2t=3*A2*B2/a
        th=2*(beta0*a**s)**2 * k2d/(k2d+m2t)
        print(f"  k={km}: F5/FN={meas:.5f} theory={th:.5f} (2b^2={2*(beta0*a**s)**2:.4f})", end='')
    print()
# 3) screening: dense top-hat
a=1.0; chibar=beta0/A2
x=(np.arange(N)+0.5)*dx
r2=((x-0.5)[:,None,None]**2+(x-0.5)[None,:,None]**2+(x-0.5)[None,None,:]**2)
delta=np.where(r2<0.06**2,1e5,0.0); delta-=delta.mean(); rho=1+delta
chi=np.full((N,N,N),chibar)
for it in range(6000):
    Sv,dS=S(chi,rho,a,chibar)
    res=lap(chi)-Sv
    dchi=-res/(-6/dx**2-dS)
    dchi=np.clip(dchi,-0.5*np.abs(chi),0.5*np.abs(chi))
    chi=np.maximum(chi+dchi,1e-12*chibar)
print(f"screening: chi_center/chibar = {chi[N//2,N//2,N//2]/chibar:.3e} (should be <<1)")
gchi=np.gradient(chi,dx,axis=0)
phik=-1.5*Om*np.fft.rfftn(delta)/np.where(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1]>0,kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1],1)
phi=np.fft.irfftn(phik,s=(N,N,N),axes=(0,1,2)); gphi=np.gradient(phi,dx,axis=0)
m=np.abs(gphi)>0.5*np.max(np.abs(gphi))
print(f"screening: median F5/FN near tophat = {np.median(-(1/B2)*A2*chi[m]*gchi[m]/(-gphi[m])):.4f} (unscreened 0.5)")

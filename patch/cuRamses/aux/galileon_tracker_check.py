#!/usr/bin/env python3
"""Validate the tracker cubic Galileon as coded: operator-split solve
with tracker coefficients; measured unscreened F5/FN vs the analytic
-xi/(9 beta2 E^2); also print the Geff(a) curve and screened behavior."""
import numpy as np
N=32; dx=1.0/N
i=np.arange(N); kd1=(2-2*np.cos(2*np.pi*i/N))/dx**2
def helm(b,m2):
    bk=np.fft.rfftn(b)
    den=-(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1]+m2)
    xk=np.zeros_like(bk)
    np.divide(bk,den,out=xk,where=np.abs(den)>1e-30)
    return np.fft.irfftn(xk,s=(N,N,N),axes=(0,1,2))
def sd(u):
    xx=(np.roll(u,1,0)+np.roll(u,-1,0)-2*u)/dx**2
    yy=(np.roll(u,1,1)+np.roll(u,-1,1)-2*u)/dx**2
    zz=(np.roll(u,1,2)+np.roll(u,-1,2)-2*u)/dx**2
    return xx,yy,zz

Om=0.3; xi=np.sqrt(6*(1-Om))
def coeffs(a):
    E2=0.5*(Om/a**3+np.sqrt((Om/a**3)**2+4*(1-Om)))
    hd=-1.5*(Om/a**3)/(2*E2-Om/a**3)
    b1=(xi/3)*(2*hd-1+(1-Om)/E2**2)
    b2=2*E2*b1/xi**2
    return E2,b1,b2

x=(np.arange(N)+0.5)*dx
print(f"{'a':>6} {'E2':>8} {'beta1':>9} {'beta2':>9} {'F5/FN meas':>11} {'analytic':>9}")
for a,amp in [(1.0,1e-4),(0.5,1e-4),(0.25,1e-4)]:
    E2,b1,b2=coeffs(a)
    coeff=1/(3*b1*a**4); srcfac=Om*a/b2
    delta=amp*np.sin(2*np.pi*2*x)[:,None,None]*np.ones((N,N,N)); delta-=delta.mean()
    u=np.zeros((N,N,N))
    for it in range(6):
        xx,yy,zz=sd(u); l=xx+yy+zz; T=xx**2+yy**2+zz**2
        S=srcfac*delta
        disc=1+4*coeff*(S+coeff*T)
        A=np.where(disc>0,(-1+np.sqrt(np.maximum(disc,0)))/(2*coeff),-1/(2*coeff))
        u=u+helm(A-l,0.0)
    u-=u.mean()
    # F5 = +xi/(6E2) grad u ; FN = -grad phi, lap phi = 1.5 Om a delta
    phi=helm(1.5*Om*a*delta,0.0)*-1  # helm solves (lap)x=b as x=-b/k2... careful
    # direct: lap phi = 1.5 Om a delta -> phik = -(1.5 Om a) deltak/kd2
    dk=np.fft.rfftn(delta); den=(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1])
    phik=np.zeros_like(dk)
    np.divide(-1.5*Om*a*dk,den,out=phik,where=den>0)
    phi=np.fft.irfftn(phik,s=(N,N,N),axes=(0,1,2))
    gu=np.gradient(u,dx,axis=0); gp=np.gradient(phi,dx,axis=0)
    m=np.abs(gp)>0.3*np.max(np.abs(gp))
    ratio=np.median((xi/(6*E2))*gu[m]/(-gp[m]))
    print(f"{a:6.2f} {E2:8.4f} {b1:9.4f} {b2:9.4f} {ratio:11.5f} {-xi/(9*b2*E2):9.5f}")

# screened test: strong top-hat at a=1
E2,b1,b2=coeffs(1.0); coeff=1/(3*b1); srcfac=Om/b2
r2=((x-0.5)[:,None,None]**2+(x-0.5)[None,:,None]**2+(x-0.5)[None,None,:]**2)
delta=np.where(r2<0.05**2,5e4,0.0); delta-=delta.mean()
u=np.zeros((N,N,N))
for it in range(30):
    xx,yy,zz=sd(u); l=xx+yy+zz; T=xx**2+yy**2+zz**2
    S=srcfac*delta; disc=1+4*coeff*(S+coeff*T)
    A=np.where(disc>0,(-1+np.sqrt(np.maximum(disc,0)))/(2*coeff),-1/(2*coeff))
    u=0.5*u+0.5*(u+helm(A-l,0.0))
dk=np.fft.rfftn(delta); den=(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1])
phik=np.zeros_like(dk)
np.divide(-1.5*Om*dk,den,out=phik,where=den>0)
phi=np.fft.irfftn(phik,s=(N,N,N),axes=(0,1,2))
gu=np.gradient(u,dx,axis=0); gp=np.gradient(phi,dx,axis=0)
m=np.abs(gp)>0.5*np.max(np.abs(gp))
print(f"screened top-hat (delta=5e4): median F5/FN = {np.median((xi/(6*E2))*gu[m]/(-gp[m])):.4f}  (unscreened would be {-xi/(9*b2*E2):.3f})")

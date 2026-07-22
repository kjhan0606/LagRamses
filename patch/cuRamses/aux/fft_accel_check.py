#!/usr/bin/env python3
"""
Validate the b2 FFT acceleration schemes exactly as coded in
force_fine.kjhan.f90 (level_fft_helmholtz + builders), on 32^3:

 1) f(R): 3 Newton-FFT iterations (spectral (lap-m2bar)^-1 with the
    level-mean Jacobian) vs a 6000-sweep Newton-GS reference.
 2) nDGP Vainshtein: 5 operator-split iterations (cell quadratic for
    A=lap(phi), then spectral lap^-1) vs a damped Newton-GS reference.

Same 7-point stencil, same discrete eigenvalues, same branch/clamp.
"""
import numpy as np
rng = np.random.default_rng(7)

N=32; dx=1.0/N
def lap(u):
    return (np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)
            +np.roll(u,1,2)+np.roll(u,-1,2)-6*u)/dx**2

i = np.arange(N)
kd1 = (2-2*np.cos(2*np.pi*i/N))/dx**2
KD2 = kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:]

x=(np.arange(N)+0.5)*dx
delta = 0.8*np.sin(2*np.pi*x)[:,None,None]*np.ones((N,N,N)) \
      + 0.3*np.sin(2*np.pi*4*x)[None,:,None] + 0.05*rng.standard_normal((N,N,N))
delta -= delta.mean()

def helm_solve(b, m2):
    bk = np.fft.rfftn(b)
    den = -(kd1[:,None,None]+kd1[None,:,None]+kd1[None,None,:N//2+1]+m2)
    with np.errstate(divide='ignore', invalid='ignore'):
        xk = np.where(np.abs(den)>1e-30, bk/den, 0.0)
    return np.fft.irfftn(xk, s=(N,N,N))

# ---------------- Test 1: f(R) ----------------
Om,OL=0.3,0.7; L=250.0; B=(L/2997.92458)**2; a=1.0; fR0=-1e-5; n=1; np1=n+1
Rbar=3*(Om/a**3+4*OL); Rbar0=3*(Om+4*OL); fRbar=fR0*(Rbar0/Rbar)**np1
a23B=a*a*B/3; rc=Om*B/a

def S_fr(u):
    R=Rbar0*(abs(fR0)/np.abs(u))**(1/np1)
    return a23B*(R-Rbar)-rc*delta, a23B*(-R/(np1*u))

# reference: Newton-GS (Jacobi-style sweeps, definite operator)
u=np.full((N,N,N),fRbar)
for it in range(6000):
    S,dS=S_fr(u); res=lap(u)-S
    u=u-res/(-6/dx**2-dS)
S,_=S_fr(u); ref_res=np.max(np.abs(lap(u)-S)); u_ref=u.copy()

# FFT-Newton exactly as coded: 3 iters
u=np.full((N,N,N),fRbar)
for it in range(3):
    S,dS=S_fr(u)
    r=lap(u)-S
    m2bar=max(np.mean(dS),0.0)
    u=u+helm_solve(-r, m2bar)
err=np.max(np.abs(u-u_ref))/np.max(np.abs(u_ref-fRbar))
S,_=S_fr(u); fft_res=np.max(np.abs(lap(u)-S))
print(f"f(R): |u_fft-u_ref|/|du_ref| = {err:.2e}  (residual: fft={fft_res:.2e}, ref={ref_res:.2e})")

# also 20 plain GS from cold start for contrast
u=np.full((N,N,N),fRbar)
for it in range(20):
    S,dS=S_fr(u); res=lap(u)-S
    u=u-res/(-6/dx**2-dS)
err20=np.max(np.abs(u-u_ref))/np.max(np.abs(u_ref-fRbar))
print(f"f(R): 20 plain GS-from-cold error = {err20:.2e}  (FFT wins by {err20/max(err,1e-30):.1e}x)")

# ---------------- Test 2: nDGP Vainshtein ----------------
Orc=0.25
E2=Om/a**3+OL; E=np.sqrt(E2)
Hd_H2=-1.5*Om/a**3/E2
beta=1+(E/np.sqrt(Orc))*(1+Hd_H2/3)   # 2*H*rc = E/sqrt(Orc)
coeff=1/(12*Orc*beta*a**4); srcfac=Om*a/beta
Ssrc=srcfac*delta

def second_derivs(u):
    uxx=(np.roll(u,1,0)+np.roll(u,-1,0)-2*u)/dx**2
    uyy=(np.roll(u,1,1)+np.roll(u,-1,1)-2*u)/dx**2
    uzz=(np.roll(u,1,2)+np.roll(u,-1,2)-2*u)/dx**2
    return uxx,uyy,uzz

# reference: damped Newton-GS, many sweeps
u=np.zeros((N,N,N))
for it in range(20000):
    uxx,uyy,uzz=second_derivs(u); l=uxx+uyy+uzz
    vain=l**2-(uxx**2+uyy**2+uzz**2)
    res=l+coeff*vain-Ssrc
    J=-6/dx**2+coeff*2*l*(-4/dx**2)
    du=-res/np.where(np.abs(J)>1e-30,J,1e-30)
    sc=0.5*dx**2*np.maximum(np.abs(Ssrc),1e-2*Om*a/beta)
    du=np.clip(du,-sc,sc)
    u=u+du
uxx,uyy,uzz=second_derivs(u); l=uxx+uyy+uzz
ref_res=np.max(np.abs(l+coeff*(l**2-(uxx**2+uyy**2+uzz**2))-Ssrc)); u_ref=u-u.mean()

# operator-split FFT, 5 iters, exactly as coded
u=np.zeros((N,N,N))
for it in range(5):
    uxx,uyy,uzz=second_derivs(u); l=uxx+uyy+uzz
    T=uxx**2+uyy**2+uzz**2
    disc=1+4*coeff*(Ssrc+coeff*T)
    A=np.where(disc>0,(-1+np.sqrt(np.maximum(disc,0)))/(2*coeff),-1/(2*coeff))
    u=u+helm_solve(A-l, 0.0)
u=u-u.mean()
uxx,uyy,uzz=second_derivs(u); l=uxx+uyy+uzz
fft_res=np.max(np.abs(l+coeff*(l**2-(uxx**2+uyy**2+uzz**2))-Ssrc))
gref=np.gradient(u_ref,dx,axis=0); gfft=np.gradient(u,dx,axis=0)
gerr=np.max(np.abs(gfft-gref))/np.max(np.abs(gref))
print(f"nDGP: force-field rel err (5 split iters vs 20000-sweep ref) = {gerr:.2e}")
print(f"nDGP: residuals fft={fft_res:.3e} ref={ref_res:.3e}  beta={beta:.3f} coeff={coeff:.3f}")

#!/usr/bin/env python3
"""
Validate the FIXED f(R) discrete equation exactly as now coded in
fR_gauss_seidel (force_fine.kjhan.f90):

  lapl u = +(a^2/3)*B*(R(u)-Rbar) - Om*B*(rho-1)/a,   B=(L/2997.92458)^2
  R(u) = Rbar0*(|fR0|/|u|)^(1/(n+1)),  Newton J = -6/dx^2 - (a^2/3)*B*dR/du
  F5 = +0.5*a^2/B * grad(u)

Linear theory target (discrete k):
  mu(k) = 1 + (1/3)*k2/(k2 + m2),  m2 = (a^2/3)*B*Rbar/((n+1)*|fR_bar|)
"""
import numpy as np

Om, OL = 0.3, 0.7
L = 250.0           # Mpc/h
B = (L/2997.92458)**2
N = 32
dx = 1.0/N

def run(a, fR0, n, kmode, A=1e-8):
    np1 = n+1
    Rbar  = 3*(Om/a**3 + 4*OL)
    Rbar0 = 3*(Om + 4*OL)
    fRbar = fR0*(Rbar0/Rbar)**np1
    x = (np.arange(N)+0.5)*dx
    delta1d = A*np.sin(2*np.pi*kmode*x)
    delta = np.broadcast_to(delta1d[:,None,None], (N,N,N)).copy()

    u = np.full((N,N,N), fRbar)
    a23B = a*a*B/3.0
    rc = Om*B/a
    dx2i = 1.0/dx**2
    for it in range(4000):
        lap = (np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)
               +np.roll(u,1,2)+np.roll(u,-1,2)-6*u)*dx2i
        R_of_u = Rbar0*(abs(fR0)/np.abs(u))**(1.0/np1)
        S = a23B*(R_of_u - Rbar) - rc*delta
        res = lap - S
        dRdu = -R_of_u/(np1*u)
        J = -6*dx2i - a23B*dRdu
        u = u - res/J
        if np.max(np.abs(res)) < 1e-16*max(np.max(np.abs(S)),1e-30):
            break

    # measured mu at the mode: F5/FN via the sine amplitudes
    du = np.mean(u - np.mean(u), axis=(1,2))
    amp_u = 2*np.mean(du*np.sin(2*np.pi*kmode*x))
    # discrete k^2 of the 7-point stencil
    k2d = (2 - 2*np.cos(2*np.pi*kmode*dx))/dx**2
    # Newtonian potential amplitude (same stencil): phi = -1.5 Om a A /k2d * sin
    amp_phi = -1.5*Om*a*A/k2d
    # F_N = -grad(phi), F5 = +(a^2/2B) grad(u): ratio needs the minus
    F5_over_FN = -(0.5*a*a/B)*amp_u/amp_phi
    m2 = a23B*Rbar/(np1*abs(fRbar))
    theory = (1.0/3.0)*k2d/(k2d + m2)
    return F5_over_FN, theory, np.sqrt(m2)*L/(2*np.pi*N*0+1)  # m in box^-1

print(f"{'a':>5} {'fR0':>8} {'n':>2} {'k':>3}   {'F5/FN meas':>12} {'theory':>12} {'rel.err':>9}")
for (a, fR0, n) in [(1.0,-1e-5,1),(0.5,-1e-5,1),(1.0,-1e-6,1),(0.25,-1e-5,1),(1.0,-1e-5,2)]:
    for kmode in [1,4,8]:
        got, th, _ = run(a, fR0, n, kmode)
        print(f"{a:5.2f} {fR0:8.0e} {n:2d} {kmode:3d}   {got:12.6f} {th:12.6f} {abs(got/th-1):9.2e}")

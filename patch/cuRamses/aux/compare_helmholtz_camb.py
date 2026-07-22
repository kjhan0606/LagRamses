#!/usr/bin/env python3
"""
Verify the LagRamses CPL DE-perturbation kappa2/alpha Helmholtz fallback
against CAMB, using the SAME Weyl-potential method that defines the
production de_table (patch/cuRamses/aux/generate_de_table.py):

    de_factor_CAMB(k,a) = (Weyl/delta_m)_[cs2=cs2_de] / (Weyl/delta_m)_[cs2=1]

Candidates compared:
  A (code today) : (k^2+kap^2)/(k^2+kap^2+alp),  kap^2=a^2E^2/(cs2 C^2),
                   alp=1.5 Ol f (1+w)/cs2 /C^2          [suppression]
  B (sign flip)  : (k^2+kap^2)/(k^2+kap^2-alp)           [enhancement]
  C (Sapone-Kunz): 1 + (Ode_a/Om) (1+w) / [(1-3w) + 2 cs2 C^2 k^2/(3 a^2 E^2)]

k in h/Mpc (comoving), C = c/H0 in Mpc/h = 2997.92458.
"""
import numpy as np
import camb
from camb.dark_energy import DarkEnergyFluid

C = 2997.92458  # c/H0 in Mpc/h

omega_b, omega_m, h0 = 0.0486, 0.3089, 0.6774
omega_l = 1.0 - omega_m
ns, As = 0.9667, 2.142e-9

def f_de(a, w0, wa):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E2(a, w0, wa):
    return omega_m/a**3 + omega_l*f_de(a, w0, wa)

def run_camb(w0, wa, cs2, zs, kmax):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h0*100, ombh2=omega_b*h0**2,
                       omch2=(omega_m-omega_b)*h0**2, mnu=0.0)
    pars.InitPower.set_params(ns=ns, As=As)
    pars.DarkEnergy = DarkEnergyFluid(w=w0, wa=wa, cs2=cs2)
    pars.set_matter_power(redshifts=sorted(zs, reverse=True), kmax=kmax*1.5)
    pars.WantTransfer = True
    return camb.get_results(pars)

def weyl_per_dm(results, q, z):
    ev = results.get_redshift_evolution(q, [z], ['Weyl','delta_cdm','delta_baryon'])
    weyl = ev[:,0,0]
    dm = ((omega_m-omega_b)*ev[:,0,1] + omega_b*ev[:,0,2])/omega_m
    return weyl/dm

def candidates(k, a, w0, wa, cs2):
    w  = w0 + wa*(1-a)
    fd = f_de(a, w0, wa)
    e2 = E2(a, w0, wa)
    ode_a = omega_l*fd*a**3          # in units of Omega_m-normalization used by code
    out = {}
    if cs2 > 0:
        kap2 = a*a*e2/(cs2*C*C)
        alp  = 1.5*omega_l*fd*(1+w)/(cs2*C*C)
        out['A_code'] = (k**2+kap2)/(k**2+kap2+alp)
        out['B_flip'] = (k**2+kap2)/(k**2+kap2-alp)
    else:
        out['A_code'] = np.full_like(k, 1.0/(1.0+1.5*omega_l*fd*(1+w)/(a*a*e2)))
        out['B_flip'] = np.full_like(k, 1.0/(1.0-1.5*omega_l*fd*(1+w)/(a*a*e2)))
    denom = (1-3*w) + (2*cs2*C*C*k**2)/(3*a*a*e2)
    out['C_SK'] = 1.0 + (ode_a/omega_m)*(1+w)/denom * np.ones_like(k)
    return out

def main():
    ktab = np.logspace(-3, 0.3, 40)   # h/Mpc
    q = ktab*h0                        # 1/Mpc for CAMB
    zs = [0.0, 0.5, 1.0, 2.0, 5.0]

    for (w0, wa, cs2) in [(-0.8, 0.0, 1e-4), (-0.8, 0.0, 0.0),
                          (-0.9, 0.1, 0.0), (-0.8, 0.0, 1e-2)]:
        print(f"\n================ w0={w0} wa={wa} cs2={cs2} ================")
        rt = run_camb(w0, wa, cs2, zs, ktab[-1])
        rs = run_camb(w0, wa, 1.0, zs, ktab[-1])
        for z in zs:
            a = 1.0/(1.0+z)
            fac_camb = weyl_per_dm(rt, q, z)/weyl_per_dm(rs, q, z)
            cand = candidates(ktab, a, w0, wa, cs2)
            # summarize at three k
            for kk in [1e-3, 1e-2, 1e-1]:
                i = np.argmin(np.abs(ktab-kk))
                print(f" z={z:3.1f} k={ktab[i]:8.1e}  CAMB={fac_camb[i]:8.5f}"
                      f"  A_code={cand['A_code'][i]:8.5f}"
                      f"  B_flip={cand['B_flip'][i]:8.5f}"
                      f"  C_SK={cand['C_SK'][i]:8.5f}")
            # error metrics over k in [3e-3, 1] (box-relevant scales)
            m = (ktab >= 3e-3) & (ktab <= 1.0)
            for name in ['A_code','B_flip','C_SK']:
                err = np.max(np.abs(cand[name][m]-fac_camb[m]))
                print(f"    max|{name}-CAMB| (k in [3e-3,1]) = {err:9.5f}")

if __name__ == '__main__':
    main()

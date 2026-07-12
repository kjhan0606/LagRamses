#!/usr/bin/env python3
"""
Exact linear forward model of the shared-IC protocol:
- All models read the SAME velocity files. RAMSES converts them to
  displacements with the MODEL's vfact => delta_i^M = (vL/vM) delta_L,
  while initial velocities (hence ddelta/dlna at a_i) are identical.
- vfact is computed with the code's own formulas (fpeebl integral
  form, mfac for coupled DE).
- Then evolve delta'' + (2+dlnE+fric) delta' = 1.5 Om_src(a) mu delta
  from (delta_i, f_L*delta_L) and form P_M/P_L at the run's aexp.
"""
import numpy as np, sys
sys.path.insert(0,'/tmp/desim')
from scipy.integrate import quad, solve_ivp
from linear_theory import quint_tables

om, ol = 0.3089, 0.6911
ai = 0.02

# ---------- backgrounds ----------
lnaq,usq,vsq,e2q,cmq = quint_tables(0.0)
b=0.1
lnac,usc,vsc,e2c,cmc = quint_tables(b)
u0c=usc[-1]

def fde_tab(lna,e2,rhoc_last=None,cm=None,us=None,bb=0.0):
    rhoc = (cm if cm else om)*np.exp(-3*lna-(bb*us if us is not None else 0))
    rphi = e2-rhoc
    return rphi/rphi[-1]
fq = fde_tab(lnaq,e2q,cm=cmq)
fc = fde_tab(lnac,e2c,cm=cmc,us=usc,bb=b)

x0=0.5005; C0=(2*x0-1)*np.sqrt(x0); rho0k=3*x0**2-x0
def kes_fw(a):
    targ=C0/a**3; x=max(x0,(targ/2)**(2/3))
    for _ in range(80):
        g=(2*x-1)*np.sqrt(x)-targ; gp=2*np.sqrt(x)+(2*x-1)/(2*np.sqrt(x))
        x=max(0.5+1e-13,x-g/gp)
    return (3*x*x-x)/rho0k, (x*x-x)/(3*x*x-x)

def mk_model(name):
    """returns dict: fde(a), w(a), mfac(a) [in y], E2(a), src_omfac(a), mu(a), fric(a)"""
    one=lambda a:1.0; zero=lambda a:0.0
    if name=='lcdm':
        return dict(fde=one,w=lambda a:-1.0,mf=one,
                    E2=lambda a: om/a**3+ol, omf=one, mu=one, fric=zero)
    if name=='quint':
        fdef=lambda a: np.interp(np.log(a),lnaq,fq)
        wq=lambda a: wtab(a,lnaq,e2q,vsq,usq,cmq,0.0)
        return dict(fde=fdef,w=wq,mf=one,
                    E2=lambda a: np.interp(np.log(a),lnaq,e2q),omf=one,mu=one,fric=zero)
    if name in ('cdeA','cde'):
        fdef=lambda a: np.interp(np.log(a),lnac,fc)
        wc=lambda a: wtab(a,lnac,e2c,vsc,usc,cmc,b)
        dm=lambda a: np.exp(-b*(np.interp(np.log(a),lnac,usc)-u0c))
        E2=lambda a: np.interp(np.log(a),lnac,e2c)
        if name=='cdeA':
            return dict(fde=fdef,w=wc,mf=dm,E2=E2,omf=one,mu=lambda a:1+2*b*b,fric=zero)
        return dict(fde=fdef,w=wc,mf=dm,E2=E2,omf=dm,mu=lambda a:1+2*b*b,
                    fric=lambda a:-b*np.interp(np.log(a),lnac,vsc))
    if name=='kes':
        return dict(fde=lambda a:kes_fw(a)[0],w=lambda a:kes_fw(a)[1],mf=one,
                    E2=lambda a: om/a**3+ol*kes_fw(a)[0],omf=one,mu=one,fric=zero)
    if name=='hs':
        E2=lambda a: om/a**3+ol
        return dict(fde=one,w=lambda a:-1.0,mf=one,E2=E2,omf=one,
                    mu=lambda a:1+0.2/E2(a),fric=zero)
    if name=='gal':
        xi=np.sqrt(6*(1-om))
        E2=lambda a: 0.5*(om/a**3+np.sqrt((om/a**3)**2+4*(1-om)))
        def mu(a):
            e2=E2(a); hd=-1.5*(om/a**3)/(2*e2-om/a**3)
            b1=(xi/3)*(2*hd-1+(1-om)/e2**2); b2=2*e2*b1/xi**2
            return 1-xi/(9*b2*e2)
        fdef=lambda a: 1.0/E2(a)
        def wg(a):
            e2=E2(a); hd=-1.5*(om/a**3)/(2*e2-om/a**3)
            return -1+2*hd/3
        return dict(fde=fdef,w=wg,mf=one,E2=E2,omf=one,mu=mu,fric=zero)
    if name=='dil':
        beta0,L=0.5,20.0
        xi=L/2997.92458; A2=1/(3*xi**2); s=3*om
        E2=lambda a: om/a**3+ol
        # mu(a,k): band-average later; here mu at a given comoving k [h/Mpc]
        return dict(fde=one,w=lambda a:-1.0,mf=one,E2=E2,omf=one,
                    mu=None,fric=zero,dil=(beta0,A2,s))

def wtab(a,lna,e2s,vs,us,cm,bb):
    x=np.log(a); e2=np.interp(x,lna,e2s); v=np.interp(x,lna,vs)
    rphi=e2-cm*np.exp(-3*x-bb*np.interp(x,lna,us))
    V=rphi-e2*v*v/6
    return (e2*v*v/6-V)/rphi

def vfact_code(M):
    """the code's vfact at ai: a*fpeebl*sqrt(y), y=om*mf/a+ol*fde*a^2"""
    y=lambda a: om*M['mf'](a)/a + ol*M['fde'](a)*a*a
    fy=lambda a: 1.0/y(a)**1.5
    fact=quad(fy,1e-6,ai,limit=200)[0]
    wa=M['w'](ai)
    dyda=-om*M['mf'](ai)/ai + ol*ai*ai*M['fde'](ai)*(-1-3*wa)
    fp=0.5*dyda/y(ai)-1+ai*fy(ai)/fact
    return ai*fp*np.sqrt(y(ai))

def evolve(M, a_end, mu_k=None):
    mu = mu_k if mu_k else M['mu']
    def rhs(x,yv):
        a=np.exp(x); e=1e-4
        dlnE=(np.log(M['E2'](a*np.exp(e)))-np.log(M['E2'](a*np.exp(-e))))/(4*e)
        return [yv[1], -(2+dlnE+M['fric'](a))*yv[1]
                + 1.5*(om/a**3/M['E2'](a))*M['omf'](a)*mu(a)*yv[0]]
    return rhs

# LCDM reference amplitude & growth rate at ai
L=mk_model('lcdm')
sL=solve_ivp(evolve(L,1.0),[np.log(1e-3),np.log(ai)],[1e-3,1e-3],rtol=1e-9,atol=1e-13)
dL_i, dLp_i = sL.y[0][-1], sL.y[1][-1]
fL = dLp_i/dL_i
vL = vfact_code(L)

def predict(name, a_end, mu_k=None):
    M=mk_model(name)
    vM=vfact_code(M)
    r=vL/vM
    y0=[r*dL_i, fL*dL_i]      # same velocities, rescaled displacements
    sM=solve_ivp(evolve(M,a_end,mu_k),[np.log(ai),np.log(a_end)],y0,rtol=1e-9,atol=1e-13)
    sL2=solve_ivp(evolve(L,a_end),[np.log(ai),np.log(a_end)],[dL_i,fL*dL_i],rtol=1e-9,atol=1e-13)
    return (sM.y[0][-1]/sL2.y[0][-1])**2, vM/vL

print(f"{'model':6} {'a':>6} {'pred P/P_L':>10} {'vM/vL':>7}")
for name,aend in [('quint',1.012),('kes',1.007),('hs',1.002),
                  ('cdeA',0.504),('cde',0.502),('cde',1.003),('gal',1.005)]:
    p,vr=predict(name,aend)
    print(f"{name:6} {aend:6.3f} {p:10.4f} {vr:7.4f}")
# dilaton: band-average mu over k in [0.05,0.12]
Md=mk_model('dil'); beta0,A2,s=Md['dil']
ks=np.linspace(0.05,0.12,8)
def mu_dil(a):
    m2=3*A2/(2997.92458**2*a)
    return float(np.mean(1+2*(beta0*a**s)**2*ks**2/(ks**2+m2)))
p,vr=predict('dil',1.0,mu_k=mu_dil)
print(f"{'dil':6} {1.0:6.3f} {p:10.4f} {vr:7.4f}")

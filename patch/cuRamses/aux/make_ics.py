#!/usr/bin/env python3
"""
Zel'dovich grafic ICs for the DE-model test suite.
64^3 particles, 100 Mpc/h box, z_start=49, Planck-ish LCDM.
Writes ic_velcx/y/z (standard grafic: header + n3 plane records, float32).
Same random seed for every model — model differences come from dynamics.
"""
import numpy as np, struct, os, sys
import camb

N = 64
Lbox_hMpc = 100.0          # Mpc/h
h = 0.6774
om, ol = 0.3089, 0.6911
zs = 49.0
astart = 1.0/(1.0+zs)
seed = 20260713

Lbox_Mpc = Lbox_hMpc / h   # grafic dxini is in Mpc
dx_Mpc = Lbox_Mpc / N

# --- P(k) at z=49 from CAMB (LCDM; shared ICs for all models) ---
pars = camb.CAMBparams()
pars.set_cosmology(H0=100*h, ombh2=0.0486*h**2, omch2=(om-0.0486)*h**2, mnu=0.0)
pars.InitPower.set_params(ns=0.9667, As=2.142e-9)
pars.set_matter_power(redshifts=[zs], kmax=50.0)
res = camb.get_results(pars)
kh, _, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=50.0, npoints=600)
pk = pk[0]  # (Mpc/h)^3 at z=49

# growth rate f and H at astart for velocities
Ez = np.sqrt(om/astart**3 + ol)
Hz = 100*h*Ez              # km/s/Mpc
Oma = om/astart**3/Ez**2
f_growth = Oma**0.55

# --- Gaussian field ---
rng = np.random.default_rng(seed)
kf = 2*np.pi/Lbox_hMpc     # h/Mpc
kx = np.fft.fftfreq(N, d=1.0/N)*kf
kz = np.fft.rfftfreq(N, d=1.0/N)*kf
KX, KY, KZ = np.meshgrid(kx, kx, kz, indexing='ij')
K2 = KX**2+KY**2+KZ**2
K = np.sqrt(K2); K[0,0,0]=1e-10

Pk_grid = np.interp(K.ravel(), kh, pk).reshape(K.shape)
# delta_k with correct normalization for discrete FFT: Var = P(k)/V
V = Lbox_hMpc**3
amp = np.sqrt(Pk_grid*V)/V*N**3   # so that irfftn gives delta(x) with right sigma
re = rng.standard_normal(K.shape); im = rng.standard_normal(K.shape)
dk = amp*(re+1j*im)/np.sqrt(2)
dk[0,0,0]=0.0
# displacement psi_k = i k delta_k / k^2  (comoving Mpc/h)
def disp(comp):
    psik = 1j*comp*dk/np.where(K2>0,K2,1)
    return np.fft.irfftn(psik, s=(N,N,N), axes=(0,1,2))
psix, psiy, psiz = disp(KX), disp(KY), disp(KZ)

sig = np.std(np.fft.irfftn(dk, s=(N,N,N)))
print(f"delta rms at z={zs}: {sig:.4f}  max disp: {np.max(np.abs(psix)):.3f} Mpc/h")

# velocities: v = a H f psi  (proper km/s; psi comoving Mpc = psi_hMpc/h)
vfac = astart*Hz*f_growth/h
vx, vy, vz = (vfac*psix).astype('f4'), (vfac*psiy).astype('f4'), (vfac*psiz).astype('f4')

outdir = sys.argv[1] if len(sys.argv)>1 else 'ics_64'
os.makedirs(outdir, exist_ok=True)

def write_grafic(fn, arr):
    with open(fn,'wb') as f:
        hdr = struct.pack('<3i8f', N,N,N, dx_Mpc, 0.0,0.0,0.0, astart, om, ol, 100*h)
        f.write(struct.pack('<i',len(hdr))); f.write(hdr); f.write(struct.pack('<i',len(hdr)))
        for i3 in range(N):
            plane = arr[:,:,i3].astype('<f4').tobytes(order='F')
            f.write(struct.pack('<i',len(plane))); f.write(plane); f.write(struct.pack('<i',len(plane)))

write_grafic(os.path.join(outdir,'ic_velcx'), vx)
write_grafic(os.path.join(outdir,'ic_velcy'), vy)
write_grafic(os.path.join(outdir,'ic_velcz'), vz)
print("wrote", outdir, " dxini(Mpc)=", dx_Mpc, " astart=", astart)

# header/density file needed by init_cosmo even for dmonly
delta = np.fft.irfftn(dk, s=(N,N,N), axes=(0,1,2)).astype('f4')
write_grafic(os.path.join(outdir,'ic_deltab'), delta)
print("wrote ic_deltab")

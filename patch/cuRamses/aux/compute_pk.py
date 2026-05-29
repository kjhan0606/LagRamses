#!/usr/bin/env python3
"""
compute_pk.py — Post-processing matter power spectrum P(k) and correlation
function xi(r) from RAMSES binary output data.

Deposits all matter (DM particles, star particles, gas AMR leaf cells, sinks)
onto a uniform grid at levelmin resolution, then computes P(k) via FFT and
xi(r) via Hankel transform with Gaussian window.

Output format matches the on-the-fly RAMSES dump (5 columns):
    k [h/Mpc]    P(k) [(Mpc/h)^3]    r [Mpc/h]    xi(r)    Nmodes

Usage:
    python3 compute_pk.py output_00010/
    python3 compute_pk.py output_00010/ output_00020/ --ngrid 256
    python3 compute_pk.py output_00010/ --no-gas        # particles only
    python3 compute_pk.py output_00010/ --no-particles  # gas only
"""

import os
import sys
import glob
import argparse
import time
import numpy as np


# ============================================================================
# Fortran unformatted binary reader
# ============================================================================
class FortranReader:
    """Read Fortran unformatted sequential (record-based) binary files."""

    def __init__(self, filename):
        self.f = open(filename, 'rb')
        self.filename = filename

    def _read_record_raw(self):
        head = np.fromfile(self.f, dtype=np.int32, count=1)
        if head.size == 0:
            raise EOFError(f"Unexpected EOF in {self.filename}")
        reclen = head[0]
        data = self.f.read(reclen)
        tail = np.fromfile(self.f, dtype=np.int32, count=1)[0]
        assert reclen == tail, f"Record mismatch {reclen} vs {tail} in {self.filename}"
        return data

    def read_ints(self):
        return np.frombuffer(self._read_record_raw(), dtype=np.int32).copy()

    def read_longs(self):
        return np.frombuffer(self._read_record_raw(), dtype=np.int64).copy()

    def read_reals(self):
        return np.frombuffer(self._read_record_raw(), dtype=np.float64).copy()

    def read_string(self):
        return self._read_record_raw().decode('ascii', errors='ignore').strip()

    def read_auto_int(self):
        """Read a record that could be int32 or int64 (auto-detect from size)."""
        data = self._read_record_raw()
        if len(data) % 8 == 0 and len(data) % 4 == 0:
            # Ambiguous — prefer int64 for single-element records (nstar_tot)
            if len(data) == 8:
                return np.frombuffer(data, dtype=np.int64).copy()
            elif len(data) == 4:
                return np.frombuffer(data, dtype=np.int32).copy()
        if len(data) % 8 == 0:
            return np.frombuffer(data, dtype=np.int64).copy()
        return np.frombuffer(data, dtype=np.int32).copy()

    def skip(self, n=1):
        for _ in range(n):
            head = np.fromfile(self.f, dtype=np.int32, count=1)
            if head.size == 0:
                raise EOFError(f"Unexpected EOF during skip in {self.filename}")
            self.f.seek(head[0], 1)
            np.fromfile(self.f, dtype=np.int32, count=1)

    def close(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


# ============================================================================
# Info file reader
# ============================================================================
def read_info(output_dir):
    """Read RAMSES info_NNNNN.txt."""
    info_files = sorted(glob.glob(os.path.join(output_dir, 'info_*.txt')))
    if not info_files:
        raise FileNotFoundError(f"No info_*.txt in {output_dir}")

    info = {}
    with open(info_files[0]) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('DOMAIN'):
                continue
            if '=' in line:
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    if key == 'ordering type':
                        info['ordering'] = val
                        continue
                    try:
                        if '.' in val or 'E' in val.upper():
                            info[key] = float(val)
                        else:
                            info[key] = int(val)
                    except ValueError:
                        info[key] = val
    return info


def get_nchar(output_dir):
    """Extract NNNNN string from output_NNNNN/."""
    return os.path.basename(output_dir.rstrip('/')).split('_')[-1]


# ============================================================================
# Particle reader
# ============================================================================
def read_particles(output_dir, info):
    """Read particle positions and masses from all CPU binary files.

    Returns: (positions [N,3], masses [N]) in code units.
    """
    ncpu = int(info['ncpu'])
    nchar = get_nchar(output_dir)
    ndim = int(info.get('ndim', 3))

    all_pos = []
    all_mass = []

    for icpu in range(1, ncpu + 1):
        fname = os.path.join(output_dir, f'part_{nchar}.out{icpu:05d}')
        if not os.path.exists(fname):
            continue

        with FortranReader(fname) as fr:
            fr.skip(1)             # ncpu
            fr.skip(1)             # ndim
            npart = fr.read_ints()[0]
            fr.skip(1)             # localseed
            fr.skip(1)             # nstar_tot (int64)
            fr.skip(2)             # mstar_tot, mstar_lost
            fr.skip(1)             # nsink

            if npart == 0:
                continue

            # Positions
            pos = np.zeros((npart, ndim))
            for idim in range(ndim):
                pos[:, idim] = fr.read_reals()

            # Skip velocities
            fr.skip(ndim)

            # Masses
            mass = fr.read_reals()

            all_pos.append(pos)
            all_mass.append(mass)

    if all_pos:
        return np.vstack(all_pos), np.concatenate(all_mass)
    return np.empty((0, 3)), np.empty(0)


# ============================================================================
# AMR + Hydro reader (gas leaf cells)
# ============================================================================
def read_amr_header(amr):
    """Read AMR file header, return dict of values."""
    hdr = {}
    hdr['ncpu'] = amr.read_ints()[0]
    hdr['ndim'] = amr.read_ints()[0]
    nxyz = amr.read_ints()
    hdr['nx'], hdr['ny'], hdr['nz'] = nxyz[0], nxyz[1], nxyz[2]
    hdr['nlevelmax'] = amr.read_ints()[0]
    hdr['ngridmax'] = amr.read_ints()[0]
    hdr['nboundary'] = amr.read_ints()[0]
    amr.skip(1)  # ngrid_current
    amr.skip(1)  # boxlen

    # Time variables
    amr.skip(1)  # noutput, iout, ifout
    amr.skip(2)  # tout, aout
    amr.skip(1)  # t
    amr.skip(2)  # dtold, dtnew
    amr.skip(1)  # nstep, nstep_coarse
    amr.skip(1)  # const, mass_tot_0, rho_tot
    cosmo = amr.read_reals()  # omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini
    hdr['boxlen_ini'] = cosmo[6]
    hdr['h0'] = cosmo[4]
    amr.skip(1)  # aexp, hexp, etc
    amr.skip(1)  # mass_sph

    # Level lists
    ncpu_amr = hdr['ncpu']
    nlevelmax = hdr['nlevelmax']
    amr.skip(1)  # headl
    amr.skip(1)  # taill
    raw = amr.read_ints()  # numbl(1:ncpu, 1:nlevelmax)
    hdr['numbl'] = raw.reshape((ncpu_amr, nlevelmax), order='F')
    amr.skip(1)  # numbtot

    nboundary = hdr['nboundary']
    if nboundary > 0:
        amr.skip(1)  # headb
        amr.skip(1)  # tailb
        raw = amr.read_ints()
        hdr['numbb'] = raw.reshape((nboundary, nlevelmax), order='F')
    else:
        hdr['numbb'] = np.zeros((0, nlevelmax), dtype=np.int32)

    # Free memory
    amr.skip(1)

    # Ordering
    ordering = amr.read_string()
    hdr['ordering'] = ordering
    if 'bisection' in ordering:
        amr.skip(5)
    elif 'ksection' in ordering:
        amr.skip(10)  # 3 header ints + 7 data arrays
    else:  # hilbert
        amr.skip(1)  # bound_key

    # Coarse level
    ncoarse = hdr['nx'] * hdr['ny'] * hdr['nz']
    amr.skip(3)  # son, flag1, cpu_map at coarse level

    return hdr


def read_gas_and_deposit(output_dir, info, ngrid):
    """Read AMR + hydro, deposit leaf cell gas mass onto uniform grid.

    Returns: mass_grid [ngrid^3] in code units.
    """
    ncpu = int(info['ncpu'])
    nchar = get_nchar(output_dir)
    ndim = int(info.get('ndim', 3))
    boxlen = info['boxlen']  # typically 1.0

    mass_grid = np.zeros(ngrid ** 3, dtype=np.float64)
    twotondim = 2 ** ndim     # 8
    twondim = 2 * ndim        # 6
    n_leaf_total = 0

    for icpu in range(1, ncpu + 1):
        amr_fname = os.path.join(output_dir, f'amr_{nchar}.out{icpu:05d}')
        hyd_fname = os.path.join(output_dir, f'hydro_{nchar}.out{icpu:05d}')
        if not os.path.exists(amr_fname) or not os.path.exists(hyd_fname):
            continue

        with FortranReader(amr_fname) as amr, FortranReader(hyd_fname) as hyd:
            # Read AMR header (positions file pointer past coarse level)
            hdr = read_amr_header(amr)
            nlevelmax = hdr['nlevelmax']
            nboundary = hdr['nboundary']
            ncpu_amr = hdr['ncpu']
            nx = hdr['nx']
            numbl = hdr['numbl']
            numbb = hdr['numbb']

            # Hydro header
            hyd.skip(1)  # ncpu
            nvar = hyd.read_ints()[0]
            hyd.skip(1)  # ndim
            hyd.skip(1)  # nlevelmax
            hyd.skip(1)  # nboundary
            hyd.skip(1)  # gamma

            # nvar records per child cell in hydro output
            nrec_per_child = nvar

            # ---- Level loop ----
            for ilevel in range(1, nlevelmax + 1):
                dx = boxlen / (nx * 2 ** ilevel)
                cell_vol = dx ** 3

                for ibound in range(1, ncpu_amr + nboundary + 1):
                    if ibound <= ncpu_amr:
                        ncache = numbl[ibound - 1, ilevel - 1]
                    else:
                        ncache = numbb[ibound - ncpu_amr - 1, ilevel - 1]

                    # Hydro always writes ilevel + ncache
                    hyd.skip(2)

                    if ncache > 0:
                        # === AMR: read grid positions + son, skip rest ===
                        amr.skip(3)  # ind_grid, next, prev

                        xg = np.zeros((ncache, ndim))
                        for idim in range(ndim):
                            xg[:, idim] = amr.read_reals()

                        amr.skip(1)        # father
                        amr.skip(twondim)  # nbor (6 records)

                        # Son indices for each child
                        son = np.zeros((ncache, twotondim), dtype=np.int32)
                        for ind in range(twotondim):
                            son[:, ind] = amr.read_ints()

                        amr.skip(twotondim)  # cpu_map
                        amr.skip(twotondim)  # flag1

                        # Only process cells owned by this CPU
                        is_owner = (ibound == icpu)

                        # === Hydro: read density for each child cell ===
                        for ind in range(twotondim):
                            density = hyd.read_reals()
                            hyd.skip(nrec_per_child - 1)  # skip vx, vy, vz, P, scalars

                            if not is_owner:
                                continue

                            # Cell positions
                            iz_c = ind // 4
                            iy_c = (ind - 4 * iz_c) // 2
                            ix_c = ind - 2 * iy_c - 4 * iz_c
                            xc = xg[:, 0] + (ix_c - 0.5) * dx
                            yc = xg[:, 1] + (iy_c - 0.5) * dx
                            zc = xg[:, 2] + (iz_c - 0.5) * dx

                            # Leaf mask: son == 0
                            leaf = (son[:, ind] == 0)
                            if not np.any(leaf):
                                continue
                            n_leaf_total += np.sum(leaf)

                            # Gas mass = rho * cell_volume
                            gmass = density[leaf] * cell_vol

                            # NGP deposit
                            ig = np.clip((xc[leaf] / boxlen * ngrid).astype(int),
                                         0, ngrid - 1)
                            jg = np.clip((yc[leaf] / boxlen * ngrid).astype(int),
                                         0, ngrid - 1)
                            kg = np.clip((zc[leaf] / boxlen * ngrid).astype(int),
                                         0, ngrid - 1)
                            idx = ig * ngrid * ngrid + jg * ngrid + kg
                            np.add.at(mass_grid, idx, gmass)
                    # ncache == 0: nothing in AMR, hydro already skipped

        if icpu % max(1, ncpu // 10) == 0 or icpu == ncpu:
            sys.stdout.write(f'\r  Gas: reading CPU {icpu}/{ncpu} ...')
            sys.stdout.flush()

    print(f'\r  Gas: {n_leaf_total} leaf cells from {ncpu} CPUs      ')
    return mass_grid


# ============================================================================
# NGP deposit for particles
# ============================================================================
def deposit_particles(positions, masses, ngrid, boxlen=1.0):
    """NGP mass deposit onto uniform grid. Returns flat array [ngrid^3]."""
    grid = np.zeros(ngrid ** 3, dtype=np.float64)
    if len(masses) == 0:
        return grid
    ix = np.clip((positions[:, 0] / boxlen * ngrid).astype(int), 0, ngrid - 1)
    iy = np.clip((positions[:, 1] / boxlen * ngrid).astype(int), 0, ngrid - 1)
    iz = np.clip((positions[:, 2] / boxlen * ngrid).astype(int), 0, ngrid - 1)
    idx = ix * ngrid * ngrid + iy * ngrid + iz
    np.add.at(grid, idx, masses)
    return grid


# ============================================================================
# Power spectrum and correlation function
# ============================================================================
def compute_pk_xi(mass_grid, ngrid, boxlen_phys, npart_total, nkbin=50):
    """Compute P(k) and xi(r) from deposited mass grid.

    Args:
        mass_grid: flat array of total mass per cell [ngrid^3]
        ngrid: grid resolution per side
        boxlen_phys: box size in Mpc/h
        npart_total: total particle count (for shot noise)
        nkbin: number of logarithmic bins

    Returns: dict with k, Pk, r, xi, counts, shot_noise, k_damp
    """
    twopi = 2.0 * np.pi

    # Overdensity: delta = M(x) / <M> - 1
    mean_mass = np.sum(mass_grid) / ngrid ** 3
    if mean_mass <= 0:
        raise ValueError("Total mass is zero — no data deposited")
    delta = mass_grid.reshape((ngrid, ngrid, ngrid)) / mean_mass - 1.0

    print(f"  FFT: delta mean = {np.mean(delta):.6e}, rms = {np.std(delta):.6e}")

    # R2C FFT
    delta_k = np.fft.rfftn(delta)

    # Physical wavenumbers
    k_fund = twopi / boxlen_phys  # h/Mpc
    k_nyq = k_fund * ngrid / 2.0

    # k-magnitude for each FFT mode
    kx_arr = np.fft.fftfreq(ngrid) * ngrid  # integer wavenumbers
    ky_arr = np.fft.fftfreq(ngrid) * ngrid
    kz_arr = np.arange(ngrid // 2 + 1)

    # Volume and normalization: P(k) = V/N^2 * |delta_k|^2
    V_box = boxlen_phys ** 3
    N_total = ngrid ** 3
    norm_fft = V_box / N_total ** 2

    # Weight for R2C symmetry
    # kz = 0 and kz = Nz/2 → weight 1, else weight 2
    wt_kz = np.ones(ngrid // 2 + 1)
    wt_kz[1:ngrid // 2] = 2.0

    # Logarithmic k-bins
    k_min_bin = k_fund
    k_max_bin = k_nyq
    dk_log = (np.log10(k_max_bin) - np.log10(k_min_bin)) / nkbin
    k_edges = 10.0 ** (np.log10(k_min_bin) + np.arange(nkbin + 1) * dk_log)
    k_centers = np.sqrt(k_edges[:-1] * k_edges[1:])

    pk_sum = np.zeros(nkbin)
    pk_cnt = np.zeros(nkbin)

    # Bin the power spectrum
    print("  Binning P(k) ...")
    power_spec = np.abs(delta_k) ** 2 * norm_fft

    for ikx in range(ngrid):
        kx = kx_arr[ikx]
        for iky in range(ngrid):
            ky = ky_arr[iky]
            for ikz in range(ngrid // 2 + 1):
                kz = kz_arr[ikz]
                if ikx == 0 and iky == 0 and ikz == 0:
                    continue  # skip DC
                k_mag = k_fund * np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
                if k_mag < k_min_bin or k_mag >= k_max_bin:
                    continue
                ibin = int((np.log10(k_mag) - np.log10(k_min_bin)) / dk_log)
                if ibin < 0 or ibin >= nkbin:
                    continue
                w = wt_kz[ikz]
                pk_sum[ibin] += power_spec[ikx, iky, ikz] * w
                pk_cnt[ibin] += w

    # Average P(k)
    pk_avg = np.zeros(nkbin)
    valid = pk_cnt > 0.5
    pk_avg[valid] = pk_sum[valid] / pk_cnt[valid]

    # Shot noise
    shot_noise = V_box / npart_total if npart_total > 0 else 0.0

    # --- xi(r) via Hankel transform with Gaussian window ---
    r_min = boxlen_phys / ngrid        # cell size
    r_max = boxlen_phys / 2.0          # half box
    dr_log = (np.log10(r_max) - np.log10(r_min)) / nkbin
    r_edges = 10.0 ** (np.log10(r_min) + np.arange(nkbin + 1) * dr_log)
    r_centers = np.sqrt(r_edges[:-1] * r_edges[1:])

    k_damp = k_nyq
    xi_vals = np.zeros(nkbin)

    print("  Computing xi(r) via Hankel transform ...")
    for j in range(nkbin):
        r = r_centers[j]
        xi = 0.0
        for i in range(nkbin):
            if pk_cnt[i] < 0.5:
                continue
            pk = pk_avg[i] - shot_noise
            dk = k_edges[i + 1] - k_edges[i]
            gauss_w = np.exp(-(k_centers[i] / k_damp) ** 2)
            kr = k_centers[i] * r
            sinc_kr = np.sin(kr) / kr
            xi += pk * gauss_w * k_centers[i] ** 2 * sinc_kr * dk
        xi_vals[j] = xi * 2.0 / twopi ** 2   # 1/(2*pi^2)

    return {
        'k': k_centers, 'Pk': pk_avg, 'counts': pk_cnt,
        'r': r_centers, 'xi': xi_vals,
        'shot_noise': shot_noise, 'k_damp': k_damp,
        'k_edges': k_edges, 'r_edges': r_edges,
    }


# ============================================================================
# Write output
# ============================================================================
def write_pk_file(filename, aexp, boxlen_phys, ngrid, npart, result):
    """Write P(k) + xi(r) ASCII file matching Fortran format."""
    nkbin = len(result['k'])
    with open(filename, 'w') as f:
        f.write(f"# Power spectrum at a_exp =  {aexp:15.8E}\n")
        f.write(f"# boxlen (Mpc/h) =  {boxlen_phys:15.8E}\n")
        f.write(f"# N_grid = {ngrid:10d}\n")
        f.write(f"# N_part = {npart:15d}\n")
        f.write(f"# shot_noise (Mpc/h)^3 =  {result['shot_noise']:15.8E}\n")
        f.write(f"# k_damp (h/Mpc) =  {result['k_damp']:15.8E}\n")
        f.write("# k [h/Mpc]    P(k) [(Mpc/h)^3]    r [Mpc/h]    xi(r)    Nmodes\n")

        for i in range(nkbin):
            if result['counts'][i] > 0.5:
                f.write(f"  {result['k'][i]:14.8E}  {result['Pk'][i]:14.8E}"
                        f"  {result['r'][i]:14.8E}  {result['xi'][i]:14.8E}"
                        f"  {result['counts'][i]:14.8E}\n")

    print(f"  Written: {filename}")


# ============================================================================
# Process one output directory
# ============================================================================
def process_output(output_dir, ngrid=0, nkbin=50,
                   include_particles=True, include_gas=True):
    """Full pipeline: read data → deposit → FFT → P(k) → xi(r) → write."""
    output_dir = output_dir.rstrip('/')
    nchar = get_nchar(output_dir)
    print(f"\n{'='*60}")
    print(f" Processing {output_dir}")
    print(f"{'='*60}")

    t0 = time.time()

    # Read info
    info = read_info(output_dir)
    ncpu = int(info['ncpu'])
    ndim = int(info.get('ndim', 3))
    levelmin = int(info['levelmin'])
    boxlen = info['boxlen']
    aexp = info['aexp']
    H0 = info.get('H0', info.get('h0', 67.66))  # km/s/Mpc

    # Determine grid resolution
    nx = 1  # typical for cosmological RAMSES
    if ngrid <= 0:
        ngrid = nx * 2 ** levelmin
    print(f"  levelmin={levelmin}, aexp={aexp:.6f}, z={1.0/aexp - 1:.2f}")
    print(f"  Grid: {ngrid}^3 = {ngrid**3} cells")

    # Get boxlen_ini from AMR file
    amr_fname = os.path.join(output_dir, f'amr_{nchar}.out00001')
    with FortranReader(amr_fname) as amr:
        hdr = read_amr_header(amr)
    boxlen_ini = hdr['boxlen_ini']
    nx = hdr['nx']
    print(f"  boxlen_ini = {boxlen_ini:.2f} Mpc/h, nx={nx}")

    # Recompute ngrid if auto
    if ngrid == 2 ** levelmin:
        ngrid = nx * 2 ** levelmin
        print(f"  Adjusted grid: {ngrid}^3 (nx={nx})")

    # Initialize total mass grid
    total_mass = np.zeros(ngrid ** 3, dtype=np.float64)
    npart_total = 0

    # --- Particles (DM + stars + sinks) ---
    if include_particles:
        print("  Reading particles ...")
        pos, mass = read_particles(output_dir, info)
        npart_total = len(mass)
        print(f"  Particles: {npart_total}, M_total = {np.sum(mass):.6e}")
        if npart_total > 0:
            total_mass += deposit_particles(pos, mass, ngrid, boxlen)
            del pos, mass

    # --- Gas (AMR leaf cells) ---
    if include_gas:
        gas_mass = read_gas_and_deposit(output_dir, info, ngrid)
        gas_total = np.sum(gas_mass)
        print(f"  Gas: M_total = {gas_total:.6e}")
        total_mass += gas_mass
        del gas_mass

    print(f"  Total deposited mass = {np.sum(total_mass):.6e}")

    # --- Compute P(k) and xi(r) ---
    result = compute_pk_xi(total_mass, ngrid, boxlen_ini, npart_total, nkbin)
    del total_mass

    # --- Write output ---
    pk_filename = os.path.join(output_dir, f'pk_{nchar}.dat')
    write_pk_file(pk_filename, aexp, boxlen_ini, ngrid, npart_total, result)

    dt = time.time() - t0
    print(f"  Done in {dt:.1f} s")
    return pk_filename


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description='Compute matter P(k) and xi(r) from RAMSES output')
    parser.add_argument('output_dirs', nargs='+',
                        help='RAMSES output directories (e.g., output_00010/)')
    parser.add_argument('--ngrid', type=int, default=0,
                        help='Grid resolution (default: nx*2^levelmin)')
    parser.add_argument('--nkbin', type=int, default=50,
                        help='Number of k-bins (default: 50)')
    parser.add_argument('--no-gas', action='store_true',
                        help='Exclude gas (particles only)')
    parser.add_argument('--no-particles', action='store_true',
                        help='Exclude particles (gas only)')

    args = parser.parse_args()

    for output_dir in args.output_dirs:
        if not os.path.isdir(output_dir):
            print(f"WARNING: {output_dir} not found, skipping")
            continue
        try:
            process_output(
                output_dir,
                ngrid=args.ngrid,
                nkbin=args.nkbin,
                include_particles=not args.no_particles,
                include_gas=not args.no_gas,
            )
        except Exception as e:
            print(f"ERROR processing {output_dir}: {e}")
            import traceback
            traceback.print_exc()


if __name__ == '__main__':
    main()

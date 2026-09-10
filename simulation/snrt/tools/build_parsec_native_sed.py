#!/usr/bin/env python3
"""Offline PARSEC high-mass common-population Q/E; NOT a full/NLTE SSP.

Measured cumulative-above-edge Q constrains five disjoint photon intervals.
Within them the shape is a Planck prior. Missing end times use the actual
full track L/Teff blackbody, explicitly distinct from measured Q. RAMSES
reads only positive cumulative moments and uses its own feedback IMF cells.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np
from build_parsec_pair_feedback import (PINNED, EXTRA_TRACKS, LOW_TRACKS, Z_GRIDS,
                                      PARSEC_LSUN, PRECISION_MODEL, source_name, compress)

ROOT = Path(__file__).resolve().parents[1]
PHOTONS = {
    'Z0.008_Y0.263_photons.zip': 'f6c78dc93e9aa9118dca8733983251e47040dd97cb3afa55230a6f981029a183',
    'Z0.014_Y0.273_photons.zip': '3f0caaa8f29fd7707f1538ca9fc8694bb0b5ad475db6350a789db59d4f0ec8d6',
}
EXTRA_PHOTONS = {
    'Z0.017_Y0.279_photons.zip': 'e38befb640ab2f78e75419b8acd580ee226b91be1740bab04e43388110bc51cc',
    'Z0.02_Y0.284_photons.zip': 'c0e7171fba4de5b6a3b7bccdfe2b31ced926021716e9c59d2088f47e5e6db7e8',
    'Z0.03_Y0.302_photons.zip': '4e99bb71e0043e3b7cb5698c570ed8d26af238ea889cc631f4a18085fe93676c',
}
LOW_PHOTONS = dict(zip((
    'Z1E-11_Y0.2485_photons.zip', 'Z0.0001_Y0.2487_photons.zip',
    'Z0.002_Y0.252_photons.zip', 'Z0.004_Y0.256_photons.zip',
    'Z0.006_Y0.259_photons.zip', 'Z0.01_Y0.267_photons.zip'), (
    '121923f8457f106b1fb216444e9dc3e8c7d2997bda9d83602ed80d572e13878f',
    '4fe530cbb84a55bf9201b8b9c8ffe2478683221232692fceec701fe937512883',
    '97de7c01c9114910815894a45c178aa6fa68b3da749c0bdc933c4f48d3b97196',
    '845e268e96570eeb55e2f8f92bd5c3d3d4476f3531e0f24c9f36e0b9194ad14b',
    'bdaf060fb592e6161f5aa07fbd036125666b9c59cfd3d8e4e2b586f358bdceac',
    'a20cd865960de3a99990ae7b39a8ffecab178c1e496be73c97974927d63656c5')))
Q_EDGES = np.array([11.2, 13.6, 24.6, 35.12, 54.4, np.inf])
EV = 1.602176634e-12
KB = 8.617333262145e-5
LSUN = PARSEC_LSUN
MYR = 31557600e6


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Spectrum:
    def __init__(self, order=64):
        self.x, self.w = np.polynomial.legendre.leggauss(order)
        self.lx, self.lw = np.polynomial.laguerre.laggauss(order)

    def moment(self, kt, low, high, reference=0.):
        """Positive scaled Planck photon/energy integrals, no Wien overflow.

        reference is the ORIGINAL constraint interval's lower edge; use the
        same scaling for its denominator and every native subinterval.
        Laguerre integrates the unbounded tail rather than truncating it.
        """
        if high <= low:
            return np.zeros(len(kt)), np.zeros(len(kt))
        if np.isinf(high):
            e = low + kt[:, None]*self.lx
            w = kt[:, None]*self.lw*np.exp(-(low-reference)/kt[:, None])
        else:
            e = np.broadcast_to(low+(self.x+1)*.5*(high-low), (len(kt), len(self.x)))
            w = self.w*.5*(high-low)*np.exp(-(e-reference)/kt[:, None])
        f = e**2/(-np.expm1(-e/kt[:, None]))
        return np.sum(w*f, axis=1), np.sum(w*f*e, axis=1)

    def planck(self, logl, logt, edges):
        kt = KB*10**logt
        lum = 10**logl*LSUN/EV
        scale = lum/(np.pi**4*kt**4/15)
        result = np.zeros((len(kt), 21))
        for g,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
            n,e = self.moment(kt,lo,hi)
            result[:,g] = scale*n;result[:,g+9] = scale*e
        result[:,18] = lum
        result[:,19] = scale*self.moment(kt,0.,edges[0])[1]
        result[:,20] = scale*self.moment(kt,edges[-1],np.inf)[1]
        return result

    def constrained(self, data, edges):
        kt = KB*10**data[:,3]
        lum = 10**data[:,2]*LSUN/EV
        q = 10**data[:,[9,5,6,8,7]]
        disjoint = np.column_stack((q[:,:-1]-q[:,1:],q[:,-1]))
        if not np.isfinite(disjoint).all() or np.any(disjoint<0):
            raise ValueError('nonmonotone/nonfinite measured Q; no clipping or inferred sentinel')
        result = np.zeros((len(kt),21))
        high_energy = np.zeros(len(kt))
        for j,(lo,hi) in enumerate(zip(Q_EDGES[:-1],Q_EDGES[1:])):
            den,energy = self.moment(kt,lo,hi,lo)
            if np.any(den<=0):
                raise ValueError('unresolved positive Planck constraint denominator')
            scale = disjoint[:,j]/den
            high_energy += scale*energy
            for g,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
                n,e = self.moment(kt,max(a,lo),min(b,hi),lo)
                result[:,g] += scale*n;result[:,g+9] += scale*e
            if np.isinf(hi):
                result[:,20] = scale*self.moment(kt,edges[-1],np.inf,lo)[1]
        residual = lum-high_energy
        if np.any(residual<0):
            raise ValueError('Q-constrained Planck spectrum exceeds actual Lbol; no fallback')
        den = self.moment(kt,0.,Q_EDGES[0])[1]
        scale = residual/den
        for g,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
            n,e = self.moment(kt,lo,min(hi,Q_EDGES[0]))
            result[:,g] += scale*n;result[:,g+9] += scale*e
        result[:,18] = lum
        result[:,19] = scale*self.moment(kt,0.,edges[0])[1]
        return result, high_energy/lum


def check_rates(values, edges):
    if not np.isfinite(values).all() or np.any(values<0):
        raise ValueError('nonfinite/negative reconstructed moments')
    if (np.any(values[:,9:18] < values[:,:9]*edges[:-1]*(1-1e-10)) or
            np.any(values[:,9:18] > values[:,:9]*edges[1:]*(1+1e-10))):
        raise ValueError('photon energy outside native group')
    residual = abs(values[:,9:18].sum(axis=1)+values[:,19:].sum(axis=1)-values[:,18])
    error = float(np.max(residual/values[:,18]))
    if error > 1e-10:
        raise ValueError(f'bolometric quadrature error {error}')
    return error


def intervals(age, rates):
    return np.vstack((np.zeros(21),np.diff(age)[:,None]*(rates[1:]+rates[:-1])*.5*MYR))


def build(args):
    feedback = json.loads((args.feedback/'manifest.json').read_text())
    grid = feedback.get('metallicity_grid', 'solar_pair')
    if grid not in Z_GRIDS:
        raise ValueError('unsupported matched feedback metallicity grid')
    coordinates = Z_GRIDS[grid]
    precision = grid == 'precision_eleven'
    if precision and feedback['model'] != PRECISION_MODEL:
        raise ValueError('precision grid requires its own baryonic/RATE model identity')
    z_nodes = [float(z) for z,_ in coordinates]
    inputs = {**PINNED, **PHOTONS}
    if grid != 'solar_pair':
        inputs.update(EXTRA_TRACKS);inputs.update(EXTRA_PHOTONS)
    if precision:
        inputs.update(LOW_TRACKS);inputs.update(LOW_PHOTONS)
    for name,sha in inputs.items():
        if digest(args.source_dir/name)!=sha:
            raise ValueError(f'physical source checksum mismatch: {name}')
    for name,sha in feedback['output_sha256'].items():
        if digest(args.feedback/name)!=sha:
            raise ValueError(f'feedback package changed: {name}')
    if feedback['source_domain'] != {'mass':[14,600],'z':z_nodes,'rotation':0}:
        raise ValueError('not the supported matched feedback population')
    if (len(feedback['nodes']) != 45*len(z_nodes) or
            sorted({n['z'] for n in feedback['nodes']}) != z_nodes):
        raise ValueError('incomplete matched feedback grid')
    if not 0<=args.escape_fraction<=1 or not np.isfinite(args.escape_fraction):
        raise ValueError('invalid escaped fraction')
    edges = np.loadtxt(args.group_edges)
    ledger = json.loads(args.transport_ledger.read_text())
    if (edges.shape!=(10,) or not np.array_equal(edges,ledger['group_edges_ev']) or
            digest(args.group_edges)!=ledger['group_edges_sha256']):
        raise ValueError('native/ledger edge mismatch')
    spectrum = Spectrum(args.order)
    nodes=[];stats=[];nr=0;bolerr=0.;highmax=0.
    for zs,ys in coordinates:
        with zipfile.ZipFile(args.source_dir/source_name(zs,ys,'photons')) as photons, \
             zipfile.ZipFile(args.source_dir/source_name(zs,ys,'tracks')) as tracks:
            for node in feedback['nodes']:
                if node['z']!=float(zs):continue
                m=node['mass'];death=node['age_yr']/1e6
                data=np.loadtxt(io.BytesIO(photons.read(source_name(zs,ys,'photons',m))),skiprows=3)
                track=np.loadtxt(io.BytesIO(tracks.read(source_name(zs,ys,'tracks',m))),skiprows=3)
                if not np.isfinite(data).all() or np.any(np.diff(data[:,0])<0):
                    raise ValueError('invalid photon table')
                data=data[np.r_[np.diff(data[:,0])>0,True]]
                if (data[0,0]<=0 or data[0,0]>=death or
                        (not precision and data[-1,0]>=death) or
                        abs(track[-1,1]/1e6-death)>1e-12*death):
                    raise ValueError('unexpected source timing; no unreviewed extrapolation')
                rates,ratio=spectrum.constrained(data,edges)
                highmax=max(highmax,float(ratio.max()))
                bolerr=max(bolerr,check_rates(rates,edges))
                measured_end = float(data[-1,0])
                measured_age = data[:,0]
                if measured_end > death:
                    # Restrict the existing positive piecewise-linear RATE
                    # integral to the material lifetime, without logQ or Teff
                    # extrapolation. This preserves the same group energy and
                    # bolometric convex constraints at the inserted endpoint.
                    j = int(np.searchsorted(measured_age,death))
                    w = (death-measured_age[j-1])/(measured_age[j]-measured_age[j-1])
                    endpoint = (1-w)*rates[j-1]+w*rates[j]
                    rates = np.vstack((rates[:j],endpoint))
                    measured_age = np.r_[measured_age[:j],death]
                    bolerr=max(bolerr,check_rates(rates,edges))
                age=np.r_[0.,measured_age]
                release=intervals(age,np.vstack((rates[0],rates)))
                values=np.cumsum(release,axis=0)
                # The missing late atmosphere is a distinct full-track Planck
                # closure, not a hold of the last measured five Q constraints.
                t=track[:,1]/1e6
                tail_age=np.r_[age[-1],t[t>age[-1]]]
                tail_age=np.unique(tail_age);tail_age[-1]=death
                lt=np.interp(tail_age,t,track[:,4]);ll=np.interp(tail_age,t,track[:,3])
                tail_rate=spectrum.planck(ll,lt,edges)
                bolerr=max(bolerr,check_rates(tail_rate,edges))
                tail_release=intervals(tail_age,tail_rate)
                tail=np.cumsum(tail_release,axis=0)
                age=np.r_[age,tail_age[1:]]
                values=np.vstack((values,values[-1]+tail[1:]))
                index,error=compress(age,values,1e-4)
                release=np.vstack((release,tail_release[1:]))
                # Store each compressed interval's positive extensive moments,
                # not the difference of large rounded cumulative endpoints.
                # This preserves tiny late hard-band photons after a bright
                # earlier phase, without clipping a physical spectral tail.
                compact=np.vstack((np.zeros(21),[release[a+1:b+1].sum(axis=0)
                                                 for a,b in zip(index[:-1],index[1:])]))
                age=age[index]
                nodes.append((node,age,compact*args.escape_fraction));nr+=len(age)
                stats.append({'mass':m,'z':float(zs),'source_photon_rows':len(data),'knots':len(age),
                              'first_photon_age_yr':float(data[0,0]*1e6),
                              'tail_duration_yr':float(max(0.,death-measured_end)*1e6),
                              'tail_lifetime_fraction':float(max(0.,1-measured_end/death)),
                              'tail_bolometric_fraction':float(tail[-1,18]/values[-1,18]),
                              'compression_relative_final_component_error':error,
                              **({'measured_overlap_discarded_yr':float(max(0.,measured_end-death)*1e6),
                                  'overlap_policy':'positive_QE_rate_interpolation_at_material_death'} if precision else {})})
    if nr > 500000:
        raise ValueError('matched source exceeds existing native knot capacity')
    lines=['SNRT_PARSEC_QE_INTERVAL_V1',f'{len(nodes)} {nr}',feedback['model']]
    for node,age,values in nodes:
        lines.append(f"{node['mass']:.17e} {node['z']:.17e} {age[-1]:.17e} {node['fate']} {len(age)}")
        lines.extend(' '.join(f'{x:.17e}' for x in (t,*row)) for t,row in zip(age,values))
    payload='\n'.join(lines)+'\n';source_sha=hashlib.sha256(payload.encode()).hexdigest()
    output=args.output.resolve()
    if "'" in str(output):raise ValueError('unsupported quote in NML path')
    nml=f"""! HIGH-MASS ONLY; no <14Msun SED and no full/NLTE SSP claim.
&snrt_stellar_sed
 version=4, na=0, nz=0,
 population_binding='match_feedback_high_mass_only',
 radiation_population='PARSEC_v2_nonrot_high_mass_only',
 source_sha256='{source_sha}',
 node_history_file='{output/'nodes.dat'}',
 imf_id={feedback['imf_id']}, population_id=0, binary_fraction=0,
 imf_min=0.08, imf_max=600,
 status='reference_control', approval_id='',
 interpolation='linear_cumulative_age_linear_Z',
 fraction_semantics='escaped', escape_fraction={args.escape_fraction:.17e},
 young_age_policy='hold_first_to_zero',
 spectral_tail_policy='Q5_Planck_and_track_tail_v1',
 energy_semantics='photon_number_and_energy_v1',
 transport_sha256='{digest(args.transport_ledger)}',
 edges_sha256='{digest(args.group_edges)}',
/
"""
    manifest={'model':'PARSEC_Q5_Planck_and_track_tail_v1','status':'high_mass_only_reference_comparison',
              'full_ssp':False,'exact_atmosphere_spectrum':False,'imf':feedback['imf_id'],
              'imf_denominator':[.08,600],'source_mass':[14,600],'source_z':z_nodes,
              'sources_sha256':inputs,'feedback_outputs_sha256':feedback['output_sha256'],
              'source_sha256':source_sha,'group_edges_ev':edges.tolist(),
              'nominal_Q_edges_ev':Q_EDGES[:-1].tolist(),'quadrature_order':args.order,
              'PARSEC_Lsun_erg_s':LSUN,
              'max_relative_bolometric_error':bolerr,'max_Q_constrained_energy_over_Lbol':highmax,
              'escaped_fraction_applied_once':args.escape_fraction,'young_hold':'first_photon_state_to_zero',
              'missing_tail':'actual_full_track_L_Teff_Planck;not_measured_Q;rate_jump_at_join',
              'mass_weighting':'native_feedback_shared_cells_IMF_mass_fraction_over_source_mass',
              'units':'endAgeMyr; preceding interval photons/star,eV/star; final columns bolometric,IR-tail,X-tail eV',
              'interpolation':'piecewise linear cumulative; positive interval storage avoids endpoint cancellation',
              'compression_tolerance_per_final_component':1e-4,'nodes':stats,
              'generator_sha256':digest(Path(__file__))}
    output.mkdir(parents=True,exist_ok=False)
    (output/'nodes.dat').write_text(payload)
    (output/'source.nml').write_text(nml)
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PARSEC_HIGH_MASS_QE',len(nodes),nr,'bolometric_error',bolerr,'max_high_E/L',highmax,output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--feedback',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--escape-fraction',type=float,required=True)
    p.add_argument('--order',type=int,choices=(32,64,128),default=64)
    p.add_argument('--group-edges',type=Path,default=ROOT/'config/p0_photon_group_edges_ev.txt')
    p.add_argument('--transport-ledger',type=Path,default=ROOT/'data/p4_pilot_agn_photon_ledger.json')
    build(p.parse_args())

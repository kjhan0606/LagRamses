#!/usr/bin/env python3
"""Opt-in mixed evolution comparison; NOT a complete/matched bolometric SSP.

Emit ordinary native yields plus v5 history and v2/v3 positive Q/E intervals.
Uses actual PARSEC2--12 tracks, a named13 proxy, KL16/Fishlock and Sukhbold.
No late radiation plateau; no double counting of terminal-source wind budgets.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile
import numpy as np
from build_parsec_pair_feedback import Z_GRIDS, source_name, read_ejecta, compress
from build_parsec_native_sed import Spectrum, check_rates, intervals, MYR
from generate_snrt_group_profile import read_profile
from read_karakas_lugaro2016 import (read_karakas_lugaro2016, read_fishlock2014,
    DEFAULT_MATRIX, LIFETIME_PATH, LIFETIME_SHA256, FISHLOCK_SHA256)
from build_sukhbold_native_lowmass import build_lowmass
from audit_g2_sukhbold2016_candidate import DEFAULT_ROOT as DEFAULT_SUKHBOLD_ROOT

MODEL = 'parsec_mixed_lowmass_truncated_v1'
ROOT = Path(__file__).resolve().parents[1]
MSUN = 1.98847e33
LOW_MASSES = [round(2+.2*i, 1) for i in range(40)] + [10., 11., 12., 13.]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nearest(items, value, key):
    return min(items, key=lambda r: (abs(key(r)-value), key(r)))


def terminal_sources(imf, sukhbold_root=DEFAULT_SUKHBOLD_ROOT):
    agb = read_karakas_lugaro2016(include_lifetimes=True)['records']
    agb = [r for r in agb if r['coordinate']['initial_mass_msun'] <= 7]
    agb += read_fishlock2014(lifetime_model='raiteri96_padova', include_one=True)
    # Existing explicitly selected stable-segment residual-remnant comparison.
    table, _ = build_lowmass(solar_coordinate=.02, wind_speed_km_s=1000.,
        lifetime_model='raiteri96_padova', remnant_policy='stable_segment_residual',
        max_remnant_difference_msun=.1, imf_id=imf, source_root=sukhbold_root)
    selection = json.loads(table.splitlines()[1][2:])
    return agb, selection


def endpoint(m, z, agb, sn):
    if m < 9:
        zs = sorted({r['coordinate']['metallicity_mass_fraction'] for r in agb})
        sz = min(zs, key=lambda x: (abs(x-z), x))
        choices = [r for r in agb if r['coordinate']['metallicity_mass_fraction'] == sz]
        r = nearest(choices, min(m, 7.), lambda x: x['coordinate']['initial_mass_msun'])
        sm = r['coordinate']['initial_mass_msun']; f = m/sm
        g = np.array(r['selected_ejecta']['tracked_ejected_mass_msun'])*f
        ret = r['selected_ejecta']['returned_mass_msun']*f
        rem = r['selected_ejecta']['remnant_mass_msun']*f
        kind = {'CO': 0, 'CO_Ne': 1, 'ONe': 2}[r['evolution']['core_kind']]
        if m > 7: kind = 3
        meta = dict(source='Fishlock2014' if sz == .001 else 'KL16', source_mass=sm, source_z=sz,
                    source_coordinate=r['coordinate'], overshoot_label=r.get('overshoot_label'),
                    same_mass_tie_policy='first_in_pinned_reader_order',
                    lifetime_source=r['evolution']['lifetime_source'], envelope_proxy=m > 7,
                    mass_policy='nearest_source_fractional_budget', z_policy='nearest_source_Z_proxy')
        return g, ret, rem, r['evolution']['stellar_lifetime_yr'], 2, kind, 0, 0., meta
    r = nearest(sn['models'], m, lambda x: x['mass']); f = m/r['mass']
    g = (np.array(r['wind_elements'])+r['terminal_elements'])*f
    ret = (r['wind']+r['terminal'])*f
    return g, ret, r['remnant']*f, r['age'], 3, 4, 1, r['terminal_energy']*f, dict(
        source='Sukhbold2016_stable_residual', source_mass=r['mass'], source_z=.02,
        lifetime_source='Raiteri96_Padova_solar_fit', mass_policy='nearest_source_fractional_budget',
        z_policy='solar_terminal_proxy_at_all_Z', envelope_proxy=False)


def clip_series(t, values, end):
    if end <= 0 or t[0] > end:
        raise ValueError('no source interval at selected end')
    k = np.searchsorted(t, end, side='left')
    if k < len(t) and t[k] == end:
        return t[:k+1], values[:k+1]
    if k == len(t):
        raise ValueError('attempted extrapolation beyond source interval')
    w = (end-t[k-1])/(t[k]-t[k-1])
    return np.r_[t[:k], end], np.vstack((values[:k], (1-w)*values[k-1]+w*values[k]))


def wind_split(track, initial, m, sm, terminal, gross, returned, speed):
    # Literal nonnegative RATE, no inferred chemical reservoir from missing light.
    t = track[:, 1]
    if np.any(np.diff(t) < 0) or np.any(track[:, 28] < 0):
        raise ValueError('negative RATE / reversed track age')
    keep = np.r_[np.diff(t) > 0, True]; t = t[keep]; track = track[keep]
    end = min(terminal, t[-1])
    # First source state held only across the small initial sampling interval.
    t = np.r_[0., t]; track = np.vstack((track[0], track))
    t, track = clip_series(t, track, end)
    surface = np.tile(initial, (len(t), 1))
    surface[:, :2] = track[:, 29:31]
    surface[:, 2] = track[:, 31]+track[:, 32]
    surface[:, 3] = track[:, 33]
    surface[:, 4] = track[:, 34]+track[:, 35]
    surface[:, 5:7] = track[:, 36:38]
    if np.any(surface < 0) or not np.isfinite(surface).all():
        raise ValueError('invalid surface composition')
    # Explicit comparison normalization only when printed projection exceeds1.
    # No elemental zero fill: heavy elements use the same PARSEC initial mixture.
    norm = np.maximum(1., surface.sum(axis=1))
    surface /= norm[:, None]
    untracked = 1-surface.sum(axis=1)
    if np.min(untracked) < -8*np.finfo(float).eps:
        raise ValueError('surface normalization does not close')
    # Bounded final arithmetic representative inside the normalization model.
    surface12 = np.column_stack((surface, np.maximum(untracked, 0.)))
    surface12 /= surface12.sum(axis=1)[:, None]
    rate = track[:, 28]*(m/sm)
    flux = rate[:, None]*surface12
    release = np.diff(t)[:, None]*(flux[1:]+flux[:-1])*.5
    wind = np.vstack((np.zeros(12), np.cumsum(release, axis=0)))
    budget = np.r_[gross, returned-gross.sum()]
    if np.any(budget < 0):
        raise ValueError('negative terminal-source untracked budget')
    positive = wind[-1] > 0
    lam = min(1., float(np.min(budget[positive]/wind[-1, positive]))) if positive.any() else 1.
    if lam < 1: lam = np.nextafter(lam, 0.)
    wind *= lam
    residual = budget-wind[-1]
    if np.any(residual < 0):
        raise ValueError('wind allocation exceeds terminal-source reservoir')
    if terminal > end:
        t = np.r_[t, terminal]; wind = np.vstack((wind, wind[-1]))
    index, error = compress(t, wind, 1e-4)
    return t[index], wind[index], residual, dict(wind_budget_scale=lam,
        max_surface_normalization=float(norm.max()), wind_compression_error=error,
        wind_stop_yr=end, wind_speed_km_s=speed)


def lower_radiation(data, track, terminal, edges, scale, spectrum):
    if np.any(np.diff(data[:, 0]) < 0) or not np.isfinite(data).all():
        raise ValueError('invalid photon source')
    data = data[np.r_[np.diff(data[:, 0]) > 0, True]]
    rates, ratio = spectrum.constrained(data, edges)
    check_rates(rates, edges)
    stop = min(track[-1, 1], terminal)/1e6
    t = np.r_[0., data[:, 0]]; rates = np.vstack((rates[0], rates))
    measured_end = t[-1]
    if stop <= measured_end:
        t, rates = clip_series(t, rates, stop)
        rel = intervals(t, rates)
    else:
        rel = intervals(t, rates)
        # Actual evolving full-track L/Teff only; no fabricated late plateau.
        tt = track[:, 1]/1e6
        tail = np.unique(np.r_[measured_end, tt[(tt > measured_end) & (tt < stop)], stop])
        tail_rate = spectrum.planck(np.interp(tail, tt, track[:, 3]),
                                    np.interp(tail, tt, track[:, 4]), edges)
        check_rates(tail_rate, edges)
        rel = np.vstack((rel, intervals(tail, tail_rate)[1:]))
        t = np.r_[t, tail[1:]]
    rel *= scale
    cumulative = np.cumsum(rel, axis=0)
    index, error = compress(t, cumulative, 1e-4)
    compact = np.vstack((np.zeros(2*(len(edges)-1)+3),
                         [rel[a+1:b+1].sum(axis=0) for a,b in zip(index[:-1], index[1:])]))
    check_rates(compact[1:], edges)
    # Round-trip Myr->yr must never put the cutoff one ulp AFTER death.
    return t[index], compact, dict(radiation_stop_yr=min(terminal,stop*1e6),
        omitted_after_track=True, q_compression_error=error,
        measured_Q_end_yr=measured_end*1e6, max_Q_energy_over_Lbol=float(ratio.max()))


def row(channel, m, z, age, ret=0., rem=0., energy=0., gross=None):
    g = np.zeros(11) if gross is None else np.asarray(gross)
    return np.r_[channel,m,z,age,ret,rem,energy,np.zeros(3),g,np.zeros(11)]


def load_upper_sed(path, *, group_set_id, group_edges_sha256, group_edges_ev):
    edges = np.asarray(group_edges_ev, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all() or \
            np.any(edges[1:] <= edges[:-1]):
        raise ValueError('upper radiation group edges are invalid')
    ng = edges.size - 1
    nodes = {}
    with path.open() as f:
        magic = f.readline().strip()
        if magic not in ('SNRT_PARSEC_QE_INTERVAL_V1', 'SNRT_POPULATION_QE_INTERVAL_V3'):
            raise ValueError('upper radiation must be a supported v1 or profile-bound v3 payload')
        try:
            nn, nr = map(int, f.readline().split())
        except (ValueError, TypeError) as exc:
            raise ValueError('upper radiation node/record header is malformed') from exc
        model = f.readline().strip()
        if not model:
            raise ValueError('upper radiation model identity is empty')
        if magic == 'SNRT_POPULATION_QE_INTERVAL_V3':
            identity = f.readline().split()
            if len(identity) != 2:
                raise ValueError('profile-bound upper radiation identity is malformed')
            payload_profile, payload_edges_sha = identity
            if payload_profile != group_set_id or payload_edges_sha != group_edges_sha256:
                raise ValueError('upper radiation profile ID or edge hash differs from its manifest')
        elif group_set_id != 'snrt_agn9_v1':
            raise ValueError('legacy v1 upper radiation is only valid for snrt_agn9_v1')
        if nn < 2 or nn > 1024 or nr < 2 * nn or nr > 500000:
            raise ValueError('upper radiation node/record count is outside the native contract')
        total = 0
        ncolumns = 2 * ng + 4  # age, group N/E, Lbol, IR tail, X-ray tail
        for _ in range(nn):
            header = f.readline().split()
            expected_header_fields = 5 if magic == 'SNRT_PARSEC_QE_INTERVAL_V1' else 8
            if len(header) != expected_header_fields:
                raise ValueError('upper radiation node header is malformed')
            try:
                if magic == 'SNRT_PARSEC_QE_INTERVAL_V1':
                    m, z, end, fate_value, n_value = map(float, header)
                    radiation_stop, terminal_channel, remnant_kind = end, 3, 4
                else:
                    (m, z, end, radiation_stop, terminal_channel,
                     remnant_kind, fate_value, n_value) = map(float, header)
            except ValueError as exc:
                raise ValueError('upper radiation node header is not numeric') from exc
            integer_fields = [terminal_channel, remnant_kind, fate_value, n_value]
            if (not np.isfinite([m, z, end, radiation_stop, *integer_fields]).all() or m <= 0 or
                    z < 0 or z > 1 or end <= 0 or fate_value != int(fate_value) or
                    radiation_stop <= 0 or radiation_stop > end or
                    any(value != int(value) for value in integer_fields)):
                raise ValueError('upper radiation node metadata is invalid')
            terminal_channel, remnant_kind = int(terminal_channel), int(remnant_kind)
            fate, n = int(fate_value), int(n_value)
            if magic == 'SNRT_PARSEC_QE_INTERVAL_V1':
                if not 1 <= fate <= 5:
                    raise ValueError('legacy upper radiation fate is invalid')
            elif terminal_channel == 2:
                if m >= 9 or fate != 0 or not 0 <= remnant_kind <= 3:
                    raise ValueError('upper radiation low-mass terminal metadata is inconsistent')
            elif (terminal_channel not in (3, 5) or m < 9 or not 1 <= fate <= 5 or
                  remnant_kind != 4 or ((terminal_channel == 5) != (fate in (3, 4)))):
                raise ValueError('upper radiation high-mass terminal metadata is inconsistent')
            if n < 2 or n > nr - total:
                raise ValueError('upper radiation node record count is invalid')
            rows = []
            for _ in range(n):
                try:
                    row = list(map(float, f.readline().split()))
                except ValueError as exc:
                    raise ValueError('upper radiation interval row is not numeric') from exc
                if len(row) != ncolumns:
                    raise ValueError('upper radiation interval row has the wrong group layout')
                rows.append(row)
            values = np.asarray(rows, dtype=float)
            if not np.isfinite(values).all() or np.any(values[:, 1:] < 0):
                raise ValueError('upper radiation intervals contain nonfinite or negative values')
            if values[0, 0] != 0 or np.any(values[0, 1:] != 0) or \
                    np.any(values[1:, 0] <= values[:-1, 0]) or values[-1, 0] != radiation_stop:
                raise ValueError('upper radiation age/zero-interval contract is invalid')
            photons = values[:, 1:1 + ng]
            energy = values[:, 1 + ng:1 + 2 * ng]
            lbol = values[:, 1 + 2 * ng]
            tails = values[:, 2 + 2 * ng:4 + 2 * ng]
            if (np.any(energy < photons * edges[:-1] * (1 - 1e-10)) or
                    np.any(energy > photons * edges[1:] * (1 + 1e-10))):
                raise ValueError('upper radiation energy lies outside its photon group')
            residual = np.abs(energy.sum(axis=1) + tails.sum(axis=1) - lbol)
            if np.any(residual > 1e-10 * np.maximum(np.finfo(float).tiny, lbol)):
                raise ValueError('upper radiation bolometric energy does not close')
            key = (z, m)
            if key in nodes:
                raise ValueError('upper radiation contains duplicate metallicity/mass nodes')
            nodes[key] = (end, fate, values)
            total += n
        if total != nr or f.read().strip():
            raise ValueError('upper radiation extent mismatch')
    return model,nodes


def build(args):
    output = args.output.resolve()
    if output.exists(): raise ValueError('output directory already exists')
    feedback = json.loads((args.feedback/'manifest.json').read_text())
    sed = json.loads((args.upper_sed/'manifest.json').read_text())
    for name, digest in feedback['output_sha256'].items():
        if sha(args.feedback/name) != digest: raise ValueError('changed upper material package')
    if sha(args.upper_sed/'nodes.dat') != sed['source_sha256']:
        raise ValueError('changed upper radiation package')
    if sed['feedback_outputs_sha256'] != feedback['output_sha256']:
        raise ValueError('upper radiation/material source mismatch')
    inputs = dict(feedback['inputs_sha256']); inputs.update(sed['sources_sha256'])
    for name, digest in inputs.items():
        if sha(args.source_dir/name) != digest: raise ValueError(f'changed source {name}')
    coordinates = Z_GRIDS[feedback.get('metallicity_grid','solar_pair')]
    if feedback['source_domain']['z'] != [float(z) for z,_ in coordinates]:
        raise ValueError('unexpected upper Z grid')
    if not all(np.isfinite(v) and 0 < v < 29979.2458 for v in (args.agb_wind_km_s,args.sn_wind_km_s)):
        raise ValueError('explicit positive nonrelativistic lower wind speeds required')
    imf = feedback['imf_id']; agb, sn = terminal_sources(imf,args.sukhbold_root)
    edges = np.array(sed['group_edges_ev']); ng=len(edges)-1
    group_set_id=sed.get('group_set_id','snrt_agn9_v1')
    profile_ng,profile_edges,registered_sha=read_profile(
        ROOT/'config/snrt_group_profiles_v1.json',group_set_id)
    group_edges_sha256=sed.get('group_edges_sha256',registered_sha)
    if (len(edges)!=profile_ng+1 or not np.array_equal(edges,profile_edges) or
            group_edges_sha256!=registered_sha):
        raise ValueError('upper SED profile ID, edge count, values or SHA-256 differs from registry')
    upper_model, radiation = load_upper_sed(
        args.upper_sed/'nodes.dat', group_set_id=group_set_id,
        group_edges_sha256=group_edges_sha256, group_edges_ev=edges)
    if upper_model != feedback['model']: raise ValueError('upper model identity mismatch')
    upper = np.loadtxt(args.feedback/'yields.dat')
    spectrum = Spectrum(64)
    material=[]; nodes=[]; spectra={}
    with zipfile.ZipFile(args.source_dir/'all_ejecta.zip') as ejecta:
        for zs,ys in coordinates:
            z = float(zs)
            initial,_ = read_ejecta(ejecta,f'ejecta/Z{zs}_Y{ys}_winds_ejecta.dat')
            with zipfile.ZipFile(args.source_dir/source_name(zs,ys,'tracks')) as tracks, \
                 zipfile.ZipFile(args.source_dir/source_name(zs,ys,'photons')) as photons:
                for m in LOW_MASSES:
                    sm = min(m,12.)
                    track = np.loadtxt(io.BytesIO(tracks.read(source_name(zs,ys,'tracks',sm))),skiprows=3)
                    data = np.loadtxt(io.BytesIO(photons.read(source_name(zs,ys,'photons',sm))),skiprows=3)
                    g,ret,rem,terminal,ch,kind,fate,energy,meta = endpoint(m,z,agb,sn)
                    speed = args.agb_wind_km_s if ch == 2 else args.sn_wind_km_s
                    tw,wind,residual,wm = wind_split(track,initial,m,sm,terminal,g,ret,speed)
                    tq,qe,qm = lower_radiation(data,track,terminal,edges,m/sm*sed['escaped_fraction_applied_once'],spectrum)
                    for t,v in zip(tw,wind):
                        material.append(row(1,m,z,t,v.sum(),0.,.5*v.sum()*MSUN*(speed*1e5)**2,v[:11]))
                    for c in (2,3,5):
                        material.append(row(c,m,z,0.))
                        if c == ch:
                            e = energy if ch == 3 else .5*residual.sum()*MSUN*(speed*1e5)**2
                            material.append(row(c,m,z,terminal,residual.sum(),rem,e,residual[:11]))
                        else:
                            material.append(row(c,m,z,terminal))
                    nodes.append(dict(mass=m,z=z,age_yr=terminal,terminal_channel=ch,remnant_kind=kind,fate=fate,
                        radiation_source_mass=sm,radiation_proxy=m==13.,total_return=ret,remnant=rem,
                        terminal_source=meta,**wm,**qm))
                    spectra[z,m] = np.column_stack((tq,qe))
            for original in feedback['nodes']:
                if original['z'] != z: continue
                m=original['mass']; end=original['age_yr']; fate=original['fate']
                upper_rows=upper[(upper[:,1]==m)&(upper[:,2]==z)].copy()
                # Mixed lifetime sources have no shared initial-abundance
                # convention for NET yields. Suppress diagnostics uniformly,
                # never gross material, under the declared unavailable policy.
                upper_rows[:,21:]=0.
                material.extend(upper_rows)
                material.extend((row(2,m,z,0.),row(2,m,z,end)))
                rad_end,rad_fate,values=radiation[z,m]
                if rad_fate!=fate or abs(rad_end-end/1e6)>1e-12*rad_end:
                    raise ValueError('upper radiation node mismatch')
                spectra[z,m]=values
                nodes.append(dict(mass=m,z=z,age_yr=end,radiation_stop_yr=end,
                    terminal_channel=5 if fate in (3,4) else 3,remnant_kind=4,fate=fate,
                    terminal_source=dict(source=feedback['model']),radiation_proxy=False))
            print(f'Converted Z={z:g}',flush=True)
    nodes.sort(key=lambda n:(n['z'],n['mass']))
    if len(nodes)>1024: raise ValueError('v5 node capacity exceeded')
    material=np.asarray(material)
    material=material[np.lexsort((material[:,3],material[:,0],material[:,1],material[:,2]))]
    if not np.isfinite(material).all() or np.any(material[:,4:7]<0) or np.any(material[:,10:21]<0):
        raise ValueError('invalid canonical payload')
    if np.any(material[:,10:21].sum(axis=1)>material[:,4]+1e-10*material[:,1]):
        raise ValueError('canonical gross mass exceeds return')
    lines=['&stellar_high_mass_history',' version=5',f' node_count={len(nodes)}',
        f" model_id='{MODEL}'", " wind_source_id='PARSEC_RATE_budgeted_mixed_v1'",
        " terminal_source_id='KL16_Fishlock_Sukhbold_PARSEC_mixed_v1'",
        " model_coordinates='nearest_terminal_M_Z_proxies;Q13=13over12_Q12;no_late_SED'",
        " timing_policy='wind_linear_terminal_step'", " metallicity_policy='linear_Z_cumulative_mixture'",
        " net_yield_policy='unavailable_diagnostic_zero'", " agb_release_policy='terminal_step'",
        f' input_imf_id={imf}', ' input_population_id=0, input_binary_fraction=0',
        ' input_imf_min=.08, input_imf_max=600, mass_domain_min=2, mass_domain_max=600']
    for i,n in enumerate(nodes,1):
        for key,value in (('mass_msun',n['mass']),('metallicity',n['z']),('terminal_age_yr',n['age_yr']),
                          ('radiation_stop_age_yr',n['radiation_stop_yr'])):
            lines.append(f' {key}({i})={value:.17e}')
        for key in ('terminal_channel','remnant_kind','fate'):
            lines.append(f" {'terminal_fate' if key=='fate' else key}({i})={n[key]}")
    history='\n'.join(lines+['/'])+'\n'
    nr=sum(len(v) for v in spectra.values())
    if nr>500000: raise ValueError(f'v5 radiation knot capacity exceeded: {nr}')
    candidate_profile=group_set_id!='snrt_agn9_v1'
    record_magic='SNRT_POPULATION_QE_INTERVAL_V3' if candidate_profile else 'SNRT_POPULATION_QE_INTERVAL_V2'
    lines=[record_magic,f'{len(nodes)} {nr}',MODEL]
    if candidate_profile:
        lines.append(f'{group_set_id} {group_edges_sha256}')
    for n in nodes:
        v=spectra[n['z'],n['mass']]
        lines.append(' '.join(f'{x:.17e}' for x in (n['mass'],n['z'],n['age_yr']/1e6,v[-1,0]))+
            f" {n['terminal_channel']} {n['remnant_kind']} {n['fate']} {len(v)}")
        lines.extend(' '.join(f'{x:.17e}' for x in r) for r in v)
    payload='\n'.join(lines)+'\n'; digest=hashlib.sha256(payload.encode()).hexdigest()
    template=(args.upper_sed/'source.nml').read_text()
    template=re.sub(r'^!.*\n','! Explicit mixed comparison; omitted late radiation; not full SSP.\n',template,count=1)
    template=template.replace('version=4','version=6' if candidate_profile else 'version=5')
    template=template.replace('match_feedback_high_mass_only','common_imf_mixed_evolution')
    template=template.replace('PARSEC_v2_nonrot_high_mass_only','PARSEC_mixed_lowmass_comparison')
    template=template.replace('Q5_Planck_and_track_tail_v1','zero_after_track_or_terminal_v1')
    template=re.sub(r"source_sha256='[^']*'",f"source_sha256='{digest}'",template)
    template=re.sub(r"node_history_file='[^']*'",f"node_history_file='{output/'nodes.dat'}'",template)
    config=f"""&stellar_enrichment_params
 feedback_mode='channel_resolved', fate_policy='user_selected_model_v1',
 high_mass_preset='source_consistent', high_mass_history_path='{output/'history.nml'}',
 yield_source_basis='per_star_cumulative', population_model='single_star_ssp',
 imf_id={imf}, imf_mass_min_msun=.08, imf_mass_max_msun=600, binary_fraction=0,
 channel_mass_min_msun=2,2,2,3,2, channel_mass_max_msun=600,600,600,8,600,
 use_wind=.true.,use_agb=.true.,use_snii=.true.,use_snia=.false.,use_pisn=.true.
/
"""
    output.mkdir(parents=True,exist_ok=False)
    np.savetxt(output/'yields.dat',material,fmt=['%d']+['%.17e']*31,
               header='Explicit mixed lifetime budgets; native v5 comparison only')
    (output/'history.nml').write_text(history);(output/'nodes.dat').write_text(payload)
    (output/'source.nml').write_text(template);(output/'enrichment.nml').write_text(config)
    manifest=dict(model=MODEL,full_ssp=False,same_evolution=False,group_set_id=group_set_id,
        group_edges_sha256=group_edges_sha256,group_edges_ev=edges.tolist(),
        source_domain=dict(mass=[2,600],z=[float(z) for z,_ in coordinates]),
        imf_id=imf,imf_denominator=[.08,600],input_upper=feedback['model'],upper_output_sha256=feedback['output_sha256'],
        parsec_inputs_sha256=inputs,upper_sed_manifest_sha256=sha(args.upper_sed/'manifest.json'),
        agb_sources=dict(selection_matrix_sha256=sha(DEFAULT_MATRIX),
            kl16=next(c for c in json.loads(DEFAULT_MATRIX.read_text())['candidates']
                      if c['candidate_id']=='karakas_lugaro2016_agb'),
            lifetime_path=str(LIFETIME_PATH),lifetime_sha256=LIFETIME_SHA256,
            fishlock_sha256=FISHLOCK_SHA256),sukhbold_selection=sn,nodes=nodes,
        radiation_knots=nr,material_rows=len(material),generator_sha256=sha(Path(__file__)),
        numerical_policy='positive_interval_QE;budgeted_wind;terminal_residual;untracked_reservoir_included',
        wind_mass_policy='PARSEC_RATE_integral_scaled_only_down_to_available_terminal_budget',
        late_radiation='zero;not_missing_photons_converted_to_heat',
        net_yields='unavailable_diagnostic_zero_for_entire_mixed_population',
        output_sha256={p.name:sha(p) for p in output.iterdir()})
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(nodes=len(nodes),radiation_knots=nr,material_rows=len(material),output=str(output))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--feedback',type=Path,required=True)
    p.add_argument('--upper-sed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sukhbold-root',type=Path,default=DEFAULT_SUKHBOLD_ROOT,
        help='parent of the pinned external/g2_candidates/sukhbold2016_ccsn archive')
    p.add_argument('--agb-wind-km-s',type=float,required=True)
    p.add_argument('--sn-wind-km-s',type=float,required=True)
    build(p.parse_args())

#!/usr/bin/env python3
"""Offline Verner/CHIMES photoionization data for SNRT photon groups.

This is a declared within-bin closure, not a resolved stellar/AGN spectrum:
J_nu is a power law whose photon-number-weighted energy equals the group's
configured mean. The native receiver uses the resulting HDF5 tables only.
The CHIMES tools checkout must carry the accompanying compatibility patch.
With --band-nodes, prepare separate partial-shell cross-section nodes for
the native N/E maxent building block instead; no grey tables are replaced.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

EDGES = [.01, 1., 5.6, 11.2, 13.6, 24.59, 54.42, 500., 2000., 10000.]
MEANS = [.1, 2.3664319132398464, 10.56942670578782, 12.2954112566775,
         17.662918001204492, 34.381528542610674, 106.63925876378897,
         869.6341490248457, 4023.594574013186]
MAIN_SHA = '8bde78faacd59249fad5e810cae43311ed03ef09131c62b0a0a5c39e8407cb0a'
PHOTO_GROUPS = ('photoion_fuv', 'photoion_euv', 'photoion_auger_fuv', 'photoion_auger_euv')
FIT_SHA = {
    'cross_section_fits_verner95.dat': '22d3afe50371471ce60838d788fb72433ba2352dc23894798c0157be2d14fd09',
    'cross_section_fits_verner96.dat': 'e8b7883ad2bd683d146c4f96b9598101a32933addffadc4b1927cdbde9e10f2e',
    'auger_ionisation_probabilities_kaastra93.dat': '48afae9a3bbeca9b9f10db036ea38d85981b6b9facc7e6c8133a119673db1b79',
}


def build_shell_bank(tools, main_data, output, n):
    """Native atomic/H-/H2-ionization bank, NOT a molecular dissociation bank.

    Keep each partial shell separate: FS2010 evaluated at a shell-averaged
    electron energy is not the average of its nonlinear energy partition.
    The existing grey tables and native CHIMES admission are untouched.
    """
    raw = tools/'generate_cross_sections/data'
    for name, digest in FIT_SHA.items():
        if sha(raw/name) != digest:
            raise ValueError(f'unmatched atomic fit input: {name}')
    sys.path.insert(0, str(tools))
    generator = tools/'generate_cross_sections/generate_cross_sections.py'
    spec = importlib.util.spec_from_file_location('chimes_band_fits', generator)
    fits = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fits)
    v96 = fits.read_in_V96_cross_sections(raw/'cross_section_fits_verner96.dat')
    v95, shells = fits.read_in_V95_cross_sections(raw/'cross_section_fits_verner95.dat')
    auger = fits.read_in_KM93_tables(raw/'auger_ionisation_probabilities_kaastra93.dat')
    edges = np.asarray(EDGES)
    ev = np.array([lo*np.exp(np.log(hi/lo)*np.arange(n)/(n-1))
                   for lo, hi in zip(edges[:-1], edges[1:])])
    ev[:, 0] = edges[:-1]; ev[:, -1] = edges[1:]
    evaluation = ev.copy()
    evaluation[:-1, -1] = np.nextafter(edges[1:-1], edges[:-2])
    neutral = [1, 4, 7, 15, 23, 33, 44, 57, 72, 89, 110]
    charges = [1, 2, 6, 7, 8, 10, 12, 14, 16, 20, 26]
    reactions, shell_reaction, binding, sigma = [], [], [], []

    def add(reaction, threshold, limit, probability, cross):
        if probability == 0:
            return
        if not np.isfinite([threshold, limit, probability]).all() or threshold <= 0 or not 0 < probability <= 1:
            raise ValueError('invalid partial shell')
        sample = np.zeros(ev.shape)
        valid = (evaluation >= threshold) & (evaluation <= limit)
        sample[valid] = probability*np.array([cross(e) for e in evaluation[valid]])
        if not np.isfinite(sample).all() or (sample < 0).any():
            raise ValueError('invalid cross section')
        shell_reaction.append(reaction); binding.append(threshold); sigma.append(sample)

    with h5py.File(main_data) as source:
        for group_id, group in enumerate(PHOTO_GROUPS):
            for index, (species, products) in enumerate(zip(source[group+'/reactants'][:], source[group+'/products'][:])):
                species = int(species)
                released = int(products[0]-species) if 'auger' in group else 1
                reaction = len(reactions)
                reactions.append([group_id, index, species, int(products[0]), released])
                if species in (3, 137):
                    threshold, cross = (.755, fits.sigma_Hminus) if species == 3 else (15.4, fits.sigma_H2)
                    add(reaction, threshold, v96[0][0][1], 1., cross)
                    continue
                element = next((i for i, start in enumerate(neutral) if start <= species < start+charges[i]), None)
                if element is None:
                    raise ValueError(f'unmapped photoionization reactant {species}')
                ion = species-neutral[element]
                if element < 2:
                    p = v96[element][ion]
                    add(reaction, p[0], p[1], 1., lambda e: fits.sigma(e, p[2:9], 0))
                    continue
                metal = element-2
                for shell in range(shells[metal][ion]):
                    probability = auger[metal][ion][shell][released-1]
                    if released == 1 and shell == shells[metal][ion]-1:
                        p = v96[element][ion]
                        add(reaction, p[0], p[1], probability, lambda e: fits.sigma(e, p[2:9], 0))
                    else:
                        k = shells[metal][ion]-shell-1
                        p = v95[metal][ion, 1:6, k]
                        angular = v95[-1][metal][ion][k][1]
                        add(reaction, v95[metal][ion, 0, k], edges[-1], probability,
                            lambda e: fits.sigma_V95(e, p, angular))
    output.mkdir(parents=True, exist_ok=False)
    path = output/'atomic_shells.h5'
    with h5py.File(path, 'x') as f:
        f.attrs['schema'] = 'snrt_chimes_atomic_shell_nodes_v1'
        f.attrs['main_sha256'] = MAIN_SHA
        f.attrs['molecular_dissociation_supported'] = 0
        f.attrs['units'] = 'energies eV; probability-weighted sigma cm2; zero-based reaction/species indices'
        f['edges_ev'] = edges; f['node_ev'] = ev
        f['reactions'] = np.asarray(reactions, dtype=np.int32)
        f['shell_reaction'] = np.asarray(shell_reaction, dtype=np.int32)
        f['binding_ev'] = np.asarray(binding)
        f.create_dataset('sigma_cm2', data=np.asarray(sigma), compression='gzip', shuffle=True)
    manifest = dict(schema='snrt_chimes_atomic_shell_nodes_v1', main_sha256=MAIN_SHA,
                    fit_sha256=FIT_SHA, tools_generator_sha256=sha(generator), wrapper_sha256=sha(Path(__file__)),
                    table_sha256=sha(path), nodes=n, reactions=len(reactions), partial_shells=len(binding),
                    edges_ev=EDGES, group_order=PHOTO_GROUPS,
                    limitations=['Atomic/H-/H2 ionization only; NOT the 30 empirical molecular dissociation channels or Leiden lines.',
                                 'Same Verner95/96 and KM93 branching as existing grey model; no new cascade/fluorescence energy.',
                                 'Finite node quadrature; upper endpoint uses left-limit opacity except final group.',
                                 'Not yet a live CVODE N/E receiver. Spectral+CHIMES remains rejected.'])
    with (output/'manifest.json').open('x') as f:
        json.dump(manifest, f, indent=2); f.write('\n')
    print(json.dumps(manifest, indent=2))


def primary_electron_moments(tools, main_data, table, lo, hi, slope):
    """Integrate the same partial-shell cross sections used for absorption.

    The original broad-spectrum heating integrals use only the V96 outer
    shell. In a hard-only RT bin that integral can be zero while V95 inner
    shell absorption remains nonzero; dividing two floored zeros gives 1 erg.
    Here primary electron energy is E - shell binding, not E - the lowest
    ionization potential. Cascade/Auger relaxation energy is not guessed.
    """
    generator = tools/'generate_cross_sections/generate_cross_sections.py'
    sys.path.insert(0, str(tools))
    spec = importlib.util.spec_from_file_location('chimes_cross_sections', generator)
    fits = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fits)
    raw = tools/'generate_cross_sections/data'
    v96 = fits.read_in_V96_cross_sections(raw/'cross_section_fits_verner96.dat')
    v95, shells = fits.read_in_V95_cross_sections(raw/'cross_section_fits_verner95.dat')
    auger = fits.read_in_KM93_tables(raw/'auger_ionisation_probabilities_kaastra93.dat')
    neutral = [1, 4, 7, 15, 23, 33, 44, 57, 72, 89, 110]
    charges = [1, 2, 6, 7, 8, 10, 12, 14, 16, 20, 26]
    pivot = np.log(hi if slope > 0 else lo)
    norm = quad(lambda x: np.exp(slope*(x-pivot)), np.log(lo), np.log(hi), epsabs=1e-12)[0]

    def integrals(threshold, limit, cross):
        lower, upper = max(lo, threshold), min(hi, limit)
        if lower >= upper:
            return 0., 0.
        # Scale sigma for the quadrature's absolute tolerance, undo below.
        weight = lambda x: np.exp(slope*(x-pivot))*cross(np.exp(x))*1e18
        number = quad(weight, np.log(lower), np.log(upper), epsabs=1e-12, epsrel=1e-7, limit=300)[0]
        heat = quad(lambda x: weight(x)*(np.exp(x)-threshold), np.log(lower), np.log(upper),
                    epsabs=1e-12, epsrel=1e-7, limit=300)[0]
        return number/norm*1e-18, heat/norm*1e-18*1.602176634e-12

    def moment(species, released):
        if species == 3:
            return integrals(.755, v96[0][0][1], fits.sigma_Hminus)
        if species == 137:
            return integrals(15.4, v96[0][0][1], fits.sigma_H2)
        for element, start in enumerate(neutral):
            if start <= species < start+charges[element]:
                ion = species-start
                break
        else:
            raise ValueError(f'unrecognized photoionization reactant {species}')
        if element < 2:
            pars = v96[element][ion]
            return integrals(pars[0], pars[1], lambda e: fits.sigma(e, pars[2:9], 0))
        metal = element-2
        n_shells = shells[metal][ion]
        number = heat = 0.
        for shell in range(n_shells):
            probability = auger[metal][ion][shell][released-1]
            if probability == 0:
                continue
            if released == 1 and shell == n_shells-1:
                pars = v96[element][ion]
                den, num = integrals(pars[0], pars[1], lambda e: fits.sigma(e, pars[2:9], 0))
            else:
                index = n_shells-shell-1
                threshold = v95[metal][ion, 0, index]
                pars = v95[metal][ion, 1:6, index]
                angular = v95[-1][metal][ion][index][1]
                den, num = integrals(threshold, hi, lambda e: fits.sigma_V95(e, pars, angular))
            number += probability*den
            heat += probability*num
        return number, heat

    with h5py.File(main_data) as source, h5py.File(table, 'r+') as dest:
        dest['Header/snrt_photoelectron_closure'] = np.bytes_('primary_shell_electron_v1')
        for group in ['photoion_fuv', 'photoion_euv', 'photoion_auger_fuv', 'photoion_auger_euv']:
            reactants = source[group+'/reactants'][:]
            products = source[group+'/products'][:]
            energy = np.zeros(len(reactants))
            sigma = dest[group+'/sigmaPhot'][:].astype(float)
            for i, species in enumerate(reactants):
                released = int(products[i, 0]-species) if 'auger' in group else 1
                den, num = moment(int(species), released)
                if abs(den-sigma[i]) > 3e-5*max(den, sigma[i], 1e-35):
                    raise ValueError(f'absorption/heating cross sections disagree: {group}:{i}: {den} vs {sigma[i]}')
                energy[i] = num/den if den > 0 else 0.
                if not 0 <= energy[i] <= hi*1.602176634e-12:
                    raise ValueError('photoelectron energy exceeds photon energy')
            dest['snrt_primary_photoelectron_energy_erg/'+group] = energy
            if group == 'photoion_fuv':
                dest[group+'/epsilonPhot'][...] = energy
            elif group == 'photoion_euv':
                # CHIMES RT reads only the zero-column numerator/denominator.
                # It does not use these arrays for spatial shielding in RT.
                one = np.full(dest[group+'/shieldFactor_1D'].shape, -100., dtype=np.float32)
                two = np.full(dest[group+'/shieldFactor_2D'].shape, -100., dtype=np.float32)
                two[:, 3, 0, 0] = np.log10(np.maximum(energy, 1e-100))
                two[:, 5, 0, 0] = 0.
                dest[group+'/shieldFactor_1D'][...] = one
                dest[group+'/shieldFactor_2D'][...] = two


def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def molecular_moments(directory, table, lo, hi, slope):
    """Leiden line-resolved optically thin moments; shielding stays native.

    Non-ionizing absorption includes fluorescent events. It is distinct from
    the dissociation cross section and excludes already-counted ionization.
    """
    pivot = np.log(hi if slope > 0 else lo)
    norm = quad(lambda x: np.exp(slope*(x-pivot)), np.log(lo), np.log(hi))[0]
    with h5py.File(table, 'r+') as dest:
        for molecule in ['H2', 'CO']:
            path = directory/f'{molecule}.h5'
            with h5py.File(path) as data:
                wavelength = data['wavelength'][:]
                diss = data['photodissociation'][:]
                absorption = np.maximum(diss, data['photoabsorption'][:]-data['photoionisation'][:])
            if not (np.diff(wavelength)>0).all():
                raise ValueError('unordered Leiden wavelengths')
            lower, upper = max(wavelength[0],1239.8419843320026/hi), min(wavelength[-1],1239.8419843320026/lo)
            values = []
            for cross in [absorption, diss]:
                if lower >= upper:
                    values.append(0.)
                    continue
                grid = np.concatenate(([lower],wavelength[(wavelength>lower)&(wavelength<upper)],[upper]))
                weight = np.exp(slope*(np.log(1239.8419843320026/grid)-pivot))/grid
                values.append(float(np.trapezoid(weight*np.interp(grid,wavelength,cross),grid)/norm))
            if not np.isfinite(values).all() or not 0 <= values[1] <= values[0]:
                raise ValueError('invalid molecular cross sections')
            dest['snrt_molecule_sigma/'+molecule] = values
            dest['snrt_molecule_source_sha256/'+molecule] = np.bytes_(sha(path))


def molecular_node_bank(directory, main_data, output, n):
    """Explicit coarse-cell projection, NOT point sampling of narrow lines.

    Node supports are midpoint/Voronoi energy intervals, whose widths equal
    the existing trapezoidal dE prior. Integrate the piecewise-linear Leiden
    cross section exactly over each support. This preserves its band integral
    for a flat photon spectrum, but is not line-resolved attenuation.
    """
    edges = np.asarray(EDGES)
    nodes = np.array([lo*np.exp(np.log(hi/lo)*np.arange(n)/(n-1))
                      for lo, hi in zip(edges[:-1], edges[1:])])
    nodes[:, 0] = edges[:-1]; nodes[:, -1] = edges[1:]
    bounds = np.concatenate((nodes[:, :1], .5*(nodes[:, :-1]+nodes[:, 1:]), nodes[:, -1:]), axis=1)
    width = np.diff(bounds, axis=1)
    mappings, absorption, dissociation, source_hash = [], [], [], {}
    # Empirical CHIMES channels keep their existing Habing-shape closure;
    # this is NOT a measured molecular cross-section reconstruction.
    evaluation = nodes.copy()
    evaluation[:-1, -1] = np.nextafter(edges[1:-1], edges[:-2])
    habing = np.where((evaluation >= 6.) & (evaluation < 13.6),
                      evaluation*1.60217657e-12/(3e10*5.29e-14), 0.)
    with h5py.File(main_data) as f:
        for group in ('photodissoc_group1', 'photodissoc_group2'):
            for reactant, products, rate in zip(f[group+'/reactants'][:], f[group+'/products'][:], f[group+'/rates'][:]):
                mappings.append([int(reactant), *map(int, products)])
                absorption.append(float(rate)*habing)
                dissociation.append(float(rate)*habing)
        for species, group in [('H2','H2_photodissoc'), ('CO','CO_photodissoc')]:
            path = directory/f'{species}.h5'
            source_hash[species] = sha(path)
            with h5py.File(path) as raw:
                wavelength = raw['wavelength'][:]
                if not np.all(np.isfinite(wavelength)) or not np.all(np.diff(wavelength)>0) or wavelength[0]<=0:
                    raise ValueError('invalid Leiden wavelength grid')
                energy = (1239.8419843320026/wavelength)[::-1]
                diss = raw['photodissociation'][:][::-1]
                nonion = (raw['photoabsorption'][:]-raw['photoionisation'][:])[::-1]
                nonion = np.maximum(diss, nonion)  # same admitted nonionizing channel definition
            if not np.all(np.isfinite([diss, nonion])) or np.any(diss<0):
                raise ValueError('invalid Leiden cross sections')
            def project(cross):
                slope = np.diff(cross)/np.diff(energy)
                integral = np.concatenate(([0.], np.cumsum(.5*(cross[:-1]+cross[1:])*np.diff(energy))))
                x = np.clip(bounds, energy[0], energy[-1])
                i = np.clip(np.searchsorted(energy, x, side='right')-1, 0, len(energy)-2)
                dx = x-energy[i]
                primitive = integral[i]+cross[i]*dx+.5*slope[i]*dx**2
                result = np.diff(primitive, axis=1)/width
                # Nonnegative segment integrals must not become subtraction
                # noise in transparent supports. Recompute only those locally.
                for g,j in zip(*np.where(result<0)):
                    lo, hi = max(bounds[g,j],energy[0]), min(bounds[g,j+1],energy[-1])
                    if lo>=hi:
                        result[g,j]=0.
                    else:
                        grid=np.r_[lo,energy[(energy>lo)&(energy<hi)],hi]
                        result[g,j]=np.trapezoid(np.interp(grid,energy,cross),grid)/width[g,j]
                return result
            ab, di = project(nonion), project(diss)
            # Store the two physically exclusive nonionizing outcomes. This
            # avoids subtracting nearly equal projected absorption/dissociation.
            fluorescence = project(nonion-diss)
            ab = di+fluorescence
            mappings.append([int(f[group+'/reactants'][0]), *map(int,f[group+'/products'][0])])
            absorption.append(ab); dissociation.append(di)
    target = output/'molecular_nodes.h5'
    with h5py.File(target,'x') as f:
        f.attrs['schema']='snrt_chimes_molecular_nodes_v1'
        f.attrs['main_sha256']=MAIN_SHA
        f.attrs['closure']='30 Habing-shape channels; H2/CO midpoint-support integral projection; no shielding applied'
        f['node_ev']=nodes; f['edges_ev']=edges; f['support_edges_ev']=bounds
        f['reactions']=np.asarray(mappings,dtype=np.int32)
        f['absorption_cm2']=np.asarray(absorption)
        f['dissociation_cm2']=np.asarray(dissociation)
        for species,digest in source_hash.items(): f.attrs[species+'_sha256']=digest
    manifest={'schema':'snrt_chimes_molecular_nodes_v1','main_sha256':MAIN_SHA,
              'table_sha256':sha(target),'source_sha256':source_hash,'nodes':n,'channels':len(mappings),
              'limitations':['Line-integral-preserving coarse projection, NOT line-resolved transport.',
                             'Thirty empirical channels retain Habing-shape approximation.',
                             'No shielding, heating partition or live admission supplied by this data file.']}
    with (output/'molecular_manifest.json').open('x') as f:
        json.dump(manifest,f,indent=2); f.write('\n')
    print(json.dumps(manifest,indent=2))


def slope_for_mean(lo, hi, target):
    # Shift the exponent so root finding stays bounded even for steep bins.
    def mean(slope):
        pivot = np.log(hi if slope > 0 else lo)
        den = quad(lambda x: np.exp(slope*(x-pivot)), np.log(lo), np.log(hi))[0]
        num = quad(lambda x: np.exp(x+slope*(x-pivot)), np.log(lo), np.log(hi))[0]
        return num/den
    slope = brentq(lambda s: mean(s)-target, -100., 100., xtol=1e-12)
    if abs(mean(slope)/target-1) > 1e-10:
        raise ValueError('photon energy moment did not close')
    return slope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chimes-tools', type=Path, required=True)
    parser.add_argument('--main-data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--leiden-data', type=Path)
    parser.add_argument('--band-nodes', type=int, choices=(0, 64, 128, 256), default=0,
                        help='prepare a separate shell-resolved native atomic bank; do not replace grey tables')
    args = parser.parse_args()
    tools = args.chimes_tools.resolve()
    main_data = args.main_data.resolve()
    if sha(main_data) != MAIN_SHA:
        raise ValueError('unadmitted main reaction table')
    if args.band_nodes:
        build_shell_bank(tools, main_data, args.output.resolve(), args.band_nodes)
        if args.leiden_data is not None:
            molecular_node_bank(args.leiden_data.resolve(),main_data,args.output.resolve(),args.band_nodes)
        return
    if args.leiden_data is None:
        parser.error('--leiden-data is required for the original grey tables')
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    generator = tools/'generate_cross_sections/generate_cross_sections.py'
    manifest = {'closure': 'photon_mean_matched_powerlaw_v1',
                'edges_ev': EDGES, 'means_ev': MEANS, 'main_sha256': MAIN_SHA,
                'generator_sha256': sha(generator), 'wrapper_sha256': sha(Path(__file__)),
                'photoelectron_closure': 'primary_shell_electron_v1',
                'molecule_source_sha256': {name: sha(args.leiden_data/f'{name}.h5') for name in ['H2','CO']},
                'groups': []}
    for index, (lo, hi, energy) in enumerate(zip(EDGES, EDGES[1:], MEANS), 1):
        slope = slope_for_mean(lo, hi, energy)
        spectrum = output/f'group_{index:02d}.spectrum'
        x = np.linspace(-6., 12., 181)
        with spectrum.open('x') as handle:
            np.savetxt(handle, np.column_stack((x, -20.+slope*(x-np.log10(energy/13.6)))))
        table = output/f'group_{index:02d}.hdf5'
        parameters = {'chimes_main_data_path': main_data,
                      'cross_sections_data_path': tools/'generate_cross_sections/data',
                      'spectrum_file': spectrum, 'output_file': table,
                      'spectrum_shape': 'SNRT_mean_matched_powerlaw_v1',
                      'version_date': '2026-09-09', 'rt_coupling_flag': 1,
                      'rt_HII_recombination_flag': 0, 'rt_HeII_recombination_flag': 0,
                      'rt_E_min_eV': lo, 'rt_E_max_eV': hi}
        paramfile = output/f'group_{index:02d}.params'
        with paramfile.open('x') as handle:
            for key, value in parameters.items():
                handle.write(f'{key} {value}\n')
        with (output/f'group_{index:02d}.log').open('x') as handle:
            subprocess.run([sys.executable, str(generator), str(paramfile)],
                           stdout=handle, stderr=subprocess.STDOUT, check=True)
        primary_electron_moments(tools, main_data, table, lo, hi, slope)
        molecular_moments(args.leiden_data.resolve(), table, lo, hi, slope)
        with h5py.File(table, 'r+') as data:
            data['Header/snrt_group_edges_ev'] = [lo, hi]
            data['Header/snrt_mean_energy_ev'] = energy
            data['Header/snrt_powerlaw_slope'] = slope
            def check(name, obj):
                if isinstance(obj, h5py.Dataset) and obj.dtype.kind in 'fc':
                    if not np.isfinite(obj[...]).all():
                        raise ValueError(f'nonfinite table: {index}:{name}')
            data.visititems(check)
        manifest['groups'].append({'index': index, 'path': table.name,
                                   'sha256': sha(table), 'powerlaw_slope': slope})
        print(f'group {index}: {lo:g}--{hi:g} eV, mean {energy:g} eV, slope {slope:g}', flush=True)
    with (output/'manifest.json').open('x') as handle:
        json.dump(manifest, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()

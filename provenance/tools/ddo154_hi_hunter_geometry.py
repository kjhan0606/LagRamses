#!/usr/bin/env python3
"""Preliminary DDO154 H I radial inventory, retaining blank-pixel ambiguity."""
import csv
import io
import json
import os
import re
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area

raw = Path(os.environ['RAR_HI_IMAGE'])
sparc = Path(os.environ['RAR_SPARC_ZIP'])
out = Path(os.environ['RAR_ANALYSIS_OUTPUT'])
distance_mpc = 4.04  # SPARC DDO154 rotmod distance
ra_deg, dec_deg = 193.5258333333333, 27.1505555555556  # Hunter+2012 optical centre
inc_deg, pa_deg = 65.2, 46.0  # Hunter+2012 optical geometry
with zipfile.ZipFile(sparc) as archive:
    rot = np.loadtxt(io.BytesIO(archive.read('DDO154_rotmod.dat')))
if rot.ndim != 2 or rot.shape[1] != 8 or np.any(np.diff(rot[:,0]) <= 0):
    raise ValueError('unexpected SPARC rotation curve')
with fits.open(raw) as hdus:
    image = np.squeeze(np.asarray(hdus[0].data, dtype=np.float64))
    header = hdus[0].header
    if image.shape != (1024,1024) or not np.isfinite(image).all():
        raise ValueError('unexpected image')
    if str(header.get('BUNIT','')).strip().upper() != 'JY/B*M/S':
        raise ValueError('unexpected intensity units')
    hist = '\n'.join(str(x) for x in header.get('HISTORY',[]))
    match = re.search(r'CLEAN BMAJ=\s*([\d.E+-]+)\s+BMIN=\s*([\d.E+-]+)\s+BPA=\s*([\d.E+-]+)',hist)
    if not match or 'HISTORY' in match.group(0):
        raise ValueError('missing AIPS clean beam record')
    major, minor, bpa = map(float, match.groups())
    # Convention assumed here: AIPS beam axes are degrees; verify against source documentation.
    major_as, minor_as = 3600*major, 3600*minor
    nu_ghz = float(header['RESTFREQ'])/1e9
    pixel_arcsec2 = float(proj_plane_pixel_area(WCS(header).celestial))*3600**2
    beam_arcsec2 = np.pi/(4*np.log(2))*major_as*minor_as
    yy, xx = np.indices(image.shape,dtype=float)
    ra,dec = WCS(header).celestial.pixel_to_world_values(xx,yy)
dra = (ra-ra_deg+180)%360-180
east = dra*np.cos(np.deg2rad(dec_deg))*3600
north = (dec-dec_deg)*3600
pa = np.deg2rad(pa_deg)
major_axis = east*np.sin(pa)+north*np.cos(pa)
minor_axis = east*np.cos(pa)-north*np.sin(pa)
r_kpc = np.hypot(major_axis,minor_axis/np.cos(np.deg2rad(inc_deg)))*distance_mpc*1000/206264.806247
# Jy/beam * m/s -> K km/s, optically thin N_HI -> H I mass surface density (no helium).
temperature_factor = 1.222e6/(nu_ghz**2*major_as*minor_as)/1000
sigma_factor = temperature_factor*(1.823e18/1.248e20)*np.cos(np.deg2rad(inc_deg))
valid_proxy = image != 0.0  # REMAG replaced blanks with 0; actual valid-zero pixels are indistinguishable
edges = np.arange(0,8.00001,0.25)
rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    ring = (r_kpc>=lo)&(r_kpc<hi)
    present=ring&valid_proxy
    n=int(ring.sum()); n_nonzero=int(present.sum())
    rows.append({'r_inner_kpc':float(lo),'r_outer_kpc':float(hi),'pixels':n,
                 'nonzero_pixels':n_nonzero,'nonzero_fraction':n_nonzero/n if n else None,
                 'mean_hi_sigma_zero_filled':float(image[ring].mean()*sigma_factor) if n else None,
                 'mean_hi_sigma_nonzero_only':float(image[present].mean()*sigma_factor) if n_nonzero else None})
with (out/'radial_probe.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
sum_flux_jykms = float(image.sum()*pixel_arcsec2/beam_arcsec2/1000)
result={'geometry_source':'Hunter et al. 2012 Table 6 optical','distance_mpc':distance_mpc,'centre_deg':[ra_deg,dec_deg],
        'inclination_deg':inc_deg,'position_angle_deg':pa_deg,'beam_arcsec_assumed':[major_as,minor_as],
        'beam_position_angle_deg':bpa,'pixel_area_arcsec2':pixel_arcsec2,
        'sparc_radius_max_kpc':float(rot[-1,0]),'sparc_radius_count':len(rot),
        'nonzero_image_pixels':int(valid_proxy.sum()),'zero_image_pixels':int((~valid_proxy).sum()),
        'max_nonzero_deprojected_radius_kpc':float(r_kpc[valid_proxy].max()),
        'sum_flux_jy_kms_zero_filled':sum_flux_jykms,
        'optically_thin_HI_mass_msun_zero_filled':2.356e5*distance_mpc**2*sum_flux_jykms,
        'limits':'Diagnostic only: zero pixels include replaced blanks and may include true zeros, so both radial averages are biased differently; no noise, missing flux, orientation, beam-unit or inclination uncertainty propagated. No helium, molecular gas, 3D source, gravitational force or QUMOND calculation.'}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

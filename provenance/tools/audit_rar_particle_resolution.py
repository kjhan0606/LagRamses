#!/usr/bin/env python3
"""Measure native self-force resolution changes with a fixed physical source."""
import json,re,subprocess,sys
from pathlib import Path
import numpy as np
from rar_native_particles import epoch

root=Path(sys.argv[1]);old=json.loads((root/'original_self_pair.json').read_text());cases={}
for case in sorted(root.glob('L*_*')):
 if not case.is_dir():continue
 match=re.fullmatch(r'L([456])_(newtonian|qumond)',case.name)
 if match is None:raise ValueError('unknown resolution case')
 level=int(match[1]);model=match[2]
 namelist=(case/'effective.nml').read_text()
 def scalar(name):
  values=re.findall(r'(?im)^\s*'+name+r'\s*=\s*([^!\n]+)',namelist)
  if len(values)!=1:raise ValueError('ambiguous refinement setting '+name)
  return float(values[0].strip().lower().replace('d','e'))
 exponent=scalar('exp_refine');diameter=scalar('r_refine')
 refinement='fixed_box' if exponent==10 else 'sphere'
 if refinement=='sphere' and (exponent!=2 or diameter!=.25):raise ValueError('unexpected sphere geometry')
 if refinement=='fixed_box' and diameter!=.375-2*2.**(-level):raise ValueError('unexpected fixed-box diameter')
 leaves=np.loadtxt(case/'composite_epoch_00000_leaf.tsv');faces=np.loadtxt(case/'composite_faces.tsv')
 rows,rho,report=epoch(case/'particles_epoch_00000.tsv',leaves)
 if len(rows)!=1 or report['cloud_width']!=.0625 or rows[0,3]!=.002 or not np.allclose(rows[0,4:7],[.6873,.51,.52],rtol=0,atol=1e-15):raise ValueError('physical cloud changed across meshes')
 if np.any(rows[:,7:10]!=0) or np.max(abs(leaves[:,6]-1-rho))>1e-12:raise ValueError('initial gas/particle source mismatch')
 subprocess.run([sys.executable,str(root/'verify.py'),str(root),str(root/f'verify_{case.name}.json'),
  '--case',case.name,'--leaf-file','composite_epoch_00000_leaf.tsv','--snapshots','0','--evolving','--particles'],check=True)
 l=faces[:,0].astype(int)-1;r=faces[:,1].astype(int)-1;axis=faces[:,2].astype(int)-1
 center=leaves[:,2:5];width=leaves[:,5];vol=width**3;mask=width[l]!=width[r]
 fine=leaves[:,1]==level+1
 fine_volume=float(np.sum(vol[fine]))
 if refinement=='fixed_box':
  fine_lo=center[fine]-width[fine,None]/2;fine_hi=center[fine]+width[fine,None]/2
  if np.any(fine_lo<.3125) or np.any(fine_hi>.6875) or fine_volume!=.375**3:
   raise ValueError('actual refined geometry is not the fixed physical box')
 l,r,axis,area=l[mask],r[mask],axis[mask],faces[mask,3]
 if len(l)==0:raise ValueError('refinement interface missing')
 facecenter=np.where((width[l]<width[r])[:,None],center[l],center[r])
 facecenter[np.arange(len(l)),axis]=np.mod(center[l,axis]+width[l]/2,1.)
 delta=rows[0,4:7]-facecenter;delta-=np.rint(delta)
 outside=np.maximum(abs(delta)-np.sqrt(area)[:,None]/2,0.)
 outside[np.arange(len(l)),axis]=abs(delta[np.arange(len(l)),axis])
 nearest=float(np.min(np.linalg.norm(outside,axis=1)))
 accel=rows[0,10:13];reference=np.array(old['cases'][model+'_a']['acceleration'][0])
 baseline_difference=float(np.max(abs(accel-reference)))
 if refinement=='sphere' and level==4 and baseline_difference>1e-10:raise ValueError('default width changed established force')
 fine_fraction=float(np.sum((rho*vol)[fine])/.002)
 if refinement=='fixed_box' and (abs(nearest-.0002)>1e-12 or abs(fine_fraction-.5032)>1e-12):
  raise ValueError('particle/interface geometry changed across fixed-box meshes')
 report.update(base_level=level,leaves=len(leaves),subfaces=len(faces),acceleration=accel.tolist(),
  refinement=refinement,requested_refinement_diameter=diameter,refined_volume=fine_volume,
  self_over_Gm_s2=float(np.linalg.norm(accel)/(.002/.0625**2)),
  ratio_to_original_amr_self_force=float(np.linalg.norm(accel)/np.linalg.norm(reference)),
  fine_level_mass_fraction=fine_fraction,
  nearest_refinement_subface_distance=nearest,
  baseline_acceleration_max_difference=baseline_difference)
 cases[case.name]=report
if not cases:raise ValueError('no resolution results')
for name,report in cases.items():
 model=name.split('_',1)[1];base=cases.get('L4_'+model)
 if base is not None:
  report['self_force_ratio_to_this_series_L4']=report['self_over_Gm_s2']/base['self_over_Gm_s2']
result={'schema':'rar_native_fixed_cloud_resolution_v1','cases':cases,
 'scope':'Same physical single-cloud mass, width and position on different native two-level meshes. Initial forces; one explicit evolution step without snapshots.',
 'limits':'Sphere mode has a moving rasterized boundary. Fixed-box mode verifies the actual physical refinement volume, particle-interface distance and refined mass fraction at each resolution. Neither mode by itself is a galaxy softening calibration or proof of exact momentum conservation.'}
(root/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(cases,indent=2))

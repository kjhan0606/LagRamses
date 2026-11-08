"""Match accepted thermal updates to actual nonuniform AMR cooling calls."""
import math
import numpy as np

def thermal_accounting(text,reference,steps,dt,unit_t,record_parser):
 intervals=np.full(steps,dt) if np.ndim(dt)==0 else np.asarray(dt,float)
 if intervals.shape!=(steps,) or np.any(intervals<=0) or not np.isfinite(intervals).all():raise ValueError('invalid thermal timestep sequence')
 times=np.cumsum(intervals)
 rows=record_parser(text,'RAR_COOLING_COMMIT_V2',11)
 if rows.ndim!=2 or rows.shape[1]!=11:raise ValueError('no evolving thermal updates')
 n=len(reference);lookup={int(v):i for i,v in enumerate(reference[:,0])};levels=reference[:,1].astype(int);volumes=reference[:,5]**3
 inventories={};lines=iter(text.splitlines())
 for line in lines:
  if 'RAR_GLOBAL_INVENTORY_V2' not in line:continue
  tokens=line.split('RAR_GLOBAL_INVENTORY_V2',1)[1].split()
  while len(tokens)<6:tokens.extend(next(lines).split())
  nv=int(tokens[5]);target=6+nv+12
  while len(tokens)<target:tokens.extend(next(lines).split())
  if len(tokens)!=target:raise ValueError('bad AMR inventory shape')
  label=tokens[0];v=np.array([float(x.replace('D','E')) for x in tokens[1:]])
  if not np.isfinite(v).all():raise ValueError('nonfinite AMR inventory')
  step=int(v[0]);lev=int(v[3]);key=(step,lev,label)
  if key in inventories:raise ValueError('duplicate cooling call identity')
  if not 0<=step<steps or abs(v[1]-times[step])>1e-16 or v[2]!=1 or int(v[4])!=nv:raise ValueError('inventory source clock/UV epoch mismatch')
  inventories[key]=(v[5:],nv)
 expected_calls={(k,int(lev),phase) for k in range(steps) for lev in set(levels) for phase in ('before_cooling','after_cooling')}
 if set(inventories)!=expected_calls or len(rows)!=steps*n:raise ValueError('missing/extra AMR thermal call coverage')
 increments=[];effective=[];calls=[]
 for k in range(steps):
  batch=rows[rows[:,3]==k]
  if len(batch)!=n or len(np.unique(batch[:,2]))!=n:raise ValueError('thermal leaf coverage mismatch')
  ix=np.array([lookup[int(c)] for c in batch[:,2]]);batch=batch[np.argsort(ix)]
  if not np.array_equal(batch[:,1],levels) or np.any(batch[:,0]!=1):raise ValueError('thermal leaf owner/level mismatch')
  if not np.allclose(batch[:,4],times[k],rtol=0,atol=1e-16) or np.any(batch[:,5]!=1):raise ValueError('thermal time/background mismatch')
  if not np.allclose(batch[:,6],intervals[k]*unit_t,rtol=1e-14,atol=0) or not np.array_equal(batch[:,7],volumes) or np.any(batch[:,10]<=0):raise ValueError('thermal interval/volume/unit mismatch')
  increments.append(batch[:,9]/volumes);effective.append(batch[:,8]/volumes)
  for lev in sorted(set(levels)):
   select=levels==lev;part=batch[select];before,nv=inventories[(k,int(lev),'before_cooling')];after,nv2=inventories[(k,int(lev),'after_cooling')]
   if nv!=nv2 or before[nv]!=n or after[nv]!=n:raise ValueError('whole-call owned leaf census mismatch')
   unchanged=np.r_[np.arange(4),np.arange(5,nv)]
   if not np.array_equal(before[unchanged],after[unchanged]):raise ValueError('thermal call changed mass/momentum/passives')
   rhs=math.fsum(part[:,8]);accepted=math.fsum(part[:,9]);actual=float(after[4]-before[4]);residual=actual-accepted
   scale=max(abs(before[4]),abs(after[4]),abs(accepted));tol=256*math.ulp(scale)+1e-10*abs(accepted)
   if abs(residual)>tol:raise ValueError('whole cooling call does not match accepted thermal receipt')
   calls.append({'step':k,'level':int(lev),'time':float(times[k]),'leaf_records':int(select.sum()),'owned_volume':float(volumes[select].sum()),'effective_rhs_code':rhs,'accepted_energy_code':accepted,'whole_call_energy_code':actual,'whole_call_residual_code':residual,'tolerance':tol})
 return {'accepted':increments,'effective':effective,'statistics':{'calls':calls,'records':len(rows),'heating_records':int(np.count_nonzero(rows[:,9]>0)),'cooling_records':int(np.count_nonzero(rows[:,9]<0)),'effective_rhs_code':math.fsum(rows[:,8]),'accepted_energy_code':math.fsum(rows[:,9]),'max_whole_call_residual_code':max(abs(v['whole_call_residual_code']) for v in calls)}}

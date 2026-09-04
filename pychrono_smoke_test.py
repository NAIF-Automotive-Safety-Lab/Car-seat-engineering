#!/usr/bin/env python3
import json, platform, sys, time, traceback

def v(obj):
 try: return obj.GetXYZ() if hasattr(obj,'GetXYZ') else [obj.x,obj.y,obj.z]
 except Exception: return str(obj)

def main():
 out={'timestamp_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'engine':'pychrono','host':platform.platform(),'arch':platform.machine()}
 try:
  import pychrono.core as chrono
  out['import']='PASS'
  out['version']=getattr(chrono,'__version__','UNKNOWN')
  sysm=chrono.ChSystemNSC(); out['system_created']=True
  a=chrono.ChBody(); a.SetFixed(True); a.SetMass(1); sysm.AddBody(a)
  b=chrono.ChBody(); b.SetMass(1); b.SetPos(chrono.ChVector3d(1,0,0)); sysm.AddBody(b)
  link=chrono.ChLinkLockRevolute(); link.Initialize(a,b,chrono.ChFramed(chrono.ChVector3d(0,0,0))); sysm.AddLink(link)
  out['joint_created']=True
  for _ in range(10): sysm.DoStepDynamics(1e-3)
  out['simulation_time']=sysm.GetChTime(); out['position']=v(b.GetPos()); out['velocity']=v(b.GetPosDt()); out['reaction']=str(link.GetReaction1())
  # GetReaction1 being callable and returning a wrench is the reaction extraction gate.
  out['reaction_extracted']=True
  out['status']='PASS'
 except Exception as e:
  out['status']='FAIL'; out['error']=f'{type(e).__name__}: {e}'; out['traceback']=traceback.format_exc()
 open('artifacts/pychrono_smoke.json','w').write(json.dumps(out,indent=2))
 print('PYCHRONO_IMPORT='+out.get('import','FAIL'))
 print('PYCHRONO_SMOKE='+out['status'])
 return 0 if out['status']=='PASS' else 7
if __name__=='__main__': sys.exit(main())

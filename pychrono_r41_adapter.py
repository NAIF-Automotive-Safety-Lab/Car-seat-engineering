from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
EXPECTED_SHA256='fbe6b17cdbf728a2e47963e567e12eeceb1352a36e719e7d1c55cc5f712a0a68'
EXPECTED_SOLIDS=62
ROOT=Path(__file__).resolve().parents[1]
STEP_PATH=ROOT/'geometry'/'T_OCS_V0_R4_1_SHOULDER_REPAIR.step'
BODY_MAP=ROOT/'manifests'/'V1_SR11_BODY_MAPPING.json'
JOINT_MAP=ROOT/'joints'/'V1_SR11_JOINT_MAPPING.json'
CONTACT_MAP=ROOT/'contacts'/'V1_SR11_CONTACT_MAPPING.json'
ASSUMPTIONS=ROOT/'parameters'/'V1_SR11_ASSUMPTIONS.json'

def sha256(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()

def static_validate()->dict:
    checks=[]
    if not STEP_PATH.is_file(): return {'status':'FAIL','failure':'STEP_MISSING','checks':checks}
    got=sha256(STEP_PATH); checks.append({'check':'sha256','expected':EXPECTED_SHA256,'actual':got,'pass':got==EXPECTED_SHA256})
    if got!=EXPECTED_SHA256: return {'status':'FAIL','failure':'R4.1_SHA_MISMATCH','checks':checks}
    from OCP.STEPControl import STEPControl_Reader
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopoDS import TopoDS
    r=STEPControl_Reader(); st=r.ReadFile(str(STEP_PATH)); r.TransferRoots(); shape=r.OneShape()
    solids=[]; ex=TopExp_Explorer(shape,TopAbs_SOLID)
    while ex.More(): solids.append(TopoDS.Solid_s(ex.Current())); ex.Next()
    read_pass=('RetDone' in str(st))
    checks.append({'check':'step_read','actual_status':str(st),'pass':read_pass})
    checks.append({'check':'solids','expected':EXPECTED_SOLIDS,'actual':len(solids),'pass':len(solids)==EXPECTED_SOLIDS})
    bm=json.loads(BODY_MAP.read_text())
    ids=sorted(i for b in bm['bodies'] for i in b['step_solid_indices'])
    checks.append({'check':'mapping_coverage','expected_count':62,'actual_count':len(ids),'unique_count':len(set(ids)),'pass':ids==list(range(1,63))})
    jm=json.loads(JOINT_MAP.read_text()); cm=json.loads(CONTACT_MAP.read_text()); am=json.loads(ASSUMPTIONS.read_text())
    checks.append({'check':'joint_refs','pass':all(j.get('body_A') and j.get('body_B') for j in jm['joints'])})
    checks.append({'check':'contact_refs','pass':all(c.get('body_A') and c.get('body_B') for c in cm['contacts'])})
    allowed={'UNKNOWN','ASSUMED','SOURCE_VERIFIED','DERIVED'}
    checks.append({'check':'no_silent_defaults','pass':all(x.get('status') in allowed for x in am['parameters'])})
    return {'status':'PASS' if all(x.get('pass',False) for x in checks) else 'FAIL','validation_class':'STATIC_ADAPTER_VALIDATION','solver_execution_performed':False,'checks':checks}

def build_chrono_system():
    try:
        import pychrono.core as chrono
        import pychrono.cascade as cascade
    except Exception as e:
        raise RuntimeError(f'PYCHRONO_REQUIRED_BUT_UNAVAILABLE: {type(e).__name__}: {e}')
    s=static_validate()
    if s['status']!='PASS': raise RuntimeError('STATIC_VALIDATION_FAILED')
    unknown=[p['parameter'] for p in json.loads(ASSUMPTIONS.read_text())['parameters'] if p['status']=='UNKNOWN']
    if unknown: raise RuntimeError('MBD_BUILD_BLOCKED_BY_UNRESOLVED_PARAMETERS: '+', '.join(unknown))
    doc=cascade.ChCascadeDoc()
    if not doc.LoadSTEP(str(STEP_PATH)): raise RuntimeError('CHRONO_CASCADE_STEP_LOAD_FAILED')
    system=chrono.ChSystemNSC()
    return system,doc

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--static-validate',action='store_true'); ap.add_argument('--build',action='store_true')
    a=ap.parse_args()
    if a.static_validate:
        print(json.dumps(static_validate(),indent=2)); return 0
    if a.build:
        try:
            system,doc=build_chrono_system(); print(json.dumps({'status':'CHRONO_SYSTEM_CREATED','system_type':type(system).__name__,'step_loaded':True},indent=2)); return 0
        except Exception as e:
            print(json.dumps({'status':'BLOCKED','error':str(e)},indent=2)); return 2
    return 1
if __name__=='__main__': sys.exit(main())

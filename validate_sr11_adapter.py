#!/usr/bin/env python3
import json, sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sr=root/'SR11'
REQUIRED={'BASE','RAIL_L','RAIL_R','CARRIAGE','ABSORBER_L','ABSORBER_R','SEATBACK','ROT_CTRL_L','ROT_CTRL_R','LOCK_L','LOCK_R','ENDSTOPS','HINGE','ANTI_RACK','SHOULDER_ANCHOR','SHOULDER_GUSSET','SHOULDER_ROOT','SHOULDER_BRACE','SHOULDER_FRAME_BRACKET'}

def main():
 p=root/'artifacts'/'sr11_validation.json'; p.parent.mkdir(parents=True,exist_ok=True)
 checks=[]
 bm=json.loads((sr/'manifests/V1_SR11_BODY_MAPPING.json').read_text()); bodies=bm['bodies']; ids=sorted(i for b in bodies for i in b['step_solid_indices'])
 checks.append({'name':'mapping_coverage','pass':ids==list(range(1,63)) and len(set(ids))==62,'count':len(ids),'unique':len(set(ids))})
 body_ids={b['mbd_body_id'] for b in bodies}; checks.append({'name':'required_model_elements','pass':REQUIRED.issubset(body_ids),'missing':sorted(REQUIRED-body_ids)})
 jm=json.loads((sr/'joints/V1_SR11_JOINT_MAPPING.json').read_text()); checks.append({'name':'joint_references','pass':all(j.get('body_A') and j.get('body_B') for j in jm['joints'])})
 cm=json.loads((sr/'contacts/V1_SR11_CONTACT_MAPPING.json').read_text()); checks.append({'name':'contact_references','pass':all(c.get('body_A') and c.get('body_B') for c in cm['contacts'])})
 am=json.loads((sr/'parameters/V1_SR11_ASSUMPTIONS.json').read_text()); checks.append({'name':'no_silent_defaults','pass':all(x.get('status') in {'UNKNOWN','ASSUMED','SOURCE_VERIFIED','DERIVED'} for x in am['parameters'])})
 checks.append({'name':'endstop_mapping','pass':body_ids.__contains__('ENDSTOPS') and any('ENDSTOP' in str(c.get('type','')).upper() for c in cm['contacts'])})
 checks.append({'name':'lock_mapping','pass':{'LOCK_L','LOCK_R'}.issubset(body_ids)})
 checks.append({'name':'bilateral_absorber_mapping','pass':{'ABSORBER_L','ABSORBER_R'}.issubset(body_ids)})
 # Verify adapter exposes required callable surface before allowing Gate 7 PASS.
 text=(sr/'adapter/pychrono_r41_adapter.py').read_text()
 required_calls=['static_validate','build_chrono_system']
 for name in required_calls: checks.append({'name':f'adapter_callable_{name}','pass':f'def {name}(' in text})
 result={'status':'PASS' if all(c['pass'] for c in checks) else 'FAIL','validation_class':'STATIC_SR11_ADAPTER','checks':checks,'solver_execution_performed':False}
 p.write_text(json.dumps(result,indent=2))
 for c in checks: print(c['name'].upper()+'='+('PASS' if c['pass'] else 'FAIL'))
 print('SR11_ADAPTER='+result['status'])
 return 0 if result['status']=='PASS' else 5
if __name__=='__main__': sys.exit(main())

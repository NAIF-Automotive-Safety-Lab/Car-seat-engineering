import json
from pathlib import Path
oldp=Path('/home/ubuntu/rc006-audit-work/package/01_MASTER/RC006_MASTER_EVIDENCE_REGISTER.json')
newp=Path('/home/ubuntu/rc006-repaired-audit-work/package/01_MASTER/RC006_MASTER_EVIDENCE_REGISTER.json')
old=[x for x in json.loads(oldp.read_text())['records'] if x['ID']=='VEH-INPUT']
new=[x for x in json.loads(newp.read_text())['records'] if x['ID'].startswith('VEH-')]
print('old_vehicle_records',len(old)); print('new_vehicle_records',len(new)); print('new_ids',[x['ID'] for x in new])
print('new_records')
for x in new: print(json.dumps(x,sort_keys=True))
# All old duplicate records should have identical invariant fields; compare each new to old template.
inv=['DOMAIN','CURRENT_VALUE','UNIT','VALUE_CLASS','SOURCE','SOURCE_TYPE','AUTHORITY','TRACEABILITY','APPLICABILITY','STATUS','BLOCKING','NEXT_ACTION']
template={k:old[0].get(k) for k in inv}
issues=[]
for x in new:
 for k in inv:
  if x.get(k)!=template[k]: issues.append((x['ID'],k,template[k],x.get(k)))
print('invariant_issues',issues)
print('all_values_null',all(x.get('CURRENT_VALUE') is None for x in new))
print('all_external_source_required',all(x.get('STATUS')=='EXTERNAL_SOURCE_REQUIRED' for x in new))
print('all_blocking',all(x.get('BLOCKING') is True for x in new))
